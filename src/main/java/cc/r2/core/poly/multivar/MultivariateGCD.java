package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.univar.PolynomialGCD.PolynomialGCD;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateGCD {
    private MultivariateGCD() {}

    private static void checkModGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {

        if (!(a.domain instanceof ModularDomain))
            throw new IllegalArgumentException();

    }

    @SuppressWarnings("unchecked")
    static MultivariatePolynomial<BigInteger> denseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {
        //todo choose first variable correctly

        BigInteger domainSize = a.domain.size();
        if (domainSize.isInt() && Math.min(a.degree(), b.degree()) > domainSize.intValue())
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");

        int evaluationStackLimit = domainSize.isInt() ? domainSize.intValue() : -1;
        int[] aDegrees = a.degrees(),
                bDegrees = b.degrees(),
                degreeBounds = new int[a.nVariables];
        for (int i = 0; i < a.nVariables; i++)
            degreeBounds[i] = Math.min(aDegrees[i], bDegrees[i]);

        MultivariatePolynomial<BigInteger> result = denseModularGCD(a, b, PrivateRandom.getRandom(), a.nVariables - 1, degreeBounds, evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return result;
    }

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> denseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("input is zero");

        MultivariatePolynomial<BigInteger> factory = a;
        if (a.isConstant() || b.isConstant())
            return factory.createOne();

        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            bMutablePolynomialZp gcd = PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b));
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        Domain<BigInteger> domain = factory.domain;

        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        MultivariatePolynomial<bMutablePolynomialZp>
                aConv = convertZp(a, variable),
                bConv = convertZp(b, variable);

        // content of a and b in Zp[x_k]
        bMutablePolynomialZp
                aCont = PolynomialGCD(aConv.data.values()),
                bCont = PolynomialGCD(bConv.data.values());

        // normalize a and b
        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(a, asMultivariate(aCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        a = qd[0];

        qd = divideAndRemainder(b, asMultivariate(bCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        b = qd[0];

        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = PolynomialGCD(aCont, bCont),
                lcGCD = PolynomialGCD(aConv.lc(), bConv.lc());

//        List<bMutablePolynomialZp> evalTest = new ArrayList<>();
//        for (int i = 0; i < variable; ++i)
//            evalTest.add(PolynomialGCD(convertZp(a, i).lc(), convertZp(b, i).lc()));

        int prevVarExponent = degreeBounds[variable - 1];
        Interpolation<BigInteger> interpolation = null;
        MultivariatePolynomial<BigInteger> previousInterpolation;
        Set<BigInteger> evaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // do division check (last chance)
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            BigInteger randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;

            evaluationStack.add(randomPoint);

            BigInteger lcVal = lcGCD.evaluate(randomPoint);
            if (lcVal.isZero())
                continue;

//            for (bMutablePolynomialZp lc : evalTest)
//                if (lc.evaluate(randomPoint).isZero())
//                    continue main;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<BigInteger> cVal = denseModularGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            // normalize gcd
            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            if (currExponent < prevVarExponent) {
                //better degree bound detected => start over
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
                continue;
            }

            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                continue;
            }

            previousInterpolation = interpolation.getInterpolatingPolynomial();
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] == 0
                    || (previousInterpolation != null && previousInterpolation.equals(interpolation.getInterpolatingPolynomial()))) {
                MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    private static MultivariatePolynomial<BigInteger> doDivisionCheck(
            MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, bMutablePolynomialZp contentGCD,
            Interpolation<BigInteger> interpolation, int variable) {
        if (interpolation == null)
            return null;
        MultivariatePolynomial<BigInteger> interpolated =
                fromZp(convertZp(interpolation.getInterpolatingPolynomial(), variable).primitivePart(), a.domain, variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        return interpolated.multiply(asMultivariate(contentGCD, a.nVariables, variable, a.ordering));
    }

    static MultivariatePolynomial<BigInteger> sparseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {
        //todo choose first variable correctly

        BigInteger domainSize = a.domain.size();
        if (domainSize.isInt() && Math.min(a.degree(), b.degree()) > domainSize.intValue())
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");

        int evaluationStackLimit = domainSize.isInt() ? domainSize.intValue() : -1;
        int[] aDegrees = a.degrees(),
                bDegrees = b.degrees(),
                degreeBounds = new int[a.nVariables];
        for (int i = 0; i < a.nVariables; i++)
            degreeBounds[i] = Math.min(aDegrees[i], bDegrees[i]);

        MultivariatePolynomial<BigInteger> result = sparseModularGCD(a, b, PrivateRandom.getRandom(), a.nVariables - 1, degreeBounds, evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return result;
    }


    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> sparseModularGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //for bivariate polynomials always use dense interpolation
        if (variable == 1)
            return denseModularGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit);

        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("input is zero");

        MultivariatePolynomial<BigInteger> factory = a;
        if (a.isConstant() || b.isConstant())
            return factory.createOne();

        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            bMutablePolynomialZp gcd = PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b));
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        ModularDomain domain = (ModularDomain) factory.domain;

        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        MultivariatePolynomial<bMutablePolynomialZp>
                aConv = convertZp(a, variable),
                bConv = convertZp(b, variable);

        // content of a and b in Zp[x_k]
        bMutablePolynomialZp
                aCont = PolynomialGCD(aConv.data.values()),
                bCont = PolynomialGCD(bConv.data.values());

        // normalize a and b
        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(a, asMultivariate(aCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        a = qd[0];

        qd = divideAndRemainder(b, asMultivariate(bCont, nVariables, variable, a.ordering));
        assert qd[1].isZero();
        b = qd[0];

        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = PolynomialGCD(aCont, bCont),
                lcGCD = PolynomialGCD(aConv.lc(), bConv.lc());

//        List<bMutablePolynomialZp> evalTest = new ArrayList<>();
//        for (int i = 0; i < variable; ++i)
//            evalTest.add(PolynomialGCD(convertZp(a, i).lc(), convertZp(b, i).lc()));

        int prevVarExponent = degreeBounds[variable - 1];
        Set<BigInteger> evaluationStack = new HashSet<>();
        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                return null;

            BigInteger seedPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(seedPoint))
                continue;

            evaluationStack.add(seedPoint);

            BigInteger lcVal = lcGCD.evaluate(seedPoint);
            if (lcVal.isZero())
                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<BigInteger> cVal = sparseModularGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            SparseGCDInterpolation sparseInterpolator = new SparseGCDInterpolation(domain, variable, a, b, cVal, rnd);

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            // we are applying dense interpolation for univariate skeleton coefficients
            Interpolation<BigInteger> interpolation = new Interpolation<>(variable, seedPoint, cVal);
            while (true) {
//                if (interpolation.numberOfPoints() == degreeBounds[variable])
//                    return null;
                BigInteger randomPoint = domain.randomElement(rnd);
                if (evaluationStack.contains(randomPoint))
                    continue;

                evaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc().equals(lcVal);
                interpolation.update(randomPoint, cVal);


                // do division test

                MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }


//            if (currExponent < prevVarExponent) {
//                //better degree bound detected => start over
//                interpolation = new Interpolation<>(variable, seedPoint, cVal);
//                degreeBounds[variable - 1] = prevVarExponent = currExponent;
//                continue;
//            }
//
//            if (interpolation == null) {
//                //first successful homomorphism
//                interpolation = new Interpolation<>(variable, seedPoint, cVal);
//                continue;
//            }
//
//            previousInterpolation = interpolation.getInterpolatingPolynomial();
//            interpolation.update(seedPoint, cVal);

        }
    }

    public static final class SparseGCDInterpolation {
        final Domain<BigInteger> domain;
        final int variable;
        final MultivariatePolynomial<BigInteger> a, b, initialGCD;
        final RandomGenerator rnd;
        final Set<DegreeVector> globalSkeleton;
        final TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton;

        final int[] evaluationVariables;
        final BigInteger[] sparseEvaluationPoint;
        final PrecomputedPowersHolder<BigInteger> powers;
        final int[] sparseUnivarDegrees;

        public SparseGCDInterpolation(
                Domain<BigInteger> domain,
                int variable,
                MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                MultivariatePolynomial<BigInteger> skeleton, RandomGenerator rnd) {
            this.domain = domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.initialGCD = skeleton;
            this.rnd = rnd;

            globalSkeleton = skeleton.data.keySet();
            univarSkeleton = getSkeleton(skeleton);
            sparseUnivarDegrees = univarSkeleton.keys();

            evaluationVariables = ArraysUtil.sequence(1, variable + 1); //variable inclusive
            sparseEvaluationPoint = new BigInteger[evaluationVariables.length];

            for (int i = variable - 2; i >= 0; --i)
                do {
                    sparseEvaluationPoint[i] = domain.randomElement(rnd);
                } while (domain.isZero(sparseEvaluationPoint[i]));

            powers = new PrecomputedPowersHolder<>(sparseEvaluationPoint, domain);
        }

        public MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint) {
            // variable = newPoint
            sparseEvaluationPoint[sparseEvaluationPoint.length - 1] = newPoint;
            powers.set(sparseEvaluationPoint.length - 1, newPoint);

            //build Vandermonde matrices

            TIntObjectHashMap<VandermondeBuilder> lhsMatrices = new TIntObjectHashMap<>(univarSkeleton.size());
            univarSkeleton.forEachEntry((i, p) -> {
                lhsMatrices.put(i,
                        new VandermondeBuilder(p, powers, variable - 1)
                                .buildMatrix());
                return true;
            });


            int requiredNumberOfEvaluations = globalSkeleton.size();
            TIntObjectHashMap<BigInteger[]> rhsValues = new TIntObjectHashMap<>(requiredNumberOfEvaluations);
            for (int uDeg : sparseUnivarDegrees)
                rhsValues.put(uDeg, new BigInteger[requiredNumberOfEvaluations]);

            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                int[] raiseFactors = ArraysUtil.arrayOf(i + 1, evaluationVariables.length);
                raiseFactors[raiseFactors.length - 1] = 1;
                MultivariatePolynomial<BigInteger>
                        aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors),
                        bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors);

                bMutablePolynomialZp gcdUnivar = PolynomialGCD(asUnivariateZp(aUnivar), asUnivariateZp(bUnivar)).monic();

                for (int uDeg : sparseUnivarDegrees) {
                    BigInteger[] rhs = rhsValues.get(uDeg);
                    rhs[i] = gcdUnivar.degree() < uDeg ? BigInteger.ZERO : gcdUnivar.get(uDeg);
                }
            }

            // all Vandermonde systems are ready => solve and transform to mulivar poly

            MultivariatePolynomial<BigInteger> gcdVal = initialGCD.createZero();
            for (int uDeg : sparseUnivarDegrees) {
                VandermondeBuilder lhs = lhsMatrices.get(uDeg);
                BigInteger[] rhs = rhsValues.get(uDeg);
                BigInteger[] solution = gaussianElimination(domain, lhs.matrix, Arrays.copyOf(rhs, lhs.matrix.length));

                for (int i = 0; i < lhs.skeleton.length; i++) {
                    DegreeVector degreeVector = lhs.skeleton[i].set(0, uDeg);
                    BigInteger value = solution[i];
                    gcdVal.add(degreeVector, value);
                }
            }

            return gcdVal;
        }
    }

    private static TIntObjectHashMap<MultivariatePolynomial<BigInteger>> getSkeleton(MultivariatePolynomial<BigInteger> poly) {
        TIntObjectHashMap<MultivariatePolynomial<BigInteger>> skeleton = new TIntObjectHashMap<>();
        for (Entry<DegreeVector, BigInteger> term : poly.data.entrySet()) {
            DegreeVector dv = term.getKey(), newDV = dv.setZero(0);
            MultivariatePolynomial<BigInteger> coeff = skeleton.get(dv.exponents[0]);
            if (coeff != null)
                coeff.add(newDV, term.getValue());
            else
                skeleton.put(dv.exponents[0], MultivariatePolynomial.create(
                        poly.domain, poly.ordering, new DegreeVector[]{newDV}, new BigInteger[]{term.getValue()}));
        }
        return skeleton;
    }

    private static final class VandermondeBuilder {
        private final Domain<BigInteger> domain;
        private final DegreeVector[] skeleton;
        private final BigInteger[][] matrix;
        private final PrecomputedPowersHolder<BigInteger> powers;
        private final int dimension;

        public VandermondeBuilder(MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int dimension) {
            this.domain = skeleton.domain;
            this.skeleton = skeleton.data.keySet().toArray(new DegreeVector[skeleton.size()]);
            this.powers = powers;
            this.dimension = dimension;
            this.matrix = new BigInteger[skeleton.size()][];
        }

        @SuppressWarnings("unchecked")
        public VandermondeBuilder buildMatrix() {
            for (int i = 0; i < matrix.length; i++) {
                BigInteger[] row = new BigInteger[skeleton.length];
                for (int j = 0; j < skeleton.length; j++) {
                    BigInteger tmp = BigInteger.ONE;
                    for (int k = 0; k < dimension; k++)
                        tmp = domain.multiply(tmp, powers.pow(k, (i + 1) * skeleton[j].exponents[k + 1]));
                    row[j] = tmp;
                }
                matrix[i] = row;
            }
            return this;
        }
    }

    public static BigInteger[] gaussianElimination(Domain<BigInteger> domain, BigInteger[][] A, BigInteger[] b) {
        assert A.length == b.length;
        int N = b.length;

        for (int p = 0; p < N; p++) {

            // find pivot row and swap
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (domain.compare(A[i][p], A[max][p]) > 0) {
                    max = i;
                }
            }
            BigInteger[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            BigInteger t = b[p]; b[p] = b[max]; b[max] = t;

            // singular or nearly singular
            if (domain.isZero(A[p][p])) {
                throw new RuntimeException("Matrix is singular or nearly singular: " + Arrays.deepToString(A) + " = " + Arrays.toString(b));
            }

            // pivot within A and b
            for (int i = p + 1; i < N; i++) {
                BigInteger alpha = domain.divideAndRemainder(A[i][p], A[p][p])[0];
                b[i] = domain.subtract(b[i], domain.multiply(alpha, b[p]));
                for (int j = p; j < N; j++) {
                    A[i][j] = domain.subtract(A[i][j], domain.multiply(alpha, A[p][j]));
                }
            }
        }

        // back substitution
        BigInteger[] x = new BigInteger[N];
        for (int i = N - 1; i >= 0; i--) {
            BigInteger sum = domain.getZero();
            for (int j = i + 1; j < N; j++) {
                sum = domain.add(sum, domain.multiply(A[i][j], x[j]));
            }
            x[i] = domain.divideAndRemainder(domain.subtract(b[i], sum), A[i][i])[0];
        }
        return x;

    }

}
