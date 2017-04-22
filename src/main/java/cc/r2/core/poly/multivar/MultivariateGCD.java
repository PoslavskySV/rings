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
import java.util.stream.IntStream;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.univar.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.util.ArraysUtil.negate;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateGCD {
    private MultivariateGCD() {}

    /**
     * Evaluates {@code poly} at random point chosen in such way that the skeleton of evaluated version is the same as of the
     * original {@code poly} with respect to all except {@code variable} variables
     */
    private static MultivariatePolynomial<BigInteger> evaluateAtRandomPreservingSkeleton(
            MultivariatePolynomial<BigInteger> poly, int variable, RandomGenerator rnd) {
        //desired skeleton
        Set<DegreeVector> skeleton = poly.getSkeletonExcept(variable);
        MultivariatePolynomial<BigInteger> tmp;
        do {
            BigInteger randomPoint = poly.domain.randomElement(rnd);
            tmp = poly.evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    /** calculates the inverse permutation */
    public static int[] inversePermutation(int[] permutation) {
        final int[] inv = new int[permutation.length];
        for (int i = permutation.length - 1; i >= 0; --i)
            inv[permutation[i]] = i;
        return inv;
    }

    /** structure with required input for GCD algorithms */
    private static final class GCDInput {
        /** input polynomials (with variables renamed) and earlyGCD if possible (in trivial cases) */
        final MultivariatePolynomial<BigInteger> a, b, earlyGCD;
        /** gcd degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** domain cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int evaluationStackLimit;

        GCDInput(MultivariatePolynomial<BigInteger> earlyGCD) {
            this.earlyGCD = earlyGCD;
            a = b = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
        }

        GCDInput(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                 int evaluationStackLimit, int[] degreeBounds, int[] mapping, int lastPresentVariable) {
            this.a = a;
            this.b = b;
            this.earlyGCD = null;
            this.evaluationStackLimit = evaluationStackLimit;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
        }

        /** recover initial order of variables in the result */
        MultivariatePolynomial<BigInteger> restoreVariablesOrder(MultivariatePolynomial<BigInteger> result) {
            return renameVariables(result, mapping);
        }
    }

    /** prepare input for modular GCD algorithms (Brown, Zippel, LinZip) */
    private static GCDInput preparedGCDInput(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("Input is zero.");
        if (a.isConstant() || b.isConstant())
            return new GCDInput(a.createOne());

        Domain<BigInteger> domain = a.domain;
        BigInteger domainSize = domain.size();
        if (domainSize == null)
            throw new IllegalArgumentException("Modular gcd algorithms are supported only for multivariate polynomial over finite fields.");

        if (domainSize.isInt() && Math.min(a.degree(), b.degree()) > domainSize.intValue())
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        // domain cardinality, i.e. number of possible random choices
        int evaluationStackLimit = domainSize.isInt() ? domainSize.intValue() : -1;

        int
                nVariables = a.nVariables,
                nUsedVariables = 0, // number of really present variables in gcd
                lastPresentVariable = -1, // last variable that present in both input polynomials
                aDegrees[] = a.degrees(),
                bDegrees[] = b.degrees(),
                degreeBounds[] = new int[nVariables]; // degree bounds for gcd

        // determine initial degree bounds for gcd
        for (int i = 0; i < nVariables; i++) {
            degreeBounds[i] = Math.min(aDegrees[i], bDegrees[i]);
            if (degreeBounds[i] != 0) {
                ++nUsedVariables;
                lastPresentVariable = i;
            }
        }

        if (nUsedVariables == 0)
            // gcd is constant => 1
            return new GCDInput(a.createOne());

        RandomGenerator rnd = PrivateRandom.getRandom();
        if (nUsedVariables != nVariables) {
            // some of the variables are redundant in one of the polys => can just substitute a random values for them
            for (int i = 0; i < nVariables; i++) {
                if (degreeBounds[i] == 0) {
                    if (a.degree(i) != 0) {
                        assert b.degree(i) == 0;
                        a = evaluateAtRandomPreservingSkeleton(a, i, rnd);
                    } else if (b.degree(i) != 0) {
                        assert a.degree(i) == 0;
                        b = evaluateAtRandomPreservingSkeleton(b, i, rnd);
                    }
                }
            }
        }

        if (nUsedVariables == 1)
            // switch to univariate gcd
            return new GCDInput(asMultivariate(PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b)), nVariables, lastPresentVariable, a.ordering));

        // now swap variables so that the first variable will have the maximal degree (univariate gcd is fast),
        // and all non-used variables are at the end of poly's

        int[] variables = ArraysUtil.sequence(nVariables);
        //sort in descending order
        ArraysUtil.quickSort(negate(degreeBounds), variables);
        negate(degreeBounds);//recover degreeBounds

        lastPresentVariable = 0; //recalculate lastPresentVariable
        for (; lastPresentVariable < degreeBounds.length; ++lastPresentVariable)
            if (degreeBounds[lastPresentVariable] == 0)
                break;
        --lastPresentVariable;

        a = renameVariables(a, variables);
        b = renameVariables(b, variables);

        return new GCDInput(a, b, evaluationStackLimit, degreeBounds, variables, lastPresentVariable);
    }

    /** holds primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
    private static final class ReduceContent {
        /** primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
        final MultivariatePolynomial<BigInteger> aPrimitive, bPrimitive;
        /** gcd of content and leading coefficient of a and b given as Zp[x_k][x_1 ... x_{k-1}] */
        final bMutablePolynomialZp contentGCD, lcGCD;

        public ReduceContent(MultivariatePolynomial<BigInteger> aPrimitive, MultivariatePolynomial<BigInteger> bPrimitive,
                             bMutablePolynomialZp contentGCD, bMutablePolynomialZp lcGCD) {
            this.aPrimitive = aPrimitive;
            this.bPrimitive = bPrimitive;
            this.contentGCD = contentGCD;
            this.lcGCD = lcGCD;
        }
    }

    /**
     * Primitive parts, contents and leading coefficients of {@code a} and {@code b} viewed as polynomials in
     * Zp[x_k][x_1 ... x_{k-1}] with {@code x_k = variable}
     */
    @SuppressWarnings("unchecked")
    private static ReduceContent reduce(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, int variable) {
        MultivariatePolynomial<BigInteger> factory = a;
        int nVariables = factory.nVariables;

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
        return new ReduceContent(a, b, contentGCD, lcGCD);
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Brown's algorithm with dense interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static MultivariatePolynomial<BigInteger> BrownGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {

        // prepare input and test for early termination
        GCDInput gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        MultivariatePolynomial<BigInteger> result = BrownGCD(
                gcdInput.a, gcdInput.b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");

        return gcdInput.restoreVariablesOrder(result);
    }

    /**
     * Actual implementation of dense interpolation
     *
     * @param variable             current variable (all variables {@code v > variable} are fixed so far)
     * @param degreeBounds         degree bounds for gcd
     * @param evaluationStackLimit domain cardinality
     */
    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> BrownGCD(
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

        ReduceContent reduceContent = reduce(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = reduceContent.aPrimitive;
        b = reduceContent.bPrimitive;
        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = reduceContent.contentGCD,
                lcGCD = reduceContent.lcGCD;

        Domain<BigInteger> domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //dense interpolation
        Interpolation<BigInteger> interpolation = null;
        //previous interpolation (used to detect whether update doesn't change the result)
        MultivariatePolynomial<BigInteger> previousInterpolation;
        //store points that were already used in interpolation
        Set<BigInteger> evaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // all elements of the domain are tried
                // do division check (last chance) and return
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            //pickup the next random element for variable
            BigInteger randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;
            evaluationStack.add(randomPoint);

            BigInteger lcVal = lcGCD.evaluate(randomPoint);
            if (lcVal.isZero())
                continue;

            // evaluate a and b at variable = randomPoint
            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<BigInteger> cVal = BrownGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
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

            // interpolate
            previousInterpolation = interpolation.getInterpolatingPolynomial();
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] == 0 || degreeBounds[variable] == interpolation.numberOfPoints()
                    || (previousInterpolation != null && previousInterpolation.equals(interpolation.getInterpolatingPolynomial()))) {
                MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    /** division test **/
    private static MultivariatePolynomial<BigInteger> doDivisionCheck(
            MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
            bMutablePolynomialZp contentGCD, Interpolation<BigInteger> interpolation, int variable) {
        if (interpolation == null)
            return null;
        MultivariatePolynomial<BigInteger> interpolated =
                fromZp(convertZp(interpolation.getInterpolatingPolynomial(), variable).primitivePart(), a.domain, variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        if (contentGCD == null)
            return interpolated;
        return interpolated.multiply(asMultivariate(contentGCD, a.nVariables, variable, a.ordering));
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Zippel's algorithm with sparse interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unhecked")
    public static MultivariatePolynomial<BigInteger> ZippelGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {

        // prepare input and test for early termination
        GCDInput gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        MultivariatePolynomial<BigInteger> aLC = gcdInput.a.coefficientOf(0, gcdInput.a.degree(0));
        MultivariatePolynomial<BigInteger> bLC = gcdInput.b.coefficientOf(0, gcdInput.b.degree(0));
        MultivariatePolynomial<BigInteger> lcFactor = ZippelGCD(aLC, bLC);
        if (lcFactor.isConstant()) {
            // monic case
            MultivariatePolynomial<BigInteger> result = monicZippelGCD(gcdInput.a, gcdInput.b, PrivateRandom.getRandom(),
                    gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
            if (result == null)
                throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
            return gcdInput.restoreVariablesOrder(result);
        }

        MultivariatePolynomial<BigInteger> result = nonMonicZippelGCD(lcFactor,
                gcdInput.a.clone().multiply(lcFactor), gcdInput.b.clone().multiply(lcFactor), PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        result = gcdInput.restoreVariablesOrder(result);

        MultivariatePolynomial<BigInteger> r = result;
        MultivariatePolynomial<BigInteger> content = IntStream.rangeClosed(0, result.degree(0))
                .mapToObj(i -> r.coefficientOf(0, i)).reduce(null, (u, v) -> {
                    if (u == null)
                        return v;
                    if (v == null)
                        return u;
                    return ZippelGCD(u, v);
                });


        return MultivariateReduction.divideAndRemainder(result, content)[0];
    }

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> monicZippelGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //for bivariate polynomials always use dense interpolation
        if (variable == 1)
            return BrownGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit);

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

        ReduceContent reduceContent = reduce(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = reduceContent.aPrimitive;
        b = reduceContent.bPrimitive;
        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = reduceContent.contentGCD,
                lcGCD = reduceContent.lcGCD;

        Domain<BigInteger> domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //store points that were already used in interpolation
        Set<BigInteger> globalEvaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == globalEvaluationStack.size())
                return null;

            BigInteger seedPoint = domain.randomElement(rnd);
            if (globalEvaluationStack.contains(seedPoint))
                continue;

            globalEvaluationStack.add(seedPoint);

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

            MultivariatePolynomial<BigInteger> cVal = monicZippelGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            if (currExponent < prevVarExponent) {
                //better degree bound detected
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
            }

            SparseGCDInterpolation sparseInterpolator = new SparseGCDInterpolation(domain, variable, a, b, cVal, rnd);

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            // we are applying dense interpolation for univariate skeleton coefficients
            Interpolation<BigInteger> interpolation = new Interpolation<>(variable, seedPoint, cVal);
            //previous interpolation (used to detect whether update doesn't change the result)
            MultivariatePolynomial<BigInteger> previousInterpolation;
            //local evaluation stack for points that are calculated via sparse interpolation (but not gcd evaluation) -> always same skeleton
            HashSet<BigInteger> localEvaluationStack = new HashSet<>(globalEvaluationStack);
            while (true) {
                if (interpolation.numberOfPoints() > degreeBounds[variable])
                    continue main;
                BigInteger randomPoint = domain.randomElement(rnd);
                if (localEvaluationStack.contains(randomPoint))
                    continue;

                localEvaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc().equals(lcVal);

                previousInterpolation = interpolation.getInterpolatingPolynomial();
                interpolation.update(randomPoint, cVal);

                if (degreeBounds[variable] == interpolation.numberOfPoints()
                        || previousInterpolation.equals(interpolation.getInterpolatingPolynomial())) {
                    // do division test
                    MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                    if (result != null)
                        return result;
                }
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> nonMonicZippelGCD(
            MultivariatePolynomial<BigInteger> lcGCD0,
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {
        //for bivariate polynomials always use dense interpolation
        if (variable == 1)
            return BrownGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit);

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

//        // a and b as Zp[x_k][x_1 ... x_{k-1}]
//        MultivariatePolynomial<bMutablePolynomialZp>
//                aConv = convertZp(a, variable),
//                bConv = convertZp(b, variable);
//
//        // content of a and b in Zp[x_k]
//        bMutablePolynomialZp
//                aCont = PolynomialGCD(aConv.data.values()),
//                bCont = PolynomialGCD(bConv.data.values());
//
//        // normalize a and b
//        MultivariatePolynomial<BigInteger>[] qd;
//        qd = divideAndRemainder(a, asMultivariate(aCont, nVariables, variable, a.ordering));
//        assert qd[1].isZero();
//        a = qd[0];
//
//        qd = divideAndRemainder(b, asMultivariate(bCont, nVariables, variable, a.ordering));
//        assert qd[1].isZero();
//        b = qd[0];
//
//        // gcd of Zp[x_k] content and lc
//        bMutablePolynomialZp
//                contentGCD = PolynomialGCD(aCont, bCont),
//                lcGCD = PolynomialGCD(aConv.lc(), bConv.lc());


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

//            BigInteger lcVal = lcGCD.evaluate(seedPoint);
//            if (lcVal.isZero())
//                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<BigInteger>
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint),
                    lcGCD0Val = lcGCD0.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<BigInteger> cVal = nonMonicZippelGCD(lcGCD0Val, aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            SparseGCDInterpolation sparseInterpolator = new SparseGCDInterpolation(domain, variable, a, b, cVal, rnd);

//            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
//            assert cVal.lc().equals(lcVal);

            // we are applying dense interpolation for univariate skeleton coefficients
            MultivariatePolynomial<BigInteger> previousInterpolation = null;
            Interpolation<BigInteger> interpolation = new Interpolation<>(variable, seedPoint, cVal);
            while (true) {
                BigInteger randomPoint = domain.randomElement(rnd);
                if (evaluationStack.contains(randomPoint))
                    continue;

                evaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
//                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
//                assert cVal.lc().equals(lcVal);
                interpolation.update(randomPoint, cVal);


                // do division test
                if (previousInterpolation == interpolation.getInterpolatingPolynomial()) {
                    MultivariatePolynomial<BigInteger> result = doDivisionCheck(a, b, null, interpolation, variable);
                    if (result != null)
                        return result;
                }
                previousInterpolation = interpolation.getInterpolatingPolynomial();
            }
        }
    }

    public static final class SparseGCDInterpolation {
        final Domain<BigInteger> domain;
        final int variable;
        final MultivariatePolynomial<BigInteger> lcGCD0, a, b, initialGCD;
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
            this(domain, variable, skeleton.createOne(), a, b, skeleton, rnd);
        }

        public SparseGCDInterpolation(
                Domain<BigInteger> domain,
                int variable,
                MultivariatePolynomial<BigInteger> lcGCD0,
                MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                MultivariatePolynomial<BigInteger> skeleton, RandomGenerator rnd) {
            this.domain = domain;
            this.variable = variable;
            this.lcGCD0 = lcGCD0;
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
                MultivariatePolynomial<BigInteger> normalization = lcGCD0.evaluate(powers, evaluationVariables, raiseFactors);
                assert normalization.degree() == 0;
                gcdUnivar = gcdUnivar.multiply(normalization.lc());

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
