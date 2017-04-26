package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.multivar.LinearAlgebra.SystemInfo;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

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
        /** GCD of monomial content of a and b **/
        final DegreeVector monomialGCD;

        GCDInput(MultivariatePolynomial<BigInteger> earlyGCD) {
            this.earlyGCD = earlyGCD;
            a = b = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialGCD = null;
        }

        GCDInput(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, DegreeVector monomialGCD,
                 int evaluationStackLimit, int[] degreeBounds, int[] mapping, int lastPresentVariable) {
            this.a = a;
            this.b = b;
            this.monomialGCD = monomialGCD;
            this.earlyGCD = null;
            this.evaluationStackLimit = evaluationStackLimit;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
        }

        /** recover initial order of variables in the result */
        MultivariatePolynomial<BigInteger> restoreGCD(MultivariatePolynomial<BigInteger> result) {
            return renameVariables(result.multiply(monomialGCD, a.domain.getOne()), mapping);
        }
    }

    private static MultivariatePolynomial<BigInteger> trivialGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("Input is zero.");
        if (a.isConstant() || b.isConstant())
            return a.createOne();
        if (a.size() == 1)
            return gcdWithMonomial(a.lt().getKey(), b);
        if (b.size() == 1)
            return gcdWithMonomial(b.lt().getKey(), a);
        return null;
    }

    /** prepare input for modular GCD algorithms (Brown, Zippel, LinZip) */
    private static GCDInput preparedGCDInput(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        MultivariatePolynomial<BigInteger> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return new GCDInput(trivialGCD);

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

        // find monomial GCD
        DegreeVector aMonomialContent = a.monomialContent();
        int[] exponentsGCD = b.monomialContent().exponents;
        MultivariatePolynomial.setMin(aMonomialContent, exponentsGCD);
        DegreeVector monomialGCD = new DegreeVector(exponentsGCD);

        a = a.divideOrNull(monomialGCD, a.domain.getOne());
        b = b.divideOrNull(monomialGCD, b.domain.getOne());

        assert a != null && b != null;

        return new GCDInput(a, b, monomialGCD, evaluationStackLimit, degreeBounds, variables, lastPresentVariable);
    }

    /** gcd with monomial */
    private static MultivariatePolynomial<BigInteger> gcdWithMonomial(DegreeVector monomial,
                                                                      MultivariatePolynomial<BigInteger> poly) {
        return poly.create(poly.commonContent(monomial), BigInteger.ONE);
    }

    /**
     * Primitive part and content of multivariate polynomial considered as polynomial over Zp[x_i][x_1, ..., x_N]
     */
    private static final class UnivariateContent {
        final MultivariatePolynomial<bMutablePolynomialZp> poly;
        final MultivariatePolynomial<BigInteger> primitivePart;
        final bMutablePolynomialZp content;

        public UnivariateContent(MultivariatePolynomial<bMutablePolynomialZp> poly, MultivariatePolynomial<BigInteger> primitivePart, bMutablePolynomialZp content) {
            this.poly = poly;
            this.primitivePart = primitivePart;
            this.content = content;
        }
    }

    /**
     * Returns primitive part and content of {@code poly} considered as polynomial over Zp[variable][x_1, ..., x_N]
     *
     * @param poly     the polynomial
     * @param variable the variable
     * @return primitive part and content
     */
    @SuppressWarnings("unchecked")
    private static UnivariateContent univariateContent(
            MultivariatePolynomial<BigInteger> poly, int variable) {
        //convert poly to Zp[var][x...]
        MultivariatePolynomial<bMutablePolynomialZp> conv = convertZp(poly, variable);
        //univariate content
        bMutablePolynomialZp content = PolynomialGCD(conv.data.values());
        MultivariatePolynomial<BigInteger>[]
                qd = divideAndRemainder(poly, asMultivariate(content, poly.nVariables, variable, poly.ordering));
        assert qd[1].isZero();
        return new UnivariateContent(conv, qd[0], content);
    }

    /** holds primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
    private static final class PrimitiveInput {
        /** primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
        final MultivariatePolynomial<BigInteger> aPrimitive, bPrimitive;
        /** gcd of content and leading coefficient of a and b given as Zp[x_k][x_1 ... x_{k-1}] */
        final bMutablePolynomialZp contentGCD, lcGCD;

        public PrimitiveInput(MultivariatePolynomial<BigInteger> aPrimitive, MultivariatePolynomial<BigInteger> bPrimitive,
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
    private static PrimitiveInput makePrimitive(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, int variable) {
        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        UnivariateContent
                aContent = univariateContent(a, variable),
                bContent = univariateContent(b, variable);

        a = aContent.primitivePart;
        b = bContent.primitivePart;

        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = PolynomialGCD(aContent.content, bContent.content),
                lcGCD = PolynomialGCD(aContent.poly.lc(), bContent.poly.lc());
        return new PrimitiveInput(a, b, contentGCD, lcGCD);
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

        return gcdInput.restoreGCD(result);
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

        PrimitiveInput primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

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
    @SuppressWarnings("unchecked")
    public static MultivariatePolynomial<BigInteger> ZippelGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {

        // prepare input and test for early termination
        GCDInput gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.a;
        b = gcdInput.b;

        MultivariatePolynomial<BigInteger> result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return gcdInput.restoreGCD(result);
    }


    private static final int MAX_SPARSE_INTERPOLATION_FAILS = 1000;

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger>[] multivariateCoefficients(
            MultivariatePolynomial<BigInteger> a,
            int variable) {
        MultivariatePolynomial[] mCfs = Arrays.stream(a.degrees(variable)).mapToObj(d -> a.coefficientOf(variable, d)).toArray(MultivariatePolynomial[]::new);
        int minSizeCoefficient = 0;
        int minSize = Integer.MAX_VALUE;
        for (int i = 0; i < mCfs.length; i++)
            if (mCfs[i].size() < minSize) {
                minSize = mCfs[i].size();
                minSizeCoefficient = i;
            }

        ArraysUtil.swap(mCfs, 0, minSizeCoefficient);
        return mCfs;
    }

    static MultivariatePolynomial<BigInteger> ZippelContentGCD(
            MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd, int variable, int[] degreeBounds, int evaluationStackLimit) {
        MultivariatePolynomial<BigInteger>[] aCfs = multivariateCoefficients(a, variable);
        MultivariatePolynomial<BigInteger>[] bCfs = multivariateCoefficients(b, variable);

        MultivariatePolynomial<BigInteger> contentGCD = ZippelGCD(aCfs[0], bCfs[0], rnd, variable - 1, degreeBounds.clone(), evaluationStackLimit);
        for (int i = 1; i < aCfs.length; i++)
            contentGCD = ZippelGCD(aCfs[i], contentGCD, rnd, variable - 1, degreeBounds.clone(), evaluationStackLimit);
        for (int i = 1; i < bCfs.length; i++)
            contentGCD = ZippelGCD(bCfs[i], contentGCD, rnd, variable - 1, degreeBounds.clone(), evaluationStackLimit);
        return contentGCD;
    }

    @SuppressWarnings("unchecked")
    private static MultivariatePolynomial<BigInteger> ZippelGCD(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        MultivariatePolynomial<BigInteger> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

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

        PrimitiveInput primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        bMutablePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD.multiply(asMultivariate(contentGCD, a.nVariables, variable, a.ordering));

//        System.out.println("===" + ZippelContentGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit));

        Domain<BigInteger> domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //store points that were already used in interpolation
        Set<BigInteger> globalEvaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        int failedSparseInterpolations = 0;
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

            MultivariatePolynomial<BigInteger> cVal = ZippelGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
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

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            LinZipGCDInterpolation sparseInterpolator = new LinZipGCDInterpolation(variable, a, b, cVal, rnd);

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

                lcVal = lcGCD.evaluate(randomPoint);
                if (lcVal.isZero())
                    continue;

                localEvaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
                if (cVal == null) {
                    ++failedSparseInterpolations;
                    if (failedSparseInterpolations == MAX_SPARSE_INTERPOLATION_FAILS)
                        throw new RuntimeException("Sparse interpolation failed");
                    continue main;
                }
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

    static final class LinZipGCDInterpolation {
        /** the domain */
        final Domain<BigInteger> domain;
        /** which we are evaluating */
        final int variable;
        /** initial polynomials and initial gcd (with variable substituted) */
        final MultivariatePolynomial<BigInteger> a, b;//, initialGCD;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final BigInteger[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final PrecomputedPowersHolder<BigInteger> powers;
        /** random */
        final RandomGenerator rnd;
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        final int monicScalingExponent;
        final boolean monic;


        LinZipGCDInterpolation(
                int variable, MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                MultivariatePolynomial<BigInteger> skeleton, RandomGenerator rnd) {
            this.monic =
                    a.coefficientOf(0, a.degree(0)).isConstant()
                            && b.coefficientOf(0, a.degree(0)).isConstant();
            this.domain = a.domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.rnd = rnd;

            this.globalSkeleton = skeleton.data.keySet();
            this.univarSkeleton = getSkeleton(skeleton);
            this.sparseUnivarDegrees = univarSkeleton.keys();

            this.evaluationVariables = ArraysUtil.sequence(1, variable + 1); //variable inclusive
            this.evaluationPoint = new BigInteger[evaluationVariables.length];

            //avoid zero evaluation points // TODO check the same skeleton in main variable!!!
            for (int i = variable - 2; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (domain.isZero(evaluationPoint[i]));

            this.powers = new PrecomputedPowersHolder<>(evaluationPoint, domain);

            int minimalSkeletonSize = -1, monicScalingExponent = -1;
            for (TIntObjectIterator<MultivariatePolynomial<BigInteger>> it = univarSkeleton.iterator(); it.hasNext(); ) {
                it.advance();
                MultivariatePolynomial<BigInteger> v = it.value();
                if (v.size() > minimalSkeletonSize)
                    minimalSkeletonSize = v.size();
                if (v.size() == 1)
                    monicScalingExponent = it.key();
            }
            this.requiredNumberOfEvaluations = minimalSkeletonSize;
            this.monicScalingExponent = monicScalingExponent;
        }

        /** Number of retries when raise condition occurs; then drop up with new homomorphism */
        private static final int NUMBER_OF_UNDER_DETERMINED_RETRIES = 100;

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        public MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            ZippelSystem[] systems = new ZippelSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new ZippelSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            // the last variable (the variable) is the same for all evaluations = newPoint
            raiseFactors[raiseFactors.length - 1] = 1;

            int nUnknowns = globalSkeleton.size(), nUnknownScalings = 0;
            int raiseFactor = 0;

            for (int nTries = 0; nTries < NUMBER_OF_UNDER_DETERMINED_RETRIES; ++nTries) {
                int previousFreeVars = -1, underDeterminedTries = 0;
                for (; ; ) {
                    // increment at each loop!
                    ++raiseFactor;
                    // sequential powers of evaluation point
                    Arrays.fill(raiseFactors, 0, raiseFactors.length - 1, raiseFactor);
                    // evaluate a and b to univariate and calculate gcd
                    bMutablePolynomialZp
                            aUnivar = asUnivariateZp(a.evaluate(powers, evaluationVariables, raiseFactors)),
                            bUnivar = asUnivariateZp(b.evaluate(powers, evaluationVariables, raiseFactors)),
                            gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                    assert gcdUnivar.isMonic();

                    int totalEquations = 0;
                    for (ZippelSystem system : systems) {
                        BigInteger rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs, nUnknownScalings - 1);
                        totalEquations += system.nEquations();
                    }
                    if (nUnknowns + nUnknownScalings <= totalEquations)
                        break;

                    if (underDeterminedTries > NUMBER_OF_UNDER_DETERMINED_RETRIES) {
                        System.out.println("underDeterminedTries");
                        // raise condition: new equations does not fix enough variables
                        return null;
                    }

                    int freeVars = nUnknowns + nUnknownScalings - totalEquations;
                    if (freeVars >= previousFreeVars)
                        ++underDeterminedTries;
                    else
                        underDeterminedTries = 0;

                    previousFreeVars = freeVars;

                    ++nUnknownScalings;
                }

                MultivariatePolynomial<BigInteger> result = a.createZero();
                SystemInfo info = solveLinZip(a, systems, nUnknownScalings, result);
                if (info == SystemInfo.UnderDetermined) {
                    //try to generate more equations
                    continue;
                }

                if (info == SystemInfo.Consistent)
                    //well done
                    return result;
                if (info == SystemInfo.Inconsistent)
                    //inconsistent system => unlucky homomorphism
                    return null;
            }

            // the system is still under-determined => bad evaluation homomorphism
            System.out.println("UnderDetermined");
            return null;
        }
    }


    /** Vandermonde system builder */
    private static class ZippelSystem {
        private final int univarDegree;
        /** the domain */
        private final Domain<BigInteger> domain;
        /** the skeleton */
        private final DegreeVector[] skeleton;
        /** the lhs matrix */
        private final ArrayList<BigInteger[]> matrix;
        /** the rhs values */
        private final ArrayList<BigInteger> rhs = new ArrayList<>();
        /** precomputed powers */
        private final PrecomputedPowersHolder<BigInteger> powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        private final int nVars;

        public ZippelSystem(int univarDegree, MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            this.skeleton = skeleton.data.keySet().toArray(new DegreeVector[skeleton.size()]);
            this.powers = powers;
            this.nVars = nVars;
            this.matrix = new ArrayList<>();
        }

        int nUnknownVariables() {
            return skeleton.length;
        }

        int nEquations() {
            return matrix.size();
        }

        private final ArrayList<BigInteger[]> scalingMatrix = new ArrayList<>();

        public void oneMoreEquation(BigInteger rhsVal, int scalingVar) {
            BigInteger[] row = new BigInteger[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);

            //todo: replace with single elemebt
            BigInteger[] scalingCoeffs = new BigInteger[scalingVar + 1];
            Arrays.fill(scalingCoeffs, BigInteger.ZERO);
            if (scalingVar != -1) {
                scalingCoeffs[scalingVar] = domain.negate(rhsVal);
                rhsVal = domain.getZero();
            }
            scalingMatrix.add(scalingCoeffs);
            rhs.add(rhsVal);
        }

        @Override
        public String toString() {
            return "{" + matrix.stream().map(Arrays::toString).collect(Collectors.joining(",")) + "} = " + rhs;
        }
    }

    private static SystemInfo solveLinZip(MultivariatePolynomial<BigInteger> factory, ZippelSystem[] subSystems, int nUnknownScalings, MultivariatePolynomial<BigInteger> destination) {
        ArrayList<DegreeVector> unknowns = new ArrayList<>();
        for (ZippelSystem system : subSystems)
            for (DegreeVector degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = unknowns.size() + nUnknownScalings;
        ArrayList<BigInteger[]> lhsGlobal = new ArrayList<>();
        ArrayList<BigInteger> rhsGlobal = new ArrayList<>();
        int offset = 0;
        for (ZippelSystem system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                BigInteger[] row = new BigInteger[nUnknownsTotal];
                Arrays.fill(row, BigInteger.ZERO);
                BigInteger[] subRow = system.matrix.get(j);
                BigInteger[] scalingRow = system.scalingMatrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                System.arraycopy(scalingRow, 0, row, nUnknownsMonomials, scalingRow.length);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        BigInteger[] solution = new BigInteger[nUnknownsTotal];
        SystemInfo info = LinearAlgebra.solve(factory.domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent)
            destination.add(unknowns.toArray(new DegreeVector[unknowns.size()]), Arrays.copyOf(solution, nUnknownsMonomials));
        return info;
    }


    static final class MonicGCDInterpolation {
        /** the domain */
        final Domain<BigInteger> domain;
        /** which we are evaluating */
        final int variable;
        /** initial polynomials and initial gcd (with variable substituted) */
        final MultivariatePolynomial<BigInteger> a, b;//, initialGCD;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final BigInteger[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final PrecomputedPowersHolder<BigInteger> powers;
        /** random */
        final RandomGenerator rnd;
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        final int monicScalingExponent;
        final boolean monic;


        MonicGCDInterpolation(
                int variable, MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                MultivariatePolynomial<BigInteger> skeleton, RandomGenerator rnd) {
            this.monic =
                    a.coefficientOf(0, a.degree(0)).isConstant()
                            && b.coefficientOf(0, a.degree(0)).isConstant();
            this.domain = a.domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.rnd = rnd;

            this.globalSkeleton = skeleton.data.keySet();
            this.univarSkeleton = getSkeleton(skeleton);
            this.sparseUnivarDegrees = univarSkeleton.keys();

            this.evaluationVariables = ArraysUtil.sequence(1, variable + 1); //variable inclusive
            this.evaluationPoint = new BigInteger[evaluationVariables.length];

            //avoid zero evaluation points // TODO check the same skeleton in main variable!!!
            for (int i = variable - 2; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (domain.isZero(evaluationPoint[i]));

            this.powers = new PrecomputedPowersHolder<>(evaluationPoint, domain);

            int minimalSkeletonSize = -1, monicScalingExponent = -1;
            for (TIntObjectIterator<MultivariatePolynomial<BigInteger>> it = univarSkeleton.iterator(); it.hasNext(); ) {
                it.advance();
                MultivariatePolynomial<BigInteger> v = it.value();
                if (v.size() > minimalSkeletonSize)
                    minimalSkeletonSize = v.size();
                if (v.size() == 1)
                    monicScalingExponent = it.key();
            }
            this.requiredNumberOfEvaluations = minimalSkeletonSize;
            this.monicScalingExponent = monicScalingExponent;
        }

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        public MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            VandermondeSystem[] systems = new VandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new VandermondeSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                Arrays.fill(raiseFactors, 0, raiseFactors.length - 1, i + 1);
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;
                // evaluate a and b to univariate and calculate gcd
                bMutablePolynomialZp
                        aUnivar = asUnivariateZp(a.evaluate(powers, evaluationVariables, raiseFactors)),
                        bUnivar = asUnivariateZp(b.evaluate(powers, evaluationVariables, raiseFactors)),
                        gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                assert gcdUnivar.isMonic();

                if (monicScalingExponent != -1) {
                    // single scaling factor
                    // scale the system according to it

                    if (gcdUnivar.degree() < monicScalingExponent || domain.isZero(gcdUnivar.get(monicScalingExponent)))
                        // unlucky homomorphism
                        return null;

                    BigInteger normalization = evaluateExceptFirst(domain, powers, domain.getOne(), univarSkeleton.get(monicScalingExponent).lt().getKey(), i + 1, variable - 1);
                    //normalize univariate gcd in order to reconstruct leading coefficient polynomial
                    normalization = domain.multiply(domain.reciprocal(gcdUnivar.get(monicScalingExponent)), normalization);
                    gcdUnivar = gcdUnivar.multiply(normalization);
                }

                boolean allDone = true;
                for (VandermondeSystem system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        BigInteger rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (VandermondeSystem system : systems) {
                //solve each system
                boolean solved = system.solve();
                if (!solved)
                    // system is inconsistent
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomial<BigInteger> gcdVal = a.createZero();
            for (VandermondeSystem system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || domain.isOne(system.solution[0]);
                for (int i = 0; i < system.skeleton.length; i++) {
                    DegreeVector degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    BigInteger value = system.solution[i];
                    gcdVal.add(degreeVector, value);
                }
            }

            return gcdVal;
        }
    }


    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
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

    /** Vandermonde system builder */
    private static final class VandermondeSystem {
        private final int univarDegree;
        /** the domain */
        private final Domain<BigInteger> domain;
        /** the skeleton */
        private final DegreeVector[] skeleton;
        /** the lhs matrix */
        private final ArrayList<BigInteger[]> matrix;
        /** the rhs values */
        private final ArrayList<BigInteger> rhs = new ArrayList<>();
        /** precomputed powers */
        private final PrecomputedPowersHolder<BigInteger> powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        private final int nVars;

        public VandermondeSystem(int univarDegree, MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            this.skeleton = skeleton.data.keySet().toArray(new DegreeVector[skeleton.size()]);
            this.powers = powers;
            this.nVars = nVars;
            this.matrix = new ArrayList<>();
        }

        int nUnknownVariables() {
            return skeleton.length;
        }

        int nEquations() {
            return matrix.size();
        }

        BigInteger[] solution = null;

        boolean solve() {
            solution = LinearAlgebra.solve(domain, matrix.toArray(new BigInteger[matrix.size()][]), rhs.toArray(new BigInteger[rhs.size()]));
            return true;
        }

        public VandermondeSystem oneMoreEquation(BigInteger rhsVal) {
            BigInteger[] row = new BigInteger[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);
            rhs.add(rhsVal);
            return this;
        }
    }

    private static BigInteger evaluateExceptFirst(Domain<BigInteger> domain,
                                                  PrecomputedPowersHolder<BigInteger> powers,
                                                  BigInteger coefficient,
                                                  DegreeVector skeleton,
                                                  int raiseFactor,
                                                  int nVars) {
        BigInteger tmp = coefficient;
        for (int k = 0; k < nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[1 + k]));
        return tmp;
    }

    private static boolean isVandermonde(BigInteger[][] lhs, BigInteger[] rhs, Domain<BigInteger> domain) {
        for (int i = 1; i < lhs.length; i++) {
            if (!rhs[i].equals(domain.pow(rhs[0], i + 1)))
                return false;
            for (int j = 0; j < lhs[0].length; j++) {
                if (!lhs[i][j].equals(domain.pow(lhs[0][j], i + 1)))
                    return false;
            }
        }
        return true;
    }
}
