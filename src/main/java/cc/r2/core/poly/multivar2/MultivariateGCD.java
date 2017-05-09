package cc.r2.core.poly.multivar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar2.LinearAlgebra.SystemInfo;
import cc.r2.core.poly.multivar2.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.stream.Collectors;

import static cc.r2.core.poly.multivar2.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar2.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar2.MultivariateReduction.dividesQ;
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
        final MultivariatePolynomial<BigInteger> aReduced, bReduced, earlyGCD;
        /** gcd degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** domain cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int evaluationStackLimit;
        /** GCD of monomial content of a and b **/
        final MonomialTerm<BigInteger> monomialGCD;

        GCDInput(MultivariatePolynomial<BigInteger> earlyGCD) {
            this.earlyGCD = earlyGCD;
            aReduced = bReduced = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialGCD = null;
        }

        GCDInput(MultivariatePolynomial<BigInteger> aReduced, MultivariatePolynomial<BigInteger> bReduced, MonomialTerm<BigInteger> monomialGCD,
                 int evaluationStackLimit, int[] degreeBounds, int[] mapping, int lastPresentVariable) {
            assert monomialGCD == null || aReduced.domain.isOne(monomialGCD.coefficient);
            this.aReduced = aReduced;
            this.bReduced = bReduced;
            this.monomialGCD = monomialGCD;
            this.earlyGCD = null;
            this.evaluationStackLimit = evaluationStackLimit;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
        }

        /** recover initial order of variables in the result */
        MultivariatePolynomial<BigInteger> restoreGCD(MultivariatePolynomial<BigInteger> result) {
            return renameVariables(result, mapping).multiply(monomialGCD);
        }
    }

    private static MultivariatePolynomial<BigInteger> trivialGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        if (a.isZero() || b.isZero())
            throw new IllegalArgumentException("Input is zero.");
        if (a.isConstant() || b.isConstant())
            return a.createOne();
        if (a.size() == 1)
            return gcdWithMonomial(a.lt(), b);
        if (b.size() == 1)
            return gcdWithMonomial(b.lt(), a);
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

        // find monomial GCD
        // and remove monomial content from a and b
        a = a.clone(); b = b.clone(); // prevent rewriting original data
        MonomialTerm<BigInteger> monomialGCD = reduceMonomialContent(a, b);

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
            return new GCDInput(a.create(monomialGCD));

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
            return new GCDInput(asMultivariate(PolynomialGCD(asUnivariateZp(a), asUnivariateZp(b)), nVariables, lastPresentVariable, a.ordering).multiply(monomialGCD));

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

        return new GCDInput(a, b, monomialGCD, evaluationStackLimit, degreeBounds, variables, lastPresentVariable);
    }

    /** gcd with monomial */
    private static MultivariatePolynomial<BigInteger> gcdWithMonomial(MonomialTerm<BigInteger> monomial,
                                                                      MultivariatePolynomial<BigInteger> poly) {
        return poly.create(poly.commonContent(monomial));
    }


    /**
     * Removes monomial content from a and b and returns monomial gcd
     */
    private static MonomialTerm<BigInteger> reduceMonomialContent(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {

        MonomialTerm<BigInteger> aMonomialContent = a.monomialContent();
        int[] exponentsGCD = b.monomialContent().exponents;
        MultivariatePolynomial.setMin(aMonomialContent, exponentsGCD);
        MonomialTerm<BigInteger> monomialGCD = new MonomialTerm<>(exponentsGCD, a.domain.getOne());

        a = a.divideOrNull(monomialGCD);
        b = b.divideOrNull(monomialGCD);
        assert a != null && b != null;

        return new MonomialTerm<>(exponentsGCD, a.domain.getOne());
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
        bMutablePolynomialZp content = PolynomialGCD(conv.coefficients());
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
                gcdInput.aReduced, gcdInput.bReduced, PrivateRandom.getRandom(),
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

        //check for trivial gcd
        MultivariatePolynomial<BigInteger> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        MultivariatePolynomial<BigInteger> factory = a;
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

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

//        MultivariatePolynomial<BigInteger> content = ZippelContentGCD(a, b, 0);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        MultivariatePolynomial<BigInteger> result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return gcdInput.restoreGCD(result);
    }

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
            int variable) {
        MultivariatePolynomial<BigInteger>[] aCfs = multivariateCoefficients(a, variable);
        MultivariatePolynomial<BigInteger>[] bCfs = multivariateCoefficients(b, variable);

        MultivariatePolynomial<BigInteger> contentGCD = ZippelGCD(aCfs[0], bCfs[0]);
        for (int i = 1; i < aCfs.length; i++)
            contentGCD = ZippelGCD(aCfs[i], contentGCD);
        for (int i = 1; i < bCfs.length; i++)
            contentGCD = ZippelGCD(bCfs[i], contentGCD);
        return contentGCD;
    }


    /** Maximal number of fails before switch to a new homomorphism */
    private static final int MAX_SPARSE_INTERPOLATION_FAILS = 1000;

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

//        MultivariatePolynomial<BigInteger> content = ZippelContentGCD(a, b, variable);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        MultivariatePolynomial<BigInteger> factory = a;
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

            SparseInterpolation sparseInterpolator = createInterpolation(variable, a, b, cVal, rnd);

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

    static boolean ALWAYS_LINZIP = false;

    static SparseInterpolation createInterpolation(int variable,
                                                   MultivariatePolynomial<BigInteger> a,
                                                   MultivariatePolynomial<BigInteger> b,
                                                   MultivariatePolynomial<BigInteger> skeleton,
                                                   RandomGenerator rnd) {
        boolean monic = a.coefficientOf(0, a.degree(0)).isConstant() && b.coefficientOf(0, a.degree(0)).isConstant();

        Set<DegreeVector> globalSkeleton = skeleton.data.keySet();
        TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        int[] evaluationVariables = ArraysUtil.sequence(1, variable + 1);//variable inclusive
        BigInteger[] evaluationPoint = new BigInteger[evaluationVariables.length];

        Domain<BigInteger> domain = a.domain;
        //avoid zero evaluation points // TODO check the same skeleton in main variable!!!
        for (int i = variable - 2; i >= 0; --i)
            do {
                evaluationPoint[i] = domain.randomElement(rnd);
            } while (domain.isZero(evaluationPoint[i]));

        PrecomputedPowersHolder<BigInteger> powers = new PrecomputedPowersHolder<>(evaluationPoint, domain);

        int requiredNumberOfEvaluations = -1, monicScalingExponent = -1;
        for (TIntObjectIterator<MultivariatePolynomial<BigInteger>> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            MultivariatePolynomial<BigInteger> v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
            if (v.size() == 1)
                monicScalingExponent = it.key();
        }

        if (!ALWAYS_LINZIP) {
            if (monic)
                monicScalingExponent = -1;

            if (monic || monicScalingExponent != -1)
                return new MonicInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                        evaluationVariables, evaluationPoint, powers, rnd, requiredNumberOfEvaluations, monicScalingExponent);
        }

        return new LinZipInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, rnd);
    }

    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
    private static TIntObjectHashMap<MultivariatePolynomial<BigInteger>> getSkeleton(MultivariatePolynomial<BigInteger> poly) {
        TIntObjectHashMap<MultivariatePolynomial<BigInteger>> skeleton = new TIntObjectHashMap<>();
        for (MonomialTerm<BigInteger> term : poly.data) {
            MonomialTerm<BigInteger> newDV = term.setZero(0);
            MultivariatePolynomial<BigInteger> coeff = skeleton.get(term.exponents[0]);
            if (coeff != null)
                coeff.add(newDV);
            else
                skeleton.put(term.exponents[0], MultivariatePolynomial.create(
                        poly.domain, poly.ordering, newDV));
        }
        return skeleton;
    }

    static abstract class SparseInterpolation {
        /** the domain */
        final Domain<BigInteger> domain;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final MultivariatePolynomial<BigInteger> a, b;
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

        SparseInterpolation(Domain<BigInteger> domain, int variable,
                            MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b,
                            Set<DegreeVector> globalSkeleton,
                            TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton,
                            int[] sparseUnivarDegrees, int[] evaluationVariables,
                            BigInteger[] evaluationPoint,
                            PrecomputedPowersHolder<BigInteger> powers, RandomGenerator rnd) {
            this.domain = domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.globalSkeleton = globalSkeleton;
            this.univarSkeleton = univarSkeleton;
            this.sparseUnivarDegrees = sparseUnivarDegrees;
            this.evaluationVariables = evaluationVariables;
            this.evaluationPoint = evaluationPoint;
            this.powers = powers;
            this.rnd = rnd;
        }

        abstract MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint);
    }

    static final class LinZipInterpolation extends SparseInterpolation {
        LinZipInterpolation(Domain<BigInteger> domain, int variable, MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, BigInteger[] evaluationPoint, PrecomputedPowersHolder<BigInteger> powers, RandomGenerator rnd) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
        }

        /** Number of retries when raise condition occurs; then drop up with new homomorphism */
        private static final int NUMBER_OF_UNDER_DETERMINED_RETRIES = 100;

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        @Override
        MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            LinZipSystem[] systems = new LinZipSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new LinZipSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            // the last variable (the variable) is the same for all evaluations = newPoint
            raiseFactors[raiseFactors.length - 1] = 1;

            int nUnknowns = globalSkeleton.size(), nUnknownScalings = -1;
            int raiseFactor = 0;

            for (int nTries = 0; nTries < NUMBER_OF_UNDER_DETERMINED_RETRIES; ++nTries) {
                int previousFreeVars = -1, underDeterminedTries = 0;
                for (; ; ) {
                    // increment at each loop!
                    ++nUnknownScalings;
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
                    for (LinZipSystem system : systems) {
                        BigInteger rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs, nUnknownScalings != 0);
                        totalEquations += system.nEquations();
                    }
                    if (nUnknowns + nUnknownScalings <= totalEquations)
                        break;


                    if (underDeterminedTries > NUMBER_OF_UNDER_DETERMINED_RETRIES) {
                        // raise condition: new equations does not fix enough variables
                        return null;
                    }

                    int freeVars = nUnknowns + nUnknownScalings - totalEquations;
                    if (freeVars >= previousFreeVars)
                        ++underDeterminedTries;
                    else
                        underDeterminedTries = 0;

                    previousFreeVars = freeVars;
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
            return null;
        }
    }

    private static SystemInfo solveLinZip(MultivariatePolynomial<BigInteger> factory, LinZipSystem[] subSystems, int nUnknownScalings, MultivariatePolynomial<BigInteger> destination) {
        ArrayList<MonomialTerm<BigInteger>> unknowns = new ArrayList<>();
        for (LinZipSystem system : subSystems)
            for (MonomialTerm<BigInteger> degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = nUnknownsMonomials + nUnknownScalings;
        ArrayList<BigInteger[]> lhsGlobal = new ArrayList<>();
        ArrayList<BigInteger> rhsGlobal = new ArrayList<>();
        int offset = 0;
        for (LinZipSystem system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                BigInteger[] row = new BigInteger[nUnknownsTotal];
                Arrays.fill(row, BigInteger.ZERO);
                BigInteger[] subRow = system.matrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                if (j > 0)
                    row[nUnknownsMonomials + j - 1] = system.scalingMatrix.get(j);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        BigInteger[] solution = new BigInteger[nUnknownsTotal];
        SystemInfo info = LinearAlgebra.solve(factory.domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent) {
            MonomialTerm<BigInteger>[] terms = new MonomialTerm[unknowns.size()];
            for (int i = 0; i < terms.length; i++)
                terms[i] = unknowns.get(i).setCoefficient(solution[i]);
            destination.add(terms);
        }
        return info;
    }


    static final class MonicInterpolation extends SparseInterpolation {
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** univar exponent with monomial factor that can be used for scaling */
        final int monicScalingExponent;

        MonicInterpolation(Domain<BigInteger> domain, int variable, MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<BigInteger>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, BigInteger[] evaluationPoint, PrecomputedPowersHolder<BigInteger> powers, RandomGenerator rnd, int requiredNumberOfEvaluations, int monicScalingExponent) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
            this.requiredNumberOfEvaluations = requiredNumberOfEvaluations;
            this.monicScalingExponent = monicScalingExponent;
        }

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        @Override
        MultivariatePolynomial<BigInteger> evaluate(BigInteger newPoint) {
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

                    BigInteger normalization = evaluateExceptFirst(domain, powers, domain.getOne(), univarSkeleton.get(monicScalingExponent).lt(), i + 1, variable - 1);
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
                SystemInfo info = system.solve();
                if (info != SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomial<BigInteger> gcdVal = a.createZero();
            for (VandermondeSystem system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || domain.isOne(system.solution[0]);
                for (int i = 0; i < system.skeleton.length; i++) {
                    MonomialTerm<BigInteger> degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    BigInteger value = system.solution[i];
                    gcdVal.add(degreeVector.setCoefficient(value));
                }
            }

            return gcdVal;
        }
    }

    private static abstract class LinearSystem {
        final int univarDegree;
        /** the domain */
        final Domain<BigInteger> domain;
        /** the skeleton */
        final MonomialTerm<BigInteger>[] skeleton;
        /** the lhs matrix */
        final ArrayList<BigInteger[]> matrix;
        /** the rhs values */
        final ArrayList<BigInteger> rhs = new ArrayList<>();
        /** precomputed powers */
        final PrecomputedPowersHolder<BigInteger> powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        final int nVars;

        LinearSystem(int univarDegree, MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            //todo refactor generics
            this.skeleton = skeleton.getSkeleton().toArray(new MonomialTerm[skeleton.size()]);
            this.powers = powers;
            this.nVars = nVars;
            this.matrix = new ArrayList<>();
        }

        final int nUnknownVariables() {
            return skeleton.length;
        }

        final int nEquations() {
            return matrix.size();
        }

        @Override
        public String toString() {
            return "{" + matrix.stream().map(Arrays::toString).collect(Collectors.joining(",")) + "} = " + rhs;
        }
    }

    /** Vandermonde system builder */
    private static final class LinZipSystem extends LinearSystem {
        public LinZipSystem(int univarDegree, MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        private final ArrayList<BigInteger> scalingMatrix = new ArrayList<>();

        public void oneMoreEquation(BigInteger rhsVal, boolean newScalingIntroduced) {
            BigInteger[] row = new BigInteger[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);

            if (newScalingIntroduced) {
                scalingMatrix.add(domain.negate(rhsVal));
                rhsVal = domain.getZero();
            } else
                scalingMatrix.add(domain.getZero());
            rhs.add(rhsVal);
        }

    }

    /** Vandermonde system builder */
    private static final class VandermondeSystem extends LinearSystem {
        public VandermondeSystem(int univarDegree, MultivariatePolynomial<BigInteger> skeleton, PrecomputedPowersHolder<BigInteger> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        BigInteger[] solution = null;

        SystemInfo solve() {
            if (solution == null)
                solution = new BigInteger[nUnknownVariables()];

            if (nUnknownVariables() <= 8)
                // for small systems Gaussian elimination is indeed faster
                return LinearAlgebra.solve(domain, matrix.toArray(new BigInteger[matrix.size()][]), rhs.toArray(new BigInteger[rhs.size()]), solution);

            // solve vandermonde system
            BigInteger[] vandermondeRow = matrix.get(0);
            SystemInfo info = LinearAlgebra.solveVandermondeT((ModularDomain) domain, vandermondeRow, rhs.toArray(new BigInteger[rhs.size()]), solution);
            if (info == SystemInfo.Consistent)
                for (int i = 0; i < solution.length; ++i)
                    solution[i] = domain.divideExact(solution[i], vandermondeRow[i]);

            return info;
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
                                                  MonomialTerm<BigInteger> skeleton,
                                                  int raiseFactor,
                                                  int nVars) {
        BigInteger tmp = coefficient;
        for (int k = 0; k < nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[1 + k]));
        return tmp;
    }

    private static boolean isVandermonde(List<BigInteger[]> lhs, Domain<BigInteger> domain) {
        return isVandermonde(lhs.toArray(new BigInteger[0][]), domain);
    }

    private static boolean isVandermonde(BigInteger[][] lhs, Domain<BigInteger> domain) {
        for (int i = 1; i < lhs.length; i++) {
            for (int j = 0; j < lhs[0].length; j++) {
                if (!lhs[i][j].equals(domain.pow(lhs[0][j], i + 1)))
                    return false;
            }
        }
        return true;
    }
}
