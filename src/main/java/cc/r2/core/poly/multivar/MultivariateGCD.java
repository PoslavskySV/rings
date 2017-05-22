package cc.r2.core.poly.multivar;


import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.LinearAlgebra.SystemInfo;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.multivar.MultivariateInterpolation.lInterpolation;
import cc.r2.core.poly.multivar.lMultivariatePolynomial.lPrecomputedPowersHolder;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TLongHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.stream.Collectors;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.PrecomputedPowersHolder;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.renameVariables;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.univar.UnivariateGCD.PolynomialGCD;
import static cc.r2.core.util.ArraysUtil.negate;


/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateGCD {
    private MultivariateGCD() {}

    /** calculates the inverse permutation */
    public static int[] inversePermutation(int[] permutation) {
        final int[] inv = new int[permutation.length];
        for (int i = permutation.length - 1; i >= 0; --i)
            inv[permutation[i]] = i;
        return inv;
    }

    /** structure with required input for GCD algorithms */
    private static final class GCDInput<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** input polynomials (with variables renamed) and earlyGCD if possible (in trivial cases) */
        final Poly aReduced, bReduced, earlyGCD;
        /** gcd degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** domain cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int evaluationStackLimit;
        /** GCD of monomial content of a and b **/
        final Term monomialGCD;

        GCDInput(Poly earlyGCD) {
            this.earlyGCD = earlyGCD;
            aReduced = bReduced = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialGCD = null;
        }

        GCDInput(Poly aReduced, Poly bReduced, Term monomialGCD,
                 int evaluationStackLimit, int[] degreeBounds, int[] mapping, int lastPresentVariable) {
            //assert monomialGCD == null || aReduced.domain.isOne(monomialGCD.coefficient);
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
        Poly restoreGCD(Poly result) {
            return renameVariables(result, mapping).multiply(monomialGCD);
        }
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly trivialGCD(Poly a, Poly b) {
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
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    GCDInput<Term, Poly> preparedGCDInput(Poly a, Poly b) {
        Poly trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return new GCDInput<>(trivialGCD);

        BigInteger domainSize = a.coefficientDomainCardinality();
        if (domainSize == null)
            throw new IllegalArgumentException("Modular gcd algorithms are supported only for multivariate polynomial over finite fields.");

        if (domainSize.isInt() && Math.min(a.degree(), b.degree()) > domainSize.intValue())
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        // domain cardinality, i.e. number of possible random choices
        int evaluationStackLimit = domainSize.isInt() ? domainSize.intValue() : -1;

        // find monomial GCD
        // and remove monomial content from a and b
        a = a.clone(); b = b.clone(); // prevent rewriting original data
        Term monomialGCD = reduceMonomialContent(a, b);

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
            return new GCDInput<>(a.create(monomialGCD));

        RandomGenerator rnd = PrivateRandom.getRandom();
        if (nUsedVariables != nVariables) {
            // some of the variables are redundant in one of the polys => can just substitute a random values for them
            for (int i = 0; i < nVariables; i++) {
                if (degreeBounds[i] == 0) {
                    if (a.degree(i) != 0) {
                        assert b.degree(i) == 0;
                        a = a.evaluateAtRandomPreservingSkeleton(i, rnd);
                    } else if (b.degree(i) != 0) {
                        assert a.degree(i) == 0;
                        b = b.evaluateAtRandomPreservingSkeleton(i, rnd);
                    }
                }
            }
        }

        if (nUsedVariables == 1)
        // switch to univariate gcd
        {
            @SuppressWarnings("unchecked")
            IUnivariatePolynomial iUnivar = PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            Poly poly = asMultivariate(iUnivar, nVariables, lastPresentVariable, a.ordering);
            return new GCDInput<>(poly.multiply(monomialGCD));
        }

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

        return new GCDInput<>(a, b, monomialGCD, evaluationStackLimit, degreeBounds, variables, lastPresentVariable);
    }

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly asMultivariate(IUnivariatePolynomial poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        if (poly instanceof UnivariatePolynomial)
            return (Poly) MultivariatePolynomial.asMultivariate((UnivariatePolynomial) poly, nVariables, variable, ordering);
        else if (poly instanceof lUnivariatePolynomialZp)
            return (Poly) lMultivariatePolynomial.asMultivariate((lUnivariatePolynomialZp) poly, nVariables, variable, ordering);
        else
            throw new RuntimeException();
    }

    /** gcd with monomial */
    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly gcdWithMonomial(Term monomial, Poly poly) {
        return poly.create(poly.commonContent(monomial));
    }

    /**
     * Removes monomial content from a and b and returns monomial gcd
     */
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Term reduceMonomialContent(Poly a, Poly b) {

        Term aMonomialContent = a.monomialContent();
        int[] exponentsGCD = b.monomialContent().exponents;
        AMultivariatePolynomial.setMin(aMonomialContent.exponents, exponentsGCD);
        Term monomialGCD = a.createTermWithUnitCoefficient(exponentsGCD);

        a = a.divideDegreeVectorOrNull(monomialGCD);
        b = b.divideDegreeVectorOrNull(monomialGCD);
        assert a != null && b != null;

        return monomialGCD;
    }

    /**
     * Primitive part and content of multivariate polynomial considered as polynomial over Zp[x_i][x_1, ..., x_N]
     */
    private static final class UnivariateContent<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        final MultivariatePolynomial<uPoly> poly;
        final Poly primitivePart;
        final uPoly content;

        public UnivariateContent(MultivariatePolynomial<uPoly> poly, Poly primitivePart, uPoly content) {
            this.poly = poly;
            this.primitivePart = primitivePart;
            this.content = content;
        }
    }

    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    MultivariatePolynomial<uPoly> asOverUnivariate(Poly poly, int variable) {
        if (poly instanceof MultivariatePolynomial)
            return (MultivariatePolynomial<uPoly>) ((MultivariatePolynomial) poly).asOverUnivariate(variable);
        else if (poly instanceof lMultivariatePolynomial)
            return (MultivariatePolynomial<uPoly>) ((lMultivariatePolynomial) poly).asOverUnivariate(variable);
        else
            throw new RuntimeException();
    }

    /**
     * Returns primitive part and content of {@code poly} considered as polynomial over Zp[variable][x_1, ..., x_N]
     *
     * @param poly     the polynomial
     * @param variable the variable
     * @return primitive part and content
     */
    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariateContent<Term, Poly, uPoly> univariateContent(
            Poly poly, int variable) {
        //convert poly to Zp[var][x...]

        MultivariatePolynomial<uPoly> conv = asOverUnivariate(poly, variable);
        //univariate content
        uPoly content = PolynomialGCD(conv.coefficients());
        Poly divider = asMultivariate(content, poly.nVariables, variable, poly.ordering);
        Poly[] qd = divideAndRemainder(poly, divider);
        assert qd[1].isZero();
        return new UnivariateContent(conv, qd[0], content);
    }

    /** holds primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
    private static final class PrimitiveInput<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        /** primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
        final Poly aPrimitive, bPrimitive;
        /** gcd of content and leading coefficient of a and b given as Zp[x_k][x_1 ... x_{k-1}] */
        final uPoly contentGCD, lcGCD;

        public PrimitiveInput(Poly aPrimitive, Poly bPrimitive,
                              uPoly contentGCD, uPoly lcGCD) {
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
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    PrimitiveInput<Term, Poly, uPoly> makePrimitive(Poly a, Poly b, int variable) {
        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        UnivariateContent<Term, Poly, uPoly>
                aContent = univariateContent(a, variable),
                bContent = univariateContent(b, variable);

        a = aContent.primitivePart;
        b = bContent.primitivePart;

        // gcd of Zp[x_k] content and lc
        uPoly
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
    public static <E> MultivariatePolynomial<E> BrownGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b) {

        // prepare input and test for early termination
        GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>> gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        MultivariatePolynomial<E> result = BrownGCD(
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
    private static <E> MultivariatePolynomial<E> BrownGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        MultivariatePolynomial<E> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        MultivariatePolynomial<E> factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            UnivariatePolynomial<E> gcd = PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        PrimitiveInput<MonomialTerm<E>, MultivariatePolynomial<E>, UnivariatePolynomial<E>> primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        UnivariatePolynomial<E>
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            MultivariatePolynomial<E> poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        Domain<E> domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //dense interpolation
        Interpolation<E> interpolation = null;
        //previous interpolation (used to detect whether update doesn't change the result)
        MultivariatePolynomial<E> previousInterpolation;
        //store points that were already used in interpolation
        Set<E> evaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // all elements of the domain are tried
                // do division check (last chance) and return
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            //pickup the next random element for variable
            E randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;
            evaluationStack.add(randomPoint);

            E lcVal = lcGCD.evaluate(randomPoint);
            if (domain.isZero(lcVal))
                continue;

            // evaluate a and b at variable = randomPoint
            MultivariatePolynomial<E>
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<E> cVal = BrownGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
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

            // Cache previous interpolation. NOTE: clone() is important, since the poly will
            // be modified inplace by the update() method
            previousInterpolation = interpolation.getInterpolatingPolynomial().clone();
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] <= interpolation.numberOfPoints()
                    || previousInterpolation.equals(interpolation.getInterpolatingPolynomial())) {
                MultivariatePolynomial<E> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    /** division test **/
    private static <E> MultivariatePolynomial<E> doDivisionCheck(
            MultivariatePolynomial<E> a, MultivariatePolynomial<E> b,
            UnivariatePolynomial<E> contentGCD, Interpolation<E> interpolation, int variable) {
        if (interpolation == null)
            return null;
        MultivariatePolynomial<E> interpolated =
                MultivariatePolynomial.asNormalMultivariate(interpolation.getInterpolatingPolynomial().asOverUnivariate(variable).primitivePart(), variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        if (contentGCD == null)
            return interpolated;
        MultivariatePolynomial<E> poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
        return interpolated.multiply(poly);
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Zippel's algorithm with sparse interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static <E> MultivariatePolynomial<E> ZippelGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b) {

        // prepare input and test for early termination
        GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>> gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

//        MultivariatePolynomial<E> content = ZippelContentGCD(a, b, 0);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        MultivariatePolynomial<E> result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return gcdInput.restoreGCD(result);
    }

    @SuppressWarnings("unchecked")
    private static <E> MultivariatePolynomial<E>[] multivariateCoefficients(
            MultivariatePolynomial<E> a,
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

    static <E> MultivariatePolynomial<E> ZippelContentGCD(
            MultivariatePolynomial<E> a, MultivariatePolynomial<E> b,
            int variable) {
        MultivariatePolynomial<E>[] aCfs = multivariateCoefficients(a, variable);
        MultivariatePolynomial<E>[] bCfs = multivariateCoefficients(b, variable);

        MultivariatePolynomial<E> contentGCD = ZippelGCD(aCfs[0], bCfs[0]);
        for (int i = 1; i < aCfs.length; i++)
            contentGCD = ZippelGCD(aCfs[i], contentGCD);
        for (int i = 1; i < bCfs.length; i++)
            contentGCD = ZippelGCD(bCfs[i], contentGCD);
        return contentGCD;
    }

    /** Maximal number of fails before switch to a new homomorphism */
    private static final int MAX_SPARSE_INTERPOLATION_FAILS = 1000;
    /** Maximal number of sparse interpolations after interpolation.numberOfPoints() > degreeBounds[variable]  */
    private static final int ALLOWED_OVER_INTERPOLATED_ATTEMPTS = 32;

    @SuppressWarnings("unchecked")
    private static <E> MultivariatePolynomial<E> ZippelGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        MultivariatePolynomial<E> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        //for bivariate polynomials always use dense interpolation
        if (variable == 1)
            return BrownGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit);

//        MultivariatePolynomial<E> content = ZippelContentGCD(a, b, variable);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        MultivariatePolynomial<E> factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            UnivariatePolynomial<E> gcd = PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }
        PrimitiveInput<
                MonomialTerm<E>,
                MultivariatePolynomial<E>,
                UnivariatePolynomial<E>>
                primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        UnivariatePolynomial<E>
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            MultivariatePolynomial<E> poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        Domain<E> domain = factory.domain;
        //store points that were already used in interpolation
        Set<E> globalEvaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        int failedSparseInterpolations = 0;

        int[] tmpDegreeBounds = degreeBounds.clone();
        main:
        while (true) {
            if (evaluationStackLimit == globalEvaluationStack.size())
                return null;

            E seedPoint = domain.randomElement(rnd);
            if (globalEvaluationStack.contains(seedPoint))
                continue;

            globalEvaluationStack.add(seedPoint);

            E lcVal = lcGCD.evaluate(seedPoint);
            if (domain.isZero(lcVal))
                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<E>
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<E> cVal = ZippelGCD(aVal, bVal, rnd, variable - 1, tmpDegreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > tmpDegreeBounds[variable - 1])
                //unlucky homomorphism
                continue;

            if (currExponent < tmpDegreeBounds[variable - 1]) {
                //better degree bound detected
                tmpDegreeBounds[variable - 1] = currExponent;
            }

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            SparseInterpolation sparseInterpolator = createInterpolation(variable, a, b, cVal, rnd);

            // we are applying dense interpolation for univariate skeleton coefficients
            Interpolation<E> denseInterpolation = new Interpolation<>(variable, seedPoint, cVal);
            //previous interpolation (used to detect whether update doesn't change the result)
            MultivariatePolynomial<E> previousInterpolation;
            //local evaluation stack for points that are calculated via sparse interpolation (but not gcd evaluation) -> always same skeleton
            HashSet<E> localEvaluationStack = new HashSet<>(globalEvaluationStack);
            while (true) {
                if (denseInterpolation.numberOfPoints() > tmpDegreeBounds[variable] + ALLOWED_OVER_INTERPOLATED_ATTEMPTS) {
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                E randomPoint = domain.randomElement(rnd);
                if (localEvaluationStack.contains(randomPoint))
                    continue;

                lcVal = lcGCD.evaluate(randomPoint);
                if (domain.isZero(lcVal))
                    continue;

                localEvaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
                if (cVal == null) {
                    ++failedSparseInterpolations;
                    if (failedSparseInterpolations == MAX_SPARSE_INTERPOLATION_FAILS)
                        throw new RuntimeException("Sparse interpolation failed");
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc().equals(lcVal);

                // Cache previous interpolation. NOTE: clone() is important, since the poly will
                // be modified inplace by the update() method
                previousInterpolation = denseInterpolation.getInterpolatingPolynomial().clone();
                denseInterpolation.update(randomPoint, cVal);

                // do division test
                if (tmpDegreeBounds[variable] <= denseInterpolation.numberOfPoints()
                        || previousInterpolation.equals(denseInterpolation.getInterpolatingPolynomial())) {
                    MultivariatePolynomial<E> result = doDivisionCheck(a, b, contentGCD, denseInterpolation, variable);
                    if (result != null)
                        return result;
                }
            }
        }
    }

    static boolean ALWAYS_LINZIP = false;

    static <E> SparseInterpolation<E> createInterpolation(int variable,
                                                          MultivariatePolynomial<E> a,
                                                          MultivariatePolynomial<E> b,
                                                          MultivariatePolynomial<E> skeleton,
                                                          RandomGenerator rnd) {
        boolean monic = a.coefficientOf(0, a.degree(0)).isConstant() && b.coefficientOf(0, a.degree(0)).isConstant();

        Set<DegreeVector> globalSkeleton = skeleton.terms.keySet();
        TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        Domain<E> domain = a.domain;

        int[] evaluationVariables = ArraysUtil.sequence(1, variable + 1);//variable inclusive
        E[] evaluationPoint = domain.createArray(evaluationVariables.length);

        PrecomputedPowersHolder<E> powers;
        search_for_good_evaluation_point:
        while (true) {
            //avoid zero evaluation points
            for (int i = variable - 2; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (domain.isZero(evaluationPoint[i]));

            //set temporary
            evaluationPoint[evaluationPoint.length - 1] = domain.randomElement(rnd);
            powers = new PrecomputedPowersHolder<>(evaluationPoint, domain);
            int[] raiseFactors = ArraysUtil.arrayOf(1, evaluationVariables.length);

            for (MultivariatePolynomial<E> p : Arrays.asList(a, b, skeleton))
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables, raiseFactors).getSkeleton()))
                    continue search_for_good_evaluation_point;
            break;
        }

        int requiredNumberOfEvaluations = -1, monicScalingExponent = -1;
        for (TIntObjectIterator<MultivariatePolynomial<E>> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            MultivariatePolynomial<E> v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
            if (v.size() == 1)
                monicScalingExponent = it.key();
        }

        if (!ALWAYS_LINZIP) {
            if (monic)
                monicScalingExponent = -1;

            if (monic || monicScalingExponent != -1)
                return new MonicInterpolation<>(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                        evaluationVariables, evaluationPoint, powers, rnd, requiredNumberOfEvaluations, monicScalingExponent);
        }

        return new LinZipInterpolation<>(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, rnd);
    }

    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
    private static <E> TIntObjectHashMap<MultivariatePolynomial<E>> getSkeleton(MultivariatePolynomial<E> poly) {
        TIntObjectHashMap<MultivariatePolynomial<E>> skeleton = new TIntObjectHashMap<>();
        for (MonomialTerm<E> term : poly.terms) {
            MonomialTerm<E> newDV = term.setZero(0);
            MultivariatePolynomial<E> coeff = skeleton.get(term.exponents[0]);
            if (coeff != null)
                coeff.add(newDV);
            else
                skeleton.put(term.exponents[0], MultivariatePolynomial.create(poly.nVariables,
                        poly.domain, poly.ordering, newDV));
        }
        return skeleton;
    }

    static abstract class SparseInterpolation<E> {
        /** the domain */
        final Domain<E> domain;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final MultivariatePolynomial<E> a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final E[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final PrecomputedPowersHolder<E> powers;
        /** random */
        final RandomGenerator rnd;

        SparseInterpolation(Domain<E> domain, int variable,
                            MultivariatePolynomial<E> a, MultivariatePolynomial<E> b,
                            Set<DegreeVector> globalSkeleton,
                            TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton,
                            int[] sparseUnivarDegrees, int[] evaluationVariables,
                            E[] evaluationPoint,
                            PrecomputedPowersHolder<E> powers, RandomGenerator rnd) {
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

        abstract MultivariatePolynomial<E> evaluate(E newPoint);
    }

    static final class LinZipInterpolation<E> extends SparseInterpolation<E> {
        LinZipInterpolation(Domain<E> domain, int variable, MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, E[] evaluationPoint, PrecomputedPowersHolder<E> powers, RandomGenerator rnd) {
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
        MultivariatePolynomial<E> evaluate(E newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            @SuppressWarnings("unchecked")
            LinZipSystem<E>[] systems = new LinZipSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new LinZipSystem<>(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


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
                    UnivariatePolynomial<E>
                            aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                    if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                        // unlucky main homomorphism or bad evaluation point
                        return null;

                    assert gcdUnivar.isMonic();

                    int totalEquations = 0;
                    for (LinZipSystem<E> system : systems) {
                        E rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
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

                MultivariatePolynomial<E> result = a.createZero();
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

    private static <E> SystemInfo solveLinZip(MultivariatePolynomial<E> factory, LinZipSystem<E>[] subSystems, int nUnknownScalings, MultivariatePolynomial<E> destination) {
        ArrayList<MonomialTerm<E>> unknowns = new ArrayList<>();
        for (LinZipSystem<E> system : subSystems)
            for (MonomialTerm<E> degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = nUnknownsMonomials + nUnknownScalings;
        ArrayList<E[]> lhsGlobal = new ArrayList<>();
        ArrayList<E> rhsGlobal = new ArrayList<>();
        int offset = 0;
        Domain<E> domain = factory.domain;
        for (LinZipSystem<E> system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                E[] row = domain.createZeroesArray(nUnknownsTotal);
                E[] subRow = system.matrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                if (j > 0)
                    row[nUnknownsMonomials + j - 1] = system.scalingMatrix.get(j);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        E[] solution = domain.createArray(nUnknownsTotal);
        SystemInfo info = LinearAlgebra.solve(domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent) {
            @SuppressWarnings("unchecked")
            MonomialTerm<E>[] terms = new MonomialTerm[unknowns.size()];
            for (int i = 0; i < terms.length; i++)
                terms[i] = unknowns.get(i).setCoefficient(solution[i]);
            destination.add(terms);
        }
        return info;
    }


    static final class MonicInterpolation<E> extends SparseInterpolation<E> {
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** univar exponent with monomial factor that can be used for scaling */
        final int monicScalingExponent;

        MonicInterpolation(Domain<E> domain, int variable, MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, E[] evaluationPoint, PrecomputedPowersHolder<E> powers, RandomGenerator rnd, int requiredNumberOfEvaluations, int monicScalingExponent) {
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
        MultivariatePolynomial<E> evaluate(E newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            @SuppressWarnings("unchecked")
            VandermondeSystem<E>[] systems = new VandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new VandermondeSystem<>(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                Arrays.fill(raiseFactors, 0, raiseFactors.length - 1, i + 1);
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;
                // evaluate a and b to univariate and calculate gcd
                UnivariatePolynomial<E>
                        aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                    // unlucky main homomorphism or bad evaluation point
                    return null;

                assert gcdUnivar.isMonic();

                if (monicScalingExponent != -1) {
                    // single scaling factor
                    // scale the system according to it

                    if (gcdUnivar.degree() < monicScalingExponent || domain.isZero(gcdUnivar.get(monicScalingExponent)))
                        // unlucky homomorphism
                        return null;

                    E normalization = evaluateExceptFirst(domain, powers, domain.getOne(), univarSkeleton.get(monicScalingExponent).lt(), i + 1, variable - 1);
                    //normalize univariate gcd in order to reconstruct leading coefficient polynomial
                    normalization = domain.multiply(domain.reciprocal(gcdUnivar.get(monicScalingExponent)), normalization);
                    gcdUnivar = gcdUnivar.multiply(normalization);
                }

                boolean allDone = true;
                for (VandermondeSystem<E> system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        E rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (VandermondeSystem<E> system : systems) {
                //solve each system
                SystemInfo info = system.solve();
                if (info != SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomial<E> gcdVal = a.createZero();
            for (VandermondeSystem<E> system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || domain.isOne(system.solution[0]);
                for (int i = 0; i < system.skeleton.length; i++) {
                    MonomialTerm<E> degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    E value = system.solution[i];
                    gcdVal.add(degreeVector.setCoefficient(value));
                }
            }

            return gcdVal;
        }
    }

    private static abstract class LinearSystem<E> {
        final int univarDegree;
        /** the domain */
        final Domain<E> domain;
        /** the skeleton */
        final MonomialTerm<E>[] skeleton;
        /** the lhs matrix */
        final ArrayList<E[]> matrix;
        /** the rhs values */
        final ArrayList<E> rhs = new ArrayList<>();
        /** precomputed powers */
        final PrecomputedPowersHolder<E> powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        final int nVars;

        @SuppressWarnings("unchecked")
        LinearSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
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
    private static final class LinZipSystem<E> extends LinearSystem<E> {
        public LinZipSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        private final ArrayList<E> scalingMatrix = new ArrayList<>();

        public void oneMoreEquation(E rhsVal, boolean newScalingIntroduced) {
            E[] row = domain.createArray(skeleton.length);
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
    private static final class VandermondeSystem<E> extends LinearSystem<E> {
        public VandermondeSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        E[] solution = null;

        SystemInfo solve() {
            if (solution == null)
                solution = domain.createArray(nUnknownVariables());

            if (nUnknownVariables() <= 8)
                // for small systems Gaussian elimination is indeed faster
                return LinearAlgebra.solve(domain, matrix.toArray(domain.createArray2d(matrix.size())), rhs.toArray(domain.createArray(rhs.size())), solution);

            // solve vandermonde system
            E[] vandermondeRow = matrix.get(0);
            SystemInfo info = LinearAlgebra.solveVandermondeT(domain, vandermondeRow, rhs.toArray(domain.createArray(rhs.size())), solution);
            if (info == SystemInfo.Consistent)
                for (int i = 0; i < solution.length; ++i)
                    solution[i] = domain.divideExact(solution[i], vandermondeRow[i]);

            return info;
        }

        public VandermondeSystem<E> oneMoreEquation(E rhsVal) {
            E[] row = domain.createArray(skeleton.length);
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);
            rhs.add(rhsVal);
            return this;
        }
    }

    private static <E> E evaluateExceptFirst(Domain<E> domain,
                                             PrecomputedPowersHolder<E> powers,
                                             E coefficient,
                                             MonomialTerm<E> skeleton,
                                             int raiseFactor,
                                             int nVars) {
        E tmp = coefficient;
        for (int k = 0; k < nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[1 + k]));
        return tmp;
    }

    private static <E> boolean isVandermonde(E[][] lhs, Domain<E> domain) {
        for (int i = 1; i < lhs.length; i++) {
            for (int j = 0; j < lhs[0].length; j++) {
                if (!lhs[i][j].equals(domain.pow(lhs[0][j], i + 1)))
                    return false;
            }
        }
        return true;
    }


    /* ========================================= Machine numbers ======================================== */

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Brown's algorithm with dense interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static lMultivariatePolynomial BrownGCD(
            lMultivariatePolynomial a,
            lMultivariatePolynomial b) {

        // prepare input and test for early termination
        GCDInput<lMonomialTerm, lMultivariatePolynomial> gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        lMultivariatePolynomial result = BrownGCD(
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
    private static lMultivariatePolynomial BrownGCD(
            lMultivariatePolynomial a,
            lMultivariatePolynomial b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        lMultivariatePolynomial trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        lMultivariatePolynomial factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            lUnivariatePolynomialZp gcd = PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        PrimitiveInput<lMonomialTerm, lMultivariatePolynomial, lUnivariatePolynomialZp> primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        lUnivariatePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            lMultivariatePolynomial poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        lIntegersModulo domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //dense interpolation
        lInterpolation interpolation = null;
        //previous interpolation (used to detect whether update doesn't change the result)
        lMultivariatePolynomial previousInterpolation;
        //store points that were already used in interpolation
        TLongHashSet evaluationStack = new TLongHashSet();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // all elements of the domain are tried
                // do division check (last chance) and return
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            //pickup the next random element for variable
            long randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;
            evaluationStack.add(randomPoint);

            long lcVal = lcGCD.evaluate(randomPoint);
            if (lcVal == 0)
                continue;

            // evaluate a and b at variable = randomPoint
            lMultivariatePolynomial
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            // calculate gcd of the result by the recursive call
            lMultivariatePolynomial cVal = BrownGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            // normalize gcd
            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc() == lcVal;

            if (currExponent < prevVarExponent) {
                //better degree bound detected => start over
                interpolation = new lInterpolation(variable, randomPoint, cVal);
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
                continue;
            }

            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new lInterpolation(variable, randomPoint, cVal);
                continue;
            }

            // Cache previous interpolation. NOTE: clone() is important, since the poly will
            // be modified inplace by the update() method
            previousInterpolation = interpolation.getInterpolatingPolynomial().clone();
            // interpolate
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] <= interpolation.numberOfPoints()
                    || previousInterpolation.equals(interpolation.getInterpolatingPolynomial())) {
                lMultivariatePolynomial result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    /** division test **/
    private static lMultivariatePolynomial doDivisionCheck(
            lMultivariatePolynomial a, lMultivariatePolynomial b,
            lUnivariatePolynomialZp contentGCD, lInterpolation interpolation, int variable) {
        if (interpolation == null)
            return null;
        lMultivariatePolynomial interpolated =
                lMultivariatePolynomial.asNormalMultivariate(interpolation.getInterpolatingPolynomial().asOverUnivariate(variable).primitivePart(), variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        if (contentGCD == null)
            return interpolated;
        lMultivariatePolynomial poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
        return interpolated.multiply(poly);
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Zippel's algorithm with sparse interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static lMultivariatePolynomial ZippelGCD(
            lMultivariatePolynomial a,
            lMultivariatePolynomial b) {

        // prepare input and test for early termination
        GCDInput<lMonomialTerm, lMultivariatePolynomial> gcdInput = preparedGCDInput(a, b);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

//        lMultivariatePolynomial content = ZippelContentGCD(a, b, 0);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        lMultivariatePolynomial result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            throw new IllegalArgumentException("Ground fill is too small for modular algorithm.");
        return gcdInput.restoreGCD(result);
    }

    @SuppressWarnings("unchecked")
    private static lMultivariatePolynomial[] multivariateCoefficients(
            lMultivariatePolynomial a,
            int variable) {
        lMultivariatePolynomial[] mCfs = Arrays.stream(a.degrees(variable)).mapToObj(d -> a.coefficientOf(variable, d)).toArray(lMultivariatePolynomial[]::new);
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

    static lMultivariatePolynomial ZippelContentGCD(
            lMultivariatePolynomial a, lMultivariatePolynomial b,
            int variable) {
        lMultivariatePolynomial[] aCfs = multivariateCoefficients(a, variable);
        lMultivariatePolynomial[] bCfs = multivariateCoefficients(b, variable);

        lMultivariatePolynomial contentGCD = ZippelGCD(aCfs[0], bCfs[0]);
        for (int i = 1; i < aCfs.length; i++)
            contentGCD = ZippelGCD(aCfs[i], contentGCD);
        for (int i = 1; i < bCfs.length; i++)
            contentGCD = ZippelGCD(bCfs[i], contentGCD);
        return contentGCD;
    }

    @SuppressWarnings("unchecked")
    private static lMultivariatePolynomial ZippelGCD(
            lMultivariatePolynomial a,
            lMultivariatePolynomial b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        lMultivariatePolynomial trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        //for bivariate polynomials always use dense interpolation
        if (variable == 1)
            return BrownGCD(a, b, rnd, variable, degreeBounds, evaluationStackLimit);

//        lMultivariatePolynomial content = ZippelContentGCD(a, b, variable);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        lMultivariatePolynomial factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0) {
            // univariate case
            lUnivariatePolynomialZp gcd = PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return asMultivariate(gcd, nVariables, variable, factory.ordering);
        }
        PrimitiveInput<
                lMonomialTerm,
                lMultivariatePolynomial,
                lUnivariatePolynomialZp>
                primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        lUnivariatePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            lMultivariatePolynomial poly = asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        lIntegersModulo domain = factory.domain;
        //store points that were already used in interpolation
        TLongHashSet globalEvaluationStack = new TLongHashSet();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        int failedSparseInterpolations = 0;

        int[] tmpDegreeBounds = degreeBounds.clone();
        main:
        while (true) {
            if (evaluationStackLimit == globalEvaluationStack.size())
                return null;

            long seedPoint = domain.randomElement(rnd);
            if (globalEvaluationStack.contains(seedPoint))
                continue;

            globalEvaluationStack.add(seedPoint);

            long lcVal = lcGCD.evaluate(seedPoint);
            if (lcVal == 0)
                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            lMultivariatePolynomial
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            lMultivariatePolynomial cVal = ZippelGCD(aVal, bVal, rnd, variable - 1, tmpDegreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > tmpDegreeBounds[variable - 1])
                //unlucky homomorphism
                continue;

            if (currExponent < tmpDegreeBounds[variable - 1]) {
                //better degree bound detected
                tmpDegreeBounds[variable - 1] = currExponent;
            }

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc() == lcVal;

            lSparseInterpolation sparseInterpolator = createInterpolation(variable, a, b, cVal, rnd);

            // we are applying dense interpolation for univariate skeleton coefficients
            lInterpolation denseInterpolation = new lInterpolation(variable, seedPoint, cVal);
            //previous interpolation (used to detect whether update doesn't change the result)
            lMultivariatePolynomial previousInterpolation;
            //local evaluation stack for points that are calculated via sparse interpolation (but not gcd evaluation) -> always same skeleton
            TLongHashSet localEvaluationStack = new TLongHashSet(globalEvaluationStack);
            while (true) {
                if (denseInterpolation.numberOfPoints() > tmpDegreeBounds[variable] + ALLOWED_OVER_INTERPOLATED_ATTEMPTS) {
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                long randomPoint = domain.randomElement(rnd);
                if (localEvaluationStack.contains(randomPoint))
                    continue;

                lcVal = lcGCD.evaluate(randomPoint);
                if (lcVal == 0)
                    continue;

                localEvaluationStack.add(randomPoint);

                cVal = sparseInterpolator.evaluate(randomPoint);
                if (cVal == null) {
                    ++failedSparseInterpolations;
                    if (failedSparseInterpolations == MAX_SPARSE_INTERPOLATION_FAILS)
                        throw new RuntimeException("Sparse interpolation failed");
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc() == lcVal;

                // Cache previous interpolation. NOTE: clone() is important, since the poly will
                // be modified inplace by the update() method
                previousInterpolation = denseInterpolation.getInterpolatingPolynomial().clone();
                denseInterpolation.update(randomPoint, cVal);

                // do division test
                if (tmpDegreeBounds[variable] <= denseInterpolation.numberOfPoints()
                        || previousInterpolation.equals(denseInterpolation.getInterpolatingPolynomial())) {
                    lMultivariatePolynomial result = doDivisionCheck(a, b, contentGCD, denseInterpolation, variable);
                    if (result != null)
                        return result;
                }
            }
        }
    }

    static lSparseInterpolation createInterpolation(int variable,
                                                    lMultivariatePolynomial a,
                                                    lMultivariatePolynomial b,
                                                    lMultivariatePolynomial skeleton,
                                                    RandomGenerator rnd) {
        boolean monic = a.coefficientOf(0, a.degree(0)).isConstant() && b.coefficientOf(0, a.degree(0)).isConstant();

        Set<DegreeVector> globalSkeleton = skeleton.getSkeleton();
        TIntObjectHashMap<lMultivariatePolynomial> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        lIntegersModulo domain = a.domain;

        int[] evaluationVariables = ArraysUtil.sequence(1, variable + 1);//variable inclusive
        long[] evaluationPoint = new long[evaluationVariables.length];

        lPrecomputedPowersHolder powers;
        search_for_good_evaluation_point:
        while (true) {
            //avoid zero evaluation points
            for (int i = variable - 2; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (evaluationPoint[i] == 0);

            //set temp value
            evaluationPoint[evaluationPoint.length - 1] = domain.randomElement(rnd);
            powers = new lPrecomputedPowersHolder(evaluationPoint, domain);
            int[] raiseFactors = ArraysUtil.arrayOf(1, evaluationVariables.length);

            for (lMultivariatePolynomial p : Arrays.asList(a, b, skeleton))
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables, raiseFactors).getSkeleton()))
                    continue search_for_good_evaluation_point;
            break;
        }

        int requiredNumberOfEvaluations = -1, monicScalingExponent = -1;
        for (TIntObjectIterator<lMultivariatePolynomial> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            lMultivariatePolynomial v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
            if (v.size() == 1)
                monicScalingExponent = it.key();
        }

        if (!ALWAYS_LINZIP) {
            if (monic)
                monicScalingExponent = -1;

            if (monic || monicScalingExponent != -1)
                return new lMonicInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                        evaluationVariables, evaluationPoint, powers, rnd, requiredNumberOfEvaluations, monicScalingExponent);
        }

        return new lLinZipInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, rnd);
    }

    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
    private static TIntObjectHashMap<lMultivariatePolynomial> getSkeleton(lMultivariatePolynomial poly) {
        TIntObjectHashMap<lMultivariatePolynomial> skeleton = new TIntObjectHashMap<>();
        for (lMonomialTerm term : poly.terms) {
            lMonomialTerm newDV = term.setZero(0);
            lMultivariatePolynomial coeff = skeleton.get(term.exponents[0]);
            if (coeff != null)
                coeff.add(newDV);
            else
                skeleton.put(term.exponents[0], lMultivariatePolynomial.create(poly.nVariables, poly.domain, poly.ordering, newDV));
        }
        return skeleton;
    }

    static abstract class lSparseInterpolation {
        /** the domain */
        final lIntegersModulo domain;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final lMultivariatePolynomial a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<lMultivariatePolynomial> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final long[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final lPrecomputedPowersHolder powers;
        /** random */
        final RandomGenerator rnd;

        lSparseInterpolation(lIntegersModulo domain, int variable,
                             lMultivariatePolynomial a, lMultivariatePolynomial b,
                             Set<DegreeVector> globalSkeleton,
                             TIntObjectHashMap<lMultivariatePolynomial> univarSkeleton,
                             int[] sparseUnivarDegrees, int[] evaluationVariables,
                             long[] evaluationPoint,
                             lPrecomputedPowersHolder powers, RandomGenerator rnd) {
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

        abstract lMultivariatePolynomial evaluate(long newPoint);
    }

    static final class lLinZipInterpolation extends lSparseInterpolation {
        lLinZipInterpolation(lIntegersModulo domain, int variable, lMultivariatePolynomial a,
                             lMultivariatePolynomial b, Set<DegreeVector> globalSkeleton,
                             TIntObjectHashMap<lMultivariatePolynomial> univarSkeleton, int[] sparseUnivarDegrees,
                             int[] evaluationVariables, long[] evaluationPoint, lPrecomputedPowersHolder powers,
                             RandomGenerator rnd) {
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
        lMultivariatePolynomial evaluate(long newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            @SuppressWarnings("unchecked")
            lLinZipSystem[] systems = new lLinZipSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new lLinZipSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


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
                    lUnivariatePolynomialZp
                            aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                    if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                        // unlucky main homomorphism or bad evaluation point
                        return null;

                    assert gcdUnivar.isMonic();

                    int totalEquations = 0;
                    for (lLinZipSystem system : systems) {
                        long rhs = gcdUnivar.degree() < system.univarDegree ? 0 : gcdUnivar.get(system.univarDegree);
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

                lMultivariatePolynomial result = a.createZero();
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

    private static SystemInfo solveLinZip(lMultivariatePolynomial factory, lLinZipSystem[] subSystems, int nUnknownScalings, lMultivariatePolynomial destination) {
        ArrayList<lMonomialTerm> unknowns = new ArrayList<>();
        for (lLinZipSystem system : subSystems)
            for (lMonomialTerm degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = nUnknownsMonomials + nUnknownScalings;
        ArrayList<long[]> lhsGlobal = new ArrayList<>();
        TLongArrayList rhsGlobal = new TLongArrayList();
        int offset = 0;
        lIntegersModulo domain = factory.domain;
        for (lLinZipSystem system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                long[] row = new long[nUnknownsTotal];
                long[] subRow = system.matrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                if (j > 0)
                    row[nUnknownsMonomials + j - 1] = system.scalingMatrix.get(j);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        long[] solution = new long[nUnknownsTotal];
        SystemInfo info = LinearAlgebra.solve(domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent) {
            @SuppressWarnings("unchecked")
            lMonomialTerm[] terms = new lMonomialTerm[unknowns.size()];
            for (int i = 0; i < terms.length; i++)
                terms[i] = unknowns.get(i).setCoefficient(solution[i]);
            destination.add(terms);
        }
        return info;
    }


    static final class lMonicInterpolation extends lSparseInterpolation {
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** univar exponent with monomial factor that can be used for scaling */
        final int monicScalingExponent;

        lMonicInterpolation(lIntegersModulo domain, int variable, lMultivariatePolynomial a, lMultivariatePolynomial b,
                            Set<DegreeVector> globalSkeleton, TIntObjectHashMap<lMultivariatePolynomial> univarSkeleton,
                            int[] sparseUnivarDegrees, int[] evaluationVariables, long[] evaluationPoint,
                            lPrecomputedPowersHolder powers, RandomGenerator rnd, int requiredNumberOfEvaluations,
                            int monicScalingExponent) {
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
        lMultivariatePolynomial evaluate(long newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationPoint.length - 1, newPoint);

            @SuppressWarnings("unchecked")
            lVandermondeSystem[] systems = new lVandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new lVandermondeSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                Arrays.fill(raiseFactors, 0, raiseFactors.length - 1, i + 1);
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;
                // evaluate a and b to univariate and calculate gcd
                lUnivariatePolynomialZp
                        aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        gcdUnivar = PolynomialGCD(aUnivar, bUnivar);

                if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                    // unlucky main homomorphism or bad evaluation point
                    return null;

                assert gcdUnivar.isMonic();

                if (monicScalingExponent != -1) {
                    // single scaling factor
                    // scale the system according to it

                    if (gcdUnivar.degree() < monicScalingExponent || gcdUnivar.get(monicScalingExponent) == 0)
                        // unlucky homomorphism
                        return null;

                    long normalization = evaluateExceptFirst(domain, powers, 1, univarSkeleton.get(monicScalingExponent).lt(), i + 1, variable - 1);
                    //normalize univariate gcd in order to reconstruct leading coefficient polynomial
                    normalization = domain.multiply(domain.reciprocal(gcdUnivar.get(monicScalingExponent)), normalization);
                    gcdUnivar = gcdUnivar.multiply(normalization);
                }

                boolean allDone = true;
                for (lVandermondeSystem system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        long rhs = gcdUnivar.degree() < system.univarDegree ? 0 : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (lVandermondeSystem system : systems) {
                //solve each system
                SystemInfo info = system.solve();
                if (info != SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            lMultivariatePolynomial gcdVal = a.createZero();
            for (lVandermondeSystem system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || system.solution[0] == 1;
                for (int i = 0; i < system.skeleton.length; i++) {
                    lMonomialTerm degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    long value = system.solution[i];
                    gcdVal.add(degreeVector.setCoefficient(value));
                }
            }

            return gcdVal;
        }
    }

    private static abstract class lLinearSystem {
        final int univarDegree;
        /** the domain */
        final lIntegersModulo domain;
        /** the skeleton */
        final lMonomialTerm[] skeleton;
        /** the lhs matrix */
        final ArrayList<long[]> matrix;
        /** the rhs values */
        final TLongArrayList rhs = new TLongArrayList();
        /** precomputed powers */
        final lPrecomputedPowersHolder powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        final int nVars;

        @SuppressWarnings("unchecked")
        lLinearSystem(int univarDegree, lMultivariatePolynomial skeleton, lPrecomputedPowersHolder powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            //todo refactor generics
            this.skeleton = skeleton.getSkeleton().toArray(new lMonomialTerm[skeleton.size()]);
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
    private static final class lLinZipSystem extends lLinearSystem {
        public lLinZipSystem(int univarDegree, lMultivariatePolynomial skeleton, lPrecomputedPowersHolder powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        private final TLongArrayList scalingMatrix = new TLongArrayList();

        public void oneMoreEquation(long rhsVal, boolean newScalingIntroduced) {
            long[] row = new long[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, 1, skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);

            if (newScalingIntroduced) {
                scalingMatrix.add(domain.negate(rhsVal));
                rhsVal = 0;
            } else
                scalingMatrix.add(0);
            rhs.add(rhsVal);
        }
    }

    /** Vandermonde system builder */
    private static final class lVandermondeSystem extends lLinearSystem {
        public lVandermondeSystem(int univarDegree, lMultivariatePolynomial skeleton, lPrecomputedPowersHolder powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        long[] solution = null;

        SystemInfo solve() {
            if (solution == null)
                solution = new long[nUnknownVariables()];

            if (nUnknownVariables() <= 8)
                // for small systems Gaussian elimination is indeed faster
                return LinearAlgebra.solve(domain, matrix.toArray(new long[matrix.size()][]), rhs.toArray(), solution);

            // solve vandermonde system
            long[] vandermondeRow = matrix.get(0);
            SystemInfo info = LinearAlgebra.solveVandermondeT(domain, vandermondeRow, rhs.toArray(), solution);
            if (info == SystemInfo.Consistent)
                for (int i = 0; i < solution.length; ++i)
                    solution[i] = domain.divide(solution[i], vandermondeRow[i]);

            return info;
        }

        public lVandermondeSystem oneMoreEquation(long rhsVal) {
            long[] row = new long[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, 1, skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);
            rhs.add(rhsVal);
            return this;
        }
    }

    private static long evaluateExceptFirst(lIntegersModulo domain,
                                            lPrecomputedPowersHolder powers,
                                            long coefficient,
                                            lMonomialTerm skeleton,
                                            int raiseFactor,
                                            int nVars) {
        long tmp = coefficient;
        for (int k = 0; k < nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[1 + k]));
        return tmp;
    }
}
