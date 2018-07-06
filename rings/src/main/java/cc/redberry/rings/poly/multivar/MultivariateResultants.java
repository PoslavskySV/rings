package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.ChineseRemainders.ChineseRemaindersMagic;
import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.linear.LinearSolver;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariateGCD;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.poly.univar.UnivariateResultants;
import cc.redberry.rings.primes.PrimesIterator;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TLongHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.stream.Stream;

import static cc.redberry.rings.ChineseRemainders.ChineseRemainders;
import static cc.redberry.rings.ChineseRemainders.createMagic;
import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.MachineArithmetic.*;
import static cc.redberry.rings.poly.Util.commonDenominator;
import static cc.redberry.rings.poly.multivar.Conversions64bit.*;
import static cc.redberry.rings.poly.multivar.MultivariateGCD.*;

/**
 * Polynomial resultants.
 */
public final class MultivariateResultants {
    private MultivariateResultants() {}

    /**
     * Computes discriminant of polynomial
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial> Poly Discriminant(Poly a, int variable) {
        Poly disc = (Poly) MultivariateDivision.divideExact(Resultant(a, a.derivative(variable), variable), a.lc(variable));
        return ((a.degree(variable) * (a.degree(variable) - 1) / 2) % 2 == 1) ? ((Poly) disc.negate()) : disc;
    }

    /**
     * Calculates polynomial resultant of two given polynomials with respect to specified variable
     *
     * @param a the first poly
     * @param b the second poly
     * @return the resultant
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial> Poly Resultant(Poly a, Poly b, int variable) {
        a.assertSameCoefficientRingWith(b);
        if (a.isOverFiniteField())
            return ResultantInGF(a, b, variable);
        if (a.isOverZ())
            return (Poly) ResultantInZ((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        if (Util.isOverRationals(a))
            return (Poly) ResultantInQ((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        if (a.isOverField())
            return (Poly) ZippelResultant((MultivariatePolynomial<BigInteger>) a, (MultivariatePolynomial<BigInteger>) b, variable);
        if (Util.isOverSimpleNumberField(a))
            return (Poly) ModularResultantInNumberField((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        if (Util.isOverRingOfIntegersOfSimpleNumberField(a))
            return (Poly) ModularResultantInRingOfIntegersOfNumberField((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        return tryNested(a, b, variable);
    }

    @SuppressWarnings("unchecked")
    private static <Poly extends AMultivariatePolynomial> Poly tryNested(Poly a, Poly b, int variable) {
        if (isOverUnivariate(a))
            return (Poly) ResultantOverUnivariate((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        else if (isOverUnivariateZp64(a))
            return (Poly) ResultantOverUnivariateZp64((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        else if (isOverMultivariate(a))
            return (Poly) ResultantOverMultivariate((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        else if (isOverMultivariateZp64(a))
            return (Poly) ResultantOverMultivariateZp64((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
        return (Poly) ClassicalResultant(a, b, variable);
    }

    private static <E> MultivariatePolynomial<UnivariatePolynomial<E>>
    ResultantOverUnivariate(MultivariatePolynomial<UnivariatePolynomial<E>> a,
                            MultivariatePolynomial<UnivariatePolynomial<E>> b, int variable) {
        return Resultant(
                MultivariatePolynomial.asNormalMultivariate(a, 0),
                MultivariatePolynomial.asNormalMultivariate(b, 0),
                variable)
                .asOverUnivariateEliminate(0);
    }

    private static <E> MultivariatePolynomial<MultivariatePolynomial<E>>
    ResultantOverMultivariate(MultivariatePolynomial<MultivariatePolynomial<E>> a,
                              MultivariatePolynomial<MultivariatePolynomial<E>> b, int variable) {
        int[] cfVars = ArraysUtil.sequence(a.lc().nVariables);
        int[] mainVars = ArraysUtil.sequence(a.lc().nVariables, a.lc().nVariables + a.nVariables);
        return Resultant(
                MultivariatePolynomial.asNormalMultivariate(a, cfVars, mainVars),
                MultivariatePolynomial.asNormalMultivariate(b, cfVars, mainVars),
                variable)
                .asOverMultivariateEliminate(cfVars);
    }

    private static MultivariatePolynomial<UnivariatePolynomialZp64>
    ResultantOverUnivariateZp64(MultivariatePolynomial<UnivariatePolynomialZp64> a,
                                MultivariatePolynomial<UnivariatePolynomialZp64> b,
                                int variable) {
        return Resultant(
                MultivariatePolynomialZp64.asNormalMultivariate(a, 0),
                MultivariatePolynomialZp64.asNormalMultivariate(b, 0),
                variable).asOverUnivariateEliminate(0);
    }

    private static MultivariatePolynomial<MultivariatePolynomialZp64>
    ResultantOverMultivariateZp64(MultivariatePolynomial<MultivariatePolynomialZp64> a,
                                  MultivariatePolynomial<MultivariatePolynomialZp64> b,
                                  int variable) {
        int[] cfVars = ArraysUtil.sequence(a.lc().nVariables);
        int[] mainVars = ArraysUtil.sequence(a.lc().nVariables, a.lc().nVariables + a.nVariables);
        return Resultant(
                MultivariatePolynomialZp64.asNormalMultivariate(a, cfVars, mainVars),
                MultivariatePolynomialZp64.asNormalMultivariate(b, cfVars, mainVars),
                variable).asOverMultivariateEliminate(cfVars);
    }

    static <E> MultivariatePolynomial<Rational<E>> ResultantInQ(
            MultivariatePolynomial<Rational<E>> a,
            MultivariatePolynomial<Rational<E>> b,
            int variable) {
        Util.Tuple2<MultivariatePolynomial<E>, E> aRat = Util.toCommonDenominator(a);
        Util.Tuple2<MultivariatePolynomial<E>, E> bRat = Util.toCommonDenominator(b);

        Ring<E> ring = aRat._1.ring;
        E correction = ring.multiply(
                ring.pow(aRat._2, b.degree(variable)),
                ring.pow(bRat._2, a.degree(variable)));
        return Util.asOverRationals(a.ring, Resultant(aRat._1, bRat._1, variable)).divideExact(new Rational<>(ring, correction));
    }

    /**
     * Computes polynomial resultant of two polynomials over finite field
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial>
    Poly ResultantInGF(Poly a, Poly b, int variable) {
        return (Poly) ZippelResultant(a, b, variable);
    }

    /**
     * Computes polynomial resultant of two polynomials over Z
     */
    public static MultivariatePolynomial<BigInteger> ResultantInZ(MultivariatePolynomial<BigInteger> a,
                                                                  MultivariatePolynomial<BigInteger> b,
                                                                  int variable) {
        return ModularResultantInZ(a, b, variable);
    }

    /**
     * Computes resultant via subresultant sequences
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly ClassicalResultant(Poly a, Poly b, int variable) {
        if (canConvertToZp64(a))
            return convertFromZp64(ClassicalResultant(asOverZp64(a), asOverZp64(b), variable));
        return UnivariateResultants.Resultant(a.asUnivariateEliminate(variable), b.asUnivariateEliminate(variable)).insertVariable(variable);
    }


    /* ============================================== Auxiliary methods ============================================= */

    /** structure with required input for resultant algorithms */
    static final class ResultantInput<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** input polynomials (with variables renamed and main variable dropped) */
        final Poly aReduced0, bReduced0;
        /** input polynomials (with variables renamed and main variable dropped) */
        final UnivariatePolynomial<Poly> aReduced, bReduced;
        /** earlyResultant if possible (in trivial cases) */
        final Poly earlyResultant;
        /** resultant degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** ring cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int evaluationStackLimit;
        /** Resultant part coming from GCD of monomial content of a and b **/
        final Poly monomialResultant;
        /**
         * degree of irreducible univariate polynomial used to construct field extension q^n (if the coefficient ring
         * has so small cardinality so that modular algorithm will fail)
         */
        final int finiteExtensionDegree;

        ResultantInput(Poly earlyResultant) {
            this.earlyResultant = earlyResultant;
            aReduced0 = bReduced0 = null;
            aReduced = bReduced = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialResultant = null;
            finiteExtensionDegree = -1;
        }

        ResultantInput(Poly aReduced0, Poly bReduced0,
                       Poly monomialResultant,
                       int evaluationStackLimit,
                       int[] degreeBounds, int[] mapping, int lastPresentVariable,
                       int finiteExtensionDegree) {
            //assert monomialGCD == null || aReduced.ring.isOne(monomialGCD.coefficient);
            this.aReduced0 = aReduced0;
            this.bReduced0 = bReduced0;
            this.aReduced = aReduced0.asUnivariateEliminate(0);
            this.bReduced = bReduced0.asUnivariateEliminate(0);
            this.monomialResultant = monomialResultant;
            this.earlyResultant = null;
            this.evaluationStackLimit = evaluationStackLimit;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
            this.finiteExtensionDegree = finiteExtensionDegree;
        }

        /** recover initial order of variables in the result */
        Poly restoreResultant(Poly result) {
            return AMultivariatePolynomial.renameVariables(result.insertVariable(0), mapping).multiply(monomialResultant);
        }
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    ResultantInput<Term, Poly> preparedResultantInput(Poly a, Poly b, int variable) {
        Poly trivialResultant = trivialResultant(a, b, variable);
        if (trivialResultant != null)
            return new ResultantInput<>(trivialResultant);

        BigInteger ringSize = a.coefficientRingCardinality();
        // ring cardinality, i.e. number of possible random choices
        int evaluationStackLimit = ringSize == null ? -1 : (ringSize.isInt() ? ringSize.intValue() : -1);

        // find monomial GCD
        // and remove monomial content from a and b
        a = a.clone();
        b = b.clone(); // prevent rewriting original data
        Term
                aContent = a.monomialContent(),
                bContent = b.monomialContent();
        a = a.divideOrNull(aContent);
        b = b.divideOrNull(bContent);

        if (aContent.exponents[variable] != 0 && bContent.exponents[variable] != 0)
            return new ResultantInput<>(a.createZero());

        Poly monomialResultant =
                resultantWithMonomial(aContent, b.create(bContent), variable)
                        .multiply(resultantWithMonomial(aContent, b, variable))
                        .multiply(resultantWithMonomial(a, bContent, variable));

        trivialResultant = trivialResultant(a, b, variable);
        if (trivialResultant != null)
            return new ResultantInput<>(trivialResultant.multiply(monomialResultant));

        int
                nVariables = a.nVariables,
                aDegrees[] = a.degreesRef(),
                bDegrees[] = b.degreesRef(),
                degreeBounds[] = new int[nVariables]; // degree bounds for resultant

        degreeBounds[variable] = Integer.MAX_VALUE;
        // populate initial resultant degree bounds based on Sylvester matrix
        int nUnused = 0;
        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                // skip main variable
                continue;
            // avoid potential int overflow
            degreeBounds[i] = safeToInt(safeAdd(
                    safeMultiply(aDegrees[i], bDegrees[variable]),
                    safeMultiply(bDegrees[i], aDegrees[variable])));
            if (degreeBounds[i] == 0)
                ++nUnused;
        }

        if (nUnused == (nVariables - 1)) {
            // all variables are unused => univariate resultant
            Poly t = trivialResultant(a, b, variable);
            assert t != null;
            return new ResultantInput<>(t.multiply(monomialResultant));
        }

        // adjust degree bounds with randomized substitutions and univariate images
        adjustDegreeBounds(a, b, variable, degreeBounds);

        // now swap variables so that the first variable will have the maximal degree (univariate resultant is fast),
        // and all non-used variables are at the end of poly's

        int[] variables = ArraysUtil.sequence(nVariables);
        //sort in descending order
        ArraysUtil.quickSort(ArraysUtil.negate(degreeBounds), variables);
        ArraysUtil.negate(degreeBounds);//recover degreeBounds
        degreeBounds = Arrays.copyOfRange(degreeBounds, 1, degreeBounds.length);

        int lastResVariable = 0; //recalculate lastPresentVariable
        for (; lastResVariable < degreeBounds.length; ++lastResVariable)
            if (degreeBounds[lastResVariable] == 0)
                break;
        --lastResVariable;

        // resultant variable is always the first
        a = AMultivariatePolynomial.renameVariables(a, variables);
        b = AMultivariatePolynomial.renameVariables(b, variables);

        // check whether coefficient ring cardinality is large enough
        int finiteExtensionDegree = 1;
        int cardinalityBound = safeToInt(safeMultiply(9, ArraysUtil.max(degreeBounds)));
        if (ringSize != null && ringSize.isInt() && ringSize.intValueExact() < cardinalityBound) {
            long ds = ringSize.intValueExact();
            finiteExtensionDegree = 2;
            long tmp = ds;
            for (; tmp < cardinalityBound; ++finiteExtensionDegree)
                tmp = tmp * ds;
        }

        return new ResultantInput<>(a, b,
                monomialResultant, evaluationStackLimit, degreeBounds, variables, lastResVariable, finiteExtensionDegree);
    }

    @SuppressWarnings("unchecked")
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void adjustDegreeBounds(Poly a, Poly b, int variable, int[] degreeBounds) {
        if (!a.isOverFiniteField())
            return;
        if (a.coefficientRingCharacteristic().bitLength() < 16)
            // don't do estimates for domains of small characteristic
            return;

        if (a instanceof MultivariatePolynomialZp64)
            adjustDegreeBounds((MultivariatePolynomialZp64) a, (MultivariatePolynomialZp64) b, variable, degreeBounds);
        else
            adjustDegreeBounds((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable, degreeBounds);
    }

    private static void adjustDegreeBounds(MultivariatePolynomialZp64 a,
                                           MultivariatePolynomialZp64 b,
                                           int variable,
                                           int[] degreeBounds) {
        int nVariables = a.nVariables;

        long[] subs = new long[nVariables];
        RandomGenerator rnd = PrivateRandom.getRandom();
        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            subs[i] = a.ring.randomNonZeroElement(rnd);
        }

        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            // all vars but i-th and variable
            int[] vars = ArraysUtil.remove(ArraysUtil.sequence(0, a.nVariables), new int[]{i, variable});
            long[] vals = ArraysUtil.remove(subs, new int[]{i, variable});
            MultivariatePolynomialZp64
                    mua = a.evaluate(vars, vals),
                    mub = b.evaluate(vars, vals);

            if (!mua.getSkeleton().equals(a.getSkeleton(variable, i)))
                continue;
            if (!mub.getSkeleton().equals(b.getSkeleton(variable, i)))
                continue;

            UnivariatePolynomial<UnivariatePolynomialZp64>
                    ua = mua.asOverUnivariateEliminate(i).asUnivariate(),
                    ub = mub.asOverUnivariateEliminate(i).asUnivariate();

            degreeBounds[i] = Math.min(degreeBounds[i], UnivariateResultants.Resultant(ua, ub).degree());
        }
    }

    private static <E> void adjustDegreeBounds(MultivariatePolynomial<E> a,
                                               MultivariatePolynomial<E> b,
                                               int variable,
                                               int[] degreeBounds) {
        int nVariables = a.nVariables;

        E[] subs = a.ring.createArray(nVariables);
        RandomGenerator rnd = PrivateRandom.getRandom();
        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            subs[i] = a.ring.randomNonZeroElement(rnd);
        }

        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            // all vars but i-th and variable
            int[] vars = ArraysUtil.remove(ArraysUtil.sequence(0, a.nVariables), new int[]{i, variable});
            E[] vals = ArraysUtil.remove(subs, new int[]{i, variable});
            MultivariatePolynomial<E>
                    mua = a.evaluate(vars, vals),
                    mub = b.evaluate(vars, vals);

            if (!mua.getSkeleton().equals(a.getSkeleton(variable, i)))
                continue;
            if (!mub.getSkeleton().equals(b.getSkeleton(variable, i)))
                continue;

            UnivariatePolynomial<UnivariatePolynomial<E>>
                    ua = mua.asOverUnivariateEliminate(i).asUnivariate(),
                    ub = mub.asOverUnivariateEliminate(i).asUnivariate();

            degreeBounds[i] = Math.min(degreeBounds[i], UnivariateResultants.Resultant(ua, ub).degree());
        }
    }

    @SuppressWarnings("unchecked")
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly trivialResultant(Poly a, Poly b, int variable) {
        if (a == b || a.isZero() || b.isZero() || a.equals(b))
            return a.createZero();
        if (a.degree(variable) == 0)
            return PolynomialMethods.polyPow(a, b.degree(variable));
        if (b.degree(variable) == 0)
            return PolynomialMethods.polyPow(b, a.degree(variable));
        if (a.size() == 1)
            return resultantWithMonomial(a.lt(), b, variable);
        if (b.size() == 1)
            return resultantWithMonomial(a, b.lt(), variable);
        if (a.univariateVariable() == variable && b.univariateVariable() == variable)
            return (Poly) AMultivariatePolynomial.asMultivariate(
                    UnivariateResultants.ResultantAsPoly(a.asUnivariate(), b.asUnivariate()),
                    a.nVariables, variable, a.ordering);
        return null;
    }

    /** resultant with monomial */
    @SuppressWarnings("unchecked")
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly resultantWithMonomial(Term monomial, Poly poly, int variable) {
        int varExponent = monomial.exponents[variable];
        Poly cFactor = PolynomialMethods.polyPow(poly.create(monomial.set(variable, 0)), poly.degree(variable));
        Poly xFactor = PolynomialMethods.polyPow(poly.evaluateAtZero(variable), varExponent);
        return cFactor.multiply(xFactor);
    }

    /** resultant with monomial */
    @SuppressWarnings("unchecked")
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly resultantWithMonomial(Poly poly, Term monomial, int variable) {
        Poly r = resultantWithMonomial(monomial, poly, variable);
        if (poly.degree(variable) % 2 == 1 && monomial.exponents[variable] % 2 == 1)
            r.negate();
        return r;
    }

    /**
     * Bivariate resultant always reduced to univariate case
     */
    static MultivariatePolynomialZp64 bivariateResultantZp64(UnivariatePolynomial<MultivariatePolynomialZp64> a,
                                                             UnivariatePolynomial<MultivariatePolynomialZp64> b) {
        MultivariatePolynomialZp64 factory = a.lc();
        IntegersZp64 ring = factory.ring;
        UnivariateRing<UnivariatePolynomialZp64> uRing = UnivariateRingZp64(ring);
        UnivariatePolynomial<UnivariatePolynomialZp64>
                aUni = a.mapCoefficients(uRing, MultivariatePolynomialZp64::asUnivariate),
                bUni = b.mapCoefficients(uRing, MultivariatePolynomialZp64::asUnivariate);
        return UnivariateResultants.Resultant(aUni, bUni).asMultivariate(factory.ordering).setNVariables(factory.nVariables);
    }

    /**
     * Bivariate resultant always reduced to univariate case
     */
    static <E> MultivariatePolynomial<E> bivariateResultantE(UnivariatePolynomial<MultivariatePolynomial<E>> a,
                                                             UnivariatePolynomial<MultivariatePolynomial<E>> b) {
        MultivariatePolynomial<E> factory = a.lc();
        Ring<E> ring = factory.ring;
        UnivariateRing<UnivariatePolynomial<E>> uRing = UnivariateRing(ring);
        UnivariatePolynomial<UnivariatePolynomial<E>>
                aUni = a.mapCoefficients(uRing, MultivariatePolynomial::asUnivariate),
                bUni = b.mapCoefficients(uRing, MultivariatePolynomial::asUnivariate);
        return UnivariateResultants.Resultant(aUni, bUni).asMultivariate(factory.ordering).setNVariables(factory.nVariables);
    }

    /**
     * Bivariate resultant always reduced to univariate case
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly bivariateResultant(UnivariatePolynomial<Poly> a,
                            UnivariatePolynomial<Poly> b) {
        if (a.lc() instanceof MultivariatePolynomialZp64)
            return (Poly) bivariateResultantZp64((UnivariatePolynomial) a, (UnivariatePolynomial) b);
        else
            return (Poly) bivariateResultantE((UnivariatePolynomial) a, (UnivariatePolynomial) b);
    }

    /* =========================================== In small characteristic ========================================== */

    /**
     * Resultant in small characteristic
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly ResultantInSmallCharacteristic(Poly a, Poly b, int variable) {
        // use "naive" algorithm for now.
        return ClassicalResultant(a, b, variable);
    }

    /**
     * Resultant in small characteristic
     */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly ResultantInSmallCharacteristic(UnivariatePolynomial<Poly> a, UnivariatePolynomial<Poly> b) {
        // use "naive" algorithm for now.
        return UnivariateResultants.Resultant(a, b);
    }


    /* ============================================== Modular resultant ============================================= */

    /**
     * Modular algorithm with Zippel sparse interpolation for resultant over Z
     */
    public static MultivariatePolynomial<BigInteger> ModularResultantInZ(MultivariatePolynomial<BigInteger> a,
                                                                         MultivariatePolynomial<BigInteger> b,
                                                                         int variable) {

        ResultantInput<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        return resInput.restoreResultant(ModularResultantInZ0(resInput.aReduced0, resInput.bReduced0));
    }

    /**
     * Number of updates in modular algorithm which doesn't change the result so we consider the obtained answer as the
     * correct one
     */
    private static final int N_UNCHANGED_INTERPOLATIONS = 5;

    private static BigInteger centralMultinomialCoefficient(int n, int d) {
        int q = n / d, r = n % d;
        return Z.factorial(n)
                .divideExact(Z.factorial(q).pow(d - r))
                .divideExact(Z.factorial(q + 1).pow(r));
    }

    /**
     * Modular algorithm with Zippel sparse interpolation for resultant over Z
     */
    static MultivariatePolynomial<BigInteger> ModularResultantInZ0(MultivariatePolynomial<BigInteger> a,
                                                                   MultivariatePolynomial<BigInteger> b) {

        // coefficient bound (not quite optimistic:)
        BigInteger
                aMax = a.maxAbsCoefficient(),
                bMax = b.maxAbsCoefficient();
        BigInteger bound2 = Z.getOne()
                .multiply(aMax.pow(b.degree()))                                // a coefficients
                .multiply(centralMultinomialCoefficient(a.size(), b.degree())) // a multiplication
                .multiply(bMax.pow(a.degree()))                                // b coefficients
                .multiply(centralMultinomialCoefficient(b.size(), a.degree())) // b multiplication
                .multiply(Z.factorial(a.degree() + b.degree()))          // overall determinant (Leibniz formula)
                .shiftLeft(1);                                                 // symmetric Zp form

        // choose better prime for start
        long startingPrime;
        if (Math.max(aMax.bitLength(), bMax.bitLength()) < 128)
            startingPrime = 1L << 30;
        else
            startingPrime = 1L << 60;

        PrimesIterator primesLoop = new PrimesIterator(startingPrime - (1 << 12));
        RandomGenerator random = PrivateRandom.getRandom();
        main_loop:
        while (true) {
            // prepare the skeleton
            long basePrime = primesLoop.take();
            assert basePrime != -1 : "long overflow";

            IntegersZp64 ring = Zp64(basePrime);
            // reduce Z -> Zp
            MultivariatePolynomialZp64
                    aMod = a.mapCoefficients(ring, c -> c.mod(basePrime).longValue()),
                    bMod = b.mapCoefficients(ring, c -> c.mod(basePrime).longValue());
            if (!aMod.sameSkeletonQ(a) || !bMod.sameSkeletonQ(b))
                continue;

            // the base image
            // accumulator to update coefficients via Chineese remainding
            MultivariatePolynomialZp64 base = ResultantInGF(aMod, bMod, 0).dropVariable(0);
            MultivariatePolynomialZp64 skeleton = base;
            MultivariatePolynomial<BigInteger> bBase = base.toBigPoly();

            // cache the previous base
            MultivariatePolynomial<BigInteger> previousBase = null;
            // number of times interpolation did not change the result
            int nUnchangedInterpolations = 0;

            BigInteger bBasePrime = Z.valueOf(basePrime);
            // over all primes
            while (true) {
                long prime = primesLoop.take();
                BigInteger bPrime = Z.valueOf(prime);
                ring = Zp64(prime);

                // reduce Z -> Zp
                aMod = a.mapCoefficients(ring, c -> c.mod(prime).longValue());
                bMod = b.mapCoefficients(ring, c -> c.mod(prime).longValue());
                if (!aMod.sameSkeletonQ(a) || !bMod.sameSkeletonQ(b))
                    continue;

                // calculate new resultant using previously calculated skeleton via sparse interpolation
                MultivariatePolynomialZp64 modularResultant = interpolateResultant(aMod, bMod, skeleton, random);
                if (modularResultant == null) {
                    // interpolation failed => assumed form is wrong => start over
                    continue main_loop;
                }

                // something wrong, start over
                if (!modularResultant.sameSkeletonQ(bBase)) {
                    base = modularResultant;
                    skeleton = base;
                    bBase = modularResultant.toBigPoly();
                    bBasePrime = bPrime;
                    previousBase = null;
                    nUnchangedInterpolations = 0;
                    continue;
                }

                //lifting
                BigInteger newBasePrime = bBasePrime.multiply(bPrime);
                PairedIterator<
                        Monomial<BigInteger>, MultivariatePolynomial<BigInteger>,
                        MonomialZp64, MultivariatePolynomialZp64
                        > iterator = new PairedIterator<>(bBase, modularResultant);
                ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, bBasePrime, bPrime);
                while (iterator.hasNext()) {
                    iterator.advance();

                    Monomial<BigInteger> baseTerm = iterator.aTerm;
                    MonomialZp64 imageTerm = iterator.bTerm;

                    if (baseTerm.coefficient.isZero())
                        // term is absent in the base
                        continue;

                    if (imageTerm.coefficient == 0) {
                        // term is absent in the modularResultant => remove it from the base
                        // bBase.subtract(baseTerm);
                        iterator.aIterator.remove();
                        continue;
                    }

                    long oth = imageTerm.coefficient;

                    // update base term
                    BigInteger newCoeff = ChineseRemainders(Z, magic, baseTerm.coefficient, BigInteger.valueOf(oth));
                    bBase.put(baseTerm.setCoefficient(newCoeff));
                }

                bBase = bBase.setRingUnsafe(new IntegersZp(newBasePrime));
                bBasePrime = newBasePrime;

                // two trials didn't change the result, probably we are done
                MultivariatePolynomial<BigInteger> candidate = MultivariatePolynomial.asPolyZSymmetric(bBase);
                if (previousBase != null && candidate.equals(previousBase))
                    ++nUnchangedInterpolations;
                else
                    nUnchangedInterpolations = 0;

                if (nUnchangedInterpolations >= N_UNCHANGED_INTERPOLATIONS || bBasePrime.compareTo(bound2) > 0)
                    return candidate;

                previousBase = candidate;
            }
        }
    }

    static MultivariatePolynomialZp64 interpolateResultant(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b, MultivariatePolynomialZp64 skeleton, RandomGenerator rnd) {
        a.assertSameCoefficientRingWith(b);
        skeleton = skeleton.setRingUnsafe(a.ring);
        if (a.nVariables == 2) {
            return bivariateResultant(
                    a.asUnivariateEliminate(0),
                    b.asUnivariateEliminate(0));
        }
        SparseInterpolationZp64 interpolation = createInterpolation(-1,
                a.asUnivariateEliminate(0),
                b.asUnivariateEliminate(0),
                skeleton, 1, rnd);
        if (interpolation == null)
            return null;
        MultivariatePolynomialZp64 res = interpolation.evaluate();
        if (res == null)
            return null;

        return res;
    }

    /* ================================== Modular resultant in simple number fields ================================ */

    /** if all coefficients are simple numbers */
    private static <E> MultivariatePolynomial<UnivariatePolynomial<E>>
    TrivialResultantInExtension(MultivariatePolynomial<UnivariatePolynomial<E>> a,
                                MultivariatePolynomial<UnivariatePolynomial<E>> b,
                                int variable) {
        AlgebraicNumberField<UnivariatePolynomial<E>> ring
                = (AlgebraicNumberField<UnivariatePolynomial<E>>) a.ring;

        if (!a.stream().allMatch(ring::isInTheBaseField)
                || !b.stream().allMatch(ring::isInTheBaseField))
            return null;

        Ring<E> cfRing = ring.getMinimalPolynomial().ring;
        MultivariatePolynomial<E>
                ar = a.mapCoefficients(cfRing, UnivariatePolynomial::cc),
                br = b.mapCoefficients(cfRing, UnivariatePolynomial::cc);
        return Resultant(ar, br, variable)
                .mapCoefficients(ring, cf -> UnivariatePolynomial.constant(cfRing, cf));
    }

    /**
     * Modular resultant in simple number field
     */
    public static MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> ModularResultantInNumberField(
            MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a,
            MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b,
            int variable) {

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> res = TrivialResultantInExtension(a, b, variable);
        if (res != null)
            return res;

        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> numberField = (AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>>) a.ring;
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly = numberField.getMinimalPolynomial();

        assert numberField.isField();

        a = a.clone();
        b = b.clone();

        // reduce problem to the case with integer monic minimal polynomial
        if (minimalPoly.stream().allMatch(Rational::isIntegral)) {
            // minimal poly is already monic & integer

            UnivariatePolynomial<BigInteger> minimalPolyZ = minimalPoly.mapCoefficients(Z, Rational::numerator);
            AlgebraicNumberField<UnivariatePolynomial<BigInteger>> numberFieldZ = new AlgebraicNumberField<>(minimalPolyZ);

            BigInteger
                    aDen = removeDenominators(a),
                    bDen = removeDenominators(b),
                    den = aDen.pow(b.degree(variable)).multiply(bDen.pow(a.degree(variable)));

            assert a.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral));
            assert b.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral));

            return ModularResultantInRingOfIntegersOfNumberField(
                    a.mapCoefficients(numberFieldZ, cf -> cf.mapCoefficients(Z, Rational::numerator)),
                    b.mapCoefficients(numberFieldZ, cf -> cf.mapCoefficients(Z, Rational::numerator)),
                    variable)
                    .mapCoefficients(numberField, cf -> cf.mapCoefficients(Q, r -> Q.mk(r, den)));
        } else {
            // replace s -> s / lc(minPoly)
            BigInteger minPolyLeadCoeff = commonDenominator(minimalPoly);
            Rational<BigInteger>
                    scale = new Rational<>(Z, Z.getOne(), minPolyLeadCoeff),
                    scaleReciprocal = scale.reciprocal();

            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>>
                    scaledNumberField = new AlgebraicNumberField<>(minimalPoly.scale(scale).monic());
            return ModularResultantInNumberField(
                    a.mapCoefficients(scaledNumberField, cf -> cf.scale(scale)),
                    b.mapCoefficients(scaledNumberField, cf -> cf.scale(scale)),
                    variable)
                    .mapCoefficients(numberField, cf -> cf.scale(scaleReciprocal));
        }
    }

    static BigInteger removeDenominators(MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a) {
        BigInteger denominator = Z.lcm(() -> a.stream().map(Util::commonDenominator).iterator());
        a.multiply(a.ring.valueOfBigInteger(denominator));
        return denominator;
    }

    /**
     * Modular algorithm with Zippel sparse interpolation for resultant over rings of integers
     */
    public static MultivariatePolynomial<UnivariatePolynomial<BigInteger>>
    ModularResultantInRingOfIntegersOfNumberField(MultivariatePolynomial<UnivariatePolynomial<BigInteger>> a,
                                                  MultivariatePolynomial<UnivariatePolynomial<BigInteger>> b,
                                                  int variable) {

        ResultantInput<Monomial<UnivariatePolynomial<BigInteger>>, MultivariatePolynomial<UnivariatePolynomial<BigInteger>>> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        return resInput.restoreResultant(ModularResultantInRingOfIntegersOfNumberField0(resInput.aReduced0, resInput.bReduced0));
    }

    /**
     * Modular algorithm with Zippel sparse interpolation for resultant over rings of integers
     */
    private static MultivariatePolynomial<UnivariatePolynomial<BigInteger>>
    ModularResultantInRingOfIntegersOfNumberField0(MultivariatePolynomial<UnivariatePolynomial<BigInteger>> a,
                                                   MultivariatePolynomial<UnivariatePolynomial<BigInteger>> b) {

        MultivariatePolynomial<UnivariatePolynomial<BigInteger>> r = TrivialResultantInExtension(a, b, 0);
        if (r != null)
            return r;

        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> numberField = (AlgebraicNumberField<UnivariatePolynomial<BigInteger>>) a.ring;
        UnivariatePolynomial<BigInteger> minimalPoly = numberField.getMinimalPolynomial();

        BigInteger
                aMax = a.stream().map(UnivariatePolynomial::maxAbsCoefficient).max(Rings.Z).orElse(BigInteger.ONE),
                bMax = b.stream().map(UnivariatePolynomial::maxAbsCoefficient).max(Rings.Z).orElse(BigInteger.ONE),
                mMax = minimalPoly.maxAbsCoefficient();

//        // bound on the value of resultant coefficients
//        assert a.degree() == a.degree(0) && b.degree() == b.degree(0);
//        BigInteger bound2 = Z.getOne()
//                .multiply(UnivariateResultants.polyPowNumFieldCfBound(aMax, mMax, minimalPoly.degree(), b.degree()))  // a coefficients
//                .multiply(centralMultinomialCoefficient(a.size(), b.degree())) // a multiplication
//                .multiply(UnivariateResultants.polyPowNumFieldCfBound(bMax, mMax, minimalPoly.degree(), a.degree()))  // b coefficients
//                .multiply(centralMultinomialCoefficient(b.size(), a.degree())) // b multiplication
//                .multiply(Z.factorial(a.degree() + b.degree()))          // overall determinant (Leibniz formula)
//                .shiftLeft(1);                                                // symmetric Zp form

        // choose better prime for start
        long startingPrime;
        if (Math.max(aMax.bitLength(), bMax.bitLength()) < 128)
            startingPrime = 1L << 30;
        else
            startingPrime = 1L << 60;

        UnivariateRing<UnivariatePolynomial<BigInteger>> auxRing = UnivariateRing(Z);

        PrimesIterator primesLoop = new PrimesIterator(startingPrime - (1 << 12));
        RandomGenerator random = PrivateRandom.getRandom();
        main_loop:
        while (true) {
            // prepare the skeleton
            long basePrime = primesLoop.take();
            assert basePrime != -1 : "long overflow";

            final IntegersZp64 baseRing = Zp64(basePrime);

            UnivariatePolynomialZp64 minimalPolyMod = UnivariatePolynomial.asOverZp64(minimalPoly, baseRing);
            FiniteField<UnivariatePolynomialZp64> numberFieldMod = new FiniteField<>(minimalPolyMod);

            // reduce Z -> Zp
            MultivariatePolynomial<UnivariatePolynomialZp64>
                    aMod = a.mapCoefficients(numberFieldMod, c -> UnivariatePolynomial.asOverZp64(c, baseRing)),
                    bMod = b.mapCoefficients(numberFieldMod, c -> UnivariatePolynomial.asOverZp64(c, baseRing));
            if (!aMod.sameSkeletonQ(a) || !bMod.sameSkeletonQ(b))
                continue;

            // the base image
            // accumulator to update coefficients via Chineese remainding
            MultivariatePolynomial<UnivariatePolynomialZp64> base;
            try {
                base = ResultantInGF(aMod, bMod, 0).dropVariable(0);
            } catch (Throwable t) {continue;}// bad base prime

            MultivariatePolynomial<UnivariatePolynomial<BigInteger>> bBase = base.mapCoefficients(auxRing, cf -> cf.asPolyZ(false).toBigPoly());

            // cache the previous base
            MultivariatePolynomial<UnivariatePolynomial<BigInteger>> previousBase = null;
            // number of times interpolation did not change the result
            int nUnchangedInterpolations = 0;

            BigInteger bBasePrime = Z.valueOf(basePrime);
            // over all primes
            while (true) {
                long prime = primesLoop.take();
                BigInteger bPrime = Z.valueOf(prime);
                IntegersZp64 ring = Zp64(prime);

                minimalPolyMod = UnivariatePolynomial.asOverZp64(minimalPoly, ring);
                numberFieldMod = new FiniteField<>(minimalPolyMod);

                // reduce Z -> Zp
                aMod = a.mapCoefficients(numberFieldMod, c -> UnivariatePolynomial.asOverZp64(c, ring));
                bMod = b.mapCoefficients(numberFieldMod, c -> UnivariatePolynomial.asOverZp64(c, ring));
                if (!aMod.sameSkeletonQ(a) || !bMod.sameSkeletonQ(b))
                    continue;

                // calculate new resultant using previously calculated skeleton via sparse interpolation
                MultivariatePolynomial<UnivariatePolynomialZp64> modularResultant;
                try {
                    modularResultant = interpolateResultant(aMod, bMod, base.mapCoefficients(numberFieldMod, cf -> cf.setModulusUnsafe(ring)), random);
                } catch (Throwable t) {
                    continue;
                }

                if (modularResultant == null) {
                    // interpolation failed => assumed form is wrong => start over
                    continue;
                }

                // something wrong, start over
                if (!modularResultant.sameSkeletonQ(bBase)) {
                    base = modularResultant;
                    bBase = modularResultant.mapCoefficients(auxRing, cf -> cf.asPolyZ(false).toBigPoly());
                    bBasePrime = bPrime;
                    previousBase = null;
                    nUnchangedInterpolations = 0;
                    continue;
                }

                //lifting
                BigInteger newBasePrime = bBasePrime.multiply(bPrime);
                PairedIterator<
                        Monomial<UnivariatePolynomial<BigInteger>>, MultivariatePolynomial<UnivariatePolynomial<BigInteger>>,
                        Monomial<UnivariatePolynomialZp64>, MultivariatePolynomial<UnivariatePolynomialZp64>
                        > iterator = new PairedIterator<>(bBase, modularResultant);
                ChineseRemaindersMagic<BigInteger> magic = createMagic(Z, bBasePrime, bPrime);
                while (iterator.hasNext()) {
                    iterator.advance();

                    Monomial<UnivariatePolynomial<BigInteger>> baseTerm = iterator.aTerm;
                    Monomial<UnivariatePolynomialZp64> imageTerm = iterator.bTerm;

                    if (baseTerm.coefficient.isZero())
                        // term is absent in the base
                        continue;

                    if (imageTerm.coefficient.isZero()) {
                        // term is absent in the modularResultant => remove it from the base
                        // bBase.subtract(baseTerm);
                        iterator.aIterator.remove();
                        continue;
                    }

                    UnivariateGCD.updateCRT(magic, baseTerm.coefficient, imageTerm.coefficient);
                }

                bBasePrime = newBasePrime;

                IntegersZp crtRing = Zp(bBasePrime);
                // two trials didn't change the result, probably we are done
                MultivariatePolynomial<UnivariatePolynomial<BigInteger>> candidate
                        = bBase.mapCoefficients(numberField, cf -> UnivariatePolynomial.asPolyZSymmetric(cf.setRingUnsafe(crtRing)));
                if (previousBase != null && candidate.equals(previousBase))
                    ++nUnchangedInterpolations;
                else
                    nUnchangedInterpolations = 0;

                if (nUnchangedInterpolations >= N_UNCHANGED_INTERPOLATIONS /* || bBasePrime.compareTo(bound2) > 0 */)
                    return candidate;

                previousBase = candidate;
            }
        }
    }

    static MultivariatePolynomial<UnivariatePolynomialZp64>
    interpolateResultant(MultivariatePolynomial<UnivariatePolynomialZp64> a,
                         MultivariatePolynomial<UnivariatePolynomialZp64> b,
                         MultivariatePolynomial<UnivariatePolynomialZp64> skeleton,
                         RandomGenerator rnd) {
        a.assertSameCoefficientRingWith(b);
        if (a.nVariables == 2) {
            return bivariateResultant(
                    a.asUnivariateEliminate(0),
                    b.asUnivariateEliminate(0));
        }

        skeleton = skeleton.setRingUnsafe(a.ring);
        SparseInterpolationE<UnivariatePolynomialZp64> interpolation = createInterpolation(-1,
                a.asUnivariateEliminate(0),
                b.asUnivariateEliminate(0),
                skeleton, 1, rnd);
        if (interpolation == null)
            return null;
        MultivariatePolynomial<UnivariatePolynomialZp64> res = interpolation.evaluate();
        if (res == null)
            return null;

        return res;
    }

    /* ============================================== Brown algorithm ============================================= */

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly BrownResultant(Poly a, Poly b, int variable) {
        if (a instanceof MultivariatePolynomialZp64)
            return (Poly) BrownResultant((MultivariatePolynomialZp64) a, (MultivariatePolynomialZp64) b, variable);
        else
            return (Poly) BrownResultant((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
    }

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    public static MultivariatePolynomialZp64 BrownResultant(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b, int variable) {
        ResultantInput<MonomialZp64, MultivariatePolynomialZp64> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        MultivariatePolynomialZp64 result = BrownResultantZp64(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        return resInput.restoreResultant(result);
    }

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    static MultivariatePolynomialZp64 BrownResultantZp64(UnivariatePolynomial<MultivariatePolynomialZp64> a,
                                                         UnivariatePolynomial<MultivariatePolynomialZp64> b,
                                                         int[] degreeBounds,
                                                         int variable) {
        if (variable == 0)
            return bivariateResultantZp64(a, b);

        MultivariateRing<MultivariatePolynomialZp64> mRing = (MultivariateRing<MultivariatePolynomialZp64>) a.ring;
        MultivariatePolynomialZp64 factory = mRing.factory();
        IntegersZp64 ring = factory.ring;

        //dense interpolation
        MultivariateInterpolation.InterpolationZp64 interpolation = null;
        //store points that were already used in interpolation
        TLongHashSet evaluationStack = new TLongHashSet();
        RandomGenerator rnd = PrivateRandom.getRandom();
        while (true) {
            if (evaluationStack.size() == ring.modulus)
                // all elements of the ring are tried
                return null;

            long v;
            do { v = ring.randomElement(rnd); }
            while (evaluationStack.contains(v));
            final long randomPoint = v;
            evaluationStack.add(randomPoint);

            MultivariateRing<MultivariatePolynomialZp64> imageRing = mRing.dropVariable();
            UnivariatePolynomial<MultivariatePolynomialZp64>
                    aMod = a.mapCoefficients(imageRing, cf -> cf.eliminate(variable, randomPoint)),
                    bMod = b.mapCoefficients(imageRing, cf -> cf.eliminate(variable, randomPoint));

            if (aMod.degree() != a.degree() || bMod.degree() != b.degree())
                continue;

            MultivariatePolynomialZp64 modResultant = BrownResultantZp64(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);
            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new MultivariateInterpolation.InterpolationZp64(variable, randomPoint, modResultant);
                continue;
            }

            // update interpolation
            interpolation.update(randomPoint, modResultant);

            if (interpolation.numberOfPoints() > degreeBounds[variable])
                return interpolation.getInterpolatingPolynomial();
        }
    }

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    public static <E> MultivariatePolynomial<E> BrownResultant(MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, int variable) {
        if (canConvertToZp64(a))
            return convertFromZp64(BrownResultant(asOverZp64(a), asOverZp64(b), variable));

        ResultantInput<Monomial<E>, MultivariatePolynomial<E>> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        MultivariatePolynomial<E> result = BrownResultantE(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        return resInput.restoreResultant(result);
    }

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    static <E> MultivariatePolynomial<E> BrownResultantE(UnivariatePolynomial<MultivariatePolynomial<E>> a,
                                                         UnivariatePolynomial<MultivariatePolynomial<E>> b,
                                                         int[] degreeBounds,
                                                         int variable) {
        if (variable == 0)
            return bivariateResultantE(a, b);

        MultivariateRing<MultivariatePolynomial<E>> mRing = (MultivariateRing<MultivariatePolynomial<E>>) a.ring;
        MultivariatePolynomial<E> factory = mRing.factory();
        Ring<E> ring = factory.ring;

        //dense interpolation
        MultivariateInterpolation.Interpolation<E> interpolation = null;
        //store points that were already used in interpolation
        Set<E> evaluationStack = new HashSet<>();
        RandomGenerator rnd = PrivateRandom.getRandom();
        while (true) {
            if (ring.cardinality().isInt() && evaluationStack.size() == ring.cardinality().intValue())
                // all elements of the ring are tried
                return null;

            E v;
            do { v = randomElement(ring, rnd); }
            while (evaluationStack.contains(v));
            final E randomPoint = v;
            evaluationStack.add(randomPoint);

            MultivariateRing<MultivariatePolynomial<E>> imageRing = mRing.dropVariable();
            UnivariatePolynomial<MultivariatePolynomial<E>>
                    aMod = a.mapCoefficients(imageRing, cf -> cf.eliminate(variable, randomPoint)),
                    bMod = b.mapCoefficients(imageRing, cf -> cf.eliminate(variable, randomPoint));

            if (aMod.degree() != a.degree() || bMod.degree() != b.degree())
                continue;

            MultivariatePolynomial<E> modResultant = BrownResultantE(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);
            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new MultivariateInterpolation.Interpolation<>(variable, randomPoint, modResultant);
                continue;
            }

            // update interpolation
            interpolation.update(randomPoint, modResultant);

            if (interpolation.numberOfPoints() > degreeBounds[variable])
                return interpolation.getInterpolatingPolynomial();
        }
    }

    /* ============================================== Zippel algorithm ============================================= */

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly ZippelResultant(Poly a, Poly b, int variable) {
        if (a instanceof MultivariatePolynomialZp64)
            return (Poly) ZippelResultant((MultivariatePolynomialZp64) a, (MultivariatePolynomialZp64) b, variable);
        else
            return (Poly) ZippelResultant((MultivariatePolynomial) a, (MultivariatePolynomial) b, variable);
    }

    /**
     * Zippel's algorithm for resultant with sparse interpolation
     */
    public static MultivariatePolynomialZp64 ZippelResultant(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b, int variable) {
        ResultantInput<MonomialZp64, MultivariatePolynomialZp64> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        MultivariatePolynomialZp64 result = ZippelResultantZp64(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        return resInput.restoreResultant(result);
    }

    /**
     * Zippel's algorithm for resultant with sparse interpolation
     */
    static MultivariatePolynomialZp64 ZippelResultantZp64(UnivariatePolynomial<MultivariatePolynomialZp64> a,
                                                          UnivariatePolynomial<MultivariatePolynomialZp64> b,
                                                          int[] degreeBounds,
                                                          int variable) {
        if (variable == 0)
            return bivariateResultantZp64(a, b);

        MultivariateRing<MultivariatePolynomialZp64> mRing = (MultivariateRing<MultivariatePolynomialZp64>) a.ring;
        MultivariatePolynomialZp64 factory = mRing.factory();
        IntegersZp64 ring = factory.ring;

        //dense interpolation
        MultivariateInterpolation.InterpolationZp64 denseInterpolation;
        //sparse interpolation
        SparseInterpolationZp64 sparseInterpolation;
        //store points that were already used in interpolation
        TLongHashSet globalEvaluationStack = new TLongHashSet();
        RandomGenerator rnd = PrivateRandom.getRandom();
        main:
        while (true) {
            if (globalEvaluationStack.size() == ring.modulus)
                // all elements of the ring are tried
                return null;

            long v;
            do { v = ring.randomElement(rnd); }
            while (globalEvaluationStack.contains(v));
            final long seedRandomPoint = v;
            globalEvaluationStack.add(seedRandomPoint);

            MultivariateRing<MultivariatePolynomialZp64> imageRing = mRing.dropVariable();
            UnivariatePolynomial<MultivariatePolynomialZp64>
                    aMod = a.mapCoefficients(imageRing, cf -> cf.eliminate(variable, seedRandomPoint)),
                    bMod = b.mapCoefficients(imageRing, cf -> cf.eliminate(variable, seedRandomPoint));

            if (aMod.degree() != a.degree() || bMod.degree() != b.degree())
                continue;

            // more checks
            for (int i = 0; i <= a.degree(); i++) {
                Set<DegreeVector> iniSkeleton = a.get(i).dropVariable(variable).getSkeleton();
                Set<DegreeVector> modSkeleton = aMod.get(i).getSkeleton();
                if (!iniSkeleton.equals(modSkeleton))
                    continue main;
            }
            for (int i = 0; i <= b.degree(); i++) {
                Set<DegreeVector> iniSkeleton = b.get(i).dropVariable(variable).getSkeleton();
                Set<DegreeVector> modSkeleton = bMod.get(i).getSkeleton();
                if (!iniSkeleton.equals(modSkeleton))
                    continue main;
            }

            // base evaluation
            MultivariatePolynomialZp64 baseResultant = ZippelResultantZp64(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);

            denseInterpolation = new MultivariateInterpolation.InterpolationZp64(variable, seedRandomPoint, baseResultant);
            sparseInterpolation = createInterpolation(variable, a, b, baseResultant, degreeBounds[variable], rnd);
            //local evaluation stack for points that are calculated via sparse interpolation (but not resultant evaluation) -> always same skeleton
            TLongHashSet localEvaluationStack = new TLongHashSet(globalEvaluationStack);
            while (true) {
                if (localEvaluationStack.size() == ring.modulus)
                    // all elements of the ring are tried
                    continue main;

                do { v = ring.randomElement(rnd); }
                while (localEvaluationStack.contains(v));
                final long randomPoint = v;
                localEvaluationStack.add(randomPoint);

                MultivariatePolynomialZp64 modResultant = sparseInterpolation.evaluate(randomPoint);
                if (modResultant == null)
                    continue;

                // update dense interpolation
                denseInterpolation.update(randomPoint, modResultant);

                if (denseInterpolation.numberOfPoints() > degreeBounds[variable])
                    return denseInterpolation.getInterpolatingPolynomial();
            }
        }
    }

    /**
     * Zippel's algorithm for resultant with sparse interpolation
     */
    public static <E> MultivariatePolynomial<E> ZippelResultant(MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, int variable) {
        if (canConvertToZp64(a))
            return convertFromZp64(ZippelResultant(asOverZp64(a), asOverZp64(b), variable));

        ResultantInput<Monomial<E>, MultivariatePolynomial<E>> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        MultivariatePolynomial<E> result = ZippelResultantE(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return resInput.restoreResultant(ResultantInSmallCharacteristic(resInput.aReduced, resInput.bReduced));

        return resInput.restoreResultant(result);
    }

    /**
     * Zippel's algorithm for resultant with sparse interpolation
     */
    static <E> MultivariatePolynomial<E> ZippelResultantE(UnivariatePolynomial<MultivariatePolynomial<E>> a,
                                                          UnivariatePolynomial<MultivariatePolynomial<E>> b,
                                                          int[] degreeBounds,
                                                          int variable) {
        if (variable == 0)
            return bivariateResultantE(a, b);

        MultivariateRing<MultivariatePolynomial<E>> mRing = (MultivariateRing<MultivariatePolynomial<E>>) a.ring;
        MultivariatePolynomial<E> factory = mRing.factory();
        Ring<E> ring = factory.ring;

        //dense interpolation
        MultivariateInterpolation.Interpolation<E> denseInterpolation;
        //sparse interpolation
        SparseInterpolationE<E> sparseInterpolation;
        //store points that were already used in interpolation
        Set<E> globalEvaluationStack = new HashSet<>();
        RandomGenerator rnd = PrivateRandom.getRandom();
        main:
        while (true) {
            if (ring.cardinality().isInt() && globalEvaluationStack.size() == ring.cardinality().intValueExact())
                // all elements of the ring are tried
                return null;

            E v;
            do { v = randomElement(ring, rnd); }
            while (globalEvaluationStack.contains(v));
            final E seedRandomPoint = v;
            globalEvaluationStack.add(seedRandomPoint);

            MultivariateRing<MultivariatePolynomial<E>> imageRing = mRing.dropVariable();
            UnivariatePolynomial<MultivariatePolynomial<E>>
                    aMod = a.mapCoefficients(imageRing, cf -> cf.eliminate(variable, seedRandomPoint)),
                    bMod = b.mapCoefficients(imageRing, cf -> cf.eliminate(variable, seedRandomPoint));

            if (aMod.degree() != a.degree() || bMod.degree() != b.degree())
                continue;

            // more checks
            for (int i = 0; i <= a.degree(); i++) {
                Set<DegreeVector> iniSkeleton = a.get(i).dropVariable(variable).getSkeleton();
                Set<DegreeVector> modSkeleton = aMod.get(i).getSkeleton();
                if (!iniSkeleton.equals(modSkeleton))
                    continue main;
            }
            for (int i = 0; i <= b.degree(); i++) {
                Set<DegreeVector> iniSkeleton = b.get(i).dropVariable(variable).getSkeleton();
                Set<DegreeVector> modSkeleton = bMod.get(i).getSkeleton();
                if (!iniSkeleton.equals(modSkeleton))
                    continue main;
            }

            // base evaluation
            MultivariatePolynomial<E> baseResultant = ZippelResultantE(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);

            denseInterpolation = new MultivariateInterpolation.Interpolation<>(variable, seedRandomPoint, baseResultant);
            sparseInterpolation = createInterpolation(variable, a, b, baseResultant, degreeBounds[variable], rnd);
            //local evaluation stack for points that are calculated via sparse interpolation (but not resultant evaluation) -> always same skeleton
            Set<E> localEvaluationStack = new HashSet<>(globalEvaluationStack);
            while (true) {
                if (ring.cardinality().isInt() && localEvaluationStack.size() == ring.cardinality().intValueExact())
                    // all elements of the ring are tried
                    continue main;

                do { v = randomElement(ring, rnd); }
                while (localEvaluationStack.contains(v));
                final E randomPoint = v;
                localEvaluationStack.add(randomPoint);

                MultivariatePolynomial<E> modResultant = sparseInterpolation.evaluate(randomPoint);
                if (modResultant == null)
                    continue;

                // update dense interpolation
                denseInterpolation.update(randomPoint, modResultant);

                if (denseInterpolation.numberOfPoints() > degreeBounds[variable])
                    return denseInterpolation.getInterpolatingPolynomial();
            }
        }
    }

    static SparseInterpolationZp64 createInterpolation(int variable,
                                                       UnivariatePolynomial<MultivariatePolynomialZp64> a,
                                                       UnivariatePolynomial<MultivariatePolynomialZp64> b,
                                                       MultivariatePolynomialZp64 skeleton,
                                                       int expectedNumberOfEvaluations,
                                                       RandomGenerator rnd) {
        MultivariatePolynomialZp64 factory = a.lc();
        assert factory.nVariables > 1;
        skeleton = skeleton.clone().setAllCoefficientsToUnit();

        Set<DegreeVector> globalSkeleton = skeleton.getSkeleton();
        TIntObjectHashMap<MultivariatePolynomialZp64> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        IntegersZp64 ring = factory.ring;

        int lastVariable = variable == -1 ? factory.nVariables - 1 : variable;
        int[] evaluationVariables = ArraysUtil.sequence(1, lastVariable + 1);//variable inclusive
        long[] evaluationPoint = new long[evaluationVariables.length];

        MultivariatePolynomialZp64.lPrecomputedPowersHolder powers;
        int fails = 0;
        search_for_good_evaluation_point:
        while (true) {
            if (fails >= MAX_FAILED_SUBSTITUTIONS)
                return null;
            //avoid zero evaluation points
            for (int i = lastVariable - 1; i >= 0; --i)
                do {
                    evaluationPoint[i] = ring.randomElement(rnd);
                } while (evaluationPoint[i] == 0);

            powers = mkPrecomputedPowers(a, b, evaluationVariables, evaluationPoint);

            Iterator<MultivariatePolynomialZp64> it = Stream.concat(Stream.concat(a.stream(), b.stream()), Stream.of(skeleton)).iterator();
            while (it.hasNext()) {
                MultivariatePolynomialZp64 p = it.next();
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables).getSkeleton())) {
                    ++fails;
                    continue search_for_good_evaluation_point;
                }
            }
            break;
        }

        int requiredNumberOfEvaluations = -1;
        for (TIntObjectIterator<MultivariatePolynomialZp64> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            MultivariatePolynomialZp64 v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
        }

        return new SparseInterpolationZp64(ring, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations, requiredNumberOfEvaluations, rnd);
    }

    static MultivariatePolynomialZp64.lPrecomputedPowersHolder mkPrecomputedPowers(
            UnivariatePolynomial<MultivariatePolynomialZp64> a, UnivariatePolynomial<MultivariatePolynomialZp64> b,
            int[] evaluationVariables, long[] evaluationPoint) {
        MultivariatePolynomialZp64 factory = a.lc();
        int[] degrees = null;
        for (int i = 0; i <= a.degree(); i++) {
            if (degrees == null)
                degrees = a.get(i).degreesRef();
            else
                degrees = ArraysUtil.max(degrees, a.get(i).degreesRef());
        }
        assert degrees != null;
        for (int i = 0; i <= b.degree(); i++)
            degrees = ArraysUtil.max(degrees, b.get(i).degreesRef());

        MultivariatePolynomialZp64.lPrecomputedPowers[] pp = new MultivariatePolynomialZp64.lPrecomputedPowers[factory.nVariables];
        for (int i = 0; i < evaluationVariables.length; ++i)
            pp[evaluationVariables[i]] = new MultivariatePolynomialZp64.lPrecomputedPowers(
                    Math.min(degrees[evaluationVariables[i]], MultivariatePolynomialZp64.MAX_POWERS_CACHE_SIZE),
                    evaluationPoint[i], factory.ring);
        return new MultivariatePolynomialZp64.lPrecomputedPowersHolder(factory.ring, pp);
    }

    static final class SparseInterpolationZp64 {
        /** the ring */
        final IntegersZp64 ring;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final UnivariatePolynomial<MultivariatePolynomialZp64> a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomialZp64> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /**
         * variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ...
         * variable] }
         */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables}; values {@code V} for variables {@code [1, 2
         * ... variable-1]} are fixed and successive powers {@code V, V^2, V^3 ...} will be used to form Vandermonde
         * matrix
         */
        final long[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final MultivariatePolynomialZp64.lPrecomputedPowersHolder powers;
        /** Efficient structures for evaluating polynomials for interpolation */
        final MultivariateGCD.ZippelEvaluationsZp64[] aEvals, bEvals;
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** random */
        final RandomGenerator rnd;
        /** random */
        final MultivariatePolynomialZp64 factory;

        SparseInterpolationZp64(IntegersZp64 ring,
                                int variable,
                                UnivariatePolynomial<MultivariatePolynomialZp64> a,
                                UnivariatePolynomial<MultivariatePolynomialZp64> b,
                                Set<DegreeVector> globalSkeleton,
                                TIntObjectHashMap<MultivariatePolynomialZp64> univarSkeleton,
                                int[] sparseUnivarDegrees,
                                int[] evaluationVariables,
                                long[] evaluationPoint,
                                MultivariatePolynomialZp64.lPrecomputedPowersHolder powers,
                                int expectedNumberOfEvaluations,
                                int requiredNumberOfEvaluations,
                                RandomGenerator rnd) {
            this.ring = ring;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.globalSkeleton = globalSkeleton;
            this.univarSkeleton = univarSkeleton;
            this.sparseUnivarDegrees = sparseUnivarDegrees;
            this.evaluationPoint = evaluationPoint;
            this.aEvals = new MultivariateGCD.ZippelEvaluationsZp64[a.degree() + 1];
            this.bEvals = new MultivariateGCD.ZippelEvaluationsZp64[b.degree() + 1];
            for (int i = 0; i < aEvals.length; ++i)
                aEvals[i] = MultivariateGCD.createEvaluations(a.get(i), evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations);
            for (int i = 0; i < bEvals.length; ++i)
                bEvals[i] = MultivariateGCD.createEvaluations(b.get(i), evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations);

            this.evaluationVariables = evaluationVariables;
            this.powers = powers;
            this.requiredNumberOfEvaluations = requiredNumberOfEvaluations;
            this.rnd = rnd;
            this.factory = a.lc();
        }

        /**
         * Returns interpolating resultant
         */
        public final MultivariatePolynomialZp64 evaluate() {
            return evaluate(evaluationPoint[evaluationPoint.length - 1]);
        }

        /**
         * Returns interpolating resultant at specified point
         *
         * @param newPoint evaluation point
         * @return resultant at {@code variable = newPoint}
         */
        public final MultivariatePolynomialZp64 evaluate(long newPoint) {
            // constant is constant
            if (globalSkeleton.size() == 1)
                return factory.create(((MonomialZp64) globalSkeleton.iterator().next()).setCoefficient(1));
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationVariables[evaluationVariables.length - 1], newPoint);
            return evaluate0(newPoint);
        }

        private MultivariatePolynomialZp64 evaluate0(long newPoint) {
            @SuppressWarnings("unchecked")
            MultivariateGCD.lVandermondeSystem[] systems = new MultivariateGCD.lVandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new MultivariateGCD.lVandermondeSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? factory.nVariables - 1 : variable - 1);

            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                int raiseFactor = i + 1;

                long lastVarValue = newPoint;
                if (variable == -1)
                    lastVarValue = ring.powMod(lastVarValue, raiseFactor);

                // evaluate a and b to bivariate and calculate resultant
                UnivariatePolynomial<UnivariatePolynomialZp64>
                        aBivar = UnivariatePolynomial.zero(UnivariateRingZp64(ring)),
                        bBivar = UnivariatePolynomial.zero(UnivariateRingZp64(ring));
                for (int j = 0; j < aEvals.length; ++j) {
                    aBivar.set(j, aEvals[j].evaluate(raiseFactor, lastVarValue));
                    if (aBivar.get(j).degree() != a.get(j).degree(0))
                        return null;
                }
                for (int j = 0; j < bEvals.length; ++j) {
                    bBivar.set(j, bEvals[j].evaluate(raiseFactor, lastVarValue));
                    if (bBivar.get(j).degree() != b.get(j).degree(0))
                        return null;
                }

                UnivariatePolynomialZp64 resUnivar = UnivariateResultants.Resultant(aBivar, bBivar);

                if (!univarSkeleton.keySet().containsAll(resUnivar.exponents()))
                    // univariate resultant contains terms that are not present in the skeleton
                    // again unlucky main homomorphism
                    return null;

                boolean allDone = true;
                for (MultivariateGCD.lVandermondeSystem system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        long rhs = resUnivar.degree() < system.univarDegree ? 0 : resUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        if (system.nEquations() < system.nUnknownVariables())
                            allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (MultivariateGCD.lVandermondeSystem system : systems) {
                //solve each system
                LinearSolver.SystemInfo info = system.solve();
                if (info != LinearSolver.SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomialZp64 resVal = factory.createZero();
            for (MultivariateGCD.lVandermondeSystem system : systems) {
                for (int i = 0; i < system.skeleton.length; i++) {
                    MonomialZp64 degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    long value = system.solution[i];
                    resVal.add(degreeVector.setCoefficient(value));
                }
            }

            return resVal;
        }
    }

    static <E> SparseInterpolationE<E> createInterpolation(int variable,
                                                           UnivariatePolynomial<MultivariatePolynomial<E>> a,
                                                           UnivariatePolynomial<MultivariatePolynomial<E>> b,
                                                           MultivariatePolynomial<E> skeleton,
                                                           int expectedNumberOfEvaluations,
                                                           RandomGenerator rnd) {
        MultivariatePolynomial<E> factory = a.lc();
        assert factory.nVariables > 1;
        skeleton = skeleton.clone().setAllCoefficientsToUnit();

        Set<DegreeVector> globalSkeleton = skeleton.getSkeleton();
        TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        Ring<E> ring = factory.ring;

        int lastVariable = variable == -1 ? factory.nVariables - 1 : variable;
        int[] evaluationVariables = ArraysUtil.sequence(1, lastVariable + 1);//variable inclusive
        E[] evaluationPoint = ring.createArray(evaluationVariables.length);

        MultivariatePolynomial.PrecomputedPowersHolder<E> powers;
        int fails = 0;
        search_for_good_evaluation_point:
        while (true) {
            if (fails >= MAX_FAILED_SUBSTITUTIONS)
                return null;
            //avoid zero evaluation points
            for (int i = lastVariable - 1; i >= 0; --i)
                do {
                    evaluationPoint[i] = randomElement(ring, rnd);
                } while (ring.isZero(evaluationPoint[i]));

            powers = mkPrecomputedPowers(a, b, evaluationVariables, evaluationPoint);

            Iterator<MultivariatePolynomial<E>> it = Stream.concat(Stream.concat(a.stream(), b.stream()), Stream.of(skeleton)).iterator();
            while (it.hasNext()) {
                MultivariatePolynomial<E> p = it.next();
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables).getSkeleton())) {
                    ++fails;
                    continue search_for_good_evaluation_point;
                }
            }
            break;
        }

        int requiredNumberOfEvaluations = -1;
        for (TIntObjectIterator<MultivariatePolynomial<E>> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            MultivariatePolynomial<E> v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
        }

        return new SparseInterpolationE<>(ring, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations, requiredNumberOfEvaluations, rnd);
    }

    static <E> MultivariatePolynomial.PrecomputedPowersHolder<E> mkPrecomputedPowers(
            UnivariatePolynomial<MultivariatePolynomial<E>> a, UnivariatePolynomial<MultivariatePolynomial<E>> b,
            int[] evaluationVariables, E[] evaluationPoint) {
        MultivariatePolynomial<E> factory = a.lc();
        int[] degrees = null;
        for (int i = 0; i <= a.degree(); i++) {
            if (degrees == null)
                degrees = a.get(i).degreesRef();
            else
                degrees = ArraysUtil.max(degrees, a.get(i).degreesRef());
        }
        assert degrees != null;
        for (int i = 0; i <= b.degree(); i++)
            degrees = ArraysUtil.max(degrees, b.get(i).degreesRef());

        @SuppressWarnings("unchecked")
        MultivariatePolynomial.PrecomputedPowers<E>[] pp = new MultivariatePolynomial.PrecomputedPowers[factory.nVariables];
        for (int i = 0; i < evaluationVariables.length; ++i)
            pp[evaluationVariables[i]] = new MultivariatePolynomial.PrecomputedPowers<>(
                    Math.min(degrees[evaluationVariables[i]], MultivariatePolynomialZp64.MAX_POWERS_CACHE_SIZE),
                    evaluationPoint[i], factory.ring);
        return new MultivariatePolynomial.PrecomputedPowersHolder<>(factory.ring, pp);
    }

    static final class SparseInterpolationE<E> {
        /** the ring */
        final Ring<E> ring;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final UnivariatePolynomial<MultivariatePolynomial<E>> a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /**
         * variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ...
         * variable] }
         */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables}; values {@code V} for variables {@code [1, 2
         * ... variable-1]} are fixed and successive powers {@code V, V^2, V^3 ...} will be used to form Vandermonde
         * matrix
         */
        final E[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final MultivariatePolynomial.PrecomputedPowersHolder<E> powers;
        /** Efficient structures for evaluating polynomials for interpolation */
        final MultivariateGCD.ZippelEvaluations<E>[] aEvals, bEvals;
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** random */
        final RandomGenerator rnd;
        /** random */
        final MultivariatePolynomial<E> factory;

        @SuppressWarnings("unchecked")
        SparseInterpolationE(Ring<E> ring,
                             int variable,
                             UnivariatePolynomial<MultivariatePolynomial<E>> a,
                             UnivariatePolynomial<MultivariatePolynomial<E>> b,
                             Set<DegreeVector> globalSkeleton,
                             TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton,
                             int[] sparseUnivarDegrees,
                             int[] evaluationVariables,
                             E[] evaluationPoint,
                             MultivariatePolynomial.PrecomputedPowersHolder<E> powers,
                             int expectedNumberOfEvaluations,
                             int requiredNumberOfEvaluations,
                             RandomGenerator rnd) {
            this.ring = ring;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.globalSkeleton = globalSkeleton;
            this.univarSkeleton = univarSkeleton;
            this.sparseUnivarDegrees = sparseUnivarDegrees;
            this.evaluationPoint = evaluationPoint;
            this.aEvals = new MultivariateGCD.ZippelEvaluations[a.degree() + 1];
            this.bEvals = new MultivariateGCD.ZippelEvaluations[b.degree() + 1];
            for (int i = 0; i < aEvals.length; ++i) {
                aEvals[i] = MultivariateGCD.createEvaluations(a.get(i), evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations);
            }
            for (int i = 0; i < bEvals.length; ++i)
                bEvals[i] = MultivariateGCD.createEvaluations(b.get(i), evaluationVariables, evaluationPoint, powers, expectedNumberOfEvaluations);

            this.evaluationVariables = evaluationVariables;
            this.powers = powers;
            this.requiredNumberOfEvaluations = requiredNumberOfEvaluations;
            this.rnd = rnd;
            this.factory = a.lc();
        }

        /**
         * Returns interpolating resultant
         */
        public final MultivariatePolynomial<E> evaluate() {
            return evaluate(evaluationPoint[evaluationPoint.length - 1]);
        }

        /**
         * Returns interpolating resultant at specified point
         *
         * @param newPoint evaluation point
         * @return resultant at {@code variable = newPoint}
         */
        public final MultivariatePolynomial<E> evaluate(E newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationVariables[evaluationVariables.length - 1], newPoint);
            return evaluate0(newPoint);
        }

        private MultivariatePolynomial<E> evaluate0(E newPoint) {
            @SuppressWarnings("unchecked")
            MultivariateGCD.VandermondeSystem<E>[] systems = new MultivariateGCD.VandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new MultivariateGCD.VandermondeSystem<>(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? factory.nVariables - 1 : variable - 1);

            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                int raiseFactor = i + 1;

                E lastVarValue = newPoint;
                if (variable == -1)
                    lastVarValue = ring.pow(lastVarValue, raiseFactor);

                // evaluate a and b to bivariate and calculate resultant
                UnivariatePolynomial<UnivariatePolynomial<E>>
                        aBivar = UnivariatePolynomial.zero(UnivariateRing(ring)),
                        bBivar = UnivariatePolynomial.zero(UnivariateRing(ring));
                for (int j = 0; j < aEvals.length; ++j) {
                    aBivar.set(j, aEvals[j].evaluate(raiseFactor, lastVarValue));
                    if (aBivar.get(j).degree() != a.get(j).degree(0))
                        return null;
                }
                for (int j = 0; j < bEvals.length; ++j) {
                    bBivar.set(j, bEvals[j].evaluate(raiseFactor, lastVarValue));
                    if (bBivar.get(j).degree() != b.get(j).degree(0))
                        return null;
                }

                UnivariatePolynomial<E> resUnivar = UnivariateResultants.Resultant(aBivar, bBivar);

                if (!univarSkeleton.keySet().containsAll(resUnivar.exponents()))
                    // univariate resultant contains terms that are not present in the skeleton
                    // again unlucky main homomorphism
                    return null;

                boolean allDone = true;
                for (MultivariateGCD.VandermondeSystem<E> system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        E rhs = resUnivar.degree() < system.univarDegree ? ring.getZero() : resUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        if (system.nEquations() < system.nUnknownVariables())
                            allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (MultivariateGCD.VandermondeSystem<E> system : systems) {
                //solve each system
                LinearSolver.SystemInfo info = system.solve();
                if (info != LinearSolver.SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomial<E> resVal = factory.createZero();
            for (MultivariateGCD.VandermondeSystem<E> system : systems) {
                for (int i = 0; i < system.skeleton.length; i++) {
                    Monomial<E> degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    E value = system.solution[i];
                    resVal.add(degreeVector.setCoefficient(value));
                }
            }

            return resVal;
        }
    }
}
