package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.linear.LinearSolver;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.poly.univar.UnivariateResultants;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TLongHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;
import java.util.stream.Stream;

import static cc.redberry.rings.Rings.UnivariateRing;
import static cc.redberry.rings.Rings.UnivariateRingZp64;
import static cc.redberry.rings.poly.MachineArithmetic.*;
import static cc.redberry.rings.poly.multivar.MultivariateGCD.*;

/**
 * Polynomial resultants.
 */
public final class MultivariateResultants {
    private MultivariateResultants() {}

    /**
     * A plain method to compute resultant which just delegates to the univariate case
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly PlainResultant(Poly a, Poly b, int variable) {
        return UnivariateResultants.Resultant(a.asUnivariateEliminate(variable), b.asUnivariateEliminate(variable)).insertVariable(variable);
    }

    /** structure with required input for resultant algorithms */
    static final class ResultantInput<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
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
            aReduced = bReduced = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialResultant = null;
            finiteExtensionDegree = -1;
        }

        ResultantInput(UnivariatePolynomial<Poly> aReduced,
                       UnivariatePolynomial<Poly> bReduced,
                       Poly monomialResultant,
                       int evaluationStackLimit,
                       int[] degreeBounds, int[] mapping, int lastPresentVariable,
                       int finiteExtensionDegree) {
            //assert monomialGCD == null || aReduced.ring.isOne(monomialGCD.coefficient);
            this.aReduced = aReduced;
            this.bReduced = bReduced;
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
                degreeBounds[] = new int[nVariables]; // degree bounds for gcd

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

        return new ResultantInput<>(a.asUnivariateEliminate(0), b.asUnivariateEliminate(0),
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
        Ring<E> ring = a.ring;
        E[] subs = ring.createZeroesArray(a.nVariables);
        RandomGenerator rnd = PrivateRandom.getRandom();
        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            subs[i] = ring.randomNonZeroElement(rnd);
        }

        for (int i = 0; i < nVariables; i++) {
            if (i == variable)
                continue;
            // all vars but i-th and variable
            int[] vars = ArraysUtil.remove(ArraysUtil.sequence(0, a.nVariables), new int[]{i, variable});
            E[] vals = ArraysUtil.remove(subs, new int[]{i, variable});
            MultivariatePolynomial<E>
                    mua = a.eliminate(vars, vals),
                    mub = b.eliminate(vars, vals);

            if (!mua.getSkeleton().equals(a.getSkeletonDrop(variable, i)))
                continue;
            if (!mub.getSkeleton().equals(b.getSkeletonDrop(variable, i)))
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
        return UnivariateResultants.Resultant(aUni, bUni).asMultivariate().setNVariables(factory.nVariables);
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
        return UnivariateResultants.Resultant(aUni, bUni).asMultivariate().setNVariables(factory.nVariables);
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

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    public static MultivariatePolynomialZp64 BrownResultant(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b, int variable) {
        ResultantInput<MonomialZp64, MultivariatePolynomialZp64> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return null;

        MultivariatePolynomialZp64 result = BrownResultant(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return null;

        return resInput.restoreResultant(result);
    }

    /**
     * Brown's algorithm for resultant with dense interpolation
     */
    static MultivariatePolynomialZp64 BrownResultant(UnivariatePolynomial<MultivariatePolynomialZp64> a,
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

            MultivariatePolynomialZp64 modResultant = BrownResultant(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);
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
     * Zippel's algorithm for resultant with sparse interpolation
     */
    public static MultivariatePolynomialZp64 ZippelResultant(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b, int variable) {
        ResultantInput<MonomialZp64, MultivariatePolynomialZp64> resInput = preparedResultantInput(a, b, variable);
        if (resInput.earlyResultant != null)
            return resInput.earlyResultant;

        if (resInput.finiteExtensionDegree > 1)
            return null;

        MultivariatePolynomialZp64 result = ZippelResultant(
                resInput.aReduced, resInput.bReduced,
                resInput.degreeBounds, resInput.lastPresentVariable);
        if (result == null)
            return null;

        return resInput.restoreResultant(result);
    }

    /**
     * Zippel's algorithm for resultant with sparse interpolation
     */
    static MultivariatePolynomialZp64 ZippelResultant(UnivariatePolynomial<MultivariatePolynomialZp64> a,
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
            MultivariatePolynomialZp64 baseResultant = ZippelResultant(aMod, bMod, degreeBounds, variable - 1).insertVariable(variable);

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

                UnivariatePolynomialZp64 gcdUnivar = UnivariateResultants.Resultant(aBivar, bBivar);

                if (!univarSkeleton.keySet().containsAll(gcdUnivar.exponents()))
                    // univariate gcd contains terms that are not present in the skeleton
                    // again unlucky main homomorphism
                    return null;

                boolean allDone = true;
                for (MultivariateGCD.lVandermondeSystem system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        long rhs = gcdUnivar.degree() < system.univarDegree ? 0 : gcdUnivar.get(system.univarDegree);
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
}
