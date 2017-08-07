package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.*;
import cc.r2.core.poly.multivar.HenselLifting.Evaluation;
import cc.r2.core.poly.multivar.HenselLifting.IEvaluation;
import cc.r2.core.poly.multivar.HenselLifting.lEvaluation;
import cc.r2.core.poly.univar.*;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TLongHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static cc.r2.core.poly.CommonPolynomialsArithmetics.polyPow;
import static cc.r2.core.poly.multivar.AMultivariatePolynomial.asMultivariate;
import static cc.r2.core.poly.multivar.AMultivariatePolynomial.renameVariables;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideExact;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateFactorization {
    private MultivariateFactorization() {}


    /* ============================================== Auxiliary methods ============================================= */

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorUnivariate(Poly poly) {
        int uVar = poly.univariateVariable();
        FactorDecomposition<? extends IUnivariatePolynomial>
                uFactors = Factorization.factor(poly.asUnivariate());
        return uFactors.map(u -> (Poly) asMultivariate(u, poly.nVariables, uVar, poly.ordering));
    }

    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorToPrimitive(Poly poly) {
        if (poly.isEffectiveUnivariate())
            return FactorDecomposition.singleFactor(poly);
        FactorDecomposition<Poly> result = FactorDecomposition.empty(poly);
        for (int i = 0; i < poly.nVariables; i++) {
            Poly factor = poly.asUnivariate(i).content();
            result.addFactor(factor, 1);
            poly = divideExact(poly, factor);
        }
        result.addFactor(poly, 1);
        return result;
    }

    private static int[] add(int[] array, int value) {
        int[] res = new int[array.length];
        for (int i = 0; i < array.length; i++)
            res[i] = array[i] = value;
        return res;
    }

    private interface FactorizationAlgorithm<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        FactorDecomposition<Poly> factor(Poly poly, boolean switchToExtensionField);
    }


    /** Number of univariate factorizations performed with different evaluation homomorphisms before doing Hensel lifting **/
    private static final long UNIVARIATE_FACTORIZATION_ATTEMPTS = 3;

    /** starting extension field exponent */
    private static final int EXTENSION_FIELD_EXPONENT = 3;

    /**
     * Factors primitive, square-free bivariate polynomial over Zp switching to extension field
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    private static FactorDecomposition<lMultivariatePolynomialZp> factorInExtensionField(
            lMultivariatePolynomialZp poly,
            FactorizationAlgorithm<MonomialTerm<lUnivariatePolynomialZp>, MultivariatePolynomial<lUnivariatePolynomialZp>> algorithm) {

        lIntegersModulo domain = poly.domain;
        int startingDegree = EXTENSION_FIELD_EXPONENT;
        while (true) {
            FiniteField<lUnivariatePolynomialZp> extensionField = new FiniteField<>(
                    IrreduciblePolynomials.randomIrreduciblePolynomial(
                            domain.modulus, startingDegree++, PrivateRandom.getRandom()));

            FactorDecomposition<lMultivariatePolynomialZp> result =
                    factorInExtensionField(poly, extensionField, algorithm);

            if (result != null)
                return result;
        }
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp switching to extension field
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    private static FactorDecomposition<lMultivariatePolynomialZp> factorInExtensionField(
            lMultivariatePolynomialZp poly,
            FiniteField<lUnivariatePolynomialZp> extensionField,
            FactorizationAlgorithm<MonomialTerm<lUnivariatePolynomialZp>, MultivariatePolynomial<lUnivariatePolynomialZp>> algorithm) {

        FactorDecomposition<MultivariatePolynomial<lUnivariatePolynomialZp>> factorization
                = algorithm.factor(poly.mapCoefficients(extensionField, extensionField::valueOf), false);
        if (factorization == null)
            // too small extension
            return null;
        if (!factorization.constantFactor.cc().isConstant())
            return null;

        FactorDecomposition<lMultivariatePolynomialZp> result = FactorDecomposition.constantFactor(poly.createConstant(factorization.constantFactor.cc().cc()));
        for (int i = 0; i < factorization.size(); i++) {
            if (!factorization.get(i).stream().allMatch(p -> p.isConstant()))
                return null;
            result.addFactor(factorization.get(i).mapCoefficients(poly.domain, lUnivariatePolynomialZp::cc), factorization.getExponent(i));
        }
        return result;
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp switching to extension field
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    private static <E>
    FactorDecomposition<MultivariatePolynomial<E>> factorInExtensionField(
            MultivariatePolynomial<E> poly,
            FactorizationAlgorithm<MonomialTerm<UnivariatePolynomial<E>>, MultivariatePolynomial<UnivariatePolynomial<E>>> algorithm) {
        Domain<E> domain = poly.domain;

        int startingDegree = EXTENSION_FIELD_EXPONENT;
        while (true) {
            FiniteField<UnivariatePolynomial<E>> extensionField = new FiniteField<>(
                    IrreduciblePolynomials.randomIrreduciblePolynomial(
                            domain, startingDegree++, PrivateRandom.getRandom()));

            FactorDecomposition<MultivariatePolynomial<E>> result =
                    factorInExtensionField(poly, extensionField, algorithm);

            if (result != null)
                return result;
        }
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp switching to extension field
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    private static <E>
    FactorDecomposition<MultivariatePolynomial<E>> factorInExtensionField(
            MultivariatePolynomial<E> poly,
            FiniteField<UnivariatePolynomial<E>> extensionField,
            FactorizationAlgorithm<MonomialTerm<UnivariatePolynomial<E>>, MultivariatePolynomial<UnivariatePolynomial<E>>> algorithm) {
        FactorDecomposition<MultivariatePolynomial<UnivariatePolynomial<E>>> factorization
                = algorithm.factor(poly.mapCoefficients(extensionField, c -> UnivariatePolynomial.constant(poly.domain, c)), false);
        if (factorization == null)
            // too small extension
            return null;
        if (!factorization.constantFactor.cc().isConstant())
            return null;

        FactorDecomposition<MultivariatePolynomial<E>> result = FactorDecomposition.constantFactor(poly.createConstant(factorization.constantFactor.cc().cc()));
        for (int i = 0; i < factorization.size(); i++) {
            if (!factorization.get(i).stream().allMatch(UnivariatePolynomial::isConstant))
                return null;
            result.addFactor(factorization.get(i).mapCoefficients(poly.domain, UnivariatePolynomial::cc), factorization.getExponent(i));
        }
        return result;
    }

    /* ================================= Bivariate factorization over finite fields ================================= */

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static FactorDecomposition<lMultivariatePolynomialZp>
    bivariateDenseFactorSquareFree(lMultivariatePolynomialZp poly) {
        return bivariateDenseFactorSquareFree(poly, true);
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly                   primitive, square-free bivariate polynomial over Zp
     * @param switchToExtensionField whether to switch to extension field if domain cardinality is too small
     * @return factor decomposition
     */
    static FactorDecomposition<lMultivariatePolynomialZp>
    bivariateDenseFactorSquareFree(lMultivariatePolynomialZp poly, boolean switchToExtensionField) {
        assert poly.nUsedVariables() <= 2 && IntStream.range(2, poly.nVariables).allMatch(i -> poly.degree(i) == 0) : poly;

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        lMultivariatePolynomialZp reducedPoly = poly;
        int[] degreeBounds = reducedPoly.degrees();

        // use main variable with maximal degree
        boolean swapVariables = false;
        if (degreeBounds[1] > degreeBounds[0]) {
            swapVariables = true;
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            ArraysUtil.swap(degreeBounds, 0, 1);
        }

        lMultivariatePolynomialZp xDerivative = reducedPoly.derivative(0);
        if (xDerivative.isZero()) {
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            swapVariables = !swapVariables;
            xDerivative = reducedPoly.derivative(0);
        }

        lMultivariatePolynomialZp dGCD = MultivariateGCD.PolynomialGCD(xDerivative, reducedPoly);
        if (!dGCD.isConstant()) {
            FactorDecomposition<lMultivariatePolynomialZp>
                    gcdFactorization = bivariateDenseFactorSquareFree(dGCD, switchToExtensionField),
                    restFactorization = bivariateDenseFactorSquareFree(divideExact(reducedPoly, dGCD), switchToExtensionField);

            if (gcdFactorization == null || restFactorization == null) {
                assert !switchToExtensionField;
                return null;
            }

            gcdFactorization.addAll(restFactorization);
            if (swapVariables)
                swap(gcdFactorization);

            return gcdFactorization;
        }

        lIntegersModulo domain = reducedPoly.domain;
        // degree in main variable
        int degree = reducedPoly.degree(0);
        // substitution value for second variable
        long ySubstitution = -1;
        // univariate factorization
        FactorDecomposition<lUnivariatePolynomialZp> uFactorization = null;

        // number of univariate factorizations tried
        int univariateFactorizations = 0;
        boolean tryZeroFirst = true;

        TLongHashSet evaluationStack = new TLongHashSet();
        RandomGenerator random = PrivateRandom.getRandom();
        while (univariateFactorizations < UNIVARIATE_FACTORIZATION_ATTEMPTS) {
            if (evaluationStack.size() == domain.modulus)
                if (switchToExtensionField)
                    // switch to extension field
                    return factorInExtensionField(poly, MultivariateFactorization::bivariateDenseFactorSquareFree);
                else
                    return null;

            long substitution;
            if (tryZeroFirst) {
                // first try to substitute 0 for second variable, then use random values
                substitution = 0;
                tryZeroFirst = false;
            } else
                do {
                    substitution = domain.randomElement(random);
                } while (evaluationStack.contains(substitution));
            evaluationStack.add(substitution);

            lMultivariatePolynomialZp image = reducedPoly.evaluate(1, substitution);
            if (image.degree() != degree)
                // unlucky substitution
                continue;

            if (image.cc() == 0)
                // c.c. must not be zero since input is primitive
                // => unlucky substitution
                continue;

            lUnivariatePolynomialZp uImage = image.asUnivariate();
            if (!SquareFreeFactorization.isSquareFree(uImage))
                // ensure that univariate image is also square free
                continue;

            FactorDecomposition<lUnivariatePolynomialZp> factorization = Factorization.factor(uImage);
            if (factorization.size() == 1)
                // irreducible polynomial
                return FactorDecomposition.singleFactor(poly);


            if (uFactorization == null || factorization.size() < uFactorization.size()) {
                // better univariate factorization found
                uFactorization = factorization;
                ySubstitution = substitution;
            }

            //if (ySubstitution == 0)
            //   break;

            ++univariateFactorizations;
        }

        assert ySubstitution != -1;
        assert uFactorization.factors.stream().allMatch(lUnivariatePolynomialZp::isMonic);

        // univariate factors are calculated
        @SuppressWarnings("unchecked")
        List<lUnivariatePolynomialZp> factorList = uFactorization.factors;

        // we don't precompute correct leading coefficients of bivariate factors
        // instead, we add the l.c. of the product to a list of lifting factors
        // in order to obtain correct factorization with monic factors mod (y - y0)^l
        // and then perform l.c. correction at the recombination stage

        long[] evals = new long[poly.nVariables - 1];
        evals[0] = ySubstitution;
        lEvaluation evaluation = new lEvaluation(poly.nVariables, evals, domain, reducedPoly.ordering);
        lMultivariatePolynomialZp lc = reducedPoly.lc(0);
        if (!lc.isConstant()) {
            // add lc to lifting factors
            lUnivariatePolynomialZp ulc = evaluation.evaluateFrom(lc, 1).asUnivariate();
            assert ulc.isConstant();
            factorList.add(0, ulc);
        } else
            factorList.get(0).multiply(lc.cc());

        // final factors to lift
        lUnivariatePolynomialZp[] factors = factorList.toArray(new lUnivariatePolynomialZp[factorList.size()]);

        // lift univariate factorization
        int liftDegree = reducedPoly.degree(1) + 1;

        Domain<lUnivariatePolynomialZp> uDomain = new UnivariatePolynomials<>(factors[0]);
        // series expansion around y = y0 for initial poly
        UnivariatePolynomial<lUnivariatePolynomialZp> baseSeries =
                HenselLifting.seriesExpansionDense(uDomain, reducedPoly, 1, evaluation);

        // lifted factors (each factor represented as series around y = y0)
        UnivariatePolynomial<lUnivariatePolynomialZp>[] lifted =
                HenselLifting.bivariateLiftDense(baseSeries, factors, liftDegree);

        if (!lc.isConstant())
            // drop auxiliary l.c. from factors
            lifted = Arrays.copyOfRange(lifted, 1, lifted.length);

        // factors are lifted => do recombination
        FactorDecomposition<lMultivariatePolynomialZp> result = denseBivariateRecombination(reducedPoly, baseSeries, lifted, evaluation, liftDegree);

        if (swapVariables)
            // reconstruct original variables order
            swap(result);

        return result;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void swap(FactorDecomposition<Poly> factorDecomposition) {
        for (int i = 0; i < factorDecomposition.factors.size(); i++)
            factorDecomposition.factors.set(i, AMultivariatePolynomial.swapVariables(factorDecomposition.get(i), 0, 1));
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static <E> FactorDecomposition<MultivariatePolynomial<E>>
    bivariateDenseFactorSquareFree(MultivariatePolynomial<E> poly) {
        return bivariateDenseFactorSquareFree(poly, true);
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly                   primitive, square-free bivariate polynomial over Zp
     * @param switchToExtensionField whether to switch to extension field if domain cardinality is too small
     * @return factor decomposition
     */
    static <E> FactorDecomposition<MultivariatePolynomial<E>>
    bivariateDenseFactorSquareFree(MultivariatePolynomial<E> poly, boolean switchToExtensionField) {
        assert poly.nUsedVariables() <= 2 && IntStream.range(2, poly.nVariables).allMatch(i -> poly.degree(i) == 0);

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        MultivariatePolynomial<E> reducedPoly = poly;
        int[] degreeBounds = reducedPoly.degrees();

        // use main variable with maximal degree
        boolean swapVariables = false;
        if (degreeBounds[1] > degreeBounds[0]) {
            swapVariables = true;
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            ArraysUtil.swap(degreeBounds, 0, 1);
        }

        MultivariatePolynomial<E> xDerivative = reducedPoly.derivative(0);
        if (xDerivative.isZero()) {
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            swapVariables = !swapVariables;
            xDerivative = reducedPoly.derivative(0);
        }

        MultivariatePolynomial<E> dGCD = MultivariateGCD.PolynomialGCD(xDerivative, reducedPoly);
        if (!dGCD.isConstant()) {
            FactorDecomposition<MultivariatePolynomial<E>>
                    gcdFactorization = bivariateDenseFactorSquareFree(dGCD, switchToExtensionField),
                    restFactorization = bivariateDenseFactorSquareFree(divideExact(reducedPoly, dGCD), switchToExtensionField);

            if (gcdFactorization == null || restFactorization == null) {
                assert !switchToExtensionField;
                return null;
            }

            gcdFactorization.addAll(restFactorization);
            if (swapVariables)
                swap(gcdFactorization);

            return gcdFactorization;
        }

        Domain<E> domain = reducedPoly.domain;
        // degree in main variable
        int degree = reducedPoly.degree(0);
        // substitution value for second variable
        E ySubstitution = null;
        // univariate factorization
        FactorDecomposition<UnivariatePolynomial<E>> uFactorization = null;

        // number of univariate factorizations tried
        int univariateFactorizations = 0;
        boolean tryZeroFirst = true;
        HashSet<E> evaluationStack = new HashSet<>();
        while (univariateFactorizations < UNIVARIATE_FACTORIZATION_ATTEMPTS) {
            if (domain.cardinality().isInt() && domain.cardinality().intValueExact() == evaluationStack.size())
                if (switchToExtensionField)
                    // switch to extension field
                    return factorInExtensionField(poly, MultivariateFactorization::bivariateDenseFactorSquareFree);
                else
                    return null;

            E substitution;
            if (tryZeroFirst) {
                // first try to substitute 0 for second variable, then use random values
                substitution = domain.getZero();
                tryZeroFirst = false;
            } else
                do {
                    substitution = domain.randomElement(PrivateRandom.getRandom());
                } while (evaluationStack.contains(substitution));
            evaluationStack.add(substitution);

            MultivariatePolynomial<E> image = reducedPoly.evaluate(1, substitution);
            if (image.degree() != degree)
                // unlucky substitution
                continue;

            if (domain.isZero(image.cc()))
                // c.c. must not be zero since input is primitive
                // => unlucky substitution
                continue;

            UnivariatePolynomial<E> uImage = image.asUnivariate();
            if (!SquareFreeFactorization.isSquareFree(uImage))
                // ensure that univariate image is also square free
                continue;

            FactorDecomposition<UnivariatePolynomial<E>> factorization = Factorization.factor(uImage);
            if (factorization.size() == 1)
                // irreducible polynomial
                return FactorDecomposition.singleFactor(poly);


            if (uFactorization == null || factorization.size() < uFactorization.size()) {
                // better univariate factorization found
                uFactorization = factorization;
                ySubstitution = substitution;
            }

            //if (ySubstitution == 0)
            //   break;

            ++univariateFactorizations;
        }

        assert ySubstitution != null;

        // univariate factors are calculated
        @SuppressWarnings("unchecked")
        List<UnivariatePolynomial<E>> factorList = uFactorization.factors;

        // we don't precompute correct leading coefficients of bivariate factors
        // instead, we add the l.c. of the product to a list of lifting factors
        // in order to obtain correct factorization with monic factors mod (y - y0)^l
        // and then perform l.c. correction at the recombination stage

        E[] evals = domain.createZeroesArray(poly.nVariables - 1);
        evals[0] = ySubstitution;
        Evaluation<E> evaluation = new Evaluation<>(poly.nVariables, evals, domain, reducedPoly.ordering);
        MultivariatePolynomial<E> lc = reducedPoly.lc(0);
        if (!lc.isConstant()) {
            // add lc to lifting factors
            UnivariatePolynomial<E> ulc = evaluation.evaluateFrom(lc, 1).asUnivariate();
            assert ulc.isConstant();
            factorList.add(0, ulc);
        } else
            factorList.get(0).multiply(lc.cc());

        // final factors to lift
        @SuppressWarnings("unchecked")
        UnivariatePolynomial<E>[] factors = factorList.toArray(new UnivariatePolynomial[factorList.size()]);

        // lift univariate factorization
        int liftDegree = reducedPoly.degree(1) + 1;

        Domain<UnivariatePolynomial<E>> uDomain = new UnivariatePolynomials<>(factors[0]);
        // series expansion around y = y0 for initial poly
        UnivariatePolynomial<UnivariatePolynomial<E>> baseSeries =
                HenselLifting.seriesExpansionDense(uDomain, reducedPoly, 1, evaluation);

        // lifted factors (each factor represented as series around y = y0)
        UnivariatePolynomial<UnivariatePolynomial<E>>[] lifted =
                HenselLifting.bivariateLiftDense(baseSeries, factors, liftDegree);

        if (!lc.isConstant())
            // drop auxiliary l.c. from factors
            lifted = Arrays.copyOfRange(lifted, 1, factors.length);

        // factors are lifted => do recombination
        FactorDecomposition<MultivariatePolynomial<E>> result = denseBivariateRecombination(reducedPoly, baseSeries, lifted, evaluation, liftDegree);

        if (swapVariables)
            // reconstruct original variables order
            for (int i = 0; i < result.factors.size(); i++)
                result.factors.set(i, AMultivariatePolynomial.swapVariables(result.get(i), 0, 1));

        return result;
    }

    /** cache of references **/
    private static int[][] naturalSequenceRefCache = new int[32][];

    private static int[] createSeq(int n) {
        int[] r = new int[n];
        for (int i = 0; i < n; i++)
            r[i] = i;
        return r;
    }

    /** returns sequence of natural numbers */
    private static int[] naturalSequenceRef(int n) {
        if (n >= naturalSequenceRefCache.length)
            return createSeq(n);
        if (naturalSequenceRefCache[n] != null)
            return naturalSequenceRefCache[n];
        return naturalSequenceRefCache[n] = createSeq(n);
    }

    /** select elements by their positions */
    private static int[] select(int[] data, int[] positions) {
        int[] r = new int[positions.length];
        int i = 0;
        for (int p : positions)
            r[i++] = data[p];
        return r;
    }

    /**
     * Naive dense recombination for bivariate factors
     *
     * @param factory        multivariate polynomial factory
     * @param poly           series around y = y0 for base polynomial (which we attempt to factor)
     * @param modularFactors univariate factors mod (y-y0)^liftDegree
     * @param evaluation     evaluation ideal
     * @param liftDegree     lifting degree (ideal power)
     * @return true factorization
     */
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    FactorDecomposition<Poly> denseBivariateRecombination(
            Poly factory,
            UnivariatePolynomial<uPoly> poly,
            UnivariatePolynomial<uPoly>[] modularFactors,
            IEvaluation<Term, Poly> evaluation,
            int liftDegree) {

        int[] modIndexes = naturalSequenceRef(modularFactors.length);
        FactorDecomposition<Poly> trueFactors = FactorDecomposition.empty(factory);
        UnivariatePolynomial<uPoly> fRest = poly;
        int s = 1;

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                UnivariatePolynomial<uPoly> mFactor = lcInSeries(fRest);
                for (int i : indexes)
                    // todo:
                    // implement IUnivariatePolynomial#multiplyLow(int)
                    // and replace truncate(int) with multiplyLow(int)
                    mFactor = mFactor.multiply(modularFactors[i]).truncate(liftDegree - 1);

                // get primitive part in first variable (remove R[y] content)
                UnivariatePolynomial<uPoly> factor =
                        changeDenseRepresentation(
                                changeDenseRepresentation(mFactor).primitivePart());

                UnivariatePolynomial<uPoly>[] qd = DivisionWithRemainder.divideAndRemainder(fRest, factor, true);
                if (qd != null && qd[1].isZero()) {
                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                    trueFactors.addFactor(HenselLifting.denseSeriesToPoly(factory, factor, 1, evaluation), 1);
                    fRest = qd[0];
                    continue factor_combinations;
                }

            }
            ++s;
        }

        if (!fRest.isConstant() || !fRest.cc().isConstant())
            trueFactors.addFactor(HenselLifting.denseSeriesToPoly(factory, fRest, 1, evaluation), 1);

        return trueFactors.monic();
    }

    /**
     * Factors primitive, square-free bivariate polynomial over Z
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static FactorDecomposition<MultivariatePolynomial<BigInteger>>
    bivariateDenseFactorSquareFreeZ(MultivariatePolynomial<BigInteger> poly) {
        assert poly.nUsedVariables() <= 2 && IntStream.range(2, poly.nVariables).allMatch(i -> poly.degree(i) == 0);

        if (poly.isEffectiveUnivariate())
            return factorUnivariate(poly);

        MultivariatePolynomial<BigInteger> reducedPoly = poly;
        int[] degreeBounds = reducedPoly.degrees();

        // use main variable with maximal degree
        boolean swapVariables = false;
        if (degreeBounds[1] > degreeBounds[0]) {
            swapVariables = true;
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            ArraysUtil.swap(degreeBounds, 0, 1);
        }

        MultivariatePolynomial<BigInteger> xDerivative = reducedPoly.derivative(0);
        if (xDerivative.isZero()) {
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            swapVariables = !swapVariables;
            xDerivative = reducedPoly.derivative(0);
        }

        MultivariatePolynomial<BigInteger> dGCD = MultivariateGCD.PolynomialGCD(xDerivative, reducedPoly);
        if (!dGCD.isConstant()) {
            FactorDecomposition<MultivariatePolynomial<BigInteger>>
                    gcdFactorization = bivariateDenseFactorSquareFreeZ(dGCD),
                    restFactorization = bivariateDenseFactorSquareFreeZ(divideExact(reducedPoly, dGCD));

            gcdFactorization.addAll(restFactorization);
            if (swapVariables)
                swap(gcdFactorization);

            return gcdFactorization;
        }

        // degree in main variable
        int degree = reducedPoly.degree(0);
        // substitution value for second variable
        BigInteger ySubstitution = null;
        // univariate factorization
        FactorDecomposition<UnivariatePolynomial<BigInteger>> uFactorization = null;

        // number of univariate factorizations tried
        int univariateFactorizations = 0;
        boolean tryZeroFirst = true;
        UnivariatePolynomial<BigInteger> uImage = null;
        while (univariateFactorizations < UNIVARIATE_FACTORIZATION_ATTEMPTS) {
            BigInteger substitution;
            if (tryZeroFirst) {
                // first try to substitute 0 for second variable, then use random values
                substitution = BigInteger.ZERO;
                tryZeroFirst = false;
            } else
                substitution = BigInteger.valueOf(PrivateRandom.getRandom().nextInt());

            MultivariatePolynomial<BigInteger> image = reducedPoly.evaluate(1, substitution);
            if (image.degree() != degree)
                // unlucky substitution
                continue;

            if (image.cc().isZero())
                // c.c. must not be zero since input is primitive
                // => unlucky substitution
                continue;

            uImage = image.asUnivariate();
            if (!SquareFreeFactorization.isSquareFree(uImage))
                // ensure that univariate image is also square free
                continue;

            FactorDecomposition<UnivariatePolynomial<BigInteger>> factorization = Factorization.factor(uImage);
            if (factorization.size() == 1)
                // irreducible polynomial
                return FactorDecomposition.singleFactor(poly);


            if (uFactorization == null || factorization.size() < uFactorization.size()) {
                // better univariate factorization found
                uFactorization = factorization;
                ySubstitution = substitution;
            }

            ++univariateFactorizations;
        }

        // univariate factors are calculated
        assert ySubstitution != null;

        // choose appropriate prime modulus
        int basePrime = 1 << 22;
        BigInteger bBasePrime;
        while (true) {
            basePrime = SmallPrimes.nextPrime(basePrime);
            bBasePrime = BigInteger.valueOf(basePrime);
            if (uImage.lc().mod(bBasePrime).isZero() || uImage.cc().mod(bBasePrime).isZero())
                continue;

            IntegersModulo moduloDomain = new IntegersModulo(bBasePrime);
            // ensure that univariate factors are still co-prime
            if (!CommonPolynomialsArithmetics.coprimeQ(uFactorization.map(f -> f.setDomain(moduloDomain)).factors))
                continue;

            break;
        }

        // chose prime**k which exceeds the coefficient bound

        BigInteger bound2 = coefficientsBound(poly).shiftLeft(1);
        BigInteger modulus = bBasePrime;
        while (modulus.compareTo(bound2) < 0)
            modulus = modulus.multiply(bBasePrime);
        IntegersModulo zpDomain = new IntegersModulo(modulus);

        @SuppressWarnings("unchecked")
        List<UnivariatePolynomial<BigInteger>> factorsListZp = uFactorization.map(f -> f.setDomain(zpDomain)).monic().factors;
        MultivariatePolynomial<BigInteger>
                baseZp = reducedPoly.setDomain(zpDomain),
                lcZp = baseZp.lc(0);
        baseZp = baseZp.divideOrNull(lcZp.evaluate(1, ySubstitution).lc());
        assert baseZp != null;

        // we don't precompute correct leading coefficients of bivariate factors
        // instead, we add the l.c. of the product to a list of lifting factors
        // in order to obtain correct factorization with monic factors mod (y - y0)^l
        // and then perform l.c. correction at the recombination stage

        Evaluation<BigInteger> evaluation = new Evaluation<>(poly.nVariables, new BigInteger[]{ySubstitution}, zpDomain, baseZp.ordering);
        if (!lcZp.isConstant()) {
            // add lc to lifting factors
            assert evaluation.evaluateFrom(lcZp, 1).isConstant();
            factorsListZp.add(0, factorsListZp.get(0).createOne());
        }

        // final factors to lift
        @SuppressWarnings("unchecked")
        UnivariatePolynomial<BigInteger>[] factorsZp = factorsListZp.toArray(new UnivariatePolynomial[factorsListZp.size()]);

        // lift univariate factorization
        int liftDegree = baseZp.degree(1) + 1;

        // series expansion around y = y0 for initial poly
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> baseSeriesZp =
                HenselLifting.seriesExpansionDense(UnivariatePolynomials.overZp(modulus), baseZp, 1, evaluation);

        // lifted factors (each factor represented as series around y = y0)
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>>[] liftedZp =
                HenselLifting.bivariateLiftDense(baseSeriesZp, factorsZp, liftDegree);

        if (!lcZp.isConstant())
            // drop auxiliary l.c. from factors
            liftedZp = Arrays.copyOfRange(liftedZp, 1, factorsZp.length);

        // factors are lifted => do recombination
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> baseSeriesZ =
                seriesExpansionDenseZ(reducedPoly, ySubstitution);
        FactorDecomposition<MultivariatePolynomial<BigInteger>> result = denseBivariateRecombinationZ(
                reducedPoly, baseZp, baseSeriesZ, liftedZp, evaluation, ySubstitution, zpDomain, liftDegree);

//        MultivariatePolynomial<BigInteger> tmp = baseZp.clone();
//        FactorDecomposition<MultivariatePolynomial<BigInteger>> result = denseBivariateRecombinationZ(
//                reducedPoly, Arrays.stream(liftedZp).map(p -> HenselLifting.denseSeriesToPoly(tmp, p, 1, evaluation)).toArray(MultivariatePolynomial[]::new), evaluation, liftDegree);

        if (swapVariables)
            // reconstruct original variables order
            for (int i = 0; i < result.factors.size(); i++)
                result.factors.set(i, AMultivariatePolynomial.swapVariables(result.get(i), 0, 1));

        return result;
    }

    static BigInteger coefficientsBound(MultivariatePolynomial<BigInteger> poly) {
        BigInteger uniformNorm = BigInteger.ZERO;
        for (BigInteger c : poly.coefficients())
            if (c.compareTo(uniformNorm) > 0)
                uniformNorm = c;

        int[] degrees = poly.degrees();
        int degreeSum = 0;
        BigInteger bound = BigInteger.ONE;
        for (int d : degrees) {
            degreeSum += d;
            bound = bound.multiply(BigInteger.valueOf(d).increment());
        }
        bound = bound.divide(BigInteger.ONE.shiftLeft(degrees.length));
        bound = BigIntegerArithmetics.sqrtCeil(bound);
        bound = bound.multiply(BigInteger.ONE.shiftLeft(degreeSum));
        bound = bound.multiply(uniformNorm);

        return bound;
    }

//    /**
//     * Naive dense recombination for bivariate factors
//     *
//     * @param base           multivariate polynomial factory (over Z)
//     * @param poly           series around y = y0 for base polynomial (which we attempt to factor) (over Z)
//     * @param modularFactors univariate factors mod (y-y0)^liftDegree (over Z/p)
//     * @param liftDegree     lifting degree (ideal power)
//     * @return true factorization
//     */
//    static FactorDecomposition<MultivariatePolynomial<BigInteger>> denseBivariateRecombinationZ(
//            MultivariatePolynomial<BigInteger> base,
//            MultivariatePolynomial<BigInteger>[] modularFactors,
//            Evaluation<BigInteger> evaluation,
//            int liftDegree) {
//
//
//        int[] modIndexes = naturalSequenceRef(modularFactors.length);
//        FactorDecomposition<MultivariatePolynomial<BigInteger>> trueFactors = FactorDecomposition.empty(base);
//        MultivariatePolynomial<BigInteger> fRest = base;
//        int s = 1;
//
//        factor_combinations:
//        while (2 * s <= modIndexes.length) {
//            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
//            for (int[] combination : combinations) {
//                int[] indexes = select(modIndexes, combination);
//
//                MultivariatePolynomial<BigInteger> factor = fRest.lc(0).setDomainFrom(modularFactors[0]);
//
//                for (int i : indexes)
//                    // todo:
//                    // implement IUnivariatePolynomial#multiplyLow(int)
//                    // and replace truncate(int) with multiplyLow(int)
//                    factor = evaluation.modImage(factor.multiply(modularFactors[i]), 1, liftDegree);//.truncate(liftDegree - 1);
//
//                factor = MultivariatePolynomial.asPolyZSymmetric(factor);
//
//                // get primitive part in first variable (remove R[y] content)
//                factor = factor.primitivePart(1);
//
//                MultivariatePolynomial<BigInteger>[] qd = MultivariateReduction.divideAndRemainder(fRest, factor);
//                if (qd != null && qd[1].isZero()) {
//                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
//                    trueFactors.addFactor(factor, 1);
//                    fRest = qd[0];
//                    continue factor_combinations;
//                }
//            }
//            ++s;
//        }
//
//        if (!fRest.isConstant() || !fRest.isConstant())
//            trueFactors.addFactor(fRest, 1);
//
//        return trueFactors;
//    }

    static FactorDecomposition<MultivariatePolynomial<BigInteger>> denseBivariateRecombinationZ(
            MultivariatePolynomial<BigInteger> baseZ,
            MultivariatePolynomial<BigInteger> baseZp,
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> baseSeriesZ,
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>>[] modularFactorsZp,
            Evaluation<BigInteger> evaluation,
            BigInteger ySubstitution,
            Domain<BigInteger> modulus,
            int liftDegree) {

        int[] modIndexes = naturalSequenceRef(modularFactorsZp.length);
        FactorDecomposition<MultivariatePolynomial<BigInteger>> trueFactors = FactorDecomposition.empty(baseZ);
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> fRest = baseSeriesZ;
        int s = 1;

        MultivariatePolynomial.USubstitution<BigInteger> lPowersZ = new MultivariatePolynomial.USubstitution<>(
                UnivariatePolynomial.create(Integers.Integers, ySubstitution.negate(), BigInteger.ONE),
                1, 2, baseZ.ordering);

        UnivariatePolynomials<UnivariatePolynomial<BigInteger>> moduloDomain = UnivariatePolynomials.overDomain(modulus);

        assert baseZ.equals(denseSeriesToPolyZ(baseZ, baseSeriesZ, lPowersZ));

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                UnivariatePolynomial<UnivariatePolynomial<BigInteger>> factor = lcInSeries(fRest).setDomain(moduloDomain);

                for (int i : indexes)
                    // todo:
                    // implement IUnivariatePolynomial#multiplyLow(int)
                    // and replace truncate(int) with multiplyLow(int)
                    factor = factor.multiply(modularFactorsZp[i]).truncate(liftDegree - 1);

                factor = seriesExpansionDenseZ(MultivariatePolynomial.asPolyZSymmetric(HenselLifting.denseSeriesToPoly(baseZp, factor, 1, evaluation)).primitivePart(1), ySubstitution);
                UnivariatePolynomial<UnivariatePolynomial<BigInteger>>[] qd = DivisionWithRemainder.divideAndRemainder(fRest, factor, true);
                if (qd != null && qd[1].isZero()) {
                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                    trueFactors.addFactor(denseSeriesToPolyZ(baseZ, factor, lPowersZ), 1);
                    fRest = qd[0];
                    continue factor_combinations;
                }
            }
            ++s;
        }

        if (!fRest.isConstant() || !fRest.cc().isConstant())
            if (trueFactors.size() == 0)
                trueFactors.addFactor(baseZ, 1);
            else
                trueFactors.addFactor(denseSeriesToPolyZ(baseZ, fRest, lPowersZ), 1);

        return trueFactors;
    }

    private static UnivariatePolynomial<UnivariatePolynomial<BigInteger>> seriesExpansionDenseZ(
            MultivariatePolynomial<BigInteger> poly,
            BigInteger ySubstitution) {
        int degree = poly.degree(1);
        UnivariatePolynomial<BigInteger>[] coefficients = UnivariatePolynomials.POLYNOMIALS_OVER_Z.createArray(degree + 1);
        for (int i = 0; i <= degree; i++)
            coefficients[i] = poly.seriesCoefficient(1, i).evaluate(1, ySubstitution).asUnivariate();
        return UnivariatePolynomial.create(UnivariatePolynomials.POLYNOMIALS_OVER_Z, coefficients);
    }

    private static MultivariatePolynomial<BigInteger> denseSeriesToPolyZ(
            MultivariatePolynomial<BigInteger> factory,
            UnivariatePolynomial<UnivariatePolynomial<BigInteger>> series,
            MultivariatePolynomial.USubstitution<BigInteger> linearPowers) {
        MultivariatePolynomial<BigInteger> result = factory.createZero();
        for (int i = 0; i <= series.degree(); i++) {
            MultivariatePolynomial<BigInteger> mPoly = AMultivariatePolynomial.asMultivariate(series.get(i), factory.nVariables, 0, factory.ordering);
            result = result.add(mPoly.multiply(linearPowers.pow(i)));
        }
        return result;
    }

    /** Given poly as R[x][y] transform it to R[y][x] */
    private static <uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomial<uPoly> changeDenseRepresentation(UnivariatePolynomial<uPoly> poly) {
        int xDegree = -1;
        for (int i = 0; i <= poly.degree(); i++)
            xDegree = Math.max(xDegree, poly.get(i).degree());

        UnivariatePolynomial<uPoly> result = poly.createZero();
        for (int i = 0; i <= xDegree; i++)
            result.set(i, coefficientInSeries(i, poly));
        return result;
    }

    /** Given poly as R[x][y] returns coefficient of x^xDegree which is R[y] */
    private static <uPoly extends IUnivariatePolynomial<uPoly>> uPoly
    coefficientInSeries(int xDegree, UnivariatePolynomial<uPoly> poly) {
        Domain<uPoly> domain = poly.domain;
        uPoly result = domain.getZero();
        for (int i = 0; i <= poly.degree(); i++)
            result.setFrom(i, poly.get(i), xDegree);
        return result;
    }

    /**
     * Given poly as R[x][y] returns leading coefficient of x which is R[y] viewed as R[x][y]
     * (with all coefficients constant)
     */
    private static <uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomial<uPoly> lcInSeries(UnivariatePolynomial<uPoly> poly) {
        UnivariatePolynomial<uPoly> result = poly.createZero();
        int xDegree = -1;
        for (int i = 0; i <= poly.degree(); i++)
            xDegree = Math.max(xDegree, poly.get(i).degree());

        for (int i = 0; i <= poly.degree(); i++)
            result.set(i, poly.get(i).getAsPoly(xDegree));
        return result;
    }

    /**
     * Factors primitive, square-free bivariate polynomial
     *
     * @param poly                   primitive, square-free bivariate polynomial over Zp
     * @param switchToExtensionField whether to switch to extension field if domain cardinality is too small
     * @return factor decomposition
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly>
    bivariateDenseFactorSquareFree(Poly poly, boolean switchToExtensionField) {
        if (poly instanceof lMultivariatePolynomialZp)
            return (FactorDecomposition<Poly>) bivariateDenseFactorSquareFree((lMultivariatePolynomialZp) poly, switchToExtensionField);
        else
            return (FactorDecomposition<Poly>) bivariateDenseFactorSquareFree((MultivariatePolynomial) poly, switchToExtensionField);
    }

    /* ================================ Multivariate factorization over finite fields ================================ */

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorInExtensionFieldGeneric(Poly poly, FactorizationAlgorithm<Term, Poly> algorithm) {
        if (poly instanceof lMultivariatePolynomialZp)
            return (FactorDecomposition<Poly>) factorInExtensionField((lMultivariatePolynomialZp) poly, (FactorizationAlgorithm<MonomialTerm<lUnivariatePolynomialZp>, MultivariatePolynomial<lUnivariatePolynomialZp>>) algorithm);
        else if (poly instanceof MultivariatePolynomial)
            return (FactorDecomposition<Poly>) factorInExtensionField((MultivariatePolynomial) poly, (FactorizationAlgorithm) algorithm);
        else
            throw new RuntimeException();
    }

    static final class OrderByDegrees<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** input polynomial (with variables renamed) */
        final Poly ordered;
        /** factors degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, variablesSorted, variablesMapping;
        /** number of variables in poly */
        final int nVariables;

        OrderByDegrees(Poly ordered, int[] degreeBounds, int[] variablesSorted, int nVariables) {
            this.ordered = ordered;
            this.degreeBounds = degreeBounds;
            this.variablesSorted = variablesSorted;
            this.variablesMapping = MultivariateGCD.inversePermutation(variablesSorted);
            this.nVariables = nVariables;
        }

        /** recover initial order of variables in the result */
        Poly restoreOrder(Poly factor) {
            return renameVariables(
                    factor.setNVariables(nVariables), variablesMapping);
        }

        Poly order(Poly factor) {
            return renameVariables(factor, variablesSorted);
        }
    }

    /**
     * @param poly             the poly
     * @param reduceNVariables whether to drop unused vars (making poly.nVariables smaller)
     * @param mainVariable     the main variable (will be x1)
     */
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    OrderByDegrees<Term, Poly> orderByDegrees(Poly poly, boolean reduceNVariables, int mainVariable) {
        int
                nVariables = poly.nVariables,
                degreeBounds[] = poly.degrees(); // degree bounds for lifting

        // swap variables so that the first variable will have the maximal degree,
        // and all non-used variables are at the end of poly

        int[] variables = ArraysUtil.sequence(nVariables);
        if (mainVariable != -1) {
            int mainDegree = degreeBounds[mainVariable];
            degreeBounds[mainVariable] = Integer.MAX_VALUE;
            //sort in descending order (NOTE: use stable sorting algorithm!!!)
            ArraysUtil.insertionSort(ArraysUtil.negate(degreeBounds), variables);
            //recover degreeBounds
            ArraysUtil.negate(degreeBounds);
            degreeBounds[ArraysUtil.firstIndexOf(mainVariable, variables)] = mainDegree;
        } else {
            //sort in descending order (NOTE: use stable sorting algorithm!!!)
            ArraysUtil.insertionSort(ArraysUtil.negate(degreeBounds), variables);
            ArraysUtil.negate(degreeBounds);//recover degreeBounds

            // chose the main variable in such way that the derivative
            // with respect to the main variable is not zero (avoid p-power)
            int i = 0;
            for (; i < variables.length; i++)
                if (!isPPower(poly, variables[i]))
                    break;

            if (i > 0) {
                ArraysUtil.swap(variables, 0, i);
                ArraysUtil.swap(degreeBounds, 0, i);
            }
        }

        int lastPresentVariable;
        if (reduceNVariables) {
            lastPresentVariable = 0; //recalculate lastPresentVariable
            for (; lastPresentVariable < degreeBounds.length; ++lastPresentVariable)
                if (degreeBounds[lastPresentVariable] == 0)
                    break;
            --lastPresentVariable;
        } else
            lastPresentVariable = nVariables - 1;

        poly = renameVariables(poly, variables)
                .setNVariables(lastPresentVariable + 1);

        return new OrderByDegrees<>(poly, degreeBounds, variables, nVariables);
    }

    /**
     * Factor multivariate polynomial over finite field
     *
     * @param polynomial the polynomial
     * @return factor decomposition
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorInField(final Poly polynomial) {
        if (!polynomial.isOverFiniteField())
            throw new IllegalArgumentException();

        if (polynomial.isEffectiveUnivariate())
            return factorUnivariate(polynomial);

        FactorDecomposition<Poly>
                // square-free decomposition
                sqf = MultivariateSquareFreeFactorization.SquareFreeFactorization(polynomial),
                // the result
                res = FactorDecomposition.constantFactor(sqf.constantFactor);
        for (int i = 0; i < sqf.size(); i++) {
            Poly factor = sqf.get(i);
            // factor into primitive polynomials
            FactorDecomposition<Poly> primitiveFactors = factorToPrimitive(factor);
            res.addConstantFactor(primitiveFactors.constantFactor);
            for (Poly primitiveFactor : primitiveFactors) {
                // factor each primitive polynomial
                FactorDecomposition<Poly> pFactors = factorPrimitive(primitiveFactor);
                res.addConstantFactor(pFactors.constantFactor);
                for (Poly pFactor : pFactors)
                    res.addFactor(pFactor, sqf.getExponent(i));
            }
        }
        return res;
    }

    /**
     * Factor primitive square-free multivariate polynomial over finite field
     *
     * @param polynomial the primitive square-free polynomial
     * @return factor decomposition
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorPrimitive(final Poly polynomial) {
        return factorPrimitive(polynomial, true);
    }

    /**
     * Factor primitive square-free multivariate polynomial over finite field
     *
     * @param polynomial             the primitive square-free polynomial
     * @param switchToExtensionField whether to switch to the extension field
     * @return factor decomposition
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorPrimitive(
            final Poly polynomial,
            boolean switchToExtensionField) {

        if (polynomial.isEffectiveUnivariate())
            return factorUnivariate(polynomial);

        // order the polynomial by degrees
        OrderByDegrees<Term, Poly> input = orderByDegrees(polynomial, true, -1);
        FactorDecomposition<Poly> decomposition = factorPrimitive0(input.ordered, switchToExtensionField);
        if (decomposition == null)
            return null;
        return decomposition.map(input::restoreOrder);
    }

    /** number of attempts to factor in base domain before switching to extension */
    private static final int N_FAILS_BEFORE_SWITCH_TO_EXTENSION = 32;

    static final class LeadingCoefficientData<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {

        // the following data represented as F[x2,x3,...,xN] (i.e. with x1 dropped, variables shifted)

        /** the original leading coefficient */
        final Poly lc;
        /** its square-free decomposition */
        final FactorDecomposition<Poly> lcSqFreeDecomposition;
        /** square-free part of l.c. */
        final Poly lcSqFreePart;
        /**
         * square-free part of l.c. divided into content in x_i and prim. part in x_i for different i
         * (sorted in order of decreasing content degrees, to provide optimal l.c. lifts)
         */
        final SplitContent<Term, Poly>[] lcSplits;
        /**
         * Whether the algorithm allows to fully reconstruct the leading coefficient. Leading coefficients can't be
         * fully reconstructed in some rare cases in the domains of very small characteritics, when l.c. factor
         * decomposition contains p-th powers in some variable (so that any univariate image is not square free)
         */
        final boolean fullyReconstructable;

        @SuppressWarnings("unchecked")
        LeadingCoefficientData(Poly lc) {
            lc = lc.dropVariable(0, true);
            this.lc = lc;
            this.lcSqFreeDecomposition = MultivariateSquareFreeFactorization.SquareFreeFactorization(lc);
            this.lcSqFreePart = lcSqFreeDecomposition.squareFreePart();

            ArrayList<SplitContent<Term, Poly>> splits = new ArrayList<>();
            // split sq.-free part of l.c. in the max degree variable
            Poly content = lcSqFreePart;
            main:
            while (!content.isConstant()) {
                // if there is some non trivial content, additional bivariate evaluations will be necessary

                int[] cDegrees = content.degrees();
                int[] variables = ArraysUtil.sequence(0, cDegrees.length);
                ArraysUtil.insertionSort(ArraysUtil.negate(cDegrees), variables);

                int iMax = 0;
                for (; iMax < variables.length; iMax++) {
                    int maxDegreeVariable = variables[iMax];
                    if (cDegrees[iMax] == 0)
                        break main;
                    Poly pContent = lcSqFreePart.contentExcept(maxDegreeVariable);
                    Poly primitivePart = MultivariateReduction.divideExact(lcSqFreePart, pContent);
                    if (containsPPower(primitivePart, maxDegreeVariable))
                        // ppPart if a p-power in main variable (e.g. (a^3*b + c) for characteristics 3)
                        // => any univariate image will not be square-free,  so we are not able to
                        // reconstruct the l.c. using this bivariate factorization and just skip it
                        continue;

                    splits.add(new SplitContent<>(maxDegreeVariable, pContent, primitivePart));
                    content = content.contentExcept(maxDegreeVariable);
                    break;
                }

            }

            this.lcSplits = splits.toArray(new SplitContent[splits.size()]);
            this.fullyReconstructable = content.isConstant();

            // assert that for domains of large characteristics, l.c. is always can be fully reconstructed
            assert lc.coefficientDomainCharacteristics().shiftRight(16).signum() == 0 || fullyReconstructable;
        }
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean isPPower(Poly p, int variable) {
        BigInteger characteristics = p.coefficientDomainCharacteristics();
        if (!characteristics.isInt())
            return false;
        int modulus = characteristics.intValueExact();
        if (modulus > p.degree())
            return false;
        for (Term term : p)
            if (term.exponents[variable] % modulus != 0)
                return false;
        return true;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean containsPPower(Poly p, int variable) {
        BigInteger characteristics = p.coefficientDomainCharacteristics();
        return characteristics.isInt()
                && characteristics.intValueExact() <= p.degree()
                && !MultivariateGCD.PolynomialGCD(p, p.derivative(variable)).isConstant();
    }

    static final class SplitContent<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        final int variable;
        final Poly content, primitivePart;
        final OrderByDegrees<Term, Poly> ppOrdered;

//        SplitContent(int variable, Poly poly) {
//            this.variable = variable;
//            this.content = poly.contentExcept(variable);
//            this.primitivePart = MultivariateReduction.divideExact(poly, content);
//            this.ppOrdered = orderByDegrees(primitivePart, false, variable);
//        }

        SplitContent(int variable, Poly content, Poly primitivePart) {
            this.variable = variable;
            this.content = content;
            this.primitivePart = primitivePart;
            this.ppOrdered = orderByDegrees(primitivePart, false, variable);
            assert !containsPPower(primitivePart, variable);
        }
    }

    /**
     * The main factorization algorithm in finite fields
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorPrimitive0(
            final Poly poly,
            boolean switchToExtensionField) {

        // assert that poly is at least bivariate
        assert poly.nUsedVariables() >= 2;
        // assert that degrees of variables are in the descending order
        assert poly.degree(1) > 0 && poly.degree(0) > 0;

        if (poly.nUsedVariables() == 2)
            // bivariate case
            return bivariateDenseFactorSquareFree(poly, switchToExtensionField);

        Poly xDerivative = poly.derivative(0);
        assert !xDerivative.isZero();

        Poly dGCD = MultivariateGCD.PolynomialGCD(xDerivative, poly);
        if (!dGCD.isConstant()) {
            FactorDecomposition<Poly>
                    gcdFactorization = factorPrimitive(dGCD, switchToExtensionField),
                    restFactorization = factorPrimitive(divideExact(poly, dGCD), switchToExtensionField);

            if (gcdFactorization == null || restFactorization == null) {
                assert !switchToExtensionField;
                return null;
            }

            return gcdFactorization.addAll(restFactorization);
        }

        // whether domain cardinality is less than 1024
        boolean isSmallCardinality = poly.coefficientDomainCardinality().bitLength() <= 10;

        // the leading coefficient
        Poly lc = poly.lc(0);
        LeadingCoefficientData<Term, Poly> lcData = new LeadingCoefficientData<>(lc);

        IEvaluationLoop<Term, Poly> evaluations = getEvaluations(poly);
        // number of attempts to find a suitable evaluation point
        int nAttempts = 0;
        // maximal number of bivariate factors
        int nBivariateFactors = Integer.MAX_VALUE;
        main:
        while (true) {
            // choose next evaluation
            IEvaluation<Term, Poly> evaluation = evaluations.next();

            if (evaluation == null || (nAttempts++ > N_FAILS_BEFORE_SWITCH_TO_EXTENSION && isSmallCardinality)) {
                // switch to field extension
                if (switchToExtensionField)
                    return factorInExtensionFieldGeneric(poly, MultivariateFactorization::factorPrimitive0);

                return null;
            }

            // check that evaluation does not change the rest degrees
            Poly[] images = poly.arrayNewInstance(poly.nVariables - 1);
            for (int i = 0; i < images.length; i++) {
                int variable = poly.nVariables - i - 1;
                images[i] = evaluation.evaluate(i == 0 ? poly : images[i - 1], variable);
                if (images[i].degree(variable - 1) != poly.degree(variable - 1))
                    continue main;
            }

            Poly
                    bivariateImage = images[images.length - 2],
                    univariateImage = images[images.length - 1];

            assert bivariateImage.degree(0) == poly.degree(0) && bivariateImage.degree(1) == poly.degree(1);

            // check that l.c. of bivariate image has same degree in second variable
            // as the original l.c.
            if (lc.degree(1) != bivariateImage.lc(0).degree(1))
                continue;

            // check that univariate image is also square-free
            if (!SquareFreeFactorization.isSquareFree(univariateImage.asUnivariate()))
                continue;

            // check that bivariate image is also primitive
            if (!bivariateImage.contentUnivariate(1).isConstant())
                continue;

            // factor bivariate image
            FactorDecomposition<Poly> biFactorsMain =
                    bivariateDenseFactorSquareFree(bivariateImage, false);
            if (biFactorsMain == null)
                if (switchToExtensionField)
                    return factorInExtensionFieldGeneric(poly, MultivariateFactorization::factorPrimitive0);
                else
                    return null;

            if (biFactorsMain.size() == 1)
                return FactorDecomposition.singleFactor(poly);

            if (biFactorsMain.size() > nBivariateFactors)
                // bad evaluation
                continue;

            nBivariateFactors = biFactorsMain.size();

            // array of bivariate factors for lifting
            // (polynomials in F[x1, x2])
            Poly[] biFactorsArrayMain;
            if (!lc.isConstant()) { // <= leading coefficients reconstruction

                // bring main bivariate factorization in canonical order
                // (required for one-to-one correspondence between different bivariate factorizations)
                toCanonicalSort(biFactorsMain, evaluation);
                biFactorsArrayMain = biFactorsMain.factors.toArray(poly.arrayNewInstance(biFactorsMain.size()));

                // the rest of l.c. (lc/lcFactors), will be constant at the end
                Poly lcRest = lc.clone();
                // the true leading coefficients (to be calculated)
                Poly[] lcFactors = poly.arrayNewInstance(biFactorsMain.size());
                // initialize lcFactors with constants (correct ones!)
                for (int i = 0; i < lcFactors.length; i++) {
                    lcFactors[i] = evaluation.evaluateFrom(biFactorsArrayMain[i].lc(0), 1);
                    lcRest = lcRest.divideByLC(lcFactors[i]);
                }

                // we perform additional bivariate factorizations in F[x1, x_i] for i = (3,..., N) (in special order)
                lc_reconstruction:
                for (int i = 0; i < lcData.lcSplits.length && !lcRest.isConstant(); i++) {
                    SplitContent<Term, Poly> lcSplit = lcData.lcSplits[i];
                    // x_i -- the variable to leave unevaluated in addition to the main variable x_1
                    // (+1 required to obtain indexing as in the original poly)
                    int freeVariable = 1 + lcSplit.variable;

                    IEvaluation<Term, Poly>
                            // original evaluation with shuffled variables
                            iEvaluation = evaluation.renameVariables(lcSplit.ppOrdered.variablesSorted),
                            // the evaluation for l.c. (x1 dropped)
                            ilcEvaluation = iEvaluation.dropVariable(1);


                    // target for lifting
                    Poly ppPart = lcSplit.ppOrdered.ordered;
                    if (!SquareFreeFactorization.isSquareFree(ilcEvaluation.evaluateFrom(ppPart, 1).asUnivariate())) {
                        // univariate image may be non square-free for two reasons:
                        //
                        // 1. bad evaluation (common for fields of large characteristics)
                        //    => we can try another evaluation
                        //
                        // 2. ppPart is a p-power in main variable (e.g. (a^3*b + c) for characteristics 3)
                        //    => any evaluation will lead to non square-free univariate image,
                        //    so we are not able to fully reconstruct the l.c.
                        //    --- this is not the case since we have filtered such splits in lcData,
                        //        (see #LeadingCoefficientData) so we just assert this case
                        assert !containsPPower(lcSplit.primitivePart, lcSplit.variable);

                        // <= try another evaluation
                        continue main;
                    }

                    // bivariate factors in F[x1, x_i]
                    FactorDecomposition<Poly> biFactors;
                    if (freeVariable == 1)
                        biFactors = biFactorsMain;
                    else {
                        // good image must have
                        //  - the same degree in the main variable
                        //  - the same degree in freeVariable
                        //  - l.c. with the same degree in freeVariable
                        Poly biImage = evaluation.evaluateFromExcept(poly, 1, freeVariable);
                        if (biImage.degree(0) != poly.degree(0)
                                || biImage.degree(freeVariable) != poly.degree(freeVariable)
                                || biImage.lc(0).degree(freeVariable) != lc.degree(freeVariable))
                            continue main;

                        // bivariate factors in F[x1, x_i]
                        biFactors = bivariateDenseFactorSquareFree(
                                orderByDegrees(biImage, false, -1).ordered, false);
                        if (biFactors == null)
                            if (switchToExtensionField)
                                return factorInExtensionFieldGeneric(poly, MultivariateFactorization::factorPrimitive0);
                            else
                                return null;
                        // bring in one-to-one correspondence with biFactorsMain
                        toCanonicalSort(biFactors, iEvaluation);
                    }

                    if (biFactors.size() != biFactorsMain.size()) {
                        // number of factors should be the same since polynomial is primitive
                        // => bad evaluation occurred
                        nBivariateFactors = Math.min(biFactors.size(), biFactorsMain.size());
                        continue main;
                    }

                    assert biFactors
                            .map(p -> iEvaluation.evaluateFrom(p, 1).asUnivariate()).monic()
                            .equals(biFactorsMain.map(p -> evaluation.evaluateFrom(p, 1).asUnivariate()).monic());

                    // square-free decomposition of the leading coefficients of bivariate factors
                    FactorDecomposition[] ulcFactors = (FactorDecomposition[])
                            biFactors.factors.stream()
                                    .map(f -> SquareFreeFactorization.SquareFreeFactorization(f.lc(0).asUnivariate()))
                                    .toArray(FactorDecomposition[]::new);

                    // move to GCD-free basis of sq.-f. decomposition (univariate, because fast)
                    GCDFreeBasis(ulcFactors);

                    // map to multivariate factors for further Hensel lifting
                    FactorDecomposition<Poly>[]
                            ilcFactors = Arrays.stream(ulcFactors)
                            .map(decomposition -> decomposition.map(p -> (Poly)
                                    asMultivariate((IUnivariatePolynomial) p, poly.nVariables - 1, 0, poly.ordering)))
                            .toArray(FactorDecomposition[]::new);

                    // pick unique factors from lc decompositions (complete square-free )
                    Set<Poly> ilcFactorsSet = Arrays.stream(ilcFactors)
                            .flatMap(FactorDecomposition::streamWithoutConstant)
                            .collect(Collectors.toSet());
                    Poly[] ilcFactorsSqFree = ilcFactorsSet
                            .toArray(poly.arrayNewInstance(ilcFactorsSet.size()));

                    assert ilcFactorsSqFree.length > 0;
                    assert Arrays.stream(ilcFactorsSqFree).noneMatch(Poly::isConstant);

                    // the sum of degrees of all unique factors in univariate gcd-free decomposition
                    // must be equal to the degree of primitive part we want to lift to
                    assert Arrays.stream(ilcFactorsSqFree).mapToInt(AMultivariatePolynomial::degree).reduce(0, (a, b) -> a + b)
                            == ppPart.degree(0);
//                    if (totalUDegree != ppPart.degree(0)) {
//                        assert !SquareFreeFactorization.isSquareFree(ilcEvaluation.evaluateFrom(ppPart, 1).asUnivariate());
//                        // univariate image is not square-free two reasons possible:
//                        // 1. bad evaluation (common for fields of large characteristics)
//                        //    => we can try another evaluation
//                        // 2. ppPart if a p-power in main variable (e.g. (a^3*b + c) for characteristics 3)
//                        //    => any evaluation will lead to non square-free univariate image,
//                        //    so we are not able to fully reconstruct the l.c.
//
//                        if (containsPPower(ppPart, 0))
//                            // <= we are not possible to reconstruct l.c. fully
//                            continue lc_reconstruction;
//                        else
//                            // <= try another evaluation, otherwise
//                            continue main;
//                    }

                    // we need to correct lcSqFreePrimitive (obtain correct numerical l.c.)
                    Poly ppPartLC = ilcEvaluation.evaluateFrom(ppPart.lc(0), 1);
                    Poly realLC = Arrays.stream(ilcFactorsSqFree)
                            .map(Poly::lcAsPoly)
                            .reduce(ilcFactorsSqFree[0].createOne(), Poly::multiply);

                    assert ppPartLC.isConstant();
                    assert realLC.isConstant();

                    Poly base = ppPart.clone().multiplyByLC(realLC.divideByLC(ppPartLC));
                    if (ilcFactorsSqFree.length == 1)
                        ilcFactorsSqFree[0].set(base);
                    else
                        // <= lifting leading coefficients
                        HenselLifting.multivariateLiftAutomaticLC(base, ilcFactorsSqFree, ilcEvaluation);

                    //assert multiply(ilcFactorsSqFree).monic().equals(base.clone().monic());

                    // l.c. has content in x2
                    for (int jFactor = 0; jFactor < lcFactors.length; jFactor++) {
                        Poly obtainedLcFactor = renameVariables(
                                ilcFactors[jFactor].toPolynomial(), lcSplit.ppOrdered.variablesMapping)
                                .insertVariable(0);
                        Poly commonPart = MultivariateGCD.PolynomialGCD(obtainedLcFactor, lcFactors[jFactor]);
                        Poly addon = MultivariateReduction.divideExact(obtainedLcFactor, commonPart);
                        // make addon monic when evaluated with evaluation
                        addon = addon.divideByLC(evaluation.evaluateFrom(addon, 1));
                        lcFactors[jFactor] = lcFactors[jFactor].multiply(addon);
                        lcRest = divideExact(lcRest, addon);
                    }
                }


                if (lcRest.isConstant()) {
                    // <= here we must be in _most_ cases

                    Poly base;
                    if (lcRest.isOne())
                        base = poly.clone().divideByLC(biFactorsMain.constantFactor);
                    else {
                        base = poly.clone();
                        base.divideByLC(lcRest);
                    }
                    HenselLifting.multivariateLift0(base, biFactorsArrayMain, lcFactors, evaluation, base.degrees(), 2);
                } else {
                    // <= very rare event (very small characteristics)

                    assert !lcData.fullyReconstructable;

                    // Poly lcCorrection = evaluation.evaluateFrom(lcRest, 2);
                    for (int i = 0; i < biFactorsMain.size(); i++) {
                        assert biFactorsArrayMain[i].lt().exponents[0] == biFactorsArrayMain[i].degree(0);

                        lcFactors[i].multiply(lcRest);
                        Poly correction = divideExact(evaluation.evaluateFrom(lcFactors[i], 2), biFactorsArrayMain[i].lc(0));
                        biFactorsArrayMain[i].multiply(correction);
                    }

                    Poly base = poly.clone().multiply(polyPow(lcRest, biFactorsMain.size() - 1, true));

                    assert IntStream.range(0, biFactorsMain.size()).allMatch(i -> biFactorsArrayMain[i].lc(0).equals(evaluation.evaluateFrom(lcFactors[i], 2)));
                    HenselLifting.multivariateLift0(base, biFactorsArrayMain, lcFactors, evaluation, base.degrees(), 2);

                    for (Poly factor : biFactorsArrayMain)
                        factor.set(HenselLifting.primitivePart(factor));
                }

            } else {
                Poly base;
                if (biFactorsMain.constantFactor.isOne())
                    base = poly;
                else {
                    base = poly.clone();
                    base.divideByLC(biFactorsMain.constantFactor);
                }

                biFactorsArrayMain = biFactorsMain.factors.toArray(poly.arrayNewInstance(biFactorsMain.size()));
                HenselLifting.multivariateLift0(base, biFactorsArrayMain, null, evaluation, poly.degrees(), 2);
            }

            FactorDecomposition<Poly> factorization
                    = FactorDecomposition.create(Arrays.asList(biFactorsArrayMain))
                    .monic()
                    .setConstantFactor(poly.lcAsPoly());

            Poly
                    lcNumeric = factorization.factors.stream().reduce(factorization.constantFactor.clone(), (a, b) -> a.lcAsPoly().multiply(b.lcAsPoly())),
                    ccNumeric = factorization.factors.stream().reduce(factorization.constantFactor.clone(), (a, b) -> a.ccAsPoly().multiply(b.ccAsPoly()));
            if (!lcNumeric.equals(poly.lcAsPoly()) || !ccNumeric.equals(poly.ccAsPoly()) || !factorization.toPolynomial().equals(poly)) {
                // bad bivariate factorization => recombination required
                // instead of recombination we try again with another evaluation
                // searching for good enough bivariate factorization
                nBivariateFactors = factorization.size() - 1;
                continue;
            }
            return factorization;
        }
    }

    /**
     * Brings bivariate factors in canonical order
     *
     * @param biFactors  some bivariate factorization
     * @param evaluation the evaluation point
     */
    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    void toCanonicalSort(FactorDecomposition<Poly> biFactors,
                         IEvaluation<Term, Poly> evaluation) {
        assert biFactors.exponents.sum() == biFactors.size();

        uPoly[] uFactorsArray = biFactors.map(p -> (uPoly) evaluation.evaluateFrom(p, 1).asUnivariate())
                .monic().factorsArrayWithoutLC();
        Poly[] biFactorsArray = biFactors.factorsArrayWithoutLC();
        ArraysUtil.quickSort(uFactorsArray, biFactorsArray);

        biFactors.factors.clear();
        biFactors.factors.addAll(Arrays.asList(biFactorsArray));
    }

    private static <Poly extends IGeneralPolynomial<Poly>> Poly multiply(Poly... p) {
        return p[0].createOne().multiply(p);
    }

    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    IEvaluationLoop<Term, Poly> getEvaluations(Poly factory) {
        if (factory instanceof lMultivariatePolynomialZp)
            return (IEvaluationLoop<Term, Poly>) new lEvaluationLoop((lMultivariatePolynomialZp) factory);
        else
            return (IEvaluationLoop<Term, Poly>) new EvaluationLoop((MultivariatePolynomial) factory);
    }

    interface IEvaluationLoop<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        IEvaluation<Term, Poly> next();
    }

    /** number of attempts to generate unique evaluation before switching to extension field */
    private static final int N_DIFF_EVALUATIONS_FAIL = 32;

    static final class lEvaluationLoop implements IEvaluationLoop<lMonomialTerm, lMultivariatePolynomialZp> {
        final lMultivariatePolynomialZp factory;
        final RandomGenerator rnd = PrivateRandom.getRandom();
        final TreeSet<long[]> tried = new TreeSet<>(ArraysUtil.COMPARATOR_LONG);

        lEvaluationLoop(lMultivariatePolynomialZp factory) {
            this.factory = factory;
        }

        @Override
        public lEvaluation next() {
            long[] point = new long[factory.nVariables - 1];
            int tries = 0;
            do {
                if (tries > N_DIFF_EVALUATIONS_FAIL)
                    return null;
                for (int i = 0; i < point.length; i++)
                    point[i] = factory.domain.randomElement(rnd);
                ++tries;
            } while (tried.contains(point));

            tried.add(point);
            return new lEvaluation(factory.nVariables, point, factory.domain, factory.ordering);
        }
    }

    static final class EvaluationLoop<E> implements IEvaluationLoop<MonomialTerm<E>, MultivariatePolynomial<E>> {
        final MultivariatePolynomial<E> factory;
        final RandomGenerator rnd = PrivateRandom.getRandom();
        final HashSet<ArrayRef<E>> tried = new HashSet<>();

        EvaluationLoop(MultivariatePolynomial<E> factory) {
            this.factory = factory;
        }

        @Override
        public Evaluation<E> next() {
            E[] point = factory.domain.createArray(factory.nVariables - 1);
            ArrayRef<E> array = new ArrayRef<>(point);
            int tries = 0;
            do {
                if (tries > N_DIFF_EVALUATIONS_FAIL)
                    return null;
                for (int i = 0; i < point.length; i++)
                    point[i] = factory.domain.randomElement(rnd);
                ++tries;
            } while (tried.contains(array));

            tried.add(array);
            return new Evaluation<>(factory.nVariables, point, factory.domain, factory.ordering);
        }

        private static class ArrayRef<T> {
            final T[] data;

            ArrayRef(T[] data) {this.data = data;}

            @Override
            public boolean equals(Object o) {
                if (this == o) return true;
                return !(o == null || getClass() != o.getClass())
                        && Arrays.equals(data, ((ArrayRef<?>) o).data);
            }

            @Override
            public int hashCode() {
                return Arrays.hashCode(data);
            }
        }
    }


    static <Poly extends IGeneralPolynomial<Poly>>
    void GCDFreeBasis(FactorDecomposition<Poly>[] decompositions) {
        ArrayList<FactorRef<Poly>> allFactors = new ArrayList<>();
        for (FactorDecomposition<Poly> decomposition : decompositions)
            for (int j = 0; j < decomposition.size(); j++)
                allFactors.add(new FactorRef<>(decomposition, j));

        for (int i = 0; i < allFactors.size() - 1; i++) {
            for (int j = i + 1; j < allFactors.size(); j++) {
                FactorRef<Poly>
                        a = allFactors.get(i),
                        b = allFactors.get(j);

                Poly gcd = CommonPolynomialsArithmetics.PolynomialGCD(a.factor(), b.factor());
                if (gcd.isConstant())
                    continue;

                Poly
                        aReduced = CommonPolynomialsArithmetics.PolynomialDivideAndRemainder(a.factor(), gcd)[0],
                        bReduced = CommonPolynomialsArithmetics.PolynomialDivideAndRemainder(b.factor(), gcd)[0];

                if (bReduced.isConstant())
                    allFactors.remove(j);

                a.update(aReduced, gcd);
                b.update(bReduced, gcd);

                FactorRef<Poly> gcdRef = new FactorRef<>();
                gcdRef.decompositions.addAll(a.decompositions);
                gcdRef.indexes.addAll(a.indexes);
                gcdRef.decompositions.addAll(b.decompositions);
                gcdRef.indexes.addAll(b.indexes);

                allFactors.add(gcdRef);
            }
        }

        Arrays.stream(decompositions).forEach(MultivariateFactorization::normalizeGCDFreeDecomposition);
    }

    private static <Poly extends IGeneralPolynomial<Poly>>
    void normalizeGCDFreeDecomposition(FactorDecomposition<Poly> decomposition) {
        main:
        for (int i = decomposition.factors.size() - 1; i >= 0; --i) {
            Poly factor = decomposition.get(i).clone();
            decomposition.addConstantFactor(polyPow(factor.lcAsPoly(), decomposition.getExponent(i), false));
            factor.monic();

            if (factor.isOne()) {
                decomposition.factors.remove(i);
                decomposition.exponents.removeAt(i);
                continue;
            }

            decomposition.factors.set(i, factor);

            for (int j = i + 1; j < decomposition.size(); j++) {
                if (decomposition.get(j).equals(factor)) {
                    decomposition.exponents.set(i, decomposition.exponents.get(j) + decomposition.exponents.get(i));
                    decomposition.factors.remove(j);
                    decomposition.exponents.removeAt(j);
                    continue main;
                }
            }
        }
    }

    private static final class FactorRef<Poly extends IGeneralPolynomial<Poly>> {
        final List<FactorDecomposition<Poly>> decompositions;
        final TIntArrayList indexes;

        public FactorRef() {
            this.decompositions = new ArrayList<>();
            this.indexes = new TIntArrayList();
        }

        public FactorRef(FactorDecomposition<Poly> decomposition, int index) {
            this.decompositions = new ArrayList<>();
            this.indexes = new TIntArrayList();
            decompositions.add(decomposition);
            indexes.add(index);
        }

        Poly factor() {return decompositions.get(0).get(indexes.get(0));}

        void update(Poly reduced, Poly gcd) {
            // add gcd to all required decompositions
            for (int i = 0; i < decompositions.size(); i++) {
                FactorDecomposition<Poly> decomposition = decompositions.get(i);
                decomposition.factors.set(indexes.get(i), reduced); // <- just in case
                decomposition.addFactor(gcd, decomposition.getExponent(indexes.get(i)));
            }
        }
    }
}
