package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
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

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateFactorization {
    private MultivariateFactorization() {}


    /* ============================================== Auxiliary methods ============================================= */

    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> fromUnivariate(Poly factory,
                                             int variable,
                                             FactorDecomposition<? extends IUnivariatePolynomial> uDecomposition) {
        return uDecomposition.map(f -> (Poly) asMultivariate(f, factory.nVariables, variable, factory.ordering));
    }

    /** calculates the inverse permutation */
    private static int[] inversePermutation(int[] permutation) {
        final int[] inv = new int[permutation.length];
        for (int i = permutation.length - 1; i >= 0; --i)
            inv[permutation[i]] = i;
        return inv;
    }

    /** structure with required input for GCD algorithms */
    private static final class FactorizationInput<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** input polynomial (with variables renamed) */
        final Poly reduced;
        /** gcd degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** domain cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int domainCardinality;
        /**
         * degree of irreducible univariate polynomial used to construct field extension q^n
         * (if the coefficient domain has so small cardinality so that modular algorithm will fail)
         */
        final int finiteExtensionDegree;

        FactorizationInput(Poly reduced,
                           int domainCardinality, int[] degreeBounds, int[] mapping, int lastPresentVariable,
                           int finiteExtensionDegree) {
            //assert monomialGCD == null || aReduced.domain.isOne(monomialGCD.coefficient);
            this.reduced = reduced;
            this.domainCardinality = domainCardinality;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
            this.finiteExtensionDegree = finiteExtensionDegree;
        }

        /** recover initial order of variables in the result */
        Poly restoreFactor(Poly factor) {
            return AMultivariatePolynomial.renameVariables(factor, mapping);
        }
    }

    /** prepare input for modular GCD algorithms (Brown, Zippel, LinZip) */
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorizationInput<Term, Poly> preparedFactorizationInput(Poly poly) {

        BigInteger domainSize = poly.coefficientDomainCardinality();
        // domain cardinality, i.e. number of possible random choices
        int domainCardinality = domainSize == null ? -1 : (domainSize.isInt() ? domainSize.intValue() : -1);

        // find monomial GCD
        // and remove monomial content from a and b
        poly = poly.clone();

        int
                nVariables = poly.nVariables,
                lastPresentVariable = -1, // last variable that present in both input polynomials
                degreeBounds[] = poly.degrees(); // degree bounds for lifting


        // now swap variables so that the first variable will have the maximal degree (univariate gcd is fast),
        // and all non-used variables are at the end of poly's

        int[] variables = ArraysUtil.sequence(nVariables);
        //sort in descending order
        ArraysUtil.quickSort(ArraysUtil.negate(degreeBounds), variables);
        ArraysUtil.negate(degreeBounds);//recover degreeBounds

        lastPresentVariable = 0; //recalculate lastPresentVariable
        for (; lastPresentVariable < degreeBounds.length; ++lastPresentVariable)
            if (degreeBounds[lastPresentVariable] == 0)
                break;
        --lastPresentVariable;

        poly = AMultivariatePolynomial.renameVariables(poly, variables);

        // check whether coefficient domain cardinality is large enough
        int finiteExtensionDegree = 1;
        int cardinalityBound = 5 * ArraysUtil.max(degreeBounds);
        if (domainSize != null && domainSize.isInt() && domainSize.intValueExact() < cardinalityBound) {
            long ds = domainSize.intValueExact();
            finiteExtensionDegree = 2;
            long tmp = ds;
            for (; tmp < cardinalityBound; ++finiteExtensionDegree)
                tmp = tmp * ds;
        }
        return new FactorizationInput<>(poly, domainCardinality, degreeBounds, variables, lastPresentVariable, finiteExtensionDegree);
    }

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorDecomposition<Poly> factorUnivariate(Poly poly) {
        int uVar = poly.univariateVariable();
        FactorDecomposition<? extends IUnivariatePolynomial>
                uFactors = Factorization.factor(poly.asUnivariate());
        return uFactors.map(u -> (Poly) asMultivariate(u, poly.nVariables, uVar, poly.ordering));
    }

    /* ================================= Bivariate factorization over finite fields ================================= */


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
    static FactorDecomposition<lMultivariatePolynomialZp>
    bivariateDenseFactorSquareFreeSmallCardinality(lMultivariatePolynomialZp poly) {
        lIntegersModulo domain = poly.domain;

        int startingDegree = EXTENSION_FIELD_EXPONENT;
        while (true) {
            FiniteField<lUnivariatePolynomialZp> extensionField = new FiniteField<>(
                    IrreduciblePolynomials.randomIrreduciblePolynomial(
                            domain.modulus, startingDegree++, PrivateRandom.getRandom()));

            FactorDecomposition<lMultivariatePolynomialZp> result =
                    bivariateDenseFactorSquareFreeSmallCardinality(poly, extensionField);

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
    static FactorDecomposition<lMultivariatePolynomialZp>
    bivariateDenseFactorSquareFreeSmallCardinality(lMultivariatePolynomialZp poly, FiniteField<lUnivariatePolynomialZp> extensionField) {
        FactorDecomposition<MultivariatePolynomial<lUnivariatePolynomialZp>> factorization
                = bivariateDenseFactorSquareFree(poly.mapCoefficients(extensionField, extensionField::valueOf), false);
        if (factorization == null)
            // too small extension
            return null;
        if (!factorization.constantFactor.cc().isConstant())
            return null;

        FactorDecomposition<lMultivariatePolynomialZp> result = FactorDecomposition.constantFactor(poly.createConstant(factorization.constantFactor.cc().cc()));
        for (int i = 0; i < factorization.size(); i++) {
            if (!factorization.get(i).stream().allMatch(p -> p.isConstant()))
                return null;
            result.addFactor(factorization.get(i).mapCoefficients(poly.domain, p -> p.cc()), factorization.getExponent(i));
        }
        return result;
    }

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
                    restFactorization = bivariateDenseFactorSquareFree(MultivariateReduction.divideExact(reducedPoly, dGCD), switchToExtensionField);

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
                    return bivariateDenseFactorSquareFreeSmallCardinality(poly);
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
     * Factors primitive, square-free bivariate polynomial over Zp switching to extension field
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static <E> FactorDecomposition<MultivariatePolynomial<E>>
    bivariateDenseFactorSquareFreeSmallCardinality(MultivariatePolynomial<E> poly) {
        Domain<E> domain = poly.domain;

        int startingDegree = EXTENSION_FIELD_EXPONENT;
        while (true) {
            FiniteField<UnivariatePolynomial<E>> extensionField = new FiniteField<>(
                    IrreduciblePolynomials.randomIrreduciblePolynomial(
                            domain, startingDegree++, PrivateRandom.getRandom()));

            FactorDecomposition<MultivariatePolynomial<E>> result =
                    bivariateDenseFactorSquareFreeSmallCardinality(poly, extensionField);

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
    static <E> FactorDecomposition<MultivariatePolynomial<E>>
    bivariateDenseFactorSquareFreeSmallCardinality(MultivariatePolynomial<E> poly, FiniteField<UnivariatePolynomial<E>> extensionField) {
        FactorDecomposition<MultivariatePolynomial<UnivariatePolynomial<E>>> factorization
                = bivariateDenseFactorSquareFree(poly.mapCoefficients(extensionField, c -> UnivariatePolynomial.constant(poly.domain, c)), false);
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
                    restFactorization = bivariateDenseFactorSquareFree(MultivariateReduction.divideExact(reducedPoly, dGCD), switchToExtensionField);

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
                    return bivariateDenseFactorSquareFreeSmallCardinality(poly);
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


    /* ================================ Multivariate factorization over finite fields ================================ */

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

    static final class FactorRef<Poly extends IGeneralPolynomial<Poly>> {
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

    private static final int N_FAILS_BEFORE_SWITCH_TO_EXTENSION = 64;

    @SuppressWarnings("unchecked")
    static FactorDecomposition<lMultivariatePolynomialZp>
    factorSquareFree(final lMultivariatePolynomialZp polynomial) {
        return factorSquareFree(polynomial, true);
    }

    @SuppressWarnings("unchecked")
    static FactorDecomposition<lMultivariatePolynomialZp>
    factorSquareFree(final lMultivariatePolynomialZp polynomial,
                     boolean switchToExtensionField) {
        if (polynomial.isEffectiveUnivariate())
            return factorUnivariate(polynomial);
        if (polynomial.nUsedVariables() == 2)
            return bivariateDenseFactorSquareFree(polynomial, switchToExtensionField);

        lMultivariatePolynomialZp xDerivative = polynomial.derivative(0);
        assert !xDerivative.isZero();

        lMultivariatePolynomialZp dGCD = MultivariateGCD.PolynomialGCD(xDerivative, polynomial);
        if (!dGCD.isConstant()) {
            FactorDecomposition<lMultivariatePolynomialZp>
                    gcdFactorization = factorSquareFree(dGCD, switchToExtensionField),
                    restFactorization = factorSquareFree(MultivariateReduction.divideExact(polynomial, dGCD), switchToExtensionField);

            if (gcdFactorization == null || restFactorization == null) {
                assert !switchToExtensionField;
                return null;
            }

            return gcdFactorization.addAll(restFactorization);
        }

        // the leading coefficient
        lMultivariatePolynomialZp lc = polynomial.lc(0);
        // square-free decomposition of l.c.
        FactorDecomposition<lMultivariatePolynomialZp>
                lcDecomposition = MultivariateSquareFreeFactorization.SquareFreeFactorization(lc);
        // square-free part of l.c. where the first variable is
        // dropped (i.e. R[x2, x3, ... xN] -> R[x1, x2, ... x(N-1)]),
        // in order to avoid redundant main variable in further algorithms
        lMultivariatePolynomialZp
                lcSqFreePart = lcDecomposition.squareFreePart()
                .dropVariable(0, true);

        // content of l.c.
        lMultivariatePolynomialZp
                lcSqFreeContent = lcSqFreePart.contentExcept(0),
                lcSqFreePrimitive = MultivariateReduction.divideExact(lcSqFreePart, lcSqFreeContent);

        lEvaluationLoop evaluations = new lEvaluationLoop(polynomial);
        // number of attempts to find a suitable evaluation point
        int nAttempts = 0;
        // maximal number of bivariate factors
        int nBivariateFactors = Integer.MAX_VALUE;
        main:
        while (true) {
            if (nAttempts++ > N_FAILS_BEFORE_SWITCH_TO_EXTENSION) {
                // switch to field extension
//                if (!switchToExtensionField)
//                    return null;

                throw new RuntimeException();
            }

            // choose next evaluation
            lEvaluation evaluation = evaluations.next();

            // check that evaluation does not change the rest degrees
            lMultivariatePolynomialZp[] images = new lMultivariatePolynomialZp[polynomial.nVariables - 1];
            for (int i = 0; i < images.length; i++) {
                int variable = polynomial.nVariables - i - 1;
                images[i] = evaluation.evaluate(i == 0 ? polynomial : images[i - 1], variable);
                if (images[i].degree(variable - 1) != polynomial.degree(variable - 1))
                    continue main;
            }

            lMultivariatePolynomialZp
                    bivariateImage = images[images.length - 2],
                    univariateImage = images[images.length - 1];

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
            // FIXME: search for best second variable (for lc reconstruction) ???
            FactorDecomposition<lMultivariatePolynomialZp>
                    bivariateFactors = bivariateDenseFactorSquareFree(bivariateImage);

            if (bivariateFactors.size() == 1)
                return FactorDecomposition.singleFactor(polynomial);

            if (bivariateFactors.size() > nBivariateFactors)
                // bad evaluation
                continue;

            lMultivariatePolynomialZp[] biFactors =
                    bivariateFactors.factors.toArray(new lMultivariatePolynomialZp[bivariateFactors.size()]);

            if (!lc.isConstant()) { // <= leading coefficients reconstruction

                lMultivariatePolynomialZp lcRest = lc;
                lMultivariatePolynomialZp[] factorsLC = new lMultivariatePolynomialZp[bivariateFactors.size()];

                if (lc.isMonomial() || lcSqFreePrimitive.isMonomial()) {
                    if (lc.univariateVariable() == 1) {
                        // if l.c. depends only on second variable =>
                        // bivariate factors already have correct leading coefficients
                        // and we don't need to manually impose l.c.

                        lcRest = polynomial.createOne();
                        factorsLC = null;
                    } else
                        // try to minimize the imposed l.c.
                        for (int i = 0; i < biFactors.length; i++) {
                            factorsLC[i] = biFactors[i].lc(0);
                            assert factorsLC[i].isMonomial();
                            lcRest = MultivariateReduction.divideExact(lcRest, factorsLC[i]);
                        }
                } else if (lcSqFreePrimitive.isConstant()) {
                    // no way to gather any info about l.c. of factors
                    // nothing to do (this should be rare case...)
                    for (int i = 0; i < factorsLC.length; i++)
                        factorsLC[i] = polynomial.createOne();
                } else {
                    // evaluation for lc (x1 dropped)
                    lEvaluation evaluation2 = evaluation.dropVariable(1);

                    // square free decomposition of evaluated l.c.
                    lUnivariatePolynomialZp ulcSqFreePart =
                            evaluation2.evaluateFrom(lcSqFreePart, 1).asUnivariate();

                    // check whether the univariate image of square-free part of l.c. is also square-free
                    if (!SquareFreeFactorization.isSquareFree(ulcSqFreePart)) {

                        // extraneous l.c. factors occurred
                        // => no way to reconstruct true leading coefficients

                        // FIXME: may be just set commonFactor to lc^#factors ?
                        continue;
                    }

                    // square-free decomposition of the leading coefficients of bivariate factors
                    FactorDecomposition<lUnivariatePolynomialZp>[] ulcFactors =
                            bivariateFactors.factors
                                    .stream()
                                    .map(f -> SquareFreeFactorization.SquareFreeFactorization(f.lc(0).asUnivariate()))
                                    .toArray(FactorDecomposition[]::new);

                    // move to GCD-free basis of sq.-f. decomposition (univariate, because fast)
                    GCDFreeBasis(ulcFactors);

                    // map to multivariate factors for further Hensel lifting
                    FactorDecomposition<lMultivariatePolynomialZp>[]
                            lcFactors = Arrays.stream(ulcFactors)
                            .map(decomposition ->
                                    decomposition.map(p -> (lMultivariatePolynomialZp)
                                            asMultivariate(p, polynomial.nVariables - 1, 0, polynomial.ordering)))
                            .toArray(FactorDecomposition[]::new);

                    // now we lift the l.c. decomposition
                    Set<lMultivariatePolynomialZp> sqFreeLCContent = Arrays.stream(lcFactors)
                            .flatMap(FactorDecomposition::streamWithoutConstant)
                            .collect(Collectors.toSet());
                    lMultivariatePolynomialZp[] lcContent = sqFreeLCContent
                            .toArray(new lMultivariatePolynomialZp[sqFreeLCContent.size()]);

                    assert lcContent.length > 0;
                    assert Arrays.stream(lcContent).noneMatch(lMultivariatePolynomialZp::isConstant);
//                        for (int i = 0; i < factorsLC.length; i++)
//                            factorsLC[i] = polynomial.createOne();


                    // we need to correct lcSqFreePrimitive (obtain correct numerical l.c.)
                    lMultivariatePolynomialZp lcInBase = evaluation2.evaluateFrom(lcSqFreePrimitive.lc(0), 1);
                    assert lcInBase.isConstant();

                    lMultivariatePolynomialZp lcReal = Arrays.stream(lcContent)
                            .map(lMultivariatePolynomialZp::lcAsPoly)
                            .reduce(lcContent[0].createOne(), lMultivariatePolynomialZp::multiply);
                    assert lcReal.isConstant();

                    lMultivariatePolynomialZp base = lcSqFreePrimitive.clone().multiplyByLC(lcReal.divideByLC(lcInBase));
                    if (lcContent.length == 1)
                        lcContent[0].set(base);
                    else
                        // <= lifting leading coefficients
                        HenselLifting.multivariateLiftAutomaticLC(base, lcContent, evaluation2);


                    for (int i = 0; i < factorsLC.length; i++) {
                        lMultivariatePolynomialZp factorLC = lcFactors[i].toPolynomial().insertVariable(0);
                        lMultivariatePolynomialZp flc = evaluation.evaluateFrom(factorLC, 1);
                        lMultivariatePolynomialZp alc = evaluation.evaluateFrom(biFactors[i].lc(0), 1);
                        factorsLC[i] = factorLC.multiplyByLC(alc.divideByLC(flc));
                        lcRest = MultivariateReduction.divideExact(lcRest, factorsLC[i]);
                    }
                }

                if (lcRest.isConstant()) {
                    lMultivariatePolynomialZp base;
                    if (lcRest.isOne())
                        base = polynomial.clone().divideByLC(bivariateFactors.constantFactor);
                    else {
                        base = polynomial.clone();
                        base.divideByLC(lcRest);
                    }

                    HenselLifting.multivariateLift0(base, biFactors, factorsLC, evaluation, polynomial.degrees(), 2);
                } else {
                    assert factorsLC != null;

                    lMultivariatePolynomialZp base = polynomial.clone();
//                    lMultivariatePolynomialZp uniRealLC = evaluation.evaluateFrom(lc, 2);
//                    lMultivariatePolynomialZp uniActualLC = multiply(evaluation.evaluateFrom(factorsLC, 2));
//
//                    assert uniRealLC.degree() == uniActualLC.degree() && MultivariateReduction.dividesQ(uniRealLC, uniActualLC);

                    lMultivariatePolynomialZp lcCorrection = evaluation.evaluateFrom(lcRest, 2);

                    for (lMultivariatePolynomialZp factor : biFactors) {
                        assert factor.lt().exponents[0] == factor.degree(0);
                        factor.multiply(lcCorrection);
                    }

                    base = base.multiply(polyPow(lcRest, biFactors.length - 1, true));

                    for (lMultivariatePolynomialZp factorLC : factorsLC)
                        factorLC.multiply(lcRest);

                    HenselLifting.multivariateLift0(base, biFactors, factorsLC, evaluation, base.degrees(), 2);

                    for (lMultivariatePolynomialZp factor : biFactors)
                        factor.set(HenselLifting.primitivePart(factor));
                }
            } else {
                lMultivariatePolynomialZp base;
                if (bivariateFactors.constantFactor.isOne())
                    base = polynomial;
                else {
                    base = polynomial.clone();
                    base.divideByLC(bivariateFactors.constantFactor);
                }
                HenselLifting.multivariateLift0(base, biFactors, null, evaluation, polynomial.degrees(), 2);
            }

            FactorDecomposition<lMultivariatePolynomialZp> factorization
                    = FactorDecomposition.create(Arrays.asList(biFactors))
                    .monic()
                    .setConstantFactor(polynomial.lcAsPoly());

            // FIXME: do fast check first (lc and cc)
            if (!factorization.toPolynomial().equals(polynomial)) {
                // bad bivariate factorization => recombination required
                // instead of recombination we try again with another evaluation
                // searching for good enough bivariate factorization
                nBivariateFactors = factorization.size() - 1;
                continue;
            }
            return factorization;
        }
    }

    private static <Poly extends IGeneralPolynomial<Poly>> Poly multiply(Poly... p) {
        return p[0].createOne().multiply(p);
    }

    interface IEvaluationLoop<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        IEvaluation<Term, Poly> next();
    }

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
            do {
                for (int i = 0; i < point.length; i++)
                    point[i] = factory.domain.randomElement(rnd);
            } while (tried.contains(point));

            tried.add(point);
            return new lEvaluation(factory.nVariables, point, factory.domain, factory.ordering);
        }
    }
}
