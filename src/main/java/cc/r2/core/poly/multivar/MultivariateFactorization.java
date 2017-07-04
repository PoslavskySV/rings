package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.UnivariatePolynomials;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.Evaluation;
import cc.r2.core.poly.multivar.HenselLifting.IEvaluation;
import cc.r2.core.poly.multivar.HenselLifting.lEvaluation;
import cc.r2.core.poly.univar.*;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly.multivar.AMultivariatePolynomial.asMultivariate;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateFactorization {
    private MultivariateFactorization() {}


    /* ============================================== Auxiliary methods ============================================= */

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

    /* ====================================== Factorization over finite fields ====================================== */


    /** Number of univariate factorizations performed with different evaluation homomorphisms before doing Hensel lifting **/
    private static final long UNIVARIATE_FACTORIZATION_ATTEMPTS = 3;

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static FactorDecomposition<lMultivariatePolynomialZp>
    bivariateDenseFactorSquareFree(lMultivariatePolynomialZp poly) {
        assert poly.nVariables == 2;

        if (poly.isEffectiveUnivariate()) {
            int uVar = poly.univariateVariable();
            FactorDecomposition<lUnivariatePolynomialZp>
                    uFactors = Factorization.factor(poly.asUnivariate());
            FactorDecomposition<lMultivariatePolynomialZp> result = FactorDecomposition.empty(poly);
            for (int i = 0; i < uFactors.size(); i++)
                result.addFactor(
                        asMultivariate(uFactors.get(i), poly.nVariables, uVar, poly.ordering),
                        uFactors.getExponent(i));
            return result;
        }

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
                    gcdFactorization = bivariateDenseFactorSquareFree(dGCD),
                    restFactorization = bivariateDenseFactorSquareFree(MultivariateReduction.divideExact(reducedPoly, dGCD));

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
        while (univariateFactorizations < UNIVARIATE_FACTORIZATION_ATTEMPTS) {
            // first try to substitute 0 for second variable, then use random values
            long substitution = tryZeroFirst ? 0 : domain.randomElement(PrivateRandom.getRandom());
            tryZeroFirst = false;

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

        lEvaluation evaluation = new lEvaluation(2, new long[]{ySubstitution}, domain, reducedPoly.ordering);
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
                HenselLifting.bivariateLift(baseSeries, factors, liftDegree);

        if (!lc.isConstant())
            // drop auxiliary l.c. from factors
            lifted = Arrays.copyOfRange(lifted, 1, factors.length);

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
        assert poly.nVariables == 2;

        if (poly.isEffectiveUnivariate()) {
            int uVar = poly.univariateVariable();
            FactorDecomposition<UnivariatePolynomial<E>>
                    uFactors = Factorization.factor(poly.asUnivariate());
            FactorDecomposition<MultivariatePolynomial<E>> result = FactorDecomposition.empty(poly);
            for (int i = 0; i < uFactors.size(); i++)
                result.addFactor(
                        asMultivariate(uFactors.get(i), poly.nVariables, uVar, poly.ordering),
                        uFactors.getExponent(i));
            return result;
        }

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
                    gcdFactorization = bivariateDenseFactorSquareFree(dGCD),
                    restFactorization = bivariateDenseFactorSquareFree(MultivariateReduction.divideExact(reducedPoly, dGCD));

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
        while (univariateFactorizations < UNIVARIATE_FACTORIZATION_ATTEMPTS) {
            // first try to substitute 0 for second variable, then use random values
            E substitution = tryZeroFirst ? domain.getZero() : domain.randomElement(PrivateRandom.getRandom());
            tryZeroFirst = false;

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

        Evaluation<E> evaluation = new Evaluation<>(2, domain.createArray(ySubstitution), domain, reducedPoly.ordering);
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
                HenselLifting.bivariateLift(baseSeries, factors, liftDegree);

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

        if (!fRest.isConstant())
            trueFactors.addFactor(HenselLifting.denseSeriesToPoly(factory, fRest, 1, evaluation), 1);

        return trueFactors;
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
}
