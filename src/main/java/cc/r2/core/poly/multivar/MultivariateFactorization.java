package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.UnivariatePolynomials;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.lEvaluation;
import cc.r2.core.poly.univar.*;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

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
    static FactorDecomposition<lMultivariatePolynomialZp> bivariateDenseFactorSquareFree(lMultivariatePolynomialZp poly) {
        assert poly.nVariables == 2;

        lMultivariatePolynomialZp reducedPoly = poly;
        int[] degreeBounds = reducedPoly.degrees();

        // use main variable with maximal degree
        boolean swapVariables = false;
        if (degreeBounds[1] > degreeBounds[0]) {
            swapVariables = true;
            reducedPoly = AMultivariatePolynomial.swapVariables(reducedPoly, 0, 1);
            ArraysUtil.swap(degreeBounds, 0, 1);
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
    static FactorDecomposition<lMultivariatePolynomialZp> denseBivariateRecombination(
            lMultivariatePolynomialZp factory,
            UnivariatePolynomial<lUnivariatePolynomialZp> poly,
            UnivariatePolynomial<lUnivariatePolynomialZp>[] modularFactors,
            lEvaluation evaluation,
            int liftDegree) {

        int[] modIndexes = naturalSequenceRef(modularFactors.length);
        FactorDecomposition<lMultivariatePolynomialZp> trueFactors = FactorDecomposition.empty(factory);
        UnivariatePolynomial<lUnivariatePolynomialZp> fRest = poly;
        int s = 1;

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                UnivariatePolynomial<lUnivariatePolynomialZp> mFactor = lcInSeries(fRest);
                for (int i : indexes)
                    // todo:
                    // implement IUnivariatePolynomial#multiplyLow(int)
                    // and replace truncate(int) with multiplyLow(int)
                    mFactor = mFactor.multiply(modularFactors[i]).truncate(liftDegree - 1);

                // get primitive part in first variable (remove R[y] content)
                UnivariatePolynomial<lUnivariatePolynomialZp> factor =
                        changeDenseRepresentation(
                                changeDenseRepresentation(mFactor).primitivePart());

                UnivariatePolynomial<lUnivariatePolynomialZp>[] qd = DivisionWithRemainder.divideAndRemainder(fRest, factor, true);
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
    private static UnivariatePolynomial<lUnivariatePolynomialZp>
    changeDenseRepresentation(UnivariatePolynomial<lUnivariatePolynomialZp> poly) {
        int xDegree = -1;
        for (int i = 0; i <= poly.degree(); i++)
            xDegree = Math.max(xDegree, poly.get(i).degree());

        UnivariatePolynomial<lUnivariatePolynomialZp> result = poly.createZero();
        for (int i = 0; i <= xDegree; i++)
            result.set(i, coefficientInSeries(i, poly));
        return result;
    }

    /** Given poly as R[x][y] returns coefficient of x^xDegree which is R[y] */
    private static lUnivariatePolynomialZp
    coefficientInSeries(int xDegree, UnivariatePolynomial<lUnivariatePolynomialZp> poly) {
        Domain<lUnivariatePolynomialZp> domain = poly.domain;
        lUnivariatePolynomialZp result = domain.getZero();
        for (int i = 0; i <= poly.degree(); i++)
            result.set(i, poly.get(i).get(xDegree));
        return result;
    }

    /**
     * Given poly as R[x][y] returns leading coefficient of x which is R[y] viewed as R[x][y]
     * (with all coefficients constant)
     */
    private static UnivariatePolynomial<lUnivariatePolynomialZp>
    lcInSeries(UnivariatePolynomial<lUnivariatePolynomialZp> poly) {
        Domain<lUnivariatePolynomialZp> domain = poly.domain;
        UnivariatePolynomial<lUnivariatePolynomialZp> result = poly.createZero();
        int xDegree = -1;
        for (int i = 0; i <= poly.degree(); i++)
            xDegree = Math.max(xDegree, poly.get(i).degree());

        for (int i = 0; i <= poly.degree(); i++)
            result.set(i, domain.valueOf(poly.get(i).get(xDegree)));
        return result;
    }


//    /** naive recombination for factors */
//    static FactorDecomposition<lMultivariatePolynomialZp> reconstructFactors(
//            lMultivariatePolynomialZp poly,
//            lMultivariatePolynomialZp[] modularFactors,
//            lEvaluation evaluation,
//            int liftDegree) {
//
//        int[] modIndexes = naturalSequenceRef(modularFactors.length);
//        FactorDecomposition<lMultivariatePolynomialZp> trueFactors = FactorDecomposition.empty(poly);
//        lMultivariatePolynomialZp fRest = poly;
////        AllProductsCache<lMultivariatePolynomialZp> prods = new AllProductsCache<>(modularFactors);
//        int s = 1;
//
//        factor_combinations:
//        while (2 * s <= modIndexes.length) {
//            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
//            for (int[] combination : combinations) {
//                int[] indexes = select(modIndexes, combination);
//
//                lMultivariatePolynomialZp mFactor = fRest.lc(0);
//                for (int i : indexes)
//                    mFactor = mFactor.multiply(modularFactors[i]);
////                lMultivariatePolynomialZp mFactor = prods.multiply(indexes).clone().multiply(fRest.lc(0));
//                mFactor = evaluation.modImage(mFactor, 1, liftDegree);
//                lMultivariatePolynomialZp factor = mFactor.primitivePart(1);
//
//                lMultivariatePolynomialZp[] qd = MultivariateReduction.divideAndRemainder(fRest, factor);
//                if (qd[1].isZero()) {
//                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
//                    trueFactors.addFactor(factor, 1);
//                    fRest = qd[0];
//                    continue factor_combinations;
//                }
//
//            }
//            ++s;
//        }
//
//        if (!fRest.isConstant())
//            trueFactors.addFactor(fRest, 1);
//
//        return trueFactors;
//    }

//    static FactorDecomposition<lMultivariatePolynomialZp> reconstructFactors(
//            lMultivariatePolynomialZp poly,
//            lMultivariatePolynomialZp[] modularFactors,
//            lEvaluation evaluation,
//            int liftDegree) {
//
//        ++liftDegree;
//        int[] modIndexes = naturalSequenceRef(modularFactors.length);
//        FactorDecomposition<lMultivariatePolynomialZp> trueFactors = FactorDecomposition.empty(poly);
//        lMultivariatePolynomialZp fRest = poly;
//        int s = 1;
//
//        factor_combinations:
//        while (2 * s <= modIndexes.length) {
//            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
//            lMultivariatePolynomialZp fRestLC = null;
//            for (int[] combination : combinations) {
//                int[] indexes = select(modIndexes, combination);
//
//                lMultivariatePolynomialZp mFactor = fRest.lc(0);
//                for (int i : indexes)
//                    mFactor = mFactor.multiply(modularFactors[i]);
//                mFactor = evaluation.modImage(mFactor, 1, liftDegree);
//
//                int[] restIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
//                lMultivariatePolynomialZp mRest = fRest.lc(0);
//                for (int i : restIndexes)
//                    mRest = mRest.multiply(modularFactors[i]);
//                mRest = evaluation.modImage(mRest, 1, liftDegree);
//
//                if (fRestLC == null)
//                    fRestLC = fRest.clone().multiply(fRest.lc(0));
//
//
//                if (mFactor.clone().multiply(mRest).equals(fRestLC)) {
//
//                    modIndexes = restIndexes;
//                    trueFactors.addFactor(mFactor.primitivePart(1), 1);
//                    fRest = mRest.primitivePart(1);
//                    continue factor_combinations;
//                }
//            }
//            ++s;
//        }
//
//        if (!fRest.isConstant())
//            trueFactors.addFactor(fRest, 1);
//
//        return trueFactors;
//    }

    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    FactorDecomposition<Poly> asMultivariate(
            FactorDecomposition<uPoly> uFactorization,
            int nVariables, int variable,
            Comparator<DegreeVector> ordering) {
        Poly constantFactor = AMultivariatePolynomial.asMultivariate(uFactorization.constantFactor, nVariables, variable, ordering);
        List<Poly> mFactors = uFactorization.factors.stream()
                .map(p -> AMultivariatePolynomial.<Term, Poly>asMultivariate(p, nVariables, variable, ordering))
                .collect(Collectors.toList());
        return FactorDecomposition.create(
                constantFactor,
                mFactors,
                uFactorization.exponents
        );
    }
}
