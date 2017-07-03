package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.AllProductsCache;
import cc.r2.core.poly.multivar.HenselLifting.lEvaluation;
import cc.r2.core.poly.univar.Factorization;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.SquareFreeFactorization;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
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
    private static final long UNIVARIATE_FACTORIZATION_ATTEMPTS = 23;

    static long UNI = 0, LIFT = 0, RECO = 0;

    /**
     * Factors primitive, square-free bivariate polynomial over Zp
     *
     * @param poly primitive, square-free bivariate polynomial over Zp
     * @return factor decomposition
     */
    static FactorDecomposition<lMultivariatePolynomialZp> bivariateFactorSquareFree(lMultivariatePolynomialZp poly) {
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

            long start = System.nanoTime();
            FactorDecomposition<lUnivariatePolynomialZp> factorization = Factorization.factor(uImage);
            UNI += System.nanoTime() - start;
            if (factorization.size() == 1)
                // irreducible polynomial
                return FactorDecomposition.singleFactor(poly);


            if (uFactorization == null || factorization.size() < uFactorization.size()) {
                // better univariate factorization found
                uFactorization = factorization;
                ySubstitution = substitution;
            }

//            if (ySubstitution == 0)
//                break;

            ++univariateFactorizations;
        }

        assert ySubstitution != -1;
        assert uFactorization.factors.stream().allMatch(lUnivariatePolynomialZp::isMonic);

        // univariate factors are calculated
        @SuppressWarnings("unchecked")
        List<lMultivariatePolynomialZp> factorList = AMultivariatePolynomial.asMultivariate(
                (List) uFactorization.factors, reducedPoly.nVariables, 0, reducedPoly.ordering);

        // we don't precompute correct leading coefficients of bivariate factors
        // instead, we add the l.c. of the product to a list of lifting factors
        // in order to obtain correct factorization with monic factors mod (y - y0)^l
        // and then perform l.c. correction at the recombination stage

        lMultivariatePolynomialZp lc = reducedPoly.lc(0);
        if (!lc.isConstant())
            // add lc to lifting factors
            factorList.add(0, lc.clone());
        else
            factorList.get(0).multiply(lc);

        // final factors to lift
        lMultivariatePolynomialZp[] factors = factorList.toArray(new lMultivariatePolynomialZp[factorList.size()]);

        // lift univariate factorization
        int liftDegree = reducedPoly.degree(1) + 1;
        lEvaluation evaluation = new lEvaluation(2, new long[]{ySubstitution}, domain, reducedPoly.ordering);
//        System.out.println(ySubstitution);
        long start = System.nanoTime();
        HenselLifting.bivariateLift(reducedPoly, factors, evaluation, liftDegree);
        LIFT += System.nanoTime() - start;

        assert reducedPoly.equals(evaluation.modImage(poly.createOne().multiply(factors), 1, liftDegree));

//        System.out.println(factors.length);
        if (!lc.isConstant())
            // drop auxiliary l.c. from factors
            factors = Arrays.copyOfRange(factors, 1, factors.length);

        // factors are lifted => do recombination
        start = System.nanoTime();
        FactorDecomposition<lMultivariatePolynomialZp> result = reconstructFactors(reducedPoly, factors, evaluation, liftDegree);
        RECO += System.nanoTime() - start;

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

    /** naive recombination for factors */
    static FactorDecomposition<lMultivariatePolynomialZp> reconstructFactors(
            lMultivariatePolynomialZp poly,
            lMultivariatePolynomialZp[] modularFactors,
            lEvaluation evaluation,
            int liftDegree) {

        int[] modIndexes = naturalSequenceRef(modularFactors.length);
        FactorDecomposition<lMultivariatePolynomialZp> trueFactors = FactorDecomposition.empty(poly);
        lMultivariatePolynomialZp fRest = poly;
//        AllProductsCache<lMultivariatePolynomialZp> prods = new AllProductsCache<>(modularFactors);
        int s = 1;

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                lMultivariatePolynomialZp mFactor = fRest.lc(0);
                for (int i : indexes)
                    mFactor = mFactor.multiply(modularFactors[i]);
//                lMultivariatePolynomialZp mFactor = prods.multiply(indexes).clone().multiply(fRest.lc(0));
                mFactor = evaluation.modImage(mFactor, 1, liftDegree);
                lMultivariatePolynomialZp factor = mFactor.primitivePart(1);

                lMultivariatePolynomialZp[] qd = MultivariateReduction.divideAndRemainder(fRest, factor);
                if (qd[1].isZero()) {
                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                    trueFactors.addFactor(factor, 1);
                    fRest = qd[0];
                    continue factor_combinations;
                }

            }
            ++s;
        }

        if (!fRest.isConstant())
            trueFactors.addFactor(fRest, 1);

        return trueFactors;
    }

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
