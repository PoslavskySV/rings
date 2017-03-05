package cc.r2.core.poly2;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.util.ArraysUtil;

import static cc.r2.core.poly2.DistinctDegreeFactorization.DistinctDegreeFactorization;
import static cc.r2.core.poly2.EqualDegreeFactorization.CantorZassenhaus;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;

/**
 * Factorization of univariate polynomials with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Factorization {
    private Factorization() {}


    /** x^n * poly */
    private static final class FactorMonomial<T extends MutablePolynomialAbstract<T>> {
        final T theRest, monomial;

        FactorMonomial(T theRest, T monomial) {
            this.theRest = theRest;
            this.monomial = monomial;
        }
    }

    /** factor out common monomial term (x^n) */
    private static <T extends MutablePolynomialAbstract<T>> FactorMonomial<T> factorOutMonomial(T poly) {
        int i = 0;
        while (poly.data[i] == 0) ++i;
        assert i < poly.data.length;

        if (i == 0)
            return new FactorMonomial<>(poly, poly.createOne());
        return new FactorMonomial<>(poly.clone().shiftLeft(i), poly.createMonomial(1, i));
    }


    /** early check for trivial cases */
    private static <T extends MutablePolynomialAbstract<T>> FactorDecomposition<T> earlyFactorizationChecks(T poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return FactorDecomposition.oneFactor(poly, 1);

        return null;
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static FactorDecomposition<MutablePolynomialMod> factor(MutablePolynomialMod poly) {
        FactorDecomposition<MutablePolynomialMod> result = earlyFactorizationChecks(poly);
        if (result != null)
            return result;

        FactorMonomial<MutablePolynomialMod> base = factorOutMonomial(poly);

        result = new FactorDecomposition<>();
        result.addFactor(base.monomial, 1);

        //do square-free factorization
        FactorDecomposition<MutablePolynomialMod> sqf = SquareFreeFactorization(base.theRest);
        for (int i = 0; i < sqf.size(); ++i) {
            //for each square-free factor
            MutablePolynomialMod sqfFactor = sqf.get(i);
            int sqfExponent = sqf.getExponent(i);

            //do distinct-degree factorization
            FactorDecomposition<MutablePolynomialMod> ddf = DistinctDegreeFactorization(sqfFactor);
            for (int j = 0; j < ddf.size(); ++j) {
                //for each distinct-degree factor
                MutablePolynomialMod ddfFactor = ddf.get(j);
                int ddfExponent = ddf.getExponent(j);

                //do equal-degree factorization
                FactorDecomposition<MutablePolynomialMod> edf = CantorZassenhaus(ddfFactor, ddfExponent);
                for (MutablePolynomialMod irreducibleFactor : edf.factors)
                    //put final irreducible factor into the result
                    result.addFactor(irreducibleFactor.monic(), sqfExponent);
            }
        }

        return result.setNumericFactor(poly.lc());
    }

    /** cache of references **/
    private static int[][] naturalSequenceRefCache = new int[32][];

    /** returns sequence of natural numbers */
    private static int[] naturalSequenceRef(int n) {
        if (naturalSequenceRefCache[n] != null)
            return naturalSequenceRefCache[n];
        int[] r = new int[n];
        for (int i = 0; i < n; i++)
            r[i] = i;
        return naturalSequenceRefCache[n] = r;
    }

    private static int[] select(int[] data, int[] positions) {
        int[] r = new int[positions.length];
        int i = 0;
        for (int p : positions)
            r[i++] = data[p];
        return r;
    }

    private static long lcMultiply(FactorDecomposition<MutablePolynomialMod> factors, int[] indexes) {
        long factor = 1;
        for (int i : indexes)
            factor = LongArithmetics.safeMultiply(factor, factors.get(i).symMod(factors.get(i).lc()));
        long overall = 1;
        for (MutablePolynomialMod f : factors)
            overall = LongArithmetics.safeMultiply(overall, f.symMod(f.lc()));

        return 0;
    }

    static FactorDecomposition<MutablePolynomialZ> reconstructFactorsZ(
            MutablePolynomialZ poly,
            FactorDecomposition<MutablePolynomialMod> modularFactors) {

        if (modularFactors.isTrivial())
            return FactorDecomposition.oneFactor(poly, 1);

        MutablePolynomialMod factory = modularFactors.get(0);

        int[] modIndexes = naturalSequenceRef(modularFactors.size());
        FactorDecomposition<MutablePolynomialZ> trueFactors = new FactorDecomposition<>();
        MutablePolynomialZ fRest = poly;
        int s = 1;

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                MutablePolynomialMod mFactor = factory.createConstant(fRest.lc());
                for (int i : indexes)
                    mFactor = mFactor.multiply(modularFactors.get(i));
                MutablePolynomialZ factor = mFactor.normalSymmetricForm().primitivePart();

                if (fRest.lc() % factor.lc() != 0)
                    continue;

                MutablePolynomialMod mRest = factory.createConstant(fRest.lc() / factor.lc());
                int[] restIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                for (int i : restIndexes)
                    mRest = mRest.multiply(modularFactors.get(i));
                MutablePolynomialZ rest = mRest.normalSymmetricForm().primitivePart();
                boolean fac;
                try {
                    fac = rest.clone().multiply(factor).equals(fRest);
                } catch (ArithmeticException e) {
                    continue;
                }

                if (fac) {
                    modIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                    trueFactors.addFactor(factor, 1);
                    fRest = rest.primitivePart();
                    continue factor_combinations;
                }

            }
            ++s;
        }

        if (!fRest.isConstant())
            trueFactors.addFactor(fRest, 1);

        return trueFactors;
    }

    static FactorDecomposition<MutablePolynomialZ> factorBigPrime(MutablePolynomialZ poly) {
        assert poly.content() == 1;
        assert poly.lc() > 0;

        long lc = poly.lc();
        double bound = Math.sqrt(poly.degree + 1) * Math.pow(2.0, poly.degree) * poly.normMax() * Math.abs(lc);
        assert 4 * bound < LongModularArithmetics.MAX_SUPPORTED_MODULUS;

        long modulus = (long) (2 * bound);
        MutablePolynomialMod moduloImage;
        do {
            modulus = BigPrimes.nextPrime(modulus);
            moduloImage = poly.modulus(modulus, true);
        } while (!SquareFreeFactorization.isSquareFree(moduloImage));

        FactorDecomposition<MutablePolynomialMod> modularFactors = factor(moduloImage.monic());
        return reconstructFactorsZ(poly, modularFactors);
    }
}
