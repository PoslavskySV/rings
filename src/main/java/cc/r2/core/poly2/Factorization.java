package cc.r2.core.poly2;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static cc.r2.core.poly2.DistinctDegreeFactorization.DistinctDegreeFactorization;
import static cc.r2.core.poly2.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly2.EqualDegreeFactorization.CantorZassenhaus;
import static cc.r2.core.poly2.LongArithmetics.fits32bitWord;
import static cc.r2.core.poly2.LongArithmetics.safeMultiply;
import static cc.r2.core.poly2.PolynomialGCD.ExtendedEuclid;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;

/**
 * Factorization of univariate polynomials with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Factorization {
    private Factorization() {}

    /* ************************** Factorization in Zp[x] ************************** */

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
    private static lFactorDecomposition<MutablePolynomialMod> earlyFactorizationChecks(MutablePolynomialMod poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return lFactorDecomposition.oneFactor(poly.isMonic() ? poly : poly.clone().monic(), poly.lc());

        return null;
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static lFactorDecomposition<MutablePolynomialMod> factor(MutablePolynomialMod poly) {
        lFactorDecomposition<MutablePolynomialMod> result = earlyFactorizationChecks(poly);
        if (result != null)
            return result;

        FactorMonomial<MutablePolynomialMod> base = factorOutMonomial(poly);

        result = new lFactorDecomposition<>();
        result.addFactor(base.monomial, 1);

        //do square-free factorization
        lFactorDecomposition<MutablePolynomialMod> sqf = SquareFreeFactorization(base.theRest);
        for (int i = 0; i < sqf.size(); ++i) {
            //for each square-free factor
            MutablePolynomialMod sqfFactor = sqf.get(i);
            int sqfExponent = sqf.getExponent(i);

            //do distinct-degree factorization
            lFactorDecomposition<MutablePolynomialMod> ddf = DistinctDegreeFactorization(sqfFactor);
            for (int j = 0; j < ddf.size(); ++j) {
                //for each distinct-degree factor
                MutablePolynomialMod ddfFactor = ddf.get(j);
                int ddfExponent = ddf.getExponent(j);

                //do equal-degree factorization
                lFactorDecomposition<MutablePolynomialMod> edf = CantorZassenhaus(ddfFactor, ddfExponent);
                for (MutablePolynomialMod irreducibleFactor : edf.factors)
                    //put final irreducible factor into the result
                    result.addFactor(irreducibleFactor.monic(), sqfExponent);
            }
        }

        return result.setNumericFactor(poly.lc());
    }


    /* ************************** Factorization in Z[x] ************************** */


    /**
     * Applies Hensel lifting to factorization. As result factorization in Z/p[x] will be lifted to factorization
     * in Z/p^(2^k)[x] where k is {@code nIterations}.
     *
     * @param modulus        the initial modulus
     * @param nIterations    number of single step Hensel iterations (so the resulting modulus is {@code modulus^{2^nIterations}}
     * @param basePoly       initial polynomial in Z[x]
     * @param modularFactors modular factorization of {@code basePoly}
     * @return factorization modulo {@code modulus^(2^nIterations)}
     */
    public static lFactorDecomposition<MutablePolynomialMod> liftFactorization(long modulus, int nIterations,
                                                                              MutablePolynomialZ basePoly,
                                                                              lFactorDecomposition<MutablePolynomialMod> modularFactors) {
        return new lFactorDecomposition<>(
                liftFactorization(modulus, nIterations, basePoly, modularFactors.factors),
                modularFactors.exponents, modularFactors.factor);
    }

    /**
     * Applies Hensel lifting to factorization. As result factorization in Z/p[x] will be lifted to factorization
     * in Z/p^(2^k)[x] where k is {@code nIterations}.
     *
     * @param modulus        the initial modulus
     * @param nIterations    number of single step Hensel iterations (so the resulting modulus is {@code modulus^{2^nIterations}}
     * @param basePoly       initial polynomial in Z[x]
     * @param modularFactors modular factorization of {@code basePoly}
     * @return factorization modulo {@code modulus^(2^nIterations)}
     */
    static List<MutablePolynomialMod> liftFactorization(long modulus, int nIterations,
                                                        MutablePolynomialZ basePoly,
                                                        List<MutablePolynomialMod> modularFactors) {
        assert modularFactors.size() > 0 : String.format("no factors: %s  mod %s", basePoly.toStringForCopy(), modulus);
        return liftFactorization0(modulus, LongArithmetics.safePow(modulus, 1L << nIterations),
                nIterations, basePoly, modularFactors);
    }

    /** actual multifactor Hensel lifting implementation **/
    static List<MutablePolynomialMod> liftFactorization0(long modulus, long newModulus, int nIterations,
                                                         MutablePolynomialZ poly,
                                                         List<MutablePolynomialMod> modularFactors) {
        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        assert nIterations > 0;

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(newModulus, true).monic());

        MutablePolynomialMod factory = modularFactors.get(0);
        MutablePolynomialMod
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        HenselData hensel = createHenselInput(modulus, poly, aFactor, bFactor);
        assertHenselLift(hensel);

        for (i = 0; i < nIterations; ++i) {
            hensel = hensel.liftQuadratic();
            assertHenselLift(hensel);
        }

        aFactor = hensel.aFactor;
        bFactor = hensel.bFactor;

        ArrayList<MutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, aFactor.normalSymmetricForm(), modularFactors.subList(0, nHalf)));
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, bFactor.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size())));
        return result;
    }

    /** creates liftable quintet */
    static HenselData createHenselInput(long modulus,
                                        MutablePolynomialZ poly,
                                        MutablePolynomialMod aFactor,
                                        MutablePolynomialMod bFactor) {
        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new HenselData(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static MutablePolynomialMod[] monicExtendedEuclid(MutablePolynomialMod a, MutablePolynomialMod b) {
        MutablePolynomialMod[] xgcd = ExtendedEuclid(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        assert xgcd[0].isConstant() : "bad xgcd: " + Arrays.toString(xgcd) + " for xgcd(" + a + ", " + b + ")";
        long inv = LongArithmetics.modInverse(xgcd[0].lc(), xgcd[0].modulus);
        for (MutablePolynomialMod e : xgcd)
            e.multiply(inv);

        return xgcd;
    }

    /** assertion for correct Hensel structure */
    private static void assertHenselLift(HenselData lift) {
        assert lift.base.modulus(lift.modulus, true).equals(lift.aFactor.clone().multiply(lift.bFactor))
                : lift.base + " != (" + lift.aFactor + ") * ( " + lift.bFactor + ")" + " mod " + lift.modulus;
        assert (lift.aCoFactor == null && lift.bCoFactor == null)
                || lift.aFactor.clone().multiply(lift.aCoFactor).add(lift.bFactor.clone().multiply(lift.bCoFactor)).isOne();
    }

    /** data used in Hensel lifting **/
    static final class HenselData {
        /** the modulus */
        final long modulus;
        /** initial Z[x] poly **/
        final MutablePolynomialZ base;
        /** it's two factors **/
        final MutablePolynomialMod aFactor, bFactor;
        /** xgcd coefficients **/
        final MutablePolynomialMod aCoFactor, bCoFactor;

        private HenselData(long modulus,
                           MutablePolynomialZ base,
                           MutablePolynomialMod aFactor, MutablePolynomialMod bFactor,
                           MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {

            assert base.lc() % modulus != 0
                    && base.degree == aFactor.degree + bFactor.degree
                    && (bCoFactor == null || bCoFactor.degree < aFactor.degree)
                    && (aCoFactor == null || aCoFactor.degree < bFactor.degree)
                    : "Bad Hensel input: " + base + " = (" + aFactor + ") (" + bFactor + ") " +
                    "xgcd: " + aCoFactor + ", " + bCoFactor;

            this.modulus = modulus;
            this.base = base;
            this.aFactor = aFactor;
            this.bFactor = bFactor;
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
        }

        /** single-step quadratic Hensel lifting */
        HenselData liftQuadratic() {
            //switch modulus
            long newModulus = safeMultiply(modulus, modulus);
            MutablePolynomialMod
                    aFactor = this.aFactor.setModulusUnsafe(newModulus),
                    bFactor = this.bFactor.setModulusUnsafe(newModulus),
                    aCoFactor = this.aCoFactor.setModulusUnsafe(newModulus),
                    bCoFactor = this.bCoFactor.setModulusUnsafe(newModulus);

            MutablePolynomialMod e = base.modulus(newModulus, true)
                    .subtract(aFactor.clone().multiply(bFactor));

            MutablePolynomialMod[] qr = divideAndRemainder(
                    aCoFactor.clone().multiply(e),
                    bFactor, false);
            MutablePolynomialMod q = qr[0], r = qr[1];

            MutablePolynomialMod aFactorNew = aFactor.clone()
                    .add(bCoFactor.clone().multiply(e))
                    .add(aFactor.clone().multiply(q));

            MutablePolynomialMod bFactorNew = bFactor.clone()
                    .add(r);

            MutablePolynomialMod b = aCoFactor.clone().multiply(aFactorNew)
                    .add(bCoFactor.clone().multiply(bFactorNew))
                    .decrement();

            MutablePolynomialMod[] cd = divideAndRemainder(
                    aCoFactor.clone().multiply(b),
                    bFactorNew, false);
            MutablePolynomialMod c = cd[0], d = cd[1];

            MutablePolynomialMod aCoFactorNew = aCoFactor.clone().subtract(d);
            MutablePolynomialMod bCoFactorNew = bCoFactor.clone()
                    .subtract(bCoFactor.clone().multiply(b))
                    .subtract(c.clone().multiply(aFactorNew));

            return new HenselData(newModulus, base, aFactorNew, bFactorNew, aCoFactorNew, bCoFactorNew);
        }

        @Override
        public String toString() {
            return "modulus=" + modulus +
                    "\n poly=" + base +
                    "\n aFactor=" + aFactor +
                    "\n bFactor=" + bFactor +
                    "\n aCoFactor=" + aCoFactor +
                    "\n bCoFactor=" + bCoFactor;
        }
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

    /** select elements by their positions */
    private static int[] select(int[] data, int[] positions) {
        int[] r = new int[positions.length];
        int i = 0;
        for (int p : positions)
            r[i++] = data[p];
        return r;
    }

    static lFactorDecomposition<MutablePolynomialZ> reconstructFactorsZ(
            MutablePolynomialZ poly,
            lFactorDecomposition<MutablePolynomialMod> modularFactors) {

        if (modularFactors.isTrivial())
            return lFactorDecomposition.oneFactor(poly, 1);

        MutablePolynomialMod factory = modularFactors.get(0);

        int[] modIndexes = naturalSequenceRef(modularFactors.size());
        lFactorDecomposition<MutablePolynomialZ> trueFactors = new lFactorDecomposition<>();
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

                if (fRest.lc() % factor.lc() != 0 || fRest.cc() % factor.cc() != 0)
                    continue;

                MutablePolynomialMod mRest = factory.createConstant(fRest.lc() / factor.lc());
                int[] restIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                for (int i : restIndexes)
                    mRest = mRest.multiply(modularFactors.get(i));
                MutablePolynomialZ rest = mRest.normalSymmetricForm().primitivePart();

                if (safeMultiply(factor.lc(), rest.lc()) != fRest.lc()
                        || safeMultiply(factor.cc(), rest.cc()) != fRest.cc())
                    continue;
                if (rest.clone().multiplyUnsafe(factor).equals(fRest)) {
                    modIndexes = restIndexes;
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

    private static final double
            MAX_PRIME_GAP = 382,
            MIGNOTTE_MAX_DOUBLE_32 = (2.0 * Integer.MAX_VALUE) - 10 * MAX_PRIME_GAP,
            MIGNOTTE_MAX_DOUBLE_64 = MIGNOTTE_MAX_DOUBLE_32 * MIGNOTTE_MAX_DOUBLE_32;

    private static final int
            LOWER_RND_MODULUS_BOUND = 1 << 24,
            UPPER_RND_MODULUS_BOUND = 1 << 30;

    /** determines the lower bound for the possible modulus for Zp trials */
    private static int findModulusLowerBound(double bound) {
        double bound2 = 2.0 * bound;

        long infinum;
        if (bound2 < MIGNOTTE_MAX_DOUBLE_32) {
            // we can use single 32-bit modulus
            infinum = (long) bound2;
        } else if (bound2 < MIGNOTTE_MAX_DOUBLE_64) {
            // we still can use machine-size words
            // two options possible:
            // 1) use 64-bit "(long) bound2" as modulus (which is more than 32-bit and thus slow)
            // 2) use 32-bit modulus and perform a single Hensel step
            // we use 2) option

            infinum = (long) Math.sqrt(bound);
        } else {
            // coefficient bound is large -> we anyway need several Hensel steps
            // so we just pick 32-bit prime at random
            infinum = new RandomDataGenerator(GlobalRandom.getRandom())
                    .nextLong(LOWER_RND_MODULUS_BOUND, UPPER_RND_MODULUS_BOUND);
        }

        assert infinum < MIGNOTTE_MAX_DOUBLE_32;
        return (int) infinum;
    }

    private static int next32BitPrime(int val) {
        if (val < 0) {
            long l = BigPrimes.nextPrime(Integer.toUnsignedLong(val));
            assert fits32bitWord(l);
            return (int) l;
        } else
            return SmallPrimes.nextPrime(val);
    }

    static lFactorDecomposition<MutablePolynomialZ> factorSquareFree(MutablePolynomialZ poly) {
        assert poly.content() == 1;
        assert poly.lc() > 0;

        long lc = poly.lc();
        double bound = Math.sqrt(poly.degree + 1) * Math.pow(2.0, poly.degree) * poly.normMax() * Math.abs(lc);

        // choose prime at random
        int trial32Modulus = findModulusLowerBound(bound) - 1;
        long modulus;
        MutablePolynomialMod moduloImage;
        do {
            trial32Modulus = next32BitPrime(trial32Modulus + 1);
            modulus = Integer.toUnsignedLong(trial32Modulus);
            moduloImage = poly.modulus(modulus, true);
        } while (!SquareFreeFactorization.isSquareFree(moduloImage));

        // do modular factorization
        lFactorDecomposition<MutablePolynomialMod> modularFactors = factor(moduloImage.monic());

        // do Hensel lifting
        // determine number of Hensel steps
        int henselIterations = 0;
        long liftedModulus = modulus;
        while (liftedModulus < bound) {
            liftedModulus = LongArithmetics.safeMultiply(liftedModulus, liftedModulus);
            ++henselIterations;
        }
        // actual lift
        if (henselIterations > 0)
            modularFactors = liftFactorization(modulus, henselIterations, poly, modularFactors);

        //reconstruct true factors
        return reconstructFactorsZ(poly, modularFactors);
    }

    static lFactorDecomposition<MutablePolynomialZ> factorBigPrime(MutablePolynomialZ poly) {
        assert poly.content() == 1;
        assert poly.lc() > 0;

        long lc = poly.lc();
        double bound = Math.sqrt(poly.degree + 1) * Math.pow(2.0, poly.degree) * poly.normMax() * Math.abs(lc);
        assert 4 * bound < MutablePolynomialMod.MAX_SUPPORTED_MODULUS;

        long modulus = (long) (2 * bound);
        MutablePolynomialMod moduloImage;
        do {
            modulus = BigPrimes.nextPrime(modulus);
            moduloImage = poly.modulus(modulus, true);
        } while (!SquareFreeFactorization.isSquareFree(moduloImage));

        lFactorDecomposition<MutablePolynomialMod> modularFactors = factor(moduloImage.monic());
        return reconstructFactorsZ(poly, modularFactors);
    }
}
