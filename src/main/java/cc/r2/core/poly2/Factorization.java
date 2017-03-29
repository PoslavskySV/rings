package cc.r2.core.poly2;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.util.ArraysUtil;

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
    private static final class FactorMonomial<T extends IMutablePolynomial<T>> {
        final T theRest, monomial;

        FactorMonomial(T theRest, T monomial) {
            this.theRest = theRest;
            this.monomial = monomial;
        }
    }

    /** factor out common monomial term (x^n) */
    private static <T extends IMutablePolynomial<T>> FactorMonomial<T> factorOutMonomial(T poly) {
        int i = poly.firstNonZeroCoefficientPosition();

        if (i == 0)
            return new FactorMonomial<>(poly, poly.createOne());
        return new FactorMonomial<>(poly.clone().shiftLeft(i), poly.createMonomial(i));
    }

    /** early check for trivial cases */
    private static lFactorDecomposition<MutablePolynomialMod> earlyFactorizationChecks(MutablePolynomialMod poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return lFactorDecomposition.oneFactor(poly.isMonic() ? poly : poly.clone().monic(), poly.lc());

        return null;
    }

    /** early check for trivial cases */
    private static bFactorDecomposition<bMutablePolynomialMod> earlyFactorizationChecks(bMutablePolynomialMod poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return bFactorDecomposition.oneFactor(poly.isMonic() ? poly : poly.clone().monic(), poly.lc());

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
        result = new lFactorDecomposition<>();
        factor(poly, result);
        return result.setNumericFactor(poly.lc());
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static bFactorDecomposition<bMutablePolynomialMod> factor(bMutablePolynomialMod poly) {
        bFactorDecomposition<bMutablePolynomialMod> result = earlyFactorizationChecks(poly);
        if (result != null)
            return result;
        result = new bFactorDecomposition<>();
        factor(poly, result);
        return result.setNumericFactor(poly.lc());
    }

    private static <T extends IMutablePolynomialZp<T>> void factor(T poly, FactorDecomposition<T> result) {
        FactorMonomial<T> base = factorOutMonomial(poly);
        result.addFactor(base.monomial, 1);

        //do square-free factorization
        FactorDecomposition<T> sqf = SquareFreeFactorization(base.theRest);
        for (int i = 0; i < sqf.size(); ++i) {
            //for each square-free factor
            T sqfFactor = sqf.get(i);
            int sqfExponent = sqf.getExponent(i);

            //do distinct-degree factorization
            FactorDecomposition<T> ddf = DistinctDegreeFactorization(sqfFactor);
            for (int j = 0; j < ddf.size(); ++j) {
                //for each distinct-degree factor
                T ddfFactor = ddf.get(j);
                int ddfExponent = ddf.getExponent(j);

                //do equal-degree factorization
                FactorDecomposition<T> edf = CantorZassenhaus(ddfFactor, ddfExponent);
                for (T irreducibleFactor : edf.factors)
                    //put final irreducible factor into the result
                    result.addFactor(irreducibleFactor.monic(), sqfExponent);
            }
        }
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
        assert nIterations < 32;
        return liftFactorization0(modulus, LongArithmetics.safePow(modulus, 1L << nIterations),
                nIterations, basePoly, modularFactors);
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
    public static bFactorDecomposition<bMutablePolynomialMod> liftFactorization(BigInteger modulus, int nIterations,
                                                                                bMutablePolynomialZ basePoly,
                                                                                bFactorDecomposition<bMutablePolynomialMod> modularFactors) {
        return new bFactorDecomposition<>(
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
    static List<bMutablePolynomialMod> liftFactorization(BigInteger modulus, int nIterations,
                                                         bMutablePolynomialZ basePoly,
                                                         List<bMutablePolynomialMod> modularFactors) {
        assert modularFactors.size() > 0 : String.format("no factors: %s  mod %s", basePoly.toStringForCopy(), modulus);
        return liftFactorization0(modulus, BigIntegerArithmetics.pow(modulus, BigInteger.ONE.shiftLeft(nIterations)),
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

        lHenselData hensel = createHenselInput(modulus, poly, aFactor, bFactor);
        assertHenselLift(hensel);

        for (i = 0; i < nIterations; ++i) {
            hensel.liftQuadratic();
            assertHenselLift(hensel);
        }

        aFactor = hensel.aFactor;
        bFactor = hensel.bFactor;

        ArrayList<MutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, aFactor.normalSymmetricForm(), modularFactors.subList(0, nHalf)));
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, bFactor.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size())));
        return result;
    }

    /** actual multifactor Hensel lifting implementation **/
    static List<bMutablePolynomialMod> liftFactorization0(BigInteger modulus, BigInteger newModulus, int nIterations,
                                                          bMutablePolynomialZ poly,
                                                          List<bMutablePolynomialMod> modularFactors) {
        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        assert nIterations > 0;

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(newModulus, true).monic());

        bMutablePolynomialMod factory = modularFactors.get(0);
        bMutablePolynomialMod
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        bHenselData hensel = createHenselInput(modulus, poly, aFactor, bFactor);
        assertHenselLift(hensel);

        for (i = 0; i < nIterations; ++i) {
            hensel.liftQuadratic();
            assertHenselLift(hensel);
        }

        aFactor = hensel.aFactor;
        bFactor = hensel.bFactor;

        ArrayList<bMutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, aFactor.normalSymmetricForm(), modularFactors.subList(0, nHalf)));
        result.addAll(liftFactorization0(modulus, newModulus, nIterations, bFactor.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size())));
        return result;
    }

//    interface HenselFactory<PolyZ extends IMutablePolynomialZ<PolyZ>,
//            PolyZp extends IMutablePolynomialZp<PolyZp>> {
//
//        HenselDataA<PolyZp> create(PolyZ poly, PolyZp aFactor, PolyZp bFactor);
//    }
//
//    static HenselFactory<MutablePolynomialZ, MutablePolynomialMod> lFactory(final long modulus) {
//        return (poly, aFactor, bFactor) -> createHenselInput(modulus, poly, aFactor, bFactor);
//    }
//
//    static HenselFactory<bMutablePolynomialZ, bMutablePolynomialMod> bFactory(final BigInteger modulus) {
//        return (poly, aFactor, bFactor) -> createHenselInput(modulus, poly, aFactor, bFactor);
//    }

    /** creates liftable quintet */
    static lHenselData createHenselInput(long modulus,
                                         MutablePolynomialZ poly,
                                         MutablePolynomialMod aFactor,
                                         MutablePolynomialMod bFactor) {
        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lHenselData(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** creates liftable quintet */
    static bHenselData createHenselInput(BigInteger modulus,
                                         bMutablePolynomialZ poly,
                                         bMutablePolynomialMod aFactor,
                                         bMutablePolynomialMod bFactor) {
        bMutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bHenselData(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static <PolyZp extends IMutablePolynomialZp<PolyZp>> PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
        PolyZp[] xgcd = ExtendedEuclid(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        assert xgcd[0].isConstant() : "bad xgcd: " + Arrays.toString(xgcd) + " for xgcd(" + a + ", " + b + ")";

        //normalize: x * a + y * b = 1
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[0].monic();

        return xgcd;
    }

    /** assertion for correct Hensel structure */
    private static void assertHenselLift(lHenselData lift) {
        assert lift.base.modulus(lift.modulus, true).equals(lift.aFactor.clone().multiply(lift.bFactor))
                : lift.base + " != (" + lift.aFactor + ") * ( " + lift.bFactor + ")" + " mod " + lift.modulus;
        assert (lift.aCoFactor == null && lift.bCoFactor == null)
                || lift.aFactor.clone().multiply(lift.aCoFactor).add(lift.bFactor.clone().multiply(lift.bCoFactor)).isOne();
    }

    /** assertion for correct Hensel structure */
    private static void assertHenselLift(bHenselData lift) {
        assert lift.base.modulus(lift.modulus, true).equals(lift.aFactor.clone().multiply(lift.bFactor))
                : lift.base + " != (" + lift.aFactor + ") * ( " + lift.bFactor + ")" + " mod " + lift.modulus;
        assert (lift.aCoFactor == null && lift.bCoFactor == null)
                || lift.aFactor.clone().multiply(lift.aCoFactor).add(lift.bFactor.clone().multiply(lift.bCoFactor)).isOne();
    }

    /** data used in Hensel lifting **/
    private static abstract class HenselDataA<PolyZp extends IMutablePolynomialZp<PolyZp>> {
        /** two factors of the initial Z[x] poly **/
        PolyZp aFactor, bFactor;
        /** xgcd coefficients **/
        PolyZp aCoFactor, bCoFactor;

        public HenselDataA(PolyZp aFactor, PolyZp bFactor, PolyZp aCoFactor, PolyZp bCoFactor) {
            this.aFactor = aFactor;
            this.bFactor = bFactor;
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
        }

        abstract void liftQuadratic();

        void henselStep0(PolyZp baseMod) {
            PolyZp e = baseMod.subtract(aFactor.clone().multiply(bFactor));

            PolyZp[] qr = divideAndRemainder(
                    aCoFactor.clone().multiply(e),
                    bFactor, false);
            PolyZp q = qr[0], r = qr[1];

            PolyZp aFactorNew = aFactor.clone()
                    .add(bCoFactor.clone().multiply(e))
                    .add(aFactor.clone().multiply(q));

            PolyZp bFactorNew = bFactor.clone().add(r);

            PolyZp b = aCoFactor.clone().multiply(aFactorNew)
                    .add(bCoFactor.clone().multiply(bFactorNew))
                    .decrement();

            PolyZp[] cd = divideAndRemainder(
                    aCoFactor.clone().multiply(b),
                    bFactorNew, false);
            PolyZp c = cd[0], d = cd[1];

            PolyZp aCoFactorNew = aCoFactor.clone().subtract(d);
            PolyZp bCoFactorNew = bCoFactor.clone()
                    .subtract(bCoFactor.clone().multiply(b))
                    .subtract(c.clone().multiply(aFactorNew));

            aFactor = aFactorNew; aCoFactor = aCoFactorNew;
            bFactor = bFactorNew; bCoFactor = bCoFactorNew;
        }
    }

    /** lift for machine-size polys **/
    static final class lHenselData extends HenselDataA<MutablePolynomialMod> {
        /** the modulus */
        long modulus;
        /** initial Z[x] poly **/
        final MutablePolynomialZ base;

        public lHenselData(long modulus, MutablePolynomialZ base, MutablePolynomialMod aFactor, MutablePolynomialMod bFactor, MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        void liftQuadratic() {
            //switch modulus
            modulus = safeMultiply(modulus, modulus);
            aFactor = aFactor.setModulusUnsafe(modulus);
            bFactor = bFactor.setModulusUnsafe(modulus);
            aCoFactor = aCoFactor.setModulusUnsafe(modulus);
            bCoFactor = bCoFactor.setModulusUnsafe(modulus);
            MutablePolynomialMod baseMod = base.modulus(modulus, true);
            henselStep0(baseMod);
        }
    }

    /** lift for BigInteger polys **/
    static final class bHenselData extends HenselDataA<bMutablePolynomialMod> {
        /** the modulus */
        BigInteger modulus;
        /** initial Z[x] poly **/
        final bMutablePolynomialZ base;

        public bHenselData(BigInteger modulus, bMutablePolynomialZ base, bMutablePolynomialMod aFactor, bMutablePolynomialMod bFactor, bMutablePolynomialMod aCoFactor, bMutablePolynomialMod bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        void liftQuadratic() {
            //switch modulus
            modulus = modulus.multiply(modulus);
            aFactor = aFactor.setModulusUnsafe(modulus);
            bFactor = bFactor.setModulusUnsafe(modulus);
            aCoFactor = aCoFactor.setModulusUnsafe(modulus);
            bCoFactor = bCoFactor.setModulusUnsafe(modulus);
            bMutablePolynomialMod baseMod = base.modulus(modulus, true);
            henselStep0(baseMod);
        }
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

    /** reconstruct true factors by enumerating all combinations */
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

    /** reconstruct true factors by enumerating all combinations */
    static bFactorDecomposition<bMutablePolynomialZ> reconstructFactorsZ(
            bMutablePolynomialZ poly,
            bFactorDecomposition<bMutablePolynomialMod> modularFactors) {

        if (modularFactors.isTrivial())
            return bFactorDecomposition.oneFactor(poly, BigInteger.ONE);

        bMutablePolynomialMod factory = modularFactors.get(0);

        int[] modIndexes = naturalSequenceRef(modularFactors.size());
        bFactorDecomposition<bMutablePolynomialZ> trueFactors = new bFactorDecomposition<>();
        bMutablePolynomialZ fRest = poly;
        int s = 1;

        factor_combinations:
        while (2 * s <= modIndexes.length) {
            IntCombinationsGenerator combinations = new IntCombinationsGenerator(modIndexes.length, s);
            for (int[] combination : combinations) {
                int[] indexes = select(modIndexes, combination);

                bMutablePolynomialMod mFactor = factory.createConstant(fRest.lc());
                for (int i : indexes)
                    mFactor = mFactor.multiply(modularFactors.get(i));
                bMutablePolynomialZ factor = mFactor.normalSymmetricForm().primitivePart();

                if (!fRest.lc().remainder(factor.lc()).isZero() || !fRest.cc().remainder(factor.cc()).isZero())
                    continue;

                bMutablePolynomialMod mRest = factory.createConstant(fRest.lc().divide(factor.lc()));
                int[] restIndexes = ArraysUtil.intSetDifference(modIndexes, indexes);
                for (int i : restIndexes)
                    mRest = mRest.multiply(modularFactors.get(i));
                bMutablePolynomialZ rest = mRest.normalSymmetricForm().primitivePart();

                if (!factor.lc().multiply(rest.lc()).equals(fRest.lc())
                        || !factor.cc().multiply(rest.cc()).equals(fRest.cc()))
                    continue;
                if (rest.clone().multiply(factor).equals(fRest)) {
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

    private static int randomModulusInf() {
        return LOWER_RND_MODULUS_BOUND + GlobalRandom.getRandom().nextInt(UPPER_RND_MODULUS_BOUND - LOWER_RND_MODULUS_BOUND);
    }

    private static int next32BitPrime(int val) {
        if (val < 0) {
            long l = BigPrimes.nextPrime(Integer.toUnsignedLong(val));
            assert fits32bitWord(l);
            return (int) l;
        } else
            return SmallPrimes.nextPrime(val);
    }

    /**
     * Factors primitive square-free polynomial using Hensel lifting
     *
     * @param poly Z[x] primitive square-free polynomial
     */
    static bFactorDecomposition<bMutablePolynomialZ> factorSquareFree(bMutablePolynomialZ poly) {
        assert poly.content().isOne();
        assert poly.lc().signum() > 0;

        BigInteger lc = poly.lc();
        BigInteger bound2 = BigInteger.TWO.multiply(poly.mignotteBound().multiply(lc.abs()));
        if (bound2.compareTo(LongArithmetics.b_MAX_SUPPORTED_MODULUS) < 0)
            return lFactorDecomposition.convert(factorSquareFree(poly.toLong()));

        // choose prime at random
        int trial32Modulus = randomModulusInf() - 1;
        long modulus;
        MutablePolynomialMod moduloImage;
        do {
            trial32Modulus = next32BitPrime(trial32Modulus + 1);
            modulus = Integer.toUnsignedLong(trial32Modulus);
            moduloImage = poly.modulus(modulus).toLong();
        } while (!SquareFreeFactorization.isSquareFree(moduloImage));

        // do modular factorization
        bFactorDecomposition<bMutablePolynomialMod> modularFactors = lFactorDecomposition.convert(factor(moduloImage.monic()));

        // do Hensel lifting
        // determine number of Hensel steps
        int henselIterations = 0;
        BigInteger bModulus = BigInteger.valueOf(modulus);
        BigInteger liftedModulus = bModulus;
        while (liftedModulus.compareTo(bound2) < 0) {
            liftedModulus = liftedModulus.multiply(liftedModulus);
            ++henselIterations;
        }
        // actual lift
        if (henselIterations > 0)
            modularFactors = liftFactorization(bModulus, henselIterations, poly, modularFactors);

        //reconstruct true factors
        return reconstructFactorsZ(poly, modularFactors);
    }


    /** determines the lower bound for the possible modulus for Zp trials */
    private static int chooseModulusLowerBound(double bound2) {
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

            infinum = (long) Math.sqrt(bound2);
        } else {
            // coefficient bound is large -> we anyway need several Hensel steps
            // so we just pick 32-bit prime at random and must use BigInteger arithmetics

            throw new IllegalArgumentException();
//            infinum = randomModulusInf();
        }

        assert infinum < MIGNOTTE_MAX_DOUBLE_32;
        return (int) infinum;
    }

    /**
     * Factors square-free polynomial using Hensel lifting
     *
     * @param poly Z[x] square-free polynomial
     */
    static lFactorDecomposition<MutablePolynomialZ> factorSquareFree(MutablePolynomialZ poly) {
        assert poly.content() == 1;
        assert poly.lc() > 0;

        long lc = poly.lc();
        double bound2 = 2.0 * poly.mignotteBound() * Math.abs(lc);
        // choose prime at random
        int trial32Modulus = chooseModulusLowerBound(bound2) - 1;
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
        while (liftedModulus < bound2) {
            liftedModulus = LongArithmetics.safeMultiply(liftedModulus, liftedModulus);
            ++henselIterations;
        }
        // actual lift
        if (henselIterations > 0)
            modularFactors = liftFactorization(modulus, henselIterations, poly, modularFactors);

        //reconstruct true factors
        return reconstructFactorsZ(poly, modularFactors);
    }

    @SuppressWarnings("unchecked")
    private static <PolyZ extends IMutablePolynomialZ<PolyZ>> FactorDecomposition<PolyZ> factorSquareFree(PolyZ poly) {
        if (poly instanceof MutablePolynomialZ)
            return (FactorDecomposition<PolyZ>) factorSquareFree((MutablePolynomialZ) poly);
        else
            return (FactorDecomposition<PolyZ>) factorSquareFree((bMutablePolynomialZ) poly);
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static lFactorDecomposition<MutablePolynomialZ> factor(MutablePolynomialZ poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return lFactorDecomposition.oneFactor(
                    poly.isMonic() ? poly : poly.clone().primitivePart(),
                    poly.isMonic() ? 1 : poly.content());

        lFactorDecomposition<MutablePolynomialZ> result = new lFactorDecomposition<>();
        long content = poly.content();
        factor(poly.clone().divideOrNull(content), result);
        return result.setNumericFactor(content);
    }

    /**
     * Factors polynomial in Zp[x].
     *
     * @param poly the polynomial
     * @return factor decomposition
     */
    public static bFactorDecomposition<bMutablePolynomialZ> factor(bMutablePolynomialZ poly) {
        if (poly.degree <= 1 || poly.isMonomial())
            return bFactorDecomposition.oneFactor(
                    poly.isMonic() ? poly : poly.clone().primitivePart(),
                    poly.isMonic() ? BigInteger.ONE : poly.content());

        bFactorDecomposition<bMutablePolynomialZ> result = new bFactorDecomposition<>();
        BigInteger content = poly.content();
        factor(poly.clone().divideOrNull(content), result);
        return result.setNumericFactor(content);
    }

    private static <T extends IMutablePolynomialZ<T>> void factor(T poly, FactorDecomposition<T> result) {
        FactorMonomial<T> base = factorOutMonomial(poly);
        result.addFactor(base.monomial, 1);

        //do square-free factorization
        FactorDecomposition<T> sqf = SquareFreeFactorization(base.theRest);
        for (int i = 0; i < sqf.size(); ++i) {
            //for each square-free factor
            T sqfFactor = sqf.get(i);
            int sqfExponent = sqf.getExponent(i);

            //do distinct-degree factorization
            FactorDecomposition<T> cz = factorSquareFree(sqfFactor);
            //do equal-degree factorization
            for (T irreducibleFactor : cz.factors)
                //put final irreducible factor into the result
                result.addFactor(irreducibleFactor, sqfExponent);
        }
    }
}
