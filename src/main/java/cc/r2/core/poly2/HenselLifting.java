package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly2.DivisionWithRemainder.InverseModMonomial;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static cc.r2.core.poly2.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.LongArithmetics.safeMultiply;
import static cc.r2.core.poly2.PolynomialGCD.ExtendedEuclid;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class HenselLifting {
    private HenselLifting() {}

    /* ************************************ Factory methods ************************************ */

    /** creates liftable quintet */
    public static lQuadraticLift createQuadraticLift(long modulus,
                                                     MutablePolynomialZ poly,
                                                     MutablePolynomialMod aFactor,
                                                     MutablePolynomialMod bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.mod(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());
        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lQuadraticLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** creates liftable quintet */
    public static bQuadraticLift createQuadraticLift(BigInteger modulus,
                                                     bMutablePolynomialZ poly,
                                                     bMutablePolynomialMod aFactor,
                                                     bMutablePolynomialMod bFactor) {
        bFactor = ensureMonic(bFactor);
        if (!aFactor.mod(poly.lc()).equals(aFactor.lc()))
            aFactor = aFactor.clone().monic(poly.lc());
        bMutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bQuadraticLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** creates liftable quintet */
    public static bQuadraticLift createQuadraticLift(BigInteger modulus,
                                                     bMutablePolynomialZ poly,
                                                     MutablePolynomialMod aFactor,
                                                     MutablePolynomialMod bFactor) {
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(modulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);
        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bQuadraticLift(modulus, poly, aFactor.toBigPoly(), bFactor.toBigPoly(), xgcd[1].toBigPoly(), xgcd[2].toBigPoly());
    }


    public static lLinearLift createLinearLift(BigInteger modulus,
                                               MutablePolynomialZ poly,
                                               MutablePolynomialMod aFactor,
                                               MutablePolynomialMod bFactor) {
        return createLinearLift(modulus.longValueExact(), poly, aFactor, bFactor);
    }

    public static bLinearLift createLinearLift(BigInteger modulus,
                                               bMutablePolynomialZ poly,
                                               MutablePolynomialMod aFactor,
                                               MutablePolynomialMod bFactor) {
        return createLinearLift(modulus.longValueExact(), poly, aFactor, bFactor);
    }

    /** creates liftable quintet */
    public static lLinearLift createLinearLift(long modulus,
                                               MutablePolynomialZ poly,
                                               MutablePolynomialMod aFactor,
                                               MutablePolynomialMod bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.mod(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());

        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lLinearLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** creates liftable quintet */
    public static bLinearLift createLinearLift(long modulus,
                                               bMutablePolynomialZ poly,
                                               MutablePolynomialMod aFactor,
                                               MutablePolynomialMod bFactor) {
        BigInteger bModulus = BigInteger.valueOf(modulus);
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(bModulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);

        MutablePolynomialMod[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bLinearLift(bModulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    static <PolyZp extends IMutablePolynomialZp<PolyZp>> PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
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

    private static <PolyZp extends IMutablePolynomialZp<PolyZp>> PolyZp ensureMonic(PolyZp p) {
        return p.isMonic() ? p : p.clone().monic();
    }

    private static long[] nIterations(long modulus, long desiredBound, boolean quadratic) {
        int nIterations = 0;
        long tmp = modulus;
        while (tmp < desiredBound) {
            tmp = LongArithmetics.safeMultiply(tmp, quadratic ? tmp : modulus);
            ++nIterations;
        }
        return new long[]{nIterations, tmp};
    }

    /** Lift for word-sized polynomials **/
    public static List<MutablePolynomialMod> liftFactorization(long modulus,
                                                               long desiredBound,
                                                               MutablePolynomialZ poly,
                                                               List<MutablePolynomialMod> modularFactors,
                                                               boolean quadraticLift) {
        long[] im = nIterations(modulus, desiredBound, quadraticLift);
        return liftFactorization(modulus, im[1], (int) im[0], poly, modularFactors, quadraticLift);
    }

    /** Lift for word-sized polynomials **/
    public static List<MutablePolynomialMod> liftFactorization(long modulus,
                                                               long finalModulus,
                                                               int nIterations,
                                                               MutablePolynomialZ poly,
                                                               List<MutablePolynomialMod> modularFactors,
                                                               boolean quadraticLift) {
        assert nIterations > 0;

        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(finalModulus, true).monic());

        MutablePolynomialMod factory = modularFactors.get(0);
        MutablePolynomialMod
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        LiftableQuintet<MutablePolynomialMod> hensel = quadraticLift
                ? createQuadraticLift(modulus, poly, aFactor, bFactor)
                : createLinearLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        MutablePolynomialMod aFactorRaised = hensel.aFactorMod();
        MutablePolynomialMod bFactorRaised = hensel.bFactorMod();

        ArrayList<MutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, aFactorRaised.normalSymmetricForm(), modularFactors.subList(0, nHalf), quadraticLift));
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, bFactorRaised.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size()), quadraticLift));
        return result;
    }

    public interface LiftFactory<PolyZp extends IMutablePolynomialZp<PolyZp>> {
        LiftableQuintet<bMutablePolynomialMod> createLift(BigInteger modulus, bMutablePolynomialZ polyZ, PolyZp aFactor, PolyZp bFactor);
    }

    /** actual multifactor Hensel lifting implementation **/
    public static <PolyZp extends IMutablePolynomialZp<PolyZp>>
    List<bMutablePolynomialMod> liftFactorization0(BigInteger modulus,
                                                   BigInteger finalModulus,
                                                   int nIterations,
                                                   bMutablePolynomialZ poly,
                                                   List<PolyZp> modularFactors,
                                                   LiftFactory<PolyZp> liftFactory) {
        assert nIterations > 0;

        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(finalModulus, true).monic());

        PolyZp factory = modularFactors.get(0);
        PolyZp
                aFactor = factory.createOne(),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        LiftableQuintet<bMutablePolynomialMod> hensel = liftFactory.createLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        bMutablePolynomialMod aFactorRaised = hensel.aFactorMod();
        bMutablePolynomialMod bFactorRaised = hensel.bFactorMod();

        ArrayList<bMutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorization0(modulus, finalModulus, nIterations, aFactorRaised.normalSymmetricForm(), modularFactors.subList(0, nHalf), liftFactory));
        result.addAll(liftFactorization0(modulus, finalModulus, nIterations, bFactorRaised.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size()), liftFactory));
        return result;
    }

    static final class LiftingInfo {
        final int nIterations;
        final BigInteger finalModulus;

        private LiftingInfo(int nIterations, BigInteger finalModulus) {
            this.nIterations = nIterations;
            this.finalModulus = finalModulus;
        }
    }

    static LiftingInfo nIterations(BigInteger modulus, BigInteger desiredBound, boolean quadratic) {
        int nIterations = 0;
        BigInteger finalModulus = modulus;
        while (finalModulus.compareTo(desiredBound) < 0) {
            finalModulus = finalModulus.multiply(quadratic ? finalModulus : modulus);
            ++nIterations;
        }
        return new LiftingInfo(nIterations, finalModulus);
    }

    public static List<bMutablePolynomialMod> liftFactorization(BigInteger modulus,
                                                                BigInteger desiredBound,
                                                                bMutablePolynomialZ poly,
                                                                List<MutablePolynomialMod> modularFactors,
                                                                boolean quadratic) {
        LiftingInfo im = nIterations(modulus, desiredBound, quadratic);
        if (im.nIterations == 0)
            return modularFactors.stream().map(MutablePolynomialMod::toBigPoly).collect(Collectors.toList());
        LiftFactory<MutablePolynomialMod> factory = quadratic ? HenselLifting::createQuadraticLift : HenselLifting::createLinearLift;
        return liftFactorization0(modulus, im.finalModulus, im.nIterations, poly, modularFactors, factory);
    }

    public static List<bMutablePolynomialMod> liftFactorizationQuadratic(BigInteger modulus,
                                                                         BigInteger desiredBound,
                                                                         bMutablePolynomialZ poly,
                                                                         List<bMutablePolynomialMod> modularFactors) {
        LiftingInfo im = nIterations(modulus, desiredBound, true);
        return liftFactorization0(modulus, im.finalModulus, im.nIterations, poly, modularFactors, HenselLifting::createQuadraticLift);
    }

    /** actual multifactor Hensel lifting implementation **/
    public static List<bMutablePolynomialMod> liftFactorizationAdaptive(long initialModulus,
                                                                        BigInteger desiredBound,
                                                                        bMutablePolynomialZ poly,
                                                                        List<MutablePolynomialMod> modularFactors) {
        return liftFactorizationAdaptive(poly, modularFactors, new AdaptiveLift(initialModulus, desiredBound));
    }

    /** actual multifactor Hensel lifting implementation **/
    private static List<bMutablePolynomialMod> liftFactorizationAdaptive(bMutablePolynomialZ poly,
                                                                         List<MutablePolynomialMod> modularFactors,
                                                                         AdaptiveLift lifter) {
        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(lifter.finalModulus, true).monic());

        MutablePolynomialMod factory = modularFactors.get(0);
        MutablePolynomialMod
                aFactor = factory.createOne(),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        bMutablePolynomialMod[] lifted = lifter.lift(poly, aFactor, bFactor);
        bMutablePolynomialMod aFactorRaised = lifted[0];
        bMutablePolynomialMod bFactorRaised = lifted[1];

        ArrayList<bMutablePolynomialMod> result = new ArrayList<>();
        result.addAll(liftFactorizationAdaptive(aFactorRaised.normalSymmetricForm(), modularFactors.subList(0, nHalf), lifter));
        result.addAll(liftFactorizationAdaptive(bFactorRaised.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size()), lifter));
        return result;
    }

    private static final int SWITCH_TO_QUADRATIC_LIFT = 64;

    private static final class AdaptiveLift {
        final long initialModulus;
        final BigInteger finalModulus;
        final int nLinearIterations, nQuadraticIterations;

        public AdaptiveLift(long initialModulus, BigInteger desiredBound) {
            this.initialModulus = initialModulus;
            LiftingInfo
                    nLinearIterations = nIterations(BigInteger.valueOf(initialModulus), desiredBound, false);

            if (nLinearIterations.nIterations < SWITCH_TO_QUADRATIC_LIFT) {
                this.nLinearIterations = nLinearIterations.nIterations;
                this.nQuadraticIterations = -1;
                this.finalModulus = nLinearIterations.finalModulus;
            } else {
                LiftingInfo nQuadraticIterations = nIterations(BigInteger.valueOf(initialModulus), desiredBound, true);
                this.nLinearIterations = -1;
                this.nQuadraticIterations = nQuadraticIterations.nIterations;
                this.finalModulus = nQuadraticIterations.finalModulus;
            }
        }

        bMutablePolynomialMod[] lift(bMutablePolynomialZ poly, MutablePolynomialMod a, MutablePolynomialMod b) {
            boolean quadratic = nLinearIterations == -1;
            LiftableQuintet<bMutablePolynomialMod> lift =
                    quadratic
                            ? createQuadraticLift(BigInteger.valueOf(initialModulus), poly, a.toBigPoly(), b.toBigPoly())
                            : createLinearLift(initialModulus, poly, a, b);
            lift.lift(quadratic ? nQuadraticIterations : nLinearIterations);
            return new bMutablePolynomialMod[]{lift.aFactorMod(), lift.bFactorMod()};
        }
    }

    private static <T extends IMutablePolynomialZp<T>> void assertHenselLift(LiftableQuintet<T> lift) {
        assert lift.polyMod().equals(lift.aFactorMod().clone().multiply(lift.bFactorMod())) : lift.toString();
        assert (lift.aCoFactorMod() == null && lift.bCoFactorMod() == null) ||
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod()))
                        .isOne() :
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod())) + "  --- " + ((bMutablePolynomialMod) lift.aFactorMod()).modulus;
    }

    /* ************************************ Quadratic lifts ************************************ */

    /** data used in Hensel lifting **/
    public static abstract class QuadraticLiftAbstract<PolyZp extends IMutablePolynomialZp<PolyZp>>
            implements LiftableQuintet<PolyZp> {
        /** two factors of the initial Z[x] poly **/
        public PolyZp aFactor, bFactor;
        /** xgcd coefficients **/
        public PolyZp aCoFactor, bCoFactor;

        public QuadraticLiftAbstract(PolyZp aFactor, PolyZp bFactor, PolyZp aCoFactor, PolyZp bCoFactor) {
            this.aFactor = aFactor;
            this.bFactor = bFactor;
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
        }

        @Override
        public PolyZp aFactorMod() {return aFactor;}

        @Override
        public PolyZp bFactorMod() {return bFactor;}

        @Override
        public PolyZp aCoFactorMod() {return aCoFactor;}

        @Override
        public PolyZp bCoFactorMod() {return bCoFactor;}

        abstract void prepare();

        @Override
        public final void lift() {
            prepare();
            henselStep0(polyMod());
        }

        @Override
        public final void liftLast() {
            prepare();
            henselLastStep0(polyMod());
        }

        private void henselStep0(PolyZp baseMod) {
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

            PolyZp aCoFactorNew = aCoFactor.subtract(d);
            PolyZp bCoFactorNew = bCoFactor
                    .subtract(bCoFactor.clone().multiply(b))
                    .subtract(c.clone().multiply(aFactorNew));

            assert aFactorNew.degree() == aFactor.degree() : String.format("%s > %s", aFactorNew.degree(), aFactor.degree());
            assert bFactorNew.degree() == bFactor.degree() : String.format("%s > %s", bFactorNew.degree(), bFactor.degree());

            aFactor = aFactorNew; aCoFactor = aCoFactorNew;
            bFactor = bFactorNew; bCoFactor = bCoFactorNew;

            assert bFactor.isMonic();
            assert aCoFactor.degree() < bFactor.degree();
            assert bCoFactor.degree() < aFactor.degree();
        }

        private void henselLastStep0(PolyZp baseMod) {
            PolyZp e = baseMod.subtract(aFactor.clone().multiply(bFactor));

            PolyZp[] qr = divideAndRemainder(
                    aCoFactor.multiply(e),
                    bFactor, false);
            PolyZp q = qr[0], r = qr[1];

            PolyZp aFactorNew = aFactor
                    .add(bCoFactor.multiply(e))
                    .add(aFactor.clone().multiply(q));

            PolyZp bFactorNew = bFactor.add(r);

            aFactor = aFactorNew; aCoFactor = null;
            bFactor = bFactorNew; bCoFactor = null;

            assert bFactor.isMonic();
        }

        @Override
        public String toString() {
            return new StringBuilder()
                    .append("aFactor   = " + aFactor + "\n")
                    .append("bFactor   = " + bFactor + "\n")
                    .append("aCoFactor = " + aCoFactor + "\n")
                    .append("bCoFactor = " + bCoFactor + "\n")
                    .toString();
        }
    }

    /** lift for machine-size polys **/
    public static final class lQuadraticLift extends QuadraticLiftAbstract<MutablePolynomialMod> {
        /** the modulus */
        public long modulus;
        /** initial Z[x] poly **/
        public final MutablePolynomialZ base;

        public lQuadraticLift(long modulus, MutablePolynomialZ base, MutablePolynomialMod aFactor, MutablePolynomialMod bFactor, MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        public MutablePolynomialMod polyMod() {
            return base.modulus(modulus, true);
        }

        @Override
        void prepare() {
            modulus = safeMultiply(modulus, modulus);
            aFactor = aFactor.setModulusUnsafe(modulus);
            bFactor = bFactor.setModulusUnsafe(modulus);
            aCoFactor = aCoFactor.setModulusUnsafe(modulus);
            bCoFactor = bCoFactor.setModulusUnsafe(modulus);
        }

        @Override
        public String toString() {
            return new StringBuilder()
                    .append(super.toString())
                    .append("base      = " + base + "\n")
                    .append("modulus   = " + modulus)
                    .toString();
        }
    }

    /** lift for BigInteger polys **/
    public static final class bQuadraticLift extends QuadraticLiftAbstract<bMutablePolynomialMod> {
        /** the modulus */
        public BigInteger modulus;
        /** initial Z[x] poly **/
        public final bMutablePolynomialZ base;

        public bQuadraticLift(BigInteger modulus, bMutablePolynomialZ base, bMutablePolynomialMod aFactor, bMutablePolynomialMod bFactor, bMutablePolynomialMod aCoFactor, bMutablePolynomialMod bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        public bMutablePolynomialMod polyMod() {
            return base.modulus(modulus, true);
        }

        @Override
        void prepare() {
            modulus = modulus.multiply(modulus);
            aFactor = aFactor.setModulusUnsafe(modulus);
            bFactor = bFactor.setModulusUnsafe(modulus);
            aCoFactor = aCoFactor.setModulusUnsafe(modulus);
            bCoFactor = bCoFactor.setModulusUnsafe(modulus);
        }

        @Override
        public String toString() {
            return new StringBuilder()
                    .append(super.toString())
                    .append("base      = " + base + "\n")
                    .append("modulus   = " + modulus)
                    .toString();
        }
    }

    /* ************************************ Linear lifts ************************************ */

    private static class LinearLiftAbstract<PolyZ extends IMutablePolynomialZ<PolyZ>> {
        final PolyZ poly;
        PolyZ aFactor, bFactor, aCoFactor, bCoFactor;
        final MutablePolynomialMod aFactorMod, aFactorModMonic, bFactorMod, aCoFactorMod, bCoFactorMod;
        final InverseModMonomial<MutablePolynomialMod> aFactorModMonicInv, bFactorModInv;

        public LinearLiftAbstract(PolyZ poly,
                                  PolyZ aFactor, PolyZ bFactor, PolyZ aCoFactor, PolyZ bCoFactor,
                                  MutablePolynomialMod aFactorMod, MutablePolynomialMod aFactorModMonic, MutablePolynomialMod bFactorMod,
                                  MutablePolynomialMod aCoFactorMod, MutablePolynomialMod bCoFactorMod) {
            assert bFactor.isMonic();
            this.poly = poly;
            this.aFactor = aFactor;
            this.bFactor = bFactor;
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
            this.aFactorMod = aFactorMod;
            this.aFactorModMonic = aFactorModMonic;
            this.bFactorMod = bFactorMod;
            this.aCoFactorMod = aCoFactorMod;
            this.bCoFactorMod = bCoFactorMod;
            this.aFactorModMonicInv = fastDivisionPreConditioning(aFactorModMonic);
            this.bFactorModInv = fastDivisionPreConditioning(bFactorMod);
        }

        protected MutablePolynomialMod aAdd, bAdd;

        final void calculateFactorsDiff(MutablePolynomialMod diff) {
            aAdd = diff.clone();
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);
            aAdd = aAdd.multiply(bCoFactorMod);
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);

            bAdd = diff.clone();
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
            bAdd = bAdd.multiply(aCoFactorMod);
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
        }

        final void calculateCoFactorsDiff(MutablePolynomialMod diff) {
            aAdd = diff.clone();
            aAdd = PolynomialArithmetics.polyMod(aAdd, bFactorMod, bFactorModInv, false);
            aAdd = aAdd.multiply(aCoFactorMod);
            aAdd = PolynomialArithmetics.polyMod(aAdd, bFactorMod, bFactorModInv, false);

            bAdd = diff.clone();
            bAdd = PolynomialArithmetics.polyMod(bAdd, aFactorModMonic, aFactorModMonicInv, false);
            bAdd = bAdd.multiply(bCoFactorMod);
            bAdd = PolynomialArithmetics.polyMod(bAdd, aFactorModMonic, aFactorModMonicInv, false);
        }
    }

    public static final class lLinearLift
            extends LinearLiftAbstract<MutablePolynomialZ>
            implements LiftableQuintet<MutablePolynomialMod> {
        final long initialModulus;
        public long modulus;

        private lLinearLift(long modulus, MutablePolynomialZ poly,
                            MutablePolynomialMod aFactor, MutablePolynomialMod aFactorMonic, MutablePolynomialMod bFactor,
                            MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).normalForm(true).multiply(poly.lc()), bFactor.normalForm(true),
                    aCoFactor.normalForm(true), bCoFactor.normalForm(true),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialModulus = modulus;
            this.modulus = modulus;
        }

        private lLinearLift(long modulus, MutablePolynomialZ poly,
                            MutablePolynomialMod aFactor, MutablePolynomialMod bFactor,
                            MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public MutablePolynomialMod polyMod() {return poly.modulus(modulus);}

        @Override
        public MutablePolynomialMod aFactorMod() {return aFactor.modulus(modulus);}

        @Override
        public MutablePolynomialMod bFactorMod() {return bFactor.modulus(modulus);}

        @Override
        public MutablePolynomialMod aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.modulus(modulus);}

        @Override
        public MutablePolynomialMod bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.modulus(modulus);}

        private void liftFactors() {
            MutablePolynomialMod factorsDiff = poly.clone().subtract(aFactor.clone().multiply(bFactor))
                    .divideOrNull(modulus)
                    .modulus(initialModulus);

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.normalForm(false).multiply(modulus));
            bFactor = bFactor.add(bAdd.normalForm(false).multiply(modulus));
        }

        private void liftCoFactors() {
            MutablePolynomialMod coFactorsDiff = aCoFactor.clone().multiply(aFactor)
                    .add(bCoFactor.clone().multiply(bFactor))
                    .decrement()
                    .negate()
                    .divideOrNull(modulus)
                    .modulus(initialModulus);

            calculateCoFactorsDiff(coFactorsDiff);
            aCoFactor = aCoFactor.add(aAdd.normalForm(false).multiply(modulus));
            bCoFactor = bCoFactor.add(bAdd.normalForm(false).multiply(modulus));
        }

        @Override
        public void lift() {
            liftFactors();
            liftCoFactors();
            modulus = LongArithmetics.safeMultiply(modulus, initialModulus);
        }

        @Override
        public void liftLast() {
            liftFactors();
            modulus = LongArithmetics.safeMultiply(modulus, initialModulus);
            aCoFactor = bCoFactor = null;
        }
    }

    public static final class bLinearLift
            extends LinearLiftAbstract<bMutablePolynomialZ>
            implements LiftableQuintet<bMutablePolynomialMod> {
        final BigInteger initialModulus;
        public BigInteger modulus;

        public bLinearLift(bMutablePolynomialZ poly,
                           bMutablePolynomialZ aFactor,
                           bMutablePolynomialZ bFactor,
                           bMutablePolynomialZ aCoFactor,
                           bMutablePolynomialZ bCoFactor,
                           MutablePolynomialMod aFactorMod,
                           MutablePolynomialMod aFactorModMonic,
                           MutablePolynomialMod bFactorMod,
                           MutablePolynomialMod aCoFactorMod,
                           MutablePolynomialMod bCoFactorMod,
                           BigInteger initialModulus,
                           BigInteger modulus) {
            super(poly, aFactor, bFactor, aCoFactor, bCoFactor, aFactorMod, aFactorModMonic, bFactorMod, aCoFactorMod, bCoFactorMod);
            this.initialModulus = initialModulus;
            this.modulus = modulus;
        }

        private bLinearLift(BigInteger modulus, bMutablePolynomialZ poly,
                            MutablePolynomialMod aFactor, MutablePolynomialMod aFactorMonic, MutablePolynomialMod bFactor,
                            MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).normalForm(true).toBigPoly().multiply(poly.lc()), bFactor.normalForm(false).toBigPoly(),
                    aCoFactor.normalForm(false).toBigPoly(), bCoFactor.normalForm(false).toBigPoly(),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialModulus = modulus;
            this.modulus = modulus;
            assert modulus.isLong();
        }

        private bLinearLift(BigInteger modulus, bMutablePolynomialZ poly,
                            MutablePolynomialMod aFactor, MutablePolynomialMod bFactor,
                            MutablePolynomialMod aCoFactor, MutablePolynomialMod bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public bMutablePolynomialMod polyMod() {return poly.modulus(modulus);}

        @Override
        public bMutablePolynomialMod aFactorMod() {return aFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialMod bFactorMod() {return bFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialMod aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialMod bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.modulus(modulus);}

        private void liftFactors() {
            MutablePolynomialMod factorsDiff = poly.clone().subtract(aFactor.clone().multiply(bFactor))
                    .divideOrNull(modulus)
                    .modulus(initialModulus)
                    .toLong();

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.normalForm(false).toBigPoly().multiply(modulus));
            bFactor = bFactor.add(bAdd.normalForm(false).toBigPoly().multiply(modulus));
        }

        private void liftCoFactors() {
            MutablePolynomialMod coFactorsDiff = aCoFactor.clone().multiply(aFactor).add(bCoFactor.clone().multiply(bFactor)).decrement().negate()
                    .divideOrNull(modulus)
                    .modulus(initialModulus)
                    .toLong();

            calculateCoFactorsDiff(coFactorsDiff);

            aCoFactor = aCoFactor.add(aAdd.normalForm(false).toBigPoly().multiply(modulus));
            bCoFactor = bCoFactor.add(bAdd.normalForm(false).toBigPoly().multiply(modulus));
        }

        @Override
        public void lift() {
            liftFactors();
            liftCoFactors();
            modulus = modulus.multiply(initialModulus);
        }

        @Override
        public void liftLast() {
            liftFactors();
            modulus = modulus.multiply(initialModulus);
            aCoFactor = bCoFactor = null;
        }
    }

    public interface LiftableQuintet<PolyZp extends IMutablePolynomialZp<PolyZp>> {
        PolyZp polyMod();

        PolyZp aFactorMod();

        PolyZp bFactorMod();

        PolyZp aCoFactorMod();

        PolyZp bCoFactorMod();

        void lift();

        void liftLast();

        default void lift(int nIterations) {
            for (int i = 0; i < nIterations - 1; ++i)
                lift();
            liftLast();
        }

        default void liftWithCoFactors(int nIterations) {
            for (int i = 0; i < nIterations; ++i)
                lift();
        }
    }
}
