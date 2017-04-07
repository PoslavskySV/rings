package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static cc.r2.core.poly.univar.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly.univar.LongArithmetics.safeMultiply;

/**
 * Methods for univariate Hensel lifting.
 * <p>
 * <i>Implementation notes.</i>
 * Two methods for Hensel lift are implemented: quadratic and linear. For {@code N} iterations quadratic lift will
 * lift to p<sup>2^N</sup> while linear just to p<sup>N</sup>. While quadratic lift converges much
 * faster, it works with BigIntegers in all intermediate steps, so each step is quite expensive. Linear lift is
 * implemented so that it starts with machine-word modulus, and perform all hard intermediate calculations with
 * machine-word arithmetics, converting to BigIntegers only a few times. In this way, a single step of linear lift is
 * very cheap, but the convergence is worse. The actual lifting used in factorization switches between linear
 * and quadratic lift in order to obtain the best trade-off.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class HenselLifting {
    private HenselLifting() {}

    /**
     * Liftable quintet. Output specifications is as follows:
     * <p>
     * <pre>
     * polyMod = aFactor * bFactor mod modulus
     * 1 = aFactor * aCoFactor + bFactor * bCoFactor mod modulus
     * </pre>
     * where {@coode modulus} is the modulus obtained by lifting
     *
     * @param <PolyZp> Zp[x] polynomial type
     */
    public interface LiftableQuintet<PolyZp extends IMutablePolynomialZp<PolyZp>> {
        /**
         * Returns initial Z[x] polynomial modulo lifted modulus
         *
         * @return initial Z[x] polynomial modulo lifted modulus
         */
        PolyZp polyMod();

        /**
         * Returns first factor lifted
         *
         * @return first factor lifted
         */
        PolyZp aFactorMod();

        /**
         * Returns second factor lifted
         *
         * @return second factor lifted
         */
        PolyZp bFactorMod();

        /**
         * Returns first co-factor lifted
         *
         * @return first co-factor lifted
         */
        PolyZp aCoFactorMod();

        /**
         * Returns second co-factor lifted
         *
         * @return second co-factor lifted
         */
        PolyZp bCoFactorMod();

        /**
         * Performs single lift step.
         */
        void lift();

        /**
         * Performs single lift step but don't lift co-factors (xgcd coefficients).
         */
        void liftLast();

        /**
         * Lifts {@code nIterations} times. Co-factor will be lost on the last step.
         *
         * @param nIterations number of lift iterations
         */
        default void lift(int nIterations) {
            for (int i = 0; i < nIterations - 1; ++i)
                lift();
            liftLast();
        }

        /**
         * Lifts {@code nIterations} times.
         *
         * @param nIterations number of lift iterations
         */
        default void liftWithCoFactors(int nIterations) {
            for (int i = 0; i < nIterations; ++i)
                lift();
        }
    }

    /* ************************************ Factory methods ************************************ */

    /**
     * Creates quadratic Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return quadratic Hensel lift
     */
    public static lQuadraticLift createQuadraticLift(long modulus,
                                                     lMutablePolynomialZ poly,
                                                     lMutablePolynomialZp aFactor,
                                                     lMutablePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.mod(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());
        lMutablePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lQuadraticLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /**
     * Creates quadratic Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return quadratic Hensel lift
     */
    public static bQuadraticLift createQuadraticLift(BigInteger modulus,
                                                     bMutablePolynomialZ poly,
                                                     bMutablePolynomialZp aFactor,
                                                     bMutablePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        if (!aFactor.mod(poly.lc()).equals(aFactor.lc()))
            aFactor = aFactor.clone().monic(poly.lc());
        bMutablePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bQuadraticLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /**
     * Creates quadratic Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return quadratic Hensel lift
     */
    public static bQuadraticLift createQuadraticLift(BigInteger modulus,
                                                     bMutablePolynomialZ poly,
                                                     lMutablePolynomialZp aFactor,
                                                     lMutablePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(modulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);
        lMutablePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bQuadraticLift(modulus, poly, aFactor.toBigPoly(), bFactor.toBigPoly(), xgcd[1].toBigPoly(), xgcd[2].toBigPoly());
    }

    /**
     * Creates linear Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return linear Hensel lift
     */
    public static lLinearLift createLinearLift(BigInteger modulus,
                                               lMutablePolynomialZ poly,
                                               lMutablePolynomialZp aFactor,
                                               lMutablePolynomialZp bFactor) {
        return createLinearLift(modulus.longValueExact(), poly, aFactor, bFactor);
    }

    /**
     * Creates linear Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return linear Hensel lift
     */
    public static bLinearLift createLinearLift(BigInteger modulus,
                                               bMutablePolynomialZ poly,
                                               lMutablePolynomialZp aFactor,
                                               lMutablePolynomialZp bFactor) {
        return createLinearLift(modulus.longValueExact(), poly, aFactor, bFactor);
    }

    /**
     * Creates linear Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return linear Hensel lift
     */
    public static lLinearLift createLinearLift(long modulus,
                                               lMutablePolynomialZ poly,
                                               lMutablePolynomialZp aFactor,
                                               lMutablePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.mod(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());

        lMutablePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lLinearLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /**
     * Creates linear Hensel lift.
     *
     * @param modulus the initial modulus
     * @param poly    Z[x] polynomial
     * @param aFactor first factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @param bFactor second factor of {@code poly} that {@code poly = aFactor * bFactor mod modulus}
     * @return linear Hensel lift
     */
    public static bLinearLift createLinearLift(long modulus,
                                               bMutablePolynomialZ poly,
                                               lMutablePolynomialZp aFactor,
                                               lMutablePolynomialZp bFactor) {
        BigInteger bModulus = BigInteger.valueOf(modulus);
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(bModulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);

        lMutablePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bLinearLift(bModulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static <PolyZp extends IMutablePolynomialZp<PolyZp>> PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
        PolyZp[] xgcd = PolynomialGCD.ExtendedEuclid(a, b);
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

    /**
     * Lifts modular factorization until {@code modulus} will overcome {@code desiredBound}.
     *
     * @param modulus        initial modulus so that {@code modularFactors} are true factors of {@code poly mod modulus}
     * @param desiredBound   desired modulus
     * @param poly           initial Z[x] polynomial
     * @param modularFactors factorization of {@code poly.modulus(modulus)}
     * @param quadratic      whether to use quadratic of linear lift
     * @return factorization of {@code poly.modulus(finalModulus) } with some {@code finalModulus} greater than {@code desiredBound}
     */
    public static List<lMutablePolynomialZp> liftFactorization(long modulus,
                                                               long desiredBound,
                                                               lMutablePolynomialZ poly,
                                                               List<lMutablePolynomialZp> modularFactors,
                                                               boolean quadratic) {
        long[] im = nIterations(modulus, desiredBound, quadratic);
        return liftFactorization(modulus, im[1], (int) im[0], poly, modularFactors, quadratic);
    }

    /**
     * Lifts modular factorization {@code nIterations} times using whether linear or quadratic lifting.
     *
     * @param modulus        initial modulus so that {@code modularFactors} are true factors of {@code poly mod modulus}
     * @param finalModulus   final modulus that will be obtained after lifting
     * @param nIterations    number of lifting steps to do
     * @param poly           initial Z[x] polynomial
     * @param modularFactors factorization of {@code poly.modulus(modulus)}
     * @param quadratic      whether to use quadratic of linear lift
     * @return factorization of {@code poly.modulus(finalModulus) }
     */
    public static List<lMutablePolynomialZp> liftFactorization(long modulus,
                                                               long finalModulus,
                                                               int nIterations,
                                                               lMutablePolynomialZ poly,
                                                               List<lMutablePolynomialZp> modularFactors,
                                                               boolean quadratic) {
        assert nIterations > 0;

        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(finalModulus, true).monic());

        lMutablePolynomialZp factory = modularFactors.get(0);
        lMutablePolynomialZp
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        LiftableQuintet<lMutablePolynomialZp> hensel = quadratic
                ? createQuadraticLift(modulus, poly, aFactor, bFactor)
                : createLinearLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        lMutablePolynomialZp aFactorRaised = hensel.aFactorMod();
        lMutablePolynomialZp bFactorRaised = hensel.bFactorMod();

        ArrayList<lMutablePolynomialZp> result = new ArrayList<>();
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, aFactorRaised.normalSymmetricForm(), modularFactors.subList(0, nHalf), quadratic));
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, bFactorRaised.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size()), quadratic));
        return result;
    }

    interface LiftFactory<PolyZp extends IMutablePolynomialZp<PolyZp>> {
        LiftableQuintet<bMutablePolynomialZp> createLift(BigInteger modulus, bMutablePolynomialZ polyZ, PolyZp aFactor, PolyZp bFactor);
    }

    /** actual multifactor Hensel lifting implementation **/
    static <PolyZp extends IMutablePolynomialZp<PolyZp>>
    List<bMutablePolynomialZp> liftFactorization0(BigInteger modulus,
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

        LiftableQuintet<bMutablePolynomialZp> hensel = liftFactory.createLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        bMutablePolynomialZp aFactorRaised = hensel.aFactorMod();
        bMutablePolynomialZp bFactorRaised = hensel.bFactorMod();

        ArrayList<bMutablePolynomialZp> result = new ArrayList<>();
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

    /**
     * Lifts modular factorization until {@code modulus} will overcome {@code desiredBound}. <b>Note:</b> if
     * {@code quadratic == false} modulus must fit 64-bit.
     *
     * @param modulus        initial modulus so that {@code modularFactors} are true factors of {@code poly mod modulus}
     * @param desiredBound   desired modulus
     * @param poly           initial Z[x] polynomial
     * @param modularFactors factorization of {@code poly.modulus(modulus)}
     * @param quadratic      whether to use quadratic of linear lift
     * @return factorization of {@code poly.modulus(finalModulus) } with some {@code finalModulus} greater than {@code desiredBound}
     */
    public static List<bMutablePolynomialZp> liftFactorization(BigInteger modulus,
                                                               BigInteger desiredBound,
                                                               bMutablePolynomialZ poly,
                                                               List<lMutablePolynomialZp> modularFactors,
                                                               boolean quadratic) {
        if (!quadratic && !modulus.isLong())
            throw new IllegalArgumentException("Only max 64-bit modulus for linear lift allowed.");
        LiftingInfo im = nIterations(modulus, desiredBound, quadratic);
        if (im.nIterations == 0)
            return modularFactors.stream().map(lMutablePolynomialZp::toBigPoly).collect(Collectors.toList());
        LiftFactory<lMutablePolynomialZp> factory = quadratic ? HenselLifting::createQuadraticLift : HenselLifting::createLinearLift;
        return liftFactorization0(modulus, im.finalModulus, im.nIterations, poly, modularFactors, factory);
    }

    /**
     * Lifts modular factorization until {@code modulus} will overcome {@code desiredBound}.
     *
     * @param modulus        initial modulus so that {@code modularFactors} are true factors of {@code poly mod modulus}
     * @param desiredBound   desired modulus
     * @param poly           initial Z[x] polynomial
     * @param modularFactors factorization of {@code poly.modulus(modulus)}
     * @return factorization of {@code poly.modulus(finalModulus) } with some {@code finalModulus} greater than {@code desiredBound}
     */
    public static List<bMutablePolynomialZp> liftFactorizationQuadratic(BigInteger modulus,
                                                                        BigInteger desiredBound,
                                                                        bMutablePolynomialZ poly,
                                                                        List<bMutablePolynomialZp> modularFactors) {
        LiftingInfo im = nIterations(modulus, desiredBound, true);
        return liftFactorization0(modulus, im.finalModulus, im.nIterations, poly, modularFactors, HenselLifting::createQuadraticLift);
    }

    /**
     * Lifts modular factorization until {@code modulus} will overcome {@code desiredBound}. <i>Implementation note:</i>
     * method will switch between linear and quadratic lift depending on the required lifting iterations.
     *
     * @param modulus        initial modulus so that {@code modularFactors} are true factors of {@code poly mod modulus}
     * @param desiredBound   desired modulus
     * @param poly           initial Z[x] polynomial
     * @param modularFactors factorization of {@code poly.modulus(modulus)}
     * @return factorization of {@code poly.modulus(finalModulus) } with some {@code finalModulus} greater than {@code desiredBound}
     */
    public static List<bMutablePolynomialZp> liftFactorization(BigInteger modulus,
                                                               BigInteger desiredBound,
                                                               bMutablePolynomialZ poly,
                                                               List<lMutablePolynomialZp> modularFactors) {
        return liftFactorization(poly, modularFactors, new AdaptiveLift(modulus, desiredBound));
    }

    /** actual multifactor Hensel lifting implementation **/
    private static List<bMutablePolynomialZp> liftFactorization(bMutablePolynomialZ poly,
                                                                List<lMutablePolynomialZp> modularFactors,
                                                                AdaptiveLift lifter) {
        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(lifter.finalModulus, true).monic());

        lMutablePolynomialZp factory = modularFactors.get(0);
        lMutablePolynomialZp
                aFactor = factory.createOne(),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        bMutablePolynomialZp[] lifted = lifter.lift(poly, aFactor, bFactor);
        bMutablePolynomialZp aFactorRaised = lifted[0];
        bMutablePolynomialZp bFactorRaised = lifted[1];

        ArrayList<bMutablePolynomialZp> result = new ArrayList<>();
        result.addAll(liftFactorization(aFactorRaised.normalSymmetricForm(), modularFactors.subList(0, nHalf), lifter));
        result.addAll(liftFactorization(bFactorRaised.normalSymmetricForm(), modularFactors.subList(nHalf, modularFactors.size()), lifter));
        return result;
    }

    private static final int SWITCH_TO_QUADRATIC_LIFT = 64;

    private static final class AdaptiveLift {
        final BigInteger initialModulus;
        final BigInteger finalModulus;
        final int nLinearIterations, nQuadraticIterations;

        public AdaptiveLift(BigInteger initialModulus, BigInteger desiredBound) {
            this.initialModulus = initialModulus;
            LiftingInfo
                    nLinearIterations = nIterations(initialModulus, desiredBound, false);

            if (nLinearIterations.nIterations < SWITCH_TO_QUADRATIC_LIFT) {
                this.nLinearIterations = nLinearIterations.nIterations;
                this.nQuadraticIterations = -1;
                this.finalModulus = nLinearIterations.finalModulus;
            } else {
                LiftingInfo nQuadraticIterations = nIterations(initialModulus, desiredBound, true);
                this.nLinearIterations = -1;
                this.nQuadraticIterations = nQuadraticIterations.nIterations;
                this.finalModulus = nQuadraticIterations.finalModulus;
            }
        }

        bMutablePolynomialZp[] lift(bMutablePolynomialZ poly, lMutablePolynomialZp a, lMutablePolynomialZp b) {
            boolean quadratic = nLinearIterations == -1;
            LiftableQuintet<bMutablePolynomialZp> lift =
                    quadratic
                            ? createQuadraticLift(initialModulus, poly, a.toBigPoly(), b.toBigPoly())
                            : createLinearLift(initialModulus, poly, a, b);
            lift.lift(quadratic ? nQuadraticIterations : nLinearIterations);
            return new bMutablePolynomialZp[]{lift.aFactorMod(), lift.bFactorMod()};
        }
    }

    private static <T extends IMutablePolynomialZp<T>> void assertHenselLift(LiftableQuintet<T> lift) {
        assert lift.polyMod().equals(lift.aFactorMod().clone().multiply(lift.bFactorMod())) : lift.toString();
        assert (lift.aCoFactorMod() == null && lift.bCoFactorMod() == null) ||
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod()))
                        .isOne() :
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod())) + "  --- " + ((bMutablePolynomialZp) lift.aFactorMod()).modulus;
    }

    /* ************************************ Quadratic lifts ************************************ */

    /** data used in Hensel lifting **/
    static abstract class QuadraticLiftAbstract<PolyZp extends IMutablePolynomialZp<PolyZp>>
            implements LiftableQuintet<PolyZp> {
        /** Two factors of the initial Z[x] poly **/
        protected PolyZp aFactor, bFactor;
        /** xgcd coefficients **/
        protected PolyZp aCoFactor, bCoFactor;

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

            PolyZp[] qr = DivisionWithRemainder.divideAndRemainder(
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

            PolyZp[] cd = DivisionWithRemainder.divideAndRemainder(
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

            PolyZp[] qr = DivisionWithRemainder.divideAndRemainder(
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
    }

    /**
     * Quadratic Hensel lift for machine word arithmetics. On each {@link #lift()} operation modulus is raised
     * as {@code modulus = modulus * modulus}.
     */
    public static final class lQuadraticLift extends QuadraticLiftAbstract<lMutablePolynomialZp> {
        /** The modulus */
        public long modulus;
        /** Initial Z[x] poly **/
        public final lMutablePolynomialZ base;

        public lQuadraticLift(long modulus, lMutablePolynomialZ base, lMutablePolynomialZp aFactor, lMutablePolynomialZp bFactor, lMutablePolynomialZp aCoFactor, lMutablePolynomialZp bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        public lMutablePolynomialZp polyMod() {
            return base.modulus(modulus, true);
        }

        @Override
        void prepare() {
            modulus = LongArithmetics.safeMultiply(modulus, modulus);
            aFactor = aFactor.setModulusUnsafe(modulus);
            bFactor = bFactor.setModulusUnsafe(modulus);
            aCoFactor = aCoFactor.setModulusUnsafe(modulus);
            bCoFactor = bCoFactor.setModulusUnsafe(modulus);
        }
    }

    /**
     * Quadratic Hensel lift for BigIntegers arithmetics. On each {@link #lift()} operation modulus is raised
     * as {@code modulus = modulus * modulus}.
     */
    public static final class bQuadraticLift extends QuadraticLiftAbstract<bMutablePolynomialZp> {
        /** The modulus */
        public BigInteger modulus;
        /** Initial Z[x] poly **/
        public final bMutablePolynomialZ base;

        public bQuadraticLift(BigInteger modulus, bMutablePolynomialZ base, bMutablePolynomialZp aFactor, bMutablePolynomialZp bFactor, bMutablePolynomialZp aCoFactor, bMutablePolynomialZp bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        public bMutablePolynomialZp polyMod() {
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
    }

    /* ************************************ Linear lifts ************************************ */

    private static class LinearLiftAbstract<PolyZ extends IMutablePolynomialZ<PolyZ>> {
        /** initial Z[x] poly */
        final PolyZ poly;
        /** lifted polynomials */
        PolyZ aFactor, bFactor, aCoFactor, bCoFactor;
        /** initial modular data */
        final lMutablePolynomialZp aFactorMod, aFactorModMonic, bFactorMod, aCoFactorMod, bCoFactorMod;
        /** precomputed inverses */
        final DivisionWithRemainder.InverseModMonomial<lMutablePolynomialZp> aFactorModMonicInv, bFactorModInv;

        public LinearLiftAbstract(PolyZ poly,
                                  PolyZ aFactor, PolyZ bFactor, PolyZ aCoFactor, PolyZ bCoFactor,
                                  lMutablePolynomialZp aFactorMod, lMutablePolynomialZp aFactorModMonic, lMutablePolynomialZp bFactorMod,
                                  lMutablePolynomialZp aCoFactorMod, lMutablePolynomialZp bCoFactorMod) {
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
            this.aFactorModMonicInv = DivisionWithRemainder.fastDivisionPreConditioning(aFactorModMonic);
            this.bFactorModInv = DivisionWithRemainder.fastDivisionPreConditioning(bFactorMod);
        }

        protected lMutablePolynomialZp aAdd, bAdd;

        final void calculateFactorsDiff(lMutablePolynomialZp diff) {
            aAdd = diff.clone();
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);
            aAdd = aAdd.multiply(bCoFactorMod);
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);

            bAdd = diff.clone();
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
            bAdd = bAdd.multiply(aCoFactorMod);
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
        }

        final void calculateCoFactorsDiff(lMutablePolynomialZp diff) {
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

    /**
     * Linear Hensel lift for machine word arithmetics. Linear Hensel lift always starts from the machine-sized modulus;
     * on each {@link #lift()} operation modulus is raised as {@code modulus = modulus * initialModulus}.
     */
    public static final class lLinearLift
            extends LinearLiftAbstract<lMutablePolynomialZ>
            implements LiftableQuintet<lMutablePolynomialZp> {
        /** The initial modulus */
        public final long initialModulus;
        /** The modulus */
        public long modulus;

        private lLinearLift(long modulus, lMutablePolynomialZ poly,
                            lMutablePolynomialZp aFactor, lMutablePolynomialZp aFactorMonic, lMutablePolynomialZp bFactor,
                            lMutablePolynomialZp aCoFactor, lMutablePolynomialZp bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).normalForm(true).multiply(poly.lc()), bFactor.normalForm(true),
                    aCoFactor.normalForm(true), bCoFactor.normalForm(true),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialModulus = modulus;
            this.modulus = modulus;
        }

        private lLinearLift(long modulus, lMutablePolynomialZ poly,
                            lMutablePolynomialZp aFactor, lMutablePolynomialZp bFactor,
                            lMutablePolynomialZp aCoFactor, lMutablePolynomialZp bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public lMutablePolynomialZp polyMod() {return poly.modulus(modulus);}

        @Override
        public lMutablePolynomialZp aFactorMod() {return aFactor.modulus(modulus);}

        @Override
        public lMutablePolynomialZp bFactorMod() {return bFactor.modulus(modulus);}

        @Override
        public lMutablePolynomialZp aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.modulus(modulus);}

        @Override
        public lMutablePolynomialZp bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.modulus(modulus);}

        private void liftFactors() {
            lMutablePolynomialZp factorsDiff = poly.clone().subtract(aFactor.clone().multiply(bFactor))
                    .divideOrNull(modulus)
                    .modulus(initialModulus);

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.normalForm(false).multiply(modulus));
            bFactor = bFactor.add(bAdd.normalForm(false).multiply(modulus));
        }

        private void liftCoFactors() {
            lMutablePolynomialZp coFactorsDiff = aCoFactor.clone().multiply(aFactor)
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

    /**
     * Linear Hensel lift for BigIntegers arithmetics. Linear Hensel lift always starts from the machine-sized modulus;
     * on each {@link #lift()} operation modulus is raised as {@code modulus = modulus * initialModulus}.
     */
    public static final class bLinearLift
            extends LinearLiftAbstract<bMutablePolynomialZ>
            implements LiftableQuintet<bMutablePolynomialZp> {
        /** The initial modulus (less than 64-bit) */
        public final BigInteger initialModulus;
        /** The modulus */
        public BigInteger modulus;

        private bLinearLift(BigInteger modulus, bMutablePolynomialZ poly,
                            lMutablePolynomialZp aFactor, lMutablePolynomialZp aFactorMonic, lMutablePolynomialZp bFactor,
                            lMutablePolynomialZp aCoFactor, lMutablePolynomialZp bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).normalForm(true).toBigPoly().multiply(poly.lc()), bFactor.normalForm(false).toBigPoly(),
                    aCoFactor.normalForm(false).toBigPoly(), bCoFactor.normalForm(false).toBigPoly(),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialModulus = modulus;
            this.modulus = modulus;
            assert modulus.isLong();
        }

        private bLinearLift(BigInteger modulus, bMutablePolynomialZ poly,
                            lMutablePolynomialZp aFactor, lMutablePolynomialZp bFactor,
                            lMutablePolynomialZp aCoFactor, lMutablePolynomialZp bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public bMutablePolynomialZp polyMod() {return poly.modulus(modulus);}

        @Override
        public bMutablePolynomialZp aFactorMod() {return aFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialZp bFactorMod() {return bFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialZp aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.modulus(modulus);}

        @Override
        public bMutablePolynomialZp bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.modulus(modulus);}

        private void liftFactors() {
            lMutablePolynomialZp factorsDiff = poly.clone().subtract(aFactor.clone().multiply(bFactor))
                    .divideOrNull(modulus)
                    .modulus(initialModulus)
                    .toLong();

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.normalForm(false).toBigPoly().multiply(modulus));
            bFactor = bFactor.add(bAdd.normalForm(false).toBigPoly().multiply(modulus));
        }

        private void liftCoFactors() {
            lMutablePolynomialZp coFactorsDiff = aCoFactor.clone().multiply(aFactor).add(bCoFactor.clone().multiply(bFactor)).decrement().negate()
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
}
