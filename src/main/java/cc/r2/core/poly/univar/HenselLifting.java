package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Integers;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.LongArithmetics;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

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
    public interface LiftableQuintet<PolyZp extends IUnivariatePolynomial<PolyZp>> {
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
                                                     lUnivariatePolynomialZ poly,
                                                     lUnivariatePolynomialZp aFactor,
                                                     lUnivariatePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.domain.modulus(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());
        lUnivariatePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new lQuadraticLift(modulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    private static void ensureIntegersDomain(UnivariatePolynomial<BigInteger> poly) {
        if (poly.domain != Integers.Integers)
            throw new IllegalArgumentException("Not an integers domain domain");
    }

    private static void ensureModularDomain(UnivariatePolynomial<BigInteger> poly) {
        if (!(poly.domain instanceof IntegersModulo))
            throw new IllegalArgumentException("Not a modular domain");
    }

    private static void ensureInputCorrect(UnivariatePolynomial<BigInteger> poly, UnivariatePolynomial<BigInteger> aFactor, UnivariatePolynomial<BigInteger> bFactor) {
        ensureIntegersDomain(poly);
        ensureModularDomain(aFactor);
        ensureModularDomain(bFactor);
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
                                                     UnivariatePolynomial<BigInteger> poly,
                                                     UnivariatePolynomial<BigInteger> aFactor,
                                                     UnivariatePolynomial<BigInteger> bFactor) {
        ensureInputCorrect(poly, aFactor, bFactor);
        bFactor = ensureMonic(bFactor);
        IntegersModulo domain = (IntegersModulo) aFactor.domain;
        if (!domain.valueOf(poly.lc()).equals(aFactor.lc()))
            aFactor = aFactor.clone().monic(poly.lc());
        UnivariatePolynomial<BigInteger>[] xgcd = monicExtendedEuclid(aFactor, bFactor);
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
                                                     UnivariatePolynomial<BigInteger> poly,
                                                     lUnivariatePolynomialZp aFactor,
                                                     lUnivariatePolynomialZp bFactor) {
        ensureIntegersDomain(poly);
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(modulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);
        lUnivariatePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
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
                                               lUnivariatePolynomialZ poly,
                                               lUnivariatePolynomialZp aFactor,
                                               lUnivariatePolynomialZp bFactor) {
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
                                               UnivariatePolynomial<BigInteger> poly,
                                               lUnivariatePolynomialZp aFactor,
                                               lUnivariatePolynomialZp bFactor) {
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
                                               lUnivariatePolynomialZ poly,
                                               lUnivariatePolynomialZp aFactor,
                                               lUnivariatePolynomialZp bFactor) {
        bFactor = ensureMonic(bFactor);
        if (aFactor.domain.modulus(poly.lc()) != aFactor.lc())
            aFactor = aFactor.clone().monic(poly.lc());

        lUnivariatePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
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
                                               UnivariatePolynomial<BigInteger> poly,
                                               lUnivariatePolynomialZp aFactor,
                                               lUnivariatePolynomialZp bFactor) {
        ensureIntegersDomain(poly);
        BigInteger bModulus = BigInteger.valueOf(modulus);
        bFactor = ensureMonic(bFactor);
        long lc = poly.lc().mod(bModulus).longValueExact();
        if (lc != aFactor.lc())
            aFactor = aFactor.clone().monic(lc);

        lUnivariatePolynomialZp[] xgcd = monicExtendedEuclid(aFactor, bFactor);
        return new bLinearLift(bModulus, poly, aFactor, bFactor, xgcd[1], xgcd[2]);
    }

    /** runs xgcd for coprime polynomials ensuring that gcd is 1 (not another constant) */
    private static <PolyZp extends IUnivariatePolynomial<PolyZp>> PolyZp[] monicExtendedEuclid(PolyZp a, PolyZp b) {
        PolyZp[] xgcd = UnivariateGCD.ExtendedEuclidGCD(a, b);
        if (xgcd[0].isOne())
            return xgcd;

        assert xgcd[0].isConstant() : "bad xgcd: " + Arrays.toString(xgcd) + " for xgcd(" + a + ", " + b + ")";

        //normalize: x * a + y * b = 1
        xgcd[2].divideByLC(xgcd[0]);
        xgcd[1].divideByLC(xgcd[0]);
        xgcd[0].monic();

        return xgcd;
    }

    private static <PolyZp extends IUnivariatePolynomial<PolyZp>> PolyZp ensureMonic(PolyZp p) {
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
    public static List<lUnivariatePolynomialZp> liftFactorization(long modulus,
                                                                  long desiredBound,
                                                                  lUnivariatePolynomialZ poly,
                                                                  List<lUnivariatePolynomialZp> modularFactors,
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
    public static List<lUnivariatePolynomialZp> liftFactorization(long modulus,
                                                                  long finalModulus,
                                                                  int nIterations,
                                                                  lUnivariatePolynomialZ poly,
                                                                  List<lUnivariatePolynomialZp> modularFactors,
                                                                  boolean quadratic) {
        assert nIterations > 0;

        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.modulus(finalModulus, true).monic());

        lUnivariatePolynomialZp factory = modularFactors.get(0);
        lUnivariatePolynomialZp
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        LiftableQuintet<lUnivariatePolynomialZp> hensel = quadratic
                ? createQuadraticLift(modulus, poly, aFactor, bFactor)
                : createLinearLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        lUnivariatePolynomialZp aFactorRaised = hensel.aFactorMod();
        lUnivariatePolynomialZp bFactorRaised = hensel.bFactorMod();

        ArrayList<lUnivariatePolynomialZp> result = new ArrayList<>();
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, aFactorRaised.asPolyZSymmetric(), modularFactors.subList(0, nHalf), quadratic));
        result.addAll(liftFactorization(modulus, finalModulus, nIterations, bFactorRaised.asPolyZSymmetric(), modularFactors.subList(nHalf, modularFactors.size()), quadratic));
        return result;
    }

    interface LiftFactory<PolyZp extends IUnivariatePolynomial<PolyZp>> {
        LiftableQuintet<UnivariatePolynomial<BigInteger>> createLift(BigInteger modulus, UnivariatePolynomial<BigInteger> polyZ, PolyZp aFactor, PolyZp bFactor);
    }

    /** actual multifactor Hensel lifting implementation **/
    static <PolyZp extends IUnivariatePolynomial<PolyZp>>
    List<UnivariatePolynomial<BigInteger>> liftFactorization0(BigInteger modulus,
                                                              BigInteger finalModulus,
                                                              int nIterations,
                                                              UnivariatePolynomial<BigInteger> poly,
                                                              List<PolyZp> modularFactors,
                                                              LiftFactory<PolyZp> liftFactory) {
        assert nIterations > 0;

        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.setDomain(new IntegersModulo(finalModulus)).monic());

        PolyZp factory = modularFactors.get(0);
        PolyZp
                aFactor = factory.createOne(),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        LiftableQuintet<UnivariatePolynomial<BigInteger>> hensel = liftFactory.createLift(modulus, poly, aFactor, bFactor);
        hensel.lift(nIterations);

        UnivariatePolynomial<BigInteger> aFactorRaised = hensel.aFactorMod();
        UnivariatePolynomial<BigInteger> bFactorRaised = hensel.bFactorMod();

        ArrayList<UnivariatePolynomial<BigInteger>> result = new ArrayList<>();
        result.addAll(liftFactorization0(modulus, finalModulus, nIterations, UnivariatePolynomial.asPolyZSymmetric(aFactorRaised), modularFactors.subList(0, nHalf), liftFactory));
        result.addAll(liftFactorization0(modulus, finalModulus, nIterations, UnivariatePolynomial.asPolyZSymmetric(bFactorRaised), modularFactors.subList(nHalf, modularFactors.size()), liftFactory));
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
    public static List<UnivariatePolynomial<BigInteger>> liftFactorization(BigInteger modulus,
                                                                           BigInteger desiredBound,
                                                                           UnivariatePolynomial<BigInteger> poly,
                                                                           List<lUnivariatePolynomialZp> modularFactors,
                                                                           boolean quadratic) {
        if (!quadratic && !modulus.isLong())
            throw new IllegalArgumentException("Only max 64-bit modulus for linear lift allowed.");
        LiftingInfo im = nIterations(modulus, desiredBound, quadratic);
        if (im.nIterations == 0)
            return modularFactors.stream().map(lUnivariatePolynomialZp::toBigPoly).collect(Collectors.toList());
        LiftFactory<lUnivariatePolynomialZp> factory = quadratic ? HenselLifting::createQuadraticLift : HenselLifting::createLinearLift;
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
    public static List<UnivariatePolynomial<BigInteger>> liftFactorizationQuadratic(BigInteger modulus,
                                                                                    BigInteger desiredBound,
                                                                                    UnivariatePolynomial<BigInteger> poly,
                                                                                    List<UnivariatePolynomial<BigInteger>> modularFactors) {
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
    public static List<UnivariatePolynomial<BigInteger>> liftFactorization(BigInteger modulus,
                                                                           BigInteger desiredBound,
                                                                           UnivariatePolynomial<BigInteger> poly,
                                                                           List<lUnivariatePolynomialZp> modularFactors) {
        return liftFactorization(poly, modularFactors, new AdaptiveLift(modulus, desiredBound));
    }

    /** actual multifactor Hensel lifting implementation **/
    private static List<UnivariatePolynomial<BigInteger>> liftFactorization(UnivariatePolynomial<BigInteger> poly,
                                                                            List<lUnivariatePolynomialZp> modularFactors,
                                                                            AdaptiveLift lifter) {
        // for the future:
        // recursion may be replaced with precomputed binary tree
        // for now the major part of execution time (~99%) is spent in actual lifting step, so irrelevant

        if (modularFactors.size() == 1)
            return Collections.singletonList(poly.setDomain(new IntegersModulo(lifter.finalModulus)).monic());

        lUnivariatePolynomialZp factory = modularFactors.get(0);
        lUnivariatePolynomialZp
                aFactor = factory.createOne(),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));

        UnivariatePolynomial<BigInteger>[] lifted = lifter.lift(poly, aFactor, bFactor);
        UnivariatePolynomial<BigInteger> aFactorRaised = lifted[0];
        UnivariatePolynomial<BigInteger> bFactorRaised = lifted[1];

        ArrayList<UnivariatePolynomial<BigInteger>> result = new ArrayList<>();
        result.addAll(liftFactorization(UnivariatePolynomial.asPolyZSymmetric(aFactorRaised), modularFactors.subList(0, nHalf), lifter));
        result.addAll(liftFactorization(UnivariatePolynomial.asPolyZSymmetric(bFactorRaised), modularFactors.subList(nHalf, modularFactors.size()), lifter));
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

        @SuppressWarnings("unchecked")
        UnivariatePolynomial<BigInteger>[] lift(UnivariatePolynomial<BigInteger> poly, lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
            boolean quadratic = nLinearIterations == -1;
            LiftableQuintet<UnivariatePolynomial<BigInteger>> lift =
                    quadratic
                            ? createQuadraticLift(initialModulus, poly, a.toBigPoly(), b.toBigPoly())
                            : createLinearLift(initialModulus, poly, a, b);
            lift.lift(quadratic ? nQuadraticIterations : nLinearIterations);
            return new UnivariatePolynomial[]{lift.aFactorMod(), lift.bFactorMod()};
        }
    }

    private static <T extends IUnivariatePolynomial<T>> void assertHenselLift(LiftableQuintet<T> lift) {
        assert lift.polyMod().equals(lift.aFactorMod().clone().multiply(lift.bFactorMod())) : lift.toString();
        assert (lift.aCoFactorMod() == null && lift.bCoFactorMod() == null) ||
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod()))
                        .isOne() :
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod())) + "  --- " + ((UnivariatePolynomial<BigInteger>) lift.aFactorMod()).domain;
    }

    /* ************************************ Quadratic lifts ************************************ */

    /** data used in Hensel lifting **/
    static abstract class QuadraticLiftAbstract<PolyZp extends IUnivariatePolynomial<PolyZp>>
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
    public static final class lQuadraticLift extends QuadraticLiftAbstract<lUnivariatePolynomialZp> {
        /** The modulus */
        public long modulus;
        /** Initial Z[x] poly **/
        public final lUnivariatePolynomialZ base;

        public lQuadraticLift(long modulus, lUnivariatePolynomialZ base, lUnivariatePolynomialZp aFactor, lUnivariatePolynomialZp bFactor, lUnivariatePolynomialZp aCoFactor, lUnivariatePolynomialZp bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.modulus = modulus;
            this.base = base;
        }

        @Override
        public lUnivariatePolynomialZp polyMod() {
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
    public static final class bQuadraticLift extends QuadraticLiftAbstract<UnivariatePolynomial<BigInteger>> {
        /** The modulus */
        public IntegersModulo domain;
        /** Initial Z[x] poly **/
        public final UnivariatePolynomial<BigInteger> base;

        public bQuadraticLift(BigInteger modulus, UnivariatePolynomial<BigInteger> base, UnivariatePolynomial<BigInteger> aFactor, UnivariatePolynomial<BigInteger> bFactor, UnivariatePolynomial<BigInteger> aCoFactor, UnivariatePolynomial<BigInteger> bCoFactor) {
            super(aFactor, bFactor, aCoFactor, bCoFactor);
            this.domain = new IntegersModulo(modulus);
            this.base = base;
        }

        @Override
        public UnivariatePolynomial<BigInteger> polyMod() {
            return base.setDomain(domain);
        }

        @Override
        void prepare() {
            domain = new IntegersModulo(domain.modulus.multiply(domain.modulus));
            aFactor = aFactor.setDomainUnsafe(domain);
            bFactor = bFactor.setDomainUnsafe(domain);
            aCoFactor = aCoFactor.setDomainUnsafe(domain);
            bCoFactor = bCoFactor.setDomainUnsafe(domain);
        }
    }

    /* ************************************ Linear lifts ************************************ */

    private static class LinearLiftAbstract<PolyZ extends IUnivariatePolynomial<PolyZ>> {
        /** initial Z[x] poly */
        final PolyZ poly;
        /** lifted polynomials */
        PolyZ aFactor, bFactor, aCoFactor, bCoFactor;
        /** initial modular data */
        final lUnivariatePolynomialZp aFactorMod, aFactorModMonic, bFactorMod, aCoFactorMod, bCoFactorMod;
        /** precomputed inverses */
        final DivisionWithRemainder.InverseModMonomial<lUnivariatePolynomialZp> aFactorModMonicInv, bFactorModInv;

        public LinearLiftAbstract(PolyZ poly,
                                  PolyZ aFactor, PolyZ bFactor, PolyZ aCoFactor, PolyZ bCoFactor,
                                  lUnivariatePolynomialZp aFactorMod, lUnivariatePolynomialZp aFactorModMonic, lUnivariatePolynomialZp bFactorMod,
                                  lUnivariatePolynomialZp aCoFactorMod, lUnivariatePolynomialZp bCoFactorMod) {
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

        protected lUnivariatePolynomialZp aAdd, bAdd;

        final void calculateFactorsDiff(lUnivariatePolynomialZp diff) {
            aAdd = diff.clone();
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);
            aAdd = aAdd.multiply(bCoFactorMod);
            aAdd = PolynomialArithmetics.polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);

            bAdd = diff.clone();
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
            bAdd = bAdd.multiply(aCoFactorMod);
            bAdd = PolynomialArithmetics.polyMod(bAdd, bFactorMod, bFactorModInv, false);
        }

        final void calculateCoFactorsDiff(lUnivariatePolynomialZp diff) {
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
            extends LinearLiftAbstract<lUnivariatePolynomialZ>
            implements LiftableQuintet<lUnivariatePolynomialZp> {
        /** The initial modulus */
        public final long initialModulus;
        /** The modulus */
        public long modulus;

        private lLinearLift(long modulus, lUnivariatePolynomialZ poly,
                            lUnivariatePolynomialZp aFactor, lUnivariatePolynomialZp aFactorMonic, lUnivariatePolynomialZp bFactor,
                            lUnivariatePolynomialZp aCoFactor, lUnivariatePolynomialZp bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).asPolyZ(true).multiply(poly.lc()), bFactor.asPolyZ(true),
                    aCoFactor.asPolyZ(true), bCoFactor.asPolyZ(true),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialModulus = modulus;
            this.modulus = modulus;
        }

        private lLinearLift(long modulus, lUnivariatePolynomialZ poly,
                            lUnivariatePolynomialZp aFactor, lUnivariatePolynomialZp bFactor,
                            lUnivariatePolynomialZp aCoFactor, lUnivariatePolynomialZp bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public lUnivariatePolynomialZp polyMod() {return poly.modulus(modulus);}

        @Override
        public lUnivariatePolynomialZp aFactorMod() {return aFactor.modulus(modulus);}

        @Override
        public lUnivariatePolynomialZp bFactorMod() {return bFactor.modulus(modulus);}

        @Override
        public lUnivariatePolynomialZp aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.modulus(modulus);}

        @Override
        public lUnivariatePolynomialZp bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.modulus(modulus);}

        private void liftFactors() {
            lUnivariatePolynomialZp factorsDiff = poly.clone().subtract(aFactor.clone().multiply(bFactor))
                    .divideOrNull(modulus)
                    .modulus(initialModulus);

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.asPolyZ(false).multiply(modulus));
            bFactor = bFactor.add(bAdd.asPolyZ(false).multiply(modulus));
        }

        private void liftCoFactors() {
            lUnivariatePolynomialZp coFactorsDiff = aCoFactor.clone().multiply(aFactor)
                    .add(bCoFactor.clone().multiply(bFactor))
                    .decrement()
                    .negate()
                    .divideOrNull(modulus)
                    .modulus(initialModulus);

            calculateCoFactorsDiff(coFactorsDiff);
            aCoFactor = aCoFactor.add(aAdd.asPolyZ(false).multiply(modulus));
            bCoFactor = bCoFactor.add(bAdd.asPolyZ(false).multiply(modulus));
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
            extends LinearLiftAbstract<UnivariatePolynomial<BigInteger>>
            implements LiftableQuintet<UnivariatePolynomial<BigInteger>> {
        /** The initial modulus (less than 64-bit) */
        public final IntegersModulo initialDomain;
        /** The modulus */
        public IntegersModulo domain;

        private bLinearLift(BigInteger modulus, UnivariatePolynomial<BigInteger> poly,
                            lUnivariatePolynomialZp aFactor, lUnivariatePolynomialZp aFactorMonic, lUnivariatePolynomialZp bFactor,
                            lUnivariatePolynomialZp aCoFactor, lUnivariatePolynomialZp bCoFactor) {
            super(poly,
                    ensureMonic(aFactor).asPolyZ(true).toBigPoly().multiply(poly.lc()), bFactor.asPolyZ(false).toBigPoly(),
                    aCoFactor.asPolyZ(false).toBigPoly(), bCoFactor.asPolyZ(false).toBigPoly(),
                    aFactor, aFactorMonic, bFactor, aCoFactor, bCoFactor);
            this.initialDomain = new IntegersModulo(modulus);
            this.domain = new IntegersModulo(modulus);
            assert modulus.isLong();
        }

        private bLinearLift(BigInteger modulus, UnivariatePolynomial<BigInteger> poly,
                            lUnivariatePolynomialZp aFactor, lUnivariatePolynomialZp bFactor,
                            lUnivariatePolynomialZp aCoFactor, lUnivariatePolynomialZp bCoFactor) {
            this(modulus, poly, aFactor, aFactor.clone().monic(), bFactor, aCoFactor, bCoFactor);
        }

        @Override
        public UnivariatePolynomial<BigInteger> polyMod() {return poly.setDomain(domain);}

        @Override
        public UnivariatePolynomial<BigInteger> aFactorMod() {return aFactor.setDomain(domain);}

        @Override
        public UnivariatePolynomial<BigInteger> bFactorMod() {return bFactor.setDomain(domain);}

        @Override
        public UnivariatePolynomial<BigInteger> aCoFactorMod() {return aCoFactor == null ? null : aCoFactor.setDomain(domain);}

        @Override
        public UnivariatePolynomial<BigInteger> bCoFactorMod() {return bCoFactor == null ? null : bCoFactor.setDomain(domain);}

        private void liftFactors() {
            lUnivariatePolynomialZp factorsDiff = UnivariatePolynomial.asLongPolyZp(
                    poly.clone().subtract(aFactor.clone().multiply(bFactor))
                            .divideOrNull(domain.modulus)
                            .setDomain(initialDomain));

            calculateFactorsDiff(factorsDiff);
            aFactor = aFactor.add(aAdd.asPolyZ(false).toBigPoly().multiply(domain.modulus));
            bFactor = bFactor.add(bAdd.asPolyZ(false).toBigPoly().multiply(domain.modulus));
        }

        private void liftCoFactors() {
            lUnivariatePolynomialZp coFactorsDiff = UnivariatePolynomial.asLongPolyZp(
                    aCoFactor.clone().multiply(aFactor).add(bCoFactor.clone().multiply(bFactor)).decrement().negate()
                            .divideOrNull(domain.modulus)
                            .setDomain(initialDomain));

            calculateCoFactorsDiff(coFactorsDiff);

            aCoFactor = aCoFactor.add(aAdd.asPolyZ(false).toBigPoly().multiply(domain.modulus));
            bCoFactor = bCoFactor.add(bAdd.asPolyZ(false).toBigPoly().multiply(domain.modulus));
        }

        @Override
        public void lift() {
            liftFactors();
            liftCoFactors();
            domain = new IntegersModulo(domain.modulus.multiply(initialDomain.modulus));
        }

        @Override
        public void liftLast() {
            liftFactors();
            domain = new IntegersModulo(domain.modulus.multiply(initialDomain.modulus));
            aCoFactor = bCoFactor = null;
        }
    }
}
