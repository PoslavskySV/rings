package cc.redberry.rings.poly.univar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.PolynomialFactorDecomposition;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.multivar.MultivariateSquareFreeFactorization;

import java.util.Arrays;

import static cc.redberry.rings.poly.univar.Conversions64bit.asOverZp64;
import static cc.redberry.rings.poly.univar.Conversions64bit.canConvertToZp64;


/**
 * Square-free factorization of univariate polynomials over Z and Zp.
 *
 * @since 1.0
 */
public final class UnivariateSquareFreeFactorization {
    private UnivariateSquareFreeFactorization() {}

    /**
     * Returns {@code true} if {@code poly} is square-free and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if {@code poly} is square-free and {@code false} otherwise
     */
    public static <T extends IUnivariatePolynomial<T>> boolean isSquareFree(T poly) {
        return UnivariateGCD.PolynomialGCD(poly, poly.derivative()).isConstant();
    }

    /**
     * Performs square-free factorization of a {@code poly}.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("unchecked")
    public static <T extends IUnivariatePolynomial<T>> PolynomialFactorDecomposition<T> SquareFreeFactorization(T poly) {
        if (UnivariateFactorization.isOverMultivariate(poly))
            return (PolynomialFactorDecomposition<T>) UnivariateFactorization.FactorOverMultivariate((UnivariatePolynomial) poly, MultivariateSquareFreeFactorization::SquareFreeFactorization);
        else if (UnivariateFactorization.isOverUnivariate(poly))
            return (PolynomialFactorDecomposition<T>) UnivariateFactorization.FactorOverUnivariate((UnivariatePolynomial) poly, MultivariateSquareFreeFactorization::SquareFreeFactorization);
        else if (poly.coefficientRingCharacteristic().isZero())
            return SquareFreeFactorizationYunZeroCharacteristics(poly);
        else
            return SquareFreeFactorizationMusser(poly);
    }

    /**
     * Returns square-free part of the {@code poly}
     *
     * @param poly the polynomial
     * @return square free part
     */
    public static <T extends IUnivariatePolynomial<T>> T SquareFreePart(T poly) {
        return SquareFreeFactorization(poly).factors.stream().filter(x -> !x.isMonomial()).reduce(poly.createOne(), IUnivariatePolynomial<T>::multiply);
    }

    /**
     * Performs square-free factorization of a {@code poly} which coefficient ring has zero characteristic using Yun's
     * algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> PolynomialFactorDecomposition<Poly>
    SquareFreeFactorizationYunZeroCharacteristics(Poly poly) {
        if (!poly.coefficientRingCharacteristic().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        if (poly.isConstant())
            return PolynomialFactorDecomposition.of(poly);

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree() && poly.isZeroAt(exponent)) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun0(poly);

        Poly expFree = poly.getRange(exponent, poly.degree() + 1);
        PolynomialFactorDecomposition<Poly> fd = SquareFreeFactorizationYun0(expFree);
        fd.addFactor(poly.createMonomial(1), exponent);
        return fd;
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    static <Poly extends IUnivariatePolynomial<Poly>> PolynomialFactorDecomposition<Poly> SquareFreeFactorizationYun0(Poly poly) {
        if (poly.isConstant())
            return PolynomialFactorDecomposition.of(poly);

        Poly content = poly.contentAsPoly();
        if (poly.signumOfLC() < 0)
            content = content.negate();

        poly = poly.clone().divideByLC(content);
        if (poly.degree() <= 1)
            return PolynomialFactorDecomposition.of(content, poly);

        PolynomialFactorDecomposition<Poly> factorization = PolynomialFactorDecomposition.of(content);
        SquareFreeFactorizationYun0(poly, factorization);
        return factorization;
    }

    private static <Poly extends IUnivariatePolynomial<Poly>> void SquareFreeFactorizationYun0(Poly poly, PolynomialFactorDecomposition<Poly> factorization) {
        Poly derivative = poly.derivative(), gcd = UnivariateGCD.PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        Poly quot = UnivariateDivision.divideAndRemainder(poly, gcd, false)[0], // safely destroy (cloned) poly (not used further)
                dQuot = UnivariateDivision.divideAndRemainder(derivative, gcd, false)[0]; // safely destroy (cloned) derivative (not used further)

        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            dQuot = dQuot.subtract(quot.derivative());
            Poly factor = UnivariateGCD.PolynomialGCD(quot, dQuot);
            quot = UnivariateDivision.divideAndRemainder(quot, factor, false)[0]; // can destroy quot in divideAndRemainder
            dQuot = UnivariateDivision.divideAndRemainder(dQuot, factor, false)[0]; // can destroy dQuot in divideAndRemainder
            if (!factor.isOne())
                factorization.addFactor(factor, i);
        }
    }

    /**
     * Performs square-free factorization of a poly which coefficient ring has zero characteristic using Musser's
     * algorithm.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static <Poly extends IUnivariatePolynomial<Poly>> PolynomialFactorDecomposition<Poly>
    SquareFreeFactorizationMusserZeroCharacteristics(Poly poly) {
        if (!poly.coefficientRingCharacteristic().isZero())
            throw new IllegalArgumentException("Characteristics 0 expected");

        if (poly.isConstant())
            return PolynomialFactorDecomposition.of(poly);

        Poly content = poly.contentAsPoly();
        if (poly.signumOfLC() < 0)
            content = content.negate();

        poly = poly.clone().divideByLC(content);
        if (poly.degree() <= 1)
            return PolynomialFactorDecomposition.of(content, poly);

        PolynomialFactorDecomposition<Poly> factorization = PolynomialFactorDecomposition.of(content);
        SquareFreeFactorizationMusserZeroCharacteristics0(poly, factorization);
        return factorization;
    }

    private static <Poly extends IUnivariatePolynomial<Poly>> void
    SquareFreeFactorizationMusserZeroCharacteristics0(Poly poly, PolynomialFactorDecomposition<Poly> factorization) {
        Poly derivative = poly.derivative(), gcd = UnivariateGCD.PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        Poly quot = UnivariateDivision.divideAndRemainder(poly, gcd, false)[0]; // safely destroy (cloned) poly
        int i = 0;
        while (true) {
            ++i;
            Poly nextQuot = UnivariateGCD.PolynomialGCD(gcd, quot);
            gcd = UnivariateDivision.divideAndRemainder(gcd, nextQuot, false)[0]; // safely destroy gcd (reassigned)
            Poly factor = UnivariateDivision.divideAndRemainder(quot, nextQuot, false)[0]; // safely destroy quot (reassigned further)
            if (!factor.isConstant())
                factorization.addFactor(factor, i);
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
    }

    /**
     * Performs square-free factorization of a {@code poly} using Musser's algorithm (both zero and non-zero
     * characteristic of coefficient ring allowed).
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> PolynomialFactorDecomposition<Poly> SquareFreeFactorizationMusser(Poly poly) {
        if (canConvertToZp64(poly))
            return SquareFreeFactorizationMusser(asOverZp64(poly)).mapTo(Conversions64bit::convert);

        poly = poly.clone();
        Poly lc = poly.lcAsPoly();
        //make poly monic
        poly = poly.monic();

        if (poly.isConstant())
            return PolynomialFactorDecomposition.of(lc);

        if (poly.degree() <= 1)
            return PolynomialFactorDecomposition.of(lc, poly);

        PolynomialFactorDecomposition<Poly> factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree() && poly.isZeroAt(exponent)) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly);
        else {
            Poly expFree = poly.getRange(exponent, poly.degree() + 1);
            factorization = SquareFreeFactorizationMusser0(expFree);
            factorization.addFactor(poly.createMonomial(1), exponent);
        }

        return factorization.setUnit(lc);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static <Poly extends IUnivariatePolynomial<Poly>> PolynomialFactorDecomposition<Poly> SquareFreeFactorizationMusser0(Poly poly) {
        poly.monic();
        if (poly.isConstant())
            return PolynomialFactorDecomposition.of(poly);

        if (poly.degree() <= 1)
            return PolynomialFactorDecomposition.of(poly);

        Poly derivative = poly.derivative();
        if (!derivative.isZero()) {
            Poly gcd = UnivariateGCD.PolynomialGCD(poly, derivative);
            if (gcd.isConstant())
                return PolynomialFactorDecomposition.of(poly);
            Poly quot = UnivariateDivision.divideAndRemainder(poly, gcd, false)[0]; // can safely destroy poly (not used further)

            PolynomialFactorDecomposition<Poly> result = PolynomialFactorDecomposition.of(poly.createOne());
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                Poly nextQuot = UnivariateGCD.PolynomialGCD(gcd, quot);
                Poly factor = UnivariateDivision.divideAndRemainder(quot, nextQuot, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant())
                    result.addFactor(factor.monic(), i);

                gcd = UnivariateDivision.divideAndRemainder(gcd, nextQuot, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd);
                PolynomialFactorDecomposition<Poly> gcdFactorization = SquareFreeFactorizationMusser0(gcd);
                gcdFactorization.raiseExponents(poly.coefficientRingCharacteristic().intValueExact());
                result.addAll(gcdFactorization);
                return result;
            } else
                return result;
        } else {
            Poly pRoot = pRoot(poly);
            PolynomialFactorDecomposition<Poly> fd = SquareFreeFactorizationMusser0(pRoot);
            fd.raiseExponents(poly.coefficientRingCharacteristic().intValueExact());
            return fd.setUnit(poly.createOne());
        }
    }

    /** p-th root of poly */
    @SuppressWarnings("unchecked")
    private static <Poly extends IUnivariatePolynomial<Poly>> Poly pRoot(Poly poly) {
        if (poly instanceof UnivariatePolynomialZp64)
            return (Poly) pRoot((UnivariatePolynomialZp64) poly);
        else if (poly instanceof UnivariatePolynomial)
            return (Poly) pRoot((UnivariatePolynomial) poly);
        else
            throw new RuntimeException(poly.getClass().toString());
    }

    /** p-th root of poly */
    private static UnivariatePolynomialZp64 pRoot(UnivariatePolynomialZp64 poly) {
        if (poly.ring.modulus > Integer.MAX_VALUE)
            throw new IllegalArgumentException("Too big modulus: " + poly.ring.modulus);
        int modulus = MachineArithmetic.safeToInt(poly.ring.modulus);
        assert poly.degree % modulus == 0;
        assert !poly.ring.isPerfectPower(); // just in case
        long[] rootData = new long[poly.degree / modulus + 1];
        Arrays.fill(rootData, 0);
        for (int i = poly.degree; i >= 0; --i)
            if (poly.data[i] != 0) {
                assert i % modulus == 0;
                rootData[i / modulus] = poly.data[i];
            }
        return poly.createFromArray(rootData);
    }

    /** p-th root of poly */
    private static <E> UnivariatePolynomial<E> pRoot(UnivariatePolynomial<E> poly) {
        if (!poly.coefficientRingCharacteristic().isInt())
            throw new IllegalArgumentException("Infinite or too large ring: " + poly.ring);
        Ring<E> ring = poly.ring;
        // p^(m -1) used for computing p-th root of elements
        BigInteger inverseFactor = ring.cardinality().divide(ring.characteristic());
        int modulus = poly.coefficientRingCharacteristic().intValueExact();
        assert poly.degree % modulus == 0;
        E[] rootData = poly.ring.createZeroesArray(poly.degree / modulus + 1);
        for (int i = poly.degree; i >= 0; --i)
            if (!poly.ring.isZero(poly.data[i])) {
                assert i % modulus == 0;
                rootData[i / modulus] = ring.pow(poly.data[i], inverseFactor); // pRoot(poly.data[i], ring);
            }
        return poly.createFromArray(rootData);
    }
}