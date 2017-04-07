package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;

import java.util.Arrays;

import static cc.r2.core.poly.univar.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly.univar.PolynomialGCD.PolynomialGCD;


/**
 * Square-free factorization of univariate polynomials over Z and Zp with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class SquareFreeFactorization {
    private SquareFreeFactorization() {}

    /**
     * Performs square-free factorization of a {@code poly} in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static lFactorDecomposition<lMutablePolynomialZ> SquareFreeFactorization(lMutablePolynomialZ poly) {
        if (poly.isConstant())
            return lFactorDecomposition.oneFactor(poly.lc());

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun(poly);

        lMutablePolynomialZ expFree = lMutablePolynomialZ.create(Arrays.copyOfRange(poly.data, exponent, poly.degree + 1));
        lFactorDecomposition<lMutablePolynomialZ> fd = SquareFreeFactorizationYun(expFree);
        fd.addFactor(lMutablePolynomialZ.monomial(1, exponent), 1);
        return fd;
    }

    /**
     * Performs square-free factorization of a {@code poly} in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static bFactorDecomposition<bMutablePolynomialZ> SquareFreeFactorization(bMutablePolynomialZ poly) {
        if (poly.isConstant())
            return bFactorDecomposition.oneFactor(poly.lc());

        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent].isZero()) { ++exponent; }
        if (exponent == 0)
            return SquareFreeFactorizationYun(poly);

        bMutablePolynomialZ expFree = poly.getRange(exponent, poly.degree + 1);
        bFactorDecomposition<bMutablePolynomialZ> fd = SquareFreeFactorizationYun(expFree);
        fd.addFactor(bMutablePolynomialZ.monomial(BigInteger.ONE, exponent), 1);
        return fd;
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static lFactorDecomposition<lMutablePolynomialZ> SquareFreeFactorizationYun(lMutablePolynomialZ poly) {
        if (poly.isConstant())
            return lFactorDecomposition.oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return lFactorDecomposition.oneFactor(poly, content);

        lFactorDecomposition<lMutablePolynomialZ> factorization = new lFactorDecomposition<>();
        SquareFreeFactorizationYun(poly, factorization);
        return factorization.setNumericFactor(content);
    }

    /**
     * Performs square-free factorization of a poly using Yun's algorithm in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static bFactorDecomposition<bMutablePolynomialZ> SquareFreeFactorizationYun(bMutablePolynomialZ poly) {
        if (poly.isConstant())
            return bFactorDecomposition.oneFactor(poly.lc());

        BigInteger content = poly.content();
        if (poly.lc().signum() < 0)
            content = content.negate();

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return bFactorDecomposition.oneFactor(poly, content);

        bFactorDecomposition<bMutablePolynomialZ> factorization = new bFactorDecomposition<>();
        SquareFreeFactorizationYun(poly, factorization);
        return factorization.setNumericFactor(content);
    }

    private static <T extends IMutablePolynomialZ<T>> void SquareFreeFactorizationYun(T poly, FactorDecomposition<T> factorization) {
        T derivative = poly.derivative(), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        T quot = divideAndRemainder(poly, gcd, false)[0], // safely destroy (cloned) poly (not used further)
                dQuot = divideAndRemainder(derivative, gcd, false)[0]; // safely destroy (cloned) derivative (not used further)

        int i = 0;
        while (!quot.isConstant()) {
            ++i;
            dQuot = dQuot.subtract(quot.derivative());
            T factor = PolynomialGCD(quot, dQuot);
            quot = divideAndRemainder(quot, factor, false)[0]; // can destroy quot in divideAndRemainder
            dQuot = divideAndRemainder(dQuot, factor, false)[0]; // can destroy dQuot in divideAndRemainder
            if (!factor.isOne())
                factorization.addFactor(factor, i);
        }
    }

    /**
     * Performs square-free factorization of a poly using Musser's algorithm in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static lFactorDecomposition<lMutablePolynomialZ> SquareFreeFactorizationMusser(lMutablePolynomialZ poly) {
        if (poly.isConstant())
            return lFactorDecomposition.oneFactor(poly.lc());

        long content = poly.content();
        if (poly.lc() < 0)
            content = -content;

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return lFactorDecomposition.oneFactor(poly, content);

        lFactorDecomposition<lMutablePolynomialZ> factorization = new lFactorDecomposition<>();
        SquareFreeFactorizationMusser(poly, factorization);
        return factorization.setNumericFactor(content);
    }

    /**
     * Performs square-free factorization of a poly using Musser's algorithm in Z[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("ConstantConditions")
    public static bFactorDecomposition<bMutablePolynomialZ> SquareFreeFactorizationMusser(bMutablePolynomialZ poly) {
        if (poly.isConstant())
            return bFactorDecomposition.oneFactor(poly.lc());

        BigInteger content = poly.content();
        if (poly.lc().signum() < 0)
            content = content.negate();

        poly = poly.clone().divideOrNull(content);
        if (poly.degree <= 1)
            return bFactorDecomposition.oneFactor(poly, content);

        bFactorDecomposition<bMutablePolynomialZ> factorization = new bFactorDecomposition<>();
        SquareFreeFactorizationMusser(poly, factorization);
        return factorization.setNumericFactor(content);
    }

    private static <T extends IMutablePolynomialZ<T>> void SquareFreeFactorizationMusser(T poly, FactorDecomposition<T> factorization) {
        T derivative = poly.derivative(), gcd = PolynomialGCD(poly, derivative);
        if (gcd.isConstant()) {
            factorization.addFactor(poly, 1);
            return;
        }

        T quot = divideAndRemainder(poly, gcd, false)[0]; // safely destroy (cloned) poly
        int i = 0;
        while (true) {
            ++i;
            T nextQuot = PolynomialGCD(gcd, quot);
            gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // safely destroy gcd (reassigned)
            T factor = divideAndRemainder(quot, nextQuot, false)[0]; // safely destroy quot (reassigned further)
            if (!factor.isConstant())
                factorization.addFactor(factor, i);
            if (nextQuot.isConstant())
                break;
            quot = nextQuot;
        }
    }

    /**
     * Performs square-free factorization of a {@code poly} in Zp[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static lFactorDecomposition<lMutablePolynomialZp> SquareFreeFactorization(lMutablePolynomialZp poly) {
        poly = poly.clone();
        long lc = poly.lc();
        //make poly monic
        poly = poly.monic();

        if (poly.isConstant())
            return lFactorDecomposition.oneFactor(lc);

        if (poly.degree <= 1)
            return lFactorDecomposition.oneFactor(poly, lc);

        lFactorDecomposition<lMutablePolynomialZp> factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent] == 0) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly);
        else {
            lMutablePolynomialZp expFree = poly.getRange(exponent, poly.degree + 1);
            factorization = SquareFreeFactorizationMusser0(expFree);
            factorization.addFactor(lMutablePolynomialZp.createMonomial(poly.modulus, 1, exponent), 1);
        }

        return factorization.setNumericFactor(lc);
    }

    /**
     * Performs square-free factorization of a {@code poly} in Zp[x].
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    public static bFactorDecomposition<bMutablePolynomialZp> SquareFreeFactorization(bMutablePolynomialZp poly) {
        if (poly.isLong())
            return lFactorDecomposition.convert(SquareFreeFactorizationMusser0(poly.toLong()));

//        if (modulus >= Integer.MAX_VALUE)// <- just to be on the safe side)
//            throw new IllegalArgumentException();

        poly = poly.clone();
        BigInteger lc = poly.lc();
        //make poly monic
        poly = poly.monic();

        if (poly.isConstant())
            return bFactorDecomposition.oneFactor(lc);

        if (poly.degree <= 1)
            return bFactorDecomposition.oneFactor(poly, lc);

        bFactorDecomposition<bMutablePolynomialZp> factorization;
        // x^2 + x^3 -> x^2 (1 + x)
        int exponent = 0;
        while (exponent <= poly.degree && poly.data[exponent].isZero()) { ++exponent; }
        if (exponent == 0)
            factorization = SquareFreeFactorizationMusser0(poly);
        else {
            bMutablePolynomialZp expFree = poly.getRange(exponent, poly.degree + 1);
            factorization = SquareFreeFactorizationMusser0(expFree);
            factorization.addFactor(bMutablePolynomialZp.createMonomial(poly.modulus, BigInteger.ONE, exponent), 1);
        }

        return factorization.setNumericFactor(lc);
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static lFactorDecomposition<lMutablePolynomialZp> SquareFreeFactorizationMusser0(lMutablePolynomialZp poly) {
        poly.monic();
        if (poly.isConstant())
            return lFactorDecomposition.oneFactor(poly.lc());

        if (poly.degree <= 1)
            return lFactorDecomposition.oneFactor(poly, 1);

        lMutablePolynomialZp derivative = poly.derivative();
        if (!derivative.isZero()) {
            lMutablePolynomialZp gcd = PolynomialGCD(poly, derivative);
            if (gcd.isConstant())
                return lFactorDecomposition.oneFactor(poly, 1);
            lMutablePolynomialZp quot = divideAndRemainder(poly, gcd, false)[0]; // can safely destroy poly (not used further)

            lFactorDecomposition<lMutablePolynomialZp> result = new lFactorDecomposition<>();
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                lMutablePolynomialZp nextQuot = PolynomialGCD(gcd, quot);
                lMutablePolynomialZp factor = divideAndRemainder(quot, nextQuot, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant())
                    result.addFactor(factor.monic(), i);

                gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant()) {
                gcd = pRoot(gcd);
                lFactorDecomposition<lMutablePolynomialZp> gcdFactorization = SquareFreeFactorizationMusser0(gcd);
                gcdFactorization.raiseExponents(poly.modulus);
                result.addAll(gcdFactorization);
                return result;
            } else
                return result;
        } else {
            lMutablePolynomialZp pRoot = pRoot(poly);
            lFactorDecomposition<lMutablePolynomialZp> fd = SquareFreeFactorizationMusser0(pRoot);
            fd.raiseExponents(poly.modulus);
            return fd.setNumericFactor(1);
        }
    }

    /** {@code poly} will be destroyed */
    @SuppressWarnings("ConstantConditions")
    private static bFactorDecomposition<bMutablePolynomialZp> SquareFreeFactorizationMusser0(bMutablePolynomialZp poly) {
        poly.monic();
        if (poly.isConstant())
            return bFactorDecomposition.oneFactor(poly.lc());

        if (poly.degree <= 1)
            return bFactorDecomposition.oneFactor(poly, BigInteger.ONE);

        bMutablePolynomialZp derivative = poly.derivative();
        if (!derivative.isZero()) {
            bMutablePolynomialZp gcd = PolynomialGCD(poly, derivative);
            if (gcd.isConstant())
                return bFactorDecomposition.oneFactor(poly, BigInteger.ONE);
            bMutablePolynomialZp quot = divideAndRemainder(poly, gcd, false)[0]; // can safely destroy poly (not used further)

            bFactorDecomposition<bMutablePolynomialZp> result = new bFactorDecomposition<>();
            int i = 0;
            //if (!quot.isConstant())
            while (true) {
                ++i;
                bMutablePolynomialZp nextQuot = PolynomialGCD(gcd, quot);
                bMutablePolynomialZp factor = divideAndRemainder(quot, nextQuot, false)[0]; // can safely destroy quot (reassigned further)
                if (!factor.isConstant())
                    result.addFactor(factor.monic(), i);

                gcd = divideAndRemainder(gcd, nextQuot, false)[0]; // can safely destroy gcd
                if (nextQuot.isConstant())
                    break;
                quot = nextQuot;
            }
            if (!gcd.isConstant())
                throw new IllegalStateException();
            return result;
        } else
            throw new IllegalStateException();
    }

    /** p-th root of poly */
    private static lMutablePolynomialZp pRoot(lMutablePolynomialZp poly) {
        long modulus = poly.modulus;
        assert poly.degree % modulus == 0;
        long[] rootData = new long[poly.degree / LongArithmetics.safeToInt(modulus) + 1];
        Arrays.fill(rootData, 0);
        for (int i = poly.degree; i >= 0; --i)
            if (poly.data[i] != 0) {
                assert i % modulus == 0;
                rootData[i / (int) modulus] = poly.data[i];
            }
        return poly.createFromArray(rootData);
    }

    /**
     * Returns {@code true} if {@code poly} is square-free and {@code false} otherwise
     *
     * @param poly the polynomial
     * @return {@code true} if {@code poly} is square-free and {@code false} otherwise
     */
    public static <T extends IMutablePolynomial<T>> boolean isSquareFree(T poly) {
        return PolynomialGCD(poly, poly.derivative()).isConstant();
    }

    /**
     * Performs square-free factorization of a {@code poly}.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("unchecked")
    public static <T extends lMutablePolynomialAbstract<T>> lFactorDecomposition<T> SquareFreeFactorization(T poly) {
        if (poly instanceof lMutablePolynomialZ)
            return (lFactorDecomposition<T>) SquareFreeFactorization((lMutablePolynomialZ) poly);
        else
            return (lFactorDecomposition<T>) SquareFreeFactorization((lMutablePolynomialZp) poly);
    }

    /**
     * Performs square-free factorization of a {@code poly}.
     *
     * @param poly the polynomial
     * @return square-free decomposition
     */
    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> FactorDecomposition<T> SquareFreeFactorization(T poly) {
        if (poly instanceof lMutablePolynomialZ)
            return (FactorDecomposition<T>) SquareFreeFactorization((lMutablePolynomialZ) poly);
        else if (poly instanceof lMutablePolynomialZp)
            return (FactorDecomposition<T>) SquareFreeFactorization((lMutablePolynomialZp) poly);
        else if (poly instanceof bMutablePolynomialZ)
            return (FactorDecomposition<T>) SquareFreeFactorization((bMutablePolynomialZ) poly);
        else
            return (FactorDecomposition<T>) SquareFreeFactorization((bMutablePolynomialZp) poly);
    }

    /**
     * Returns square-free part of the {@code poly}
     *
     * @param poly the polynomial
     * @return square free part
     */
    public static <T extends IMutablePolynomial<T>> T SquareFreePart(T poly) {
        return SquareFreeFactorization(poly).factors.stream().filter(x -> !x.isMonomial()).reduce(poly.createOne(), IMutablePolynomial<T>::multiply);
    }
}