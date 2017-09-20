package cc.redberry.rings.bigint;

import static cc.redberry.rings.bigint.BigInteger.ONE;

/**
 * @since 1.0
 */
public final class BigIntegerUtil {
    private BigIntegerUtil() {}

    public static BigInteger max(BigInteger a, BigInteger b) {
        return a.compareTo(b) > 0 ? a : b;
    }

    public static BigInteger abs(BigInteger a) {
        return a.abs();
    }

    public static BigInteger gcd(BigInteger a, BigInteger b) {
        return a.gcd(b);
    }

    /**
     * Returns the greatest common an array of longs
     *
     * @param integers array of longs
     * @param from     from position (inclusive)
     * @param to       to position (exclusive)
     * @return greatest common divisor of array
     */
    public static BigInteger gcd(final BigInteger[] integers, int from, int to) {
        if (integers.length < 2)
            throw new IllegalArgumentException();
        BigInteger gcd = gcd(integers[from], integers[from + 1]);
        if (gcd.isOne())
            return gcd;
        for (int i = from + 2; i < to; i++) {
            gcd = gcd(integers[i], gcd);
            if (gcd.isOne())
                return gcd;
        }
        return gcd;
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e}
     * @throws ArithmeticException if the result overflows a long
     */
    public static BigInteger pow(final long base, long exponent) {
        return pow(BigInteger.valueOf(base), exponent);
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e}
     * @throws ArithmeticException if the result overflows a long
     */
    public static BigInteger pow(final BigInteger base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        BigInteger result = ONE;
        BigInteger k2p = base;
        for (; ; ) {
            if ((exponent & 1) != 0)
                result = result.multiply(k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p);
        }
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e}
     * @throws ArithmeticException if the result overflows a long
     */
    public static BigInteger pow(final BigInteger base, BigInteger exponent) {
        if (exponent.signum() < 0)
            throw new IllegalArgumentException();

        BigInteger result = ONE;
        BigInteger k2p = base;
        for (; ; ) {
            if (exponent.testBit(0))
                result = result.multiply(k2p);
            exponent = exponent.shiftRight(1);
            if (exponent.isZero())
                return result;
            k2p = k2p.multiply(k2p);
        }
    }

    /**
     * Returns floor square root of {@code val}
     *
     * @param val the number
     * @return floor square root
     */
    public static BigInteger sqrtFloor(BigInteger val)
            throws IllegalArgumentException {
        if (val.signum() < 0)
            throw new IllegalArgumentException("Negative argument.");
        if (val.isZero() || val.isOne())
            return val;

        BigInteger y;
        // starting with y = x / 2 avoids magnitude issues with x squared
        for (y = val.shiftRight(1);
             y.compareTo(val.divide(y)) > 0;
             y = ((val.divide(y)).add(y)).shiftRight(1)) {}

        return y;
    }

    /**
     * Returns ceil square root of {@code val}
     *
     * @param val the number
     * @return floor square root
     */
    public static BigInteger sqrtCeil(BigInteger val)
            throws IllegalArgumentException {
        if (val.signum() < 0)
            throw new IllegalArgumentException("Negative argument.");
        if (val.isZero() || val.isOne())
            return val;

        BigInteger y;
        // starting with y = x / 2 avoids magnitude issues with x squared
        for (y = val.shiftRight(1);
             y.compareTo(val.divide(y)) > 0;
             y = ((val.divide(y)).add(y)).shiftRight(1)) {}

        if (val.compareTo(y.multiply(y)) == 0)
            return y;
        else
            return y.add(BigInteger.ONE);

    }

    /**
     * Tests whether {@code n} is a perfect power {@code n == a^b} and returns {@code {a, b}} if so and {@code null}
     * otherwise
     *
     * @param n the number
     * @return array {@code {a, b}} so that {@code n = a^b} or {@code null} is {@code n} is not a perfect power
     */
    public static BigInteger[] perfectPowerDecomposition(BigInteger n) {
        if (n.signum() < 0) {
            n = n.negate();
            BigInteger[] ipp = perfectPowerDecomposition(n);
            if (ipp == null)
                return null;
            if (ipp[1].testBit(0))
                return null;
            ipp[0] = ipp[0].negate();
            return ipp;
        }

        if (n.bitCount() == 1)
            return new BigInteger[]{BigInteger.TWO, BigInteger.valueOf(n.bitLength() - 1)};

        int lgn = 1 + n.bitLength();
        for (int b = 2; b < lgn; b++) {
            //b lg a = lg n
            BigInteger lowa = BigInteger.ONE;
            BigInteger higha = BigInteger.ONE.shiftLeft(lgn / b + 1);
            while (lowa.compareTo(higha.decrement()) < 0) {
                BigInteger mida = (lowa.add(higha)).shiftRight(1);
                BigInteger ab = pow(mida, b);
                if (ab.compareTo(n) > 0)
                    higha = mida;
                else if (ab.compareTo(n) < 0)
                    lowa = mida;
                else {
                    BigInteger[] ipp = perfectPowerDecomposition(mida);
                    if (ipp != null)
                        return new BigInteger[]{ipp[0], ipp[1].multiply(b)};
                    else
                        return new BigInteger[]{mida, BigInteger.valueOf(b)};
                }
            }
        }
        return null;
    }
}
