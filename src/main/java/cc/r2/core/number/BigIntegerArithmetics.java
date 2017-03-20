package cc.r2.core.number;

import static cc.r2.core.number.BigInteger.ONE;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class BigIntegerArithmetics {
    public BigIntegerArithmetics() {}

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
    public static BigInteger pow(final BigInteger base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        BigInteger result = ONE;
        BigInteger k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
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

    /* ************************ mock methods for @Specialization ************************ */

    public static BigInteger safeAdd(BigInteger a, BigInteger b) {
        return a.add(b);
    }

    public static BigInteger safeSubtract(BigInteger a, BigInteger b) {
        return a.subtract(b);
    }

    public static BigInteger safeMultiply(BigInteger a, BigInteger b) {
        return a.multiply(b);
    }

    public static BigInteger safePow(BigInteger a, long exp) {
        return pow(a, exp);
    }

    public static BigInteger modInverse(BigInteger a, BigInteger mod) {
        return a.modInverse(mod);
    }

    public static BigInteger mod(BigInteger a, BigInteger mod) {
        return a.mod(mod);
    }
}
