package cc.r2.core.number;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class BigIntegerMath {
    public BigIntegerMath() {}

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
}
