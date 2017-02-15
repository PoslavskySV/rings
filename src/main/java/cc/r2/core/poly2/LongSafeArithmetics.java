package cc.r2.core.poly2;

/**
 * Created by poslavsky on 15/02/2017.
 */
public final class LongSafeArithmetics {
    private LongSafeArithmetics() {}

    /**
     * Delegates to {@link Math#multiplyExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long multiply(long x, long y) {
        return Math.multiplyExact(x, y);
    }

    /**
     * Delegates to {@link Math#multiplyExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long multiply(long x, long y, long z) {
        return Math.multiplyExact(Math.multiplyExact(x, y), z);
    }

    /**
     * Delegates to {@link Math#addExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long add(long x, long y) {
        return Math.addExact(x, y);
    }

    /**
     * Delegates to {@link Math#subtractExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long subtract(long a, long b) {
        return Math.subtractExact(a, b);
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static long pow(final long base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        long result = 1L;
        long k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = multiply(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiply(k2p, k2p);
        }
    }

}
