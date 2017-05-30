package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class CommonPolynomialsArithmetics {
    private CommonPolynomialsArithmetics() {}

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e}
     */
    public static <T extends IGeneralPolynomial<T>> T polyPow(final T base, long exponent, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        T result = base.createOne();
        T k2p = copy ? base.clone() : base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = result.multiply(k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p);
        }
    }
}
