package cc.r2.core.poly;

import gnu.trove.map.hash.TIntObjectHashMap;

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
        if (exponent == 1)
            return copy ? base.clone() : base;
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

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @param cache    cache to store all intermediate powers
     * @return {@code base} in a power of {@code e}
     */
    public static <T extends IGeneralPolynomial<T>> T polyPow(final T base, int exponent, boolean copy,
                                                              TIntObjectHashMap<T> cache) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 1)
            return copy ? base.clone() : base;

        T cached = cache.get(exponent);
        if (cached != null)
            return cached.clone();

        T result = base.createOne();
        T k2p = copy ? base.clone() : base;
        int rExp = 0, kExp = 1;
        for (; ; ) {
            if ((exponent&1) != 0)
                cache.put(rExp += kExp, result.multiply(k2p).clone());
            exponent = exponent >> 1;
            if (exponent == 0) {
                cache.put(rExp, result);
                return result;
            }
            cache.put(kExp *= 2, k2p.square().clone());
        }
    }
}
