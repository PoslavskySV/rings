package cc.r2.core.polynomial;

import static cc.r2.core.polynomial.LongArithmetics.multiply;

final class SmallPolynomialArithmetics {
    @SuppressWarnings("ConstantConditions")
    public static MutableLongPoly mod(MutableLongPoly a, MutableLongPoly polyModulus, long modulus) {
        if (a.degree < polyModulus.degree)
            return a.modulus(modulus);
        return SmallPolynomials.divideAndRemainder(a, polyModulus, modulus)[1];
    }

    public static MutableLongPoly multiplyMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus) {
        return mod(a.clone().multiply(b, modulus), polyModulus, modulus);
    }

    public static MutableLongPoly addMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus) {
        return mod(a.clone().add(b, modulus), polyModulus, modulus);
    }

    public static MutableLongPoly subtractMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus) {
        return mod(a.clone().subtract(b, modulus), polyModulus, modulus);
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static MutableLongPoly pow(final MutableLongPoly base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = base;
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
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static MutableLongPoly powMod(final MutableLongPoly base, long exponent, long modulus) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = result.multiply(k2p, modulus);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p, modulus);
        }
    }

    public static MutableLongPoly derivative(MutableLongPoly poly) {
        if (poly.isConstant())
            return MutableLongPoly.zero();
        long[] data = new long[poly.degree];
        for (int i = poly.degree; i > 0; --i)
            data[i - 1] = multiply(poly.data[i], i);
        return MutableLongPoly.create(data);
    }

    public static MutableLongPoly derivative(MutableLongPoly poly, long modulus) {
        if (poly.isConstant())
            return MutableLongPoly.zero();
        long[] data = new long[poly.degree];
        for (int i = poly.degree; i > 0; --i)
            data[i - 1] = LongArithmetics.multiplyMod(poly.data[i], i, modulus);
        return MutableLongPoly.create(data);
    }
}
