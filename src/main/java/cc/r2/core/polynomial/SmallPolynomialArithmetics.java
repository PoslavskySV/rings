package cc.r2.core.polynomial;

import cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.InverseModMonomial;

import static cc.r2.core.polynomial.LongArithmetics.multiply;
import static cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.remainder;
import static cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.remainderFastWithSwitch;

final class SmallPolynomialArithmetics {
    @SuppressWarnings("ConstantConditions")
    public static MutableLongPoly polyMod(MutableLongPoly a, MutableLongPoly polyModulus, long modulus, boolean copy) {
        return remainder(a, polyModulus, modulus, copy);
    }

    @SuppressWarnings("ConstantConditions")
    public static MutableLongPoly polyMod(MutableLongPoly a, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return remainderFastWithSwitch(a, polyModulus, invMod, modulus, copy);
    }

    public static MutableLongPoly polyMultiplyMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b, modulus), polyModulus, modulus, false);
    }

    public static MutableLongPoly polyMultiplyMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b, modulus), polyModulus, invMod, modulus, false);
    }

    public static MutableLongPoly polyAddMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).add(b, modulus), polyModulus, modulus, false);
    }

    public static MutableLongPoly polyAddMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).add(b, modulus), polyModulus, invMod, modulus, false);
    }

    public static MutableLongPoly polySubtractMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).subtract(b, modulus), polyModulus, modulus, false);
    }

    public static MutableLongPoly polySubtractMod(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).subtract(b, modulus), polyModulus, invMod, modulus, false);
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static MutableLongPoly polyPow(final MutableLongPoly base, long exponent, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = copy ? base.clone() : base;
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
    public static MutableLongPoly polyPowMod(final MutableLongPoly base, long exponent, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        if (exponent == modulus) // fast exponent
            return powMod0(base, (int) exponent);

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = (copy ? base.clone() : base).modulus(modulus);
        for (; ; ) {
            if ((exponent&1) != 0)
                result = result.multiply(k2p, modulus);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = k2p.multiply(k2p, modulus);
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
    public static MutableLongPoly polyPowMod(final MutableLongPoly base, long exponent, MutableLongPoly polyModulus, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return MutableLongPoly.one();

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = polyMod(base, polyModulus, modulus, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = polyMod(result.multiply(k2p, modulus), polyModulus, modulus, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = polyMod(k2p.multiply(k2p, modulus), polyModulus, modulus, false);
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
    public static MutableLongPoly polyPowMod(final MutableLongPoly base, long exponent, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return MutableLongPoly.one();

        MutableLongPoly result = MutableLongPoly.one();
        MutableLongPoly k2p = polyMod(base, polyModulus, invMod, modulus, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = polyMod(result.multiply(k2p, modulus), polyModulus, invMod, modulus, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = polyMod(k2p.multiply(k2p, modulus), polyModulus, invMod, modulus, false);
        }
    }

    static MutableLongPoly powMod0(MutableLongPoly base, int modulus) {
        long[] result = new long[base.degree * modulus + 1];
        for (int i = base.degree; i >= 0; --i)
            result[i * modulus] = LongArithmetics.mod(base.data[i], modulus);
        return MutableLongPoly.create(result);
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
