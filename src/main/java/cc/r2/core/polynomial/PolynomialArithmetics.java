package cc.r2.core.polynomial;

import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;

import static cc.r2.core.polynomial.LongArithmetics.multiply;
import static cc.r2.core.polynomial.DivideAndRemainder.remainder;
import static cc.r2.core.polynomial.DivideAndRemainder.remainderFastWithSwitch;

final class PolynomialArithmetics {
    @SuppressWarnings("ConstantConditions")
    public static MutablePolynomial polyMod(MutablePolynomial a, MutablePolynomial polyModulus, long modulus, boolean copy) {
        return remainder(a, polyModulus, modulus, copy);
    }

    @SuppressWarnings("ConstantConditions")
    public static MutablePolynomial polyMod(MutablePolynomial a, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return remainderFastWithSwitch(a, polyModulus, invMod, modulus, copy);
    }

    public static MutablePolynomial polyMultiplyMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b, modulus), polyModulus, modulus, false);
    }

    public static MutablePolynomial polyMultiplyMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b, modulus), polyModulus, invMod, modulus, false);
    }

    public static MutablePolynomial polyAddMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).add(b, modulus), polyModulus, modulus, false);
    }

    public static MutablePolynomial polyAddMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).add(b, modulus), polyModulus, invMod, modulus, false);
    }

    public static MutablePolynomial polySubtractMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, long modulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).subtract(b, modulus), polyModulus, modulus, false);
    }

    public static MutablePolynomial polySubtractMod(MutablePolynomial a, MutablePolynomial b, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
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
    public static MutablePolynomial polyPow(final MutablePolynomial base, long exponent, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        MutablePolynomial result = MutablePolynomial.one();
        MutablePolynomial k2p = copy ? base.clone() : base;
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
    public static MutablePolynomial polyPowMod(final MutablePolynomial base, long exponent, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        if (exponent == modulus) // fast exponent
            return powMod0(base, (int) exponent);

        MutablePolynomial result = MutablePolynomial.one();
        MutablePolynomial k2p = (copy ? base.clone() : base).modulus(modulus);
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
    public static MutablePolynomial polyPowMod(final MutablePolynomial base, long exponent, MutablePolynomial polyModulus, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return MutablePolynomial.one();

        MutablePolynomial result = MutablePolynomial.one();
        MutablePolynomial k2p = polyMod(base, polyModulus, modulus, copy); // this will copy the base
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
    public static MutablePolynomial polyPowMod(final MutablePolynomial base, long exponent, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return MutablePolynomial.one();

        MutablePolynomial result = MutablePolynomial.one();
        MutablePolynomial k2p = polyMod(base, polyModulus, invMod, modulus, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = polyMod(result.multiply(k2p, modulus), polyModulus, invMod, modulus, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = polyMod(k2p.multiply(k2p, modulus), polyModulus, invMod, modulus, false);
        }
    }

    static MutablePolynomial powMod0(MutablePolynomial base, int modulus) {
        long[] result = new long[base.degree * modulus + 1];
        for (int i = base.degree; i >= 0; --i)
            result[i * modulus] = LongArithmetics.mod(base.data[i], modulus);
        return MutablePolynomial.create(result);
    }
}
