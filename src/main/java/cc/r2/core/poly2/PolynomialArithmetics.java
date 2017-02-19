package cc.r2.core.poly2;


import cc.r2.core.poly2.DivideAndRemainder.InverseModMonomial;

import static cc.r2.core.poly2.DivideAndRemainder.remainder;
import static cc.r2.core.poly2.DivideAndRemainder.remainderFastWithSwitch;

/**
 * Created by poslavsky on 15/02/2017.
 */
public final class PolynomialArithmetics {
    private PolynomialArithmetics() {}

    @SuppressWarnings("ConstantConditions")
    public static MutablePolynomialMod polyMod(MutablePolynomialMod a, MutablePolynomialMod polyModulus, boolean copy) {
        return remainder(a, polyModulus, copy);
    }


    public static MutablePolynomialMod polyMultiplyMod(MutablePolynomialMod a, MutablePolynomialMod b, MutablePolynomialMod polyModulus, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b), polyModulus, false);
    }

    @SuppressWarnings("ConstantConditions")
    public static MutablePolynomialMod polyMod(MutablePolynomialMod a, MutablePolynomialMod polyModulus, InverseModMonomial invMod, boolean copy) {
        return remainderFastWithSwitch(a, polyModulus, invMod, copy);
    }

    public static MutablePolynomialMod polyMultiplyMod(MutablePolynomialMod a, MutablePolynomialMod b, MutablePolynomialMod polyModulus, InverseModMonomial invMod, boolean copy) {
        return polyMod((copy ? a.clone() : a).multiply(b), polyModulus, invMod, false);
    }


    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static <T extends MutablePolynomialAbstract<T>> T polyPow(final T base, long exponent, boolean copy) {
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

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if the result overflows a long
     */
    public static MutablePolynomialMod polyPowMod(final MutablePolynomialMod base, long exponent, MutablePolynomialMod polyModulus, InverseModMonomial invMod, boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return base.createOne();

        MutablePolynomialMod result = base.createOne();
        MutablePolynomialMod k2p = polyMod(base, polyModulus, invMod, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = polyMod(result.multiply(k2p), polyModulus, invMod, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = polyMod(k2p.multiply(k2p), polyModulus, invMod, false);
        }
    }
}
