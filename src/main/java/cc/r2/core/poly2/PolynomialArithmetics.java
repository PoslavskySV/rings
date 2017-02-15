package cc.r2.core.poly2;

import static cc.r2.core.poly2.DivideAndRemainder.remainder;

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
}
