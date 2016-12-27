package cc.r2.core.polynomial;

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
}
