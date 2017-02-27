package cc.r2.core.poly2;


import cc.r2.core.poly2.DivisionWithRemainder.*;

import static cc.r2.core.poly2.DivisionWithRemainder.*;

/**
 * Helper methods for polynomial arithmetics.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class PolynomialArithmetics {
    private PolynomialArithmetics() {}

    /**
     * Returns the remainder of dividing {@code dividend} by {@code polyModulus}.
     *
     * @param dividend    the polynomial
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code dividend}; if not, the result will be placed directly to
     *                    {@code dividend} and the original {@code dividend} data will be lost
     * @return {@code dividend % polyModulus}
     */
    public static MutablePolynomialMod polyMod(MutablePolynomialMod dividend,
                                               MutablePolynomialMod polyModulus, boolean copy) {
        return remainder(dividend, polyModulus, copy);
    }


    /**
     * Returns the remainder of dividing the product {@code (m1 * m2)} by {@code polyModulus}.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 * m2) % polyModulus}
     */
    public static MutablePolynomialMod polyMultiplyMod(MutablePolynomialMod m1, MutablePolynomialMod m2,
                                                       MutablePolynomialMod polyModulus, boolean copy) {
        return polyMod((copy ? m1.clone() : m1).multiply(m2), polyModulus, false);
    }

    /**
     * Returns the remainder of dividing {@code dividend} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param dividend    the polynomial
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)} )})
     * @param copy        whether to clone {@code dividend}; if not, the result will be placed directly to
     *                    {@code dividend} and the original {@code dividend} data will be lost
     * @return {@code dividend % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod polyMod(MutablePolynomialMod dividend, MutablePolynomialMod polyModulus,
                                               InverseModMonomial invMod, boolean copy) {
        return remainderFast(dividend, polyModulus, invMod, copy);
    }

    /**
     * Returns the remainder of dividing the product {@code (m1 * m2)} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)} )})
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 * m2) % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod polyMultiplyMod(MutablePolynomialMod m1, MutablePolynomialMod m2,
                                                       MutablePolynomialMod polyModulus, InverseModMonomial invMod,
                                                       boolean copy) {
        return polyMod((copy ? m1.clone() : m1).multiply(m2), polyModulus, invMod, false);
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e}
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
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code polyModulus}
     *
     * @param base        the base
     * @param exponent    the non-negative exponent
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e} modulo {@code polyModulus}
     */
    public static MutablePolynomialMod polyPowMod(final MutablePolynomialMod base, long exponent,
                                                  MutablePolynomialMod polyModulus,
                                                  boolean copy) {
        return polyPowMod(base, exponent, polyModulus, fastDivisionPreConditioning(polyModulus), copy);
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code polyModulus}
     *
     * @param base        the base
     * @param exponent    the non-negative exponent
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)} )})
     * @param copy        whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e} modulo {@code polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod polyPowMod(final MutablePolynomialMod base, long exponent,
                                                  MutablePolynomialMod polyModulus, InverseModMonomial invMod,
                                                  boolean copy) {
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


    /** switch between plain and log2 algorithms */
    private static final long MONOMIAL_MOD_EXPONENT_THRESHOLD = 64;

    /**
     * Creates {@code x^exponent mod polyModulus}.
     *
     * @param exponent    the monomial exponent
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(MutablePolynomialMod)} )})
     * @return {@code x^exponent mod polyModulus}
     */
    public static MutablePolynomialMod createMonomialMod(long exponent,
                                                         MutablePolynomialMod polyModulus,
                                                         InverseModMonomial invMod) {
        if (exponent < 0)
            throw new IllegalArgumentException("Negative exponent: " + exponent);

        if (exponent == 0)
            return polyModulus.createOne();

        if (exponent < MONOMIAL_MOD_EXPONENT_THRESHOLD)
            return smallMonomial(exponent, polyModulus, invMod);
        else
            return largeMonomial(exponent, polyModulus, invMod);
    }

    /** plain create and reduce */
    static MutablePolynomialMod smallMonomial(long exponent, MutablePolynomialMod polyModulus, InverseModMonomial invMod) {
        return PolynomialArithmetics.polyMod(MutablePolynomialMod.createMonomial(polyModulus.modulus, 1, LongArithmetics.safeToInt(exponent)), polyModulus, invMod, false);
    }

    /** repeated squaring */
    static MutablePolynomialMod largeMonomial(long exponent, MutablePolynomialMod polyModulus, InverseModMonomial invMod) {
        MutablePolynomialMod base = PolynomialArithmetics.polyMod(
                MutablePolynomialMod.createMonomial(polyModulus.modulus, 1, LongArithmetics.safeToInt(MONOMIAL_MOD_EXPONENT_THRESHOLD)),
                polyModulus, invMod, false);

        MutablePolynomialMod result = base.clone();
        long exp = MONOMIAL_MOD_EXPONENT_THRESHOLD;
        for (; ; ) {
            if (exp + exp > exponent)
                break;
            result = PolynomialArithmetics.polyMultiplyMod(result, result, polyModulus, invMod, false);
            exp += exp;
        }

        MutablePolynomialMod rest = createMonomialMod(exponent - exp, polyModulus, invMod);
        return PolynomialArithmetics.polyMultiplyMod(result, rest, polyModulus, invMod, false);
    }
}
