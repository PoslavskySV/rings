package cc.r2.core.poly.univar2;


import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.LongArithmetics;

import static cc.r2.core.poly.univar2.DivisionWithRemainder.remainder;

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
    public static <T extends IUnivariatePolynomial<T>> T polyMod(T dividend, T polyModulus, boolean copy) {
        return remainder(dividend, polyModulus, copy);
    }

    /**
     * Returns the remainder of dividing {@code dividend} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param dividend    the polynomial
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code dividend}; if not, the result will be placed directly to
     *                    {@code dividend} and the original {@code dividend} data will be lost
     * @return {@code dividend % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyMod(T dividend, T polyModulus,
                                                                 DivisionWithRemainder.InverseModMonomial<T> invMod, boolean copy) {
        return DivisionWithRemainder.remainderFast(dividend, polyModulus, invMod, copy);
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
    public static <T extends IUnivariatePolynomial<T>> T polyMultiplyMod(T m1, T m2, T polyModulus, boolean copy) {
        return polyMod((copy ? m1.clone() : m1).multiply(m2), polyModulus, false);
    }

    /**
     * Returns the remainder of dividing the product {@code (m1 * m2)} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 * m2) % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyMultiplyMod(T m1, T m2,
                                                                         T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                         boolean copy) {
        return polyMod((copy ? m1.clone() : m1).multiply(m2), polyModulus, invMod, false);
    }

    /**
     * Returns the remainder of dividing the sum {@code (m1 + m2)} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 + m2) % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyAddMod(T m1, T m2,
                                                                    T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                    boolean copy) {
        return polyMod((copy ? m1.clone() : m1).add(m2), polyModulus, invMod, false);
    }

    /**
     * Returns the remainder of dividing the sum {@code (m1 + m2)} by {@code polyModulus}.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 + m2) % polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T polyAddMod(T m1, T m2, T polyModulus, boolean copy) {
        return polyMod((copy ? m1.clone() : m1).add(m2), polyModulus, false);
    }

    /**
     * Returns the remainder of dividing the difference {@code (m1 - m2)} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 - m2) % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polySubtractMod(T m1, T m2,
                                                                         T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                         boolean copy) {
        return polyMod((copy ? m1.clone() : m1).subtract(m2), polyModulus, invMod, false);
    }

    /**
     * Returns the remainder of dividing the difference {@code (m1 - m2)} by {@code polyModulus}.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 - m2) % polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T polySubtractMod(T m1, T m2, T polyModulus, boolean copy) {
        return polyMod((copy ? m1.clone() : m1).subtract(m2), polyModulus, false);
    }

    /**
     * Returns the remainder of dividing the negated poly {@code -m1} by {@code polyModulus} using fast algorithm for
     * pre-conditioned modulus.
     *
     * @param m1          the polynomial
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (-m1) % polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyNegateMod(T m1,
                                                                       T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                       boolean copy) {
        return polyMod((copy ? m1.clone() : m1).negate(), polyModulus, invMod, false);
    }

    /**
     * Returns the remainder of dividing the negated poly {@code -m1} by {@code polyModulus}.
     *
     * @param m1          the polynomial
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (-m1) % polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T polyNegateMod(T m1, T polyModulus, boolean copy) {
        return polyMod((copy ? m1.clone() : m1).negate(), polyModulus, false);
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param copy     whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e}
     */
    public static <T extends IUnivariatePolynomial<T>> T polyPow(final T base, long exponent, boolean copy) {
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
    public static <T extends IUnivariatePolynomial<T>> T polyPowMod(final T base, long exponent,
                                                                    T polyModulus,
                                                                    boolean copy) {
        return polyPowMod(base, exponent, polyModulus, DivisionWithRemainder.fastDivisionPreConditioning(polyModulus), copy);
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code polyModulus}
     *
     * @param base        the base
     * @param exponent    the non-negative exponent
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e} modulo {@code polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyPowMod(final T base, long exponent,
                                                                    T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                    boolean copy) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return base.createOne();

        T result = base.createOne();
        T k2p = polyMod(base, polyModulus, invMod, copy); // this will copy the base
        for (; ; ) {
            if ((exponent&1) != 0)
                result = polyMod(result.multiply(k2p), polyModulus, invMod, false);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = polyMod(k2p.multiply(k2p), polyModulus, invMod, false);
        }
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code polyModulus}
     *
     * @param base        the base
     * @param exponent    the non-negative exponent
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @param copy        whether to clone {@code base}; if not the data of {@code base} will be lost
     * @return {@code base} in a power of {@code e} modulo {@code polyModulus}
     * @see DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)
     */
    public static <T extends IUnivariatePolynomial<T>> T polyPowMod(final T base, BigInteger exponent,
                                                                    T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod,
                                                                    boolean copy) {
        if (exponent.signum() < 0)
            throw new IllegalArgumentException();
        if (exponent.isZero())
            return base.createOne();


        T result = base.createOne();
        T k2p = polyMod(base, polyModulus, invMod, copy); // this will copy the base
        for (; ; ) {
            if (exponent.testBit(0))
                result = polyMod(result.multiply(k2p), polyModulus, invMod, false);
            exponent = exponent.shiftRight(1);
            if (exponent.isZero())
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
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @return {@code x^exponent mod polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T createMonomialMod(long exponent,
                                                                           T polyModulus,
                                                                           DivisionWithRemainder.InverseModMonomial<T> invMod) {
        if (exponent < 0)
            throw new IllegalArgumentException("Negative exponent: " + exponent);

        if (exponent == 0)
            return polyModulus.createOne();

        if (exponent < MONOMIAL_MOD_EXPONENT_THRESHOLD)
            return smallMonomial(exponent, polyModulus, invMod);
        else
            return largeMonomial(exponent, polyModulus, invMod);
    }

    /**
     * Creates {@code x^exponent mod polyModulus}.
     *
     * @param exponent    the monomial exponent
     * @param polyModulus the modulus
     * @param invMod      pre-conditioned modulus ({@link DivisionWithRemainder#fastDivisionPreConditioning(IUnivariatePolynomial)} )})
     * @return {@code x^exponent mod polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T createMonomialMod(BigInteger exponent,
                                                                           T polyModulus,
                                                                           DivisionWithRemainder.InverseModMonomial<T> invMod) {
        if (exponent.signum() < 0)
            throw new IllegalArgumentException("Negative exponent: " + exponent);

        if (exponent.isZero())
            return polyModulus.createOne();

        if (exponent.isLong())
            return createMonomialMod(exponent.longValueExact(), polyModulus, invMod);
        else
            return largeMonomial(exponent, polyModulus, invMod);
    }

    /** plain create and reduce */
    static <T extends IUnivariatePolynomial<T>> T smallMonomial(long exponent, T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod) {
        return PolynomialArithmetics.polyMod(polyModulus.createMonomial(LongArithmetics.safeToInt(exponent)), polyModulus, invMod, false);
    }

    /** repeated squaring */
    static <T extends IUnivariatePolynomial<T>> T largeMonomial(long exponent, T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod) {
        T base = PolynomialArithmetics.polyMod(
                polyModulus.createMonomial(LongArithmetics.safeToInt(MONOMIAL_MOD_EXPONENT_THRESHOLD)),
                polyModulus, invMod, false);

        T result = base.clone();
        long exp = MONOMIAL_MOD_EXPONENT_THRESHOLD;
        for (; ; ) {
            if (exp + exp > exponent)
                break;
            result = PolynomialArithmetics.polyMultiplyMod(result, result, polyModulus, invMod, false);
            exp += exp;
        }

        T rest = createMonomialMod(exponent - exp, polyModulus, invMod);
        return PolynomialArithmetics.polyMultiplyMod(result, rest, polyModulus, invMod, false);
    }

    /** repeated squaring */
    static <T extends IUnivariatePolynomial<T>> T largeMonomial(BigInteger exponent, T polyModulus, DivisionWithRemainder.InverseModMonomial<T> invMod) {
        T base = PolynomialArithmetics.polyMod(
                polyModulus.createMonomial(LongArithmetics.safeToInt(MONOMIAL_MOD_EXPONENT_THRESHOLD)),
                polyModulus, invMod, false);

        T result = base.clone();
        BigInteger exp = BigInteger.valueOf(MONOMIAL_MOD_EXPONENT_THRESHOLD);
        for (; ; ) {
            if (exp.shiftLeft(1).compareTo(exponent) > 0)
                break;
            result = PolynomialArithmetics.polyMultiplyMod(result, result, polyModulus, invMod, false);
            exp = exp.shiftLeft(1);
        }

        T rest = createMonomialMod(exponent.subtract(exp), polyModulus, invMod);
        return PolynomialArithmetics.polyMultiplyMod(result, rest, polyModulus, invMod, false);
    }
}
