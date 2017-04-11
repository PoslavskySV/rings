package cc.r2.core.poly.univar;


import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.ArrayList;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.ZERO;
import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Algorithms for division with remainder.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class DivisionWithRemainder {
    private DivisionWithRemainder() {}

    /* **************************************** Common methods  *************************************** */

    /**
     * Returns the remainder of {@code dividend} divided by monomial {@code x^xDegree}
     *
     * @param dividend the dividend
     * @param xDegree  monomial degree
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static <T extends IMutablePolynomial<T>> T remainderMonomial(T dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).truncate(xDegree - 1);
    }

    private static void checkZeroDivider(IMutablePolynomial p) {
        if (p.isZero())
            throw new ArithmeticException("divide by zero");
    }

    /* ************************************ Machine-precision division in Z[x]  ************************************ */

    /**
     * Returns {@code {quotient, remainder}} or {@code null} if the division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    public static lMutablePolynomialZ[] divideAndRemainder(final lMutablePolynomialZ dividend,
                                                           final lMutablePolynomialZ divider,
                                                           boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), lMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            lMutablePolynomialZ div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new lMutablePolynomialZ[]{div, lMutablePolynomialZ.zero()};
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
        return divideAndRemainderGeneral0(dividend, divider, 1, copy);
    }

    /**
     * Returns quotient and remainder using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static lMutablePolynomialZ[] pseudoDivideAndRemainder(final lMutablePolynomialZ dividend,
                                                                 final lMutablePolynomialZ divider,
                                                                 final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), lMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        long factor = LongArithmetics.safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new lMutablePolynomialZ[]{(copy ? dividend.clone() : dividend).multiply(factor / dividend.lc()), lMutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderGeneral0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static lMutablePolynomialZ[] divideAndRemainderGeneral0(final lMutablePolynomialZ dividend,
                                                            final lMutablePolynomialZ divider,
                                                            final long dividendRaiseFactor,
                                                            final boolean copy) {
        assert dividend.degree >= divider.degree;
        if (divider.lc() == 1 && dividendRaiseFactor == 1)
            return divideAndRemainderGeneralMonic(dividend, divider, copy);

        lMutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        long[] quotient = new long[dividend.degree - divider.degree + 1];


        Magic magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc())
                    return null;

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), remainder};
    }

    /** Plain school implementation */
    private static lMutablePolynomialZ[] divideAndRemainderGeneralMonic(final lMutablePolynomialZ dividend,
                                                                        final lMutablePolynomialZ divider,
                                                                        final boolean copy) {
        assert divider.lc() == 1;

        lMutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend);
        long[] quotient = new long[dividend.degree - divider.degree + 1];
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.lc();
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }
        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), remainder};
    }

    /**
     * Returns quotient and remainder using adaptive pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    static lMutablePolynomialZ[] pseudoDivideAndRemainderAdaptive(final lMutablePolynomialZ dividend,
                                                                  final lMutablePolynomialZ divider,
                                                                  final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), lMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new lMutablePolynomialZ[]{copy ? dividend.clone() : dividend, lMutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static lMutablePolynomialZ[] pseudoDivideAndRemainderAdaptive0(final lMutablePolynomialZ dividend,
                                                                   final lMutablePolynomialZ divider,
                                                                   final boolean copy) {
        assert dividend.degree >= divider.degree;

        lMutablePolynomialZ remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        Magic magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc()) {
                    long gcd = LongArithmetics.gcd(remainder.lc(), divider.lc());
                    long factor = divider.lc() / gcd;
                    remainder.multiply(factor);
                    for (int j = i + 1; j < quotient.length; ++j)
                        quotient[j] = LongArithmetics.safeMultiply(quotient[j], factor);
                    quot = divideSignedFast(remainder.lc(), magic);
                }

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] pseudoDivideAndRemainderLinearDividerAdaptive(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc(), factor = 1;
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        Magic magic = magicSigned(lc);
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = LongArithmetics.safeAdd(LongArithmetics.safeMultiply(res, cc), LongArithmetics.safeMultiply(factor, tmp));
            if (i == 0) break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res) {
                long gcd = LongArithmetics.gcd(res, lc), f = lc / gcd;
                factor = LongArithmetics.safeMultiply(factor, f);
                res = LongArithmetics.safeMultiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = LongArithmetics.safeMultiply(quotient[j], f);
                quot = divideSignedFast(res, magic);
            }
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), lMutablePolynomialZ.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] divideAndRemainderLinearDivider(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] pseudoDivideAndRemainderLinearDivider(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, LongArithmetics.safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] divideAndRemainderLinearDivider0(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, long raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc();
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        Magic magic = magicSigned(lc);
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = LongArithmetics.safeAdd(LongArithmetics.safeMultiply(res, cc), LongArithmetics.safeMultiply(raiseFactor, tmp));
            if (i == 0)
                break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res)
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), lMutablePolynomialZ.create(res)};
    }

    /**
     * Returns remainder of dividing {@code dividend} by {@code divider} or {@code null} if division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder or {@code null} if division is not possible
     */
    public static lMutablePolynomialZ remainder(final lMutablePolynomialZ dividend,
                                                final lMutablePolynomialZ divider,
                                                final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return lMutablePolynomialZ.zero();
        if (divider.degree == 1) {
            if (divider.cc() % divider.lc() != 0)
                return null;
            return lMutablePolynomialZ.create(dividend.evaluate(-divider.cc() / divider.lc()));
        }
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static lMutablePolynomialZ remainder0(final lMutablePolynomialZ dividend,
                                          final lMutablePolynomialZ divider,
                                          final boolean copy) {
        assert dividend.degree >= divider.degree;

        lMutablePolynomialZ remainder = copy ? dividend.clone() : dividend;
        Magic magic = magicSigned(remainder.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc())
                    return null;
                remainder.subtract(divider, quot, i);
            }
        return remainder;
    }

    /* ************************************ Multi-precision division in Z[x]  ************************************ */

    /**
     * Returns {@code {quotient, remainder}} or {@code null} if the division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    public static bMutablePolynomialZ[] divideAndRemainder(final bMutablePolynomialZ dividend,
                                                           final bMutablePolynomialZ divider,
                                                           boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), bMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            bMutablePolynomialZ div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new bMutablePolynomialZ[]{div, bMutablePolynomialZ.zero()};
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
        return divideAndRemainderGeneral0(dividend, divider, ONE, copy);
    }

    /**
     * Returns quotient and remainder using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static bMutablePolynomialZ[] pseudoDivideAndRemainder(final bMutablePolynomialZ dividend,
                                                                 final bMutablePolynomialZ divider,
                                                                 final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), bMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        BigInteger factor = BigIntegerArithmetics.safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new bMutablePolynomialZ[]{(copy ? dividend.clone() : dividend).multiply(factor.divide(dividend.lc())), bMutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderGeneral0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static bMutablePolynomialZ[] divideAndRemainderGeneral0(final bMutablePolynomialZ dividend,
                                                            final bMutablePolynomialZ divider,
                                                            final BigInteger dividendRaiseFactor,
                                                            final boolean copy) {
        assert dividend.degree >= divider.degree;
        if (divider.lc().isOne() && dividendRaiseFactor.isOne())
            return divideAndRemainderGeneralMonic(dividend, divider, copy);

        bMutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        BigInteger[] quotient = new BigInteger[dividend.degree - divider.degree + 1];


        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                BigInteger quot = remainder.lc().divide(divider.lc());
                if (!quot.multiply(divider.lc()).equals(remainder.lc()))
                    return null;

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = ZERO;
        }

        return new bMutablePolynomialZ[]{bMutablePolynomialZ.create(quotient), remainder};
    }

    /** Plain school implementation */
    private static bMutablePolynomialZ[] divideAndRemainderGeneralMonic(final bMutablePolynomialZ dividend,
                                                                        final bMutablePolynomialZ divider,
                                                                        final boolean copy) {
        assert divider.lc().isOne();

        bMutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend);
        BigInteger[] quotient = new BigInteger[dividend.degree - divider.degree + 1];
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.lc();
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = ZERO;
        }
        return new bMutablePolynomialZ[]{bMutablePolynomialZ.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] divideAndRemainderLinearDivider(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, ONE, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] pseudoDivideAndRemainderLinearDivider(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, BigIntegerArithmetics.safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] divideAndRemainderLinearDivider0(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, BigInteger raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        BigInteger cc = divider.cc().negate(), lc = divider.lc();
        BigInteger[] quotient = copy ? new BigInteger[dividend.degree] : dividend.data;
        BigInteger res = ZERO;
        for (int i = dividend.degree; ; --i) {
            BigInteger tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = BigIntegerArithmetics.safeAdd(BigIntegerArithmetics.safeMultiply(res, cc), BigIntegerArithmetics.safeMultiply(raiseFactor, tmp));
            if (i == 0)
                break;
            BigInteger quot = res.divide(lc);
            if (!quot.multiply(lc).equals(res))
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = ZERO;
        return new bMutablePolynomialZ[]{bMutablePolynomialZ.create(quotient), bMutablePolynomialZ.create(res)};
    }

    /**
     * Returns remainder of dividing {@code dividend} by {@code divider} or {@code null} if division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder or {@code null} if division is not possible
     */
    public static bMutablePolynomialZ remainder(final bMutablePolynomialZ dividend,
                                                final bMutablePolynomialZ divider,
                                                final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return bMutablePolynomialZ.zero();
        if (divider.degree == 1) {
            if (!(divider.cc().mod(divider.lc())).isZero())
                return null;
            return bMutablePolynomialZ.create(dividend.evaluate(divider.cc().negate().divide(divider.lc())));
        }
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static bMutablePolynomialZ remainder0(final bMutablePolynomialZ dividend,
                                          final bMutablePolynomialZ divider,
                                          final boolean copy) {
        assert dividend.degree >= divider.degree;

        bMutablePolynomialZ remainder = copy ? dividend.clone() : dividend;
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i) {
                BigInteger quot = remainder.lc().divide(remainder.lc());
                if (!quot.multiply(divider.lc()).equals(remainder.lc()))
                    return null;
                remainder.subtract(divider, quot, i);
            }
        return remainder;
    }


    /* ************************************ Machine-precision division in Zp[x]  ************************************ */


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(IMutablePolynomialZp<?> dividend,
                                                IMutablePolynomialZp<?> divider) {
        // practical benchmarks show that without pre-conditioning,
        // classical division is always faster or at least the same fast
        return true;
    }

    /** early checks for division */
    private static lMutablePolynomialZp[] earlyDivideAndRemainderChecks(final lMutablePolynomialZp dividend,
                                                                        final lMutablePolynomialZp divider,
                                                                        final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lMutablePolynomialZp[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new lMutablePolynomialZp[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new lMutablePolynomialZp[]{(copy ? dividend.clone() : dividend).multiply(LongArithmetics.modInverse(divider.lc(), dividend.modulus)), dividend.createZero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, copy);
        return null;
    }

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static lMutablePolynomialZp[] divideAndRemainder(final lMutablePolynomialZp dividend,
                                                            final lMutablePolynomialZp divider,
                                                            final boolean copy) {
        lMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy);
    }

    /**
     * Classical algorithm for division with remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static lMutablePolynomialZp[] divideAndRemainderClassic(final lMutablePolynomialZp dividend,
                                                                   final lMutablePolynomialZp divider,
                                                                   final boolean copy) {
        lMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static lMutablePolynomialZp[] divideAndRemainderClassic0(final lMutablePolynomialZp dividend,
                                                             final lMutablePolynomialZp divider,
                                                             final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        lMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        long lcInverse = LongArithmetics.modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.multiplyMod(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new lMutablePolynomialZp[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZp[] divideAndRemainderLinearDividerModulus(lMutablePolynomialZp dividend, lMutablePolynomialZp divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatibleModulus(divider);

        //apply Horner's method

        long cc = LongArithmetics.mod(-divider.cc(), dividend.modulus);
        long lcInverse = LongArithmetics.modInverse(divider.lc(), dividend.modulus);

        if (divider.lc() != 1)
            cc = dividend.multiplyMod(cc, lcInverse);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = dividend.multiplyMod(res, lcInverse);
            res = dividend.addMod(dividend.multiplyMod(res, cc), tmp);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lMutablePolynomialZp[]{dividend.createFromArray(quotient), dividend.createFromArray(new long[]{res})};
    }

    /* ************************************ Multi-precision division in Zp[x]  ************************************ */

    /** early checks for division */
    private static bMutablePolynomialZp[] earlyDivideAndRemainderChecks(final bMutablePolynomialZp dividend,
                                                                        final bMutablePolynomialZp divider,
                                                                        final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new bMutablePolynomialZp[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new bMutablePolynomialZp[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new bMutablePolynomialZp[]{(copy ? dividend.clone() : dividend).multiply(BigIntegerArithmetics.modInverse(divider.lc(), dividend.modulus)), dividend.createZero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, copy);
        return null;
    }

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static bMutablePolynomialZp[] divideAndRemainder(final bMutablePolynomialZp dividend,
                                                            final bMutablePolynomialZp divider,
                                                            final boolean copy) {
        bMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy);
    }

    /**
     * Classical algorithm for division with remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static bMutablePolynomialZp[] divideAndRemainderClassic(final bMutablePolynomialZp dividend,
                                                                   final bMutablePolynomialZp divider,
                                                                   final boolean copy) {
        bMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static bMutablePolynomialZp[] divideAndRemainderClassic0(final bMutablePolynomialZp dividend,
                                                             final bMutablePolynomialZp divider,
                                                             final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        bMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        BigInteger[] quotient = new BigInteger[dividend.degree - divider.degree + 1];

        BigInteger lcInverse = BigIntegerArithmetics.modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.multiplyMod(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = ZERO;
        }

        return new bMutablePolynomialZp[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZp[] divideAndRemainderLinearDividerModulus(bMutablePolynomialZp dividend, bMutablePolynomialZp divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatibleModulus(divider);

        //apply Horner's method

        BigInteger cc = BigIntegerArithmetics.mod(divider.cc().negate(), dividend.modulus);
        BigInteger lcInverse = BigIntegerArithmetics.modInverse(divider.lc(), dividend.modulus);

        if (!divider.lc().isOne())
            cc = dividend.multiplyMod(cc, lcInverse);

        BigInteger[] quotient = copy ? new BigInteger[dividend.degree] : dividend.data;
        BigInteger res = ZERO;
        for (int i = dividend.degree; i >= 0; --i) {
            BigInteger tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = dividend.multiplyMod(res, lcInverse);
            res = dividend.addMod(dividend.multiplyMod(res, cc), tmp);
        }
        if (!copy) quotient[dividend.degree] = ZERO;
        return new bMutablePolynomialZp[]{dividend.createFromArray(quotient), dividend.createFromArray(new BigInteger[]{res})};
    }


    /* ********************************** Fast division algorithm ********************************** */

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomial<T extends IMutablePolynomial<T>> {
        final T poly;

        private InverseModMonomial(T poly) {
            if (!poly.isUnitCC())
                throw new IllegalArgumentException("Smallest coefficient is not a unit: " + poly);
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<T> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public T getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() >= r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            T gPrev = currentSize == 0 ? poly.createOne() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                T tmp = gPrev.clone().multiply(2).subtract(gPrev.clone().square().multiply(poly));
                inverses.add(gPrev = remainderMonomial(tmp, 1 << i, false));
            }
            return gPrev;
        }
    }

    /**
     * Prepares {@code rev(divider)^(-1) mod x^i } for fast division.
     *
     * @param divider the divider
     */
    public static <T extends IMutablePolynomial<T>> InverseModMonomial<T> fastDivisionPreConditioning(T divider) {
        if (!divider.isMonic())
            throw new IllegalArgumentException("Only monic polynomials allowed. Input: " + divider);
        return new InverseModMonomial<>(divider.clone().reverse());
    }

    /** fast division implementation */
    public static <T extends IMutablePolynomial<T>> T[] divideAndRemainderFast0(T dividend, T divider,
                                                                                InverseModMonomial<T> invRevMod,
                                                                                boolean copy) {
        int m = dividend.degree() - divider.degree();
        T q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1)), m + 1, false).reverse();
        if (q.degree() < m)
            q.shiftRight(m - q.degree());
        T r = (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(q));
        return dividend.arrayNewInstance(q, r);
    }

    /* ********************************* Machine-precision fast division in Zp[x]  ******************************** */

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static lMutablePolynomialZp[] divideAndRemainderFast(lMutablePolynomialZp dividend,
                                                                lMutablePolynomialZp divider,
                                                                boolean copy) {
        lMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, copy);
    }

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static lMutablePolynomialZp[] divideAndRemainderFast(lMutablePolynomialZp dividend,
                                                                lMutablePolynomialZp divider,
                                                                InverseModMonomial<lMutablePolynomialZp> invMod,
                                                                boolean copy) {
        lMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invMod, copy);
    }

    static lMutablePolynomialZp[] divideAndRemainderFast0(lMutablePolynomialZp dividend,
                                                          lMutablePolynomialZp divider,
                                                          boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        long lc = divider.lc();
        long lcInv = LongArithmetics.modInverse(lc, dividend.modulus);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        lMutablePolynomialZp[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    /* ********************************* Multi-precision fast division in Zp[x]  ******************************** */

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static bMutablePolynomialZp[] divideAndRemainderFast(bMutablePolynomialZp dividend,
                                                                bMutablePolynomialZp divider,
                                                                boolean copy) {
        bMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, copy);
    }

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static bMutablePolynomialZp[] divideAndRemainderFast(bMutablePolynomialZp dividend,
                                                                bMutablePolynomialZp divider,
                                                                InverseModMonomial<bMutablePolynomialZp> invMod,
                                                                boolean copy) {
        bMutablePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invMod, copy);
    }

    static bMutablePolynomialZp[] divideAndRemainderFast0(bMutablePolynomialZp dividend,
                                                          bMutablePolynomialZp divider,
                                                          boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.lc().isOne())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        BigInteger lc = divider.lc();
        BigInteger lcInv = BigIntegerArithmetics.modInverse(lc, dividend.modulus);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        bMutablePolynomialZp[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    /* ********************************* Remainder in Zp[x]  ******************************** */

    /** fast division checks */
    private static lMutablePolynomialZp earlyRemainderChecks(final lMutablePolynomialZp dividend,
                                                             final lMutablePolynomialZp divider,
                                                             final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1)
            return dividend.createFromArray(new long[]{
                    dividend.evaluate(
                            dividend.multiplyMod(LongArithmetics.mod(-divider.cc(), dividend.modulus), LongArithmetics.modInverse(divider.lc(), dividend.modulus)))
            });
        return null;
    }

    /**
     * Returns remainder of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static lMutablePolynomialZp remainder(final lMutablePolynomialZp dividend,
                                                 final lMutablePolynomialZp divider,
                                                 final boolean copy) {
        lMutablePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static lMutablePolynomialZp remainderClassical0(final lMutablePolynomialZp dividend,
                                                    final lMutablePolynomialZp divider,
                                                    final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        lMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long lcInverse = LongArithmetics.modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, remainder.multiplyMod(remainder.lc(), lcInverse), i);

        return remainder;
    }

    /**
     * Fast remainder using Newton's iteration with switch to classical remainder for small polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static lMutablePolynomialZp remainderFast(final lMutablePolynomialZp dividend,
                                                     final lMutablePolynomialZp divider,
                                                     final InverseModMonomial<lMutablePolynomialZp> invMod,
                                                     final boolean copy) {
        lMutablePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[1];
    }

    /**
     * Returns quotient of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static lMutablePolynomialZp quotient(final lMutablePolynomialZp dividend,
                                                final lMutablePolynomialZp divider,
                                                final boolean copy) {
        lMutablePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic(dividend, divider, copy)[0];

        return divideAndRemainderFast0(dividend, divider, copy)[0];
    }

    /**
     * Fast quotient using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static lMutablePolynomialZp quotientFast(final lMutablePolynomialZp dividend,
                                                    final lMutablePolynomialZp divider,
                                                    final InverseModMonomial<lMutablePolynomialZp> invMod,
                                                    final boolean copy) {
        lMutablePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[0];
    }


    /** fast division checks */
    private static bMutablePolynomialZp earlyRemainderChecks(final bMutablePolynomialZp dividend,
                                                             final bMutablePolynomialZp divider,
                                                             final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1)
            return dividend.createFromArray(new BigInteger[]{
                    dividend.evaluate(
                            dividend.multiplyMod(BigIntegerArithmetics.mod(divider.cc().negate(), dividend.modulus), BigIntegerArithmetics.modInverse(divider.lc(), dividend.modulus)))
            });
        return null;
    }

    /**
     * Returns remainder of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static bMutablePolynomialZp remainder(final bMutablePolynomialZp dividend,
                                                 final bMutablePolynomialZp divider,
                                                 final boolean copy) {
        bMutablePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static bMutablePolynomialZp remainderClassical0(final bMutablePolynomialZp dividend,
                                                    final bMutablePolynomialZp divider,
                                                    final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        bMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        BigInteger lcInverse = BigIntegerArithmetics.modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, remainder.multiplyMod(remainder.lc(), lcInverse), i);

        return remainder;
    }

    /**
     * Fast remainder using Newton's iteration with switch to classical remainder for small polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static bMutablePolynomialZp remainderFast(final bMutablePolynomialZp dividend,
                                                     final bMutablePolynomialZp divider,
                                                     final InverseModMonomial<bMutablePolynomialZp> invMod,
                                                     final boolean copy) {
        bMutablePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[1];
    }

    /**
     * Returns quotient of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static bMutablePolynomialZp quotient(final bMutablePolynomialZp dividend,
                                                final bMutablePolynomialZp divider,
                                                final boolean copy) {
        bMutablePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic(dividend, divider, copy)[0];

        return divideAndRemainderFast0(dividend, divider, copy)[0];
    }

    /**
     * Fast quotient using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static bMutablePolynomialZp quotientFast(final bMutablePolynomialZp dividend,
                                                    final bMutablePolynomialZp divider,
                                                    final InverseModMonomial<bMutablePolynomialZp> invMod,
                                                    final boolean copy) {
        bMutablePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[0];
    }


    /* ********************************** Common conversion ********************************** */

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomialZ<T>> T[] pseudoDivideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T[]) pseudoDivideAndRemainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else
            return (T[]) pseudoDivideAndRemainder((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T[] divideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T[]) divideAndRemainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T[]) divideAndRemainder((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof bMutablePolynomialZ)
            return (T[]) divideAndRemainder((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T[]) divideAndRemainder((bMutablePolynomialZp) dividend, (bMutablePolynomialZp) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T remainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T) remainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T) remainder((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof bMutablePolynomialZ)
            return (T) remainder((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T) remainder((bMutablePolynomialZp) dividend, (bMutablePolynomialZp) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T quotient(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T) quotient((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T) quotient((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof bMutablePolynomialZ)
            return (T) quotient((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T) quotient((bMutablePolynomialZp) dividend, (bMutablePolynomialZp) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomialZp<T>> T remainderFast(final T dividend,
                                                                      final T divider,
                                                                      final InverseModMonomial<T> invMod,
                                                                      final boolean copy) {
        if (dividend instanceof lMutablePolynomialZp)
            return (T) remainderFast((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, (InverseModMonomial<lMutablePolynomialZp>) invMod, copy);
        else
            return (T) remainderFast((bMutablePolynomialZp) dividend, (bMutablePolynomialZp) divider, (InverseModMonomial<bMutablePolynomialZp>) invMod, copy);
    }
}
