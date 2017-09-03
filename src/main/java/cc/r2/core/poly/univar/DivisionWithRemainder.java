package cc.r2.core.poly.univar;


import cc.r2.core.poly.Domain;
import cc.r2.core.poly.LongArithmetics;
import cc.r2.core.poly.lIntegersModulo;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.ArrayList;

import static cc.r2.core.poly.LongArithmetics.*;
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
    public static <T extends IUnivariatePolynomial<T>> T remainderMonomial(T dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).truncate(xDegree - 1);
    }

    private static void checkZeroDivider(IUnivariatePolynomial p) {
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
    public static lUnivariatePolynomialZ[] divideAndRemainder(final lUnivariatePolynomialZ dividend,
                                                              final lUnivariatePolynomialZ divider,
                                                              boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), lUnivariatePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            lUnivariatePolynomialZ div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new lUnivariatePolynomialZ[]{div, lUnivariatePolynomialZ.zero()};
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
        return divideAndRemainderClassic0(dividend, divider, 1, copy);
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
    public static lUnivariatePolynomialZ[] pseudoDivideAndRemainder(final lUnivariatePolynomialZ dividend,
                                                                    final lUnivariatePolynomialZ divider,
                                                                    final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), lUnivariatePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        long factor = safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new lUnivariatePolynomialZ[]{(copy ? dividend.clone() : dividend).multiply(factor / divider.lc()), lUnivariatePolynomialZ.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderClassic0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static lUnivariatePolynomialZ[] divideAndRemainderClassic0(final lUnivariatePolynomialZ dividend,
                                                               final lUnivariatePolynomialZ divider,
                                                               final long dividendRaiseFactor,
                                                               final boolean copy) {
        assert dividend.degree >= divider.degree;

        if (divider.lc() == 1 && dividendRaiseFactor == 1)
            return divideAndRemainderClassicMonic(dividend, divider, copy);

        lUnivariatePolynomialZ
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

        return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.create(quotient), remainder};
    }

    /** Plain school implementation */
    private static lUnivariatePolynomialZ[] divideAndRemainderClassicMonic(final lUnivariatePolynomialZ dividend,
                                                                           final lUnivariatePolynomialZ divider,
                                                                           final boolean copy) {
        assert divider.lc() == 1;

        lUnivariatePolynomialZ
                remainder = (copy ? dividend.clone() : dividend);
        long[] quotient = new long[dividend.degree - divider.degree + 1];
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.lc();
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }
        return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.create(quotient), remainder};
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
    static lUnivariatePolynomialZ[] pseudoDivideAndRemainderAdaptive(final lUnivariatePolynomialZ dividend,
                                                                     final lUnivariatePolynomialZ divider,
                                                                     final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), lUnivariatePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new lUnivariatePolynomialZ[]{copy ? dividend.clone() : dividend, lUnivariatePolynomialZ.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static lUnivariatePolynomialZ[] pseudoDivideAndRemainderAdaptive0(final lUnivariatePolynomialZ dividend,
                                                                      final lUnivariatePolynomialZ divider,
                                                                      final boolean copy) {
        assert dividend.degree >= divider.degree;

        lUnivariatePolynomialZ remainder = copy ? dividend.clone() : dividend;
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
                        quotient[j] = safeMultiply(quotient[j], factor);
                    quot = divideSignedFast(remainder.lc(), magic);
                }

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lUnivariatePolynomialZ[] pseudoDivideAndRemainderLinearDividerAdaptive(lUnivariatePolynomialZ dividend, lUnivariatePolynomialZ divider, boolean copy) {
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
            res = safeAdd(safeMultiply(res, cc), safeMultiply(factor, tmp));
            if (i == 0) break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res) {
                long gcd = LongArithmetics.gcd(res, lc), f = lc / gcd;
                factor = safeMultiply(factor, f);
                res = safeMultiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = safeMultiply(quotient[j], f);
                quot = divideSignedFast(res, magic);
            }
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.create(quotient), lUnivariatePolynomialZ.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lUnivariatePolynomialZ[] divideAndRemainderLinearDivider(lUnivariatePolynomialZ dividend, lUnivariatePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lUnivariatePolynomialZ[] pseudoDivideAndRemainderLinearDivider(lUnivariatePolynomialZ dividend, lUnivariatePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lUnivariatePolynomialZ[] divideAndRemainderLinearDivider0(lUnivariatePolynomialZ dividend, lUnivariatePolynomialZ divider, long raiseFactor, boolean copy) {
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
            res = safeAdd(safeMultiply(res, cc), safeMultiply(raiseFactor, tmp));
            if (i == 0)
                break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res)
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lUnivariatePolynomialZ[]{lUnivariatePolynomialZ.create(quotient), lUnivariatePolynomialZ.create(res)};
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
    public static lUnivariatePolynomialZ remainder(final lUnivariatePolynomialZ dividend,
                                                   final lUnivariatePolynomialZ divider,
                                                   final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return lUnivariatePolynomialZ.zero();
        if (divider.degree == 1)
            if (divider.cc() % divider.lc() == 0)
                return lUnivariatePolynomialZ.create(dividend.evaluate(-divider.cc() / divider.lc()));
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static lUnivariatePolynomialZ remainder0(final lUnivariatePolynomialZ dividend,
                                             final lUnivariatePolynomialZ divider,
                                             final boolean copy) {
        assert dividend.degree >= divider.degree;

        lUnivariatePolynomialZ remainder = copy ? dividend.clone() : dividend;
        Magic magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc())
                    return null;
                remainder.subtract(divider, quot, i);
            }
        return remainder;
    }

    /**
     * Returns quotient {@code dividend/ divider} or null if exact division o
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static lUnivariatePolynomialZ quotient(final lUnivariatePolynomialZ dividend,
                                                  final lUnivariatePolynomialZ divider,
                                                  final boolean copy) {
        lUnivariatePolynomialZ[] qd = divideAndRemainder(dividend, divider, copy);
        if (qd == null)
            return null;

        return qd[0];
    }

    /* ************************************ Machine-precision division in Zp[x]  ************************************ */


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(IUnivariatePolynomial dividend,
                                                IUnivariatePolynomial divider) {
        // practical benchmarks show that without pre-conditioning,
        // classical division is always faster or at least the same fast
        return true;
    }

    /** early checks for division */
    private static lUnivariatePolynomialZp[] earlyDivideAndRemainderChecks(final lUnivariatePolynomialZp dividend,
                                                                           final lUnivariatePolynomialZp divider,
                                                                           final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lUnivariatePolynomialZp[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new lUnivariatePolynomialZp[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new lUnivariatePolynomialZp[]{(copy ? dividend.clone() : dividend).divide(divider.lc()), dividend.createZero()};
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
    public static lUnivariatePolynomialZp[] divideAndRemainder(final lUnivariatePolynomialZp dividend,
                                                               final lUnivariatePolynomialZp divider,
                                                               final boolean copy) {
        lUnivariatePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static lUnivariatePolynomialZp[] divideAndRemainderClassic(final lUnivariatePolynomialZp dividend,
                                                                      final lUnivariatePolynomialZp divider,
                                                                      final boolean copy) {
        lUnivariatePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static lUnivariatePolynomialZp[] divideAndRemainderClassic0(final lUnivariatePolynomialZp dividend,
                                                                final lUnivariatePolynomialZp divider,
                                                                final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkSameDomainWith(divider);

        lUnivariatePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        long lcInverse = dividend.domain.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.domain.multiply(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new lUnivariatePolynomialZp[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lUnivariatePolynomialZp[] divideAndRemainderLinearDividerModulus(lUnivariatePolynomialZp dividend, lUnivariatePolynomialZp divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkSameDomainWith(divider);

        //apply Horner's method

        long cc = dividend.domain.negate(divider.cc());
        long lcInverse = dividend.domain.reciprocal(divider.lc());

        if (divider.lc() != 1)
            cc = dividend.domain.multiply(cc, lcInverse);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = dividend.domain.multiply(res, lcInverse);
            res = dividend.domain.add(dividend.domain.multiply(res, cc), tmp);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lUnivariatePolynomialZp[]{dividend.createFromArray(quotient), dividend.createFromArray(new long[]{res})};
    }

    /* ************************************ Multi-precision division ************************************ */

    /** early checks for division */
    private static <E> UnivariatePolynomial<E>[] earlyDivideAndRemainderChecks(final UnivariatePolynomial<E> dividend,
                                                                               final UnivariatePolynomial<E> divider,
                                                                               final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return dividend.arrayNewInstance(dividend.createZero(), dividend.createZero());
        if (dividend.degree < divider.degree)
            return dividend.arrayNewInstance(dividend.createZero(), copy ? dividend.clone() : dividend);
        if (divider.degree == 0) {
            UnivariatePolynomial<E> div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return dividend.arrayNewInstance(div, dividend.createZero());
        } if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
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
    public static <E> UnivariatePolynomial<E>[] divideAndRemainder(final UnivariatePolynomial<E> dividend,
                                                                   final UnivariatePolynomial<E> divider,
                                                                   final boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic0(dividend, divider, dividend.domain.getOne(), copy);

        return divideAndRemainderFast0(dividend, divider, copy);
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
    public static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainder(final UnivariatePolynomial<E> dividend,
                                                                         final UnivariatePolynomial<E> divider,
                                                                         final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return dividend.arrayNewInstance(dividend.createZero(), dividend.createZero());
        if (dividend.degree < divider.degree)
            return dividend.arrayNewInstance(dividend.createZero(), copy ? dividend.clone() : dividend);
        E factor = dividend.domain.pow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return dividend.arrayNewInstance((copy ? dividend.clone() : dividend).multiply(dividend.domain.divideExact(factor, divider.lc())), dividend.createZero());
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderClassic0(dividend, divider, factor, copy);
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
    public static <E> UnivariatePolynomial<E>[] divideAndRemainderClassic(final UnivariatePolynomial<E> dividend,
                                                                          final UnivariatePolynomial<E> divider,
                                                                          final boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, dividend.domain.getOne(), copy);
    }

    /** Plain school implementation */
    static <E> UnivariatePolynomial<E>[] divideAndRemainderClassic0(final UnivariatePolynomial<E> dividend,
                                                                    final UnivariatePolynomial<E> divider,
                                                                    final E dividendRaiseFactor,
                                                                    final boolean copy) {
        assert dividend.degree >= divider.degree;

        Domain<E> domain = dividend.domain;
        UnivariatePolynomial<E>
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        E[] quotient = domain.createArray(dividend.degree - divider.degree + 1);

        E lcMultiplier, lcDivider;
        if (domain.isField()) {
            lcMultiplier = domain.reciprocal(divider.lc());
            lcDivider = domain.getOne();
        } else {
            lcMultiplier = domain.getOne();
            lcDivider = divider.lc();
        }

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                E quot = domain.divideOrNull(domain.multiply(remainder.lc(), lcMultiplier), lcDivider);
                if (quot == null)
                    return null;

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = domain.getZero();
        }

        return dividend.arrayNewInstance(dividend.createFromArray(quotient), remainder);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] divideAndRemainderLinearDivider(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.domain.getOne(), copy);
    }


    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainderLinearDivider(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.domain.pow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] divideAndRemainderLinearDivider0(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, E raiseFactor, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkSameDomainWith(divider);

        //apply Horner's method

        Domain<E> domain = dividend.domain;
        E cc = domain.negate(divider.cc()), lcDivider, lcMultiplier;
        if (domain.isField()) {
            lcMultiplier = domain.reciprocal(divider.lc());
            lcDivider = domain.getOne();
        } else {
            lcMultiplier = domain.getOne();
            lcDivider = divider.lc();
        }

        E[] quotient = copy ? domain.createArray(dividend.degree) : dividend.data;
        E res = domain.getZero();
        for (int i = dividend.degree; ; --i) {
            E tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = domain.copy(res);
            res = domain.addMutable(domain.multiplyMutable(res, cc), domain.multiply(raiseFactor, tmp));
            if (i == 0)
                break;
            E quot = domain.divideOrNull(domain.multiply(res, lcMultiplier), lcDivider);
            if (quot == null)
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = domain.getZero();
        return dividend.arrayNewInstance(dividend.createFromArray(quotient), dividend.createConstant(res));
    }


    /* ********************************** Fast division algorithm ********************************** */

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomial<Poly extends IUnivariatePolynomial<Poly>>
            implements java.io.Serializable {
        final Poly poly;

        private InverseModMonomial(Poly poly) {
            if (!poly.isUnitCC())
                throw new IllegalArgumentException("Smallest coefficient is not a unit: " + poly);
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<Poly> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public Poly getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() >= r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            Poly gPrev = currentSize == 0 ? poly.createOne() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                Poly tmp = gPrev.clone().multiply(2).subtract(gPrev.clone().square().multiply(poly));
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
    public static <Poly extends IUnivariatePolynomial<Poly>> InverseModMonomial<Poly> fastDivisionPreConditioning(Poly divider) {
        if (!divider.isMonic())
            throw new IllegalArgumentException("Only monic polynomials allowed. Input: " + divider);
        return new InverseModMonomial<>(divider.clone().reverse());
    }

    /**
     * Prepares {@code rev(divider)^(-1) mod x^i } for fast division.
     *
     * @param divider the divider
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> InverseModMonomial<Poly> fastDivisionPreConditioningWithLCCorrection(Poly divider) {
        return new InverseModMonomial<>(divider.clone().monic().reverse());
    }

    /** fast division implementation */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] divideAndRemainderFast0(Poly dividend, Poly divider,
                                                                                            InverseModMonomial<Poly> invRevMod,
                                                                                            boolean copy) {
        int m = dividend.degree() - divider.degree();
        Poly q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1)), m + 1, false).reverse();
        if (q.degree() < m)
            q.shiftRight(m - q.degree());
        Poly r = (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(q));
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
    public static lUnivariatePolynomialZp[] divideAndRemainderFast(lUnivariatePolynomialZp dividend,
                                                                   lUnivariatePolynomialZp divider,
                                                                   boolean copy) {
        lUnivariatePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static lUnivariatePolynomialZp[] divideAndRemainderFast(lUnivariatePolynomialZp dividend,
                                                                   lUnivariatePolynomialZp divider,
                                                                   InverseModMonomial<lUnivariatePolynomialZp> invMod,
                                                                   boolean copy) {
        lUnivariatePolynomialZp[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy);
    }

    static lUnivariatePolynomialZp[] divideAndRemainderFastCorrectLC(lUnivariatePolynomialZp dividend,
                                                                     lUnivariatePolynomialZp divider,
                                                                     InverseModMonomial<lUnivariatePolynomialZp> invMod,
                                                                     boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, invMod, copy);

        long lc = divider.lc();
        long lcInv = divider.domain.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        lUnivariatePolynomialZp[] result = divideAndRemainderFast0(dividend, divider, invMod, copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    static lUnivariatePolynomialZp[] divideAndRemainderFast0(lUnivariatePolynomialZp dividend,
                                                             lUnivariatePolynomialZp divider,
                                                             boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        long lc = divider.lc();
        long lcInv = divider.domain.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        lUnivariatePolynomialZp[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
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
    public static <E> UnivariatePolynomial<E>[] divideAndRemainderFast(UnivariatePolynomial<E> dividend,
                                                                       UnivariatePolynomial<E> divider,
                                                                       boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static <E> UnivariatePolynomial<E>[] divideAndRemainderFast(UnivariatePolynomial<E> dividend,
                                                                       UnivariatePolynomial<E> divider,
                                                                       InverseModMonomial<UnivariatePolynomial<E>> invMod,
                                                                       boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy);
    }

    static <E> UnivariatePolynomial<E>[] divideAndRemainderFastCorrectLC(UnivariatePolynomial<E> dividend,
                                                                         UnivariatePolynomial<E> divider,
                                                                         InverseModMonomial<UnivariatePolynomial<E>> invMod,
                                                                         boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, invMod, copy);

        E lc = divider.lc();
        E lcInv = dividend.domain.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        UnivariatePolynomial<E>[] result = divideAndRemainderFast0(dividend, divider, invMod, copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    static <E> UnivariatePolynomial<E>[] divideAndRemainderFast0(UnivariatePolynomial<E> dividend,
                                                                 UnivariatePolynomial<E> divider,
                                                                 boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        E lc = divider.lc();
        E lcInv = dividend.domain.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        UnivariatePolynomial<E>[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    /* ********************************* Machine-precision remainders ******************************** */

    /** fast division checks */
    private static lUnivariatePolynomialZp earlyRemainderChecks(final lUnivariatePolynomialZp dividend,
                                                                final lUnivariatePolynomialZp divider,
                                                                final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1) {
            lIntegersModulo domain = dividend.domain;
            return dividend.createFromArray(new long[]{
                    dividend.evaluate(
                            domain.multiply(domain.negate(divider.cc()), domain.reciprocal(divider.lc())))
            });
        }
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
    public static lUnivariatePolynomialZp remainder(final lUnivariatePolynomialZp dividend,
                                                    final lUnivariatePolynomialZp divider,
                                                    final boolean copy) {
        lUnivariatePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static lUnivariatePolynomialZp remainderClassical0(final lUnivariatePolynomialZp dividend,
                                                       final lUnivariatePolynomialZp divider,
                                                       final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkSameDomainWith(divider);

        lUnivariatePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long lcInverse = dividend.domain.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, remainder.domain.multiply(remainder.lc(), lcInverse), i);

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
    public static lUnivariatePolynomialZp remainderFast(final lUnivariatePolynomialZp dividend,
                                                        final lUnivariatePolynomialZp divider,
                                                        final InverseModMonomial<lUnivariatePolynomialZp> invMod,
                                                        final boolean copy) {
        lUnivariatePolynomialZp rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy)[1];
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
    public static lUnivariatePolynomialZp quotient(final lUnivariatePolynomialZp dividend,
                                                   final lUnivariatePolynomialZp divider,
                                                   final boolean copy) {
        lUnivariatePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static lUnivariatePolynomialZp quotientFast(final lUnivariatePolynomialZp dividend,
                                                       final lUnivariatePolynomialZp divider,
                                                       final InverseModMonomial<lUnivariatePolynomialZp> invMod,
                                                       final boolean copy) {
        lUnivariatePolynomialZp[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy)[0];
    }

    /* ********************************* Multi-precision remainders ******************************** */


    /** fast division checks */
    private static <E> UnivariatePolynomial<E> earlyRemainderChecks(final UnivariatePolynomial<E> dividend,
                                                                    final UnivariatePolynomial<E> divider,
                                                                    final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1) {
            E p = dividend.domain.divideOrNull(dividend.domain.negate(divider.cc()), divider.lc());
            if (p == null)
                return null;
            return dividend.createConstant(dividend.evaluate(p));
        }
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
    public static <E> UnivariatePolynomial<E> remainder(final UnivariatePolynomial<E> dividend,
                                                        final UnivariatePolynomial<E> divider,
                                                        final boolean copy) {
        UnivariatePolynomial<E> rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static <E> UnivariatePolynomial<E> remainderClassical0(final UnivariatePolynomial<E> dividend,
                                                           final UnivariatePolynomial<E> divider,
                                                           final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkSameDomainWith(divider);

        UnivariatePolynomial<E> remainder = copy ? dividend.clone() : dividend;
        Domain<E> domain = dividend.domain;
        if (domain.isField()) {
            E lcInverse = domain.reciprocal(divider.lc());
            for (int i = dividend.degree - divider.degree; i >= 0; --i)
                if (remainder.degree == divider.degree + i)
                    remainder.subtract(divider, domain.multiply(remainder.lc(), lcInverse), i);
        } else {
            for (int i = dividend.degree - divider.degree; i >= 0; --i)
                if (remainder.degree == divider.degree + i) {
                    E quot = domain.divideOrNull(remainder.lc(), divider.lc());
                    if (quot == null)
                        return null;
                    remainder.subtract(divider, quot, i);
                }
        }

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
    public static <E> UnivariatePolynomial<E> remainderFast(final UnivariatePolynomial<E> dividend,
                                                            final UnivariatePolynomial<E> divider,
                                                            final InverseModMonomial<UnivariatePolynomial<E>> invMod,
                                                            final boolean copy) {
        UnivariatePolynomial<E> rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy)[1];
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
    public static <E> UnivariatePolynomial<E> quotient(final UnivariatePolynomial<E> dividend,
                                                       final UnivariatePolynomial<E> divider,
                                                       final boolean copy) {
        UnivariatePolynomial<E>[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static <E> UnivariatePolynomial<E> quotientFast(final UnivariatePolynomial<E> dividend,
                                                           final UnivariatePolynomial<E> divider,
                                                           final InverseModMonomial<UnivariatePolynomial<E>> invMod,
                                                           final boolean copy) {
        UnivariatePolynomial<E>[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy)[0];
    }


    /* ********************************** Common conversion ********************************** */

    /**
     * Returns quotient and remainder of {@code dividend / divider} using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] pseudoDivideAndRemainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof lUnivariatePolynomialZ)
            return (Poly[]) pseudoDivideAndRemainder((lUnivariatePolynomialZ) dividend, (lUnivariatePolynomialZ) divider, copy);
        if (dividend instanceof lUnivariatePolynomialZp)
            return (Poly[]) divideAndRemainder((lUnivariatePolynomialZp) dividend, (lUnivariatePolynomialZp) divider, copy);
        else if (dividend instanceof UnivariatePolynomial) {
            if (dividend.isOverField())
                return (Poly[]) divideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
            else
                return (Poly[]) pseudoDivideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        } else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Returns {@code {quotient, remainder}} of {@code dividend / divider} or {@code null} if the division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] divideAndRemainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof lUnivariatePolynomialZ)
            return (Poly[]) divideAndRemainder((lUnivariatePolynomialZ) dividend, (lUnivariatePolynomialZ) divider, copy);
        else if (dividend instanceof lUnivariatePolynomialZp)
            return (Poly[]) divideAndRemainder((lUnivariatePolynomialZp) dividend, (lUnivariatePolynomialZp) divider, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly[]) divideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Divides {@code dividend} by {@code divider} or throws {@code ArithmeticException} if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider}
     * @throws ArithmeticException if exact division is not possible
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly divideExact(Poly dividend, Poly divider, boolean copy) {
        Poly[] qr = divideAndRemainder(dividend, divider, copy);
        if (!qr[1].isZero())
            throw new ArithmeticException("Not divisible: (" + dividend + ") / (" + divider + ")");
        return qr[0];
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
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly remainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof lUnivariatePolynomialZ)
            return (Poly) remainder((lUnivariatePolynomialZ) dividend, (lUnivariatePolynomialZ) divider, copy);
        else if (dividend instanceof lUnivariatePolynomialZp)
            return (Poly) remainder((lUnivariatePolynomialZp) dividend, (lUnivariatePolynomialZp) divider, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly) remainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Returns quotient {@code dividend/ divider} or null if exact division o
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly quotient(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof lUnivariatePolynomialZ)
            return (Poly) quotient((lUnivariatePolynomialZ) dividend, (lUnivariatePolynomialZ) divider, copy);
        else if (dividend instanceof lUnivariatePolynomialZp)
            return (Poly) quotient((lUnivariatePolynomialZp) dividend, (lUnivariatePolynomialZp) divider, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly) quotient((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
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
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly remainderFast(final Poly dividend,
                                                                                final Poly divider,
                                                                                final InverseModMonomial<Poly> invMod,
                                                                                final boolean copy) {
        if (dividend instanceof lUnivariatePolynomialZp)
            return (Poly) remainderFast((lUnivariatePolynomialZp) dividend, (lUnivariatePolynomialZp) divider, (InverseModMonomial<lUnivariatePolynomialZp>) invMod, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly) remainderFast((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, (InverseModMonomial) invMod, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }
}
