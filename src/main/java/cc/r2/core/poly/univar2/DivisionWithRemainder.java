package cc.r2.core.poly.univar2;


import cc.r2.core.poly.Domain;
import cc.r2.core.poly.LongArithmetics;
import cc.r2.core.poly.lModularDomain;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.ArrayList;

import static cc.r2.core.poly.LongArithmetics.safeAdd;
import static cc.r2.core.poly.LongArithmetics.safeMultiply;
import static cc.r2.core.poly.LongArithmetics.safePow;
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
    public static lMutablePolynomialZ[] pseudoDivideAndRemainder(final lMutablePolynomialZ dividend,
                                                                 final lMutablePolynomialZ divider,
                                                                 final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), lMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new lMutablePolynomialZ[]{lMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        long factor = safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new lMutablePolynomialZ[]{(copy ? dividend.clone() : dividend).multiply(factor / divider.lc()), lMutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderClassic0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static lMutablePolynomialZ[] divideAndRemainderClassic0(final lMutablePolynomialZ dividend,
                                                            final lMutablePolynomialZ divider,
                                                            final long dividendRaiseFactor,
                                                            final boolean copy) {
        assert dividend.degree >= divider.degree;

        if (divider.lc() == 1 && dividendRaiseFactor == 1)
            return divideAndRemainderClassicMonic(dividend, divider, copy);

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
    private static lMutablePolynomialZ[] divideAndRemainderClassicMonic(final lMutablePolynomialZ dividend,
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
                        quotient[j] = safeMultiply(quotient[j], factor);
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
        return new lMutablePolynomialZ[]{lMutablePolynomialZ.create(quotient), lMutablePolynomialZ.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] divideAndRemainderLinearDivider(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZ[] pseudoDivideAndRemainderLinearDivider(lMutablePolynomialZ dividend, lMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, safePow(divider.lc(), dividend.degree), copy);
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
            res = safeAdd(safeMultiply(res, cc), safeMultiply(raiseFactor, tmp));
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


    /* ************************************ Machine-precision division in Zp[x]  ************************************ */


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(IMutablePolynomial<?> dividend,
                                                IMutablePolynomial<?> divider) {
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
            return new lMutablePolynomialZp[]{(copy ? dividend.clone() : dividend).divide(divider.lc()), dividend.createZero()};
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
        dividend.checkCompatible(divider);

        lMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        long lcInverse = dividend.domain.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.domain.multiplyMod(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new lMutablePolynomialZp[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static lMutablePolynomialZp[] divideAndRemainderLinearDividerModulus(lMutablePolynomialZp dividend, lMutablePolynomialZp divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatible(divider);

        //apply Horner's method

        long cc = dividend.domain.negateMod(divider.cc());
        long lcInverse = dividend.domain.reciprocal(divider.lc());

        if (divider.lc() != 1)
            cc = dividend.domain.multiplyMod(cc, lcInverse);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = dividend.domain.multiplyMod(res, lcInverse);
            res = dividend.domain.addMod(dividend.domain.multiplyMod(res, cc), tmp);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new lMutablePolynomialZp[]{dividend.createFromArray(quotient), dividend.createFromArray(new long[]{res})};
    }

    /* ************************************ Multi-precision division ************************************ */

    /** early checks for division */
    private static <E> gMutablePolynomial<E>[] earlyDivideAndRemainderChecks(final gMutablePolynomial<E> dividend,
                                                                             final gMutablePolynomial<E> divider,
                                                                             final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return dividend.arrayNewInstance(dividend.createZero(), dividend.createZero());
        if (dividend.degree < divider.degree)
            return dividend.arrayNewInstance(dividend.createZero(), copy ? dividend.clone() : dividend);
        if (divider.degree == 0) {
            gMutablePolynomial<E> div = copy ? dividend.clone() : dividend;
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
    public static <E> gMutablePolynomial<E>[] divideAndRemainder(final gMutablePolynomial<E> dividend,
                                                                 final gMutablePolynomial<E> divider,
                                                                 final boolean copy) {
        gMutablePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static <E> gMutablePolynomial<E>[] pseudoDivideAndRemainder(final gMutablePolynomial<E> dividend,
                                                                       final gMutablePolynomial<E> divider,
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
    public static <E> gMutablePolynomial<E>[] divideAndRemainderClassic(final gMutablePolynomial<E> dividend,
                                                                        final gMutablePolynomial<E> divider,
                                                                        final boolean copy) {
        gMutablePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, dividend.domain.getOne(), copy);
    }

    /** Plain school implementation */
    static <E> gMutablePolynomial<E>[] divideAndRemainderClassic0(final gMutablePolynomial<E> dividend,
                                                                  final gMutablePolynomial<E> divider,
                                                                  final E dividendRaiseFactor,
                                                                  final boolean copy) {
        assert dividend.degree >= divider.degree;

        Domain<E> domain = dividend.domain;
        gMutablePolynomial<E>
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
    static <E> gMutablePolynomial<E>[] divideAndRemainderLinearDivider(gMutablePolynomial<E> dividend, gMutablePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.domain.getOne(), copy);
    }


    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> gMutablePolynomial<E>[] pseudoDivideAndRemainderLinearDivider(gMutablePolynomial<E> dividend, gMutablePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.domain.pow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> gMutablePolynomial<E>[] divideAndRemainderLinearDivider0(gMutablePolynomial<E> dividend, gMutablePolynomial<E> divider, E raiseFactor, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatible(divider);

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
                quotient[i] = res;
            res = domain.add(domain.multiply(res, cc), domain.multiply(raiseFactor, tmp));
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
        long lcInv = dividend.domain.reciprocal(lc);
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
    public static <E> gMutablePolynomial<E>[] divideAndRemainderFast(gMutablePolynomial<E> dividend,
                                                                     gMutablePolynomial<E> divider,
                                                                     boolean copy) {
        gMutablePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static <E> gMutablePolynomial<E>[] divideAndRemainderFast(gMutablePolynomial<E> dividend,
                                                                     gMutablePolynomial<E> divider,
                                                                     InverseModMonomial<gMutablePolynomial<E>> invMod,
                                                                     boolean copy) {
        gMutablePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invMod, copy);
    }

    static <E> gMutablePolynomial<E>[] divideAndRemainderFast0(gMutablePolynomial<E> dividend,
                                                               gMutablePolynomial<E> divider,
                                                               boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        E lc = divider.lc();
        E lcInv = dividend.domain.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        gMutablePolynomial<E>[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    /* ********************************* Machine-precision remainders ******************************** */

    /** fast division checks */
    private static lMutablePolynomialZp earlyRemainderChecks(final lMutablePolynomialZp dividend,
                                                             final lMutablePolynomialZp divider,
                                                             final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1) {
            lModularDomain domain = dividend.domain;
            return dividend.createFromArray(new long[]{
                    dividend.evaluate(
                            domain.multiplyMod(domain.negateMod(divider.cc()), domain.reciprocal(divider.lc())))
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
        dividend.checkCompatible(divider);

        lMutablePolynomialZp remainder = copy ? dividend.clone() : dividend;
        long lcInverse = dividend.domain.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, remainder.domain.multiplyMod(remainder.lc(), lcInverse), i);

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

    /* ********************************* Multi-precision remainders ******************************** */


    /** fast division checks */
    private static <E> gMutablePolynomial<E> earlyRemainderChecks(final gMutablePolynomial<E> dividend,
                                                                  final gMutablePolynomial<E> divider,
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
    public static <E> gMutablePolynomial<E> remainder(final gMutablePolynomial<E> dividend,
                                                      final gMutablePolynomial<E> divider,
                                                      final boolean copy) {
        gMutablePolynomial<E> rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static <E> gMutablePolynomial<E> remainderClassical0(final gMutablePolynomial<E> dividend,
                                                         final gMutablePolynomial<E> divider,
                                                         final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatible(divider);

        gMutablePolynomial<E> remainder = copy ? dividend.clone() : dividend;
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
    public static <E> gMutablePolynomial<E> remainderFast(final gMutablePolynomial<E> dividend,
                                                          final gMutablePolynomial<E> divider,
                                                          final InverseModMonomial<gMutablePolynomial<E>> invMod,
                                                          final boolean copy) {
        gMutablePolynomial<E> rem = earlyRemainderChecks(dividend, divider, copy);
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
    public static <E> gMutablePolynomial<E> quotient(final gMutablePolynomial<E> dividend,
                                                     final gMutablePolynomial<E> divider,
                                                     final boolean copy) {
        gMutablePolynomial<E>[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static <E> gMutablePolynomial<E> quotientFast(final gMutablePolynomial<E> dividend,
                                                         final gMutablePolynomial<E> divider,
                                                         final InverseModMonomial<gMutablePolynomial<E>> invMod,
                                                         final boolean copy) {
        gMutablePolynomial<E>[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[0];
    }


    /* ********************************** Common conversion ********************************** */

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T[] pseudoDivideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T[]) pseudoDivideAndRemainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        if (dividend instanceof lMutablePolynomialZp)
            return (T[]) divideAndRemainder((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof gMutablePolynomial) {
            if (dividend.isOverField())
                return (T[]) divideAndRemainder((gMutablePolynomial) dividend, (gMutablePolynomial) divider, copy);
            else
                return (T[]) pseudoDivideAndRemainder((gMutablePolynomial) dividend, (gMutablePolynomial) divider, copy);
        } else
            throw new RuntimeException(dividend.getClass().toString());
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T[] divideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T[]) divideAndRemainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T[]) divideAndRemainder((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof gMutablePolynomial)
            return (T[]) divideAndRemainder((gMutablePolynomial) dividend, (gMutablePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T remainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T) remainder((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T) remainder((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof gMutablePolynomial)
            return (T) remainder((gMutablePolynomial) dividend, (gMutablePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T quotient(T dividend, T divider, boolean copy) {
        if (dividend instanceof lMutablePolynomialZ)
            return (T) quotient((lMutablePolynomialZ) dividend, (lMutablePolynomialZ) divider, copy);
        else if (dividend instanceof lMutablePolynomialZp)
            return (T) quotient((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, copy);
        else if (dividend instanceof gMutablePolynomial)
            return (T) quotient((gMutablePolynomial) dividend, (gMutablePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    @SuppressWarnings("unchecked")
    public static <T extends IMutablePolynomial<T>> T remainderFast(final T dividend,
                                                                    final T divider,
                                                                    final InverseModMonomial<T> invMod,
                                                                    final boolean copy) {
        if (dividend instanceof lMutablePolynomialZp)
            return (T) remainderFast((lMutablePolynomialZp) dividend, (lMutablePolynomialZp) divider, (InverseModMonomial<lMutablePolynomialZp>) invMod, copy);
        else if (dividend instanceof gMutablePolynomial)
            return (T) remainderFast((gMutablePolynomial) dividend, (gMutablePolynomial) divider, (InverseModMonomial) invMod, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }
}
