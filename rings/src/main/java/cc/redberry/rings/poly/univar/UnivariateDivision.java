package cc.redberry.rings.poly.univar;


import cc.redberry.libdivide4j.FastDivision.Magic;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.MachineArithmetic;

import java.util.ArrayList;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Division with remainder of univariate polynomials.
 *
 * @since 1.0
 */
public final class UnivariateDivision {
    private UnivariateDivision() {}

    /* **************************************** Common methods  *************************************** */

    /**
     * Returns the remainder of {@code dividend} and monomial {@code x^xDegree}
     *
     * @param dividend the dividend
     * @param xDegree  monomial degree
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    public static UnivariatePolynomialZ64[] divideAndRemainder(final UnivariatePolynomialZ64 dividend,
                                                               final UnivariatePolynomialZ64 divider,
                                                               boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), UnivariatePolynomialZ64.zero()};
        if (dividend.degree < divider.degree)
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            UnivariatePolynomialZ64 div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new UnivariatePolynomialZ64[]{div, UnivariatePolynomialZ64.zero()};
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static UnivariatePolynomialZ64[] pseudoDivideAndRemainder(final UnivariatePolynomialZ64 dividend,
                                                                     final UnivariatePolynomialZ64 divider,
                                                                     final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), UnivariatePolynomialZ64.zero()};
        if (dividend.degree < divider.degree)
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), copy ? dividend.clone() : dividend};
        long factor = MachineArithmetic.safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new UnivariatePolynomialZ64[]{(copy ? dividend.clone() : dividend).multiply(factor / divider.lc()), UnivariatePolynomialZ64.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderClassic0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static UnivariatePolynomialZ64[] divideAndRemainderClassic0(final UnivariatePolynomialZ64 dividend,
                                                                final UnivariatePolynomialZ64 divider,
                                                                final long dividendRaiseFactor,
                                                                final boolean copy) {
        assert dividend.degree >= divider.degree;

        if (divider.lc() == 1 && dividendRaiseFactor == 1)
            return divideAndRemainderClassicMonic(dividend, divider, copy);

        UnivariatePolynomialZ64
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

        return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.create(quotient), remainder};
    }

    /** Plain school implementation */
    private static UnivariatePolynomialZ64[] divideAndRemainderClassicMonic(final UnivariatePolynomialZ64 dividend,
                                                                            final UnivariatePolynomialZ64 divider,
                                                                            final boolean copy) {
        assert divider.lc() == 1;

        UnivariatePolynomialZ64
                remainder = (copy ? dividend.clone() : dividend);
        long[] quotient = new long[dividend.degree - divider.degree + 1];
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.lc();
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }
        return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.create(quotient), remainder};
    }

    /**
     * Returns quotient and remainder using adaptive pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    static UnivariatePolynomialZ64[] pseudoDivideAndRemainderAdaptive(final UnivariatePolynomialZ64 dividend,
                                                                      final UnivariatePolynomialZ64 divider,
                                                                      final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), UnivariatePolynomialZ64.zero()};
        if (dividend.degree < divider.degree)
            return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new UnivariatePolynomialZ64[]{copy ? dividend.clone() : dividend, UnivariatePolynomialZ64.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static UnivariatePolynomialZ64[] pseudoDivideAndRemainderAdaptive0(final UnivariatePolynomialZ64 dividend,
                                                                       final UnivariatePolynomialZ64 divider,
                                                                       final boolean copy) {
        assert dividend.degree >= divider.degree;

        UnivariatePolynomialZ64 remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        Magic magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc()) {
                    long gcd = MachineArithmetic.gcd(remainder.lc(), divider.lc());
                    long factor = divider.lc() / gcd;
                    remainder.multiply(factor);
                    for (int j = i + 1; j < quotient.length; ++j)
                        quotient[j] = MachineArithmetic.safeMultiply(quotient[j], factor);
                    quot = divideSignedFast(remainder.lc(), magic);
                }

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static UnivariatePolynomialZ64[] pseudoDivideAndRemainderLinearDividerAdaptive(UnivariatePolynomialZ64 dividend, UnivariatePolynomialZ64 divider, boolean copy) {
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
            res = MachineArithmetic.safeAdd(MachineArithmetic.safeMultiply(res, cc), MachineArithmetic.safeMultiply(factor, tmp));
            if (i == 0) break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res) {
                long gcd = MachineArithmetic.gcd(res, lc), f = lc / gcd;
                factor = MachineArithmetic.safeMultiply(factor, f);
                res = MachineArithmetic.safeMultiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = MachineArithmetic.safeMultiply(quotient[j], f);
                quot = divideSignedFast(res, magic);
            }
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.create(quotient), UnivariatePolynomialZ64.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static UnivariatePolynomialZ64[] divideAndRemainderLinearDivider(UnivariatePolynomialZ64 dividend, UnivariatePolynomialZ64 divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static UnivariatePolynomialZ64[] pseudoDivideAndRemainderLinearDivider(UnivariatePolynomialZ64 dividend, UnivariatePolynomialZ64 divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, MachineArithmetic.safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static UnivariatePolynomialZ64[] divideAndRemainderLinearDivider0(UnivariatePolynomialZ64 dividend, UnivariatePolynomialZ64 divider, long raiseFactor, boolean copy) {
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
            res = MachineArithmetic.safeAdd(MachineArithmetic.safeMultiply(res, cc), MachineArithmetic.safeMultiply(raiseFactor, tmp));
            if (i == 0)
                break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res)
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new UnivariatePolynomialZ64[]{UnivariatePolynomialZ64.create(quotient), UnivariatePolynomialZ64.create(res)};
    }

    /**
     * Returns remainder of {@code dividend} and {@code divider} or {@code null} if division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the remainder or {@code null} if division is not possible
     */
    public static UnivariatePolynomialZ64 remainder(final UnivariatePolynomialZ64 dividend,
                                                    final UnivariatePolynomialZ64 divider,
                                                    final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return UnivariatePolynomialZ64.zero();
        if (divider.degree == 1)
            if (divider.cc() % divider.lc() == 0)
                return UnivariatePolynomialZ64.create(dividend.evaluate(-divider.cc() / divider.lc()));
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static UnivariatePolynomialZ64 remainder0(final UnivariatePolynomialZ64 dividend,
                                              final UnivariatePolynomialZ64 divider,
                                              final boolean copy) {
        assert dividend.degree >= divider.degree;

        UnivariatePolynomialZ64 remainder = copy ? dividend.clone() : dividend;
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
     * Returns quotient {@code dividend/ divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static UnivariatePolynomialZ64 quotient(final UnivariatePolynomialZ64 dividend,
                                                   final UnivariatePolynomialZ64 divider,
                                                   final boolean copy) {
        UnivariatePolynomialZ64[] qd = divideAndRemainder(dividend, divider, copy);
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
    private static UnivariatePolynomialZp64[] earlyDivideAndRemainderChecks(final UnivariatePolynomialZp64 dividend,
                                                                            final UnivariatePolynomialZp64 divider,
                                                                            final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new UnivariatePolynomialZp64[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new UnivariatePolynomialZp64[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new UnivariatePolynomialZp64[]{(copy ? dividend.clone() : dividend).divide(divider.lc()), dividend.createZero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, copy);
        return null;
    }

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static UnivariatePolynomialZp64[] divideAndRemainder(final UnivariatePolynomialZp64 dividend,
                                                                final UnivariatePolynomialZp64 divider,
                                                                final boolean copy) {
        UnivariatePolynomialZp64[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static UnivariatePolynomialZp64[] divideAndRemainderClassic(final UnivariatePolynomialZp64 dividend,
                                                                       final UnivariatePolynomialZp64 divider,
                                                                       final boolean copy) {
        UnivariatePolynomialZp64[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static UnivariatePolynomialZp64[] divideAndRemainderClassic0(final UnivariatePolynomialZp64 dividend,
                                                                 final UnivariatePolynomialZp64 divider,
                                                                 final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.assertSameCoefficientRingWith(divider);

        UnivariatePolynomialZp64 remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        long lcInverse = dividend.ring.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.ring.multiply(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new UnivariatePolynomialZp64[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static UnivariatePolynomialZp64[] divideAndRemainderLinearDividerModulus(UnivariatePolynomialZp64 dividend, UnivariatePolynomialZp64 divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.assertSameCoefficientRingWith(divider);

        //apply Horner's method

        long cc = dividend.ring.negate(divider.cc());
        long lcInverse = dividend.ring.reciprocal(divider.lc());

        if (divider.lc() != 1)
            cc = dividend.ring.multiply(cc, lcInverse);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = dividend.ring.multiply(res, lcInverse);
            res = dividend.ring.add(dividend.ring.multiply(res, cc), tmp);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new UnivariatePolynomialZp64[]{dividend.createFromArray(quotient), dividend.createFromArray(new long[]{res})};
    }

    /* ************************************ Multi-precision division ************************************ */

    /** early checks for division */
    private static <E> UnivariatePolynomial<E>[] earlyDivideAndRemainderChecks(final UnivariatePolynomial<E> dividend,
                                                                               final UnivariatePolynomial<E> divider,
                                                                               final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return dividend.createArray(dividend.createZero(), dividend.createZero());
        if (dividend.degree < divider.degree)
            return dividend.createArray(dividend.createZero(), copy ? dividend.clone() : dividend);
        if (divider.degree == 0) {
            UnivariatePolynomial<E> div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return dividend.createArray(div, dividend.createZero());
        }
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider(dividend, divider, copy);
        return null;
    }

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static <E> UnivariatePolynomial<E>[] divideAndRemainder(final UnivariatePolynomial<E> dividend,
                                                                   final UnivariatePolynomial<E> divider,
                                                                   final boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic0(dividend, divider, dividend.ring.getOne(), copy);

        return divideAndRemainderFast0(dividend, divider, copy);
    }

    /**
     * Returns quotient and remainder using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainder(final UnivariatePolynomial<E> dividend,
                                                                         final UnivariatePolynomial<E> divider,
                                                                         final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return dividend.createArray(dividend.createZero(), dividend.createZero());
        if (dividend.degree < divider.degree)
            return dividend.createArray(dividend.createZero(), copy ? dividend.clone() : dividend);
        E factor = dividend.ring.pow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return dividend.createArray((copy ? dividend.clone() : dividend).multiply(dividend.ring.divideExact(factor, divider.lc())), dividend.createZero());
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderClassic0(dividend, divider, factor, copy);
    }

    /**
     * Classical algorithm for division with remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static <E> UnivariatePolynomial<E>[] divideAndRemainderClassic(final UnivariatePolynomial<E> dividend,
                                                                          final UnivariatePolynomial<E> divider,
                                                                          final boolean copy) {
        UnivariatePolynomial<E>[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, dividend.ring.getOne(), copy);
    }

    /** Plain school implementation */
    static <E> UnivariatePolynomial<E>[] divideAndRemainderClassic0(final UnivariatePolynomial<E> dividend,
                                                                    final UnivariatePolynomial<E> divider,
                                                                    final E dividendRaiseFactor,
                                                                    final boolean copy) {
        assert dividend.degree >= divider.degree;

        Ring<E> ring = dividend.ring;
        UnivariatePolynomial<E>
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        E[] quotient = ring.createArray(dividend.degree - divider.degree + 1);

        E lcMultiplier, lcDivider;
        if (ring.isField()) {
            lcMultiplier = ring.reciprocal(divider.lc());
            lcDivider = ring.getOne();
        } else {
            lcMultiplier = ring.getOne();
            lcDivider = divider.lc();
        }

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                E quot = ring.divideOrNull(ring.multiply(remainder.lc(), lcMultiplier), lcDivider);
                if (quot == null)
                    return null;

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = ring.getZero();
        }

        return dividend.createArray(dividend.createFromArray(quotient), remainder);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] divideAndRemainderLinearDivider(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.ring.getOne(), copy);
    }


    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainderLinearDivider(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, dividend.ring.pow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static <E> UnivariatePolynomial<E>[] divideAndRemainderLinearDivider0(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, E raiseFactor, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.assertSameCoefficientRingWith(divider);

        //apply Horner's method

        Ring<E> ring = dividend.ring;
        E cc = ring.negate(divider.cc()), lcDivider, lcMultiplier;
        if (ring.isField()) {
            lcMultiplier = ring.reciprocal(divider.lc());
            lcDivider = ring.getOne();
        } else {
            lcMultiplier = ring.getOne();
            lcDivider = divider.lc();
        }

        E[] quotient = copy ? ring.createArray(dividend.degree) : dividend.data;
        E res = ring.getZero();
        for (int i = dividend.degree; ; --i) {
            E tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = ring.copy(res);
            res = ring.addMutable(ring.multiplyMutable(res, cc), ring.multiply(raiseFactor, tmp));
            if (i == 0)
                break;
            E quot = ring.divideOrNull(ring.multiply(res, lcMultiplier), lcDivider);
            if (quot == null)
                return null;
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = ring.getZero();
        return dividend.createArray(dividend.createFromArray(quotient), dividend.createConstant(res));
    }


    /**
     * Returns quotient and remainder using adaptive pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainderAdaptive(final UnivariatePolynomial<E> dividend,
                                                                          final UnivariatePolynomial<E> divider,
                                                                          final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return new UnivariatePolynomial[]{UnivariatePolynomial.zero(dividend.ring), UnivariatePolynomial.zero(dividend.ring)};
        if (dividend.degree < divider.degree)
            return new UnivariatePolynomial[]{UnivariatePolynomial.zero(dividend.ring), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new UnivariatePolynomial[]{copy ? dividend.clone() : dividend, UnivariatePolynomial.zero(dividend.ring)};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainderAdaptive0(final UnivariatePolynomial<E> dividend,
                                                                           final UnivariatePolynomial<E> divider,
                                                                           final boolean copy) {
        assert dividend.degree >= divider.degree;

        Ring<E> ring = dividend.ring;
        UnivariatePolynomial<E> remainder = copy ? dividend.clone() : dividend;
        E[] quotient = ring.createArray(dividend.degree - divider.degree + 1);

        E dlc = divider.lc();
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                E quot = ring.divideOrNull(remainder.lc(), dlc);
                if (quot == null) {
                    E gcd = ring.gcd(remainder.lc(), divider.lc());
                    E factor = ring.divideExact(divider.lc(), gcd);
                    remainder.multiply(factor);
                    for (int j = i + 1; j < quotient.length; ++j)
                        quotient[j] = ring.multiply(quotient[j], factor);
                    quot = ring.divideExact(remainder.lc(), dlc);
                }

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = ring.getZero();
        }

        return new UnivariatePolynomial[]{UnivariatePolynomial.create(ring, quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E>[] pseudoDivideAndRemainderLinearDividerAdaptive(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method
        Ring<E> ring = dividend.ring;
        E cc = ring.negate(divider.cc()), lc = divider.lc(), factor = ring.getOne();
        E[] quotient = copy ? ring.createArray(dividend.degree) : dividend.data;
        E res = ring.getZero();
        for (int i = dividend.degree; ; --i) {
            E tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = ring.add(ring.multiply(res, cc), ring.multiply(factor, tmp));
            if (i == 0) break;
            E quot = ring.divideOrNull(res, lc);
            if (quot == null) {
                E gcd = ring.gcd(res, lc), f = ring.divideExact(lc, gcd);
                factor = ring.multiply(factor, f);
                res = ring.multiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = ring.multiply(quotient[j], f);
                quot = ring.divideExact(res, lc);
            }
            res = quot;
        }
        if (!copy) quotient[dividend.degree] = ring.getZero();
        return new UnivariatePolynomial[]{UnivariatePolynomial.create(ring, quotient), UnivariatePolynomial.create(ring, res)};
    }

    /**
     * Returns quotient and remainder using adaptive pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E> pseudoRemainderAdaptive(final UnivariatePolynomial<E> dividend,
                                                               final UnivariatePolynomial<E> divider,
                                                               final boolean copy) {
        checkZeroDivider(divider);
        if (dividend.isZero())
            return UnivariatePolynomial.zero(dividend.ring);
        if (dividend.degree < divider.degree)
            return copy ? dividend.clone() : dividend;
        if (divider.degree == 0)
            return UnivariatePolynomial.zero(dividend.ring);
        if (divider.degree == 1)
            return pseudoRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E> pseudoRemainderAdaptive0(final UnivariatePolynomial<E> dividend,
                                                                final UnivariatePolynomial<E> divider,
                                                                final boolean copy) {
        assert dividend.degree >= divider.degree;

        Ring<E> ring = dividend.ring;
        UnivariatePolynomial<E> remainder = copy ? dividend.clone() : dividend;

        E dlc = divider.lc();
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                E quot = ring.divideOrNull(remainder.lc(), dlc);
                if (quot == null) {
                    E gcd = ring.gcd(remainder.lc(), divider.lc());
                    E factor = ring.divideExact(divider.lc(), gcd);
                    remainder.multiply(factor);
                    quot = ring.divideExact(remainder.lc(), dlc);
                }

                remainder.subtract(divider, quot, i);

            }
        }

        return remainder;
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    @SuppressWarnings("unchecked")
    static <E> UnivariatePolynomial<E> pseudoRemainderLinearDividerAdaptive(UnivariatePolynomial<E> dividend, UnivariatePolynomial<E> divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method
        Ring<E> ring = dividend.ring;
        E cc = ring.negate(divider.cc()), lc = divider.lc(), factor = ring.getOne();
        E res = ring.getZero();
        for (int i = dividend.degree; ; --i) {
            E tmp = dividend.data[i];
            res = ring.add(ring.multiply(res, cc), ring.multiply(factor, tmp));
            if (i == 0) break;
            E quot = ring.divideOrNull(res, lc);
            if (quot == null) {
                E gcd = ring.gcd(res, lc), f = ring.divideExact(lc, gcd);
                factor = ring.multiply(factor, f);
                res = ring.multiply(res, f);
                quot = ring.divideExact(res, lc);
            }
            res = quot;
        }
        return UnivariatePolynomial.create(ring, res);
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
        private static final long serialVersionUID = 1L;
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
        return dividend.createArray(q, r);
    }

    /* ********************************* Machine-precision fast division in Zp[x]  ******************************** */

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static UnivariatePolynomialZp64[] divideAndRemainderFast(UnivariatePolynomialZp64 dividend,
                                                                    UnivariatePolynomialZp64 divider,
                                                                    boolean copy) {
        UnivariatePolynomialZp64[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static UnivariatePolynomialZp64[] divideAndRemainderFast(UnivariatePolynomialZp64 dividend,
                                                                    UnivariatePolynomialZp64 divider,
                                                                    InverseModMonomial<UnivariatePolynomialZp64> invMod,
                                                                    boolean copy) {
        UnivariatePolynomialZp64[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy);
    }

    static UnivariatePolynomialZp64[] divideAndRemainderFastCorrectLC(UnivariatePolynomialZp64 dividend,
                                                                      UnivariatePolynomialZp64 divider,
                                                                      InverseModMonomial<UnivariatePolynomialZp64> invMod,
                                                                      boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, invMod, copy);

        long lc = divider.lc();
        long lcInv = divider.ring.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        UnivariatePolynomialZp64[] result = divideAndRemainderFast0(dividend, divider, invMod, copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    static UnivariatePolynomialZp64[] divideAndRemainderFast0(UnivariatePolynomialZp64 dividend,
                                                              UnivariatePolynomialZp64 divider,
                                                              boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.isMonic())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        long lc = divider.lc();
        long lcInv = divider.ring.reciprocal(lc);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        UnivariatePolynomialZp64[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
        E lcInv = dividend.ring.reciprocal(lc);
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
        E lcInv = dividend.ring.reciprocal(lc);
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
    private static UnivariatePolynomialZp64 earlyRemainderChecks(final UnivariatePolynomialZp64 dividend,
                                                                 final UnivariatePolynomialZp64 divider,
                                                                 final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1) {
            IntegersZp64 ring = dividend.ring;
            return dividend.createFromArray(new long[]{
                    dividend.evaluate(
                            ring.multiply(ring.negate(divider.cc()), ring.reciprocal(divider.lc())))
            });
        }
        return null;
    }

    /**
     * Returns remainder of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static UnivariatePolynomialZp64 remainder(final UnivariatePolynomialZp64 dividend,
                                                     final UnivariatePolynomialZp64 divider,
                                                     final boolean copy) {
        UnivariatePolynomialZp64 rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static UnivariatePolynomialZp64 remainderClassical0(final UnivariatePolynomialZp64 dividend,
                                                        final UnivariatePolynomialZp64 divider,
                                                        final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.assertSameCoefficientRingWith(divider);

        UnivariatePolynomialZp64 remainder = copy ? dividend.clone() : dividend;
        long lcInverse = dividend.ring.reciprocal(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, remainder.ring.multiply(remainder.lc(), lcInverse), i);

        return remainder;
    }

    /**
     * Fast remainder using Newton's iteration with switch to classical remainder for small polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   pre-conditioned divider ({@code fastDivisionPreConditioning(divider)})
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static UnivariatePolynomialZp64 remainderFast(final UnivariatePolynomialZp64 dividend,
                                                         final UnivariatePolynomialZp64 divider,
                                                         final InverseModMonomial<UnivariatePolynomialZp64> invMod,
                                                         final boolean copy) {
        UnivariatePolynomialZp64 rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        return divideAndRemainderFastCorrectLC(dividend, divider, invMod, copy)[1];
    }

    /**
     * Returns quotient of dividing {@code dividend} by {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static UnivariatePolynomialZp64 quotient(final UnivariatePolynomialZp64 dividend,
                                                    final UnivariatePolynomialZp64 divider,
                                                    final boolean copy) {
        UnivariatePolynomialZp64[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static UnivariatePolynomialZp64 quotientFast(final UnivariatePolynomialZp64 dividend,
                                                        final UnivariatePolynomialZp64 divider,
                                                        final InverseModMonomial<UnivariatePolynomialZp64> invMod,
                                                        final boolean copy) {
        UnivariatePolynomialZp64[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
            E p = dividend.ring.divideOrNull(dividend.ring.negate(divider.cc()), divider.lc());
            if (p == null)
                return null;
            return dividend.createConstant(dividend.evaluate(p));
        }
        return null;
    }

    /**
     * Returns remainder of {@code dividend} and {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
        dividend.assertSameCoefficientRingWith(divider);

        UnivariatePolynomial<E> remainder = copy ? dividend.clone() : dividend;
        Ring<E> ring = dividend.ring;
        if (ring.isField()) {
            E lcInverse = ring.reciprocal(divider.lc());
            for (int i = dividend.degree - divider.degree; i >= 0; --i)
                if (remainder.degree == divider.degree + i)
                    remainder.subtract(divider, ring.multiply(remainder.lc(), lcInverse), i);
        } else {
            for (int i = dividend.degree - divider.degree; i >= 0; --i)
                if (remainder.degree == divider.degree + i) {
                    E quot = ring.divideOrNull(remainder.lc(), divider.lc());
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
     * Returns quotient of {@code dividend} and {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
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
     * Returns quotient and remainder of {@code dividend} and {@code divider} using pseudo division.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] pseudoDivideAndRemainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof UnivariatePolynomialZ64)
            return (Poly[]) pseudoDivideAndRemainder((UnivariatePolynomialZ64) dividend, (UnivariatePolynomialZ64) divider, copy);
        if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly[]) divideAndRemainder((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, copy);
        else if (dividend instanceof UnivariatePolynomial) {
            if (dividend.isOverField())
                return (Poly[]) divideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
            else
                return (Poly[]) pseudoDivideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        } else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Returns {@code {quotient, remainder}} of {@code dividend} and {@code divider} or {@code null} if the division is
     * not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] divideAndRemainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof UnivariatePolynomialZ64)
            return (Poly[]) divideAndRemainder((UnivariatePolynomialZ64) dividend, (UnivariatePolynomialZ64) divider, copy);
        else if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly[]) divideAndRemainder((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly[]) divideAndRemainder((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Returns {@code {quotient, remainder}} of {@code dividend} and {@code divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @param invMod   precomputed Newton inverses
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] divideAndRemainderFast(Poly dividend, Poly divider, InverseModMonomial<Poly> invMod, boolean copy) {
        if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly[]) divideAndRemainderFast((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, (InverseModMonomial) invMod, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly[]) divideAndRemainderFast((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, (InverseModMonomial) invMod, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Divides {@code dividend} by {@code divider} or throws {@code ArithmeticException} if exact division is not
     * possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider}
     * @throws ArithmeticException if exact division is not possible
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly divideExact(Poly dividend, Poly divider, boolean copy) {
        Poly[] qr = divideAndRemainder(dividend, divider, copy);
        if (qr == null || !qr[1].isZero())
            throw new ArithmeticException("Not divisible: (" + dividend + ") / (" + divider + ")");
        return qr[0];
    }

    /**
     * Divides {@code dividend} by {@code divider} or returns {@code null} if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider} or {@code null} if exact division is not possible
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly divideOrNull(Poly dividend, Poly divider, boolean copy) {
        Poly[] qr = divideAndRemainder(dividend, divider, copy);
        if (qr == null || !qr[1].isZero())
            return null;
        return qr[0];
    }

    /**
     * Returns remainder of {@code dividend} and {@code divider}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly remainder(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof UnivariatePolynomialZ64)
            return (Poly) remainder((UnivariatePolynomialZ64) dividend, (UnivariatePolynomialZ64) divider, copy);
        else if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly) remainder((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly quotient(Poly dividend, Poly divider, boolean copy) {
        if (dividend instanceof UnivariatePolynomialZ64)
            return (Poly) quotient((UnivariatePolynomialZ64) dividend, (UnivariatePolynomialZ64) divider, copy);
        else if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly) quotient((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, copy);
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
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to {@code
     *                 dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly remainderFast(final Poly dividend,
                                                                                final Poly divider,
                                                                                final InverseModMonomial<Poly> invMod,
                                                                                final boolean copy) {
        if (dividend instanceof UnivariatePolynomialZp64)
            return (Poly) remainderFast((UnivariatePolynomialZp64) dividend, (UnivariatePolynomialZp64) divider, (InverseModMonomial<UnivariatePolynomialZp64>) invMod, copy);
        else if (dividend instanceof UnivariatePolynomial)
            return (Poly) remainderFast((UnivariatePolynomial) dividend, (UnivariatePolynomial) divider, (InverseModMonomial) invMod, copy);
        else
            throw new RuntimeException(dividend.getClass().toString());
    }

    /**
     * Gives an upper bound on the coefficients of remainder of division of {@code dividend} by {@code divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return upper bound on the coefficients of remainder
     */
    public static <E> E remainderCoefficientBound(UnivariatePolynomial<E> dividend,
                                                  UnivariatePolynomial<E> divider) {
        if (divider.degree < dividend.degree)
            return dividend.maxAbsCoefficient();
        Ring<E> ring = dividend.ring;
        // see e.g. http://www.csd.uwo.ca/~moreno//AM583/Lectures/Newton2Hensel.html/node13.html
        return ring.multiply(dividend.maxAbsCoefficient(),
                ring.pow(ring.increment(ring.quotient(divider.maxAbsCoefficient(), divider.lc())),
                        dividend.degree - divider.degree + 1));
    }
}
