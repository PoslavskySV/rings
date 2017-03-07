package cc.r2.core.poly2;

import cc.r2.core.poly2.LongModularArithmetics.MagicDivider;

import java.util.ArrayList;

import static cc.r2.core.poly2.LongArithmetics.*;
import static cc.r2.core.poly2.LongModularArithmetics.divideSignedFast;
import static cc.r2.core.poly2.LongModularArithmetics.magicSigned;

/**
 * Algorithms for division with remainder.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class DivisionWithRemainder {
    private DivisionWithRemainder() {}

    /**
     * Returns {@code {quotient, remainder}} or {@code null} if the division is not possible.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder} or {@code null} if the division is not possible
     */
    public static MutablePolynomialZ[] divideAndRemainder(final MutablePolynomialZ dividend,
                                                          final MutablePolynomialZ divider,
                                                          boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), MutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            MutablePolynomialZ div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new MutablePolynomialZ[]{div, MutablePolynomialZ.zero()};
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
    public static MutablePolynomialZ[] pseudoDivideAndRemainder(final MutablePolynomialZ dividend,
                                                                final MutablePolynomialZ divider,
                                                                final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), MutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        long factor = safePow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new MutablePolynomialZ[]{(copy ? dividend.clone() : dividend).multiply(factor / dividend.lc()), MutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderGeneral0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static MutablePolynomialZ[] divideAndRemainderGeneral0(final MutablePolynomialZ dividend,
                                                           final MutablePolynomialZ divider,
                                                           final long dividendRaiseFactor,
                                                           final boolean copy) {
        assert dividend.degree >= divider.degree;
        if (divider.lc() == 1 && dividendRaiseFactor == 1)
            return divideAndRemainderGeneralMonic(dividend, divider, copy);

        MutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        long[] quotient = new long[dividend.degree - divider.degree + 1];


        MagicDivider magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc())
                    return null;

                quotient[i] = quot;
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new MutablePolynomialZ[]{MutablePolynomialZ.create(quotient), remainder};
    }

    /** Plain school implementation */
    private static MutablePolynomialZ[] divideAndRemainderGeneralMonic(final MutablePolynomialZ dividend,
                                                                       final MutablePolynomialZ divider,
                                                                       final boolean copy) {
        assert divider.lc() == 1;

        MutablePolynomialZ
                remainder = (copy ? dividend.clone() : dividend);
        long[] quotient = new long[dividend.degree - divider.degree + 1];
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.lc();
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }
        return new MutablePolynomialZ[]{MutablePolynomialZ.create(quotient), remainder};
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
    static MutablePolynomialZ[] pseudoDivideAndRemainderAdaptive(final MutablePolynomialZ dividend,
                                                                 final MutablePolynomialZ divider,
                                                                 final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), MutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomialZ[]{MutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new MutablePolynomialZ[]{copy ? dividend.clone() : dividend, MutablePolynomialZ.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static MutablePolynomialZ[] pseudoDivideAndRemainderAdaptive0(final MutablePolynomialZ dividend,
                                                                  final MutablePolynomialZ divider,
                                                                  final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomialZ remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        MagicDivider magic = magicSigned(divider.lc());
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                long quot = divideSignedFast(remainder.lc(), magic);
                if (quot * divider.lc() != remainder.lc()) {
                    long gcd = gcd(remainder.lc(), divider.lc());
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

        return new MutablePolynomialZ[]{MutablePolynomialZ.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomialZ[] pseudoDivideAndRemainderLinearDividerAdaptive(MutablePolynomialZ dividend, MutablePolynomialZ divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc(), factor = 1;
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        MagicDivider magic = magicSigned(lc);
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = safeAdd(safeMultiply(res, cc), safeMultiply(factor, tmp));
            if (i == 0) break;
            long quot = divideSignedFast(res, magic);
            if (quot * lc != res) {
                long gcd = gcd(res, lc), f = lc / gcd;
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
        return new MutablePolynomialZ[]{MutablePolynomialZ.create(quotient), MutablePolynomialZ.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomialZ[] divideAndRemainderLinearDivider(MutablePolynomialZ dividend, MutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomialZ[] pseudoDivideAndRemainderLinearDivider(MutablePolynomialZ dividend, MutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomialZ[] divideAndRemainderLinearDivider0(MutablePolynomialZ dividend, MutablePolynomialZ divider, long raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc();
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        MagicDivider magic = magicSigned(lc);
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
        return new MutablePolynomialZ[]{MutablePolynomialZ.create(quotient), MutablePolynomialZ.create(res)};
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
    public static MutablePolynomialZ remainder(final MutablePolynomialZ dividend,
                                               final MutablePolynomialZ divider,
                                               final boolean copy) {
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return MutablePolynomialZ.zero();
        if (divider.degree == 1) {
            if (divider.cc() % divider.lc() != 0)
                return null;
            return MutablePolynomialZ.create(dividend.evaluate(-divider.cc() / divider.lc()));
        }
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static MutablePolynomialZ remainder0(final MutablePolynomialZ dividend,
                                         final MutablePolynomialZ divider,
                                         final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomialZ remainder = copy ? dividend.clone() : dividend;
        MagicDivider magic = magicSigned(remainder.lc());
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
     * Returns the remainder of {@code dividend} divided by monomial {@code x^xDegree}
     *
     * @param dividend the dividend
     * @param xDegree  monomial degree
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static <T extends MutablePolynomialAbstract<T>> T remainderMonomial(T dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).truncate(xDegree - 1);
    }


    /* ********************************** Modular methods ********************************** */


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(MutablePolynomialMod dividend,
                                                MutablePolynomialMod divider) {
        // practical benchmarks show that without pre-conditioning,
        // classical division is always faster or at least the same fast
        return true;
    }

    /** early checks for division */
    private static MutablePolynomialMod[] earlyDivideAndRemainderChecks(final MutablePolynomialMod dividend,
                                                                        final MutablePolynomialMod divider,
                                                                        final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomialMod[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomialMod[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new MutablePolynomialMod[]{(copy ? dividend.clone() : dividend).multiply(modInverse(divider.lc(), dividend.modulus)), dividend.createZero()};
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
    public static MutablePolynomialMod[] divideAndRemainder(final MutablePolynomialMod dividend,
                                                            final MutablePolynomialMod divider,
                                                            final boolean copy) {
        MutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static MutablePolynomialMod[] divideAndRemainderClassic(final MutablePolynomialMod dividend,
                                                                   final MutablePolynomialMod divider,
                                                                   final boolean copy) {
        MutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static MutablePolynomialMod[] divideAndRemainderClassic0(final MutablePolynomialMod dividend,
                                                             final MutablePolynomialMod divider,
                                                             final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        MutablePolynomialMod remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        long lcInverse = modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.multiplyMod(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new MutablePolynomialMod[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomialMod[] divideAndRemainderLinearDividerModulus(MutablePolynomialMod dividend, MutablePolynomialMod divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatibleModulus(divider);

        //apply Horner's method

        long cc = mod(-divider.cc(), dividend.modulus);
        long lcInverse = modInverse(divider.lc(), dividend.modulus);

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
        return new MutablePolynomialMod[]{dividend.createFromArray(quotient), dividend.createFromArray(new long[]{res})};
    }

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomial {
        final MutablePolynomialMod poly;

        private InverseModMonomial(MutablePolynomialMod poly) {
            if (poly.cc() != 1)
                throw new IllegalArgumentException("Smallest coefficient is not a unit: " + poly);
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<MutablePolynomialMod> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public MutablePolynomialMod getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() >= r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            MutablePolynomialMod gPrev = currentSize == 0 ? poly.createOne() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                MutablePolynomialMod tmp = gPrev.clone().multiply(2).subtract(gPrev.clone().square().multiply(poly));
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
    public static InverseModMonomial fastDivisionPreConditioning(MutablePolynomialMod divider) {
        if (divider.lc() != 1)
            throw new IllegalArgumentException("Only monic polynomials allowed. Input: " + divider);
        return new InverseModMonomial(divider.clone().reverse());
    }

    /**
     * Fast algorithm for division with remainder using Newton's iteration.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutablePolynomialMod[] divideAndRemainderFast(MutablePolynomialMod dividend,
                                                                MutablePolynomialMod divider,
                                                                boolean copy) {
        MutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static MutablePolynomialMod[] divideAndRemainderFast(MutablePolynomialMod dividend,
                                                                MutablePolynomialMod divider,
                                                                InverseModMonomial invMod,
                                                                boolean copy) {
        MutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invMod, copy);
    }

    static MutablePolynomialMod[] divideAndRemainderFast0(MutablePolynomialMod dividend,
                                                          MutablePolynomialMod divider,
                                                          boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.lc() == 1)
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        long lc = divider.lc();
        long lcInv = modInverse(lc, dividend.modulus);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        MutablePolynomialMod[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    static MutablePolynomialMod[] divideAndRemainderFast0(MutablePolynomialMod dividend,
                                                          MutablePolynomialMod divider,
                                                          InverseModMonomial invRevMod,
                                                          boolean copy) {
        dividend.checkCompatibleModulus(divider);
        int m = dividend.degree - divider.degree;
        MutablePolynomialMod q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1)), m + 1, false).reverse();
        if (q.degree < m)
            q.shiftRight(m - q.degree);
        return new MutablePolynomialMod[]{q, (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(q))};
    }

    /** fast division checks */
    private static MutablePolynomialMod earlyRemainderChecks(final MutablePolynomialMod dividend,
                                                             final MutablePolynomialMod divider,
                                                             final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1)
            return dividend.createFromArray(new long[]{
                    dividend.evaluate(
                            dividend.multiplyMod(mod(-divider.cc(), dividend.modulus), modInverse(divider.lc(), dividend.modulus)))
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
    public static MutablePolynomialMod remainder(final MutablePolynomialMod dividend,
                                                 final MutablePolynomialMod divider,
                                                 final boolean copy) {
        MutablePolynomialMod rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static MutablePolynomialMod remainderClassical0(final MutablePolynomialMod dividend,
                                                    final MutablePolynomialMod divider,
                                                    final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        MutablePolynomialMod remainder = copy ? dividend.clone() : dividend;
        long lcInverse = modInverse(divider.lc(), dividend.modulus);
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
    public static MutablePolynomialMod remainderFast(final MutablePolynomialMod dividend,
                                                     final MutablePolynomialMod divider,
                                                     final InverseModMonomial invMod,
                                                     final boolean copy) {
        MutablePolynomialMod rem = earlyRemainderChecks(dividend, divider, copy);
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
    public static MutablePolynomialMod quotient(final MutablePolynomialMod dividend,
                                                final MutablePolynomialMod divider,
                                                final boolean copy) {
        MutablePolynomialMod[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static MutablePolynomialMod quotientFast(final MutablePolynomialMod dividend,
                                                    final MutablePolynomialMod divider,
                                                    final InverseModMonomial invMod,
                                                    final boolean copy) {
        MutablePolynomialMod[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[0];
    }


    /* ********************************** Common conversion ********************************** */

    @SuppressWarnings("unchecked")
    public static <T extends MutablePolynomialAbstract<T>> T[] divideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof MutablePolynomialZ)
            return (T[]) divideAndRemainder((MutablePolynomialZ) dividend, (MutablePolynomialZ) divider, copy);
        else
            return (T[]) divideAndRemainder((MutablePolynomialMod) dividend, (MutablePolynomialMod) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends MutablePolynomialAbstract<T>> T remainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof MutablePolynomialZ)
            return (T) remainder((MutablePolynomialZ) dividend, (MutablePolynomialZ) divider, copy);
        else
            return (T) remainder((MutablePolynomialMod) dividend, (MutablePolynomialMod) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends MutablePolynomialAbstract<T>> T quotient(T dividend, T divider, boolean copy) {
        if (dividend instanceof MutablePolynomialZ)
            return (T) quotient((MutablePolynomialZ) dividend, (MutablePolynomialZ) divider, copy);
        else
            return (T) quotient((MutablePolynomialMod) dividend, (MutablePolynomialMod) divider, copy);
    }
}
