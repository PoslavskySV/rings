package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;

import java.util.ArrayList;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.number.BigIntegerArithmetics.*;

/**
 * Algorithms for division with remainder.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class bDivisionWithRemainder {
    private bDivisionWithRemainder() {}

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
        if (dividend.isZero())
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), bMutablePolynomialZ.zero()};
        if (dividend.degree < divider.degree)
            return new bMutablePolynomialZ[]{bMutablePolynomialZ.zero(), copy ? dividend.clone() : dividend};
        BigInteger factor = safePow(divider.lc(), dividend.degree - divider.degree + 1);
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

            } else quotient[i] = 0;
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
            } else quotient[i] = 0;
        }
        return new bMutablePolynomialZ[]{bMutablePolynomialZ.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] divideAndRemainderLinearDivider(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, ONE, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] pseudoDivideAndRemainderLinearDivider(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, safePow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialZ[] divideAndRemainderLinearDivider0(bMutablePolynomialZ dividend, bMutablePolynomialZ divider, BigInteger raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        BigInteger cc = -divider.cc(), lc = divider.lc();
        BigInteger[] quotient = copy ? new BigInteger[dividend.degree] : dividend.data;
        BigInteger res = ZERO;
        for (int i = dividend.degree; ; --i) {
            BigInteger tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = safeAdd(safeMultiply(res, cc), safeMultiply(raiseFactor, tmp));
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

    /**
     * Returns the remainder of {@code dividend} divided by monomial {@code x^xDegree}
     *
     * @param dividend the dividend
     * @param xDegree  monomial degree
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static <T extends bMutablePolynomialAbstract<T>> T remainderMonomial(T dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).truncate(xDegree - 1);
    }


    /* ********************************** Modular methods ********************************** */


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(bMutablePolynomialMod dividend,
                                                bMutablePolynomialMod divider) {
        // practical benchmarks show that without pre-conditioning,
        // classical division is always faster or at least the same fast
        return true;
    }

    /** early checks for division */
    private static bMutablePolynomialMod[] earlyDivideAndRemainderChecks(final bMutablePolynomialMod dividend,
                                                                         final bMutablePolynomialMod divider,
                                                                         final boolean copy) {
        if (dividend.isZero())
            return new bMutablePolynomialMod[]{dividend.createZero(), dividend.createZero()};
        if (dividend.degree < divider.degree)
            return new bMutablePolynomialMod[]{dividend.createZero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new bMutablePolynomialMod[]{(copy ? dividend.clone() : dividend).multiply(modInverse(divider.lc(), dividend.modulus)), dividend.createZero()};
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
    public static bMutablePolynomialMod[] divideAndRemainder(final bMutablePolynomialMod dividend,
                                                             final bMutablePolynomialMod divider,
                                                             final boolean copy) {
        bMutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static bMutablePolynomialMod[] divideAndRemainderClassic(final bMutablePolynomialMod dividend,
                                                                    final bMutablePolynomialMod divider,
                                                                    final boolean copy) {
        bMutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static bMutablePolynomialMod[] divideAndRemainderClassic0(final bMutablePolynomialMod dividend,
                                                              final bMutablePolynomialMod divider,
                                                              final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        bMutablePolynomialMod remainder = copy ? dividend.clone() : dividend;
        BigInteger[] quotient = new BigInteger[dividend.degree - divider.degree + 1];

        BigInteger lcInverse = modInverse(divider.lc(), dividend.modulus);
        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = remainder.multiplyMod(remainder.lc(), lcInverse);
                remainder.subtract(divider, quotient[i], i);
            } else quotient[i] = 0;
        }

        return new bMutablePolynomialMod[]{dividend.createFromArray(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static bMutablePolynomialMod[] divideAndRemainderLinearDividerModulus(bMutablePolynomialMod dividend, bMutablePolynomialMod divider, boolean copy) {
        assert divider.degree == 1;
        assert dividend.degree > 0;
        dividend.checkCompatibleModulus(divider);

        //apply Horner's method

        BigInteger cc = mod(divider.cc().negate(), dividend.modulus);
        BigInteger lcInverse = modInverse(divider.lc(), dividend.modulus);

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
        return new bMutablePolynomialMod[]{dividend.createFromArray(quotient), dividend.createFromArray(new BigInteger[]{res})};
    }

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomial {
        final bMutablePolynomialMod poly;

        private InverseModMonomial(bMutablePolynomialMod poly) {
            if (!poly.cc().isOne())
                throw new IllegalArgumentException("Smallest coefficient is not a unit: " + poly);
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<bMutablePolynomialMod> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public bMutablePolynomialMod getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() >= r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            bMutablePolynomialMod gPrev = currentSize == 0 ? poly.createOne() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                bMutablePolynomialMod tmp = gPrev.clone().multiply(TWO).subtract(gPrev.clone().square().multiply(poly));
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
    public static InverseModMonomial fastDivisionPreConditioning(bMutablePolynomialMod divider) {
        if (!divider.lc().isOne())
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
    public static bMutablePolynomialMod[] divideAndRemainderFast(bMutablePolynomialMod dividend,
                                                                 bMutablePolynomialMod divider,
                                                                 boolean copy) {
        bMutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static bMutablePolynomialMod[] divideAndRemainderFast(bMutablePolynomialMod dividend,
                                                                 bMutablePolynomialMod divider,
                                                                 InverseModMonomial invMod,
                                                                 boolean copy) {
        bMutablePolynomialMod[] r = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invMod, copy);
    }

    static bMutablePolynomialMod[] divideAndRemainderFast0(bMutablePolynomialMod dividend,
                                                           bMutablePolynomialMod divider,
                                                           boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.lc().isOne())
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);

        BigInteger lc = divider.lc();
        BigInteger lcInv = modInverse(lc, dividend.modulus);
        // make the divisor monic
        divider.multiply(lcInv);
        // perform fast arithmetic with monic divisor
        bMutablePolynomialMod[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider), copy);
        // reconstruct divisor's lc
        divider.multiply(lc);
        // reconstruct actual quotient
        result[0].multiply(lcInv);
        return result;
    }

    static bMutablePolynomialMod[] divideAndRemainderFast0(bMutablePolynomialMod dividend,
                                                           bMutablePolynomialMod divider,
                                                           InverseModMonomial invRevMod,
                                                           boolean copy) {
        dividend.checkCompatibleModulus(divider);
        int m = dividend.degree - divider.degree;
        bMutablePolynomialMod q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1)), m + 1, false).reverse();
        if (q.degree < m)
            q.shiftRight(m - q.degree);
        return new bMutablePolynomialMod[]{q, (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(q))};
    }

    /** fast division checks */
    private static bMutablePolynomialMod earlyRemainderChecks(final bMutablePolynomialMod dividend,
                                                              final bMutablePolynomialMod divider,
                                                              final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend);
        if (divider.degree == 0)
            return dividend.createZero();
        if (divider.degree == 1)
            return dividend.createFromArray(new BigInteger[]{
                    dividend.evaluate(
                            dividend.multiplyMod(mod(divider.cc().negate(), dividend.modulus), modInverse(divider.lc(), dividend.modulus)))
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
    public static bMutablePolynomialMod remainder(final bMutablePolynomialMod dividend,
                                                  final bMutablePolynomialMod divider,
                                                  final boolean copy) {
        bMutablePolynomialMod rem = earlyRemainderChecks(dividend, divider, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, copy);

        return divideAndRemainderFast0(dividend, divider, copy)[1];
    }

    /** Plain school implementation */
    static bMutablePolynomialMod remainderClassical0(final bMutablePolynomialMod dividend,
                                                     final bMutablePolynomialMod divider,
                                                     final boolean copy) {
        assert dividend.degree >= divider.degree;
        dividend.checkCompatibleModulus(divider);

        bMutablePolynomialMod remainder = copy ? dividend.clone() : dividend;
        BigInteger lcInverse = modInverse(divider.lc(), dividend.modulus);
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
    public static bMutablePolynomialMod remainderFast(final bMutablePolynomialMod dividend,
                                                      final bMutablePolynomialMod divider,
                                                      final InverseModMonomial invMod,
                                                      final boolean copy) {
        bMutablePolynomialMod rem = earlyRemainderChecks(dividend, divider, copy);
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
    public static bMutablePolynomialMod quotient(final bMutablePolynomialMod dividend,
                                                 final bMutablePolynomialMod divider,
                                                 final boolean copy) {
        bMutablePolynomialMod[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
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
    public static bMutablePolynomialMod quotientFast(final bMutablePolynomialMod dividend,
                                                     final bMutablePolynomialMod divider,
                                                     final InverseModMonomial invMod,
                                                     final boolean copy) {
        bMutablePolynomialMod[] qd = earlyDivideAndRemainderChecks(dividend, divider, copy);
        if (qd != null)
            return qd[0];

        return divideAndRemainderFast0(dividend, divider, invMod, copy)[0];
    }


    /* ********************************** Common conversion ********************************** */

    @SuppressWarnings("unchecked")
    public static <T extends bMutablePolynomialAbstract<T>> T[] divideAndRemainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof bMutablePolynomialZ)
            return (T[]) divideAndRemainder((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T[]) divideAndRemainder((bMutablePolynomialMod) dividend, (bMutablePolynomialMod) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends bMutablePolynomialAbstract<T>> T remainder(T dividend, T divider, boolean copy) {
        if (dividend instanceof bMutablePolynomialZ)
            return (T) remainder((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T) remainder((bMutablePolynomialMod) dividend, (bMutablePolynomialMod) divider, copy);
    }

    @SuppressWarnings("unchecked")
    public static <T extends bMutablePolynomialAbstract<T>> T quotient(T dividend, T divider, boolean copy) {
        if (dividend instanceof bMutablePolynomialZ)
            return (T) quotient((bMutablePolynomialZ) dividend, (bMutablePolynomialZ) divider, copy);
        else
            return (T) quotient((bMutablePolynomialMod) dividend, (bMutablePolynomialMod) divider, copy);
    }
}
