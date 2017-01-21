package cc.r2.core.polynomial;

import java.util.ArrayList;

import static cc.r2.core.polynomial.LongArithmetics.*;

public final class DivideAndRemainder {
    private DivideAndRemainder() {}

    /**
     * Returns quotient and remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutablePolynomial[] divideAndRemainder(final MutablePolynomial dividend,
                                                         final MutablePolynomial divider,
                                                         boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomial[]{MutablePolynomial.zero(), MutablePolynomial.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomial[]{MutablePolynomial.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0) {
            MutablePolynomial div = copy ? dividend.clone() : dividend;
            div = div.divideOrNull(divider.lc());
            if (div == null) return null;
            return new MutablePolynomial[]{div, MutablePolynomial.zero()};
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
    public static MutablePolynomial[] pseudoDivideAndRemainder(final MutablePolynomial dividend,
                                                               final MutablePolynomial divider,
                                                               final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomial[]{MutablePolynomial.zero(), MutablePolynomial.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomial[]{MutablePolynomial.zero(), copy ? dividend.clone() : dividend};
        long factor = pow(divider.lc(), dividend.degree - divider.degree + 1);
        if (divider.degree == 0)
            return new MutablePolynomial[]{(copy ? dividend.clone() : dividend).multiply(factor / dividend.lc()), MutablePolynomial.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDivider0(dividend, divider, factor, copy);
        return divideAndRemainderGeneral0(dividend, divider, factor, copy);
    }

    /** Plain school implementation */
    static MutablePolynomial[] divideAndRemainderGeneral0(final MutablePolynomial dividend,
                                                          final MutablePolynomial divider,
                                                          final long dividendRaiseFactor,
                                                          final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomial
                remainder = (copy ? dividend.clone() : dividend).multiply(dividendRaiseFactor);
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0)
                    return null;

                quotient[i] = remainder.lc() / divider.lc();
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new MutablePolynomial[]{MutablePolynomial.create(quotient), remainder};
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
    static MutablePolynomial[] pseudoDivideAndRemainderAdaptive(final MutablePolynomial dividend,
                                                                final MutablePolynomial divider,
                                                                final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomial[]{MutablePolynomial.zero(), MutablePolynomial.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomial[]{MutablePolynomial.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new MutablePolynomial[]{copy ? dividend.clone() : dividend, MutablePolynomial.zero()};
        if (divider.degree == 1)
            return pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, copy);
        return pseudoDivideAndRemainderAdaptive0(dividend, divider, copy);
    }

    /** general implementation */
    static MutablePolynomial[] pseudoDivideAndRemainderAdaptive0(final MutablePolynomial dividend,
                                                                 final MutablePolynomial divider,
                                                                 final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomial remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0) {
                    long gcd = gcd(remainder.lc(), divider.lc());
                    long factor = divider.lc() / gcd;
                    remainder.multiply(factor);
                    for (int j = i + 1; j < quotient.length; ++j)
                        quotient[j] = multiply(quotient[j], factor);
                }

                quotient[i] = remainder.lc() / divider.lc();
                remainder.subtract(divider, quotient[i], i);

            } else quotient[i] = 0;
        }

        return new MutablePolynomial[]{MutablePolynomial.create(quotient), remainder};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomial[] pseudoDivideAndRemainderLinearDividerAdaptive(MutablePolynomial dividend, MutablePolynomial divider, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc(), factor = 1;
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(factor, tmp));
            if (i == 0) break;
            if (res % lc != 0) {
                long gcd = gcd(res, lc), f = lc / gcd;
                factor = multiply(factor, f);
                res = multiply(res, f);
                if (i != dividend.degree)
                    for (int j = quotient.length - 1; j >= i; --j)
                        quotient[j] = multiply(quotient[j], f);
            }
            res = res / lc;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutablePolynomial[]{MutablePolynomial.create(quotient), MutablePolynomial.create(res)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomial[] divideAndRemainderLinearDivider(MutablePolynomial dividend, MutablePolynomial divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, 1, copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomial[] pseudoDivideAndRemainderLinearDivider(MutablePolynomial dividend, MutablePolynomial divider, boolean copy) {
        return divideAndRemainderLinearDivider0(dividend, divider, pow(divider.lc(), dividend.degree), copy);
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomial[] divideAndRemainderLinearDivider0(MutablePolynomial dividend, MutablePolynomial divider, long raiseFactor, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = -divider.cc(), lc = divider.lc();
        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; ; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = res;
            res = add(multiply(res, cc), multiply(raiseFactor, tmp));
            if (i == 0) break;
            if (res % lc != 0) return null;
            res = res / lc;
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutablePolynomial[]{MutablePolynomial.create(quotient), MutablePolynomial.create(res)};
    }

    /**
     * Returns the remainder of {@code dividend} divided by {@code divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return the remainder
     */
    public static MutablePolynomial remainder(final MutablePolynomial dividend,
                                              final MutablePolynomial divider,
                                              final boolean copy) {
        if (dividend.degree < divider.degree)
            return dividend;
        if (divider.degree == 0)
            return MutablePolynomial.zero();
        if (divider.degree == 1) {
            if (divider.cc() % divider.lc() != 0)
                return null;
            return MutablePolynomial.create(dividend.evaluate(-divider.cc() / divider.lc()));
        }
        return remainder0(dividend, divider, copy);
    }

    /** Plain school implementation */
    static MutablePolynomial remainder0(final MutablePolynomial dividend,
                                        final MutablePolynomial divider,
                                        final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomial remainder = copy ? dividend.clone() : dividend;
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i) {
                if (remainder.lc() % divider.lc() != 0)
                    return null;
                remainder.subtract(divider, remainder.lc() / divider.lc(), i);
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
    public static MutablePolynomial remainderMonomial(MutablePolynomial dividend, int xDegree, boolean copy) {
        return (copy ? dividend.clone() : dividend).cut(xDegree - 1);
    }



    /* ********************************** Modular methods ********************************** */

    /** when to switch to fast division using Newton iterations */
    static final int
            /** divider.lc() == 1 (no need to additionally make monic etc) **/
            DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD_MONIC = 36,
            DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD = 60;
    /** deg(dividend) - deg(divider) when to switch to fast division using Newton iterations */
    static final int
            /** divider.lc() == 1 (no need to additionally make monic etc) **/
            DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD_MONIC = 6,
            DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD = 10,
            DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD_SMALL_DEG = 20;
    ;


    /** when to switch between classical and Newton's */
    private static boolean useClassicalDivision(MutablePolynomial dividend,
                                                MutablePolynomial divider) {
        assert dividend.degree >= divider.degree;

        if (dividend.degree < DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD_MONIC) // <- very small degrees
            // we always use plain division for very small polynomials
            return true;

        if (divider.isMonic())
            // we use plain division when deg(dividend) is almost equal to deg(divisor)
            return dividend.degree - divider.degree < DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD_MONIC;
        else
            // we use plain division for small polynomials
            if (dividend.degree < DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD) // <- small degrees
                //for small polynomial we require larger degree diff between divider and dividend
                return dividend.degree - divider.degree < DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD_SMALL_DEG;
            else // <- large degrees
                // we use plain division when deg(dividend) is almost equal to deg(divisor)
                return dividend.degree - divider.degree < DIVIDEND_DIVIDER_DEGREE_DIFF_FAST_DIVISION_THRESHOLD;
    }

    /** early checks for division */
    private static MutablePolynomial[] earlyDivideAndRemainderChecks(final MutablePolynomial dividend,
                                                                     final MutablePolynomial divider,
                                                                     final long modulus,
                                                                     final boolean copy) {
        if (dividend.isZero())
            return new MutablePolynomial[]{MutablePolynomial.zero(), MutablePolynomial.zero()};
        if (dividend.degree < divider.degree)
            return new MutablePolynomial[]{MutablePolynomial.zero(), copy ? dividend.clone() : dividend};
        if (divider.degree == 0)
            return new MutablePolynomial[]{(copy ? dividend.clone() : dividend).multiply(modInverse(divider.lc(), modulus), modulus), MutablePolynomial.zero()};
        if (divider.degree == 1)
            return divideAndRemainderLinearDividerModulus(dividend, divider, modulus, copy);
        return null;
    }

    /**
     * Returns quotient and remainder modulo {@code modulus}. Classical division algorithm and Newton's iterations used
     * depending on the degrees of input polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutablePolynomial[] divideAndRemainder(final MutablePolynomial dividend,
                                                         final MutablePolynomial divider,
                                                         final long modulus,
                                                         final boolean copy) {
        MutablePolynomial[] r = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (r != null)
            return r;

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic0(dividend, divider, modulus, copy);

        return divideAndRemainderFast0(dividend, divider, modulus, copy);
    }

    /**
     * Returns quotient and remainder modulo {@code modulus}.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutablePolynomial[] divideAndRemainderClassic(final MutablePolynomial dividend,
                                                                final MutablePolynomial divider,
                                                                final long modulus,
                                                                final boolean copy) {
        MutablePolynomial[] r = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (r != null)
            return r;
        return divideAndRemainderClassic0(dividend, divider, modulus, copy);
    }

    /** Plain school implementation */
    static MutablePolynomial[] divideAndRemainderClassic0(final MutablePolynomial dividend,
                                                          final MutablePolynomial divider,
                                                          final long modulus,
                                                          final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomial remainder = copy ? dividend.clone() : dividend;
        long[] quotient = new long[dividend.degree - divider.degree + 1];

        for (int i = dividend.degree - divider.degree; i >= 0; --i) {
            if (remainder.degree == divider.degree + i) {
                quotient[i] = divideMod(remainder.lc(), divider.lc(), modulus);
                remainder.subtract(divider, quotient[i], i, modulus);
            } else quotient[i] = 0;
        }

        return new MutablePolynomial[]{MutablePolynomial.create(quotient).modulus(modulus), remainder.modulus(modulus)};
    }

    /** Fast division with remainder for divider of the form f(x) = x - u **/
    static MutablePolynomial[] divideAndRemainderLinearDividerModulus(MutablePolynomial dividend, MutablePolynomial divider, long modulus, boolean copy) {
        assert divider.degree == 1;

        //apply Horner's method

        long cc = mod(-divider.cc(), modulus);
        long lcInverse = modInverse(divider.lc(), modulus);

        if (divider.lc() != 1)
            cc = mod(multiply(cc, lcInverse), modulus);

        long[] quotient = copy ? new long[dividend.degree] : dividend.data;
        long res = 0;
        for (int i = dividend.degree; i >= 0; --i) {
            long tmp = dividend.data[i];
            if (i != dividend.degree)
                quotient[i] = mod(multiply(res, lcInverse), modulus);
            res = addMod(multiply(res, cc), tmp, modulus);
        }
        if (!copy) quotient[dividend.degree] = 0;
        return new MutablePolynomial[]{MutablePolynomial.create(quotient), MutablePolynomial.create(res)};
    }

    /* that is [log2] */
    static int log2(int l) {
        if (l <= 0)
            throw new IllegalArgumentException();
        return 33 - Integer.numberOfLeadingZeros(l - 1);
    }

    /** Holds {@code poly^(-1) mod x^i } */
    public static final class InverseModMonomial {
        final long modulus;
        final MutablePolynomial poly;

        public InverseModMonomial(MutablePolynomial poly, long modulus) {
            if (poly.cc() != 1)
                throw new IllegalArgumentException("Smallest coefficient is not a unit: " + poly);
            this.modulus = modulus;
            this.poly = poly;
        }

        /** the inverses */
        private final ArrayList<MutablePolynomial> inverses = new ArrayList<>();

        /**
         * Returns {@code poly^(-1) mod x^xDegree }. Newton iterations are inside.
         *
         * @param xDegree monomial degree
         * @return {@code poly^(-1) mod x^xDegree }
         */
        public MutablePolynomial getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() > r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            MutablePolynomial gPrev = currentSize == 0 ? MutablePolynomial.one() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                MutablePolynomial tmp = gPrev.clone().multiply(2, modulus).subtract(gPrev.clone().square(modulus).multiply(poly, modulus), modulus);
                inverses.add(gPrev = remainderMonomial(tmp, 1 << i, false));
            }
            return gPrev;
        }
    }

    /**
     * Prepares {@code rev(divider)^(-1) mod x^i } for fast division.
     *
     * @param divider the divider
     * @param modulus the modulus
     */
    public static InverseModMonomial fastDivisionPreConditioning(MutablePolynomial divider, long modulus) {
        if (divider.lc() != 1)
            throw new IllegalArgumentException("Only monic polynomials allowed. Input: " + divider);
        return new InverseModMonomial(divider.clone().reverse(), modulus);
    }

    /**
     * Fast divide and remainder algorithm using Newton iterations.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return
     */
    public static MutablePolynomial[] divideAndRemainderFast(MutablePolynomial dividend,
                                                             MutablePolynomial divider,
                                                             long modulus,
                                                             boolean copy) {
        MutablePolynomial[] r = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, modulus, copy);
    }

    /**
     * Fast divide and remainder algorithm using Newton iterations.
     *
     * @param dividend  the dividend
     * @param divider   the divider
     * @param invRevMod {@code divider^(-1) mod x^i }
     * @param modulus   the modulus
     * @param copy      whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                  {@code dividend} and {@code dividend} data will be lost
     * @return {quotient, remainder}
     */
    public static MutablePolynomial[] divideAndRemainderFast(MutablePolynomial dividend,
                                                             MutablePolynomial divider,
                                                             InverseModMonomial invRevMod,
                                                             long modulus,
                                                             boolean copy) {
        MutablePolynomial[] r = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (r != null)
            return r;
        return divideAndRemainderFast0(dividend, divider, invRevMod, modulus, copy);
    }

    static MutablePolynomial[] divideAndRemainderFast0(MutablePolynomial dividend,
                                                       MutablePolynomial divider,
                                                       long modulus,
                                                       boolean copy) {
        // if the divider can be directly inverted modulo x^i
        if (divider.lc() == 1)
            return divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider, modulus), modulus, copy);

        long lc = divider.lc();
        long lcInv = LongArithmetics.modInverse(lc, modulus);
        // make the divisor monic
        divider.multiply(lcInv, modulus);
        // perform fast arithmetic with monic divisor
        MutablePolynomial[] result = divideAndRemainderFast0(dividend, divider, fastDivisionPreConditioning(divider, modulus), modulus, copy);
        // reconstruct divisor's lc
        divider.multiply(lc, modulus);
        // reconstruct actual quotient
        result[0].multiply(lcInv, modulus);
        return result;
    }

    static MutablePolynomial[] divideAndRemainderFast0(MutablePolynomial dividend,
                                                       MutablePolynomial divider,
                                                       InverseModMonomial invRevMod,
                                                       long modulus,
                                                       boolean copy) {
        int m = dividend.degree - divider.degree;
        MutablePolynomial q = remainderMonomial(dividend.clone().reverse().multiply(invRevMod.getInverse(m + 1), modulus), m + 1, false).reverse();
        if (q.degree < m)
            q.shiftRight(m - q.degree);
        return new MutablePolynomial[]{q, (copy ? dividend.clone() : dividend).subtract(divider.clone().multiply(q, modulus), modulus)};
    }

    /** fast division checks */
    private static MutablePolynomial earlyRemainderChecks(final MutablePolynomial dividend,
                                                          final MutablePolynomial divider,
                                                          final long modulus,
                                                          final boolean copy) {
        if (dividend.degree < divider.degree)
            return (copy ? dividend.clone() : dividend).modulus(modulus);
        if (divider.degree == 0)
            return MutablePolynomial.zero();
        if (divider.degree == 1)
            return MutablePolynomial.create(dividend.evaluate(multiplyMod(-divider.cc(), modInverse(divider.lc(), modulus), modulus), modulus));
        return null;
    }

    /**
     * Returns the remainder of {@code dividend} divided by {@code divider} modulo {@code modulus}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  prime modulus
     * @return the remainder
     */
    public static MutablePolynomial remainder(final MutablePolynomial dividend,
                                              final MutablePolynomial divider,
                                              final long modulus,
                                              final boolean copy) {
        MutablePolynomial rem = earlyRemainderChecks(dividend, divider, modulus, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, modulus, copy);

        return divideAndRemainderFast0(dividend, divider, modulus, copy)[1];
    }

    /** Plain school implementation */
    static MutablePolynomial remainderClassical0(final MutablePolynomial dividend,
                                                 final MutablePolynomial divider,
                                                 final long modulus,
                                                 final boolean copy) {
        assert dividend.degree >= divider.degree;

        MutablePolynomial remainder = copy ? dividend.clone() : dividend;
        for (int i = dividend.degree - divider.degree; i >= 0; --i)
            if (remainder.degree == divider.degree + i)
                remainder.subtract(divider, divideMod(remainder.lc(), divider.lc(), modulus), i, modulus);

        return remainder.modulus(modulus);
    }

    /**
     * Fast remainder using Newton iterations with switch to classical remainder for small polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   {@code divider^(-1) mod x^i }
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the remainder
     */
    public static MutablePolynomial remainderFastWithSwitch(final MutablePolynomial dividend,
                                                            final MutablePolynomial divider,
                                                            final InverseModMonomial invMod,
                                                            final long modulus,
                                                            final boolean copy) {
        MutablePolynomial rem = earlyRemainderChecks(dividend, divider, modulus, copy);
        if (rem != null)
            return rem;

        if (useClassicalDivision(dividend, divider))
            return remainderClassical0(dividend, divider, modulus, copy);

        return divideAndRemainderFast0(dividend, divider, invMod, modulus, copy)[1];
    }

    /**
     * Returns the quotient of {@code dividend} divided by {@code divider} modulo {@code modulus}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param modulus  prime modulus
     * @return the quotient
     */
    public static MutablePolynomial quotient(final MutablePolynomial dividend,
                                             final MutablePolynomial divider,
                                             final long modulus,
                                             final boolean copy) {
        MutablePolynomial[] qd = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (qd != null)
            return qd[0];

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic(dividend, divider, modulus, copy)[0];

        return divideAndRemainderFast0(dividend, divider, modulus, copy)[0];
    }

    /**
     * Fast quotient using Newton iterations with switch to classical remainder for small polynomials.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @param invMod   {@code divider^(-1) mod x^i }
     * @param modulus  the modulus
     * @param copy     whether to clone {@code dividend}; if not, the remainder will be placed directly to
     *                 {@code dividend} and {@code dividend} data will be lost
     * @return the quotient
     */
    public static MutablePolynomial quotientFastWithSwitch(final MutablePolynomial dividend,
                                                           final MutablePolynomial divider,
                                                           final InverseModMonomial invMod,
                                                           final long modulus,
                                                           final boolean copy) {
        MutablePolynomial[] qd = earlyDivideAndRemainderChecks(dividend, divider, modulus, copy);
        if (qd != null)
            return qd[0];

        if (useClassicalDivision(dividend, divider))
            return divideAndRemainderClassic(dividend, divider, modulus, copy)[0];

        return divideAndRemainderFast0(dividend, divider, invMod, modulus, copy)[0];
    }
}
