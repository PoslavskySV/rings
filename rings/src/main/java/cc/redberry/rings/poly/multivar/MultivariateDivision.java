package cc.redberry.rings.poly.multivar;

import java.util.Collection;

/**
 * Division with remainder of multivariate polynomials (multivariate reduction).
 *
 * @since 1.0
 */
public final class MultivariateDivision {
    private MultivariateDivision() {}

    /**
     * Performs multivariate division with remainder. The resulting array of quotients and remainder (last element of
     * the returned array) satisfies {@code dividend = quotient_1 * divider_1 + quotient_2 * divider_2 + ... + remainder
     * }.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder in the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] divideAndRemainder(Poly dividend, Poly... dividers) {
        Poly[] quotients = dividend.createArray(dividers.length + 1);
        int i = 0;
        int constDivider = -1;
        for (; i < dividers.length; i++) {
            if (dividers[i].isZero())
                throw new ArithmeticException("divide by zero");
            if (dividers[i].isConstant())
                constDivider = i;
            quotients[i] = dividend.createZero();
        }
        quotients[i] = dividend.createZero();

        if (constDivider != -1) {
            if (dividers[constDivider].isOne()) {
                quotients[constDivider] = dividend.clone();
                return quotients;
            }
            Poly dd = dividend.clone().divideByLC(dividers[constDivider]);
            if (dd != null) {
                quotients[constDivider] = dd;
                return quotients;
            }
        }

        dividend = dividend.clone();
        Poly remainder = quotients[quotients.length - 1];
        while (!dividend.isZero()) {
            Term ltDiv = null;
            for (i = 0; i < dividers.length; i++) {
                ltDiv = dividend.divideOrNull(dividend.lt(), dividers[i].lt());
                if (ltDiv != null)
                    break;
            }
            if (ltDiv != null) {
                quotients[i] = quotients[i].add(ltDiv);
                dividend = dividend.subtract(dividend.create(ltDiv).multiply(dividers[i]));
            } else {
                remainder = remainder.add(dividend.lt());
                dividend = dividend.subtractLt();
            }
        }
        return quotients;
    }

    /**
     * Performs multivariate division with remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return array of quotient and remainder
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        Poly[] array = divider.createArray(1);
        array[0] = divider;
        return divideAndRemainder(dividend, array);
    }

    /**
     * Performs multivariate division with remainder and returns the remainder.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Poly... dividers) {
        Poly[] r = divideAndRemainder(dividend, dividers);
        return r[r.length - 1];
    }

    /**
     * Performs multivariate division with remainder and rerurns the remainder.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Collection<Poly> dividers) {
        Poly[] r = divideAndRemainder(dividend, dividers.toArray(dividend.createArray(dividers.size())));
        return r[r.length - 1];
    }

    /**
     * Performs multivariate division with remainder and rerurns the remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Poly divider) {
        Poly[] r = divideAndRemainder(dividend, divider);
        return r[r.length - 1];
    }

    /**
     * Divides {@code dividend} by {@code divider} or throws exception if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider}
     * @throws ArithmeticException if exact division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly divideExact(Poly dividend, Poly divider) {
        Poly[] qd = divideAndRemainder(dividend, divider);
        if (qd == null || !qd[1].isZero())
            throw new ArithmeticException("not divisible: " + dividend + " / " + divider);
        return qd[0];
    }

    /**
     * Divides {@code dividend} by {@code divider} or returns null if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider} or null if exact division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly divideOrNull(Poly dividend, Poly divider) {
        Poly[] qd = divideAndRemainder(dividend, divider);
        if (qd == null || !qd[1].isZero())
            return null;
        return qd[0];
    }

    /**
     * Tests whether {@code divisor} is a divisor of {@code poly}
     *
     * @param dividend the polynomial
     * @param divider  the divisor to check
     * @return whether {@code divisor} is a divisor of {@code poly}
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean dividesQ(Poly dividend, Poly divider) {
        if (divider.isOne())
            return true;
        dividend = dividend.clone();
        if (divider.isConstant())
            return dividend.divideByLC(divider) != null;
        int[]
                dividendDegrees = dividend.degrees(),
                dividerDegrees = divider.degrees();
        for (int i = 0; i < dividendDegrees.length; i++)
            if (dividendDegrees[i] < dividerDegrees[i])
                return false;

        while (!dividend.isZero()) {
            Term ltDiv = dividend.divideOrNull(dividend.lt(), divider.lt());
            if (ltDiv == null)
                return false;
            dividend = dividend.subtract(divider.clone().multiply(ltDiv));
        }
        return true;
    }
}
