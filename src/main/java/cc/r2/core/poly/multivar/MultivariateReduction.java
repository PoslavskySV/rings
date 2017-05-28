package cc.r2.core.poly.multivar;

import java.util.Collection;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateReduction {
    private MultivariateReduction() {}

    /**
     * Performs multivariate division with remainder. The resulting array of quotients and remainder (last element of
     * the returned array) satisfies {@code dividend = quotient_1 * divider_1 + quotient_2 * divider_2 + ... + remainder }.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] divideAndRemainder(Poly dividend, Poly...  dividers) {
        Poly[] quotients = dividend.arrayNewInstance(dividers.length + 1);
        int i = 0;
        for (; i < dividers.length; i++) {
            if (dividers[i].isZero())
                throw new ArithmeticException("divide by zero");
            quotients[i] = dividend.createZero();
        }
        quotients[i] = dividend.createZero();

        Poly remainder = quotients[quotients.length - 1];
        dividend = dividend.clone();
        while (!dividend.isZero()) {
            Term ltDiv = null;
            for (i = 0; i < dividers.length; i++) {
                ltDiv = dividend.divideOrNull(dividend.lt(), dividers[i].lt());
                if (ltDiv != null)
                    break;
            }
            if (ltDiv != null) {
                quotients[i] = quotients[i].add(ltDiv);
                dividend = dividend.subtract(dividers[i].clone().multiply(ltDiv));
            } else {
                remainder = remainder.add(dividend.lt());
                dividend = dividend.subtractLt();
            }
        }
        return quotients;
    }

    /**
     * Performs multivariate division with remainder and rerurns the remainder.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Poly[] dividers) {
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
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Collection<Poly> dividers) {
        Poly[] r = divideAndRemainder(dividend, dividers.toArray(dividend.arrayNewInstance(dividers.size())));
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
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
            throw new ArithmeticException("not divisible");
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
        dividend = dividend.clone();
        while (!dividend.isZero()) {
            Term ltDiv = dividend.divideOrNull(dividend.lt(), divider.lt());
            if (ltDiv == null)
                return false;
            dividend = dividend.subtract(divider.clone().multiply(ltDiv));
        }
        return true;
    }
}
