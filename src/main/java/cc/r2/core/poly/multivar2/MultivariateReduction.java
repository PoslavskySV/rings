package cc.r2.core.poly.multivar2;

import cc.r2.core.poly.generics.Domain;

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
    public static <E> MultivariatePolynomial<E>[] divideAndRemainder(MultivariatePolynomial<E> dividend, MultivariatePolynomial<E>... dividers) {
        MultivariatePolynomial<E>[] quotients = new MultivariatePolynomial[dividers.length + 1];
        int i = 0;
        for (; i < dividers.length; i++) {
            if (dividers[i].isZero())
                throw new ArithmeticException("divide by zero");
            quotients[i] = dividend.createZero();
        }
        quotients[i] = dividend.createZero();

        MultivariatePolynomial<E> remainder = quotients[quotients.length - 1];
        dividend = dividend.clone();
        while (!dividend.isZero()) {
            MonomialTerm<E> ltDiv = null;
            for (i = 0; i < dividers.length; i++) {
                ltDiv = divide(dividend.domain, dividend.lt(), dividers[i].lt());
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
     * Divides {@code dividend} by {@code divider} or throws exception if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider}
     * @throws ArithmeticException if exact division is not possible
     */
    @SuppressWarnings("unchecked")
    public static <E> MultivariatePolynomial<E> divideExact(MultivariatePolynomial<E> dividend, MultivariatePolynomial<E> divider) {
        MultivariatePolynomial<E>[] qd = divideAndRemainder(dividend, divider);
        if (qd == null || !qd[1].isZero())
            throw new ArithmeticException("not divisible");
        return qd[0];
    }

    /**
     * Tests whether {@code divisor} is a divisor of {@code poly}
     *
     * @param poly    the polynomial
     * @param divisor the divisor to check
     * @return whether {@code divisor} is a divisor of {@code poly}
     */
    @SuppressWarnings("unchecked")
    public static <E> boolean dividesQ(MultivariatePolynomial<E> poly, MultivariatePolynomial<E> divisor) {
        poly = poly.clone();
        while (!poly.isZero()) {
            MonomialTerm<E> ltDiv = divide(poly.domain, poly.lt(), divisor.lt());
            if (ltDiv == null)
                return false;
            poly = poly.subtract(divisor.clone().multiply(ltDiv));
        }
        return true;
    }

    /** Monomial division */
    private static <E> MonomialTerm<E> divide(Domain<E> domain,
                                              MonomialTerm<E> ltDividend,
                                              MonomialTerm<E> ltDivider) {
        E[] qr = domain.divideAndRemainder(ltDividend.coefficient, ltDivider.coefficient);
        if (!domain.isZero(qr[1]))
            return null;
        return ltDividend.divide(ltDivider, qr[0]);
    }
}
