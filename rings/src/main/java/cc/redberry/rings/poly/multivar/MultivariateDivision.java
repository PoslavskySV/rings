package cc.redberry.rings.poly.multivar;

import java.util.Arrays;
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

        IMonomialAlgebra<Term> mAlgebra = dividend.monomialAlgebra;
        // cache leading terms
        Term[] dividersLTs = Arrays.stream(dividers).map(Poly::lt).toArray(mAlgebra::createArray);
        dividend = dividend.clone();
        Poly remainder = quotients[quotients.length - 1];
        while (!dividend.isZero()) {
            Term ltDiv = null;
            Term lt = dividend.lt();
            for (i = 0; i < dividers.length; i++) {
                ltDiv = mAlgebra.divideOrNull(lt, dividersLTs[i]);
                if (ltDiv != null)
                    break;
            }
            if (ltDiv != null) {
                quotients[i] = quotients[i].add(ltDiv);
                dividend = dividend.subtract(dividend.create(ltDiv).multiply(dividers[i]));
            } else {
                remainder = remainder.add(lt);
                dividend = dividend.subtractLt();
            }
        }
        return quotients;
    }

    /**
     * Performs multivariate division with remainder and returns the remainder.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return the remainder
     */
    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Poly... dividers) {
        int i = 0;
        int constDivider = -1;
        for (; i < dividers.length; i++) {
            if (dividers[i].isZero())
                throw new ArithmeticException("divide by zero");
            if (dividers[i].isConstant())
                constDivider = i;
        }

        if (constDivider != -1) {
            if (dividers[constDivider].isOne())
                return dividend.createZero();

            Poly dd = dividend.clone().divideByLC(dividers[constDivider]);
            if (dd != null)
                return dividend.createZero();
        }

        IMonomialAlgebra<Term> mAlgebra = dividend.monomialAlgebra;
        // cache leading terms
        Term[] dividersLTs = Arrays.stream(dividers).map(Poly::lt).toArray(mAlgebra::createArray);
        dividend = dividend.clone();
        Poly remainder = dividend.createZero();
        while (!dividend.isZero()) {
            Term ltDiv = null;
            Term lt = dividend.lt();
            for (i = 0; i < dividersLTs.length; ++i) {
                ltDiv = mAlgebra.divideOrNull(lt, dividersLTs[i]);
                if (ltDiv != null)
                    break;
            }
            if (ltDiv != null)
                dividend = dividend.subtract(dividend.create(ltDiv).multiply(dividers[i]));
            else {
                remainder = remainder.add(lt);
                dividend = dividend.subtractLt();
            }
        }
        return remainder;
    }

    /**
     * Performs multivariate division with remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return array of quotient and remainder
     */
    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        Poly[] array = divider.createArray(1);
        array[0] = divider;
        return divideAndRemainder(dividend, array);
    }

    /**
     * Performs multivariate division with remainder and rerurns the remainder.
     *
     * @param dividend the dividend
     * @param dividers the dividers
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Collection<Poly> dividers) {
        return remainder(dividend, dividers.toArray(dividend.createArray(dividers.size())));
    }

    /**
     * Performs multivariate division with remainder and rerurns the remainder.
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return array of quotients and remainder at the last position
     */
    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder(Poly dividend, Poly divider) {
        Poly[] array = divider.createArray(1);
        array[0] = divider;
        return remainder(dividend, array);
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

        IMonomialAlgebra<Term> mAlgebra = dividend.monomialAlgebra;
        while (!dividend.isZero()) {
            Term ltDiv = mAlgebra.divideOrNull(dividend.lt(), divider.lt());
            if (ltDiv == null)
                return false;
            dividend = dividend.subtract(divider.clone().multiply(ltDiv));
        }
        return true;
    }

    /**
     * Tests whether there is nontrivial quotient {@code dividend / divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return whether {@code divisor} is a divisor of {@code poly}
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean nontrivialQuotientQ(Poly dividend, Poly divider) {
        Term lt = divider.lt();
        for (Term term : dividend)
            if (term.dvDivisibleBy(lt))
                return true;
        return false;
    }
}
