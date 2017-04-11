package cc.r2.core.poly.multivar;

import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.multivar.MultivariatePolynomial.DegreeVector;

import java.util.Map.Entry;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class DivisionWithRemainderMultivariate {
    private DivisionWithRemainderMultivariate() {}

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
            Entry<DegreeVector, E> ltDiv = null;
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

    /** Monomial division */
    private static <E> Entry<DegreeVector, E> divide(Domain<E> domain,
                                                     Entry<DegreeVector, E> ltDividend,
                                                     Entry<DegreeVector, E> ltDivider) {
        E[] qr = domain.divideAndRemainder(ltDividend.getValue(), ltDivider.getValue());
        if (!domain.isZero(qr[1]))
            return null;
        int[] dividendV = ltDividend.getKey().exponents, dividerV = ltDivider.getKey().exponents;
        int[] exponents = new int[dividendV.length];
        for (int i = 0; i < exponents.length; i++) {
            exponents[i] = dividendV[i] - dividerV[i];
            if (exponents[i] < 0)
                return null;
        }
        return new MultivariatePolynomial.EntryImpl<>(new DegreeVector(exponents), qr[0]);
    }

}
