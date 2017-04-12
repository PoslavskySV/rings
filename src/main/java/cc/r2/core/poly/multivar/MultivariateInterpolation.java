package cc.r2.core.poly.multivar;

import cc.r2.core.poly.generics.Domain;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateInterpolation {
    private MultivariateInterpolation() {}

    private static void checkInput(Object[] points, Object[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Newton's mixed radix iterations.
     *
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    @SuppressWarnings("unchecked")
    public static <E> MultivariatePolynomial<E> interpolateNewton(int variable, E[] points, MultivariatePolynomial<E>[] values) {
        checkInput(points, values);
        int length = points.length;

        // Newton's representation
        MultivariatePolynomial<E>[] mixedRadix = new MultivariatePolynomial[length];
        mixedRadix[0] = values[0].clone();
        MultivariatePolynomial<E> lins = values[0].createOne();
        MultivariatePolynomial<E> poly = mixedRadix[0].clone();

        Domain<E> domain = poly.domain;
        for (int k = 1; k < length; ++k) {
            E reciprocal = domain.subtract(points[k], points[0]);
            MultivariatePolynomial<E> accumulator = mixedRadix[0];
            for (int i = 1; i < k; ++i) {
                accumulator = accumulator.add(mixedRadix[i].clone().multiply(reciprocal));
                reciprocal = domain.multiply(reciprocal, domain.subtract(points[k], points[i]));
            }
            mixedRadix[k] = values[k].clone().subtract(accumulator).multiply(domain.reciprocal(reciprocal));

            lins = lins.multiply(lins.createLinear(variable, domain.negate(points[k - 1]), domain.getOne()));
            poly = poly.add(lins.clone().multiply(mixedRadix[k]));
        }
        return poly;
    }
}
