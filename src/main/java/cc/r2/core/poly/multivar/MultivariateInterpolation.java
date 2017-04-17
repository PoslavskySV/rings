package cc.r2.core.poly.multivar;

import cc.r2.core.poly.generics.Domain;

import java.util.ArrayList;
import java.util.List;

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

//    /**
//     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
//     * method uses Newton's mixed radix iterations.
//     *
//     * @param points evaluation points
//     * @param values corresponding polynomial values
//     * @return the interpolating polynomial
//     */
//    @SuppressWarnings("unchecked")
//    public static <E> Interpolation<E> interpolation(int variable, E[] points, MultivariatePolynomial<E>[] values) {
//        checkInput(points, values);
//        int length = points.length;
//
//        // Newton's representation
//        MultivariatePolynomial<E>[] mixedRadix = new MultivariatePolynomial[length];
//        mixedRadix[0] = values[0].clone();
//        MultivariatePolynomial<E> lins = values[0].createOne();
//        MultivariatePolynomial<E> poly = mixedRadix[0].clone();
//
//        Domain<E> domain = poly.domain;
//        for (int k = 1; k < length; ++k) {
//            E reciprocal = domain.subtract(points[k], points[0]);
//            MultivariatePolynomial<E> accumulator = mixedRadix[0];
//            for (int i = 1; i < k; ++i) {
//                accumulator = accumulator.add(mixedRadix[i].clone().multiply(reciprocal));
//                reciprocal = domain.multiply(reciprocal, domain.subtract(points[k], points[i]));
//            }
//            mixedRadix[k] = values[k].clone().subtract(accumulator).multiply(domain.reciprocal(reciprocal));
//
//            lins = lins.multiply(lins.createLinear(variable, domain.negate(points[k - 1]), domain.getOne()));
//            poly = poly.add(lins.clone().multiply(mixedRadix[k]));
//        }
//        return new Interpolation<>(variable, points, values, mixedRadix, lins, poly);
//    }

    static final class Interpolation<E> {
        final int variable;
        final List<E> points;
        final List<MultivariatePolynomial<E>> values;
        final List<MultivariatePolynomial<E>> mixedRadix;
        final MultivariatePolynomial<E> lins;
        final MultivariatePolynomial<E> poly;

        public Interpolation(int variable, E point, MultivariatePolynomial<E> value) {
            this.variable = variable;
            this.points = new ArrayList<>();
            this.values = new ArrayList<>();
            this.mixedRadix = new ArrayList<>();
            this.lins = value.createOne();
            this.poly = value.clone();

            points.add(point);
            values.add(value);
            mixedRadix.add(value);
        }

        public void update(E point, MultivariatePolynomial<E> value) {
            Domain<E> domain = poly.domain;
            E reciprocal = domain.subtract(point, points.get(0));
            MultivariatePolynomial<E> accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = accumulator.add(mixedRadix.get(i).clone().multiply(reciprocal));
                reciprocal = domain.multiply(reciprocal, domain.subtract(point, points.get(i)));
            }
            mixedRadix.add(value.clone().subtract(accumulator).multiply(domain.reciprocal(reciprocal)));

            lins.multiply(lins.createLinear(variable, domain.negate(points.get(points.size() - 1)), domain.getOne()));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
        }
    }

    public static <E> MultivariatePolynomial<E> interpolateNewton(int variable, E[] points,
                                                                  MultivariatePolynomial<E>[] values,
                                                                  MultivariatePolynomial<E> result) {
        checkInput(points, values);
        int length = points.length;

        // Newton's representation
        MultivariatePolynomial<E>[] mixedRadix = new MultivariatePolynomial[length];
        mixedRadix[0] = values[0].clone();
        MultivariatePolynomial<E> lins = values[0].createOne();

        Domain<E> domain = result.domain;
        for (int k = 1; k < length; ++k) {
            E reciprocal = domain.subtract(points[k], points[0]);
            MultivariatePolynomial<E> accumulator = mixedRadix[0];
            for (int i = 1; i < k; ++i) {
                accumulator = accumulator.add(mixedRadix[i].clone().multiply(reciprocal));
                reciprocal = domain.multiply(reciprocal, domain.subtract(points[k], points[i]));
            }
            mixedRadix[k] = values[k].clone().subtract(accumulator).multiply(domain.reciprocal(reciprocal));

            lins = lins.multiply(lins.createLinear(variable, domain.negate(points[k - 1]), domain.getOne()));
            result = result.add(lins.clone().multiply(mixedRadix[k]));
        }
        return result;
    }
}
