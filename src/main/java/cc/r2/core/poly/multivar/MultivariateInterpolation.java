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
        Interpolation<E> interpolation = new Interpolation<>(variable, points[0], values[0]);
        for (int i = 1; i < points.length; i++)
            interpolation.update(points[i], values[i]);
        return interpolation.getInterpolatingPolynomial();
//        int length = points.length;
//
//        // Newton's representation
//        MultivariatePolynomial<E>[] mixedRadix = new MultivariatePolynomial[length];
//        mixedRadix[0] = values[0].clone();
//        MultivariatePolynomial<E> lins = values[0].createOne();
//        MultivariatePolynomial<E> poly = values[0].clone();
//
//        Domain<E> domain = poly.domain;
//        for (int k = 1; k < length; ++k) {
//            E reciprocal = domain.subtract(points[k], points[0]);
//            MultivariatePolynomial<E> accumulator = mixedRadix[0].clone();
//            for (int i = 1; i < k; ++i) {
//                accumulator = accumulator.add(mixedRadix[i].clone().multiply(reciprocal));
//                reciprocal = domain.multiply(reciprocal, domain.subtract(points[k], points[i]));
//            }
//            mixedRadix[k] = values[k].clone().subtract(accumulator).multiply(domain.reciprocal(reciprocal));
//
//            lins = lins.multiply(lins.createLinear(variable, domain.negate(points[k - 1]), domain.getOne()));
//            poly = poly.add(lins.clone().multiply(mixedRadix[k]));
//        }
//        return poly;
    }

    /**
     * Updatable Newton interpolation
     */
    public static final class Interpolation<E> {
        /** variable */
        private final int variable;
        /** list of evaluation points */
        private final List<E> points;
        /** list of values at points */
        private final List<MultivariatePolynomial<E>> values;
        /** mixed radix form of interpolating polynomial */
        private final List<MultivariatePolynomial<E>> mixedRadix;
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final MultivariatePolynomial<E> lins;
        /** resulting interpolating polynomial */
        private final MultivariatePolynomial<E> poly;
        /** domain */
        private final Domain<E> domain;

        /**
         * Start new interpolation with {@code interpolation[variable = point] = value}
         *
         * @param variable interpolating variable
         * @param point    evaluation point
         * @param value    polynomial value at {@code point}
         */
        public Interpolation(int variable, E point, MultivariatePolynomial<E> value) {
            this.variable = variable;
            this.points = new ArrayList<>();
            this.values = new ArrayList<>();
            this.mixedRadix = new ArrayList<>();
            this.lins = value.createOne();
            this.poly = value.clone();
            this.domain = poly.domain;

            points.add(point);
            values.add(value);
            mixedRadix.add(value.clone());
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public void update(E point, MultivariatePolynomial<E> value) {
            E reciprocal = domain.subtract(point, points.get(0));
            MultivariatePolynomial<E> accumulator = mixedRadix.get(0).clone();
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

        /**
         * Returns variable used in the interpolation
         *
         * @return variable used in the interpolation
         */
        public int getVariable() {return variable;}

        /**
         * Returns resulting interpolating polynomial
         *
         * @return interpolating polynomial
         */
        public MultivariatePolynomial<E> getInterpolatingPolynomial() {return poly;}

        /**
         * Returns the list of evaluation points used in interpolation
         *
         * @return list of evaluation points used in interpolation
         */
        public List<E> getPoints() {return points;}

        /**
         * Returns the list of polynomial values at interpolation points
         *
         * @return the list of polynomial values at interpolation points
         */
        public List<MultivariatePolynomial<E>> getValues() {return values;}
    }
}
