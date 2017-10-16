package cc.redberry.rings.poly.multivar;


import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.PolynomialRing;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;
import java.util.List;

/**
 * Multivariate interpolation
 *
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
        return new MultivariateInterpolation.Interpolation<>(variable, values[0])
                .update(points, values)
                .getInterpolatingPolynomial();
//        int length = points.length;
//
//        // Newton's representation
//        MultivariatePolynomial<E>[] mixedRadix = new MultivariatePolynomial[length];
//        mixedRadix[0] = values[0].clone();
//        MultivariatePolynomial<E> lins = values[0].createOne();
//        MultivariatePolynomial<E> poly = values[0].clone();
//
//        Ring<E> ring = poly.ring;
//        for (int k = 1; k < length; ++k) {
//            E reciprocal = ring.subtract(points[k], points[0]);
//            MultivariatePolynomial<E> accumulator = mixedRadix[0].clone();
//            for (int i = 1; i < k; ++i) {
//                accumulator = accumulator.add(mixedRadix[i].clone().multiply(reciprocal));
//                reciprocal = ring.multiply(reciprocal, ring.subtract(points[k], points[i]));
//            }
//            mixedRadix[k] = values[k].clone().subtract(accumulator).multiply(ring.reciprocal(reciprocal));
//
//            lins = lins.multiply(lins.createLinear(variable, ring.negate(points[k - 1]), ring.getOne()));
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
        private final List<E> points = new ArrayList<>();
        /** list of values at points */
        private final List<MultivariatePolynomial<E>> values = new ArrayList<>();
        /** mixed radix form of interpolating polynomial */
        private final List<MultivariatePolynomial<E>> mixedRadix = new ArrayList<>();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final MultivariatePolynomial<E> lins;
        /** resulting interpolating polynomial */
        private final MultivariatePolynomial<E> poly;
        /** ring */
        private final Ring<E> ring;

        /**
         * Start new interpolation with {@code interpolation[variable = point] = value}
         *
         * @param variable interpolating variable
         * @param point    evaluation point
         * @param value    polynomial value at {@code point}
         */
        public Interpolation(int variable, E point, MultivariatePolynomial<E> value) {
            this.variable = variable;
            this.lins = value.createOne();
            this.poly = value.clone();
            this.ring = poly.ring;

            points.add(point);
            values.add(value);
            mixedRadix.add(value.clone());
        }

        /**
         * Start new interpolation
         *
         * @param variable interpolating variable
         * @param factory  factory polynomial
         */
        public Interpolation(int variable, MultivariatePolynomial<E> factory) {
            this.variable = variable;
            this.lins = factory.createOne();
            this.poly = factory.createOne();
            this.ring = poly.ring;
        }

        /**
         * Start new interpolation
         *
         * @param variable interpolating variable
         * @param factory  factory polynomial
         */
        public Interpolation(int variable, PolynomialRing<MultivariatePolynomial<E>> factory) {
            this(variable, factory.factory());
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public Interpolation<E> update(E point, MultivariatePolynomial<E> value) {
            if (points.isEmpty()) {
                poly.multiply(value);
                points.add(point);
                values.add(value);
                mixedRadix.add(value.clone());
                return this;
            }
            E reciprocal = ring.subtract(point, points.get(0));
            MultivariatePolynomial<E> accumulator = mixedRadix.get(0).clone();
            for (int i = 1; i < points.size(); ++i) {
                accumulator = accumulator.add(mixedRadix.get(i).clone().multiply(reciprocal));
                reciprocal = ring.multiply(reciprocal, ring.subtract(point, points.get(i)));
            }
            if (ring.isZero(reciprocal))
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            mixedRadix.add(value.clone().subtract(accumulator).multiply(ring.reciprocal(reciprocal)));

            lins.multiply(lins.createLinear(variable, ring.negate(points.get(points.size() - 1)), ring.getOne()));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
            return this;
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param points evaluation points
         * @param values polynomial values at {@code point}
         */
        public Interpolation<E> update(E[] points, MultivariatePolynomial<E>[] values) {
            for (int i = 0; i < points.length; i++)
                update(points[i], values[i]);
            return this;
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

        /**
         * Returns the number of interpolation points used
         *
         * @return number of interpolation points used
         */
        public int numberOfPoints() {return points.size();}
    }

    /**
     * Updatable Newton interpolation
     */
    public static final class InterpolationZp64 {
        /** variable */
        private final int variable;
        /** list of evaluation points */
        private final TLongArrayList points = new TLongArrayList();
        /** list of values at points */
        private final List<MultivariatePolynomialZp64> values = new ArrayList<>();
        /** mixed radix form of interpolating polynomial */
        private final List<MultivariatePolynomialZp64> mixedRadix = new ArrayList<>();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final MultivariatePolynomialZp64 lins;
        /** resulting interpolating polynomial */
        private final MultivariatePolynomialZp64 poly;
        /** ring */
        private final IntegersZp64 ring;

        /**
         * Start new interpolation with {@code interpolation[variable = point] = value}
         *
         * @param variable interpolating variable
         * @param point    evaluation point
         * @param value    polynomial value at {@code point}
         */
        public InterpolationZp64(int variable, long point, MultivariatePolynomialZp64 value) {
            this.variable = variable;
            this.lins = value.createOne();
            this.poly = value.clone();
            this.ring = poly.ring;

            points.add(point);
            values.add(value);
            mixedRadix.add(value.clone());
        }

        /**
         * Start new interpolation
         *
         * @param variable interpolating variable
         * @param factory  factory polynomial
         */
        public InterpolationZp64(int variable, MultivariatePolynomialZp64 factory) {
            this.variable = variable;
            this.lins = factory.createOne();
            this.poly = factory.createOne();
            this.ring = poly.ring;
        }

        /**
         * Start new interpolation
         *
         * @param variable interpolating variable
         * @param factory  factory polynomial
         */
        public InterpolationZp64(int variable, PolynomialRing<MultivariatePolynomialZp64> factory) {
            this(variable, factory.factory());
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public InterpolationZp64 update(long point, MultivariatePolynomialZp64 value) {
            if (points.isEmpty()) {
                poly.multiply(value);
                points.add(point);
                values.add(value);
                mixedRadix.add(value.clone());
                return this;
            }
            long reciprocal = ring.subtract(point, points.get(0));
            MultivariatePolynomialZp64 accumulator = mixedRadix.get(0).clone();
            for (int i = 1; i < points.size(); ++i) {
                accumulator = accumulator.add(mixedRadix.get(i).clone().multiply(reciprocal));
                reciprocal = ring.multiply(reciprocal, ring.subtract(point, points.get(i)));
            }
            if (reciprocal == 0)
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            mixedRadix.add(value.clone().subtract(accumulator).multiply(ring.reciprocal(reciprocal)));

            lins.multiply(lins.createLinear(variable, ring.negate(points.get(points.size() - 1)), 1));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
            return this;
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param points evaluation points
         * @param values polynomial values at {@code point}
         */
        public InterpolationZp64 update(long points[], MultivariatePolynomialZp64 values[]) {
            for (int i = 0; i < points.length; i++)
                update(points[i], values[i]);
            return this;
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
        public MultivariatePolynomialZp64 getInterpolatingPolynomial() {return poly;}

        /**
         * Returns the list of evaluation points used in interpolation
         *
         * @return list of evaluation points used in interpolation
         */
        public TLongArrayList getPoints() {return points;}

        /**
         * Returns the list of polynomial values at interpolation points
         *
         * @return the list of polynomial values at interpolation points
         */
        public List<MultivariatePolynomialZp64> getValues() {return values;}

        /**
         * Returns the number of interpolation points used
         *
         * @return number of interpolation points used
         */
        public int numberOfPoints() {return points.size();}
    }
}