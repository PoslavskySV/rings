package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

/**
 * Univariate polynomial interpolation.
 *
 * @since 1.0
 */
public final class UnivariateInterpolation {
    private UnivariateInterpolation() {}

    private static void checkInput(long[] points, long[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Lagrange's interpolation formula.
     *
     * @param modulus the modulus
     * @param points  evaluation points
     * @param values  corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static UnivariatePolynomialZp64 interpolateLagrange(long modulus, long[] points, long[] values) {
        checkInput(points, values);

        int length = points.length;
        IntegersZp64 ring = new IntegersZp64(modulus);
        UnivariatePolynomialZp64 result = UnivariatePolynomialZp64.zero(ring);
        for (int i = 0; i < length; ++i) {
            UnivariatePolynomialZp64 interpolant = UnivariatePolynomialZp64.constant(modulus, values[i]);
            for (int j = 0; j < length; ++j) {
                if (j == i)
                    continue;
                UnivariatePolynomialZp64 linear = result
                        .createLinear(ring.negate(points[j]), 1)
                        .divide(ring.subtract(points[i], points[j]));
                interpolant = interpolant.multiply(linear);
            }
            result = result.add(interpolant);
        }
        return result;
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Lagrange's interpolation formula.
     *
     * @param ring   the ring
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static <E> UnivariatePolynomial<E> interpolateLagrange(Ring<E> ring, E[] points, E[] values) {
        checkInput(points, values);

        int length = points.length;
        UnivariatePolynomial<E> result = UnivariatePolynomial.zero(ring);
        for (int i = 0; i < length; ++i) {
            UnivariatePolynomial<E> interpolant = UnivariatePolynomial.constant(ring, values[i]);
            for (int j = 0; j < length; ++j) {
                if (j == i)
                    continue;
                UnivariatePolynomial<E> linear = result
                        .createLinear(ring.negate(points[j]), ring.getOne())
                        .divideExact(ring.subtract(points[i], points[j]));
                interpolant = interpolant.multiply(linear);
            }
            result = result.add(interpolant);
        }
        return result;
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Newton's mixed radix iterations.
     *
     * @param modulus the modulus
     * @param points  evaluation points
     * @param values  corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static UnivariatePolynomialZp64 interpolateNewton(long modulus, long[] points, long[] values) {
        return interpolateNewton(new IntegersZp64(modulus), points, values);
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Newton's mixed radix iterations.
     *
     * @param ring   the ring
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static UnivariatePolynomialZp64 interpolateNewton(IntegersZp64 ring, long[] points, long[] values) {
        checkInput(points, values);
        return new InterpolationZp64(ring)
                .update(points, values)
                .getInterpolatingPolynomial();
    }

    private static void checkInput(Object[] points, Object[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Newton's mixed radix iterations.
     *
     * @param ring   the ring
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static <E> UnivariatePolynomial<E> interpolateNewton(Ring<E> ring, E[] points, E[] values) {
        checkInput(points, values);
        return new Interpolation<>(ring)
                .update(points, values)
                .getInterpolatingPolynomial();
    }

    /**
     * Updatable Newton interpolation
     */
    public static final class InterpolationZp64 {
        private final IntegersZp64 ring;
        /** list of evaluation points */
        private final TLongArrayList points = new TLongArrayList();
        /** list of values at points */
        private final TLongArrayList values = new TLongArrayList();
        /** mixed radix form of interpolating polynomial */
        private final TLongArrayList mixedRadix = new TLongArrayList();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final UnivariatePolynomialZp64 lins;
        /** resulting interpolating polynomial */
        private final UnivariatePolynomialZp64 poly;

        /**
         * Start new interpolation with {@code interpolation[point] = value}
         *
         * @param ring the ring
         */
        public InterpolationZp64(IntegersZp64 ring) {
            this.ring = ring;
            this.lins = UnivariatePolynomialZp64.one(ring);
            this.poly = UnivariatePolynomialZp64.one(ring);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public InterpolationZp64 update(long point, long value) {
            if (points.isEmpty()) {
                poly.multiply(value);
                points.add(point);
                values.add(value);
                mixedRadix.add(value);
                return this;
            }
            long reciprocal = poly.subtract(point, points.get(0));
            long accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = ring.add(accumulator, ring.multiply(mixedRadix.get(i), reciprocal));
                reciprocal = ring.multiply(reciprocal, ring.subtract(point, points.get(i)));
            }
            if (reciprocal == 0)
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = ring.reciprocal(reciprocal);
            mixedRadix.add(ring.multiply(reciprocal, ring.subtract(value, accumulator)));

            lins.multiply(lins.createLinear(ring.negate(points.get(points.size() - 1)), 1));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
            return this;
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param points evaluation points
         * @param values polynomial values at {@code points}
         */
        public InterpolationZp64 update(long[] points, long[] values) {
            for (int i = 0; i < points.length; i++)
                update(points[i], values[i]);
            return this;
        }

        /**
         * Returns resulting interpolating polynomial
         *
         * @return interpolating polynomial
         */
        public UnivariatePolynomialZp64 getInterpolatingPolynomial() {return poly;}

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
        public TLongArrayList getValues() {return values;}

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
    public static final class Interpolation<E> {
        private final Ring<E> ring;
        /** list of evaluation points */
        private final ArrayList<E> points = new ArrayList<>();
        /** list of values at points */
        private final ArrayList<E> values = new ArrayList<>();
        /** mixed radix form of interpolating polynomial */
        private final ArrayList<E> mixedRadix = new ArrayList<>();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final UnivariatePolynomial<E> lins;
        /** resulting interpolating polynomial */
        private final UnivariatePolynomial<E> poly;

        /**
         * Start new interpolation with {@code interpolation[point] = value}
         *
         * @param ring the ring
         */
        public Interpolation(Ring<E> ring) {
            this.ring = ring;
            this.lins = UnivariatePolynomial.one(ring);
            this.poly = UnivariatePolynomial.one(ring);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public Interpolation<E> update(E point, E value) {
            if (points.isEmpty()) {
                points.add(point);
                values.add(value);
                mixedRadix.add(value);
                poly.multiply(value);
                return this;
            }
            E reciprocal = ring.subtract(point, points.get(0));
            E accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = ring.add(accumulator, ring.multiply(mixedRadix.get(i), reciprocal));
                reciprocal = ring.multiply(reciprocal, ring.subtract(point, points.get(i)));
            }
            if (ring.isZero(reciprocal))
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = ring.reciprocal(reciprocal);
            mixedRadix.add(ring.multiply(reciprocal, ring.subtract(value, accumulator)));

            lins.multiply(lins.createLinear(ring.negate(points.get(points.size() - 1)), ring.getOne()));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
            return this;
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param points evaluation points
         * @param values polynomial values at {@code points}
         */
        public Interpolation<E> update(E[] points, E[] values) {
            for (int i = 0; i < points.length; i++)
                update(points[i], values[i]);
            return this;
        }

        /**
         * Returns resulting interpolating polynomial
         *
         * @return interpolating polynomial
         */
        public UnivariatePolynomial<E> getInterpolatingPolynomial() {return poly;}

        /**
         * Returns the list of evaluation points used in interpolation
         *
         * @return list of evaluation points used in interpolation
         */
        public ArrayList<E> getPoints() {return points;}

        /**
         * Returns the list of polynomial values at interpolation points
         *
         * @return the list of polynomial values at interpolation points
         */
        public ArrayList<E> getValues() {return values;}

        /**
         * Returns the number of interpolation points used
         *
         * @return number of interpolation points used
         */
        public int numberOfPoints() {return points.size();}
    }
}