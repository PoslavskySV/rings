package cc.r2.core.poly.univar;

import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersZp64;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

/**
 * Univariate polynomial interpolation.
 *
 * @author Stanislav Poslavsky
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
        IntegersZp64 domain = new IntegersZp64(modulus);
        UnivariatePolynomialZp64 result = UnivariatePolynomialZp64.zero(domain);
        for (int i = 0; i < length; ++i) {
            UnivariatePolynomialZp64 interpolant = UnivariatePolynomialZp64.constant(modulus, values[i]);
            for (int j = 0; j < length; ++j) {
                if (j == i)
                    continue;
                UnivariatePolynomialZp64 linear = result
                        .createLinear(domain.negate(points[j]), 1)
                        .divide(domain.subtract(points[i], points[j]));
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
     * @param domain the domain
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static <E> UnivariatePolynomial<E> interpolateLagrange(Domain<E> domain, E[] points, E[] values) {
        checkInput(points, values);

        int length = points.length;
        UnivariatePolynomial<E> result = UnivariatePolynomial.zero(domain);
        for (int i = 0; i < length; ++i) {
            UnivariatePolynomial<E> interpolant = UnivariatePolynomial.constant(domain, values[i]);
            for (int j = 0; j < length; ++j) {
                if (j == i)
                    continue;
                UnivariatePolynomial<E> linear = result
                        .createLinear(domain.negate(points[j]), domain.getOne())
                        .divideExact(domain.subtract(points[i], points[j]));
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
     * @param domain the domain
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static UnivariatePolynomialZp64 interpolateNewton(IntegersZp64 domain, long[] points, long[] values) {
        checkInput(points, values);
        InterpolationZp64 interpolation = new InterpolationZp64(domain);
        for (int i = 0; i < points.length; i++)
            interpolation.update(points[i], values[i]);
        return interpolation.getInterpolatingPolynomial();
    }

    private static void checkInput(Object[] points, Object[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Newton's mixed radix iterations.
     *
     * @param domain the domain
     * @param points evaluation points
     * @param values corresponding polynomial values
     * @return the interpolating polynomial
     */
    public static <E> UnivariatePolynomial<E> interpolateNewton(Domain<E> domain, E[] points, E[] values) {
        checkInput(points, values);
        Interpolation<E> interpolation = new Interpolation<>(domain);
        for (int i = 0; i < points.length; i++)
            interpolation.update(points[i], values[i]);
        return interpolation.getInterpolatingPolynomial();
    }

    /**
     * Updatable Newton interpolation
     */
    public static final class InterpolationZp64 {
        private final IntegersZp64 domain;
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
         * @param domain the domain
         */
        public InterpolationZp64(IntegersZp64 domain) {
            this.domain = domain;
            this.lins = UnivariatePolynomialZp64.one(domain);
            this.poly = UnivariatePolynomialZp64.one(domain);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public void update(long point, long value) {
            if (points.isEmpty()) {
                poly.multiply(value);
                points.add(point);
                values.add(value);
                mixedRadix.add(value);
                return;
            }
            long reciprocal = poly.subtract(point, points.get(0));
            long accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = domain.add(accumulator, domain.multiply(mixedRadix.get(i), reciprocal));
                reciprocal = domain.multiply(reciprocal, domain.subtract(point, points.get(i)));
            }
            if (reciprocal == 0)
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = domain.reciprocal(reciprocal);
            mixedRadix.add(domain.multiply(reciprocal, domain.subtract(value, accumulator)));

            lins.multiply(lins.createLinear(domain.negate(points.get(points.size() - 1)), 1));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
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
        private final Domain<E> domain;
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
         * @param domain the domain
         */
        public Interpolation(Domain<E> domain) {
            this.domain = domain;
            this.lins = UnivariatePolynomial.one(domain);
            this.poly = UnivariatePolynomial.one(domain);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public void update(E point, E value) {
            if (points.isEmpty()) {
                points.add(point);
                values.add(value);
                mixedRadix.add(value);
                poly.multiply(value);
                return;
            }
            E reciprocal = domain.subtract(point, points.get(0));
            E accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = domain.add(accumulator, domain.multiply(mixedRadix.get(i), reciprocal));
                reciprocal = domain.multiply(reciprocal, domain.subtract(point, points.get(i)));
            }
            if (domain.isZero(reciprocal))
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = domain.reciprocal(reciprocal);
            mixedRadix.add(domain.multiply(reciprocal, domain.subtract(value, accumulator)));

            lins.multiply(lins.createLinear(domain.negate(points.get(points.size() - 1)), domain.getOne()));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
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