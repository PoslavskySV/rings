package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import gnu.trove.list.array.TLongArrayList;

import java.util.ArrayList;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class PolynomialInterpolation {
    private PolynomialInterpolation() {}

    private static void checkInput(long[] points, long[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
    }

    /**
     * Constructs an interpolating polynomial which values at {@code points[i]} are exactly {@code values[i]}. This
     * method uses Lagrange's interpolation formula.
     *
     * @param points  evaluation points
     * @param values  corresponding polynomial values
     * @param modulus the modulus
     * @return the interpolating polynomial
     */
    public static lMutablePolynomialZp interpolateLagrange(long[] points, long[] values, long modulus) {
        checkInput(points, values);

        int length = points.length;
        lMutablePolynomialZp result = lMutablePolynomialZp.zero(modulus);
        for (int i = 0; i < length; ++i) {
            lMutablePolynomialZp interpolant = lMutablePolynomialZp.constant(modulus, values[i]);
            for (int j = 0; j < length; ++j) {
                if (j == i)
                    continue;
                lMutablePolynomialZp linear = result
                        .createLinear(interpolant.negateMod(points[j]), 1)
                        .divide(interpolant.subtractMod(points[i], points[j]));
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
    public static lMutablePolynomialZp interpolateNewton(long modulus, long[] points, long[] values) {
        checkInput(points, values);
        lInterpolation interpolation = new lInterpolation(modulus, points[0], values[0]);
        for (int i = 1; i < points.length; i++)
            interpolation.update(points[i], values[i]);
        return interpolation.getInterpolatingPolynomial();
//        int length = points.length;
//
//        // Newton's representation
//        long[] mixedRadix = new long[length];
//        mixedRadix[0] = values[0];
//        lMutablePolynomialZp lins = lMutablePolynomialZp.one(modulus);
//        lMutablePolynomialZp poly = lMutablePolynomialZp.constant(modulus, mixedRadix[0]);
//
//        for (int k = 1; k < length; ++k) {
//            long reciprocal = poly.subtractMod(points[k], points[0]);
//            long accumulator = mixedRadix[0];
//            for (int i = 1; i < k; ++i) {
//                accumulator = poly.addMod(accumulator, poly.multiplyMod(mixedRadix[i], reciprocal));
//                reciprocal = poly.multiplyMod(reciprocal, poly.subtractMod(points[k], points[i]));
//            }
//            reciprocal = LongArithmetics.modInverse(reciprocal, modulus);
//            mixedRadix[k] = poly.multiplyMod(reciprocal, poly.subtractMod(values[k], accumulator));
//
//            lins = lins.multiply(lins.createLinear(poly.negateMod(points[k - 1]), 1));
//            poly = poly.add(lins.clone().multiply(mixedRadix[k]));
//        }
//        return poly;
    }

    private static void checkInput(BigInteger[] points, BigInteger[] values) {
        if (points.length != values.length)
            throw new IllegalArgumentException();
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
    public static bMutablePolynomialZp interpolateNewton(BigInteger modulus, BigInteger[] points, BigInteger[] values) {
        checkInput(points, values);
        bInterpolation interpolation = new bInterpolation(modulus, points[0], values[0]);
        for (int i = 1; i < points.length; i++)
            interpolation.update(points[i], values[i]);
        return interpolation.getInterpolatingPolynomial();
    }

    /**
     * Updatable Newton interpolation
     */
    public static final class lInterpolation {
        private final long modulus;
        /** list of evaluation points */
        private final TLongArrayList points = new TLongArrayList();
        /** list of values at points */
        private final TLongArrayList values = new TLongArrayList();
        /** mixed radix form of interpolating polynomial */
        private final TLongArrayList mixedRadix = new TLongArrayList();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final lMutablePolynomialZp lins;
        /** resulting interpolating polynomial */
        private final lMutablePolynomialZp poly;

        /**
         * Start new interpolation with {@code interpolation[point] = value}
         *
         * @param modulus the modulus
         * @param point   evaluation point
         * @param value   polynomial value at {@code point}
         */
        public lInterpolation(long modulus, long point, long value) {
            this.modulus = modulus;
            this.lins = lMutablePolynomialZp.one(modulus);
            this.poly = lMutablePolynomialZp.constant(modulus, value);

            points.add(point);
            values.add(value);
            mixedRadix.add(value);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public void update(long point, long value) {
            long reciprocal = poly.subtractMod(point, points.get(0));
            long accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = poly.addMod(accumulator, poly.multiplyMod(mixedRadix.get(i), reciprocal));
                reciprocal = poly.multiplyMod(reciprocal, poly.subtractMod(point, points.get(i)));
            }
            if (reciprocal == 0)
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = LongArithmetics.modInverse(reciprocal, modulus);
            mixedRadix.add(poly.multiplyMod(reciprocal, poly.subtractMod(value, accumulator)));

            lins.multiply(lins.createLinear(poly.negateMod(points.get(points.size() - 1)), 1));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
        }

        /**
         * Returns resulting interpolating polynomial
         *
         * @return interpolating polynomial
         */
        public lMutablePolynomialZp getInterpolatingPolynomial() {return poly;}

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
    public static final class bInterpolation {
        private final BigInteger modulus;
        /** list of evaluation points */
        private final ArrayList<BigInteger> points = new ArrayList<>();
        /** list of values at points */
        private final ArrayList<BigInteger> values = new ArrayList<>();
        /** mixed radix form of interpolating polynomial */
        private final ArrayList<BigInteger> mixedRadix = new ArrayList<>();
        /** total modulus (x_i - points[0])*(x_i - points[1])*... */
        private final bMutablePolynomialZp lins;
        /** resulting interpolating polynomial */
        private final bMutablePolynomialZp poly;

        /**
         * Start new interpolation with {@code interpolation[point] = value}
         *
         * @param modulus the modulus
         * @param point   evaluation point
         * @param value   polynomial value at {@code point}
         */
        public bInterpolation(BigInteger modulus, BigInteger point, BigInteger value) {
            this.modulus = modulus;
            this.lins = bMutablePolynomialZp.one(modulus);
            this.poly = bMutablePolynomialZp.constant(modulus, value);

            points.add(point);
            values.add(value);
            mixedRadix.add(value);
        }

        /**
         * Updates interpolation, so that interpolating polynomial satisfies {@code interpolation[point] = value}
         *
         * @param point evaluation point
         * @param value polynomial value at {@code point}
         */
        public void update(BigInteger point, BigInteger value) {
            BigInteger reciprocal = poly.subtractMod(point, points.get(0));
            BigInteger accumulator = mixedRadix.get(0);
            for (int i = 1; i < points.size(); ++i) {
                accumulator = poly.addMod(accumulator, poly.multiplyMod(mixedRadix.get(i), reciprocal));
                reciprocal = poly.multiplyMod(reciprocal, poly.subtractMod(point, points.get(i)));
            }
            if (reciprocal.isZero())
                throw new IllegalArgumentException("Point " + point + " was already used in interpolation.");
            reciprocal = BigIntegerArithmetics.modInverse(reciprocal, modulus);
            mixedRadix.add(poly.multiplyMod(reciprocal, poly.subtractMod(value, accumulator)));

            lins.multiply(lins.createLinear(poly.negateMod(points.get(points.size() - 1)), BigInteger.ONE));
            poly.add(lins.clone().multiply(mixedRadix.get(mixedRadix.size() - 1)));

            points.add(point);
            values.add(value);
        }

        /**
         * Returns resulting interpolating polynomial
         *
         * @return interpolating polynomial
         */
        public bMutablePolynomialZp getInterpolatingPolynomial() {return poly;}

        /**
         * Returns the list of evaluation points used in interpolation
         *
         * @return list of evaluation points used in interpolation
         */
        public ArrayList<BigInteger> getPoints() {return points;}

        /**
         * Returns the list of polynomial values at interpolation points
         *
         * @return the list of polynomial values at interpolation points
         */
        public ArrayList<BigInteger> getValues() {return values;}

        /**
         * Returns the number of interpolation points used
         *
         * @return number of interpolation points used
         */
        public int numberOfPoints() {return points.size();}
    }
}
