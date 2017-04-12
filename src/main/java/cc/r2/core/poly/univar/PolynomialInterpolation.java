package cc.r2.core.poly.univar;

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
     * @param points  evaluation points
     * @param values  corresponding polynomial values
     * @param modulus the modulus
     * @return the interpolating polynomial
     */
    public static lMutablePolynomialZp interpolateNewton(long[] points, long[] values, long modulus) {
        checkInput(points, values);
        int length = points.length;

        // Newton's representation
        long[] mixedRadix = new long[length];
        mixedRadix[0] = values[0];
        lMutablePolynomialZp lins = lMutablePolynomialZp.one(modulus);
        lMutablePolynomialZp poly = lMutablePolynomialZp.constant(modulus, mixedRadix[0]);

        for (int k = 1; k < length; ++k) {
            long reciprocal = poly.subtractMod(points[k], points[0]);
            long accumulator = mixedRadix[0];
            for (int i = 1; i < k; ++i) {
                accumulator = poly.addMod(accumulator, poly.multiplyMod(mixedRadix[i], reciprocal));
                reciprocal = poly.multiplyMod(reciprocal, poly.subtractMod(points[k], points[i]));
            }
            reciprocal = LongArithmetics.modInverse(reciprocal, modulus);
            mixedRadix[k] = poly.multiplyMod(reciprocal, poly.subtractMod(values[k], accumulator));

            lins = lins.multiply(lins.createLinear(poly.negateMod(points[k - 1]), 1));
            poly = poly.add(lins.clone().multiply(mixedRadix[k]));
        }
        return poly;
    }
}
