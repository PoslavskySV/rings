package cc.r2.core.poly2;

import cc.r2.core.polynomial.MutablePolynomial;

import java.util.Arrays;

/**
 * Created by poslavsky on 15/02/2017.
 */
final class MutablePolynomialZ extends MutablePolynomialAbstract<MutablePolynomialZ> {
    private MutablePolynomialZ(long[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    /** copy constructor */
    private MutablePolynomialZ(long[] data, int degree) {
        this.data = data;
        this.degree = degree;
    }

    static MutablePolynomialZ create(long... data) {
        return new MutablePolynomialZ(data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    static MutablePolynomialZ createMonomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new MutablePolynomialZ(data);
    }

    static MutablePolynomialZ zero() {
        return create(0);
    }

    static MutablePolynomialZ one() {
        return create(1);
    }

    MutablePolynomialMod modulus(long modulus, boolean copy) {
        return MutablePolynomialMod.createSigned(modulus, copy ? data.clone() : data);
    }

    MutablePolynomialMod modulus(long modulus) {
        return modulus(modulus, true);
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} if some of the elements can't be exactly
     * divided by the {@code factor}
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    MutablePolynomialZ divideOrNull(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        LongModularArithmetics.MagicDivider magic = LongModularArithmetics.magicSigned(factor);
        for (int i = degree; i >= 0; --i) {
            long l = LongModularArithmetics.divideSignedFast(data[i], magic);
            if (l * factor != data[i])
                return null;
            data[i] = l;
        }
        return this;
    }

    @Override
    MutablePolynomialZ createFromArray(long[] data) {
        return new MutablePolynomialZ(data);
    }

    @Override
    MutablePolynomialZ createZero() {
        return new MutablePolynomialZ(new long[]{0}, 0);
    }

    @Override
    MutablePolynomialZ createOne() {
        return new MutablePolynomialZ(new long[]{1}, 0);
    }

    @Override
    MutablePolynomialZ derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        for (int i = degree; i > 0; --i)
            newData[i - 1] = LongArithmetics.safeMultiply(data[i], i);
        return createFromArray(newData);
    }

    @Override
    long evaluate(long point) {
        if (point == 0)
            return cc();
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = LongArithmetics.safeAdd(LongArithmetics.safeMultiply(res, point), data[i]);
        return res;
    }

    /**
     * Evaluates this poly at a give rational point {@code num/den}
     *
     * @param num point numerator
     * @param den point denominator
     * @return value at {@code num/den}
     * @throws ArithmeticException if the result is not integer
     */
    long evaluateAtRational(long num, long den) {
        if (num == 0)
            return cc();
        long res = 0;
        for (int i = degree; i >= 0; --i) {
            long x = LongArithmetics.safeMultiply(res, num);
            if (x % den != 0)
                throw new IllegalArgumentException("The answer is not integer");
            res = LongArithmetics.safeAdd(x / den, data[i]);
        }
        return res;
    }

    @Override
    MutablePolynomialZ add(MutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeAdd(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialZ addMonomial(long coefficient, int exponent) {
        if (coefficient == 0)
            return this;

        ensureCapacity(exponent);
        data[exponent] = LongArithmetics.safeAdd(data[exponent], coefficient);
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialZ addMul(MutablePolynomialZ oth, long factor) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeAdd(data[i], LongArithmetics.safeMultiply(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialZ subtract(MutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeSubtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialZ subtract(MutablePolynomialZ oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = LongArithmetics.safeSubtract(data[i], LongArithmetics.safeMultiply(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialZ negate() {
        for (int i = degree; i >= 0; --i)
            data[i] = -data[i];
        return this;
    }

    @Override
    MutablePolynomialZ multiply(long factor) {
        for (int i = degree; i >= 0; --i)
            data[i] = LongArithmetics.safeMultiply(factor, data[i]);
        return this;
    }

    @Override
    public MutablePolynomialZ clone() {
        return new MutablePolynomialZ(data.clone(), degree);
    }

    @Override
    MutablePolynomialZ multiply(MutablePolynomialZ oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        double rBound = norm1() * oth.norm1() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE)
            // we can apply fast integer arithmetic
            data = multiplyUnsafe(oth);
        else
            data = multiplySafe(oth);

        degree += oth.degree;
        fixDegree();
        return this;
    }

    /** switch algorithms */
    private long[] multiplyUnsafe(MutablePolynomialZ oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** switch algorithms */
    private long[] multiplySafe(MutablePolynomialZ oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassical(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsuba(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    @Override
    MutablePolynomialZ square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = norm1();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE)
            // we can apply fast integer arithmetic
            data = squareUnsafe();
        else
            data = squareSafe();

        degree += degree;
        fixDegree();
        return this;
    }

    /** switch algorithms */
    private long[] squareUnsafe() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalUnsafe(data, 0, degree + 1);
        else
            return squareKaratsubaUnsafe(data, 0, degree + 1);
    }

    /** switch algorithms */
    private long[] squareSafe() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassical(data, 0, degree + 1);
        else
            return squareKaratsuba(data, 0, degree + 1);
    }

    /**
     * Classical n*m multiplication algorithm
     *
     * @param result where to write the result
     * @param a      the first multiplier
     * @param aFrom  begin in a
     * @param aTo    end in a
     * @param b      the second multiplier
     * @param bFrom  begin in b
     * @param bTo    end in b
     */
    static void multiplyClassical(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassical(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = LongArithmetics.safeAdd(result[i + j], LongArithmetics.safeMultiply(c, b[bFrom + j]));
        }
    }

    /**
     * Classical n*m multiplication algorithm
     *
     * @param a     the first multiplier
     * @param aFrom begin in a
     * @param aTo   end in a
     * @param b     the second multiplier
     * @param bFrom begin in b
     * @param bTo   end in b
     * @return the result
     */
    static long[] multiplyClassical(final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        multiplyClassical(result, a, aFrom, aTo, b, bFrom, bTo);
        return result;
    }

    /**
     * Karatsuba multiplication
     *
     * @param f     the first multiplier
     * @param g     the second multiplier
     * @param fFrom begin in f
     * @param fTo   end in f
     * @param gFrom begin in g
     * @param gTo   end in g
     * @return the result
     */
    static long[] multiplyKaratsuba(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = LongArithmetics.safeMultiply(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = LongArithmetics.safeMultiply(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = LongArithmetics.safeMultiply(f[fFrom], g[gFrom]);
            result[1] = LongArithmetics.safeAdd(LongArithmetics.safeMultiply(f[fFrom], g[gFrom + 1]), LongArithmetics.safeMultiply(f[fFrom + 1], g[gFrom]));
            result[2] = LongArithmetics.safeMultiply(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if (LongArithmetics.safeMultiply(fTo - fFrom, gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyClassical(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyKaratsuba(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyKaratsuba(f, fFrom, fFrom + split, g, gFrom, gTo);
            long[] f1g = multiplyKaratsuba(f, fFrom + split, fTo, g, gFrom, gTo);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = LongArithmetics.safeAdd(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyKaratsuba(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyKaratsuba(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.safeAdd(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = LongArithmetics.safeAdd(g0_plus_g1[i - gMid], g[i]);

        long[] mid = multiplyKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.safeSubtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.safeSubtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.safeAdd(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.safeAdd(result[i + 2 * split], f1g1[i]);

        return result;
    }

    /** classical square */
    static long[] squareClassical(long[] a, int from, int to) {
        long[] x = new long[(to - from) * 2 - 1];
        squareClassical(x, a, from, to);
        return x;
    }

    /**
     * Square the poly {@code data} using classical algorithm
     *
     * @param result result destination
     * @param data   the data
     * @param from   data from
     * @param to     end point in the {@code data}
     */
    static void squareClassical(final long[] result, long[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[i + j] = LongArithmetics.safeAdd(result[i + j], LongArithmetics.safeMultiply(c, data[from + j]));
        }
    }

    /**
     * Karatsuba squaring
     *
     * @param f     the data
     * @param fFrom begin in f
     * @param fTo   end in f
     * @return the result
     */
    static long[] squareKaratsuba(final long[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return new long[0];
        if (fTo - fFrom == 1)
            return new long[]{LongArithmetics.safeMultiply(f[fFrom], f[fFrom])};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = LongArithmetics.safeMultiply(f[fFrom], f[fFrom]);
            result[1] = LongArithmetics.safeMultiply(2L, f[fFrom], f[fFrom + 1]);
            result[2] = LongArithmetics.safeMultiply(f[fFrom + 1], f[fFrom + 1]);
            return result;
        }
        //switch to classical
        if (LongArithmetics.safeMultiply(fTo - fFrom, fTo - fFrom) < KARATSUBA_THRESHOLD)
            return squareClassical(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareKaratsuba(f, fFrom, fMid);
        long[] f1g1 = squareKaratsuba(f, fMid, fTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.safeAdd(f0_plus_f1[i - fMid], f[i]);

        long[] mid = squareKaratsuba(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.safeSubtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.safeSubtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.safeAdd(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.safeAdd(result[i + 2 * split], f1g1[i]);

        return result;
    }
}
