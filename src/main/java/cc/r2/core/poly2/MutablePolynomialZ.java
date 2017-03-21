package cc.r2.core.poly2;

import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.Arrays;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Univariate polynomial over Z.
 * All operations (except where it is specifically stated) changes the content of this.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MutablePolynomialZ extends MutablePolynomialAbstract<MutablePolynomialZ> implements IMutablePolynomialZ<MutablePolynomialZ> {

    /** main constructor */
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

    /**
     * Creates Z[x] polynomial from the specified coefficients
     *
     * @param data coefficients
     * @return Z[x] polynomial
     */
    public static MutablePolynomialZ create(long... data) {
        return new MutablePolynomialZ(data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static MutablePolynomialZ monomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new MutablePolynomialZ(data);
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ createMonomial(long coefficient, int degree) {
        return monomial(coefficient, degree);
    }

    @Override
    public MutablePolynomialZ getRange(int from, int to) {
        return new MutablePolynomialZ(Arrays.copyOfRange(data, from, to));
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @return polynomial 0
     */
    public static MutablePolynomialZ zero() {
        return new MutablePolynomialZ(new long[]{0}, 0);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @return polynomial 1
     */
    public static MutablePolynomialZ one() {
        return new MutablePolynomialZ(new long[]{1}, 0);
    }

    public bMutablePolynomialZ toBigPoly() {
        return bMutablePolynomialZ.create(dataToBigIntegers());
    }

    /**
     * Returns Mignotte's bound (sqrt(n+1) * 2^n max |this|)
     *
     * @return Mignotte's bound
     */
    public double mignotteBound() {
        return Math.pow(2.0, degree) * norm2();
    }

    /**
     * Reduces polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @param copy    whether to copy the internal data or reduce inplace
     * @return Zp[x] polynomial from this
     */
    public MutablePolynomialMod modulus(long modulus, boolean copy) {
        return MutablePolynomialMod.createSigned(modulus, copy ? data.clone() : data);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @return Zp[x] polynomial from this
     */
    public MutablePolynomialMod modulus(long modulus) {
        return modulus(modulus, true);
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public MutablePolynomialZ divideOrNull(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        Magic magic = magicSigned(factor);
        for (int i = degree; i >= 0; --i) {
            long l = divideSignedFast(data[i], magic);
            if (l * factor != data[i])
                return null;
            data[i] = l;
        }
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ divideOrNullByLC(MutablePolynomialZ other) {
        return divideOrNull(other.lc());
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ createFromArray(long[] data) {
        return new MutablePolynomialZ(data);
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ createConstant(long val) {
        return new MutablePolynomialZ(new long[]{val}, 0);
    }

    @Override
    public MutablePolynomialZ[] arrayNewInstance(int length) {
        return new MutablePolynomialZ[length];
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        for (int i = degree; i > 0; --i)
            newData[i - 1] = LongArithmetics.safeMultiply(data[i], i);
        return createFromArray(newData);
    }

    /** {@inheritDoc} */
    @Override
    public long evaluate(long point) {
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
    public long evaluateAtRational(long num, long den) {
        if (num == 0)
            return cc();
        long res = 0;
        Magic magic = magicSigned(den);
        for (int i = degree; i >= 0; --i) {
            long x = LongArithmetics.safeMultiply(res, num);
            long q = divideSignedFast(x, magic);
            if (q * den != x)
                throw new IllegalArgumentException("The answer is not integer");
            res = LongArithmetics.safeAdd(q, data[i]);
        }
        return res;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ add(MutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeAdd(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ addMonomial(long coefficient, int exponent) {
        if (coefficient == 0)
            return this;

        ensureCapacity(exponent);
        data[exponent] = LongArithmetics.safeAdd(data[exponent], coefficient);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ addMul(MutablePolynomialZ oth, long factor) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeAdd(data[i], LongArithmetics.safeMultiply(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ subtract(MutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.safeSubtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ subtract(MutablePolynomialZ oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = LongArithmetics.safeSubtract(data[i], LongArithmetics.safeMultiply(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ negate() {
        for (int i = degree; i >= 0; --i)
            data[i] = -data[i];
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ multiply(long factor) {
        for (int i = degree; i >= 0; --i)
            data[i] = LongArithmetics.safeMultiply(factor, data[i]);
        return this;
    }

    public MutablePolynomialZ multiplyUnsafe(long factor) {
        for (int i = degree; i >= 0; --i)
            data[i] *= factor;
        return this;
    }

    @Override
    public MutablePolynomialZ clone() {
        return new MutablePolynomialZ(data.clone(), degree);
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ multiply(MutablePolynomialZ oth) {
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

        double rBound = normMax() * oth.normMax() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE)
            // we can apply fast integer arithmetic
            data = multiplyUnsafe0(oth);
        else
            data = multiplySafe(oth);

        degree += oth.degree;
        fixDegree();
        return this;
    }

    public MutablePolynomialZ multiplyUnsafe(MutablePolynomialZ oth) {
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
            return multiplyUnsafe(factor);
        }

        data = multiplyUnsafe0(oth);
        degree += oth.degree;
        fixDegree();
        return this;
    }

    /** switch algorithms */
    private long[] multiplyUnsafe0(MutablePolynomialZ oth) {
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

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialZ square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = normMax();
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
