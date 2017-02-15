package cc.r2.core.poly2;

//import cc.r2.core.polynomial.LongArithmetics;

import cc.r2.core.poly2.LongModularArithmetics.*;

import java.util.Arrays;

import static cc.r2.core.poly2.LongModularArithmetics.*;
import static cc.r2.core.polynomial.LongArithmetics.modInverse;

/**
 * Created by poslavsky on 26/01/2017.
 */
final class MutablePolynomialMod extends MutablePolynomialAbstract<MutablePolynomialMod> {
    /** the modulus */
    final long modulus;
    /** magic **/
    final MagicDivider magic, magic32MulMod;
    /** whether modulus less then 2^32 **/
    final boolean modulusFits32;

    /** copy constructor */
    private MutablePolynomialMod(long modulus, MagicDivider magic, MagicDivider magic32MulMod, long[] data, int degree, boolean modulusFits32) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32MulMod;
        this.data = data;
        this.degree = degree;
        this.modulusFits32 = modulusFits32;
    }

    /** private constructor */
    private MutablePolynomialMod(long modulus, MagicDivider magic, long[] data) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32ForMultiplyMod(modulus);
        this.data = data;
        this.degree = data.length - 1;
        this.modulusFits32 = Long.compareUnsigned(modulus, 1L << 32) < 0;
        fixDegree();
    }

    private static void checkModulus(long modulus) {
        if (Long.compareUnsigned(modulus, LongModularArithmetics.MAX_SUPPORTED_MODULUS) > 0)
            throw new IllegalArgumentException("Too large modulus. Max allowed is " + MAX_SUPPORTED_MODULUS);
    }

    /**
     * Creates poly with specified coefficients
     *
     * @param data coefficients
     * @return the polynomial
     */
    static MutablePolynomialMod create(long modulus, long... data) {
        MagicDivider magic = magicUnsigned(modulus);
        reduce(data, magic);
        return new MutablePolynomialMod(modulus, magic, data);
    }

    /**
     * Creates poly with specified signed coefficients
     *
     * @param data coefficients
     * @return the polynomial
     */
    static MutablePolynomialMod createSigned(long modulus, long... data) {
        MagicDivider magic = magicSigned(modulus);
        for (int i = 0; i < data.length; ++i)
            data[i] = modSignedFast(data[i], magic);
        return new MutablePolynomialMod(modulus, magicUnsigned(modulus), data);
    }

    /** reduce data mod modulus **/
    private static void reduce(long[] data, MagicDivider magic) {
        for (int i = 0; i < data.length; ++i)
            data[i] = modUnsignedFast(data[i], magic);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param modulus     the modulus
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    static MutablePolynomialMod createMonomial(long modulus, long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        MagicDivider magic = magicUnsigned(modulus);
        data[exponent] = modUnsignedFast(coefficient, magic);
        return new MutablePolynomialMod(modulus, magic, data);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @param modulus the modulus
     * @return polynomial 1
     */
    static MutablePolynomialMod one(long modulus) {
        return create(modulus, 1);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @param modulus the modulus
     * @return polynomial 0
     */
    static MutablePolynomialMod zero(long modulus) {
        return create(modulus, 0);
    }

    /** modulus operation */
    private long mod(long val) { return modUnsignedFast(val, magic);}

    /** mulMod operation */
    private long mulMod(long a, long b) {
        return modulusFits32 ? mod(a * b) : multiplyMod128Unsigned(a, b, modulus, magic32MulMod);
    }

    /** addMod operation */
    private long addMod(long a, long b) {
        long r = a + b;
        return r - modulus >= 0 ? r - modulus : r;
    }

    /** subtractMod operation */
    private long subtractMod(long a, long b) {
        long r = a - b;
        return r + ((r >> 63)&modulus);
    }

    private void checkCompatibleModulus(MutablePolynomialMod oth) {
        if (modulus != oth.modulus)
            throw new IllegalArgumentException();
    }

    @Override
    MutablePolynomialMod create(long[] data) {
        return create(this.modulus, data);
    }

    private MutablePolynomialMod constant(long val) {
        return new MutablePolynomialMod(modulus, magic, magic32MulMod, new long[]{val}, 0, modulusFits32);
    }

    @Override
    MutablePolynomialMod zero() {
        return constant(0);
    }

    @Override
    MutablePolynomialMod one() {
        return constant(1);
    }

    /**
     * Reduces this polynomial and returns its monic part (that is {@code this} multiplied by
     * its inversed leading coefficient).
     *
     * @return {@code this} as a monic polynomial
     */
    MutablePolynomialMod monic() {
        if (data[degree] == 0) // isZero()
            return this;
        if (degree == 0) {
            data[0] = 1;
            return this;
        }
        return multiply(modInverse(lc(), modulus));
    }

    @Override
    long evaluate(long point) {
        if (point == 0)
            return cc();

        point = mod(point);
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = mod(mulMod(res, point) + data[i]);
        return res;
    }


    @Override
    MutablePolynomialMod add(MutablePolynomialMod oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth);

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialMod addMonomial(long coefficient, int exponent) {
        if (coefficient == 0)
            return this;

        ensureCapacity(exponent);
        data[exponent] = addMod(data[exponent], mod(coefficient));
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialMod addMul(MutablePolynomialMod oth, long factor) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor == 0)
            return this;

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], mulMod(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialMod subtract(MutablePolynomialMod oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth).negate();

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = subtractMod(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialMod subtract(MutablePolynomialMod oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor == 0)
            return this;

        checkCompatibleModulus(oth);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtractMod(data[i], mulMod(factor, oth.data[i - exponent]));

        fixDegree();
        return this;
    }

    @Override
    MutablePolynomialMod negate() {
        for (int i = degree; i >= 0; --i)
            data[i] = modulus - data[i];
        return this;
    }

    @Override
    MutablePolynomialMod multiply(long factor) {
        factor = mod(factor);
        if (factor == 1)
            return this;

        if (factor == 0)
            return toZero();

        for (int i = degree; i >= 0; --i)
            data[i] = mulMod(data[i], factor);
        return this;
    }

    @Override
    MutablePolynomialMod derivative() {
        if (isConstant())
            return zero();
        long[] newData = new long[degree];
        if (degree < modulus)
            for (int i = degree; i > 0; --i)
                newData[i - 1] = mulMod(data[i], i);
        else {
            int i = degree;
            for (; i >= modulus; --i)
                newData[i - 1] = mulMod(data[i], mod(i));
            for (; i > 0; --i)
                newData[i - 1] = mulMod(data[i], i);
        }
        return create(newData);
    }

    /**
     * Deep copy
     */
    @Override
    public MutablePolynomialMod clone() {
        return new MutablePolynomialMod(modulus, magic, magic32MulMod, data.clone(), degree, modulusFits32);
    }

    @Override
    MutablePolynomialMod multiply(MutablePolynomialMod oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        checkCompatibleModulus(oth);
        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        double rBound = norm1() * oth.norm1() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = multiplyExact(oth);
            degree += oth.degree;
            reduce(data, magic);
            fixDegree();
        } else {
            data = multiplyMod(oth);
            degree += oth.degree;
            fixDegree();
        }
        return this;
    }


    /** switch algorithms */
    private long[] multiplyExact(MutablePolynomialMod oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** switch algorithms */
    private long[] multiplyMod(MutablePolynomialMod oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyModClassical(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyModKaratsuba(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    @Override
    MutablePolynomialMod square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = norm1();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = squareExact();
            degree += degree;
            reduce(data, magic);
            fixDegree();
        } else {
            data = squareMod();
            degree += degree;
            fixDegree();
        }
        return this;
    }

    /** switch algorithms */
    private long[] squareExact() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalUnsafe(data, 0, degree + 1);
        else
            return squareKaratsubaUnsafe(data, 0, degree + 1);
    }

    /** switch algorithms */
    private long[] squareMod() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareModClassical(data, 0, degree + 1);
        else
            return squareModKaratsuba(data, 0, degree + 1);
    }

    /* *
    *
    * Multiplication with safe modular arithmetic
    *
    * */

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
    long[] multiplyModClassical(final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        multiplyModClassical(result, a, aFrom, aTo, b, bFrom, bTo);
        return result;
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
    void multiplyModClassical(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyModClassical(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = addMod(result[i + j], mulMod(c, b[bFrom + j]));
        }
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
    long[] multiplyModKaratsuba(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = mulMod(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = mulMod(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = mulMod(f[fFrom], g[gFrom]);
            result[1] = addMod(mulMod(f[fFrom], g[gFrom + 1]), mulMod(f[fFrom + 1], g[gFrom]));
            result[2] = mulMod(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyModClassical(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyModKaratsuba(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyModKaratsuba(f, fFrom, fFrom + split, g, gFrom, gTo);
            long[] f1g = multiplyModKaratsuba(f, fFrom + split, fTo, g, gFrom, gTo);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = addMod(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyModKaratsuba(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyModKaratsuba(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = addMod(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = addMod(g0_plus_g1[i - gMid], g[i]);

        long[] mid = multiplyModKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = subtractMod(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = subtractMod(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = addMod(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = addMod(result[i + 2 * split], f1g1[i]);

        return result;
    }

    long[] squareModClassical(long[] a, int from, int to) {
        long[] x = new long[(to - from) * 2 - 1];
        squareModClassical(x, a, from, to);
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
    void squareModClassical(final long[] result, long[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[i + j] = addMod(result[i + j], mulMod(c, data[from + j]));
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
    long[] squareModKaratsuba(final long[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return new long[0];
        if (fTo - fFrom == 1)
            return new long[]{mulMod(f[fFrom], f[fFrom])};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = mulMod(f[fFrom], f[fFrom]);
            result[1] = mulMod(mulMod(2L, f[fFrom]), f[fFrom + 1]);
            result[2] = mulMod(f[fFrom + 1], f[fFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (fTo - fFrom) < KARATSUBA_THRESHOLD)
            return squareModClassical(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareModKaratsuba(f, fFrom, fMid);
        long[] f1g1 = squareModKaratsuba(f, fMid, fTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = addMod(f0_plus_f1[i - fMid], f[i]);

        long[] mid = squareModKaratsuba(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = subtractMod(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = subtractMod(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = addMod(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = addMod(result[i + 2 * split], f1g1[i]);

        return result;
    }
}
