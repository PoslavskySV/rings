package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.redberry.libdivide4j.FastDivision.*;

import java.util.Arrays;

import static cc.r2.core.poly2.LongArithmetics.fits32bitWord;
import static cc.r2.core.poly2.LongArithmetics.modInverse;
import static cc.redberry.libdivide4j.FastDivision.*;

/**
 * Univariate polynomial over Zp.
 * All operations (except where it is specifically stated) changes the content of this.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MutablePolynomialMod extends MutablePolynomialAbstract<MutablePolynomialMod> implements IMutablePolynomialZp<MutablePolynomialMod> {
    /** the modulus */
    final long modulus;
    /** magic **/
    private final Magic magic, magic32MulMod;
    /** whether modulus less then 2^32 (if so, faster mulmod available) **/
    private final boolean modulusFits32;

    /** copy constructor */
    private MutablePolynomialMod(long modulus, Magic magic, Magic magic32MulMod, long[] data, int degree, boolean modulusFits32) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32MulMod;
        this.data = data;
        this.degree = degree;
        this.modulusFits32 = modulusFits32;
    }

    /** main constructor */
    private MutablePolynomialMod(long modulus, Magic magic, long[] data) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32ForMultiplyMod(modulus);
        this.data = data;
        this.degree = data.length - 1;
        this.modulusFits32 = fits32bitWord(modulus);
        fixDegree();
    }

    /** Max supported modulus bits */
    public static final int MAX_SUPPORTED_MODULUS_BITS = 62;

    /** Max supported modulus */
    public static final long MAX_SUPPORTED_MODULUS = (1L << MAX_SUPPORTED_MODULUS_BITS) - 1L;

    private static void checkModulus(long modulus) {
        if (Long.compareUnsigned(modulus, MAX_SUPPORTED_MODULUS) > 0)
            throw new IllegalArgumentException("Too large modulus. Max allowed is " + MAX_SUPPORTED_MODULUS);
    }

    /* =========================== Factory methods =========================== */

    /**
     * Creates poly with specified coefficients represented as signed integers reducing them modulo {@code modulus}
     *
     * @param modulus the modulus
     * @param data    coefficients
     * @return the polynomial
     */
    public static MutablePolynomialMod createSigned(long modulus, long[] data) {
        reduceSigned(data, magicSigned(modulus));
        return new MutablePolynomialMod(modulus, magicUnsigned(modulus), data);
    }

    /** reduce data mod modulus **/
    private static void reduceUnsigned(long[] data, Magic magicUnsigned) {
        for (int i = 0; i < data.length; ++i)
            data[i] = modUnsignedFast(data[i], magicUnsigned);
    }

    /** reduce data mod modulus **/
    private static void reduceSigned(long[] data, Magic magicSigned) {
        for (int i = 0; i < data.length; ++i)
            data[i] = modSignedFast(data[i], magicSigned);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param modulus     the modulus
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static MutablePolynomialMod createMonomial(long modulus, long coefficient, int exponent) {
        Magic magic = magicUnsigned(modulus);
        if (coefficient >= modulus)
            coefficient = modUnsignedFast(coefficient, magic);
        else if (coefficient < 0)
            coefficient = LongArithmetics.mod(coefficient, modulus);

        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new MutablePolynomialMod(modulus, magic, data);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param modulus the modulus
     * @param value   the value
     * @return constant polynomial
     */
    public static MutablePolynomialMod constant(long modulus, long value) {
        Magic magic = magicUnsigned(modulus);
        if (value >= modulus)
            value = modUnsignedFast(value, magic);
        else if (value < 0)
            value = LongArithmetics.mod(value, modulus);

        return new MutablePolynomialMod(modulus, magic, new long[]{value});
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @param modulus the modulus
     * @return polynomial 1
     */
    public static MutablePolynomialMod one(long modulus) {
        return constant(modulus, 1);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @param modulus the modulus
     * @return polynomial 0
     */
    public static MutablePolynomialMod zero(long modulus) {
        return constant(modulus, 0);
    }


    /*=========================== Main methods ===========================*/

    /** modulus operation */
    long mod(long val) {
        return modUnsignedFast(val, magic);
    }

    /** multiplyMod operation */
    long multiplyMod(long a, long b) {
        return modulusFits32 ? mod(a * b) : multiplyMod128Unsigned(a, b, modulus, magic32MulMod);
    }

    /** addMod operation */
    long addMod(long a, long b) {
        long r = a + b;
        return r - modulus >= 0 ? r - modulus : r;
    }

    /** subtractMod operation */
    long subtractMod(long a, long b) {
        long r = a - b;
        return r + ((r >> 63)&modulus);
    }

    /** to symmetric modulus */
    long symmetricForm(long value) {
        return value <= modulus / 2 ? value : value - modulus;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod createFromArray(long[] newData) {
        reduceSigned(newData, magicSigned(modulus));
        MutablePolynomialMod r = new MutablePolynomialMod(modulus, magic, magic32MulMod, newData, newData.length - 1, modulusFits32);
        r.fixDegree();
        return r;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod createConstant(long val) {
        if (val < 0)
            val = LongArithmetics.mod(val, modulus);
        else if (val >= modulus)
            val = mod(val);
        return new MutablePolynomialMod(modulus, magic, magic32MulMod, new long[]{val}, 0, modulusFits32);
    }

    @Override
    public MutablePolynomialMod getRange(int from, int to) {
        MutablePolynomialMod r = new MutablePolynomialMod(modulus, magic, magic32MulMod, Arrays.copyOfRange(data, from, to), to - from - 1, modulusFits32);
        r.fixDegree();
        return r;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod createMonomial(long coefficient, int newDegree) {
        if (coefficient < 0)
            coefficient = LongArithmetics.mod(coefficient, modulus);
        else if (coefficient >= modulus)
            coefficient = mod(coefficient);

        long[] newData = new long[newDegree + 1];
        newData[newDegree] = coefficient;

        return new MutablePolynomialMod(modulus, magic, magic32MulMod, newData, newDegree, modulusFits32);
    }

    @Override
    public MutablePolynomialMod[] arrayNewInstance(int length) {
        return new MutablePolynomialMod[length];
    }

    @Override
    public MutablePolynomialMod[] arrayNewInstance(MutablePolynomialMod a, MutablePolynomialMod b) {
        return new MutablePolynomialMod[]{a, b};
    }

    /** does not copy the data and does not reduce the data with new modulus */
    public MutablePolynomialMod setModulusUnsafe(long newModulus) {
        return new MutablePolynomialMod(newModulus, magicUnsigned(newModulus), data);
    }

    /**
     * Creates new Zp[x] polynomial with specified modulus.
     *
     * @param newModulus the new modulus
     * @return the new Zp[x] polynomial with specified modulus
     */
    public MutablePolynomialMod setModulus(long newModulus) {
        long[] newData = data.clone();
        Magic newMagic = magicUnsigned(newModulus);
        reduceUnsigned(newData, newMagic);
        return new MutablePolynomialMod(newModulus, newMagic, newData);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this.
     *
     * @param copy whether to copy the internal data
     * @return Z[x] version of this
     */
    public MutablePolynomialZ normalForm(boolean copy) {
        return MutablePolynomialZ.create(copy ? data.clone() : data);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[x] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     */
    public MutablePolynomialZ normalSymmetricForm() {
        long[] newData = new long[degree + 1];
        for (int i = degree; i >= 0; --i)
            newData[i] = symmetricForm(data[i]);
        return MutablePolynomialZ.create(newData);
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod monic() {
        if (data[degree] == 0) // isZero()
            return this;
        if (degree == 0) {
            data[0] = 1;
            return this;
        }
        return multiply(modInverse(lc(), modulus));
    }

    @Override
    public MutablePolynomialMod divideByLC(MutablePolynomialMod other) {
        return multiply(modInverse(other.lc(), modulus));
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public MutablePolynomialMod monic(long factor) {
        return multiply(multiplyMod(mod(factor), modInverse(lc(), modulus)));
    }

    /** {@inheritDoc} */
    @Override
    public long evaluate(long point) {
        if (point == 0)
            return cc();

        point = mod(point);
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = addMod(multiplyMod(res, point), data[i]);
        return res;
    }

    @Override
    public void checkCompatibleModulus(MutablePolynomialMod oth) {
        if (modulus != oth.modulus)
            throw new IllegalArgumentException();
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod add(MutablePolynomialMod oth) {
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

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod addMonomial(long coefficient, int exponent) {
        if (coefficient == 0)
            return this;

        ensureCapacity(exponent);
        data[exponent] = addMod(data[exponent], mod(coefficient));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod addMul(MutablePolynomialMod oth, long factor) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor == 0)
            return this;

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], multiplyMod(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod subtract(MutablePolynomialMod oth) {
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

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod subtract(MutablePolynomialMod oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor == 0)
            return this;

        checkCompatibleModulus(oth);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtractMod(data[i], multiplyMod(factor, oth.data[i - exponent]));

        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod negate() {
        for (int i = degree; i >= 0; --i)
            if (data[i] != 0)
                data[i] = modulus - data[i];
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod multiply(long factor) {
        factor = mod(factor);
        if (factor == 1)
            return this;

        if (factor == 0)
            return toZero();

        for (int i = degree; i >= 0; --i)
            data[i] = multiplyMod(data[i], factor);
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        if (degree < modulus)
            for (int i = degree; i > 0; --i)
                newData[i - 1] = multiplyMod(data[i], i);
        else {
            int i = degree;
            for (; i >= modulus; --i)
                newData[i - 1] = multiplyMod(data[i], mod(i));
            for (; i > 0; --i)
                newData[i - 1] = multiplyMod(data[i], i);
        }
        return createFromArray(newData);
    }

    public bMutablePolynomialMod toBigPoly() {
        return bMutablePolynomialMod.createUnsafe(BigInteger.valueOf(modulus), dataToBigIntegers());
    }

    /** Deep copy */
    @Override
    public MutablePolynomialMod clone() {
        return new MutablePolynomialMod(modulus, magic, magic32MulMod, data.clone(), degree, modulusFits32);
    }

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod multiply(MutablePolynomialMod oth) {
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

        double rBound = normMax() * oth.normMax() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = multiplyExact(oth);
            degree += oth.degree;
            reduceUnsigned(data, magic);
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

    /** {@inheritDoc} */
    @Override
    public MutablePolynomialMod square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = normMax();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = squareExact();
            degree += degree;
            reduceUnsigned(data, magic);
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

    /**
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
                    result[i + j] = addMod(result[i + j], multiplyMod(c, b[bFrom + j]));
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
                result[i - gFrom] = multiplyMod(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = multiplyMod(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = multiplyMod(f[fFrom], g[gFrom]);
            result[1] = addMod(multiplyMod(f[fFrom], g[gFrom + 1]), multiplyMod(f[fFrom + 1], g[gFrom]));
            result[2] = multiplyMod(f[fFrom + 1], g[gFrom + 1]);
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
                    result[i + j] = addMod(result[i + j], multiplyMod(c, data[from + j]));
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
            return new long[]{multiplyMod(f[fFrom], f[fFrom])};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = multiplyMod(f[fFrom], f[fFrom]);
            result[1] = multiplyMod(multiplyMod(2L, f[fFrom]), f[fFrom + 1]);
            result[2] = multiplyMod(f[fFrom + 1], f[fFrom + 1]);
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
