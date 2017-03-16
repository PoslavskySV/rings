package cc.r2.core.poly2;

import cc.r2.core.util.ArraysUtil;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.Arrays;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Univariate polynomials over Z ({@link MutablePolynomialZ}) or Zp ({@link MutablePolynomialMod}).
 * All operations (except where it is specifically stated) changes the content of this.
 */
abstract class MutablePolynomialAbstract<T extends MutablePolynomialAbstract> implements Comparable<T> {
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    long[] data;
    /** points to the last non zero element in the data array */
    int degree;

    /**
     * Ensures that the capacity of internal storage is enough for storing polynomial of the {@code desiredDegree}.
     * The degree of {@code this} is set to {@code desiredDegree} if the latter is greater than the former.
     *
     * @param desiredDegree desired degree
     */
    void ensureCapacity(int desiredDegree) {
        if (degree < desiredDegree)
            degree = desiredDegree;

        if (data.length < (desiredDegree + 1))
            data = Arrays.copyOf(data, desiredDegree + 1);
    }

    /**
     * Removes zeroes from the end of {@code data} and adjusts the degree
     */
    void fixDegree() {
        int i = degree;
        while (i >= 0 && data[i] == 0) --i;
        if (i < 0) i = 0;

        if (i != degree) {
            degree = i;
            Arrays.fill(data, degree + 1, data.length, 0);
        }
    }

    /**
     * Returns {@code true} if this is zero
     *
     * @return whether {@code this} is zero
     */
    final boolean isZero() {return data[degree] == 0;}

    /**
     * Returns {@code true} if this is one
     *
     * @return whether {@code this} is one
     */
    final boolean isOne() {return degree == 0 && data[0] == 1;}

    /**
     * Returns {@code true} if this polynomial is monic
     *
     * @return whether {@code this} is monic
     */
    final boolean isMonic() {return lc() == 1;}

    /**
     * Returns {@code true} if this polynomial has only constant term
     *
     * @return whether {@code this} is constant
     */
    final boolean isConstant() {return degree == 0;}

    /**
     * Returns {@code true} if this polynomial has only one monomial term
     *
     * @return whether {@code this} has the form {@code c*x^i} (one term)
     */
    final boolean isMonomial() {
        for (int i = degree - 1; i >= 0; --i)
            if (data[i] != 0)
                return false;
        return true;
    }

    /**
     * Returns L1 norm of this polynomial, i.e. sum of abs coefficients
     *
     * @return L1 norm of {@code this}
     */
    final double norm1() {
        double norm = 0;
        for (int i = 0; i <= degree; ++i)
            norm += Math.abs(data[i]);
        return norm;
    }

    /**
     * Returns L2 norm of this polynomial, i.e. a square root of a sum of coefficient squares.
     *
     * @return L2 norm of {@code this}
     */
    final double norm2() {
        double norm = 0;
        for (int i = 0; i <= degree; ++i)
            norm += ((double) data[i]) * data[i];
        return Math.ceil(Math.sqrt(norm));
    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    final double normMax() {
        return (double) maxAbsCoefficient();
    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    final long maxAbsCoefficient() {
        long max = Math.abs(data[0]);
        for (int i = 1; i <= degree; ++i)
            max = Math.max(Math.abs(data[i]), max);
        return max;
    }

    /**
     * Returns the leading coefficient of the poly
     *
     * @return leading coefficient
     */
    final long lc() {return data[degree];}

    /**
     * Returns the constant coefficient of the poly
     *
     * @return constant coefficient
     */
    final long cc() {return data[0];}

    /**
     * Returns the content of the poly
     *
     * @return polynomial content
     */
    final long content() {
        if (degree == 0)
            return data[0];
        return LongArithmetics.gcd(data, 0, degree + 1);
    }

    /** casted self **/
    @SuppressWarnings("unchecked")
    private final T self = (T) this;

    /** fill content with zeroes */
    final T toZero() {
        Arrays.fill(data, 0, degree + 1, 0);
        degree = 0;
        return self;
    }

    /**
     * Sets the content of this to {@code oth}
     *
     * @param oth the polynomial
     * @return this := oth
     */
    final T set(T oth) {
        this.data = oth.data.clone();
        this.degree = oth.degree;
        return self;
    }

    /**
     * Returns the quotient {@code this / x^offset}, it is polynomial with coefficient list formed by shifting coefficients
     * of {@code this} to the left by {@code offset}.
     *
     * @param offset shift amount
     * @return the quotient {@code this / x^offset}
     */
    final T shiftLeft(int offset) {
        if (offset == 0)
            return self;
        if (offset > degree)
            return toZero();

        System.arraycopy(data, offset, data, 0, degree - offset + 1);
        Arrays.fill(data, degree - offset + 1, degree + 1, 0);
        degree = degree - offset;
        return self;
    }

    /**
     * Multiplies {@code this} by the {@code x^offset}.
     *
     * @param offset monomial exponent
     * @return {@code this * x^offset}
     */
    final T shiftRight(int offset) {
        if (offset == 0)
            return self;
        int degree = this.degree;
        ensureCapacity(offset + degree);
        System.arraycopy(data, 0, data, offset, degree + 1);
        Arrays.fill(data, 0, offset, 0);
        return self;
    }

    /**
     * Returns the remainder {@code this rem x^(newDegree + 1)}, it is polynomial with the coefficient list truncated
     * to the new degree {@code newDegree}.
     *
     * @param newDegree new degree
     * @return remainder {@code this rem x^(newDegree + 1)}
     */
    final T truncate(int newDegree) {
        if (newDegree >= degree)
            return self;
        Arrays.fill(data, newDegree + 1, degree + 1, 0);
        degree = newDegree;
        fixDegree();
        return self;
    }

    /**
     * Reverses the coefficients of this
     *
     * @return reversed polynomial
     */
    final T reverse() {
        ArraysUtil.reverse(data, 0, degree + 1);
        fixDegree();
        return self;
    }

    /**
     * Reduces poly to its primitive part
     *
     * @return primitive part (poly will be modified)
     */
    final T primitivePart() {
        long content = content();
        if (content == 1)
            return self;
        if (lc() < 0)
            content = -content;
        Magic magic = magicSigned(content);
        for (int i = degree; i >= 0; --i)
            data[i] = divideSignedFast(data[i], magic);
        return self;
    }

    /**
     * Adds 1 to this
     *
     * @return {@code this + 1}
     */
    final T increment() {
        return add(createOne());
    }

    /**
     * Subtracts 1 from this
     *
     * @return {@code this - 1}
     */
    final T decrement() {
        return subtract(createOne());
    }

    /**
     * Factory
     *
     * @param data the data
     * @return polynomial
     */
    abstract T createFromArray(long[] data);

    /**
     * Creates constant polynomial with specified value
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    abstract T createConstant(long val);

    /**
     * Creates monomial {@code coefficient * x^degree}
     *
     * @param coefficient monomial coefficient
     * @param degree      monomial degree
     * @return {@code coefficient * x^degree}
     */
    abstract T createMonomial(long coefficient, int degree);

    /**
     * Returns 0 (new instance)
     *
     * @return new instance of 0
     */
    final T createZero() {return createConstant(0);}

    /**
     * Returns 1 (new instance)
     *
     * @return new instance of 1
     */
    final T createOne() {return createConstant(1);}

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    abstract long evaluate(long point);

    /**
     * Adds {@code oth} to {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this + oth}
     */
    abstract T add(T oth);

    /**
     * Adds {@code coefficient*x^exponent} to {@code this}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code this + coefficient*x^exponent}
     */
    abstract T addMonomial(long coefficient, int exponent);

    /**
     * Adds {@code oth * factor} to {@code this}
     *
     * @param oth    the polynomial
     * @param factor the factor
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    abstract T addMul(T oth, long factor);

    /**
     * Subtracts {@code oth} from {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this - oth}
     */
    abstract T subtract(T oth);

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this}
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @return {@code this - factor * x^exponent * oth}
     */
    abstract T subtract(T oth, long factor, int exponent);

    /**
     * Negates this and returns
     *
     * @return this negated
     */
    abstract T negate();

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    abstract T multiply(long factor);

    /**
     * Sets this to {@code this * oth }
     *
     * @param oth the polynomial
     * @return {@code this * oth }
     */
    abstract T multiply(T oth);

    /**
     * Square of {@code this}
     *
     * @return {@code this * this}
     */
    abstract T square();

    /**
     * Returns the formal derivative of this poly
     *
     * @return the formal derivative
     */
    abstract T derivative();

    @Override
    public abstract T clone();

    @Override
    public final int compareTo(T o) {
        int c = Integer.compare(degree, o.degree);
        if (c != 0)
            return c;
        for (int i = degree; i >= 0; --i) {
            c = Long.compare(data[i], o.data[i]);
            if (c != 0)
                return c;
        }
        return 0;
    }

    @Override
    public final String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.length; i++) {
            if (data[i] == 0)
                continue;
            if (i != 0 && data[i] == 1) {
                if (sb.length() != 0)
                    sb.append("+");
                sb.append("x^").append(i);
            } else {
                String c = String.valueOf(data[i]);
                if (!c.startsWith("-") && sb.length() != 0)
                    sb.append("+");
                sb.append(c);
                if (i != 0)
                    sb.append("x^").append(i);
            }
        }

        if (sb.length() == 0)
            return "0";
        return sb.toString();
    }

    String toStringForCopy() {
        String s = ArraysUtil.toString(data, 0, degree + 1);
        return "create(" + s.substring(1, s.length() - 1) + ")";
    }

    @Override
    public final boolean equals(Object obj) {
        if (obj.getClass() != this.getClass())
            return false;
        MutablePolynomialAbstract oth = (MutablePolynomialAbstract) obj;
        if (degree != oth.degree)
            return false;
        for (int i = 0; i <= degree; ++i)
            if (data[i] != oth.data[i])
                return false;
        return true;
    }

    @Override
    public final int hashCode() {
        int result = 1;
        for (int i = degree; i >= 0; --i) {
            long element = data[i];
            int elementHash = (int) (element^(element >>> 32));
            result = 31 * result + elementHash;
        }
        return result;
    }

    /* *
    *
    * Exact multiplication with unsafe arithmetic
    *
    * */

    /** switch to classical multiplication */
    static final long KARATSUBA_THRESHOLD = 1024L;
    /** when use Karatsuba fast multiplication */
    static final long
            MUL_CLASSICAL_THRESHOLD = 256L * 256L,
            MUL_MOD_CLASSICAL_THRESHOLD = 128L * 128L;

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
    static void multiplyClassicalUnsafe(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassicalUnsafe(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = result[i + j] + c * b[bFrom + j];
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
    static long[] multiplyClassicalUnsafe(final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        multiplyClassicalUnsafe(result, a, aFrom, aTo, b, bFrom, bTo);
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
    static long[] multiplyKaratsubaUnsafe(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = f[fFrom] * g[i];
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = g[gFrom] * f[i];
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = f[fFrom] * g[gFrom];
            result[1] = f[fFrom] * g[gFrom + 1] + f[fFrom + 1] * g[gFrom];
            result[2] = f[fFrom + 1] * g[gFrom + 1];
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyClassicalUnsafe(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyKaratsubaUnsafe(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyKaratsubaUnsafe(f, fFrom, fFrom + split, g, gFrom, gTo);
            long[] f1g = multiplyKaratsubaUnsafe(f, fFrom + split, fTo, g, gFrom, gTo);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = result[i + split] + f1g[i];
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyKaratsubaUnsafe(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyKaratsubaUnsafe(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid] + f[i];

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = g0_plus_g1[i - gMid] + g[i];

        long[] mid = multiplyKaratsubaUnsafe(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = mid[i] - f0g0[i];
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i] - f1g1[i];


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split] + mid[i];
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split] + f1g1[i];

        return result;
    }

    /** classical square */
    static long[] squareClassicalUnsafe(long[] a, int from, int to) {
        long[] x = new long[(to - from) * 2 - 1];
        squareClassicalUnsafe(x, a, from, to);
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
    static void squareClassicalUnsafe(final long[] result, long[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[i + j] = result[i + j] + c * data[from + j];
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
    static long[] squareKaratsubaUnsafe(final long[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return new long[0];
        if (fTo - fFrom == 1)
            return new long[]{f[fFrom] * f[fFrom]};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = f[fFrom] * f[fFrom];
            result[1] = 2L * f[fFrom] * f[fFrom + 1];
            result[2] = f[fFrom + 1] * f[fFrom + 1];
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (fTo - fFrom) < KARATSUBA_THRESHOLD)
            return MutablePolynomialMod.squareClassicalUnsafe(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareKaratsubaUnsafe(f, fFrom, fMid);
        long[] f1g1 = squareKaratsubaUnsafe(f, fMid, fTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid] + f[i];

        long[] mid = squareKaratsubaUnsafe(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = mid[i] - f0g0[i];
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i] - f1g1[i];


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split] + mid[i];
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split] + f1g1[i];

        return result;
    }
}
