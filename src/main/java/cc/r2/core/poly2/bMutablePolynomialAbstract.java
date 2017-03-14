package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.ZERO;
import static cc.r2.core.number.BigIntegerMath.*;

/**
 * Univariate polynomials over Z ({@link MutablePolynomialZ}) or Zp ({@link MutablePolynomialMod}).
 * All operations (except where it is specifically stated) changes the content of this.
 */
abstract class bMutablePolynomialAbstract<T extends bMutablePolynomialAbstract> implements Comparable<T> {
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    BigInteger[] data;
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
        while (i >= 0 && data[i].isZero()) --i;
        if (i < 0) i = 0;

        if (i != degree) {
            degree = i;
            Arrays.fill(data, degree + 1, data.length, ZERO);
        }
    }

    /**
     * Returns {@code true} if this is zero
     *
     * @return whether {@code this} is zero
     */
    final boolean isZero() {return data[degree].isZero();}

    /**
     * Returns {@code true} if this is one
     *
     * @return whether {@code this} is one
     */
    final boolean isOne() {return degree == 0 && data[0].isOne();}

    /**
     * Returns {@code true} if this polynomial is monic
     *
     * @return whether {@code this} is monic
     */
    final boolean isMonic() {return lc().isOne();}

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
            if (!data[i].isZero())
                return false;
        return true;
    }

    /**
     * Returns L1 norm of this polynomial, i.e. sum of abs coefficients
     *
     * @return L1 norm of {@code this}
     */
    final BigInteger norm1() {
        BigInteger norm = 0;
        for (int i = 0; i <= degree; ++i)
            norm = norm + abs(data[i]);
        return norm;
    }

//    /**
//     * Returns L2 norm of this polynomial, i.e. a square root of a sum of coefficient squares.
//     *
//     * @return L2 norm of {@code this}
//     */
//    final BigInteger norm2() {
//        BigInteger norm = 0;
//        for (int i = 0; i <= degree; ++i)
//            norm = norm + data[i] * data[i];
//        return norm.sqrtCeil();
//    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    final BigInteger normMax() {
        return maxAbsCoefficient();
    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    final BigInteger maxAbsCoefficient() {
        BigInteger max = abs(data[0]);
        for (int i = 1; i <= degree; ++i)
            max = max(abs(data[i]), max);
        return max;
    }

    /**
     * Returns the leading coefficient of the poly
     *
     * @return leading coefficient
     */
    final BigInteger lc() {return data[degree];}

    /**
     * Returns the constant coefficient of the poly
     *
     * @return constant coefficient
     */
    final BigInteger cc() {return data[0];}

    /**
     * Returns the content of the poly
     *
     * @return polynomial content
     */
    final BigInteger content() {
        if (degree == 0)
            return data[0];
        return gcd(data, 0, degree + 1);
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
        BigInteger content = content();
        if (content.isOne())
            return self;
        if (lc().signum() < 0)
            content = -content;
        for (int i = degree; i >= 0; --i)
            data[i] = data[i] / content;
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
    abstract T createFromArray(BigInteger[] data);

    /**
     * Creates constant polynomial with specified value
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    abstract T createConstant(BigInteger val);

    /**
     * Creates monomial {@code coefficient * x^degree}
     *
     * @param coefficient monomial coefficient
     * @param degree      monomial degree
     * @return {@code coefficient * x^degree}
     */
    abstract T createMonomial(BigInteger coefficient, int degree);

    /**
     * Returns 0 (new instance)
     *
     * @return new instance of 0
     */
    final T createZero() {return createConstant(ZERO);}

    /**
     * Returns 1 (new instance)
     *
     * @return new instance of 1
     */
    final T createOne() {return createConstant(ONE);}

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    abstract BigInteger evaluate(BigInteger point);

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
    abstract T addMonomial(BigInteger coefficient, int exponent);

    /**
     * Adds {@code oth * factor} to {@code this}
     *
     * @param oth    the polynomial
     * @param factor the factor
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    abstract T addMul(T oth, BigInteger factor);

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
    abstract T subtract(T oth, BigInteger factor, int exponent);

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
    abstract T multiply(BigInteger factor);

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
            c = data[i].compareTo(o.data[i]);
            if (c != 0)
                return c;
        }
        return 0;
    }

    @Override
    public final String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.length; i++) {
            if (data[i].isZero())
                continue;
            if (i != 0 && data[i].isOne()) {
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
        bMutablePolynomialAbstract oth = (bMutablePolynomialAbstract) obj;
        if (degree != oth.degree)
            return false;
        for (int i = 0; i <= degree; ++i)
            if (!data[i].equals(oth.data[i]))
                return false;
        return true;
    }

    @Override
    public final int hashCode() {
        int result = 1;
        for (int i = degree; i >= 0; --i) {
            int elementHash = data[i].hashCode();
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
    static void multiplyClassicalUnsafe(final BigInteger[] result, final BigInteger[] a, final int aFrom, final int aTo, final BigInteger[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassicalUnsafe(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            BigInteger c = a[aFrom + i];
            if (!c.isZero())
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
    static BigInteger[] multiplyClassicalUnsafe(final BigInteger[] a, final int aFrom, final int aTo, final BigInteger[] b, final int bFrom, final int bTo) {
        BigInteger[] result = new BigInteger[aTo - aFrom + bTo - bFrom - 1];
        Arrays.fill(result, ZERO);
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
    static BigInteger[] multiplyKaratsubaUnsafe(
            final BigInteger[] f, final int fFrom, final int fTo,
            final BigInteger[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new BigInteger[0];

        // single element in f
        if (fTo - fFrom == 1) {
            BigInteger[] result = new BigInteger[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = f[fFrom] * g[i];
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            BigInteger[] result = new BigInteger[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = g[gFrom] * f[i];
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            BigInteger[] result = new BigInteger[3];
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
            BigInteger[] f0g = multiplyKaratsubaUnsafe(f, fFrom, fFrom + split, g, gFrom, gTo);
            BigInteger[] f1g = multiplyKaratsubaUnsafe(f, fFrom + split, fTo, g, gFrom, gTo);

            int oldLen = f0g.length, newLen = fTo - fFrom + gTo - gFrom - 1;
            BigInteger[] result = Arrays.copyOf(f0g, newLen);
            Arrays.fill(result, oldLen, newLen, ZERO);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = result[i + split] + f1g[i];
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        BigInteger[] f0g0 = multiplyKaratsubaUnsafe(f, fFrom, fMid, g, gFrom, gMid);
        BigInteger[] f1g1 = multiplyKaratsubaUnsafe(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        BigInteger[] f0_plus_f1 = new BigInteger[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        Arrays.fill(f0_plus_f1, fMid - fFrom, f0_plus_f1.length, ZERO);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid] + f[i];

        //g0 + g1
        BigInteger[] g0_plus_g1 = new BigInteger[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        Arrays.fill(g0_plus_g1, gMid - gFrom, g0_plus_g1.length, ZERO);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = g0_plus_g1[i - gMid] + g[i];

        BigInteger[] mid = multiplyKaratsubaUnsafe(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f0g0.length);
            Arrays.fill(mid, oldLen, mid.length, ZERO);
        }
        if (mid.length < f1g1.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f1g1.length);
            Arrays.fill(mid, oldLen, mid.length, ZERO);
        }

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = mid[i] - f0g0[i];
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i] - f1g1[i];


        int oldLen = f0g0.length;
        BigInteger[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        Arrays.fill(result, oldLen, result.length, ZERO);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split] + mid[i];
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split] + f1g1[i];

        return result;
    }

    /** classical square */
    static BigInteger[] squareClassicalUnsafe(BigInteger[] a, int from, int to) {
        BigInteger[] x = new BigInteger[(to - from) * 2 - 1];
        Arrays.fill(x, ZERO);
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
    static void squareClassicalUnsafe(final BigInteger[] result, BigInteger[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            BigInteger c = data[from + i];
            if (!c.isZero())
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
    static BigInteger[] squareKaratsubaUnsafe(final BigInteger[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return new BigInteger[0];
        if (fTo - fFrom == 1)
            return new BigInteger[]{f[fFrom] * f[fFrom]};
        if (fTo - fFrom == 2) {
            BigInteger[] result = new BigInteger[3];
            result[0] = f[fFrom] * f[fFrom];
            result[1] = BigInteger.TWO * f[fFrom] * f[fFrom + 1];
            result[2] = f[fFrom + 1] * f[fFrom + 1];
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (fTo - fFrom) < KARATSUBA_THRESHOLD)
            return squareClassicalUnsafe(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        BigInteger[] f0g0 = squareKaratsubaUnsafe(f, fFrom, fMid);
        BigInteger[] f1g1 = squareKaratsubaUnsafe(f, fMid, fTo);

        // f0 + f1
        BigInteger[] f0_plus_f1 = new BigInteger[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        Arrays.fill(f0_plus_f1, fMid - fFrom, f0_plus_f1.length, ZERO);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid] + f[i];

        BigInteger[] mid = squareKaratsubaUnsafe(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f0g0.length);
            Arrays.fill(mid, oldLen, mid.length, ZERO);
        }
        if (mid.length < f1g1.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f1g1.length);
            Arrays.fill(mid, oldLen, mid.length, ZERO);
        }


        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = mid[i] - f0g0[i];
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i] - f1g1[i];


        int oldLen = f0g0.length;
        BigInteger[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        Arrays.fill(result, oldLen, result.length, ZERO);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split] + mid[i];
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split] + f1g1[i];

        return result;
    }
}
