package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.ZERO;
import static cc.r2.core.number.BigIntegerArithmetics.*;

/**
 * Univariate polynomials over Z ({@link lMutablePolynomialZ}) or Zp ({@link lMutablePolynomialZp}).
 * All operations (except where it is specifically stated) changes the content of this.
 */
abstract class bMutablePolynomialAbstract<T extends bMutablePolynomialAbstract> implements IMutablePolynomial<T> {
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    BigInteger[] data;
    /** points to the last non zero element in the data array */
    int degree;

    /** {@inheritDoc} */
    @Override
    public final int degree() {return degree;}

    /**
     * Returns i-th element of this poly
     */
    public final BigInteger get(int i) { return data[i];}

    /**
     * Ensures that the capacity of internal storage is enough for storing polynomial of the {@code desiredDegree}.
     * The degree of {@code this} is set to {@code desiredDegree} if the latter is greater than the former.
     *
     * @param desiredDegree desired degree
     */
    void ensureCapacity(int desiredDegree) {
        if (degree < desiredDegree)
            degree = desiredDegree;

        if (data.length < (desiredDegree + 1)) {
            int oldLen = data.length;
            data = Arrays.copyOf(data, desiredDegree + 1);
            Arrays.fill(data, oldLen, data.length, ZERO);
        }
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

    /** {@inheritDoc} */
    @Override
    public final boolean isZero() {return data[degree].isZero();}

    /** {@inheritDoc} */
    @Override
    public final boolean isOne() {return degree == 0 && data[0].isOne();}

    /** {@inheritDoc} */
    @Override
    public final boolean isMonic() {return lc().isOne();}

    /** {@inheritDoc} */
    @Override
    public final boolean isUnitCC() {return cc().isOne();}

    /** {@inheritDoc} */
    @Override
    public final boolean isConstant() {return degree == 0;}

    /** {@inheritDoc} */
    @Override
    public final boolean isMonomial() {
        for (int i = degree - 1; i >= 0; --i)
            if (!data[i].isZero())
                return false;
        return true;
    }

    @Override
    public int signum() {
        return lc().signum();
    }

    @Override
    public int firstNonZeroCoefficientPosition() {
        int i = 0;
        while (data[i].isZero()) ++i;
        assert i < data.length;
        return i;
    }

    protected long[] toLong0() {
        long[] lData = new long[degree + 1];
        for (int i = degree; i >= 0; --i)
            lData[i] = data[i].longValueExact();
        return lData;
    }

    /**
     * Returns L1 norm of this polynomial, i.e. sum of abs coefficients
     *
     * @return L1 norm of {@code this}
     */
    public final BigInteger norm1() {
        BigInteger norm = ZERO;
        for (int i = 0; i <= degree; ++i)
            norm = norm.add(abs(data[i]));
        return norm;
    }

    /**
     * Returns L2 norm of this polynomial, i.e. a square root of a sum of coefficient squares.
     *
     * @return L2 norm of {@code this}
     */
    final BigInteger norm2() {
        BigInteger norm = ZERO;
        for (int i = 0; i <= degree; ++i)
            norm = norm.add(data[i].multiply(data[i]));
        return BigIntegerArithmetics.sqrtCeil(norm);
    }

    /**
     * Returns L2 norm of this polynomial, i.e. a square root of a sum of coefficient squares.
     *
     * @return L2 norm of {@code this}
     */
    final double norm2Double() {
        double norm = 0;
        for (int i = 0; i <= degree; ++i) {
            double d = data[i].doubleValue();
            norm += d * d;
        }
        return Math.sqrt(norm);
    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    public final BigInteger normMax() {
        return maxAbsCoefficient();
    }

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    public final BigInteger maxAbsCoefficient() {
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
    public final BigInteger lc() {return data[degree];}

    /**
     * Returns the constant coefficient of the poly
     *
     * @return constant coefficient
     */
    public final BigInteger cc() {return data[0];}

    /**
     * Returns the content of the poly
     *
     * @return polynomial content
     */
    public final BigInteger content() {
        if (degree == 0)
            return data[0];
        return gcd(data, 0, degree + 1);
    }

    /** casted self **/
    @SuppressWarnings("unchecked")
    private final T self = (T) this;

    /** {@inheritDoc} */
    @Override
    public final T toZero() {
        Arrays.fill(data, 0, degree + 1, ZERO);
        degree = 0;
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T set(T oth) {
        this.data = oth.data.clone();
        this.degree = oth.degree;
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T shiftLeft(int offset) {
        if (offset == 0)
            return self;
        if (offset > degree)
            return toZero();

        System.arraycopy(data, offset, data, 0, degree - offset + 1);
        Arrays.fill(data, degree - offset + 1, degree + 1, ZERO);
        degree = degree - offset;
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T shiftRight(int offset) {
        if (offset == 0)
            return self;
        int degree = this.degree;
        ensureCapacity(offset + degree);
        System.arraycopy(data, 0, data, offset, degree + 1);
        Arrays.fill(data, 0, offset, ZERO);
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T truncate(int newDegree) {
        if (newDegree >= degree)
            return self;
        Arrays.fill(data, newDegree + 1, degree + 1, ZERO);
        degree = newDegree;
        fixDegree();
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T reverse() {
        ArraysUtil.reverse(data, 0, degree + 1);
        fixDegree();
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T primitivePart() {
        BigInteger content = content();
        if (lc().signum() < 0)
            content = content.negate();
        if (content.isMinusOne())
            return negate();
        return primitivePart0(content);
    }

    @Override
    public T primitivePartSameSign() {
        return primitivePart0(content());
    }

    private T primitivePart0(BigInteger content) {
        if (content.isOne())
            return self;
        for (int i = degree; i >= 0; --i)
            data[i] = data[i].divide(content);
        return self;
    }

    /** {@inheritDoc} */
    @Override
    public final T increment() {
        return add(createOne());
    }

    /** {@inheritDoc} */
    @Override
    public final T decrement() {
        return subtract(createOne());
    }

    /** {@inheritDoc} */
    @Override
    public final T createMonomial(int degree) {return createMonomial(ONE, degree);}

    /**
     * Factory
     *
     * @param data the data
     * @return polynomial
     */
    public abstract T createFromArray(BigInteger[] data);

    /**
     * Creates constant polynomial with specified value
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    public abstract T createConstant(BigInteger val);

    /**
     * Creates monomial {@code coefficient * x^degree}
     *
     * @param coefficient monomial coefficient
     * @param degree      monomial degree
     * @return {@code coefficient * x^degree}
     */
    public abstract T createMonomial(BigInteger coefficient, int degree);

    /** {@inheritDoc} */
    @Override
    public final T createZero() {return createConstant(ZERO);}

    /** {@inheritDoc} */
    @Override
    public final T createOne() {return createConstant(ONE);}

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    public abstract BigInteger evaluate(BigInteger point);

    /**
     * Adds {@code coefficient*x^exponent} to {@code this}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code this + coefficient*x^exponent}
     */
    public abstract T addMonomial(BigInteger coefficient, int exponent);

    /**
     * Adds {@code oth * factor} to {@code this}
     *
     * @param oth    the polynomial
     * @param factor the factor
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    public abstract T addMul(T oth, BigInteger factor);

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this}
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @return {@code this - factor * x^exponent * oth}
     */
    public abstract T subtract(T oth, BigInteger factor, int exponent);

    /** {@inheritDoc} */
    @Override
    public final T multiply(long factor) {return multiply(BigInteger.valueOf(factor));}

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public abstract T multiply(BigInteger factor);

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
        String s = toStringForCopy(data, 0, degree + 1);
        return "create(" + s.substring(1, s.length() - 1) + ")";
    }

    private static String toStringForCopy(BigInteger[] a, int from, int to) {
        if (a == null)
            return "null";
        int iMax = to - 1;
        if (iMax == -1)
            return "[]";

        StringBuilder b = new StringBuilder();
        b.append('[');
        for (int i = from; ; i++) {
            b.append(toStringForCopy(a[i]));
            if (i == iMax)
                return b.append(']').toString();
            b.append(", ");
        }
    }

    private static String toStringForCopy(BigInteger b) {
        return "new BigInteger(\"" + b + "\")";
    }

    public BigInteger[] getDataReferenceUnsafe() {
        return data;
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
                    result[i + j] = result[i + j].add(c.multiply(b[bFrom + j]));
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
                result[i - gFrom] = f[fFrom].multiply(g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            BigInteger[] result = new BigInteger[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = g[gFrom].multiply(f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            BigInteger[] result = new BigInteger[3];
            //both a and b are linear
            result[0] = f[fFrom].multiply(g[gFrom]);
            result[1] = f[fFrom].multiply(g[gFrom + 1]).add(f[fFrom + 1].multiply(g[gFrom]));
            result[2] = f[fFrom + 1].multiply(g[gFrom + 1]);
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
                result[i + split] = result[i + split].add(f1g[i]);
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
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid].add(f[i]);

        //g0 + g1
        BigInteger[] g0_plus_g1 = new BigInteger[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        Arrays.fill(g0_plus_g1, gMid - gFrom, g0_plus_g1.length, ZERO);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = g0_plus_g1[i - gMid].add(g[i]);

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
            mid[i] = mid[i].subtract(f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i].subtract(f1g1[i]);


        int oldLen = f0g0.length;
        BigInteger[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        Arrays.fill(result, oldLen, result.length, ZERO);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split].add(mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split].add(f1g1[i]);

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
                    result[i + j] = result[i + j].add(c.multiply(data[from + j]));
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
            return new BigInteger[]{f[fFrom].multiply(f[fFrom])};
        if (fTo - fFrom == 2) {
            BigInteger[] result = new BigInteger[3];
            result[0] = f[fFrom].multiply(f[fFrom]);
            result[1] = BigInteger.TWO.multiply(f[fFrom]).multiply(f[fFrom + 1]);
            result[2] = f[fFrom + 1].multiply(f[fFrom + 1]);
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
            f0_plus_f1[i - fMid] = f0_plus_f1[i - fMid].add(f[i]);

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
            mid[i] = mid[i].subtract(f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = mid[i].subtract(f1g1[i]);


        int oldLen = f0g0.length;
        BigInteger[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        Arrays.fill(result, oldLen, result.length, ZERO);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = result[i + split].add(mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = result[i + 2 * split].add(f1g1[i]);

        return result;
    }
}
