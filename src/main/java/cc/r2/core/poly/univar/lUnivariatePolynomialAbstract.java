package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.MachineArithmetic;
import cc.r2.core.util.ArraysUtil;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.Arrays;
import java.util.function.LongFunction;
import java.util.stream.LongStream;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Univariate polynomials over machine integers.
 */
abstract class lUnivariatePolynomialAbstract<lPoly extends lUnivariatePolynomialAbstract<lPoly>> implements IUnivariatePolynomial<lPoly> {
    private static final long serialVersionUID = 1L;
    /** array of coefficients { x^0, x^1, ... , x^degree } */
    long[] data;
    /** points to the last non zero element in the data array */
    int degree;
    /** casted self **/
    @SuppressWarnings("unchecked")
    private final lPoly self = (lPoly) this;

    /** {@inheritDoc} */
    @Override
    public final int degree() {return degree;}

    /**
     * Returns the i-th coefficient of this poly (coefficient before x^i)
     */
    public final long get(int i) { return i > degree ? 0 : data[i];}

    /**
     * Sets i-th element of this poly with the specified value
     */
    public final lPoly set(int i, long el) {
        el = valueOf(el);
        if (el == 0) {
            if (i > degree)
                return self;
            data[i] = el;
            fixDegree();
            return self;
        }
        ensureCapacity(i);
        data[i] = el;
        fixDegree();
        return self;
    }

    /**
     * Sets hte leading coefficient of this poly with specified value
     */
    public final lPoly setLC(long lc) {
        return set(degree, lc);
    }

    /** {@inheritDoc} */
    @Override
    public final int firstNonZeroCoefficientPosition() {
        if (isZero()) return -1;
        int i = 0;
        while (data[i] == 0) ++i;
        return i;
    }

    /**
     * Returns the leading coefficient of this poly
     *
     * @return leading coefficient
     */
    public final long lc() {return data[degree];}

    @Override
    public final lPoly lcAsPoly() {return createConstant(lc());}

    @Override
    public final lPoly ccAsPoly() {return createConstant(cc());}

    @Override
    public lPoly getAsPoly(int i) {return createConstant(get(i));}

    /**
     * Returns the constant coefficient of this poly
     *
     * @return constant coefficient
     */
    public final long cc() {return data[0];}

    @Override
    public final void ensureInternalCapacity(int desiredCapacity) {
        if (data.length < desiredCapacity)
            // the rest will be filled with zeros
            data = Arrays.copyOf(data, desiredCapacity);
    }

    /**
     * Ensures that the capacity of internal storage is enough for storing polynomial of the {@code desiredDegree}.
     * The degree of {@code this} is set to {@code desiredDegree} if the latter is greater than the former.
     *
     * @param desiredDegree desired degree
     */
    final void ensureCapacity(int desiredDegree) {
        if (degree < desiredDegree)
            degree = desiredDegree;

        if (data.length < (desiredDegree + 1))
            data = Arrays.copyOf(data, desiredDegree + 1);
    }

    /**
     * Removes zeroes from the end of {@code data} and adjusts the degree
     */
    final void fixDegree() {
        int i = degree;
        while (i >= 0 && data[i] == 0) --i;
        if (i < 0) i = 0;

        if (i != degree) {
            degree = i;
            // unnecessary clearing
            // Arrays.fill(data, degree + 1, data.length, 0);
        }
    }

    /**
     * Converts this to a polynomial over BigIntegers
     *
     * @return polynomial over BigIntegers
     */
    public abstract UnivariatePolynomial<BigInteger> toBigPoly();

    /**
     * Creates new poly with the specified coefficients (over the same domain)
     *
     * @param data the data
     * @return polynomial
     */
    public abstract lPoly createFromArray(long[] data);

    @Override
    public final lPoly createMonomial(int degree) {return createMonomial(1L, degree);}

    /**
     * Creates linear polynomial of form {@code cc + x * lc}  (with the same coefficient domain)
     *
     * @param cc the  constant coefficient
     * @param lc the  leading coefficient
     * @return {@code cc + x * lc}
     */
    public final lPoly createLinear(long cc, long lc) {
        return createFromArray(new long[]{cc, lc});
    }

    /**
     * Creates monomial {@code coefficient * x^degree} (with the same coefficient domain)
     *
     * @param coefficient monomial coefficient
     * @param degree      monomial degree
     * @return {@code coefficient * x^degree}
     */
    public abstract lPoly createMonomial(long coefficient, int degree);

    /**
     * Creates constant polynomial with specified value (with the same coefficient domain)
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    public final lPoly createConstant(long val) {
        return createFromArray(new long[]{val});
    }

    @Override
    public final lPoly createZero() {return createConstant(0);}

    @Override
    public final lPoly createOne() {return createConstant(1);}

    @Override
    public final boolean isZeroAt(int i) {return i >= data.length || data[i] == 0;}

    @Override
    public final lPoly setZero(int i) {
        if (i < data.length)
            data[i] = 0;
        return self;
    }

    @Override
    public final lPoly setFrom(int indexInThis, lPoly poly, int indexInPoly) {
        ensureCapacity(indexInThis);
        data[indexInThis] = poly.get(indexInPoly);
        fixDegree();
        return self;
    }

    @Override
    public final boolean isZero() {return data[degree] == 0;}

    @Override
    public final boolean isOne() {return degree == 0 && data[0] == 1;}

    @Override
    public final boolean isMonic() {return lc() == 1;}

    @Override
    public final boolean isUnitCC() {return cc() == 1;}

    @Override
    public final boolean isConstant() {return degree == 0;}

    @Override
    public final boolean isMonomial() {
        for (int i = degree - 1; i >= 0; --i)
            if (data[i] != 0)
                return false;
        return true;
    }

    @Override
    public final int signumOfLC() {
        return Long.signum(lc());
    }

    /**
     * Returns L1 norm of this polynomial, i.e. sum of abs coefficients
     *
     * @return L1 norm of {@code this}
     */
    public final double norm1() {
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
    public final double norm2() {
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
    public final double normMax() {return (double) maxAbsCoefficient();}

    /**
     * Returns max coefficient (by absolute value) of this poly
     *
     * @return max coefficient (by absolute value)
     */
    public final long maxAbsCoefficient() {
        long max = Math.abs(data[0]);
        for (int i = 1; i <= degree; ++i)
            max = Math.max(Math.abs(data[i]), max);
        return max;
    }

    @Override
    public final lPoly toZero() {
        Arrays.fill(data, 0, degree + 1, 0);
        degree = 0;
        return self;
    }

    @Override
    public final lPoly set(lPoly oth) {
        this.data = oth.data.clone();
        this.degree = oth.degree;
        assert data.length > 0;
        return self;
    }

    @Override
    public final lPoly setAndDestroy(lPoly oth) {
        this.data = oth.data;
        oth.data = null; // destroy
        this.degree = oth.degree;
        assert data.length > 0;
        return self;
    }

    @Override
    public final lPoly setDomainFrom(lPoly poly) {
        return clone();
    }

    @Override
    public final lPoly shiftLeft(int offset) {
        if (offset == 0)
            return self;
        if (offset > degree)
            return toZero();

        System.arraycopy(data, offset, data, 0, degree - offset + 1);
        Arrays.fill(data, degree - offset + 1, degree + 1, 0);
        degree = degree - offset;
        return self;
    }

    @Override
    public final lPoly shiftRight(int offset) {
        if (offset == 0)
            return self;
        int degree = this.degree;
        ensureCapacity(offset + degree);
        System.arraycopy(data, 0, data, offset, degree + 1);
        Arrays.fill(data, 0, offset, 0);
        return self;
    }

    @Override
    public final lPoly truncate(int newDegree) {
        if (newDegree >= degree)
            return self;
        Arrays.fill(data, newDegree + 1, degree + 1, 0);
        degree = newDegree;
        fixDegree();
        return self;
    }

    @Override
    public final lPoly reverse() {
        ArraysUtil.reverse(data, 0, degree + 1);
        fixDegree();
        return self;
    }

    /**
     * Returns the content of this poly (gcd of its coefficients)
     *
     * @return polynomial content
     */
    public final long content() {
        if (degree == 0)
            return data[0];
        return MachineArithmetic.gcd(data, 0, degree + 1);
    }

    @Override
    public final lPoly contentAsPoly() {
        return createConstant(content());
    }

    @Override
    public final lPoly primitivePart() {
        if (isZero())
            return self;
        long content = content();
        if (lc() < 0)
            content = -content;
        if (content == -1)
            return negate();
        return primitivePart0(content);
    }

    @Override
    public final lPoly primitivePartSameSign() {
        return primitivePart0(content());
    }

    private lPoly primitivePart0(long content) {
        if (content == 1)
            return self;
        Magic magic = magicSigned(content);
        for (int i = degree; i >= 0; --i)
            data[i] = divideSignedFast(data[i], magic);
        return self;
    }

    /** addition in the coefficient domain **/
    abstract long add(long a, long b);

    /** subtraction in the coefficient domain **/
    abstract long subtract(long a, long b);

    /** multiplication in the coefficient domain **/
    abstract long multiply(long a, long b);

    /** negation in the coefficient domain **/
    abstract long negate(long a);

    /** convert long to element of this coefficient domain **/
    abstract long valueOf(long a);

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    public final long evaluate(long point) {
        if (point == 0)
            return cc();

        point = valueOf(point);
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = add(multiply(res, point), data[i]);
        return res;
    }

    @Override
    public final lPoly composition(lPoly value) {
        if (value.isOne())
            return this.clone();
        if (value.isZero())
            return ccAsPoly();

        lPoly result = createZero();
        for (int i = degree; i >= 0; --i)
            result = result.multiply(value).add(data[i]);
        return result;
    }

    /**
     * Shifts variables x -> x + value and returns the result (new instance)
     *
     * @param value shift amount
     * @return a copy of this with x -> x + value
     */
    public final lPoly shift(long value) {
        return composition(createLinear(value, 1));
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor};
     *
     * @param factor the factor
     * @return {@code this}
     */
    public abstract lPoly monic(long factor);

    @Override
    public lPoly monicWithLC(lPoly other) {
        return monic(other.lc());
    }

    /**
     * Add constant to this.
     *
     * @param val some number
     * @return this + val
     */
    public final lPoly add(long val) {
        data[0] = add(data[0], valueOf(val));
        fixDegree();
        return self;
    }

    @Override
    public final lPoly decrement() {
        return subtract(createOne());
    }

    @Override
    public final lPoly increment() {
        return add(createOne());
    }

    @Override
    public final lPoly add(lPoly oth) {
        if (oth.isZero())
            return self;
        if (isZero())
            return set(oth);

        assertSameDomainWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = add(data[i], oth.data[i]);
        fixDegree();
        return self;
    }

    /**
     * Adds {@code coefficient*x^exponent} to {@code this}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code this + coefficient*x^exponent}
     */
    public final lPoly addMonomial(long coefficient, int exponent) {
        if (coefficient == 0)
            return self;

        ensureCapacity(exponent);
        data[exponent] = add(data[exponent], valueOf(coefficient));
        fixDegree();
        return self;
    }

    /**
     * Adds {@code oth * factor} to {@code this}
     *
     * @param oth    the polynomial
     * @param factor the factor
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    public final lPoly addMul(lPoly oth, long factor) {
        if (oth.isZero())
            return self;

        factor = valueOf(factor);
        if (factor == 0)
            return self;

        assertSameDomainWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = add(data[i], multiply(factor, oth.data[i]));
        fixDegree();
        return self;
    }

    @Override
    public final lPoly subtract(lPoly oth) {
        if (oth.isZero())
            return self;
        if (isZero())
            return set(oth).negate();

        assertSameDomainWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = subtract(data[i], oth.data[i]);
        fixDegree();
        return self;
    }

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this}
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @return {@code this - factor * x^exponent * oth}
     */
    public final lPoly subtract(lPoly oth, long factor, int exponent) {
        if (oth.isZero())
            return self;

        factor = valueOf(factor);
        if (factor == 0)
            return self;

        assertSameDomainWith(oth);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtract(data[i], multiply(factor, oth.data[i - exponent]));

        fixDegree();
        return self;
    }

    @Override
    public final lPoly negate() {
        for (int i = degree; i >= 0; --i)
            data[i] = negate(data[i]);
        return self;
    }

    @Override
    public lPoly multiplyByLC(lPoly other) {
        return multiply(other.lc());
    }

    @Override
    public final lPoly multiply(long factor) {
        factor = valueOf(factor);
        if (factor == 1)
            return self;

        if (factor == 0)
            return toZero();

        for (int i = degree; i >= 0; --i)
            data[i] = multiply(data[i], factor);
        return self;
    }

    @Override
    public abstract lPoly clone();

    /** convert this long[] data to BigInteger[] */
    final BigInteger[] dataToBigIntegers() {
        BigInteger[] bData = new BigInteger[degree + 1];
        for (int i = degree; i >= 0; --i)
            bData[i] = BigInteger.valueOf(data[i]);
        return bData;
    }

    /** internal API >>> direct unsafe access to internal storage */
    public final long[] getDataReferenceUnsafe() {
        return data;
    }

    @Override
    public final int compareTo(lPoly o) {
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
    public String toString() {
        if (isZero())
            return "0";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.length; i++) {
            String str = termToString(i);
            if (sb.length() == 0 && str.startsWith("+"))
                str = str.substring(1);
            sb.append(str);
        }
        return sb.toString();
    }

    private String termToString(int i) {
        long el = data[i];
        if (el == 0)
            return "";
        String coefficient;
        boolean needSeparator;
        if (el == 1) {
            coefficient = "";
            needSeparator = false;
        } else if (el == -1) {
            coefficient = "-";
            needSeparator = false;
        } else {
            coefficient = Long.toString(el);
            needSeparator = true;
        }

        if (!coefficient.startsWith("-") && !coefficient.startsWith("+"))
            coefficient = "+" + coefficient;

        String m;
        if (i == 0)
            m = "";
        else
            m = ((needSeparator ? "*" : "") + "x" + (i == 1 ? "" : "^" + i));
        if (m.isEmpty())
            if (el == 1 || el == -1)
                coefficient = coefficient + "1";

        return coefficient + m;
    }

    public String toStringForCopy() {
        String s = ArraysUtil.toString(data, 0, degree + 1);
        return "of(" + s.substring(1, s.length() - 1) + ")";
    }

    /**
     * Returns a sequential {@code Stream} with coefficients of this as its source.
     *
     * @return a sequential {@code Stream} over the coefficients in this polynomial
     */
    public final LongStream stream() {
        return Arrays.stream(data, 0, degree + 1);
    }

    /**
     * Applies transformation function to this and returns the result. This method is equivalent of
     * {@code stream().mapToObj(mapper).collect(new PolynomialCollector<>(domain))}.
     *
     * @param domain domain of the new polynomial
     * @param mapper function that maps coefficients of this to coefficients of the result
     * @param <T>    result elements type
     * @return a new polynomial with the coefficients obtained from this by applying {@code mapper}
     */
    public final <T> UnivariatePolynomial<T> mapElements(Domain<T> domain, LongFunction<T> mapper) {
        return stream()
                .mapToObj(mapper)
                .collect(new UnivariatePolynomial.PolynomialCollector<T>(domain));
    }

    @Override
    public final boolean equals(Object obj) {
        if (obj.getClass() != this.getClass())
            return false;
        lUnivariatePolynomialAbstract oth = (lUnivariatePolynomialAbstract) obj;
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
            int elementHash = (int) (element ^ (element >>> 32));
            result = 31 * result + elementHash;
        }
        return result;
    }

    /* =========================== Multiplication with safe arithmetic =========================== */

    /** switch to classical multiplication */
    static final long KARATSUBA_THRESHOLD = 2048L;
    /** when use Karatsuba fast multiplication */
    static final long
            MUL_CLASSICAL_THRESHOLD = 256L * 256L,
            MUL_MOD_CLASSICAL_THRESHOLD = 128L * 128L;

    /** switch algorithms */
    final long[] multiplyUnsafe0(lPoly oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** switch algorithms */
    final long[] multiplySafe0(lPoly oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalSafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaSafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }


    /** switch algorithms */
    final long[] squareUnsafe0() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalUnsafe(data, 0, degree + 1);
        else
            return squareKaratsubaUnsafe(data, 0, degree + 1);
    }

    /** switch algorithms */
    final long[] squareSafe0() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalSafe(data, 0, degree + 1);
        else
            return squareKaratsubaSafe(data, 0, degree + 1);
    }

    /* =========================== Exact multiplication with safe arithmetics =========================== */

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
    final long[] multiplyClassicalSafe(final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        multiplyClassicalSafe(result, a, aFrom, aTo, b, bFrom, bTo);
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
    void multiplyClassicalSafe(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassicalSafe(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = add(result[i + j], multiply(c, b[bFrom + j]));
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
    long[] multiplyKaratsubaSafe(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = multiply(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = multiply(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = multiply(f[fFrom], g[gFrom]);
            result[1] = add(multiply(f[fFrom], g[gFrom + 1]), multiply(f[fFrom + 1], g[gFrom]));
            result[2] = multiply(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyClassicalSafe(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyKaratsubaSafe(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyKaratsubaSafe(f, fFrom, fFrom + split, g, gFrom, gTo);
            long[] f1g = multiplyKaratsubaSafe(f, fFrom + split, fTo, g, gFrom, gTo);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = add(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyKaratsubaSafe(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyKaratsubaSafe(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = add(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = add(g0_plus_g1[i - gMid], g[i]);

        long[] mid = multiplyKaratsubaSafe(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = subtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = subtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = add(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = add(result[i + 2 * split], f1g1[i]);

        return result;
    }

    long[] squareClassicalSafe(long[] a, int from, int to) {
        long[] x = new long[(to - from) * 2 - 1];
        squareClassicalSafe(x, a, from, to);
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
    void squareClassicalSafe(final long[] result, long[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[i + j] = add(result[i + j], multiply(c, data[from + j]));
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
    long[] squareKaratsubaSafe(final long[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return new long[0];
        if (fTo - fFrom == 1)
            return new long[]{multiply(f[fFrom], f[fFrom])};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = multiply(f[fFrom], f[fFrom]);
            result[1] = multiply(multiply(valueOf(2L), f[fFrom]), f[fFrom + 1]);
            result[2] = multiply(f[fFrom + 1], f[fFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (fTo - fFrom) < KARATSUBA_THRESHOLD)
            return squareClassicalSafe(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareKaratsubaSafe(f, fFrom, fMid);
        long[] f1g1 = squareKaratsubaSafe(f, fMid, fTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = add(f0_plus_f1[i - fMid], f[i]);

        long[] mid = squareKaratsubaSafe(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = subtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = subtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = add(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = add(result[i + 2 * split], f1g1[i]);

        return result;
    }

    /* =========================== Exact multiplication with unsafe arithmetics =========================== */

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
            return squareClassicalUnsafe(f, fFrom, fTo);


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
