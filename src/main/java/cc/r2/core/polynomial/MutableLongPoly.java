package cc.r2.core.polynomial;

import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;

import static cc.r2.core.polynomial.LongArithmetics.*;

/**
 * Intermediate structure for {@code long[]} univariate polynomial. All operations such as {@link #multiply(MutableLongPoly)},
 * {@link #add(MutableLongPoly)} etc. modify the content of this. All arithmetics operations are safe from long
 * overflow, so that in the case of overflow the {@link ArithmeticException} is always thrown. Thus, the results of all
 * math operations are garantied to be correct.
 */
final class MutableLongPoly implements Comparable<MutableLongPoly> {
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    long[] data;
    /** points to the last non zero element in the data array */
    int degree;

    /** copy constructor */
    private MutableLongPoly(long[] data, int degree) {
        this.data = data;
        this.degree = degree;
    }

    /** private constructor */
    private MutableLongPoly(long[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    /**
     * Creates poly with specified coefficients
     *
     * @param data coefficients
     * @return the polynomial
     */
    static MutableLongPoly create(long... data) {
        return new MutableLongPoly(data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    static MutableLongPoly createMonomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new MutableLongPoly(data);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @return polynomial 1
     */
    static MutableLongPoly one() {
        return new MutableLongPoly(new long[]{1}, 0);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @return polynomial 0
     */
    static MutableLongPoly zero() {
        return new MutableLongPoly(new long[]{0}, 0);
    }

    /**
     * Returns {@code true} if this is zero
     *
     * @return whether {@code this} is zero
     */
    boolean isZero() {return data[degree] == 0;}

    /**
     * Returns {@code true} if this is one
     *
     * @return whether {@code this} is one
     */
    boolean isOne() {return degree == 0 && data[0] == 1;}

    /**
     * Returns {@code true} if this polynomial is monic
     *
     * @return whether {@code this} is monic
     */
    boolean isMonic() {return lc() == 1;}

    /**
     * Returns {@code true} if this polynomial has only constant term
     *
     * @return whether {@code this} is constant
     */
    boolean isConstant() {return degree == 0;}

    /**
     * Returns {@code true} if this polynomial is has only one monomial term
     *
     * @return whether {@code this} has the form {@code c*x^i} (one term)
     */
    boolean isMonomial() {
        for (int i = degree - 1; i >= 0; --i)
            if (data[i] != 0)
                return false;
        return true;
    }

    /**
     * Returns L2 norm of this polynomial, i.e. a square root of a sum of coefficient squares.
     *
     * @return L2 norm of {@code this}
     */
    double norm() {
        double norm = 0;
        for (int i = 0; i <= degree; ++i)
            norm += ((double) data[i]) * data[i];
        return Math.ceil(Math.sqrt(norm));
    }

    /**
     * Returns L1 norm of this polynomial, i.e. the maximal absolute coefficient of {@code this}
     *
     * @return L1 norm of {@code this}
     */
    double norm1() {
        double norm = Math.abs(data[0]);
        for (int i = 1; i <= degree; ++i)
            norm = Math.max((double) Math.abs(data[i]), norm);
        return norm;
    }

    /**
     * Returns the leading coefficient of the poly
     *
     * @return leading coefficient
     */
    long lc() {return data[degree];}

    /**
     * Returns the constant coefficient of the poly
     *
     * @return constant coefficient
     */
    long cc() {return data[0];}

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    long evaluate(long point) {
        if (point == 0)
            return cc();
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = LongArithmetics.add(LongArithmetics.multiply(res, point), data[i]);
        return res;
    }

    /**
     * Evaluates this poly at a given {@code point} modulo {@code modulus} (via Horner method).
     *
     * @param point   {@code point}
     * @param modulus the modulus
     * @return value at {@code point} modulo {@code modulus}
     */
    long evaluate(long point, long modulus) {
        if (point == 0)
            return mod(cc(), modulus);
        long res = 0;
        for (int i = degree; i >= 0; --i)
            res = LongArithmetics.addMod(LongArithmetics.multiplyMod(res, point, modulus), data[i], modulus);
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
            long x = LongArithmetics.multiply(res, num);
            if (x % den != 0)
                throw new IllegalArgumentException("The answer is not integer");
            res = LongArithmetics.add(x / den, data[i]);
        }
        return res;
    }

    /**
     * Ensures that the capcity of internal storage is enough for storing polynomial of the {@code desiredDegree}.
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
     * Removes the zero tail of the polynoamial and sets the degree to its actual value.
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
     * Reduces coefficients of this polynomial modulo {@code modulus}
     *
     * @param modulus the modulus
     * @return {@code this} reduced modulo {@code modulus}
     */
    MutableLongPoly modulus(long modulus) {
        for (int i = degree; i >= 0; --i)
            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    /**
     * Reduces coefficients of this polynomial modulo {@code modulus} using the symmetric representation for Zp
     * elements, i.e. so that each coefficient is in the bound {@code -modulus / 2 <= c <= modulus / 2}.
     *
     * @param modulus the modulus
     * @return {@code this} reduced modulo {@code modulus} in the symmetric representation for Zp elements
     */
    MutableLongPoly symModulus(long modulus) {
        for (int i = degree; i >= 0; --i)
            data[i] = symMod(data[i], modulus);
        fixDegree();
        return this;
    }

    /**
     * Divides this polynomial by a {@code factor}.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor}
     * @throws ArithmeticException if some of the elements can't be exactly divided by the {@code factor}
     */
    MutableLongPoly divide(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        for (int i = degree; i >= 0; --i) {
            assert data[i] % factor == 0 : "not divisible";
            data[i] /= factor;
        }
        return this;
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} if some of the elements can't be exactly
     * divided by the {@code factor}
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    MutableLongPoly divideOrNull(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        for (int i = degree; i >= 0; --i) {
            if (data[i] % factor != 0)
                return null;
            data[i] /= factor;
        }
        return this;
    }

    /**
     * Reduces this polynomial modulo {@code modulus} and returns its monic part (that is {@code this} multiplied by
     * its inversed leading coefficient).
     *
     * @param modulus the modulus
     * @return {@code this} as a monic polynomial
     */
    MutableLongPoly monic(long modulus) {
        if (data[degree] == 0) // isZero()
            return this;
        if (degree == 0) {
            data[0] = 1;
            return this;
        }
        return multiply(modInverse(lc(), modulus), modulus);
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor, modulus)} ).
     *
     * @param factor  the factor
     * @param modulus the modulus
     * @return {@code this }
     */
    MutableLongPoly monic(long factor, long modulus) {
        return multiply(LongArithmetics.multiply(mod(factor, modulus), modInverse(lc(), modulus)), modulus);
    }

    /** fill content with zeroes */
    private MutableLongPoly toZero() {
        degree = 0;
        Arrays.fill(data, 0);
        return this;
    }

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    MutableLongPoly multiply(long factor) {
        if (factor == 1)
            return this;
        if (factor == 0)
            return toZero();
        for (int i = degree; i >= 0; --i)
            data[i] = LongArithmetics.multiply(data[i], factor);
        return this;
    }

    /**
     * Raises {@code this} by the {@code factor} modulo {@code modulus}
     *
     * @param factor  the factor
     * @param modulus the modulus
     * @return {@code} this multiplied by the {@code factor} modulo {@code modulus}
     */
    MutableLongPoly multiply(long factor, long modulus) {
        factor = mod(factor, modulus);
        if (factor == 1)
            return modulus(modulus);
        if (factor == 0)
            return toZero();
        for (int i = degree; i >= 0; --i)
            data[i] = mod(LongArithmetics.multiply(mod(data[i], modulus), factor), modulus);
        return this;
    }

    /**
     * Adds {@code oth} to {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this + oth}
     */
    MutableLongPoly add(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.add(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /**
     * Adds {@code oth} to {@code this} modulo {@code modulus}.
     *
     * @param oth     the polynomial
     * @param modulus the modulus
     * @return {@code this + oth} modulo {@code modulus}
     */
    MutableLongPoly add(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], oth.data[i], modulus);
//        for (int i = degree; i > oth.degree; --i)
//            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    /**
     * Adds {@code oth * factor} to {@code this} modulo {@code modulus}.
     *
     * @param oth     the polynomial
     * @param factor  the factor
     * @param modulus the modulus
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    MutableLongPoly addMul(MutableLongPoly oth, long factor, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], multiplyMod(factor, oth.data[i], modulus), modulus);
//        for (int i = degree; i > oth.degree; --i)
//            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    /**
     * Subtracts {@code oth} from {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this - oth}
     */
    MutableLongPoly subtract(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.subtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /**
     * Subtracts {@code oth} from {@code this} modulo {@code modulus}.
     *
     * @param oth     the polynomial
     * @param modulus the modulus
     * @return {@code this - oth} modulo {@code modulus}
     */
    MutableLongPoly subtract(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = subtractMod(data[i], oth.data[i], modulus);
//        for (int i = degree; i > oth.degree; --i)
//            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this}.
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @return {@code this - factor * x^exponent * oth}
     */
    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = LongArithmetics.subtract(data[i], LongArithmetics.multiply(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this} modulo {@code modulus}.
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @param modulus  the modulus
     * @return {@code this - factor * x^exponent * oth} modulo {@code modulus}
     */
    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree + exponent);
        for (int i = exponent - 1; i >= 0; --i)
            data[i] = mod(data[i], modulus);

        factor = mod(factor, modulus);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtractMod(data[i], LongArithmetics.multiply(factor, mod(oth.data[i - exponent], modulus)), modulus);
//        for (int i = degree, to = oth.degree + exponent; i > to; --i)
//            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    /** when use Karatsuba fast multiplication */
    private static final long
            MUL_CLASSICAL_THRESHOLD = 256L * 256L,
            MUL_MOD_CLASSICAL_THRESHOLD = 128L * 128L;

    /**
     * Multiply {@code this} and {@code oth}. Classical multiplication and
     * Karatsuba used depending on the degree of input polynomials.
     *
     * @param oth this polynomial
     * @return {@code this * oth}
     */
    MutableLongPoly multiply(MutableLongPoly oth) {
        if (this == oth)
            return square();
        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (isZero())
            return this;
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            multiply(factor);
            return this;
        }
        return multiply0(oth);
    }

    /**
     * Square {@code this}. Classical multiplication and
     * Karatsuba used depending on the degree.
     *
     * @return {@code this * this}
     */
    MutableLongPoly square() {
        if (degree == 0) {
            data[0] = LongArithmetics.multiply(data[0], data[0]);
            return this;
        }
        return square0();
    }

    /** switch algorithms */
    private MutableLongPoly multiply0(MutableLongPoly oth) {
        if (LongArithmetics.multiply(degree + 1, oth.degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassical(oth);
        else
            return multiplyKaratsuba(oth);
    }

    /** switch algorithms */
    private MutableLongPoly square0() {
        if (LongArithmetics.multiply(degree + 1, degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassical();
        else
            return squareKaratsuba();
    }

    /** Karatsuba multiplication */
    private MutableLongPoly multiplyKaratsuba(MutableLongPoly oth) {
        data = multiplyKaratsuba(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        degree = degree + oth.degree;
        fixDegree();
        return this;
    }

    /** classical multiplication */
    private MutableLongPoly multiplyClassical(MutableLongPoly oth) {
        data = multiplyClassical(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        degree = degree + oth.degree;
        assert data[degree] != 0;
        return this;
    }

    /** Karatsuba square */
    private MutableLongPoly squareKaratsuba() {
        data = squareKaratsuba(data, 0, degree + 1);
        degree = 2 * degree;
        fixDegree();
        return this;
    }

    /** classical square */
    private MutableLongPoly squareClassical() {
        data = squareClassical(data, 0, degree + 1);
        degree = 2 * degree;
        assert data[degree] != 0;
        return this;
    }

    /**
     * Multiply {@code this} and {@code oth} modulo {@code modulus}. Classical multiplication and
     * Karatsuba used depending on the degree of input polynomials.
     *
     * @param oth the polynomial
     * @return {@code this * oth}
     */
    MutableLongPoly multiply(MutableLongPoly oth, long modulus) {
        if (this == oth)
            return square(modulus);
        if (oth.degree == 0)
            return multiply(oth.data[0], modulus);
        if (isZero())
            return this;
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            multiply(factor, modulus);
            return this;
        }

        double rBound = norm1() * oth.norm1() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE) {
            //don't use modular arithmetic at all
            return multiply0(oth).modulus(modulus);
        }

        if (LongArithmetics.multiply(degree + 1, oth.degree + 1) <= MUL_MOD_CLASSICAL_THRESHOLD)
            return multiplyModClassical(oth, modulus);
        else
            return multiplyModKaratsuba(oth, modulus);
    }

    /**
     * Square {@code this} and {@code oth} modulo {@code modulus}. Classical multiplication and
     * Karatsuba used depending on the degree.
     *
     * @return {@code this * this}
     */
    MutableLongPoly square(long modulus) {
        if (degree == 0) {
            data[0] = LongArithmetics.multiplyMod(data[0], data[0], modulus);
            return this;
        }

        double norm1 = norm1();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE) {
            //don't use modular arithmetic at all
            return square().modulus(modulus);
        }

        if (LongArithmetics.multiply(degree + 1, degree + 1) <= MUL_MOD_CLASSICAL_THRESHOLD)
            return squareModClassical(modulus);
        else
            return squareModKaratsuba(modulus);
    }

    /** Karatsuba multiplication */
    private MutableLongPoly multiplyModKaratsuba(MutableLongPoly oth, long modulus) {
        data = multiplyModKaratsuba(data, 0, degree + 1, oth.data, 0, oth.degree + 1, modulus);
        degree = degree + oth.degree;
        fixDegree();
        return this;
    }

    /** classical multiplication */
    private MutableLongPoly multiplyModClassical(MutableLongPoly oth, long modulus) {
        data = multiplyModClassical(data, 0, degree + 1, oth.data, 0, oth.degree + 1, modulus);
        degree = degree + oth.degree;
        fixDegree();
        return this;
    }

    /** Karatsuba multiplication */
    private MutableLongPoly squareModKaratsuba(long modulus) {
        data = squareModKaratsuba(data, 0, degree + 1, modulus);
        degree = 2 * degree;
        fixDegree();
        return this;
    }

    /** classical multiplication */
    private MutableLongPoly squareModClassical(long modulus) {
        data = squareModClassical(data, 0, degree + 1, modulus);
        degree = 2 * degree;
        fixDegree();
        return this;
    }

    /**
     * Multiplies {@code this} by the monomial with the specified {@code exponent}.
     *
     * @param exponent monomial exponent
     * @return {@code this * x^exponent}
     */
    MutableLongPoly multiplyMonomial(int exponent) {
        if (exponent == 0)
            return this;
        int oldDegree = degree;
        ensureCapacity(oldDegree * exponent);
        System.arraycopy(data, 0, data, exponent, oldDegree + 1);
        Arrays.fill(data, 0, exponent, 0);
        return this;
    }

    @Override
    public int compareTo(MutableLongPoly o) {
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

    /**
     * Deep copy
     */
    @Override
    public MutableLongPoly clone() {
        return new MutableLongPoly(data.clone(), degree);
    }

    @Override
    public String toString() {
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
    public boolean equals(Object obj) {
        if (obj.getClass() != this.getClass())
            return false;
        MutableLongPoly oth = (MutableLongPoly) obj;
        if (degree != oth.degree)
            return false;
        for (int i = 0; i <= degree; ++i)
            if (data[i] != oth.data[i])
                return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = 1;
        for (int i = degree; i >= 0; --i) {
            long element = data[i];
            int elementHash = (int) (element^(element >>> 32));
            result = 31 * result + elementHash;
        }
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
    static void multiplyClassical(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassical(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = LongArithmetics.add(result[i + j], LongArithmetics.multiply(c, b[bFrom + j]));
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
     * Classical n*m multiplication algorithm
     *
     * @param result  where to write the result
     * @param a       the first multiplier
     * @param aFrom   begin in a
     * @param aTo     end in a
     * @param b       the second multiplier
     * @param bFrom   begin in b
     * @param bTo     end in b
     * @param modulus the modulus
     */
    static void multiplyModClassical(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo, long modulus) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyModClassical(result, b, bFrom, bTo, a, aFrom, aTo, modulus);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = mod(a[aFrom + i], modulus);
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = LongArithmetics.addMod(result[i + j], LongArithmetics.multiplyMod(c, b[bFrom + j], modulus), modulus);
        }
    }

    /**
     * Classical n*m multiplication algorithm
     *
     * @param a       the first multiplier
     * @param aFrom   begin in a
     * @param aTo     end in a
     * @param b       the second multiplier
     * @param bFrom   begin in b
     * @param bTo     end in b
     * @param modulus the modulus
     * @return the result
     */
    static long[] multiplyModClassical(final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo, final long modulus) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        multiplyModClassical(result, a, aFrom, aTo, b, bFrom, bTo, modulus);
        return result;
    }

    /** switch to classical multiplication */
    private static final long KARATSUBA_THRESHOLD = 1024L;

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
                result[i - gFrom] = LongArithmetics.multiply(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = LongArithmetics.multiply(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = LongArithmetics.multiply(f[fFrom], g[gFrom]);
            result[1] = LongArithmetics.add(LongArithmetics.multiply(f[fFrom], g[gFrom + 1]), LongArithmetics.multiply(f[fFrom + 1], g[gFrom]));
            result[2] = LongArithmetics.multiply(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if (LongArithmetics.multiply(fTo - fFrom, gTo - gFrom) < KARATSUBA_THRESHOLD)
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
                result[i + split] = LongArithmetics.add(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyKaratsuba(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyKaratsuba(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.add(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = LongArithmetics.add(g0_plus_g1[i - gMid], g[i]);

        long[] mid = multiplyKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.add(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.add(result[i + 2 * split], f1g1[i]);

        return result;
    }

    /**
     * Karatsuba multiplication
     *
     * @param f       the first multiplier
     * @param g       the second multiplier
     * @param fFrom   begin in f
     * @param fTo     end in f
     * @param gFrom   begin in g
     * @param gTo     end in g
     * @param modulus the modulus
     * @return the result
     */
    static long[] multiplyModKaratsuba(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo,
            long modulus) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = LongArithmetics.multiplyMod(f[fFrom], g[i], modulus);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = LongArithmetics.multiplyMod(g[gFrom], f[i], modulus);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = LongArithmetics.multiplyMod(f[fFrom], g[gFrom], modulus);
            result[1] = LongArithmetics.addMod(LongArithmetics.multiplyMod(f[fFrom], g[gFrom + 1], modulus), LongArithmetics.multiplyMod(f[fFrom + 1], g[gFrom], modulus), modulus);
            result[2] = LongArithmetics.multiplyMod(f[fFrom + 1], g[gFrom + 1], modulus);
            return result;
        }
        //switch to classical
        if (LongArithmetics.multiply(fTo - fFrom, gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyModClassical(g, gFrom, gTo, f, fFrom, fTo, modulus);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyModKaratsuba(g, gFrom, gTo, f, fFrom, fTo, modulus);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyModKaratsuba(f, fFrom, fFrom + split, g, gFrom, gTo, modulus);
            long[] f1g = multiplyModKaratsuba(f, fFrom + split, fTo, g, gFrom, gTo, modulus);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = LongArithmetics.addMod(result[i + split], f1g[i], modulus);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyModKaratsuba(f, fFrom, fMid, g, gFrom, gMid, modulus);
        long[] f1g1 = multiplyModKaratsuba(f, fMid, fTo, g, gMid, gTo, modulus);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.addMod(f0_plus_f1[i - fMid], f[i], modulus);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = LongArithmetics.addMod(g0_plus_g1[i - gMid], g[i], modulus);

        long[] mid = multiplyModKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length, modulus);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.subtractMod(mid[i], f0g0[i], modulus);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.subtractMod(mid[i], f1g1[i], modulus);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.addMod(result[i + split], mid[i], modulus);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.addMod(result[i + 2 * split], f1g1[i], modulus);

        return result;
    }

    /** classical square */
    static long[] squareClassical(long[] a, int from, int to) {
        long[] x = new long[(to - from) * 2 - 1];
        squareClassical(x, a, from, to);
        return x;
    }

    static long[] squareModClassical(long[] a, int from, int to, long modulus) {
        long[] x = new long[(to - from) * 2 - 1];
        squareModClassical(x, a, from, to, modulus);
        return x;
    }


//     NOT EFFECTIVE METHODS !
//    /**
//     * Square the poly {@code data} using classical algorithm
//     *
//     * @param result result destination
//     * @param data   the data
//     * @param from   data from
//     * @param to     end point in the {@code data}
//     */
//    static void squareClassical(final long[] result, long[] data, int from, int to) {
//        int degree = to - from - 1;
//        for (int i = 0, s = 2 * degree; i <= s; ++i) {
//            int jFrom = Math.max(0, i - degree);
//            int jTo = Math.min(degree, i);
//            int m = jTo - jFrom + 1;
//            int m2 = m >> 1;
//            jTo = jFrom + m2 - 1;
//            long r = 0;
//            for (int j = jFrom; j <= jTo; ++j)
//                r = LongArithmetics.add(r, LongArithmetics.multiply(data[from + j], data[from + i - j]));
//
//            r = LongArithmetics.add(r, r);
//            if ((m&1) == 1)
//                r = LongArithmetics.add(r, LongArithmetics.multiply(data[from + jTo + 1], data[from + jTo + 1]));
//            result[i] = r;
//        }
//    }

//    /**
//     * Square the poly {@code data} using classical algorithm modulo {@code modulus}
//     *
//     * @param result  result destination
//     * @param data    the data
//     * @param from    end point in the {@code data}
//     * @param modulus the modulus
//     */
//    static void squareModClassical(final long[] result, long[] data, int from, int to, long modulus) {
//        int degree = to - from - 1;
//        for (int i = 0, s = 2 * degree; i <= s; ++i) {
//            int jFrom = Math.max(0, i - degree);
//            int jTo = Math.min(degree, i);
//            int m = jTo - jFrom + 1;
//            int m2 = m >> 1;
//            jTo = jFrom + m2 - 1;
//            long r = 0;
//            for (int j = jFrom; j <= jTo; ++j)
//                r = LongArithmetics.addMod(r, LongArithmetics.multiplyMod(data[from + j], data[from + i - j], modulus), modulus);
//
//            r = LongArithmetics.addMod(r, r, modulus);
//            if ((m&1) == 1)
//                r = LongArithmetics.addMod(r, LongArithmetics.multiplyMod(data[from + jTo + 1], data[from + jTo + 1], modulus), modulus);
//            result[i] = r;
//        }
//    }

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
                    result[i + j] = LongArithmetics.add(result[i + j], LongArithmetics.multiply(c, data[from + j]));
        }
    }

    /**
     * Square the poly {@code data} using classical algorithm
     *
     * @param result result destination
     * @param data   the data
     * @param from   data from
     * @param to     end point in the {@code data}
     */
    static void squareModClassical(final long[] result, long[] data, int from, int to, long modulus) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[i + j] = LongArithmetics.addMod(result[i + j], LongArithmetics.multiplyMod(c, data[from + j], modulus), modulus);
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
            return new long[]{LongArithmetics.multiply(f[fFrom], f[fFrom])};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = LongArithmetics.multiply(f[fFrom], f[fFrom]);
            result[1] = LongArithmetics.multiply(2L, f[fFrom], f[fFrom + 1]);
            result[2] = LongArithmetics.multiply(f[fFrom + 1], f[fFrom + 1]);
            return result;
        }
        //switch to classical
        if (LongArithmetics.multiply(fTo - fFrom, fTo - fFrom) < KARATSUBA_THRESHOLD)
            return MutableLongPoly.squareClassical(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareKaratsuba(f, fFrom, fMid);
        long[] f1g1 = squareKaratsuba(f, fMid, fTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.add(f0_plus_f1[i - fMid], f[i]);

        long[] mid = squareKaratsuba(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.add(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.add(result[i + 2 * split], f1g1[i]);

        return result;
    }

    /**
     * Karatsuba squaring
     *
     * @param f     the data
     * @param fFrom begin in f
     * @param fTo   end in f
     * @return the result
     */
    static long[] squareModKaratsuba(final long[] f, final int fFrom, final int fTo, long modulus) {
        if (fFrom >= fTo)
            return new long[0];
        if (fTo - fFrom == 1)
            return new long[]{LongArithmetics.multiplyMod(f[fFrom], f[fFrom], modulus)};
        if (fTo - fFrom == 2) {
            long[] result = new long[3];
            result[0] = LongArithmetics.multiplyMod(f[fFrom], f[fFrom], modulus);
            result[1] = LongArithmetics.multiplyMod(2L, f[fFrom], f[fFrom + 1], modulus);
            result[2] = LongArithmetics.multiplyMod(f[fFrom + 1], f[fFrom + 1], modulus);
            return result;
        }
        //switch to classical
        if (LongArithmetics.multiply(fTo - fFrom, fTo - fFrom) < KARATSUBA_THRESHOLD)
            return MutableLongPoly.squareModClassical(f, fFrom, fTo, modulus);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        long[] f0g0 = squareModKaratsuba(f, fFrom, fMid, modulus);
        long[] f1g1 = squareModKaratsuba(f, fMid, fTo, modulus);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.addMod(f0_plus_f1[i - fMid], f[i], modulus);

        long[] mid = squareModKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, modulus);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.subtractMod(mid[i], f0g0[i], modulus);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.subtractMod(mid[i], f1g1[i], modulus);


        long[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.addMod(result[i + split], mid[i], modulus);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.addMod(result[i + 2 * split], f1g1[i], modulus);

        return result;
    }
}
