package cc.r2.core.polynomial;

import java.util.Arrays;

import static cc.r2.core.polynomial.LongArithmetics.*;

/**
 * Intermediate structure for long[] polynomial.
 */
final class MutableLongPoly {
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    long[] data;
    /** points to the last non zero element in the data array */
    int degree;

    MutableLongPoly(int desiredDegree) {
        data = new long[desiredDegree + 1];
        degree = desiredDegree;
    }

    MutableLongPoly(long[] data, int degree) {
        this.data = data;
        this.degree = degree;
    }

    private MutableLongPoly(long[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    static MutableLongPoly create(long... data) {
        return new MutableLongPoly(data);
    }

    static MutableLongPoly one() {
        return new MutableLongPoly(new long[]{1});
    }

    static MutableLongPoly zero() {
        return new MutableLongPoly(new long[]{0});
    }

    boolean isZero() {return degree == 0 && data[0] == 0;}

    boolean isOne() {return degree == 0 && data[0] == 1;}

    boolean isMonic() {return lc() == 1;}

    boolean isConstant() {return degree == 0;}

    double norm() {
        double norm = 0;
        for (int i = 0; i <= degree; ++i)
            norm += ((double) data[i]) * data[i];
        return Math.ceil(Math.sqrt(norm));
    }

    double norm1() {
        double norm = data[0];
        for (int i = 1; i <= degree; ++i)
            norm = Math.max((double) data[i], norm);
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
     * Evaluates this poly at a give {@code point}
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
     * Evaluates this poly at a give rational point {@code num/den}
     *
     * @param num evaluation point numerator
     * @param den evaluation point denominator
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

    void ensureCapacity(int desiredDegree) {
        if (degree < desiredDegree)
            degree = desiredDegree;

        if (data.length < (desiredDegree + 1))
            data = Arrays.copyOf(data, desiredDegree + 1);
    }

    void fixDegree() {
        int i = degree;
        while (i >= 0 && data[i] == 0) --i;
        if (i < 0) i = 0;

        if (i != degree) {
            degree = i;
            Arrays.fill(data, degree + 1, data.length, 0);
        }
    }

    MutableLongPoly modulus(long modulus) {
        for (int i = degree; i >= 0; --i)
            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly symModulus(long modulus) {
        for (int i = degree; i >= 0; --i)
            data[i] = symMod(data[i], modulus);
        fixDegree();
        return this;
    }

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

    MutableLongPoly monic(long modulus) {
        return multiply(modInverse(lc(), modulus), modulus);
    }

    MutableLongPoly monic(long factor, long modulus) {
        return multiply(LongArithmetics.multiply(mod(factor, modulus), modInverse(lc(), modulus)), modulus);
    }

    private MutableLongPoly toZero() {
        degree = 0;
        Arrays.fill(data, 0);
        return this;
    }

    MutableLongPoly multiply(long factor) {
        if (factor == 1)
            return this;
        if (factor == 0)
            return toZero();
        for (int i = degree; i >= 0; --i)
            data[i] = LongArithmetics.multiply(data[i], factor);
        return this;
    }

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

    MutableLongPoly add(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.add(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    MutableLongPoly add(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], oth.data[i], modulus);
        for (int i = degree; i > oth.degree; --i)
            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = LongArithmetics.subtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = subtractMod(data[i], oth.data[i], modulus);
        for (int i = degree; i > oth.degree; --i)
            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = LongArithmetics.subtract(data[i], LongArithmetics.multiply(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree + exponent);
        for (int i = exponent - 1; i >= 0; --i)
            data[i] = mod(data[i], modulus);

        factor = mod(factor, modulus);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtractMod(data[i], LongArithmetics.multiply(factor, mod(oth.data[i - exponent], modulus)), modulus);
        for (int i = degree, to = oth.degree + exponent; i > to; --i)
            data[i] = mod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly multiply(MutableLongPoly oth) {
        if (oth.degree == 0)
            return multiply(oth.data[0]);

        final long[] newData = new long[degree + oth.degree + 1];
        for (int j = oth.degree; j >= 0; --j) // <- upper loop over oth (safe redundant mod operations)
            for (int i = degree; i >= 0; --i)
                newData[i + j] = LongArithmetics.add(newData[i + j], LongArithmetics.multiply(data[i], oth.data[j]));

        data = newData;
        degree = degree + oth.degree;
        fixDegree();
        return this;
    }

    MutableLongPoly multiply(MutableLongPoly oth, long modulus) {
        if (oth.degree == 0)
            return multiply(oth.data[0], modulus);

        modulus(modulus);
        final long[] newData = new long[degree + oth.degree + 1];
        for (int j = oth.degree; j >= 0; --j)// <- upper loop over oth (safe redundant mod operations)
            for (int i = degree; i >= 0; --i)
                newData[i + j] = mod(LongArithmetics.add(
                        newData[i + j],
                        mod(LongArithmetics.multiply(data[i], mod(oth.data[j], modulus)), modulus))
                        , modulus);

        data = newData;
        degree = degree + oth.degree;
        fixDegree();
        return this;
    }

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
}
