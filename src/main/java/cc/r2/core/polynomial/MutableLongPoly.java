package cc.r2.core.polynomial;

import java.util.Arrays;

import static cc.r2.core.number.ArithmeticUtils.modInverse;
import static cc.r2.core.number.ArithmeticUtils.symMod;
import static java.lang.Math.*;

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
        return ceil(sqrt(norm));
    }

    double norm1() {
        double norm = data[0];
        for (int i = 1; i <= degree; ++i)
            norm = max((double) data[i], norm);
        return norm;
    }

    /**
     * Returns the leading coefficient of the poly
     *
     * @return leading coefficient
     */
    long lc() {return data[degree];}

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
        for (int i = 0; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly symModulus(long modulus) {
        for (int i = 0; i <= degree; ++i)
            data[i] = symMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly divide(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        for (int i = 0; i <= degree; ++i) {
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
        for (int i = 0; i <= degree; ++i) {
            if (data[i] % factor != 0)
                return null;
            data[i] /= factor;
        }
        return this;
    }

    MutableLongPoly monic(long modulus) {
        long inverse = modInverse(lc(), modulus);
        for (int i = 0; i <= degree; ++i)
            data[i] = floorMod(floorMod(multiplyExact(floorMod(data[i], modulus), inverse), modulus), modulus);
        return this;
    }

    MutableLongPoly monic(long factor, long modulus) {
        if (factor == 0)
            return toZero();
        if (lc() == 1)
            return multiply(factor, modulus);

        long inverse = modInverse(lc(), modulus);
        for (int i = 0; i <= degree; ++i)
            data[i] = floorMod(
                    multiplyExact(
                            floorMod(multiplyExact(floorMod(data[i], modulus), inverse), modulus),
                            floorMod(factor, modulus)),
                    modulus);
        return this;
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
        for (int i = 0; i <= degree; ++i)
            data[i] = multiplyExact(data[i], factor);
        return this;
    }

    MutableLongPoly multiply(long factor, long modulus) {
        if (factor == 1)
            return modulus(modulus);
        if (factor == 0)
            return toZero();
        for (int i = 0; i <= degree; ++i)
            data[i] = floorMod(multiplyExact(floorMod(data[i], modulus), floorMod(factor, modulus)), modulus);
        return this;
    }

    MutableLongPoly add(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = 0; i <= oth.degree; ++i)
            data[i] = addExact(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    MutableLongPoly add(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        int i = 0;
        for (; i <= oth.degree; ++i)
            data[i] = floorMod(addExact(floorMod(data[i], modulus), floorMod(oth.data[i], modulus)), modulus);
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = 0; i <= oth.degree; ++i)
            data[i] = subtractExact(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree);
        int i = 0;
        for (; i <= oth.degree; ++i)
            data[i] = floorMod(subtractExact(floorMod(data[i], modulus), floorMod(oth.data[i], modulus)), modulus);
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = exponent; i <= oth.degree + exponent; ++i)
            data[i] = subtractExact(data[i], multiplyExact(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent, long modulus) {
        if (oth.isZero())
            return modulus(modulus);

        ensureCapacity(oth.degree + exponent);
        int i = 0;
        for (; i < exponent; ++i)
            data[i] = floorMod(data[i], modulus);

        factor = floorMod(factor, modulus);
        for (; i <= oth.degree + exponent; ++i)
            data[i] = floorMod(subtractExact(
                    floorMod(data[i], modulus),
                    floorMod(multiplyExact(factor, floorMod(oth.data[i - exponent], modulus)), modulus))
                    , modulus);
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly multiply(MutableLongPoly oth) {
        if (oth.degree == 0)
            return multiply(oth.data[0]);

        final long[] newData = new long[degree + oth.degree + 1];
        for (int j = 0; j <= oth.degree; ++j) // <- upper loop over oth (safe redundant floorMod operations)
            for (int i = 0; i <= degree; ++i)
                newData[i + j] = addExact(newData[i + j], multiplyExact(data[i], oth.data[j]));

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
        for (int j = 0; j <= oth.degree; ++j)// <- upper loop over oth (safe redundant floorMod operations)
            for (int i = 0; i <= degree; ++i)
                newData[i + j] = floorMod(addExact(
                        newData[i + j],
                        floorMod(multiplyExact(data[i], floorMod(oth.data[j], modulus)), modulus))
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
