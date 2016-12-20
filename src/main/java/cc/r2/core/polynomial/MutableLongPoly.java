package cc.r2.core.polynomial;

import java.util.Arrays;

import static cc.r2.core.number.ArithmeticUtils.*;
import static java.lang.Math.floorMod;

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

    MutableLongPoly(long... data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    boolean isZero() {return degree == 0 && data[0] == 0;}

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

    MutableLongPoly multiply(long factor) {
        if (factor == 1)
            return this;
        if (factor == 0) {
            degree = 0;
            Arrays.fill(data, 0);
            return this;
        }
        for (int i = 0; i <= degree; ++i) {
            assert safeMultiply(data[i], factor) : "long overflow";
            data[i] *= factor;
        }
        return this;
    }

    MutableLongPoly divide(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        for (int i = 0; i <= degree; ++i) {
            assert data[i] % factor == 0 : "long overflow";
            data[i] /= factor;
        }
        return this;
    }

    MutableLongPoly multiply(long coefficient, long modulus) {
        if (coefficient == 0) {
            degree = 0;
            Arrays.fill(data, 0);
            return this;
        }
        for (int i = 0; i <= degree; ++i) {
            long x = floorMod(data[i], modulus);
            long y = floorMod(coefficient, modulus);
            assert safeMultiply(x, y) : "long overflow";
            data[i] = floorMod(x * y, modulus);
        }
        return this;
    }

    MutableLongPoly add(MutableLongPoly oth) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree);

        for (int i = 0; i <= oth.degree; ++i) {
            assert safeAdd(data[i], oth.data[i]) : "long overflow";
            data[i] += oth.data[i];
        }
        fixDegree();
        return this;
    }

    MutableLongPoly add(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree);

        int i = 0;
        for (; i <= oth.degree; ++i) {
            long x = floorMod(data[i], modulus);
            long y = floorMod(oth.data[i], modulus);
            assert safeAdd(x, y) : "long overflow";
            data[i] = floorMod(x + y, modulus);
        }
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree);

        for (int i = 0; i <= oth.degree; ++i) {
            assert safeSubtract(data[i], oth.data[i]) : "long overflow";
            data[i] -= oth.data[i];
        }
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long modulus) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree);

        int i = 0;
        for (; i <= oth.degree; ++i) {
            long x = floorMod(data[i], modulus);
            long y = floorMod(oth.data[i], modulus);
            assert safeSubtract(x, y) : "long overflow";
            data[i] = floorMod(x - y, modulus);
        }
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree + exponent);

        for (int i = exponent; i <= oth.degree + exponent; ++i) {
            assert safeMultiply(factor, oth.data[i - exponent]) : "long overflow";
            assert safeSubtract(data[i], factor * oth.data[i - exponent]) : "long overflow";
            data[i] -= factor * oth.data[i - exponent];
        }
        fixDegree();
        return this;
    }

    MutableLongPoly subtract(MutableLongPoly oth, long factor, int exponent, long modulus) {
        if (oth.isZero())
            return this;
        ensureCapacity(oth.degree + exponent);

        int i = 0;
        for (; i < exponent; ++i)
            data[i] = floorMod(data[i], modulus);

        factor = floorMod(factor, modulus);
        for (; i <= oth.degree + exponent; ++i) {
            long xThis = floorMod(data[i], modulus);
            long xOther = floorMod(oth.data[i - exponent], modulus);
            assert safeMultiply(xOther, factor) : "long overflow";
            xOther = floorMod(factor * xOther, modulus);
            assert safeSubtract(xThis, xOther) : "long overflow";
            data[i] = floorMod(xThis - xOther, modulus);
        }
        for (; i <= degree; ++i)
            data[i] = floorMod(data[i], modulus);
        fixDegree();
        return this;
    }

    MutableLongPoly multiply(MutableLongPoly oth) {
        if (oth.degree == 0)
            return multiply(oth.data[0]);

        final long[] thisData = Arrays.copyOf(data, degree + 1);
        ensureCapacity(degree * oth.degree);
        for (int i = 0; i <= degree; ++i)
            for (int j = 0; j <= oth.degree; ++j) {
                assert safeMultiply(thisData[i], oth.data[j]) : "long overflow";
                data[i + j] = thisData[i] * oth.data[j];
            }

        fixDegree();
        return this;
    }

    MutableLongPoly multiply(MutableLongPoly oth, long modulus) {
        if (oth.degree == 0)
            return multiply(oth.data[0]);

        modulus(modulus);
        final long[] thisData = Arrays.copyOf(data, degree + 1);

        ensureCapacity(degree * oth.degree);
        for (int i = 0; i <= degree; ++i) // <- upper loop over oth (safe redundant floorMod operations)
            for (int j = 0; j <= oth.degree; ++j) {
                assert safeMultiply(thisData[j], floorMod(oth.data[i], modulus)) : "long overflow";
                data[i + j] = floorMod(thisData[j] * floorMod(oth.data[i], modulus), modulus);
            }

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
            String c = String.valueOf(data[i]);
            if (!c.startsWith("-") && sb.length() != 0)
                sb.append("+");
            sb.append(c);
            if (i != 0)
                sb.append("x^").append(i);
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
        for (int i = 0; i <= degree; i++)
            if (data[i] != oth.data[i])
                return false;
        return true;
    }
}
