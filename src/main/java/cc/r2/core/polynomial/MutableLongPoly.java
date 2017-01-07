package cc.r2.core.polynomial;

import cc.r2.core.util.ArraysUtil;

import java.util.Arrays;

import static cc.r2.core.polynomial.LongArithmetics.*;

/**
 * Intermediate structure for long[] polynomial.
 */
final class MutableLongPoly implements Comparable<MutableLongPoly> {
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

    static MutableLongPoly createMonomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new MutableLongPoly(data);
    }

    static MutableLongPoly one() {
        return new MutableLongPoly(new long[]{1});
    }

    static MutableLongPoly zero() {
        return new MutableLongPoly(new long[]{0});
    }

    boolean isZero() {return data[degree] == 0;}

    boolean isOne() {return degree == 0 && data[0] == 1;}

    boolean isMonic() {return lc() == 1;}

    boolean isConstant() {return degree == 0;}

    boolean isMonomial() {
        for (int i = 0; i < degree; ++i)
            if (data[i] != 0)
                return false;
        return true;
    }

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
     * Evaluates this poly at a give {@code point}
     *
     * @param point {@code point}
     * @return value at {@code point}
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
        if (data[degree] == 0) // isZero()
            return this;
        if (degree == 0) {
            data[0] = 1;
            return this;
        }
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
        for (int j = oth.degree; j >= 0; --j) { // <- upper loop over oth (safe redundant mod operations)
            long c = oth.data[j];
            if (c != 0)
                for (int i = degree; i >= 0; --i)
                    newData[i + j] = LongArithmetics.add(newData[i + j], LongArithmetics.multiply(data[i], c));
        }

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

    MutableLongPoly multiply(MutableLongPoly oth, long modulus, MutableLongPoly cache) {
        if (oth.degree == 0) {
            cache.data = multiply(oth.data[0], modulus).data;
            return this;
        }

        modulus(modulus);
        cache.ensureCapacity(degree + oth.degree + 1);
        cache.degree = degree + oth.degree + 1;

        final long[] newData = cache.data;
        for (int j = oth.degree; j >= 0; --j)// <- upper loop over oth (safe redundant mod operations)
            for (int i = degree; i >= 0; --i) {
                long di = data[i];
                newData[i + j] = mod(LongArithmetics.add(
                        newData[i + j],
                        mod(LongArithmetics.multiply(di, mod(oth.data[j], modulus)), modulus))
                        , modulus);
            }

        data = newData;
        degree = degree + oth.degree;
        return this;
    }

    static MutableLongPoly create0(long... data) {
        if (data.length == 0) data = new long[1];
        return new MutableLongPoly(data);
    }

    static long[] copyOfRange(long[] data, int from, int to) {
        if (from >= to)
            return new long[0];
        return Arrays.copyOfRange(data, from, to);
    }

    static long[] copy(long[] a, int from, int to) {
        if (from >= a.length)
            return new long[0];
        if (to > a.length) return a;
        return Arrays.copyOfRange(a, from, to);
    }

    long[] multiplyKaratsuba4(long[] a, long[] b) {
        if (a.length == 0)
            return new long[0];
        if (b.length == 0)
            return new long[0];

        if (a.length == 1) {
            long[] result = new long[b.length];
            //single element in a
            for (int i = 0; i < b.length; ++i)
                result[i] = a[0] * b[i];
            return result;
        }
        if (b.length == 1) {
            long[] result = new long[a.length];
            //single element in b
            for (int i = 0; i < a.length; ++i)
                result[i] = b[0] * a[i];
            return result;
        }
        if (a.length == 2 && b.length == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = a[0] * b[0];
            result[1] = a[0] * b[1] + a[1] * b[0]; ;
            result[2] = a[1] * b[1];
            return result;
        }


        int imid = Math.max(a.length, b.length) / 2;
//        int aMid = aFrom + (aTo - aFrom) / 2, bMid = bFrom + (bTo - bFrom) / 2;
//        int aMid = aFrom + imid, bMid = bFrom + imid;

//        int nDegree = Math.max(aTo - aFrom, bTo - bFrom);
//        if (Integer.bitCount(nDegree) != 1)
//            nDegree = 1 << (32 - Integer.numberOfLeadingZeros(nDegree));
//
//        assert nDegree % 2 == 0;


        //f0*g0
        long[] f0 = copy(a, 0, imid), f1 = copy(a, imid, a.length);
        long[] g0 = copy(b, 0, imid), g1 = copy(b, imid, b.length);

        long[] f0g0 = multiplyKaratsuba4(f0, g0);
        long[] f1g1 = multiplyKaratsuba4(f1, g1);


        // f0 + f1
        if (f1.length > f0.length) {
            long[] tmp = f0;
            f0 = f1;
            f1 = tmp;
        }
        if (g1.length > g0.length) {
            long[] tmp = g0;
            g0 = g1;
            g1 = tmp;
        }

        for (int i = 0; i < f1.length; ++i)
            f0[i] += f1[i];
        for (int i = 0; i < g1.length; ++i)
            g0[i] += g1[i];

        long[] mid = multiplyKaratsuba4(f0, g0);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; i++)
            mid[i] -= f0g0[i];

        for (int i = 0; i < f1g1.length; i++)
            mid[i] -= f1g1[i];


        long[] result = Arrays.copyOf(f0g0, a.length + b.length - 1);
        for (int i = 0; i < mid.length; i++)
            result[i + imid] = LongArithmetics.add(result[i + imid], mid[i]);
        for (int i = 0; i < f1g1.length; i++)
            result[i + 2 * imid] = LongArithmetics.add(result[i + 2 * imid], f1g1[i]);

        return result;
    }


    MutableLongPoly multiplyKaratsuba4(MutableLongPoly oth) {
        return MutableLongPoly.create(multiplyKaratsuba4(data, oth.data));
    }

    static long[] multiplyNaive(final long[] a, final long[] b, final int aFrom, final int aTo, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        for (int i = 0; i < aTo - aFrom; i++) {
            for (int j = 0; j < bTo - bFrom; j++) {
                result[i + j] += a[aFrom + i] * b[bFrom + j];
            }
        }
        return result;
    }

    static long[] multiplyKaratsuba4a(final long[] a, final long[] b, final int aFrom, final int aTo, final int bFrom, final int bTo) {
        if (aFrom >= aTo)
            return new long[0];
        if (bFrom >= bTo)
            return new long[0];

        if (aTo - aFrom == 1) {
            long[] result = new long[bTo - bFrom];
            //single element in a
            for (int i = bFrom; i < bTo; ++i)
                result[i - bFrom] = a[aFrom] * b[i];
            return result;
        }
        if (bTo - bFrom == 1) {
            long[] result = new long[aTo - aFrom];
            //single element in b
            for (int i = aFrom; i < aTo; ++i)
                result[i - aFrom] = b[bFrom] * a[i];
            return result;
        }
        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = a[aFrom] * b[bFrom];
            result[1] = a[aFrom] * b[aFrom + 1] + a[aFrom + 1] * b[bFrom]; ;
            result[2] = a[aFrom + 1] * b[bFrom + 1];
            return result;
        }
        if ((aTo - aFrom) * (bTo - bFrom) < 1000)
            return multiplyNaive(a, b, aFrom, aTo, bFrom, bTo);


        int imid = (Math.max(aTo - aFrom, bTo - bFrom)) / 2;
//        System.out.println(aTo - aFrom);
//        System.out.println(bTo - bFrom);
//        System.out.println(imid);
//        int aMid = aFrom + (aTo - aFrom) / 2, bMid = bFrom + (bTo - bFrom) / 2;
//        int aMid = aFrom + imid, bMid = bFrom + imid;

//        int nDegree = Math.max(aTo - aFrom, bTo - bFrom);
//        if (Integer.bitCount(nDegree) != 1)
//            nDegree = 1 << (32 - Integer.numberOfLeadingZeros(nDegree));
//        imid = nDegree /2;
//        assert nDegree % 2 == 0;


        //f0*g0
        int aMid = Math.min(aFrom + imid, aTo);
        int bMid = Math.min(bFrom + imid, bTo);

//        long[] f0 = copy(a, 0, imid), f1 = copy(a, imid, a.length);
//        long[] g0 = copy(b, 0, imid), g1 = copy(b, imid, b.length);

        long[] f0g0 = multiplyKaratsuba4a(a, b, aFrom, aMid, bFrom, bMid);
        long[] f1g1 = multiplyKaratsuba4a(a, b, aMid, aTo, bMid, bTo);


        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(aMid - aFrom, aTo - aMid)];
        long[] g0_plus_g1 = new long[Math.max(bMid - bFrom, bTo - bMid)];

        for (int i = aFrom; i < aMid; i++)
            f0_plus_f1[i - aFrom] = a[i];
        for (int i = aMid; i < aTo; ++i)
            f0_plus_f1[i - aMid] += a[i];

        for (int i = bFrom; i < bMid; i++)
            g0_plus_g1[i - bFrom] = b[i];
        for (int i = bMid; i < bTo; ++i)
            g0_plus_g1[i - bMid] += b[i];

        long[] mid = multiplyKaratsuba4a(f0_plus_f1, g0_plus_g1, 0, f0_plus_f1.length, 0, g0_plus_g1.length);
//        System.out.println(Arrays.toString(mid));

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; i++)
            mid[i] -= f0g0[i];

        for (int i = 0; i < f1g1.length; i++)
            mid[i] -= f1g1[i];


        long[] result = Arrays.copyOf(f0g0, (aTo - aFrom) + (bTo - bFrom) - 1);
        for (int i = 0; i < mid.length; i++) {
//            if (mid[i] != 0)
            result[i + imid] += mid[i];
        }
        for (int i = 0; i < f1g1.length; i++) {
//            if(f1g1[i] != 0)
            result[i + 2 * imid] += f1g1[i];
        }

        return result;
    }

    MutableLongPoly multiplyKaratsuba4a(MutableLongPoly oth) {
        return MutableLongPoly.create(multiplyKaratsuba4a(data, oth.data, 0, degree + 1, 0, oth.degree + 1));
    }

    long[] multiplyKaratsuba3(long[] a, long[] b, int aFrom, int aTo, int bFrom, int bTo) {
//
//        aTo = Math.min(aTo, a.length);
//        bTo = Math.min(bTo, b.length);

        if (aFrom >= aTo || bFrom >= bTo)
            return new long[0];

        if (aTo - aFrom == 1) {
            long[] result = new long[bTo - bFrom];
            //single element in a
            for (int i = bFrom; i < bTo; ++i)
                result[i - bFrom] = LongArithmetics.multiply(a[aFrom], b[i]);
            return result;
        }
        if (bTo - bFrom == 1) {
            long[] result = new long[aTo - aFrom];
            //single element in b
            for (int i = aFrom; i < aTo; ++i)
                result[i - aFrom] = LongArithmetics.multiply(b[bFrom], a[i]);
            return result;
        }
        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = LongArithmetics.multiply(a[aFrom], b[bFrom]);
            result[1] = LongArithmetics.add(LongArithmetics.multiply(a[aFrom], b[bFrom + 1]), LongArithmetics.multiply(a[aFrom + 1], b[bFrom]));
            result[2] = LongArithmetics.multiply(a[aFrom + 1], b[bFrom + 1]);
            return result;
        }


        int imid = Math.max(aTo - aFrom, bTo - bFrom) / 2;
//        int aMid = aFrom + (aTo - aFrom) / 2, bMid = bFrom + (bTo - bFrom) / 2;
//        int aMid = aFrom + imid, bMid = bFrom + imid;

//        int nDegree = Math.max(aTo - aFrom, bTo - bFrom);
//        if (Integer.bitCount(nDegree) != 1)
//            nDegree = 1 << (32 - Integer.numberOfLeadingZeros(nDegree));
//
//        assert nDegree % 2 == 0;


        //f0*g0
        long[] f0g0 = multiplyKaratsuba3(a, b, aFrom, aFrom + imid, bFrom, bFrom + imid);
        //f1*g1
        long[] f1g1 = multiplyKaratsuba3(a, b, aFrom + imid, aTo, bFrom + imid, bTo);


        // a <- f0 + f1
        for (int i = aFrom + imid; i < aTo; i++)
            a[i - imid] = LongArithmetics.add(a[i - imid], a[i]);

        // b <- g0 + g1
        for (int i = bFrom + imid; i < bTo; i++)
            b[i - imid] = LongArithmetics.add(b[i - imid], b[i]);

        long[] mid = multiplyKaratsuba3(a, b, aFrom, aFrom + imid, bFrom, bFrom + imid);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; i++)
            mid[i] = LongArithmetics.subtract(mid[i], f0g0[i]);

        for (int i = 0; i < f1g1.length; i++)
            mid[i] = LongArithmetics.subtract(mid[i], f1g1[i]);


        // restore a: a <- f0 + f1
        for (int i = aFrom + imid; i < aTo; i++)
            a[i - imid] = LongArithmetics.subtract(a[i - imid], a[i]);

        for (int i = bFrom + imid; i < bTo; i++)
            b[i - imid] = LongArithmetics.subtract(b[i - imid], b[i]);


        long[] result = Arrays.copyOf(f0g0, (aTo - aFrom - 1) + (bTo - bFrom - 1) + 1);
        for (int i = 0; i < mid.length; i++)
            result[i + imid] = LongArithmetics.add(result[i + imid], mid[i]);
        for (int i = 0; i < f1g1.length; i++)
            result[i + 2 * imid] = LongArithmetics.add(result[i + 2 * imid], f1g1[i]);

        return result;
    }

    MutableLongPoly multiplyKaratsuba3(MutableLongPoly oth) {
        return MutableLongPoly.create(multiplyKaratsuba3(data, oth.data, 0, degree + 1, 0, oth.degree + 1));
    }

    private final long[] cache = new long[100000];

    int multiplyKaratsuba(int recursion, int aDegree, int bDegree, long[] a, long[] b, long[] target, int aFrom, int aTo, int bFrom, int bTo, int targetFrom) {
//        assert aTo - aFrom == bTo - bFrom;

        if (aFrom > (aDegree + 1) || bFrom > (bDegree + 1))
            return 0;

        aTo = Math.min(aTo, aDegree + 1);
        bTo = Math.min(bTo, bDegree + 1);


        if (aTo - aFrom == 1) {
            //single element in a
            for (int i = bFrom; i < bTo; ++i)
                target[i - bFrom + targetFrom] = LongArithmetics.add(target[i - bFrom + targetFrom], LongArithmetics.multiply(a[aFrom], b[i]));
            return bTo - bFrom;
        }
        if (bTo - bFrom == 1) {
            //single element in b
            for (int i = aFrom; i < aTo; ++i)
                target[i - aFrom + targetFrom] = LongArithmetics.add(target[i - aFrom + targetFrom], LongArithmetics.multiply(b[bFrom], a[i]));
            return aTo - aFrom;
        }
        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
            //both a and b are linear
            target[targetFrom] = LongArithmetics.add(target[targetFrom], LongArithmetics.multiply(a[aFrom], b[bFrom]));
            target[targetFrom + 1] = LongArithmetics.add(target[targetFrom + 1],
                    LongArithmetics.add(LongArithmetics.multiply(a[aFrom], b[bFrom + 1]), LongArithmetics.multiply(a[aFrom + 1], b[bFrom])));
            target[targetFrom + 2] = LongArithmetics.add(target[targetFrom + 2], LongArithmetics.multiply(a[aFrom + 1], b[bFrom + 1]));
            return 3;
        }


        int nDegree = Math.max(aTo - aFrom, bTo - bFrom);
        if (Integer.bitCount(nDegree) != 1)
            nDegree = 1 << (32 - Integer.numberOfLeadingZeros(nDegree));

        assert nDegree % 2 == 0;

        //f0*g0
        int f0g0 = multiplyKaratsuba(recursion + 1, aDegree, bDegree, a, b, target, aFrom, aFrom + nDegree / 2, bFrom, bFrom + nDegree / 2, targetFrom);
        //f1*g1
        int f1g1 = multiplyKaratsuba(recursion + 1, aDegree, bDegree, a, b, target, aFrom + nDegree / 2, aFrom + nDegree, bFrom + nDegree / 2, bFrom + nDegree, targetFrom + nDegree);


        // a <- f0 + f1
        for (int i = 0; i < nDegree / 2 && aFrom + nDegree / 2 + i < aTo; i++)
            a[aFrom + i] = LongArithmetics.add(a[aFrom + i], a[aFrom + nDegree / 2 + i]);

        // b <- g0 + g1
        for (int i = 0; i < nDegree / 2 && bFrom + nDegree / 2 + i < bTo; i++)
            b[bFrom + i] = LongArithmetics.add(b[bFrom + i], b[bFrom + nDegree / 2 + i]);

        int cacheFrom = recursion * nDegree;
        Arrays.fill(cache, cacheFrom, cacheFrom + nDegree, 0);
        int mid = multiplyKaratsuba(recursion + 1, aDegree, bDegree, a, b, cache, aFrom, aFrom + nDegree / 2, bFrom, bFrom + nDegree / 2, cacheFrom);

        // restore a: a <- f0 + f1
        for (int i = 0; i < nDegree / 2 && aFrom + nDegree / 2 + i < aTo; i++)
            a[aFrom + i] = LongArithmetics.subtract(a[aFrom + i], a[aFrom + nDegree / 2 + i]);

        // restore b: b <- g0 + g1
        for (int i = 0; i < nDegree / 2 && bFrom + nDegree / 2 + i < bTo; i++)
            b[bFrom + i] = LongArithmetics.subtract(b[bFrom + i], b[bFrom + nDegree / 2 + i]);


        // subtract f0g0
        for (int i = 0; i < f0g0; ++i)
            cache[cacheFrom + i] = LongArithmetics.subtract(cache[cacheFrom + i], target[targetFrom + i]);

        // subtract f1g1
        for (int i = 0; i < f1g1; ++i)
            cache[cacheFrom + i] = LongArithmetics.subtract(cache[cacheFrom + i], target[targetFrom + nDegree + i]);

        int l = Math.max(f0g0, f1g1);
        l = Math.max(l, mid);
        //write (f0+f1)(g0+g1) - f0g0 - f1g1
        for (int i = nDegree / 2; i < nDegree / 2 + l; i++) {
            target[i] = LongArithmetics.add(target[i], cache[cacheFrom + i - nDegree / 2]);
        }
        return f0g0 + f1g1;
    }


    MutableLongPoly multiplyKaratsuba2(MutableLongPoly oth) {
        long[] data = new long[Math.max(degree, oth.degree) * 4];
        multiplyKaratsuba(0, degree, oth.degree, this.data, oth.data, data, 0, degree + 1, 0, oth.degree + 1, 0);
        return MutableLongPoly.create(data);
    }

    MutableLongPoly multiplyKaratsuba(MutableLongPoly oth, long modulus) {
        if (oth.degree * degree <= 10)
            return multiply(oth, modulus);

        //thisF
        int nDegree = Math.max(degree, oth.degree);
        do {
            ++nDegree;
        } while (Integer.bitCount(nDegree) != 1);

        MutableLongPoly f0 = create0(copyOfRange(data, 0, nDegree / 2));
        MutableLongPoly f1 = create0(copyOfRange(data, nDegree / 2, data.length));

        MutableLongPoly g0 = create0(copyOfRange(oth.data, 0, nDegree / 2));
        MutableLongPoly g1 = create0(copyOfRange(oth.data, nDegree / 2, oth.data.length));

        MutableLongPoly f0g0 = f0.clone().multiplyKaratsuba(g0, modulus);
        MutableLongPoly f1g1 = f1.clone().multiplyKaratsuba(g1, modulus);
        MutableLongPoly mid = f0.add(f1, modulus).multiplyKaratsuba(g0.add(g1, modulus), modulus).subtract(f0g0, modulus).subtract(f1g1, modulus);

        long[] newData = Arrays.copyOf(f0g0.data, degree + oth.degree + 1);
        for (int i = 0; i <= mid.degree; i++) {
            if (mid.data[i] != 0)
                newData[i + nDegree / 2] = addMod(newData[i + nDegree / 2], mid.data[i], modulus);
        }

        for (int i = 0; i <= f1g1.degree; i++) {
            if (f1g1.data[i] != 0)
                newData[i + nDegree] = addMod(newData[i + nDegree], f1g1.data[i], modulus);
        }
//        System.out.println(degree + oth.degree + 1);
////        System.out.println(f0g0.degree + 1 + mid.degree + 1 + f1g1.degree + 1);
//        System.arraycopy(mid.data, 0, newData, nDegree / 2, mid.degree + 1);
//        System.arraycopy(f1g1.data, 0, newData, nDegree, f1g1.degree + 1);
        this.data = newData;
        degree = degree + oth.degree;
        return this;
    }

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

    public String toStringForCopy() {
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
}
