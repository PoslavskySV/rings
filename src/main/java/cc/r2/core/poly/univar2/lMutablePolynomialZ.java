package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.IntegersDomain;
import cc.r2.core.poly.LongArithmetics;
import cc.redberry.libdivide4j.FastDivision.Magic;

import java.util.Arrays;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lMutablePolynomialZ extends lMutablePolynomialAbstract<lMutablePolynomialZ> {

    /** main constructor */
    private lMutablePolynomialZ(long[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    /** copy constructor */
    private lMutablePolynomialZ(long[] data, int degree) {
        this.data = data;
        this.degree = degree;
    }

    /**
     * Creates Z[x] polynomial from the specified coefficients
     *
     * @param data coefficients
     * @return Z[x] polynomial
     */
    public static lMutablePolynomialZ create(long... data) {
        return new lMutablePolynomialZ(data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static lMutablePolynomialZ monomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new lMutablePolynomialZ(data);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @return polynomial 0
     */
    public static lMutablePolynomialZ zero() {
        return new lMutablePolynomialZ(new long[]{0}, 0);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @return polynomial 1
     */
    public static lMutablePolynomialZ one() {
        return new lMutablePolynomialZ(new long[]{1}, 0);
    }

    /**
     * Reduces polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @param copy    whether to copy the internal data or reduce inplace
     * @return Zp[x] polynomial from this
     */
    public lMutablePolynomialZp modulus(long modulus, boolean copy) {
        return lMutablePolynomialZp.create(modulus, copy ? data.clone() : data);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @return Zp[x] polynomial from this
     */
    public lMutablePolynomialZp modulus(long modulus) {
        return modulus(modulus, true);
    }

    lMutablePolynomialZp modulusUnsafe(long modulus) {
        return lMutablePolynomialZp.createUnsafe(modulus, data);
    }

    @Override
    public gMutablePolynomial<BigInteger> toBigPoly() {
        return gMutablePolynomial.createUnsafe(IntegersDomain.IntegersDomain, dataToBigIntegers());
    }

    /**
     * Returns Mignotte's bound (sqrt(n+1) * 2^n max |this|)
     *
     * @return Mignotte's bound
     */
    public double mignotteBound() {
        return Math.pow(2.0, degree) * norm2();
    }

    /**
     * Evaluates this poly at a give rational point {@code num/den}
     *
     * @param num point numerator
     * @param den point denominator
     * @return value at {@code num/den}
     * @throws ArithmeticException if the result is not integer
     */
    public long evaluateAtRational(long num, long den) {
        if (num == 0)
            return cc();
        long res = 0;
        Magic magic = magicSigned(den);
        for (int i = degree; i >= 0; --i) {
            long x = multiply(res, num);
            long q = divideSignedFast(x, magic);
            if (q * den != x)
                throw new IllegalArgumentException("The answer is not integer");
            res = add(q, data[i]);
        }
        return res;
    }

    @Override
    public lMutablePolynomialZ getRange(int from, int to) {
        return new lMutablePolynomialZ(Arrays.copyOfRange(data, from, to));
    }

    @Override
    public lMutablePolynomialZ[] arrayNewInstance(int length) {
        return new lMutablePolynomialZ[length];
    }

    @Override
    public void checkCompatible(lMutablePolynomialZ oth) {}

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ createFromArray(long[] data) {
        return new lMutablePolynomialZ(data);
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ createMonomial(long coefficient, int degree) {
        return monomial(coefficient, degree);
    }

    @Override
    public boolean isOverField() {return false;}

    @Override
    public boolean isOverFiniteField() {return false;}

    @Override
    public BigInteger domainCardinality() {return null;}

    @Override
    long add(long a, long b) {return LongArithmetics.safeAdd(a, b);}

    @Override
    long subtract(long a, long b) {return LongArithmetics.safeSubtract(a, b);}

    @Override
    long multiply(long a, long b) {return LongArithmetics.safeMultiply(a, b);}

    @Override
    long negate(long a) {return -a;}

    @Override
    long valueOf(long a) {return a;}

    @Override
    public lMutablePolynomialZ monic() {
        return divideOrNull(lc());
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public lMutablePolynomialZ divideOrNull(long factor) {
        if (factor == 0)
            throw new ArithmeticException("Divide by zero");
        if (factor == 1)
            return this;
        Magic magic = magicSigned(factor);
        for (int i = degree; i >= 0; --i) {
            long l = divideSignedFast(data[i], magic);
            if (l * factor != data[i])
                return null;
            data[i] = l;
        }
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ divideByLC(lMutablePolynomialZ other) {
        return divideOrNull(other.lc());
    }

    public lMutablePolynomialZ multiplyUnsafe(long factor) {
        for (int i = degree; i >= 0; --i)
            data[i] *= factor;
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ multiply(lMutablePolynomialZ oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        double rBound = normMax() * oth.normMax() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE)
            // we can apply fast integer arithmetic
            data = multiplyUnsafe0(oth);
        else
            data = multiplySafe0(oth);

        degree += oth.degree;
        fixDegree();
        return this;
    }

    public lMutablePolynomialZ multiplyUnsafe(lMutablePolynomialZ oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiplyUnsafe(factor);
        }

        data = multiplyUnsafe0(oth);
        degree += oth.degree;
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = normMax();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE)
            // we can apply fast integer arithmetic
            data = squareUnsafe0();
        else
            data = squareSafe0();

        degree += degree;
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZ derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        for (int i = degree; i > 0; --i)
            newData[i - 1] = multiply(data[i], i);
        return createFromArray(newData);
    }

    @Override
    public lMutablePolynomialZ clone() {
        return new lMutablePolynomialZ(data.clone(), degree);
    }
}
