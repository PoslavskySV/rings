package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;

import java.util.Arrays;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.ZERO;

/**
 * Univariate polynomial over Z.
 * All operations (except where it is specifically stated) changes the content of this.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class bMutablePolynomialZ extends bMutablePolynomialAbstract<bMutablePolynomialZ> implements IMutablePolynomialZ<bMutablePolynomialZ> {

    static BigInteger safeAdd(BigInteger a, BigInteger b) {return a.add(b);}

    static BigInteger safeMultiply(BigInteger a, BigInteger b) {return a.multiply(b);}

    static BigInteger safeMultiply(BigInteger a, int i) {return a.multiply(BigInteger.valueOf(i));}

    static BigInteger safeSubtract(BigInteger a, BigInteger b) {return a.subtract(b);}

    /** main constructor */
    private bMutablePolynomialZ(BigInteger[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    /** copy constructor */
    private bMutablePolynomialZ(BigInteger[] data, int degree) {
        this.data = data;
        this.degree = degree;
    }

    /**
     * Creates Z[x] polynomial from the specified coefficients
     *
     * @param data coefficients
     * @return Z[x] polynomial
     */
    static bMutablePolynomialZ create(BigInteger... data) {
        return new bMutablePolynomialZ(data);
    }

    /**
     * Creates Z[x] polynomial from the specified coefficients
     *
     * @param data coefficients
     * @return Z[x] polynomial
     */
    static bMutablePolynomialZ create(long... data) {
        return new bMutablePolynomialZ(Arrays.stream(data).mapToObj(BigInteger::valueOf).toArray(BigInteger[]::new));
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    static bMutablePolynomialZ monomial(BigInteger coefficient, int exponent) {
        BigInteger[] data = new BigInteger[exponent + 1];
        data[exponent] = coefficient;
        return new bMutablePolynomialZ(data);
    }

    @Override
    public bMutablePolynomialZ createMonomial(BigInteger coefficient, int degree) {
        return monomial(coefficient, degree);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @return polynomial 0
     */
    static bMutablePolynomialZ zero() {
        return new bMutablePolynomialZ(new BigInteger[]{ZERO}, 0);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @return polynomial 1
     */
    static bMutablePolynomialZ one() {
        return new bMutablePolynomialZ(new BigInteger[]{ONE}, 0);
    }

//    /**
//     * Returns Mignotte's bound (sqrt(n+1) * 2^n max |this|)
//     *
//     * @return Mignotte's bound
//     */
//    double mignotteBound() {
//        return Math.pow(2.0, degree) * norm2();
//    }

    /**
     * Reduces polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @param copy    whether to copy the internal data or reduce inplace
     * @return Zp[x] polynomial from this
     */
    bMutablePolynomialMod modulus(BigInteger modulus, boolean copy) {
        return bMutablePolynomialMod.create(modulus, copy ? data.clone() : data);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @return Zp[x] polynomial from this
     */
    bMutablePolynomialMod modulus(BigInteger modulus) {
        return modulus(modulus, true);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns Zp[x] result.
     *
     * @param modulus the modulus
     * @return Zp[x] polynomial from this
     */
    bMutablePolynomialMod modulus(long modulus) {
        return modulus(BigInteger.valueOf(modulus));
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    bMutablePolynomialZ divideOrNull(BigInteger factor) {
        if (factor.isZero())
            throw new ArithmeticException("Divide by zero");
        if (factor.isOne())
            return this;
        for (int i = degree; i >= 0; --i) {
            BigInteger[] qr = data[i].divideAndRemainder(factor);
            if (qr[1].isZero())
                return null;
            data[i] = qr[0];
        }
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ divideOrNullByLC(bMutablePolynomialZ other) {
        return divideOrNull(other.lc());
    }

    @Override
    public bMutablePolynomialZ createFromArray(BigInteger[] data) {
        return new bMutablePolynomialZ(data);
    }

    @Override
    public bMutablePolynomialZ createConstant(BigInteger val) {
        return new bMutablePolynomialZ(new BigInteger[]{val}, 0);
    }

    @Override
    public bMutablePolynomialZ[] arrayNewInstance(int length) {
        return new bMutablePolynomialZ[length];
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ derivative() {
        if (isConstant())
            return createZero();
        BigInteger[] newData = new BigInteger[degree];
        for (int i = degree; i > 0; --i)
            newData[i - 1] = safeMultiply(data[i], i);
        return createFromArray(newData);
    }

    /** {@inheritDoc} */
    @Override
    public BigInteger evaluate(BigInteger point) {
        if (point.isZero())
            return cc();
        BigInteger res = 0;
        for (int i = degree; i >= 0; --i)
            res = safeAdd(safeMultiply(res, point), data[i]);
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
    BigInteger evaluateAtRational(BigInteger num, BigInteger den) {
        if (num.isZero())
            return cc();
        BigInteger res = 0;
        for (int i = degree; i >= 0; --i) {
            BigInteger x = safeMultiply(res, num);
            BigInteger[] qr = x.divideAndRemainder(den);
            if (!qr[1].isZero())
                throw new IllegalArgumentException("The answer is not integer");
            res = safeAdd(qr[0], data[i]);
        }
        return res;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ add(bMutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = safeAdd(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ addMonomial(BigInteger coefficient, int exponent) {
        if (coefficient.isZero())
            return this;

        ensureCapacity(exponent);
        data[exponent] = safeAdd(data[exponent], coefficient);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ addMul(bMutablePolynomialZ oth, BigInteger factor) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = safeAdd(data[i], safeMultiply(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ subtract(bMutablePolynomialZ oth) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = safeSubtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ subtract(bMutablePolynomialZ oth, BigInteger factor, int exponent) {
        if (oth.isZero())
            return this;

        ensureCapacity(oth.degree + exponent);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = safeSubtract(data[i], safeMultiply(factor, oth.data[i - exponent]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ negate() {
        for (int i = degree; i >= 0; --i)
            data[i] = data[i].negate();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ multiply(BigInteger factor) {
        for (int i = degree; i >= 0; --i)
            data[i] = safeMultiply(factor, data[i]);
        return this;
    }

    @Override
    public bMutablePolynomialZ clone() {
        return new bMutablePolynomialZ(data.clone(), degree);
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ multiply(bMutablePolynomialZ oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            BigInteger factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        data = multiplyUnsafe0(oth);

        degree += oth.degree;
        fixDegree();
        return this;
    }

    /** switch algorithms */
    private BigInteger[] multiplyUnsafe0(bMutablePolynomialZ oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialZ square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        data = squareUnsafe();
        degree += degree;
        fixDegree();
        return this;
    }

    /** switch algorithms */
    private BigInteger[] squareUnsafe() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalUnsafe(data, 0, degree + 1);
        else
            return squareKaratsubaUnsafe(data, 0, degree + 1);
    }

    public static bMutablePolynomialZ parse(String string) {
        String[] terms = string.split("\\+");
        Object[][] data = new Object[terms.length][];
        int maxDegree = 0;
        for (int i = 0; i < terms.length; i++) {
            data[i] = parseTerm(terms[i]);
            maxDegree = Integer.max(maxDegree, (int) data[i][1]);
        }
        BigInteger[] bdata = new BigInteger[maxDegree + 1];
        Arrays.fill(bdata, ZERO);
        for (Object[] datum : data)
            bdata[(int) datum[1]] = (BigInteger) datum[0];

        return create(bdata);
    }

    private static Object[] parseTerm(String term) {
        if (!term.contains("x"))
            return new Object[]{new BigInteger(term.trim()), 0};

        term = term.replace("x", "").trim();
        if (term.isEmpty())
            return new Object[]{ONE, 1};
        String[] be = term.split("\\^");
        String cf = be[0].trim();
        if (cf.isEmpty()) cf = "1";
        if (be.length == 1)
            return new Object[]{new BigInteger(cf), 0};
        else
            return new Object[]{new BigInteger(cf), Integer.parseInt(be[1].trim())};
    }
}
