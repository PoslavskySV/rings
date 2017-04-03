package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;

import java.util.Arrays;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.ZERO;

/**
 * Univariate polynomial over Zp.
 * All operations (except where it is specifically stated) changes the content of this.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class bMutablePolynomialMod extends bMutablePolynomialAbstract<bMutablePolynomialMod>
        implements IMutablePolynomialZp<bMutablePolynomialMod> {
    /** the modulus */
    final BigInteger modulus;

    /** copy constructor */
    private bMutablePolynomialMod(BigInteger modulus, BigInteger[] data, int degree) {
        this.modulus = modulus;
        this.data = data;
        this.degree = degree;
    }

    /** main constructor */
    private bMutablePolynomialMod(BigInteger modulus, BigInteger[] data) {
        this.modulus = modulus;
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
    }

    /* =========================== Factory methods =========================== */

    /**
     * Creates poly with specified coefficients represented as signed integers reducing them modulo {@code modulus}
     *
     * @param modulus the modulus
     * @param data    coefficients
     * @return the polynomial
     */
    static bMutablePolynomialMod create(BigInteger modulus, BigInteger[] data) {
        reduce(data, modulus);
        return new bMutablePolynomialMod(modulus, data);
    }

    static bMutablePolynomialMod createUnsafe(BigInteger modulus, BigInteger[] data) {
        return new bMutablePolynomialMod(modulus, data);
    }

    /** reduce data mod modulus **/
    private static void reduce(BigInteger[] data, BigInteger modulus) {
        for (int i = 0; i < data.length; ++i)
            data[i] = data[i].mod(modulus);
    }

    private static BigInteger reduce(BigInteger val, BigInteger modulus) {
        return (val.signum() < 0 || val.compareTo(modulus) >= 0) ? val.mod(modulus) : val;
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param modulus     the modulus
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    static bMutablePolynomialMod createMonomial(BigInteger modulus, BigInteger coefficient, int exponent) {
        BigInteger[] data = new BigInteger[exponent + 1];
        data[exponent] = reduce(coefficient, modulus);
        return new bMutablePolynomialMod(modulus, data);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param modulus the modulus
     * @param value   the value
     * @return constant polynomial
     */
    static bMutablePolynomialMod constant(BigInteger modulus, BigInteger value) {
        return new bMutablePolynomialMod(modulus, new BigInteger[]{reduce(value, modulus)});
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @param modulus the modulus
     * @return polynomial 1
     */
    static bMutablePolynomialMod one(BigInteger modulus) {
        return constant(modulus, ONE);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @param modulus the modulus
     * @return polynomial 0
     */
    static bMutablePolynomialMod zero(BigInteger modulus) {
        return constant(modulus, ZERO);
    }




    /*=========================== Main methods ===========================*/


    /** modulus operation */
    BigInteger mod(BigInteger val) {
        return val.compareTo(modulus) < 0 ? val : val.mod(modulus);
    }


    /** multiplyMod operation */
    BigInteger multiplyMod(BigInteger a, BigInteger b) {
        return mod(a.multiply(b));
    }

    /** addMod operation */
    BigInteger addMod(BigInteger a, BigInteger b) {
        BigInteger r = a.add(b), rm = r.subtract(modulus);
        return rm.signum() >= 0 ? rm : r;
    }

    /** subtractMod operation */
    BigInteger subtractMod(BigInteger a, BigInteger b) {
        BigInteger r = a.subtract(b);
        return r.signum() < 0 ? r.add(modulus) : r;
    }

    /** to symmetric modulus */
    BigInteger symMod(BigInteger value) {
        return value.compareTo(modulus.shiftRight(1)) <= 0 ? value : value.subtract(modulus);
    }

    @Override
    public bMutablePolynomialMod createFromArray(BigInteger[] newData) {
        reduce(newData, modulus);
        bMutablePolynomialMod r = new bMutablePolynomialMod(modulus, newData, newData.length - 1);
        r.fixDegree();
        return r;
    }

    /** Whether this BigInteger polynomial can be converted to machine-precision polynomial (long) */
    public boolean isLong() { return modulus.compareTo(LongArithmetics.b_MAX_SUPPORTED_MODULUS) <= 0;}

    /**
     * Creates constant polynomial with specified value
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    @Override
    public bMutablePolynomialMod createConstant(BigInteger val) {
        return new bMutablePolynomialMod(modulus, new BigInteger[]{reduce(val, modulus)}, 0);
    }

    @Override
    public bMutablePolynomialMod getRange(int from, int to) {
        return new bMutablePolynomialMod(modulus, Arrays.copyOfRange(data, from, to));
    }

    @Override
    public bMutablePolynomialMod createMonomial(BigInteger coefficient, int newDegree) {
        BigInteger[] newData = new BigInteger[newDegree + 1];
        Arrays.fill(newData, ZERO);
        newData[newDegree] = reduce(coefficient, modulus);
        return new bMutablePolynomialMod(modulus, newData, newDegree);
    }

    @Override
    public bMutablePolynomialMod[] arrayNewInstance(int length) {
        return new bMutablePolynomialMod[length];
    }

    @Override
    public bMutablePolynomialMod[] arrayNewInstance(bMutablePolynomialMod a, bMutablePolynomialMod b) {
        return new bMutablePolynomialMod[]{a, b};
    }

    /** does not copy the data and does not reduce the data with new modulus */
    bMutablePolynomialMod setModulusUnsafe(BigInteger newModulus) {
        return new bMutablePolynomialMod(newModulus, data, degree);
    }

    /**
     * Creates new Zp[x] polynomial with specified modulus.
     *
     * @param newModulus the new modulus
     * @return the new Zp[x] polynomial with specified modulus
     */
    bMutablePolynomialMod setModulus(BigInteger newModulus) {
        BigInteger[] newData = data.clone();
        reduce(newData, newModulus);
        return new bMutablePolynomialMod(newModulus, newData);
    }

    @Override
    public BigInteger modulusAsBigInt() {
        return modulus;
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this.
     *
     * @param copy whether to copy the internal data
     * @return Z[x] version of this
     */
    bMutablePolynomialZ normalForm(boolean copy) {
        return bMutablePolynomialZ.create(copy ? data.clone() : data);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[x] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     */
    bMutablePolynomialZ normalSymmetricForm() {
        BigInteger[] newData = new BigInteger[degree + 1];
        for (int i = degree; i >= 0; --i)
            newData[i] = symMod(data[i]);
        return bMutablePolynomialZ.create(newData);
    }

    MutablePolynomialMod toLong() {
        if (modulus.compareTo(LongArithmetics.b_MAX_SUPPORTED_MODULUS) > 0)
            throw new ArithmeticException("modulus out of supported range");
        return MutablePolynomialZ.create(toLong0()).modulus(modulus.longValueExact(), false);
    }

    /** {@inheritDoc} */
    public bMutablePolynomialMod monic() {
        if (data[degree].isZero()) // isZero()
            return this;
        if (degree == 0) {
            data[0] = ONE;
            return this;
        }
        return multiply(lc().modInverse(modulus));
    }

    @Override
    public bMutablePolynomialMod divideByLC(bMutablePolynomialMod other) {
        return multiply(other.lc().modInverse(modulus));
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    bMutablePolynomialMod monic(BigInteger factor) {
        return multiply(multiplyMod(mod(factor), lc().modInverse(modulus)));
    }

    /** {@inheritDoc} */
    @Override
    public BigInteger evaluate(BigInteger point) {
        if (point.isZero())
            return cc();

        point = mod(point);
        BigInteger res = ZERO;
        for (int i = degree; i >= 0; --i)
            res = addMod(multiplyMod(res, point), data[i]);
        return res;
    }

    /** {@inheritDoc} */
    public void checkCompatibleModulus(bMutablePolynomialMod oth) {
        if (!modulus.equals(oth.modulus))
            throw new IllegalArgumentException();
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod add(bMutablePolynomialMod oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth);

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod addMonomial(BigInteger coefficient, int exponent) {
        if (coefficient.isZero())
            return this;

        ensureCapacity(exponent);
        data[exponent] = addMod(data[exponent], mod(coefficient));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod addMul(bMutablePolynomialMod oth, BigInteger factor) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor.isZero())
            return this;

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = addMod(data[i], multiplyMod(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod subtract(bMutablePolynomialMod oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth).negate();

        checkCompatibleModulus(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = subtractMod(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod subtract(bMutablePolynomialMod oth, BigInteger factor, int exponent) {
        if (oth.isZero())
            return this;

        factor = mod(factor);
        if (factor.isZero())
            return this;

        checkCompatibleModulus(oth);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = subtractMod(data[i], multiplyMod(factor, oth.data[i - exponent]));

        fixDegree();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod negate() {
        for (int i = degree; i >= 0; --i)
            if (!data[i].isZero())
                data[i] = modulus.subtract(data[i]);
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod multiply(BigInteger factor) {
        factor = mod(factor);
        if (factor.isOne())
            return this;

        if (factor.isZero())
            return toZero();

        for (int i = degree; i >= 0; --i)
            data[i] = multiplyMod(data[i], factor);
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod derivative() {
        if (isConstant())
            return createZero();
        BigInteger[] newData = new BigInteger[degree];
        int i = degree;
        for (; i > 0; --i)
            newData[i - 1] = multiplyMod(data[i], BigInteger.valueOf(i));
        return createFromArray(newData);
    }

    /**
     * Deep copy
     */
    @Override
    public bMutablePolynomialMod clone() {
        return new bMutablePolynomialMod(modulus, data.clone(), degree);
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod multiply(bMutablePolynomialMod oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        checkCompatibleModulus(oth);
        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            BigInteger factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        // can apply fast integer arithmetic and then reduce
        data = multiplyExact(oth);
        degree += oth.degree;
        reduce(data, modulus);
        fixDegree();

        return this;
    }


    /** switch algorithms */
    private BigInteger[] multiplyExact(bMutablePolynomialMod oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaUnsafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** {@inheritDoc} */
    @Override
    public bMutablePolynomialMod square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        // we apply fast integer arithmetic and then reduce
        data = squareExact();
        degree += degree;
        reduce(data, modulus);
        fixDegree();

        return this;
    }

    /** switch algorithms */
    private BigInteger[] squareExact() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalUnsafe(data, 0, degree + 1);
        else
            return squareKaratsubaUnsafe(data, 0, degree + 1);
    }
}
