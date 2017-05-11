package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.LongArithmetics;
import cc.r2.core.poly.ModularDomain;
import cc.r2.core.poly.lModularDomain;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lMutablePolynomialZp extends lMutablePolynomialAbstract<lMutablePolynomialZp> {
    final lModularDomain domain;

    /** copy constructor */
    private lMutablePolynomialZp(lModularDomain domain, long[] data, int degree) {
        this.domain = domain;
        this.data = data;
        this.degree = degree;
    }

    /** main constructor */
    private lMutablePolynomialZp(lModularDomain domain, long[] data) {
        this(domain, data, data.length - 1);
        fixDegree();
    }

    private static void checkModulus(long modulus) {
        if (Long.compare(modulus, LongArithmetics.MAX_SUPPORTED_MODULUS) > 0)
            throw new IllegalArgumentException("Too large modulus. Max allowed is " + LongArithmetics.MAX_SUPPORTED_MODULUS);
    }

    /* =========================== Factory methods =========================== */

    /**
     * Creates poly with specified coefficients represented as signed integers reducing them modulo {@code modulus}
     *
     * @param modulus the modulus
     * @param data    coefficients
     * @return the polynomial
     */
    public static lMutablePolynomialZp create(long modulus, long[] data) {
        lModularDomain domain = new lModularDomain(modulus);
        domain.mod(data);
        return new lMutablePolynomialZp(domain, data);
    }

    /**
     * Creates linear polynomial of form {@code cc + x * lc}
     *
     * @param cc      the  constant coefficient
     * @param lc      the  leading coefficient
     * @param modulus the modulus
     * @return {@code cc + x * lc}
     */
    public static lMutablePolynomialZp linear(long cc, long lc, long modulus) {
        return create(modulus, new long[]{cc, lc});
    }

    public static lMutablePolynomialZp createUnsafe(long modulus, long[] data) {
        return new lMutablePolynomialZp(new lModularDomain(modulus), data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param modulus     the modulus
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static lMutablePolynomialZp createMonomial(long modulus, long coefficient, int exponent) {
        lModularDomain domain = new lModularDomain(modulus);
        coefficient = domain.mod(coefficient);
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new lMutablePolynomialZp(domain, data);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param modulus the modulus
     * @param value   the value
     * @return constant polynomial
     */
    public static lMutablePolynomialZp constant(long modulus, long value) {
        lModularDomain domain = new lModularDomain(modulus);
        return new lMutablePolynomialZp(domain, new long[]{domain.mod(value)}, 0);
    }

    /**
     * Returns polynomial corresponding to math 0
     *
     * @param modulus the modulus
     * @return polynomial 0
     */
    public static lMutablePolynomialZp zero(long modulus) {
        return constant(modulus, 0L);
    }

    /**
     * Returns polynomial corresponding to math 1
     *
     * @param modulus the modulus
     * @return polynomial 1
     */
    public static lMutablePolynomialZp one(long modulus) {
        return constant(modulus, 1L);
    }

    /** does not copy the data and does not reduce the data with new modulus */
    public lMutablePolynomialZp setModulusUnsafe(long newModulus) {
        return new lMutablePolynomialZp(new lModularDomain(newModulus), data, degree);
    }

    /**
     * Creates new Zp[x] polynomial with specified modulus.
     *
     * @param newModulus the new modulus
     * @return the new Zp[x] polynomial with specified modulus
     */
    public lMutablePolynomialZp setModulus(long newModulus) {
        long[] newData = data.clone();
        lModularDomain newDomain = new lModularDomain(newModulus);
        newDomain.mod(newData);
        return new lMutablePolynomialZp(newDomain, newData);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[x] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     */
    @SuppressWarnings("unchecked")
    public lMutablePolynomialZ normalSymmetricForm() {
        long[] newData = new long[degree + 1];
        for (int i = degree; i >= 0; --i)
            newData[i] = domain.symmetricForm(data[i]);
        return lMutablePolynomialZ.create(newData);
    }

//    @Override
//    public BigInteger modulusAsBigInt() {
//        return BigInteger.valueOfUnsigned(modulus);
//    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this.
     *
     * @param copy whether to copy the internal data
     * @return Z[x] version of this
     */
    @SuppressWarnings("unchecked")
    public lMutablePolynomialZ normalForm(boolean copy) {
        return lMutablePolynomialZ.create(copy ? data.clone() : data);
    }

    @Override
    public lMutablePolynomialZp[] arrayNewInstance(int length) {
        return new lMutablePolynomialZp[length];
    }

    @Override
    public lMutablePolynomialZp[] arrayNewInstance(lMutablePolynomialZp a, lMutablePolynomialZp b) {
        return new lMutablePolynomialZp[]{a, b};
    }

    @Override
    public lMutablePolynomialZp getRange(int from, int to) {
        return new lMutablePolynomialZp(domain, Arrays.copyOfRange(data, from, to));
    }

    @Override
    public void checkCompatible(lMutablePolynomialZp oth) {
        if (domain.modulus != oth.domain.modulus)
            throw new IllegalArgumentException();
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp createFromArray(long[] newData) {
        domain.mod(newData);
        return new lMutablePolynomialZp(domain, newData);
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp createMonomial(long coefficient, int newDegree) {
        long[] newData = new long[newDegree + 1];
        newData[newDegree] = valueOf(coefficient);
        return new lMutablePolynomialZp(domain, newData, newDegree);
    }

    @Override
    public boolean isOverField() {return true;}

    @Override
    public boolean isOverFiniteField() {return true;}

    @Override
    public BigInteger domainCardinality() {
        return BigInteger.valueOf(domain.modulus);
    }

    /*=========================== Main methods ===========================*/

    @Override
    long add(long a, long b) {
        return domain.addMod(a, b);
    }

    @Override
    long subtract(long a, long b) {
        return domain.subtractMod(a, b);
    }

    @Override
    long multiply(long a, long b) {
        return domain.multiplyMod(a, b);
    }

    @Override
    long negate(long a) {
        return domain.negateMod(a);
    }

    @Override
    long valueOf(long a) {
        return domain.mod(a);
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp monic() {
        if (isMonic())
            return this;
        if (isZero())
            return this;
        if (degree == 0) {
            data[0] = 1;
            return this;
        }
        return multiply(domain.reciprocal(lc()));
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public lMutablePolynomialZp monic(long factor) {
        return multiply(multiply(valueOf(factor), domain.reciprocal(lc())));
    }

    @Override
    public lMutablePolynomialZp divideByLC(lMutablePolynomialZp other) {
        return divide(other.lc());
    }

    /**
     * Divide by specified value
     *
     * @param val the value
     * @return {@code this / val}
     */
    public lMutablePolynomialZp divide(long val) {
        return multiply(domain.reciprocal(val));
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp multiply(lMutablePolynomialZp oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        checkCompatible(oth);
        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        double rBound = normMax() * oth.normMax() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = multiplyUnsafe0(oth);
            degree += oth.degree;
            domain.mod(data);
            fixDegree();
        } else {
            data = multiplySafe0(oth);
            degree += oth.degree;
            fixDegree();
        }
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        double norm1 = normMax();
        double rBound = norm1 * norm1 * (degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = squareUnsafe0();
            degree += degree;
            domain.mod(data);
            fixDegree();
        } else {
            data = squareSafe0();
            degree += degree;
            fixDegree();
        }
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public lMutablePolynomialZp derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        if (degree < domain.modulus)
            for (int i = degree; i > 0; --i)
                newData[i - 1] = multiply(data[i], i);
        else {
            int i = degree;
            for (; i >= domain.modulus; --i)
                newData[i - 1] = multiply(data[i], valueOf(i));
            for (; i > 0; --i)
                newData[i - 1] = multiply(data[i], i);
        }
        return new lMutablePolynomialZp(domain, newData);
    }

    @Override
    public gMutablePolynomial<BigInteger> toBigPoly() {
        return gMutablePolynomial.createUnsafe(new ModularDomain(domain.modulus), dataToBigIntegers());
    }

    /** Deep copy */
    @Override
    public lMutablePolynomialZp clone() {
        return new lMutablePolynomialZp(domain, data.clone(), degree);
    }
}
