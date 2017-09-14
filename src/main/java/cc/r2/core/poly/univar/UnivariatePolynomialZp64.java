package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.IntegersZp;
import cc.r2.core.poly.IntegersZp64;
import cc.r2.core.poly.MachineArithmetic;

import java.util.Arrays;

/**
 * Univariate polynomial over Zp domain with modulus in the range of {@code [2, 2^62) } (the last value is
 * specified by {@link MachineArithmetic#MAX_SUPPORTED_MODULUS_BITS}. Fast methods from {@link IntegersZp64} are used
 * to perform all arithmetic operations.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariatePolynomialZp64 extends AUnivariatePolynomial64<UnivariatePolynomialZp64> {
    private static final long serialVersionUID = 1L;

    /** The domain */
    public final IntegersZp64 domain;

    private UnivariatePolynomialZp64(IntegersZp64 domain, long[] data, int degree) {
        this.domain = domain;
        this.data = data;
        this.degree = degree;
        assert data.length > 0;
    }

    private UnivariatePolynomialZp64(IntegersZp64 domain, long[] data) {
        this(domain, data, data.length - 1);
        fixDegree();
    }

    private static void checkModulus(long modulus) {
        if (Long.compare(modulus, MachineArithmetic.MAX_SUPPORTED_MODULUS) > 0)
            throw new IllegalArgumentException("Too large modulus. Modulus should be less than 2^" + MachineArithmetic.MAX_SUPPORTED_MODULUS_BITS);
    }

    /**
     * Parse string into polynomial
     */
    public static UnivariatePolynomialZp64 parse(long modulus, String string) {
        return UnivariatePolynomialZ64.parse(string).modulus(modulus);
    }

    /**
     * Parse string into polynomial
     */
    public static UnivariatePolynomialZp64 parse(IntegersZp64 modulus, String string) {
        return UnivariatePolynomialZ64.parse(string).modulus(modulus);
    }

    /**
     * Creates poly with specified coefficients represented as signed integers reducing them modulo {@code modulus}
     *
     * @param modulus the modulus
     * @param data    coefficients
     * @return the polynomial
     */
    public static UnivariatePolynomialZp64 create(long modulus, long[] data) {
        IntegersZp64 domain = new IntegersZp64(modulus);
        domain.modulus(data);
        return new UnivariatePolynomialZp64(domain, data);
    }

    /**
     * Creates linear polynomial of form {@code cc + x * lc}
     *
     * @param cc      the  constant coefficient
     * @param lc      the  leading coefficient
     * @param modulus the modulus
     * @return {@code cc + x * lc}
     */
    public static UnivariatePolynomialZp64 linear(long cc, long lc, long modulus) {
        return create(modulus, new long[]{cc, lc});
    }

    /** data is not reduced modulo modulus */
    public static UnivariatePolynomialZp64 createUnsafe(long modulus, long[] data) {
        return new UnivariatePolynomialZp64(new IntegersZp64(modulus), data);
    }

    /** data is not reduced modulo modulus */
    public static UnivariatePolynomialZp64 createUnsafe(IntegersZp64 domain, long[] data) {
        return new UnivariatePolynomialZp64(domain, data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param modulus     the modulus
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static UnivariatePolynomialZp64 monomial(long modulus, long coefficient, int exponent) {
        IntegersZp64 domain = new IntegersZp64(modulus);
        coefficient = domain.modulus(coefficient);
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new UnivariatePolynomialZp64(domain, data);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param modulus the modulus
     * @param value   the value
     * @return constant polynomial
     */
    public static UnivariatePolynomialZp64 constant(long modulus, long value) {
        IntegersZp64 domain = new IntegersZp64(modulus);
        return new UnivariatePolynomialZp64(domain, new long[]{domain.modulus(value)}, 0);
    }

    /**
     * Creates zero polynomial
     *
     * @param modulus the modulus
     * @return zero polynomial
     */
    public static UnivariatePolynomialZp64 zero(long modulus) {
        return constant(modulus, 0L);
    }

    /**
     * Creates zero polynomial
     *
     * @param domain the domain
     * @return zero polynomial
     */
    public static UnivariatePolynomialZp64 zero(IntegersZp64 domain) {
        return new UnivariatePolynomialZp64(domain, new long[]{0L}, 0);
    }

    /**
     * Creates unit polynomial
     *
     * @param modulus the modulus
     * @return unit polynomial
     */
    public static UnivariatePolynomialZp64 one(long modulus) {
        return constant(modulus, 1L);
    }

    /**
     * Creates unit polynomial
     *
     * @param domain the domain
     * @return unit polynomial
     */
    public static UnivariatePolynomialZp64 one(IntegersZp64 domain) {
        return new UnivariatePolynomialZp64(domain, new long[]{1L}, 0);
    }

    /** Returns the modulus */
    public long modulus() {
        return domain.modulus;
    }

    /** does not copy the data and does not reduce the data with new modulus */
    public UnivariatePolynomialZp64 setModulusUnsafe(long newModulus) {
        return new UnivariatePolynomialZp64(new IntegersZp64(newModulus), data, degree);
    }

    /**
     * Creates new Zp[x] polynomial by coping the coefficients of this and reducing them modulo new modulus.
     *
     * @param newModulus the new modulus
     * @return the copy of this reduced modulo new modulus
     */
    public UnivariatePolynomialZp64 setModulus(long newModulus) {
        long[] newData = data.clone();
        IntegersZp64 newDomain = new IntegersZp64(newModulus);
        newDomain.modulus(newData);
        return new UnivariatePolynomialZp64(newDomain, newData);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[x] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     */
    @SuppressWarnings("unchecked")
    public UnivariatePolynomialZ64 asPolyZSymmetric() {
        long[] newData = new long[degree + 1];
        for (int i = degree; i >= 0; --i)
            newData[i] = domain.symmetricForm(data[i]);
        return UnivariatePolynomialZ64.create(newData);
    }

    /**
     * Returns Z[x] polynomial formed from the coefficients of this.
     *
     * @param copy whether to copy the internal data
     * @return Z[x] version of this
     */
    @SuppressWarnings("unchecked")
    public UnivariatePolynomialZ64 asPolyZ(boolean copy) {
        return UnivariatePolynomialZ64.create(copy ? data.clone() : data);
    }

    @Override
    public UnivariatePolynomialZp64[] createArray(int length) {
        return new UnivariatePolynomialZp64[length];
    }

    @Override
    public UnivariatePolynomialZp64[] createArray(UnivariatePolynomialZp64 a, UnivariatePolynomialZp64 b) {
        return new UnivariatePolynomialZp64[]{a, b};
    }

    @Override
    public UnivariatePolynomialZp64[][] createArray2d(int length) {
        return new UnivariatePolynomialZp64[length][];
    }

    @Override
    public UnivariatePolynomialZp64[][] createArray2d(int length1, int length2) {
        return new UnivariatePolynomialZp64[length1][length2];
    }

    @Override
    public UnivariatePolynomialZp64 getRange(int from, int to) {
        return new UnivariatePolynomialZp64(domain, Arrays.copyOfRange(data, from, to));
    }

    @Override
    public boolean sameDomainWith(UnivariatePolynomialZp64 oth) {
        return domain.modulus == oth.domain.modulus;
    }

    @Override
    public UnivariatePolynomialZp64 createFromArray(long[] newData) {
        domain.modulus(newData);
        return new UnivariatePolynomialZp64(domain, newData);
    }

    @Override
    public UnivariatePolynomialZp64 createMonomial(long coefficient, int newDegree) {
        long[] newData = new long[newDegree + 1];
        newData[newDegree] = valueOf(coefficient);
        return new UnivariatePolynomialZp64(domain, newData, newDegree);
    }

    @Override
    public boolean isOverField() {return true;}

    @Override
    public boolean isOverFiniteField() {return true;}

    @Override
    public boolean isOverZ() {return false;}

    @Override
    public BigInteger coefficientDomainCardinality() {
        return BigInteger.valueOf(modulus());
    }

    @Override
    public BigInteger coefficientDomainCharacteristics() {
        return BigInteger.valueOf(modulus());
    }

    @Override
    public boolean isOverPerfectPower() {
        return domain.isPerfectPower();
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerBase() {
        return BigInteger.valueOf(domain.perfectPowerBase());
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerExponent() {
        return BigInteger.valueOf(domain.perfectPowerExponent());
    }

    @Override
    long add(long a, long b) {
        return domain.add(a, b);
    }

    @Override
    long subtract(long a, long b) {
        return domain.subtract(a, b);
    }

    @Override
    long multiply(long a, long b) {
        return domain.multiply(a, b);
    }

    @Override
    long negate(long a) {
        return domain.negate(a);
    }

    @Override
    long valueOf(long a) {
        return domain.modulus(a);
    }

    @Override
    public UnivariatePolynomialZp64 monic() {
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
    @Override
    public UnivariatePolynomialZp64 monic(long factor) {
        return multiply(multiply(valueOf(factor), domain.reciprocal(lc())));
    }

    @Override
    public UnivariatePolynomialZp64 divideByLC(UnivariatePolynomialZp64 other) {
        return divide(other.lc());
    }

    /**
     * Divide by specified value
     *
     * @param val the value
     * @return {@code this / val}
     */
    public UnivariatePolynomialZp64 divide(long val) {
        return multiply(domain.reciprocal(val));
    }

    @Override
    public UnivariatePolynomialZp64 multiplyByBigInteger(BigInteger factor) {
        return multiply(factor.mod(BigInteger.valueOf(modulus())).longValueExact());
    }

    @Override
    public UnivariatePolynomialZp64 multiply(UnivariatePolynomialZp64 oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        assertSameDomainWith(oth);
        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            long factor = data[0];
            this.set(oth);
            return multiply(factor);
        }

        double rBound = normMax() * oth.normMax() * Math.max(degree + 1, oth.degree + 1);
        if (rBound < Long.MAX_VALUE) {
            // we can apply fast integer arithmetic and then reduce
            data = multiplyUnsafe0(oth);
            degree += oth.degree;
            domain.modulus(data);
            fixDegree();
        } else {
            data = multiplySafe0(oth);
            degree += oth.degree;
            fixDegree();
        }
        return this;
    }

    @Override
    public UnivariatePolynomialZp64 square() {
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
            domain.modulus(data);
            fixDegree();
        } else {
            data = squareSafe0();
            degree += degree;
            fixDegree();
        }
        return this;
    }

    @Override
    public UnivariatePolynomialZp64 derivative() {
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
        return new UnivariatePolynomialZp64(domain, newData);
    }

    /**
     * {@inheritDoc}.
     * The domain of the result will be exactly those returned by {@code this.domain.asDomain() }
     */
    @Override
    public UnivariatePolynomial<BigInteger> toBigPoly() {
        return UnivariatePolynomial.createUnsafe(new IntegersZp(domain.modulus), dataToBigIntegers());
    }

    @Override
    public UnivariatePolynomialZp64 clone() {
        return new UnivariatePolynomialZp64(domain, data.clone(), degree);
    }

    @Override
    public UnivariatePolynomialZp64 parsePoly(String string) {
        return UnivariatePolynomialZ64.parse(string).modulus(domain);
    }

    @Override
    void multiplyClassicalSafe(long[] result, long[] a, int aFrom, int aTo, long[] b, int bFrom, int bTo) {
        if (domain.modulusFits32)
            multiplyClassicalSafeTrick(result, a, aFrom, aTo, b, bFrom, bTo);
        else
            super.multiplyClassicalSafe(result, a, aFrom, aTo, b, bFrom, bTo);
    }

    void multiplyClassicalSafeNoTrick(long[] result, long[] a, int aFrom, int aTo, long[] b, int bFrom, int bTo) {
        super.multiplyClassicalSafe(result, a, aFrom, aTo, b, bFrom, bTo);
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
    final void multiplyClassicalSafeTrick(final long[] result, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        // this trick is taken from Seyed Mohammad Mahdi Javadi PhD
        // thesis "EFFICIENT ALGORITHMS FOR COMPUTATIONS WITH SPARSE POLYNOMIALS" page 79
        // This trick works only when modulus (so as all elements in arrays) is less
        // than 2^32 (so modulus*modulus fits machine word).

        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassicalSafeTrick(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }

        long p2 = domain.modulus * domain.modulus; // this is safe multiplication
        int
                aDegree = aTo - aFrom - 1,
                bDegree = bTo - bFrom - 1,
                resultDegree = aDegree + bDegree;

        for (int i = 0; i <= resultDegree; ++i) {
            long acc = 0;
            for (int j = Math.max(0, i - bDegree), to = Math.min(i, aDegree); j <= to; ++j) {
                if (acc > 0)
                    acc = acc - p2;
                acc = acc + a[aFrom + j] * b[bFrom + i - j];
            }
            result[i] = domain.modulus(acc);
        }
    }

    @Override
    public String coefficientDomainToString() {
        return domain.toString();
    }
}
