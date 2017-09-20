package cc.redberry.rings.poly.univar;

import cc.redberry.libdivide4j.FastDivision.Magic;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MachineArithmetic;

import java.util.Arrays;

import static cc.redberry.libdivide4j.FastDivision.divideSignedFast;
import static cc.redberry.libdivide4j.FastDivision.magicSigned;

/**
 * Univariate polynomial over machine integers in range [-2^63, 2^63]. <b>NOTE:</b> this class is used in internal
 * routines for performance reasons, for usual polynomials over Z use {@link UnivariatePolynomial} over BigIntegers.
 *
 * <p> Arithmetic operations on instances of this type may cause long overflow in which case a proper {@link
 * ArithmeticException} will be thrown.
 *
 * @since 1.0
 */
public final class UnivariatePolynomialZ64 extends AUnivariatePolynomial64<UnivariatePolynomialZ64> {
    private static final long serialVersionUID = 1L;

    /** main constructor */
    private UnivariatePolynomialZ64(long[] data) {
        this.data = data;
        this.degree = data.length - 1;
        fixDegree();
        assert data.length > 0;
    }

    /** copy constructor */
    private UnivariatePolynomialZ64(long[] data, int degree) {
        this.data = data;
        this.degree = degree;
        assert data.length > 0;
    }

    /**
     * Parse string into polynomial
     */
    public static UnivariatePolynomialZ64 parse(String string) {
        return UnivariatePolynomial.asLongPolyZ(Parser.parse(Rings.Z, Rings.Z, string));
    }

    /**
     * Creates Z[x] polynomial from the specified coefficients
     *
     * @param data coefficients
     * @return Z[x] polynomial
     */
    public static UnivariatePolynomialZ64 create(long... data) {
        return new UnivariatePolynomialZ64(data);
    }

    /**
     * Creates monomial {@code coefficient * x^exponent}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code coefficient * x^exponent}
     */
    public static UnivariatePolynomialZ64 monomial(long coefficient, int exponent) {
        long[] data = new long[exponent + 1];
        data[exponent] = coefficient;
        return new UnivariatePolynomialZ64(data);
    }

    /**
     * Creates zero polynomial
     *
     * @return zero polynomial
     */
    public static UnivariatePolynomialZ64 zero() {
        return new UnivariatePolynomialZ64(new long[]{0}, 0);
    }

    /**
     * Creates unit polynomial
     *
     * @return unit polynomial
     */
    public static UnivariatePolynomialZ64 one() {
        return new UnivariatePolynomialZ64(new long[]{1}, 0);
    }

    /**
     * Returns constant with specified value
     *
     * @return constant with specified value
     */
    public static UnivariatePolynomialZ64 constant(long value) {
        return new UnivariatePolynomialZ64(new long[]{value}, 0);
    }

    /**
     * Reduces this polynomial modulo {@code modulus} and returns the result.
     *
     * @param modulus the modulus
     * @param copy    whether to copy the internal data or reduce inplace (in which case the data of this will be lost)
     * @return this modulo {@code modulus}
     */
    public UnivariatePolynomialZp64 modulus(long modulus, boolean copy) {
        return UnivariatePolynomialZp64.create(modulus, copy ? data.clone() : data);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns the result.
     *
     * @param modulus the modulus
     * @return a copy of this modulo {@code modulus}
     */
    public UnivariatePolynomialZp64 modulus(long modulus) {
        return modulus(modulus, true);
    }

    /**
     * Reduces this polynomial modulo {@code modulus} and returns the result.
     *
     * @param ring the modulus
     * @param copy whether to copy the internal data or reduce inplace (in which case the data of this will be lost)
     * @return this modulo {@code modulus}
     */
    public UnivariatePolynomialZp64 modulus(IntegersZp64 ring, boolean copy) {
        long[] data = copy ? this.data.clone() : this.data;
        for (int i = degree; i >= 0; --i)
            data[i] = ring.modulus(data[i]);
        return UnivariatePolynomialZp64.createUnsafe(ring, data);
    }

    /**
     * Reduces (copied) polynomial modulo {@code modulus} and returns the result.
     *
     * @param ring the modulus
     * @return a copy of this modulo {@code modulus}
     */
    public UnivariatePolynomialZp64 modulus(IntegersZp64 ring) {
        return modulus(ring, true);
    }

    /** internal API */
    UnivariatePolynomialZp64 modulusUnsafe(long modulus) {
        return UnivariatePolynomialZp64.createUnsafe(modulus, data);
    }

    /**
     * {@inheritDoc}. The ring of the result is {@link Rings#Z}
     */
    @Override
    public UnivariatePolynomial<BigInteger> toBigPoly() {
        return UnivariatePolynomial.createUnsafe(Rings.Z, dataToBigIntegers());
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
     * Evaluates this poly at a given rational point {@code num/den}
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
    public UnivariatePolynomialZ64 getRange(int from, int to) {
        return new UnivariatePolynomialZ64(Arrays.copyOfRange(data, from, to));
    }

    @Override
    public UnivariatePolynomialZ64[] createArray(int length) {
        return new UnivariatePolynomialZ64[length];
    }

    @Override
    public UnivariatePolynomialZ64[][] createArray2d(int length) {
        return new UnivariatePolynomialZ64[length][];
    }

    @Override
    public UnivariatePolynomialZ64[][] createArray2d(int length1, int length2) {
        return new UnivariatePolynomialZ64[length1][length2];
    }

    @Override
    public boolean sameCoefficientRingWith(UnivariatePolynomialZ64 oth) {return true;}

    @Override
    public UnivariatePolynomialZ64 createFromArray(long[] data) {
        return new UnivariatePolynomialZ64(data);
    }

    @Override
    public UnivariatePolynomialZ64 createMonomial(long coefficient, int degree) {
        return monomial(coefficient, degree);
    }

    @Override
    public boolean isOverField() {return false;}

    @Override
    public boolean isOverFiniteField() {return false;}

    @Override
    public boolean isOverZ() {return true;}

    @Override
    public BigInteger coefficientRingCardinality() {return null;}

    @Override
    public BigInteger coefficientRingCharacteristic() {
        return BigInteger.ZERO;
    }

    @Override
    public boolean isOverPerfectPower() {
        return false;
    }

    @Override
    public BigInteger coefficientRingPerfectPowerBase() {
        return null;
    }

    @Override
    public BigInteger coefficientRingPerfectPowerExponent() {
        return null;
    }

    @Override
    long add(long a, long b) {return MachineArithmetic.safeAdd(a, b);}

    @Override
    long subtract(long a, long b) {return MachineArithmetic.safeSubtract(a, b);}

    @Override
    long multiply(long a, long b) {return MachineArithmetic.safeMultiply(a, b);}

    @Override
    long negate(long a) {return -a;}

    @Override
    long valueOf(long a) {return a;}

    @Override
    public UnivariatePolynomialZ64 monic() {
        if (isZero())
            return this;
        return divideOrNull(lc());
    }

    @Override
    public UnivariatePolynomialZ64 monic(long factor) {
        long lc = lc();
        long gcd = MachineArithmetic.gcd(lc, factor);
        factor = factor / gcd;
        lc = lc / gcd;
        UnivariatePolynomialZ64 r = divideOrNull(lc);
        if (r == null)
            return null;
        return r.multiply(factor);
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of
     * the elements can't be exactly divided by the {@code factor}. NOTE: is {@code null} is returned, the content of
     * {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public UnivariatePolynomialZ64 divideOrNull(long factor) {
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

    @Override
    public UnivariatePolynomialZ64 divideByLC(UnivariatePolynomialZ64 other) {
        return divideOrNull(other.lc());
    }

    @Override
    public UnivariatePolynomialZ64 multiplyByBigInteger(BigInteger factor) {
        return multiply(factor.longValueExact());
    }

    /** internal API */
    UnivariatePolynomialZ64 multiplyUnsafe(long factor) {
        for (int i = degree; i >= 0; --i)
            data[i] *= factor;
        return this;
    }

    @Override
    public UnivariatePolynomialZ64 multiply(UnivariatePolynomialZ64 oth) {
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

    /** internal API */
    UnivariatePolynomialZ64 multiplyUnsafe(UnivariatePolynomialZ64 oth) {
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

    @Override
    public UnivariatePolynomialZ64 square() {
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

    @Override
    public UnivariatePolynomialZ64 derivative() {
        if (isConstant())
            return createZero();
        long[] newData = new long[degree];
        for (int i = degree; i > 0; --i)
            newData[i - 1] = multiply(data[i], i);
        return createFromArray(newData);
    }

    @Override
    public UnivariatePolynomialZ64 clone() {
        return new UnivariatePolynomialZ64(data.clone(), degree);
    }

    @Override
    public UnivariatePolynomialZ64 parsePoly(String string) {
        return parse(string);
    }

    @Override
    public String coefficientRingToString() {
        return "Z";
    }
}
