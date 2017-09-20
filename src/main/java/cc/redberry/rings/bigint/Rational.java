package cc.redberry.rings.bigint;

import org.apache.commons.math3.util.FastMath;

import java.io.Serializable;

/**
 * @since 1.0
 */
public final class Rational implements Comparable<Rational>, Serializable {

    /** A fraction representing "2 / 1". */
    public static final Rational TWO = new Rational(2);

    /** A fraction representing "1". */
    public static final Rational ONE = new Rational(1);

    /** A fraction representing "0". */
    public static final Rational ZERO = new Rational(0);

    /** A fraction representing "-1 / 1". */
    public static final Rational MINUS_ONE = new Rational(-1);

    /** A fraction representing "4/5". */
    public static final Rational FOUR_FIFTHS = new Rational(4, 5);

    /** A fraction representing "1/5". */
    public static final Rational ONE_FIFTH = new Rational(1, 5);

    /** A fraction representing "1/2". */
    public static final Rational ONE_HALF = new Rational(1, 2);

    /** A fraction representing "1/4". */
    public static final Rational ONE_QUARTER = new Rational(1, 4);

    /** A fraction representing "1/3". */
    public static final Rational ONE_THIRD = new Rational(1, 3);

    /** A fraction representing "3/5". */
    public static final Rational THREE_FIFTHS = new Rational(3, 5);

    /** A fraction representing "3/4". */
    public static final Rational THREE_QUARTERS = new Rational(3, 4);

    /** A fraction representing "2/5". */
    public static final Rational TWO_FIFTHS = new Rational(2, 5);

    /** A fraction representing "2/4". */
    public static final Rational TWO_QUARTERS = new Rational(2, 4);

    /** A fraction representing "2/3". */
    public static final Rational TWO_THIRDS = new Rational(2, 3);

    /** Serializable version identifier. */
    private static final long serialVersionUID = -5630213147331578515L;

    /** <code>BigInteger</code> representation of 100. */
    private static final BigInteger ONE_HUNDRED = BigInteger.valueOf(100);

    /** The numerator. */
    private final BigInteger numerator;

    /** The denominator. */
    private final BigInteger denominator;

    /**
     * <p> Create a {@link Rational} equivalent to the passed <tt>BigInteger</tt>, ie "num / 1". </p>
     *
     * @param num the numerator.
     */
    public Rational(final BigInteger num) {
        this(num, BigInteger.ONE);
    }

    /**
     * Create a {@link Rational} given the numerator and denominator as {@code BigInteger}. The {@link Rational} is
     * reduced to lowest terms.
     *
     * @param num the numerator, must not be {@code null}.
     * @param den the denominator, must not be {@code null}.
     */
    public Rational(BigInteger num, BigInteger den) {
        if (den.isZero())
            throw new ArithmeticException("divide by zero");
        if (num.isZero()) {
            numerator = BigInteger.ZERO;
            denominator = BigInteger.ONE;
        } else {

            // reduce numerator and denominator by greatest common denominator
            final BigInteger gcd = num.gcd(den);
            if (BigInteger.ONE.compareTo(gcd) < 0) {
                num = num.divide(gcd);
                den = den.divide(gcd);
            }

            // move sign to numerator
            if (den.signum() < 0) {
                num = num.negate();
                den = den.negate();
            }

            // store the values in the final fields
            numerator = num;
            denominator = den;
        }
    }

    /**
     * <p> Create a {@link Rational} equivalent to the passed <tt>int</tt>, ie "num / 1". </p>
     *
     * @param num the numerator.
     */
    public Rational(final int num) {
        this(BigInteger.valueOf(num), BigInteger.ONE);
    }

    /**
     * <p> Create a {@link Rational} given the numerator and denominator as simple <tt>int</tt>. The {@link Rational} is
     * reduced to lowest terms. </p>
     *
     * @param num the numerator.
     * @param den the denominator.
     */
    public Rational(final int num, final int den) {
        this(BigInteger.valueOf(num), BigInteger.valueOf(den));
    }

    /**
     * <p> Create a {@link Rational} equivalent to the passed long, ie "num / 1". </p>
     *
     * @param num the numerator.
     */
    public Rational(final long num) {
        this(BigInteger.valueOf(num), BigInteger.ONE);
    }

    /**
     * <p> Create a {@link Rational} given the numerator and denominator as simple <tt>long</tt>. The {@link Rational}
     * is reduced to lowest terms. </p>
     *
     * @param num the numerator.
     * @param den the denominator.
     */
    public Rational(final long num, final long den) {
        this(BigInteger.valueOf(num), BigInteger.valueOf(den));
    }

    /**
     * <p> Creates a <code>Rational</code> instance with the 2 parts of a fraction Y/Z. </p> <p> <p> Any negative signs
     * are resolved to be on the numerator. </p>
     *
     * @param numerator   the numerator, for example the three in 'three sevenths'.
     * @param denominator the denominator, for example the seven in 'three sevenths'.
     * @return a new fraction instance, with the numerator and denominator reduced.
     * @throws ArithmeticException if the denominator is <code>zero</code>.
     */
    public static Rational getReducedFraction(final int numerator,
                                              final int denominator) {
        if (numerator == 0)
            return ZERO; // normalize zero.

        return new Rational(numerator, denominator);
    }

    public boolean isZero() {
        return numerator.isZero();
    }

    public boolean isOne() {
        return numerator.isOne() && denominator.isOne();
    }

    public int signum() {
        return numerator.signum();
    }

    /**
     * <p> Returns the absolute value of this {@link Rational}. </p>
     *
     * @return the absolute value as a {@link Rational}.
     */
    public Rational abs() {
        return (BigInteger.ZERO.compareTo(numerator) <= 0) ? this : negate();
    }

    /**
     * <p> Adds the value of this fraction to the passed {@link BigInteger}, returning the result in reduced form. </p>
     *
     * @param bg the {@link BigInteger} to add, must'nt be <code>null</code>.
     * @return a <code>Rational</code> instance with the resulting values.
     */
    public Rational add(final BigInteger bg) {
        return new Rational(numerator.add(denominator.multiply(bg)), denominator);
    }

    /**
     * <p> Adds the value of this fraction to the passed <tt>integer</tt>, returning the result in reduced form. </p>
     *
     * @param i the <tt>integer</tt> to add.
     * @return a <code>Rational</code> instance with the resulting values.
     */
    public Rational add(final int i) {
        return add(BigInteger.valueOf(i));
    }

    /**
     * <p> Adds the value of this fraction to the passed <tt>long</tt>, returning the result in reduced form. </p>
     *
     * @param l the <tt>long</tt> to add.
     * @return a <code>Rational</code> instance with the resulting values.
     */
    public Rational add(final long l) {
        return add(BigInteger.valueOf(l));
    }

    /**
     * <p> Adds the value of this fraction to another, returning the result in reduced form. </p>
     *
     * @param fraction the {@link Rational} to add, must not be <code>null</code>.
     * @return a {@link Rational} instance with the resulting values.
     */
    public Rational add(final Rational fraction) {
        if (fraction.isZero())
            return this;

        BigInteger num, den;

        if (denominator.equals(fraction.denominator)) {
            num = numerator.add(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).add((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return new Rational(num, den);

    }

    /**
     * <p> Gets the fraction as a <code>BigDecimal</code>. This calculates the fraction as the numerator divided by
     * denominator. </p>
     *
     * @return the fraction as a <code>BigDecimal</code>.
     * @throws ArithmeticException if the exact quotient does not have a terminating decimal expansion.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue() {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator));
    }

    /**
     * <p> Gets the fraction as a <code>BigDecimal</code> following the passed rounding mode. This calculates the
     * fraction as the numerator divided by denominator. </p>
     *
     * @param roundingMode rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @throws IllegalArgumentException if <tt>roundingMode</tt> does not represent a valid rounding mode.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue(final int roundingMode) {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator), roundingMode);
    }

    /**
     * <p> Gets the fraction as a <code>BigDecimal</code> following the passed scale and rounding mode. This calculates
     * the fraction as the numerator divided by denominator. </p>
     *
     * @param scale        scale of the <code>BigDecimal</code> quotient to be returned. see {@link BigDecimal} for more
     *                     information.
     * @param roundingMode rounding mode to apply. see {@link BigDecimal} constants.
     * @return the fraction as a <code>BigDecimal</code>.
     * @see BigDecimal
     */
    public BigDecimal bigDecimalValue(final int scale, final int roundingMode) {
        return new BigDecimal(numerator).divide(new BigDecimal(denominator), scale, roundingMode);
    }

    /**
     * <p> Compares this object to another based on size. </p>
     *
     * @param object the object to compare to, must not be <code>null</code>.
     * @return -1 if this is less than <tt>object</tt>, +1 if this is greater than <tt>object</tt>, 0 if they are equal.
     * @see java.lang.Comparable#compareTo(java.lang.Object)
     */
    public int compareTo(final Rational object) {
        BigInteger nOd = numerator.multiply(object.denominator);
        BigInteger dOn = denominator.multiply(object.numerator);
        return nOd.compareTo(dOn);
    }

    /**
     * <p> Divide the value of this fraction by the passed {@code BigInteger}, ie {@code this * 1 / bg}, returning the
     * result in reduced form. </p>
     *
     * @param bg the {@code BigInteger} to divide by, must not be {@code null}
     * @return a {@link Rational} instance with the resulting values
     */
    public Rational divide(final BigInteger bg) {
        if (bg.isZero())
            throw new ArithmeticException("divide by zero");
        return new Rational(numerator, denominator.multiply(bg));
    }

    /**
     * <p> Divide the value of this fraction by the passed {@code int}, ie {@code this * 1 / i}, returning the result in
     * reduced form. </p>
     *
     * @param i the {@code int} to divide by
     * @return a {@link Rational} instance with the resulting values
     */
    public Rational divide(final int i) {
        return divide(BigInteger.valueOf(i));
    }

    /**
     * <p> Divide the value of this fraction by the passed {@code long}, ie {@code this * 1 / l}, returning the result
     * in reduced form. </p>
     *
     * @param l the {@code long} to divide by
     * @return a {@link Rational} instance with the resulting values
     */
    public Rational divide(final long l) {
        return divide(BigInteger.valueOf(l));
    }

    /**
     * <p> Divide the value of this fraction by another, returning the result in reduced form. </p>
     *
     * @param fraction Fraction to divide by, must not be {@code null}.
     * @return a {@link Rational} instance with the resulting values.
     */
    public Rational divide(final Rational fraction) {
        if (fraction.numerator.isZero())
            throw new ArithmeticException("divide by zero");

        return multiply(fraction.reciprocal());
    }

    /**
     * <p> Gets the fraction as a <tt>double</tt>. This calculates the fraction as the numerator divided by denominator.
     * </p>
     *
     * @return the fraction as a <tt>double</tt>
     * @see java.lang.Number#doubleValue()
     */
    public double doubleValue() {
        double result = numerator.doubleValue() / denominator.doubleValue();
        if (Double.isNaN(result)) {
            // Numerator and/or denominator must be out of range:
            // Calculate how far to shift them to put them in range.
            int shift = Math.max(numerator.bitLength(),
                    denominator.bitLength()) - FastMath.getExponent(Double.MAX_VALUE);
            result = numerator.shiftRight(shift).doubleValue() /
                    denominator.shiftRight(shift).doubleValue();
        }
        return result;
    }

    /**
     * Access the denominator as a <code>BigInteger</code>.
     *
     * @return the denominator as a <code>BigInteger</code>.
     */
    public BigInteger getDenominator() {
        return denominator;
    }

    /**
     * Access the denominator as a <tt>int</tt>.
     *
     * @return the denominator as a <tt>int</tt>.
     */
    public int getDenominatorAsInt() {
        return denominator.intValue();
    }

    /**
     * <p> Access the denominator as a <tt>long</tt>. </p>
     *
     * @return the denominator as a <tt>long</tt>.
     */
    public long getDenominatorAsLong() {
        return denominator.longValue();
    }

    /**
     * <p> Access the numerator as a <code>BigInteger</code>. </p>
     *
     * @return the numerator as a <code>BigInteger</code>.
     */
    public BigInteger getNumerator() {
        return numerator;
    }

    /**
     * <p> Access the numerator as a <tt>int</tt>. </p>
     *
     * @return the numerator as a <tt>int</tt>.
     */
    public int getNumeratorAsInt() {
        return numerator.intValue();
    }

    /**
     * <p> Access the numerator as a <tt>long</tt>. </p>
     *
     * @return the numerator as a <tt>long</tt>.
     */
    public long getNumeratorAsLong() {
        return numerator.longValue();
    }

    /**
     * <p> Multiplies the value of this fraction by the passed <code>BigInteger</code>, returning the result in reduced
     * form. </p>
     *
     * @param bg the {@code BigInteger} to multiply by.
     * @return a {@code Rational} instance with the resulting values.
     */
    public Rational multiply(final BigInteger bg) {
        return new Rational(bg.multiply(numerator), denominator);
    }

    /**
     * <p> Multiply the value of this fraction by the passed <tt>int</tt>, returning the result in reduced form. </p>
     *
     * @param i the <tt>int</tt> to multiply by.
     * @return a {@link Rational} instance with the resulting values.
     */
    public Rational multiply(final int i) {
        return multiply(BigInteger.valueOf(i));
    }

    /**
     * <p> Multiply the value of this fraction by the passed <tt>long</tt>, returning the result in reduced form. </p>
     *
     * @param l the <tt>long</tt> to multiply by.
     * @return a {@link Rational} instance with the resulting values.
     */
    public Rational multiply(final long l) {
        return multiply(BigInteger.valueOf(l));
    }

    /**
     * <p> Multiplies the value of this fraction by another, returning the result in reduced form. </p>
     *
     * @param fraction Fraction to multiply by, must not be {@code null}.
     * @return a {@link Rational} instance with the resulting values.
     */
    public Rational multiply(final Rational fraction) {
        if (numerator.isZero() || fraction.numerator.isZero())
            return ZERO;
        return new Rational(numerator.multiply(fraction.numerator),
                denominator.multiply(fraction.denominator));
    }

    /**
     * <p> Return the additive inverse of this fraction, returning the result in reduced form. </p>
     *
     * @return the negation of this fraction.
     */
    public Rational negate() {
        return new Rational(numerator.negate(), denominator);
    }

    /**
     * <p> Returns a {@code Rational} whose value is {@code (this<sup>exponent</sup>)}, returning the result in reduced
     * form. </p>
     *
     * @param exponent exponent to which this {@code Rational} is to be raised.
     * @return <tt>this<sup>exponent</sup></tt>.
     */
    public Rational pow(final int exponent) {
        if (exponent < 0)
            return new Rational(denominator.pow(-exponent), numerator.pow(-exponent));
        return new Rational(numerator.pow(exponent), denominator.pow(exponent));
    }

    /**
     * <p> Returns a <code>Rational</code> whose value is <tt>(this<sup>exponent</sup>)</tt>, returning the result in
     * reduced form. </p>
     *
     * @param exponent exponent to which this <code>Rational</code> is to be raised.
     * @return <tt>this<sup>exponent</sup></tt> as a <code>Rational</code>.
     */
    public Rational pow(final long exponent) {
        if (exponent < 0)
            return new Rational(BigIntegerUtil.pow(denominator, -exponent),
                    BigIntegerUtil.pow(numerator, -exponent));
        return new Rational(BigIntegerUtil.pow(numerator, exponent),
                BigIntegerUtil.pow(denominator, exponent));
    }

    /**
     * <p> Returns a <code>Rational</code> whose value is <tt>(this<sup>exponent</sup>)</tt>, returning the result in
     * reduced form. </p>
     *
     * @param exponent exponent to which this <code>Rational</code> is to be raised.
     * @return <tt>this<sup>exponent</sup></tt> as a <code>Rational</code>.
     */
    public Rational pow(final BigInteger exponent) {
        if (exponent.compareTo(BigInteger.ZERO) < 0) {
            final BigInteger eNeg = exponent.negate();
            return new Rational(BigIntegerUtil.pow(denominator, eNeg),
                    BigIntegerUtil.pow(numerator, eNeg));
        }
        return new Rational(BigIntegerUtil.pow(numerator, exponent),
                BigIntegerUtil.pow(denominator, exponent));
    }

    /**
     * <p> Return the multiplicative inverse of this fraction. </p>
     *
     * @return the reciprocal fraction.
     */
    public Rational reciprocal() {
        return new Rational(denominator, numerator);
    }

    /**
     * <p> Subtracts the value of an {@link BigInteger} from the value of this {@code Rational}, returning the result in
     * reduced form. </p>
     *
     * @param bg the {@link BigInteger} to subtract, cannot be {@code null}.
     * @return a {@code Rational} instance with the resulting values.
     */
    public Rational subtract(final BigInteger bg) {
        return new Rational(numerator.subtract(denominator.multiply(bg)), denominator);
    }

    /**
     * <p> Subtracts the value of an {@code integer} from the value of this {@code Rational}, returning the result in
     * reduced form. </p>
     *
     * @param i the {@code integer} to subtract.
     * @return a {@code Rational} instance with the resulting values.
     */
    public Rational subtract(final int i) {
        return subtract(BigInteger.valueOf(i));
    }

    /**
     * <p> Subtracts the value of a {@code long} from the value of this {@code Rational}, returning the result in
     * reduced form. </p>
     *
     * @param l the {@code long} to subtract.
     * @return a {@code Rational} instance with the resulting values.
     */
    public Rational subtract(final long l) {
        return subtract(BigInteger.valueOf(l));
    }

    /**
     * <p> Subtracts the value of another fraction from the value of this one, returning the result in reduced form.
     * </p>
     *
     * @param fraction {@link Rational} to subtract, must not be {@code null}.
     * @return a {@link Rational} instance with the resulting values
     */
    public Rational subtract(final Rational fraction) {
        if (fraction.isZero())
            return this;

        BigInteger num, den;
        if (denominator.equals(fraction.denominator)) {
            num = numerator.subtract(fraction.numerator);
            den = denominator;
        } else {
            num = (numerator.multiply(fraction.denominator)).subtract((fraction.numerator).multiply(denominator));
            den = denominator.multiply(fraction.denominator);
        }
        return new Rational(num, den);

    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Rational that = (Rational) o;

        if (!numerator.equals(that.numerator)) return false;
        return denominator.equals(that.denominator);
    }

    @Override
    public int hashCode() {
        int result = numerator.hashCode();
        result = 31 * result + denominator.hashCode();
        return result;
    }

    /**
     * <p> Returns the <code>String</code> representing this fraction, ie "num / dem" or just "num" if the denominator
     * is one. </p>
     *
     * @return a string representation of the fraction.
     * @see java.lang.Object#toString()
     */
    @Override
    public String toString() {
        String str = null;
        if (denominator.isOne()) {
            str = numerator.toString();
        } else if (numerator.isZero()) {
            str = "0";
        } else {
            str = numerator + "/" + denominator;
        }
        return str;
    }
}
