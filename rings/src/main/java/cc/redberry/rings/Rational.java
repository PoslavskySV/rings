package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;

import java.util.function.Function;
import java.util.stream.Stream;

/**
 * Rational expression with numerator and denominator from the ring. Objects of this class are immutable.
 *
 * @param <E> type
 * @since 1.0
 */
public final class Rational<E>
        implements Comparable<Rational<E>>,
                   java.io.Serializable {
    private static final long serialVersionUID = 1L;
    /** The ring. */
    public final Ring<E> ring;
    /** The numerator. */
    public final E numerator;
    /** The denominator. */
    public final E denominator;

    /**
     * Constructs rational with the specified numerator and denominator
     */
    public Rational(Ring<E> ring, E numerator, E denominator) {
        if (ring.isZero(denominator))
            throw new ArithmeticException("division by zero");
        this.ring = ring;
        if (!ring.isZero(numerator)) {
            E gcd = ring.gcd(numerator, denominator);
            numerator = ring.divideExact(numerator, gcd);
            denominator = ring.divideExact(denominator, gcd);

            if (ring.isUnit(denominator)) {
                numerator = ring.divideExact(numerator, denominator);
                denominator = ring.getOne();
            }
        }

        if (ring.signum(denominator) < 0) {
            denominator = ring.negate(denominator);
            numerator = ring.negate(numerator);
        }
        this.numerator = numerator;
        this.denominator = denominator;
    }

    /**
     * Constructs rational with the specified numerator and unit denominator
     */
    public Rational(Ring<E> ring, E numerator) {
        this.ring = ring;
        this.numerator = numerator;
        this.denominator = ring.getOne();
    }

    /**
     * Constructs zero
     */
    public static <E> Rational<E> zero(Ring<E> ring) {
        return new Rational<>(ring, ring.getZero());
    }

    /**
     * Constructs one
     */
    public static <E> Rational<E> one(Ring<E> ring) {
        return new Rational<>(ring, ring.getOne());
    }

    /* private constructor without gcd reduction */
    private Rational(boolean b, Ring<E> ring, E numerator, E denominator) {
        if (ring.isUnit(denominator)) {
            numerator = ring.divideExact(numerator, denominator);
            denominator = ring.getOne();
        }
        this.ring = ring;
        this.numerator = numerator;
        this.denominator = denominator;
    }

    /**
     * Tests whether the denominator is one
     */
    public boolean isIntegral() { return ring.isOne(denominator);}

    /**
     * Whether this is zero
     */
    public boolean isZero() {
        return ring.isZero(numerator);
    }

    /**
     * Whether this is one
     */
    public boolean isOne() {
        return ring.isOne(numerator) && ring.isOne(denominator);
    }

    /**
     * Signum of this rational
     */
    public int signum() {
        return ring.signum(numerator);
    }

    /**
     * Returns the absolute value of this {@link Rational}.
     *
     * @return the absolute value as a {@link Rational}.
     */
    public Rational<E> abs() {
        return (ring.signum(numerator) <= 0) ? this : negate();
    }

    /**
     * Add {@code other} to this
     */
    public Rational<E> add(E val) {
        return new Rational<>(ring, ring.add(numerator, ring.multiply(denominator, val)), denominator);
    }

    /**
     * Add {@code other} to this
     */
    public Rational<E> add(long other) {
        return add(ring.valueOf(other));
    }

    /**
     * Add {@code other} to this
     */
    public Rational<E> add(Rational<E> other) {
        if (other.isZero())
            return this;

        E num, den;
        if (denominator.equals(other.denominator)) {
            num = ring.add(numerator, other.numerator);
            den = denominator;
        } else {
            num = ring.add(ring.multiply(numerator, other.denominator), ring.multiply(denominator, other.numerator));
            den = ring.multiply(denominator, other.denominator);
        }
        return new Rational<>(ring, num, den);
    }

    @Override
    public int compareTo(Rational<E> object) {
        E nOd = ring.multiply(numerator, object.denominator);
        E dOn = ring.multiply(denominator, object.numerator);
        return ring.compare(nOd, dOn);
    }

    /**
     * Divide this by {@code other}
     */
    public Rational<E> divide(E other) {
        if (ring.isZero(other))
            throw new ArithmeticException("divide by zero");
        return new Rational<>(ring, numerator, ring.multiply(denominator, other));
    }


    /**
     * Divide this by {@code other}
     */
    public Rational<E> divide(long l) {
        return divide(ring.valueOf(l));
    }

    /**
     * Divide this by {@code other}
     */

    public Rational<E> divide(Rational other) {
        if (other.isZero())
            throw new ArithmeticException("divide by zero");
        return multiply(other.reciprocal());
    }

    /**
     * Multiply this by {@code other}
     */
    public Rational<E> multiply(E other) {
        return new Rational<>(ring, ring.multiply(other, numerator), denominator);
    }

    /**
     * Multiply this by {@code other}
     */
    public Rational<E> multiply(long other) {
        return multiply(ring.valueOf(other));
    }

    /**
     * Multiply this by {@code other}
     */
    public Rational<E> multiply(Rational<E> other) {
        if (isZero())
            return this;
        if (other.isZero())
            return other;
        return new Rational<>(ring, ring.multiply(numerator, other.numerator), ring.multiply(denominator, other.denominator));
    }

    /**
     * Negates this
     */
    public Rational<E> negate() {
        return new Rational<>(ring, ring.negate(numerator), denominator);
    }

    /**
     * Return the multiplicative inverse of this rational.
     */
    public Rational<E> reciprocal() {
        return new Rational<>(true, ring, denominator, numerator);
    }

    /**
     * Subtracts {@code other} from this
     */
    public Rational<E> subtract(E other) {
        return new Rational<>(ring, ring.subtract(numerator, ring.multiply(denominator, other)), denominator);
    }

    /**
     * Subtracts {@code other} from this
     */
    public Rational<E> subtract(final long l) {
        return subtract(ring.valueOf(l));
    }

    /**
     * Subtracts {@code other} from this
     */
    public Rational<E> subtract(final Rational<E> other) {
        if (other.isZero())
            return this;

        E num, den;
        if (denominator.equals(other.denominator)) {
            num = ring.subtract(numerator, other.numerator);
            den = denominator;
        } else {
            num = ring.subtract(ring.multiply(numerator, other.denominator), ring.multiply(other.numerator, denominator));
            den = ring.multiply(denominator, other.denominator);
        }
        return new Rational<>(ring, num, den);
    }

    /**
     * Raise this in a power {@code exponent}
     *
     * @param exponent exponent
     */
    public Rational<E> pow(int exponent) {
        return new Rational<>(ring, ring.pow(numerator, exponent), ring.pow(denominator, exponent));
    }

    /**
     * Raise this in a power {@code exponent}
     *
     * @param exponent exponent
     */
    public Rational<E> pow(long exponent) {
        return new Rational<>(ring, ring.pow(numerator, exponent), ring.pow(denominator, exponent));
    }

    /**
     * Raise this in a power {@code exponent}
     *
     * @param exponent exponent
     */
    public Rational<E> pow(BigInteger exponent) {
        return new Rational<>(ring, ring.pow(numerator, exponent), ring.pow(denominator, exponent));
    }

    /**
     * Maps rational to a new ring
     */
    public <O> Rational<O> map(Ring<O> ring, Function<E, O> function) {
        return new Rational<>(ring, function.apply(numerator), function.apply(denominator));
    }

    /**
     * Stream of numerator and denominator
     */
    public Stream<E> stream() {
        return Stream.of(numerator, denominator);
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

    @Override
    public String toString() {
        return toString((ToStringSupport<E>) null);
    }

    public String toString(ToStringSupport<E> toString) {
        String str;
        if (ring.isOne(denominator))
            str = toString(numerator, toString, false);
        else if (ring.isZero(numerator))
            str = "0";
        else
            str = toString(numerator, toString, true) + "/" + toString(denominator, toString, true);
        return str;
    }
//
//    @Override
//    public String toString(String[] variables) {
//        return toString(numerator instanceof WithVariables
//                ? (ToStringSupport<E>) (ToStringSupport.withVariables(variables))
//                : null);
//    }

    private String toString(E e, ToStringSupport<E> toString, boolean needBrackets) {
        String str;
        if (toString != null)
            str = toString.toString(e);
        else
            str = ring.toString(e);
        if (needBrackets && !(str.matches("[0-9]+") || str.matches("[a-zA-Z]+")))
            str = "(" + str + ")";
        return str;
    }
}
