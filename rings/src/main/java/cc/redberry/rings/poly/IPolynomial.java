package cc.redberry.rings.poly;

import cc.redberry.rings.WithVariables;
import cc.redberry.rings.bigint.BigInteger;

/**
 * Parent interface for all polynomials. All polynomial instances are mutable, so all structural operations except those
 * where it is stated explicitly will in general modify the instance. All arithmetic operations ({@code add(oth),
 * multiply(oth), monic()} etc.) applies to {@code this} inplace and return {@code this} reference ( so e.g. {@code
 * (poly == poly.add(other))}).
 *
 * <p><b>Note:</b> modifier operations are not synchronized.
 *
 * @param <Poly> the type of polynomial (self type)
 * @since 1.0
 */
public interface IPolynomial<Poly extends IPolynomial<Poly>>
        extends Comparable<Poly>, WithVariables, java.io.Serializable {
    /**
     * Returns whether {@code oth} and {@code this} have the same coefficient ring
     *
     * @param oth other polynomial
     * @return whether this and oth are over the same coefficient ring
     */
    boolean sameCoefficientRingWith(Poly oth);

    /**
     * Checks whether {@code oth} and {@code this} have the same coefficient ring, if not exception will be thrown
     *
     * @param oth other polynomial
     * @throws IllegalArgumentException if this and oth have different coefficient ring
     */
    default void assertSameCoefficientRingWith(Poly oth) {
        if (!sameCoefficientRingWith(oth))
            throw new IllegalArgumentException("Mixing polynomials over different coefficient rings: " + this.coefficientRingToString() + " and " + oth.coefficientRingToString());
    }

    /**
     * Set the coefficient ring from specified poly
     *
     * @param poly the polynomial
     * @return a copy of this with the coefficient ring taken from {@code poly}
     */
    Poly setCoefficientRingFrom(Poly poly);

    /**
     * Returns the degree of this polynomial
     *
     * @return the degree
     */
    int degree();

    /**
     * Returns the size of this polynomial
     *
     * @return the size
     */
    int size();

    /**
     * Returns {@code true} if this is zero
     *
     * @return whether {@code this} is zero
     */
    boolean isZero();

    /**
     * Returns {@code true} if this is one
     *
     * @return whether {@code this} is one
     */
    boolean isOne();

    /**
     * Returns {@code true} if this polynomial is monic
     *
     * @return whether {@code this} is monic
     */
    boolean isMonic();

    /**
     * Returns true if constant term is equal to one
     *
     * @return whether constant term is 1
     */
    boolean isUnitCC();

    /**
     * Returns {@code true} if this polynomial has only constant term
     *
     * @return whether {@code this} is constant
     */
    boolean isConstant();

    /**
     * Returns {@code true} if this polynomial has only one monomial term
     *
     * @return whether {@code this} has only one monomial term
     */
    boolean isMonomial();

    /**
     * Returns whether the coefficient ring of this polynomial is a field
     *
     * @return whether the coefficient ring of this polynomial is a field
     */
    boolean isOverField();

    /**
     * Returns whether the coefficient ring of this polynomial is Z
     *
     * @return whether the coefficient ring of this polynomial is Z
     */
    boolean isOverZ();

    /**
     * Returns whether the coefficient ring of this polynomial is a finite field
     *
     * @return whether the coefficient ring of this polynomial is a finite field
     */
    boolean isOverFiniteField();

    /**
     * Returns cardinality of the coefficient ring of this poly
     *
     * @return cardinality of the coefficient ring
     */
    BigInteger coefficientRingCardinality();

    /**
     * Returns characteristic of the coefficient ring of this poly
     *
     * @return characteristic of the coefficient ring
     */
    BigInteger coefficientRingCharacteristic();

    /**
     * Returns whether the {@code coefficientRingCardinality()} is a perfect power
     *
     * @return whether the {@code coefficientRingCardinality()} is a perfect power
     */
    boolean isOverPerfectPower();

    /**
     * Returns {@code base} so that {@code coefficientRingCardinality() == base^exponent} or null if cardinality is not
     * finite
     *
     * @return {@code base} so that {@code coefficientRingCardinality() == base^exponent} or null if cardinality is not
     * finite
     */
    BigInteger coefficientRingPerfectPowerBase();

    /**
     * Returns {@code exponent} so that {@code coefficientRingCardinality() == base^exponent} or null if cardinality is
     * not finite
     *
     * @return {@code exponent} so that {@code coefficientRingCardinality() == base^exponent} or null if cardinality is
     * not finite
     */
    BigInteger coefficientRingPerfectPowerExponent();

    /**
     * Sets {@code this} to its monic part (that is {@code this} divided by its leading coefficient), or returns {@code
     * null} (causing loss of internal data) if some of the elements can't be exactly divided by the {@code lc()}. NOTE:
     * if {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @return monic {@code this} or {@code null}
     */
    Poly monic();

    /**
     * Sets {@code this} to its monic part (that is {@code this} divided by its leading coefficient), or throws {@code
     * ArithmeticException} if some of the elements can't be exactly divided by the l.c.
     *
     * @return monic {@code this} or {@code null}
     * @throws ArithmeticException if some of the elements can't be exactly divided by the l.c.
     */
    default Poly monicExact() {
        Poly self = monic();
        if (self == null)
            throw new ArithmeticException("Not divisible by lc.");
        return self;
    }

    /**
     * Gives signum of the leading coefficient
     *
     * @return signum of the leading coefficient
     */
    int signumOfLC();

    /**
     * Sets this to zero
     *
     * @return this := zero
     */
    Poly toZero();

    /**
     * Sets the content of this to {@code oth}
     *
     * @param oth the polynomial
     * @return this := oth
     */
    Poly set(Poly oth);

    /**
     * Reduces poly to its primitive part (primitive part will always have positive l.c.)
     *
     * @return primitive part (poly will be modified)
     */
    Poly primitivePart();

    /**
     * Reduces poly to its primitive part, so that primitive part will have the same signum as the initial poly
     *
     * @return primitive part (poly will be modified)
     */
    Poly primitivePartSameSign();

    /**
     * Adds 1 to this
     *
     * @return {@code this + 1}
     */
    Poly increment();

    /**
     * Subtracts 1 from this
     *
     * @return {@code this - 1}
     */
    Poly decrement();

    /**
     * Returns the new instance of zero polynomial (with the same coefficient ring)
     *
     * @return new instance of 0
     */
    Poly createZero();

    /**
     * Returns the new instance of unit polynomial (with the same coefficient ring)
     *
     * @return new instance of 1
     */
    Poly createOne();

    /**
     * Creates constant polynomial with specified value
     *
     * @param value the value
     * @return constant polynomial
     */
    default Poly createConstant(long value) {
        return createOne().multiply(value);
    }

    /**
     * Adds {@code oth} to {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this + oth}
     */
    Poly add(Poly oth);


    /**
     * Adds {@code oth} to {@code this}.
     *
     * @param oth the polynomials
     * @return {@code this + oth}
     */
    @SuppressWarnings("unchecked")
    default Poly add(Poly... oth) {
        for (Poly t : oth)
            add(t);
        return (Poly) this;
    }

    /**
     * Subtracts {@code oth} from {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this - oth}
     */
    Poly subtract(Poly oth);

    /**
     * Subtracts {@code oth} from {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this - oth}
     */
    @SuppressWarnings("unchecked")
    default Poly subtract(Poly... oth) {
        for (Poly t : oth)
            subtract(t);
        return (Poly) this;
    }

    /**
     * Negates this and returns
     *
     * @return this negated
     */
    Poly negate();

    /**
     * Multiplies this by {@code oth }
     *
     * @param oth the polynomial
     * @return {@code this * oth }
     */
    Poly multiply(Poly oth);

    /**
     * Multiplies this by {@code oth }
     *
     * @param oth the polynomials
     * @return {@code this * oth }
     */
    @SuppressWarnings("unchecked")
    default Poly multiply(Poly... oth) {
        for (Poly t : oth)
            multiply(t);
        return (Poly) this;
    }

    /**
     * Multiplies this by {@code oth }
     *
     * @param oth the polynomials
     * @return {@code this * oth }
     */
    @SuppressWarnings("unchecked")
    default Poly multiply(Iterable<Poly> oth) {
        for (Poly t : oth)
            multiply(t);
        return (Poly) this;
    }

    /**
     * Multiplies this by {@code factor}
     *
     * @param factor the factor
     * @return {@code this * factor}
     */
    Poly multiply(long factor);

    /**
     * Multiplies this by {@code factor}
     *
     * @param factor the factor
     * @return {@code this * factor}
     */
    Poly multiplyByBigInteger(BigInteger factor);

    /**
     * Squares {@code this}
     *
     * @return {@code this * this}
     */
    Poly square();

    /**
     * Returns the content of this (gcd of coefficients) as a constant poly
     */
    Poly contentAsPoly();

    /**
     * Returns the leading coefficient as a constant poly
     */
    Poly lcAsPoly();

    /**
     * Returns the constant coefficient as a constant poly
     */
    Poly ccAsPoly();

    /**
     * Divides this polynomial by the leading coefficient of {@code other} or returns {@code null} (causing loss of
     * internal data) if some of the elements can't be exactly divided by the {@code other.lc()}. NOTE: if {@code null}
     * is returned, the content of {@code this} is destroyed.
     *
     * @param other the polynomial
     * @return {@code this} divided by the {@code other.lc()} or {@code null} if exact division is not possible
     */
    Poly divideByLC(Poly other);

    /**
     * Sets {@code this} to its monic part multiplied by the leading coefficient of {@code other};
     *
     * @param other other polynomial
     * @return monic part multiplied by the leading coefficient of {@code other} or null if exact division by the
     * reduced leading coefficient is not possible
     */
    Poly monicWithLC(Poly other);

    /**
     * Multiply this by the leading coefficient of {@code other}
     *
     * @param other polynomial
     * @return this * lc(other)
     */
    Poly multiplyByLC(Poly other);

    /**
     * Deep copy of this
     *
     * @return deep copy of this
     */
    Poly clone();

    /**
     * Deep copy of this (alias for {@link #clone()}, required for scala)
     *
     * @return deep copy of this
     */
    default Poly copy() { return clone(); }

    /** overcome Java generics... */
    Poly[] createArray(int length);

    /** overcome Java generics... */
    Poly[][] createArray2d(int length);

    /** overcome Java generics... */
    Poly[][] createArray2d(int length1, int length2);

    /** overcome Java generics... */
    default Poly[] createArray(Poly a, Poly b) {
        Poly[] r = createArray(2);
        r[0] = a;
        r[1] = b;
        return r;
    }

    /**
     * Parse string representation of polynomial
     *
     * @param string string
     * @return the polynomial corresponding to specified string
     */
    Poly parsePoly(String string);

    /**
     * Parse string representation of polynomial
     *
     * @param string    string
     * @param variables names of variables
     * @return the polynomial corresponding to specified string
     */
    Poly parsePoly(String string, String[] variables);

    /**
     * String representation of the coefficient ring of this
     */
    String coefficientRingToString();

    @Override
    String toString(String[] variables);
}
