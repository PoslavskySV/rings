package cc.r2.core.poly.univar;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
interface IMutablePolynomial<T extends IMutablePolynomial> extends Comparable<T> {
    /**
     * Return the degree of this
     *
     * @return the degree
     */
    int degree();

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
     * Returns true if constant term is a unit
     *
     * @return whether constant term is unit
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
     * @return whether {@code this} has the form {@code c*x^i} (one term)
     */
    boolean isMonomial();

    /**
     * Gives signum of leading coefficient
     *
     * @return signum of leading coefficient
     */
    int signum();

    /**
     * Returns position of the first non-zero coefficient, that is common monomial exponent (e.g. 2 for x^2 + x^3 + ...)
     *
     * @return position of the first non-zero coefficient
     */
    int firstNonZeroCoefficientPosition();

    /**
     * Sets this to zero
     *
     * @return this := zero
     */
    T toZero();

    /**
     * Sets the content of this to {@code oth}
     *
     * @param oth the polynomial
     * @return this := oth
     */
    T set(T oth);

    /**
     * Returns the quotient {@code this / x^offset}, it is polynomial with coefficient list formed by shifting coefficients
     * of {@code this} to the left by {@code offset}.
     *
     * @param offset shift amount
     * @return the quotient {@code this / x^offset}
     */
    T shiftLeft(int offset);

    /**
     * Multiplies {@code this} by the {@code x^offset}.
     *
     * @param offset monomial exponent
     * @return {@code this * x^offset}
     */
    T shiftRight(int offset);

    /**
     * Returns the remainder {@code this rem x^(newDegree + 1)}, it is polynomial with the coefficient list truncated
     * to the new degree {@code newDegree}.
     *
     * @param newDegree new degree
     * @return remainder {@code this rem x^(newDegree + 1)}
     */
    T truncate(int newDegree);

    /**
     * Creates polynomial formed from the coefficients of this starting from {@code from} (inclusive) to {@code to} (exclusive)
     *
     * @param from the initial index of the range to be copied, inclusive
     * @param to   the final index of the range to be copied, exclusive.
     * @return polynomial formed from the range of coefficients
     */
    T getRange(int from, int to);

    /**
     * Reverses the coefficients of this
     *
     * @return reversed polynomial
     */
    T reverse();

    /**
     * Reduces poly to its primitive part (primitive part will always have positive l.c.)
     *
     * @return primitive part (poly will be modified)
     */
    T primitivePart();

    /**
     * Reduces poly to its primitive part preserving signum of l.c.
     *
     * @return primitive part (poly will be modified)
     */
    T primitivePartSameSign();

    /**
     * Adds 1 to this
     *
     * @return {@code this + 1}
     */
    T increment();

    /**
     * Subtracts 1 from this
     *
     * @return {@code this - 1}
     */
    T decrement();

    /**
     * Returns 0 (new instance)
     *
     * @return new instance of 0
     */
    T createZero();

    /**
     * Returns 1 (new instance)
     *
     * @return new instance of 1
     */
    T createOne();

    /**
     * Creates monomial {@code x^degree}
     *
     * @param degree monomial degree
     * @return {@code coefficient * x^degree}
     */
    T createMonomial(int degree);

    /** overcome Java generics... */
    T[] arrayNewInstance(int length);

    /** overcome Java generics... */
    default T[] arrayNewInstance(T a, T b) {
        T[] r = arrayNewInstance(2);
        r[0] = a; r[1] = b;
        return r;
    }

    /**
     * Adds {@code oth} to {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this + oth}
     */
    T add(T oth);

    /**
     * Subtracts {@code oth} from {@code this}.
     *
     * @param oth the polynomial
     * @return {@code this - oth}
     */
    T subtract(T oth);


    /**
     * Negates this and returns
     *
     * @return this negated
     */
    T negate();

    /**
     * Sets this to {@code this * oth }
     *
     * @param oth the polynomial
     * @return {@code this * oth }
     */
    T multiply(T oth);

    @SuppressWarnings("unchecked")
    default T multiply(T... oth) {
        for (T t : oth)
            multiply(t);
        return (T) this;
    }

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    T multiply(long factor);

    /**
     * Square of {@code this}
     *
     * @return {@code this * this}
     */
    T square();

    /**
     * Returns the formal derivative of this poly
     *
     * @return the formal derivative
     */
    T derivative();

    /**
     * Deep copy of this
     *
     * @return deep copy of this
     */
    T clone();
}
