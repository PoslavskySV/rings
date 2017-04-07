package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface IGeneralPolynomial<T extends IGeneralPolynomial> {
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
}
