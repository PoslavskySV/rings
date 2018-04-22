package cc.redberry.rings;

/**
 * @since 1.0
 */
public interface ToStringSupport<E> {
    /**
     * Gives string representation of specified element
     *
     * @param element the element
     * @return string representation
     */
    String toString(E element);

    /**
     * {@code Object::toString}
     */
    ToStringSupport plain = Object::toString;

    /**
     * {@code Object::toString}
     */
    @SuppressWarnings("unchecked")
    static <E> ToStringSupport<E> plain() { return (ToStringSupport<E>) plain;}

    /**
     * {@code p -> p.toString(vars)}
     */
    static <P extends WithVariables> ToStringSupport<P> withVariables(String[] vars) { return p -> p.toString(vars);}
}
