package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface ToStringSupport<E> {
    String toString(E element);

    ToStringSupport plain = Object::toString;

    @SuppressWarnings("unchecked")
    static <E> ToStringSupport<E> plain() { return (ToStringSupport<E>) plain;}

    static <P extends IPolynomial<P>> ToStringSupport<P> poly(String[] vars) { return p -> p.toString(vars);}
}
