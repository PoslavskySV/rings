package cc.redberry.rings.poly;

/**
 * @since 1.0
 */
public final class Util {
    private Util() {}

    public static void ensureOverFiniteField(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverFiniteField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureOverField(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureOverZ(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverZ())
                throw new IllegalArgumentException("Polynomial over Z is expected, but got " + poly.getClass());
    }
}
