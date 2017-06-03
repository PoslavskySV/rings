package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class CommonUtils {
    private CommonUtils() {}

    public static void ensureFiniteFieldDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverFiniteField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureFieldDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureZDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverZ())
                throw new IllegalArgumentException("Polynomial over Z is expected, but got " + poly.getClass());
    }
}
