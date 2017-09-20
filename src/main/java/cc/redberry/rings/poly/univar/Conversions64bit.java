package cc.redberry.rings.poly.univar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.Util;

/**
 * @since 1.0
 */
final class Conversions64bit {
    private Conversions64bit() {}

    /**
     * whether to switch to 64 bit integer arithmetic when possible (false in tests)
     */
    static boolean SWITCH_TO_64bit = true;

    static boolean canConvertToZp64(IUnivariatePolynomial poly) {
        return SWITCH_TO_64bit && Util.canConvertToZp64(poly);
    }

    @SuppressWarnings("unchecked")
    static <T extends IUnivariatePolynomial<T>> UnivariatePolynomialZp64 asOverZp64(T poly) {
        return UnivariatePolynomial.asOverZp64((UnivariatePolynomial<BigInteger>) poly);
    }

    @SuppressWarnings("unchecked")
    static <T extends IUnivariatePolynomial<T>> T convert(UnivariatePolynomialZp64 p) {
        return (T) p.toBigPoly();
    }

    @SuppressWarnings("unchecked")
    static <T extends IUnivariatePolynomial<T>> T[] convert(T factory, UnivariatePolynomialZp64[] p) {
        T[] r = factory.createArray(p.length);
        for (int i = 0; i < p.length; i++)
            r[i] = convert(p[i]);
        return r;
    }
}
