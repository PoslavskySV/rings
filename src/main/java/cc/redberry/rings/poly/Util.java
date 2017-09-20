package cc.redberry.rings.poly;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;

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

    /**
     * Test whether poly is over Zp with modulus less then 2^63
     */
    public static boolean canConvertToZp64(IPolynomial poly) {
        Ring ring;
        if (poly instanceof UnivariatePolynomial)
            ring = ((UnivariatePolynomial) poly).ring;
        else if (poly instanceof MultivariatePolynomial)
            ring = ((MultivariatePolynomial) poly).ring;
        else
            return false;

        return ring instanceof IntegersZp && ((IntegersZp) ring).modulus.bitLength() < MachineArithmetic.MAX_SUPPORTED_MODULUS_BITS;
    }
}
