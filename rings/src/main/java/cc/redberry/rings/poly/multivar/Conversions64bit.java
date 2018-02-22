package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.Util;

import java.util.List;
import java.util.stream.Collectors;

/**
 * @since 1.0
 */
class Conversions64bit {
    private Conversions64bit() {}

    /**
     * whether to switch to 64 bit integer arithmetic when possible (false in tests)
     */
    static boolean SWITCH_TO_64bit = true;

    static boolean canConvertToZp64(AMultivariatePolynomial poly) {
        return SWITCH_TO_64bit && Util.canConvertToZp64(poly);
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    MultivariatePolynomialZp64 asOverZp64(Poly poly) {
        return MultivariatePolynomial.asOverZp64((MultivariatePolynomial<BigInteger>) poly);
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<MultivariatePolynomialZp64> asOverZp64(List<Poly> list) {
        return list.stream()
                .map(Conversions64bit::asOverZp64)
                .collect(Collectors.toList());
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly convertFromZp64(MultivariatePolynomialZp64 p) {
        return (Poly) p.toBigPoly();
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> convertFromZp64(List<MultivariatePolynomialZp64> list) {
        return (List<Poly>) list.stream()
                .map(MultivariatePolynomialZp64::toBigPoly)
                .collect(Collectors.toList());
    }

    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] convertFromZp64(Poly factory, MultivariatePolynomialZp64[] p) {
        Poly[] r = factory.createArray(p.length);
        for (int i = 0; i < p.length; i++)
            r[i] = convertFromZp64(p[i]);
        return r;
    }
}
