package cc.redberry.rings.scaladsl;

import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateDivision;

/** Ugly type erasures for interoperability with Java */
final class JavaConversions {
    private JavaConversions() {}

    @SuppressWarnings("unchecked")
    static <Poly> Poly[] divideAndRemainder(Poly dividend, Poly... dividers) {
        return (Poly[]) MultivariateDivision.divideAndRemainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) dividers);
    }

    @SuppressWarnings("unchecked")
    static <Poly> Poly remainder(Poly dividend, Poly... dividers) {
        return (Poly) MultivariateDivision.remainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) dividers);
    }

    @SuppressWarnings("unchecked")
    static <Poly> Poly swapVariables(Poly poly, int i, int j) {
        return (Poly) AMultivariatePolynomial.swapVariables((AMultivariatePolynomial) poly, i, j);
    }
}
