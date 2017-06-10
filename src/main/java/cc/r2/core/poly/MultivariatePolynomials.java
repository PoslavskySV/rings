package cc.r2.core.poly;

import cc.r2.core.poly.multivar.AMultivariatePolynomial;
import cc.r2.core.poly.multivar.MultivariateGCD;
import cc.r2.core.poly.multivar.MultivariateReduction;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomials<Poly extends AMultivariatePolynomial<?, Poly>> extends APolynomialsDomain<Poly> {
    public MultivariatePolynomials(Poly factory) {
        super(factory);
    }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        Poly[] arr = divider.arrayNewInstance(1);
        arr[0] = divider;
        return (Poly[]) MultivariateReduction.divideAndRemainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) arr);
    }

    @Override
    public Poly gcd(Poly a, Poly b) {
        return MultivariateGCD.PolynomialGCD(a, b);
    }
}
