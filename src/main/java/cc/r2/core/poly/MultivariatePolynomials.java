package cc.r2.core.poly;

import cc.r2.core.poly.multivar.AMultivariatePolynomial;
import cc.r2.core.poly.multivar.MultivariateGCD;
import cc.r2.core.poly.multivar.MultivariateReduction;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Domain of multivariate polynomials.
 *
 * @param <Poly> type of univariate polynomials
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomials<Poly extends AMultivariatePolynomial<?, Poly>> extends APolynomialsDomain<Poly> {
    private static final long serialVersionUID = 1L;

    /**
     * Creates the domain of multivariate polynomials which support operations over multivariate polynomials of the
     * type and number of variables same as of provided {@code factory} polynomial
     *
     * @param factory a representative polynomial (the exact value of {@code factory} is irrelevant) which fixes
     *                the element type of this domain and number of variables
     */
    public MultivariatePolynomials(Poly factory) { super(factory); }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        Poly[] arr = divider.createArray(1);
        arr[0] = divider;
        return (Poly[]) MultivariateReduction.divideAndRemainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) arr);
    }

    @Override
    public Poly gcd(Poly a, Poly b) {
        return MultivariateGCD.PolynomialGCD(a, b);
    }

    @Override
    public Poly gcd(Poly[] elements) {
        return MultivariateGCD.PolynomialGCD(elements);
    }

    @Override
    public Poly gcd(Iterable<Poly> elements) {
        return MultivariateGCD.PolynomialGCD(elements);
    }

    /**
     * Gives a random constant polynomial
     *
     * @param rnd the source of randomness
     * @return random constant polynomial
     */
    @Override
    public Poly randomElement(RandomGenerator rnd) {
        return super.randomElement(rnd);
    }

}
