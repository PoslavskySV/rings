package cc.redberry.rings.poly;

import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateDivision;
import cc.redberry.rings.poly.multivar.MultivariateGCD;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Ring of multivariate polynomials.
 *
 * @param <Poly> type of multivariate polynomials
 * @since 1.0
 */
public final class MultivariateRing<Poly extends AMultivariatePolynomial<?, Poly>> extends APolynomialRing<Poly> {
    private static final long serialVersionUID = 1L;

    /**
     * Creates ring of multivariate polynomials which support operations over multivariate polynomials of the type and
     * number of variables same as of provided {@code factory} polynomial
     *
     * @param factory the factory polynomial (the exact value of {@code factory} is irrelevant) which fixes the element
     *                type of this ring, coefficient ring and the number of variables
     */
    public MultivariateRing(Poly factory) { super(factory); }

    @Override
    public int nVariables() { return factory.nVariables; }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        Poly[] arr = divider.createArray(1);
        arr[0] = divider;
        return (Poly[]) MultivariateDivision.divideAndRemainder((AMultivariatePolynomial) dividend, (AMultivariatePolynomial[]) arr);
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

    @Override
    @SuppressWarnings("unchecked")
    public Poly variable(int variable) {
        return factory.createMonomial(variable, 1);
    }

    /**
     * Gives a random constant polynomial. For generating non-constant random polynomials see {@link
     * cc.redberry.rings.poly.multivar.RandomMultivariatePolynomials}
     *
     * @param rnd the source of randomness
     * @return random constant polynomial
     * @see cc.redberry.rings.poly.multivar.RandomMultivariatePolynomials
     */
    @Override
    public Poly randomElement(RandomGenerator rnd) {
        return super.randomElement(rnd);
    }
}
