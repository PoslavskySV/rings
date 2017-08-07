package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.univar.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariatePolynomials<Poly extends IUnivariatePolynomial<Poly>> extends APolynomialsDomain<Poly> {
    /**
     * Polynomials over Z (arbitrary precision integers)
     */
    public static final UnivariatePolynomials<UnivariatePolynomial<BigInteger>> POLYNOMIALS_OVER_Z = new UnivariatePolynomials<>(UnivariatePolynomial.zero(Integers.Integers));

    /**
     * Polynomials over Z (machine precision integers)
     */
    public static final UnivariatePolynomials<lUnivariatePolynomialZ> lPOLYNOMIALS_OVER_Z = new UnivariatePolynomials<>(lUnivariatePolynomialZ.zero());

    /**
     * Polynomials over specified domain
     *
     * @param domain the domain
     */
    public static <E> UnivariatePolynomials<UnivariatePolynomial<E>> overDomain(Domain<E> domain) {
        return new UnivariatePolynomials<>(UnivariatePolynomial.zero(domain));
    }

    /**
     * Polynomials over Z/p (arbitrary precision integers)
     *
     * @param modulus the modulus
     */
    public static UnivariatePolynomials<UnivariatePolynomial<BigInteger>> overZp(BigInteger modulus) {
        return new UnivariatePolynomials<>(UnivariatePolynomial.zero(new IntegersModulo(modulus)));
    }

    public UnivariatePolynomials(Poly factory) {
        super(factory);
    }

    @Override
    public Poly remainder(Poly a, Poly b) {return DivisionWithRemainder.remainder(a, b, true);}

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {return DivisionWithRemainder.divideAndRemainder(a, b, true);}

    @Override
    public Poly gcd(Poly a, Poly b) {return UnivariateGCD.PolynomialGCD(a, b);}

//    @Override
//    public Poly randomElement(RandomGenerator rnd) {
//        return RandomPolynomials.randomPoly(32, factory.);
//    }
}
