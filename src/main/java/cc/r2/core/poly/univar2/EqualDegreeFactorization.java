package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;

import static cc.r2.core.poly.univar2.Factorization.ensureFiniteFieldDomain;


/**
 * Equal-degree factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    @SuppressWarnings("unchecked")
    private static <T extends IMutablePolynomial<T>> T randomMonicPoly(T factory) {
        if (factory instanceof lMutablePolynomialZp) {
            lMutablePolynomialZp fm = (lMutablePolynomialZp) factory;
            return (T) RandomPolynomials.randomMonicPoly(fm.degree - 1, fm.domain.modulus, PrivateRandom.getRandom());
        } else {
            gMutablePolynomial fm = (gMutablePolynomial) factory;
            return (T) RandomPolynomials.randomMonicPoly(fm.degree - 1, fm.domain, PrivateRandom.getRandom());
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    static <Poly extends IMutablePolynomial<Poly>> FactorDecomposition<Poly> CantorZassenhaus(Poly input, int d) {
        ensureFiniteFieldDomain(input);
        FactorDecomposition<Poly> result = FactorDecomposition.constantFactor(input.lcAsPoly());
        CantorZassenhaus(input, d, result);
        return result;
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <T extends IMutablePolynomial<T>> void CantorZassenhaus(T input, int d, FactorDecomposition<T> result) {
        assert input.degree() % d == 0;
        int nFactors = input.degree() / d;
        if (input.degree() == 1 || nFactors == 1) {
            result.addFactor(input, 1);
            return;
        }

        assert input.degree() != d;

        T poly = input.clone();
        while (true) {

            poly = poly.monic();
            if (poly.degree() == d) {
                result.addFactor(poly, 1);
                break;
            }

            T factor;
            do {
                factor = CantorZassenhaus0(poly, d);
            } while (factor == null);

            if (factor.degree() != d)
                CantorZassenhaus(factor, d, result);
            else
                result.addFactor(factor.monic(), 1);
            poly = DivisionWithRemainder.quotient(poly, factor, false);
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <Poly extends IMutablePolynomial<Poly>> Poly CantorZassenhaus0(Poly poly, int d) {
        assert poly.isMonic();

        Poly a = randomMonicPoly(poly);
        if (a.isConstant())
            return null;

        Poly gcd1 = PolynomialGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant())
            return gcd1;

        // (modulus^d - 1) / 2
        BigInteger exponent = poly.domainCardinality().pow(d).decrement().shiftRight(1);
        Poly b = PolynomialArithmetics.polyPowMod(a, exponent, poly, DivisionWithRemainder.fastDivisionPreConditioning(poly), true);

        Poly gcd2 = PolynomialGCD.PolynomialGCD(b.decrement(), poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }
}
