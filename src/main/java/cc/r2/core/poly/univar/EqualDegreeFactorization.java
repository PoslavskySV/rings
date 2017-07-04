package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.LongArithmetics;

import static cc.r2.core.poly.CommonUtils.ensureFiniteFieldDomain;


/**
 * Equal-degree factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    @SuppressWarnings("unchecked")
    private static <T extends IUnivariatePolynomial<T>> T randomMonicPoly(T factory) {
        if (factory instanceof lUnivariatePolynomialZp) {
            lUnivariatePolynomialZp fm = (lUnivariatePolynomialZp) factory;
            return (T) RandomPolynomials.randomMonicPoly(fm.degree - 1, fm.domain.modulus, PrivateRandom.getRandom());
        } else {
            UnivariatePolynomial fm = (UnivariatePolynomial) factory;
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
    static <Poly extends IUnivariatePolynomial<Poly>> FactorDecomposition<Poly> CantorZassenhaus(Poly input, int d) {
        ensureFiniteFieldDomain(input);
        FactorDecomposition<Poly> result = FactorDecomposition.constantFactor(input.lcAsPoly());
        if (!input.coefficientDomainCardinality().testBit(0))
            //even characteristics => GF2p
            CantorZassenhaus(input, d, result, input.coefficientDomainPerfectPowerExponent().intValueExact());
        else
            CantorZassenhaus(input, d, result, -1);
        return result;
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <T extends IUnivariatePolynomial<T>> void CantorZassenhaus(T input, int d, FactorDecomposition<T> result, int pPower) {
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
                if (pPower == -1)
                    factor = CantorZassenhaus0(poly, d);
                else
                    factor = CantorZassenhausGF2p(poly, d, pPower);
            } while (factor == null);

            if (factor.degree() != d)
                CantorZassenhaus(factor, d, result, pPower);
            else
                result.addFactor(factor.monic(), 1);
            poly = DivisionWithRemainder.quotient(poly, factor, false);
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm for odd characteristics
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <Poly extends IUnivariatePolynomial<Poly>> Poly CantorZassenhaus0(Poly poly, int d) {
        assert poly.isMonic();

        Poly a = randomMonicPoly(poly);
        if (a.isConstant())
            return null;

        Poly gcd1 = UnivariateGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant())
            return gcd1;

        // (modulus^d - 1) / 2
        BigInteger exponent = poly.coefficientDomainCardinality().pow(d).decrement().shiftRight(1);
        Poly b = PolynomialArithmetics.polyPowMod(a, exponent, poly, DivisionWithRemainder.fastDivisionPreConditioning(poly), true);

        Poly gcd2 = UnivariateGCD.PolynomialGCD(b.decrement(), poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }

    /**
     * Plain Cantor-Zassenhaus algorithm in GF2p
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <Poly extends IUnivariatePolynomial<Poly>> Poly CantorZassenhausGF2p(Poly poly, int d, int pPower) {
        assert poly.isMonic();

        Poly a = randomMonicPoly(poly);
        if (a.isConstant())
            return null;

        Poly gcd1 = UnivariateGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant())
            return gcd1;

        DivisionWithRemainder.InverseModMonomial<Poly> invMod = DivisionWithRemainder.fastDivisionPreConditioning(poly);
        Poly b = tracePolyGF2(a, LongArithmetics.safeToInt(1L * pPower * d), poly, invMod);
        Poly gcd2 = UnivariateGCD.PolynomialGCD(b.decrement(), poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }

    static <Poly extends IUnivariatePolynomial<Poly>> Poly tracePolyGF2(Poly a, int m, Poly modulus, DivisionWithRemainder.InverseModMonomial<Poly> invMod) {
        Poly result = a.createZero();
        for (int i = 0; i < m; i++)
            result.add(PolynomialArithmetics.polyPowMod(a, BigInteger.valueOf(1).shiftLeft(i), modulus, invMod, true));
        return result;
    }
}
