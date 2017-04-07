package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;


/**
 * Equal-degree factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    @SuppressWarnings("unchecked")
    private static <T extends IMutablePolynomialZp<T>> T randomMonicPoly(T factory) {
        if (factory instanceof lMutablePolynomialZp) {
            lMutablePolynomialZp fm = (lMutablePolynomialZp) factory;
            return (T) RandomPolynomials.randomMonicPoly(fm.degree - 1, fm.modulus, GlobalRandom.getRandom());
        } else {
            bMutablePolynomialZp fm = (bMutablePolynomialZp) factory;
            return (T) RandomPolynomials.randomMonicPoly(fm.degree - 1, fm.modulus, GlobalRandom.getRandom());
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    static <T extends IMutablePolynomialZp<T>> T CantorZassenhaus0(T poly, int d) {
        assert poly.isMonic();

        T a = randomMonicPoly(poly);
        if (a.isConstant())
            return null;

        T gcd1 = PolynomialGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant())
            return gcd1;

        // (modulus^d - 1) / 2
        BigInteger exponent = poly.modulusAsBigInt().pow(d).decrement().shiftRight(1);
        T b = PolynomialArithmetics.polyPowMod(a, exponent, poly, DivisionWithRemainder.fastDivisionPreConditioning(poly), true);

        T gcd2 = PolynomialGCD.PolynomialGCD(b.decrement(), poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    static lFactorDecomposition<lMutablePolynomialZp> CantorZassenhaus(lMutablePolynomialZp input, int d) {
        lFactorDecomposition<lMutablePolynomialZp> result = new lFactorDecomposition<>();
        CantorZassenhaus(input, d, result);
        return result.setNumericFactor(input.lc());
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    static bFactorDecomposition<bMutablePolynomialZp> CantorZassenhaus(bMutablePolynomialZp input, int d) {
        bFactorDecomposition<bMutablePolynomialZp> result = new bFactorDecomposition<>();
        CantorZassenhaus(input, d, result);
        return result.setNumericFactor(input.lc());
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    @SuppressWarnings("unchecked")
    static <T extends IMutablePolynomialZp<T>> FactorDecomposition<T> CantorZassenhaus(T input, int d) {
        if (input instanceof lMutablePolynomialZp)
            return (FactorDecomposition<T>) CantorZassenhaus((lMutablePolynomialZp) input, d);
        else
            return (FactorDecomposition<T>) CantorZassenhaus((bMutablePolynomialZp) input, d);
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <T extends IMutablePolynomialZp<T>> void CantorZassenhaus(T input, int d, FactorDecomposition<T> result) {
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
}
