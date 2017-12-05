package cc.redberry.rings.poly.univar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.Util;
import org.apache.commons.math3.random.RandomGenerator;


/**
 * Equal-degree factorization of univariate polynomials over finite fields.
 *
 * @since 1.0
 */
public final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    @SuppressWarnings("unchecked")
    static <T extends IUnivariatePolynomial<T>> T randomMonicPoly(T factory) {
        RandomGenerator rnd = PrivateRandom.getRandom();
        int degree = Math.max(1, rnd.nextInt(2 * factory.degree() + 1));
        if (factory instanceof UnivariatePolynomialZp64) {
            UnivariatePolynomialZp64 fm = (UnivariatePolynomialZp64) factory;
            return (T) RandomUnivariatePolynomials.randomMonicPoly(degree, fm.ring.modulus, rnd);
        } else {
            UnivariatePolynomial fm = (UnivariatePolynomial) factory;
            return (T) RandomUnivariatePolynomials.randomMonicPoly(degree, fm.ring, rnd);
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm implementation
     *
     * @param input the polynomial
     * @param d     distinct degree
     * @return irreducible factor of {@code poly}
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> FactorDecomposition<Poly> CantorZassenhaus(Poly input, int d) {
        Util.ensureOverFiniteField(input);
        FactorDecomposition<Poly> result = FactorDecomposition.unit(input.lcAsPoly());
        if (!input.coefficientRingCardinality().testBit(0))
            //even characteristic => GF2p
            CantorZassenhaus(input, d, result, input.coefficientRingPerfectPowerExponent().intValueExact());
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
            poly = UnivariateDivision.quotient(poly, factor, false);
        }
    }

    /**
     * Plain Cantor-Zassenhaus algorithm for odd characteristic
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <Poly extends IUnivariatePolynomial<Poly>> Poly CantorZassenhaus0(Poly poly, int d) {
        assert poly.isMonic();

        Poly a = randomMonicPoly(poly);
        if (a.isConstant() || a.equals(poly))
            return null;

        Poly gcd1 = UnivariateGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant() && !gcd1.equals(poly))
            return gcd1;

        // (modulus^d - 1) / 2
        BigInteger exponent = poly.coefficientRingCardinality().pow(d).decrement().shiftRight(1);
        Poly b = UnivariatePolynomialArithmetic.polyPowMod(a, exponent, poly, UnivariateDivision.fastDivisionPreConditioning(poly), true);

        Poly gcd2 = UnivariateGCD.PolynomialGCD(b.decrement(), poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }

    /**
     * Plain Cantor-Zassenhaus algorithm in GF2p.
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    private static <Poly extends IUnivariatePolynomial<Poly>> Poly CantorZassenhausGF2p(Poly poly, int d, int pPower) {
        assert poly.isMonic();

        Poly a = randomMonicPoly(poly);
        if (a.isConstant() || a.equals(poly))
            return null;

        Poly gcd1 = UnivariateGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant() && !gcd1.equals(poly))
            return gcd1;

        UnivariateDivision.InverseModMonomial<Poly> invMod = UnivariateDivision.fastDivisionPreConditioning(poly);
        Poly b = tracePolyGF2(a, MachineArithmetic.safeToInt(1L * pPower * d), poly, invMod);

        Poly gcd2 = UnivariateGCD.PolynomialGCD(b, poly);
        if (!gcd2.isConstant() && !gcd2.equals(poly))
            return gcd2;

        return null;
    }

    static <Poly extends IUnivariatePolynomial<Poly>> Poly tracePolyGF2(Poly a, int m, Poly modulus, UnivariateDivision.InverseModMonomial<Poly> invMod) {
        Poly tmp = a.clone();
        Poly result = a.clone();
        for (int i = 0; i < (m - 1); i++) {
            tmp = UnivariatePolynomialArithmetic.polyMultiplyMod(tmp, tmp, modulus, invMod, false);
            result.add(tmp);
        }

        return result;
    }
}
