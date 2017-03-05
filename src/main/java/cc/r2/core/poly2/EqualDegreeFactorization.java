package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;


/**
 * Equal-degree factorization of univariate polynomials over finite fields with single-precision coefficients.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    /** thread local instance of random */
    private static final ThreadLocal<RandomGenerator> CZ_ThreadLocalRandom
            = ThreadLocal.withInitial(() -> new Well1024a(0x7f67fcad528cfae9L));

    /**
     * Plain Cantor-Zassenhaus algorithm
     *
     * @param poly the monic polynomial
     * @param d    distinct degree
     * @return irreducible factor of {@code poly}
     */
    static MutablePolynomialMod CantorZassenhaus0(MutablePolynomialMod poly, int d) {
        assert poly.lc() == 1;

        MutablePolynomialMod a = RandomPolynomials.randomMonicPoly(poly.degree - 1, poly.modulus, CZ_ThreadLocalRandom.get());
        if (a.isConstant())
            return null;

        MutablePolynomialMod gcd1 = PolynomialGCD.PolynomialGCD(a, poly);
        if (!gcd1.isConstant())
            return gcd1;

        // (modulus^d - 1) / 2
        BigInteger exponent = BigInteger.valueOfUnsigned(poly.modulus).pow(d).decrement().shiftRight(1);
        MutablePolynomialMod b = PolynomialArithmetics.polyPowMod(a, exponent, poly, fastDivisionPreConditioning(poly), true);

        MutablePolynomialMod gcd2 = PolynomialGCD.PolynomialGCD(b.decrement(), poly);
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
    static FactorDecomposition<MutablePolynomialMod> CantorZassenhaus(MutablePolynomialMod input, int d) {
        assert input.degree % d == 0;
        int nFactors = input.degree / d;
        if (input.degree == 1 || nFactors == 1)
            return FactorDecomposition.oneFactor(input, 1);

        assert input.degree != d;

        FactorDecomposition<MutablePolynomialMod> result = new FactorDecomposition<>();

        long overallFactor = input.lc();
        MutablePolynomialMod poly = input.clone();
        while (true) {

            poly = poly.monic();
            if (poly.degree == d) {
                result.addFactor(poly, 1);
                break;
            }

            MutablePolynomialMod factor;
            do {
                factor = CantorZassenhaus0(poly, d);
            } while (factor == null);

            if (factor.degree != d)
                result.addAll(CantorZassenhaus(factor, d));
            else
                result.addFactor(factor.monic(), 1);
            poly = DivisionWithRemainder.quotient(poly, factor, false);
        }

        return result.setNumericFactor(overallFactor);
    }
}
