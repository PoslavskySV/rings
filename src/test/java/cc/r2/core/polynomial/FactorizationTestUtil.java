package cc.r2.core.polynomial;

import static cc.r2.core.polynomial.PolynomialArithmetics.polyPow;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
public final class FactorizationTestUtil {
    public FactorizationTestUtil() {}

    static void assertFactorization(MutablePolynomial poly, Factorization factorization) {
        MutablePolynomial r = MutablePolynomial.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(polyPow(factorization.factors[i], factorization.exponents[i], true));
        assertEquals(poly, r);
    }

    static void assertFactorization(MutablePolynomial poly, Factorization factorization, long modulus) {
        MutablePolynomial r = MutablePolynomial.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(PolynomialArithmetics.polyPowMod(factorization.factors[i], factorization.exponents[i], modulus, true), modulus);
        assertEquals(poly.clone().modulus(modulus), r);
    }

    static void assertDistinctDegreeFactorization(MutablePolynomial poly, Factorization factorization, long modulus) {
        MutablePolynomial r = MutablePolynomial.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(factorization.factors[i], modulus);
        assertEquals(poly.clone().modulus(modulus), r);
    }
}
