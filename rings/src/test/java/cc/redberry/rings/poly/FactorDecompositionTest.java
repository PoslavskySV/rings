package cc.redberry.rings.poly;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MultivariateFactorization;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @since 1.0
 */
public class FactorDecompositionTest {

    public static <T extends IPolynomial<T>> void assertFactorization(T poly, PolynomialFactorDecomposition<T> factorization) {
        assertEquals(poly, factorization.multiply());
    }

    @Test
    public void test1() {
        MultivariatePolynomial<BigInteger> p = MultivariatePolynomial.parse("x^3 + x^2*y");
        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.Factor(p);

        MultivariatePolynomial<BigInteger> pOrig = p.clone();
        for (int i = 0; i < 3; ++i) {
            factors.canonical();
            assertEquals(pOrig, factors.multiply());
        }
    }
}