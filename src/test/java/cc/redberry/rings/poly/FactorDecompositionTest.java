package cc.redberry.rings.poly;

import static org.junit.Assert.assertEquals;

/**
 * @since 1.0
 */
public class FactorDecompositionTest {

    public static <T extends IPolynomial<T>> void assertFactorization(T poly, FactorDecomposition<T> factorization) {
        assertEquals(poly, factorization.toPolynomial());
    }
}