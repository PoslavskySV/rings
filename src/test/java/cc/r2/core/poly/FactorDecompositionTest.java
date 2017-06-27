package cc.r2.core.poly;

import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class FactorDecompositionTest {

    public static <T extends IGeneralPolynomial<T>> void assertFactorization(T poly, FactorDecomposition<T> factorization) {
        assertEquals(poly, factorization.toPolynomial());
    }
}