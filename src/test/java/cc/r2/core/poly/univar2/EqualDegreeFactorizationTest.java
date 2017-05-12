package cc.r2.core.poly.univar2;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class EqualDegreeFactorizationTest {

    @Test
    public void test1() throws Exception {
        int modulus = 6101;
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(5224, 5225, 5225, 5225, 1).modulus(modulus);
        for (int i = 0; i < 10; i++)
            Assert.assertEquals(4, EqualDegreeFactorization.CantorZassenhaus(a, 1).size());
    }

    @Test
    public void test2() throws Exception {
        int modulus = 13;
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(9, 0, 1).modulus(modulus);
        Assert.assertEquals(2, EqualDegreeFactorization.CantorZassenhaus(poly, 1).size());
    }
}