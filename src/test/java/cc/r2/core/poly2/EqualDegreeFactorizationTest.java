package cc.r2.core.poly2;

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
        MutablePolynomialMod a = MutablePolynomialZ.create(5224, 5225, 5225, 5225, 1).modulus(modulus);
        for (int i = 0; i < 10; i++)
            Assert.assertEquals(4, EqualDegreeFactorization.CantorZassenhaus(a, 1).size());
    }
}