package cc.r2.core.poly.univar2;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lUnivariatePolynomialZTest {
    @Test
    public void test1() throws Exception {
        Assert.assertEquals(-1, lUnivariatePolynomialZ.create(0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(-1, lUnivariatePolynomialZ.create(0, 0, 0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(0, lUnivariatePolynomialZ.create(1).firstNonZeroCoefficientPosition());
        Assert.assertEquals(1, lUnivariatePolynomialZ.create(0, 1).firstNonZeroCoefficientPosition());
    }
}