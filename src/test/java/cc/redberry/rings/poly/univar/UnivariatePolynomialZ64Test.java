package cc.redberry.rings.poly.univar;

import org.junit.Assert;
import org.junit.Test;

/**
 * @since 1.0
 */
public class UnivariatePolynomialZ64Test extends AUnivariateTest {
    @Test
    public void test1() throws Exception {
        Assert.assertEquals(-1, UnivariatePolynomialZ64.create(0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(-1, UnivariatePolynomialZ64.create(0, 0, 0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(0, UnivariatePolynomialZ64.create(1).firstNonZeroCoefficientPosition());
        Assert.assertEquals(1, UnivariatePolynomialZ64.create(0, 1).firstNonZeroCoefficientPosition());
    }
}