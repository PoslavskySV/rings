package cc.r2.core.poly.univar;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lUnivariatePolynomialZTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        Assert.assertEquals(-1, lUnivariatePolynomialZ.create(0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(-1, lUnivariatePolynomialZ.create(0, 0, 0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(0, lUnivariatePolynomialZ.create(1).firstNonZeroCoefficientPosition());
        Assert.assertEquals(1, lUnivariatePolynomialZ.create(0, 1).firstNonZeroCoefficientPosition());
    }
}