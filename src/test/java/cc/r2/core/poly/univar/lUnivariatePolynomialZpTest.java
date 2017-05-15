package cc.r2.core.poly.univar;

import cc.r2.core.poly.AbstractPolynomialTest;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lUnivariatePolynomialZpTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        lUnivariatePolynomialZp aL = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, 6).modulus(59);
        for (int i = 0; i < 5; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            Assert.assertTrue(check(aL));
        }
    }

    @Test
    public void test2() throws Exception {
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(3);
        Assert.assertEquals(0, factory.domain.negate(0));
        Assert.assertEquals(0, factory.negate().lc());
    }

    @Test
    public void test4() throws Exception {
        System.out.println(lUnivariatePolynomialZ.create(0).firstNonZeroCoefficientPosition());
    }

    private static boolean check(lUnivariatePolynomialZp poly) {
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] >= poly.domain.modulus)
                return false;
        }
        return true;
    }
}