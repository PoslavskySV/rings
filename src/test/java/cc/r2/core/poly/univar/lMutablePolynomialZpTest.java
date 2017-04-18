package cc.r2.core.poly.univar;

import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lMutablePolynomialZpTest {
    @Test
    public void test1() throws Exception {
        lMutablePolynomialZp aL = lMutablePolynomialZ.create(1, 2, 3, 4, 5, 6).modulus(59);
        for (int i = 0; i < 5; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            Assert.assertTrue(check(aL));
        }
    }

    @Test
    public void test2() throws Exception {
        lMutablePolynomialZp factory = lMutablePolynomialZp.zero(3);
        Assert.assertEquals(0, factory.negateMod(0));
        Assert.assertEquals(0, factory.negate().lc());
    }

    private static boolean check(lMutablePolynomialZp poly) {
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] >= poly.modulus)
                return false;
        }
        return true;
    }
}