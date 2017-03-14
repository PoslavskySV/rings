package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class bMutablePolynomialModTest {

    @Test
    public void test1() throws Exception {
        bMutablePolynomialZ aZ = bMutablePolynomialZ.create(1, 2, 3, 4, 5, 6);
        BigInteger modulus = BigInteger.valueOf(59);
        bMutablePolynomialMod a = aZ.modulus(modulus);
        MutablePolynomialMod aL = a.toLong();

        for (int i = 0; i < 5; i++) {
            a = (a.clone() * a.clone().decrement() - a.clone().derivative() + (a.clone().square())) * a.clone();
            a = a.shiftRight(2);
            aZ = (aZ.clone() * aZ.clone().decrement() - aZ.clone().derivative() + (aZ.clone().square())) * aZ.clone();

            aL = (aL.clone() * aL.clone().decrement() - aL.clone().derivative() + (aL.clone().square())) * aL.clone();
        }

        Assert.assertEquals(aL, a.toLong());
        Assert.assertEquals(a, aZ.modulus(modulus));
    }
}