package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.*;

/**
 * @since 1.0
 */
public class RationalsTest {
    @Test
    public void test1() throws Exception {
        for (int i = 19; i < 1000; i++)
            Assert.assertNotNull(Frac(UnivariateRingZp64(17)).randomElement());
    }

    @Test
    public void test2() throws Exception {
        Rational<BigInteger> rat1 = new Rational<>(Z, Z.valueOf(2), Z.valueOf(3));
        Rational<BigInteger> rat2 = new Rational<>(Z, Z.valueOf(0), Z.valueOf(13));
        Assert.assertEquals(rat2, rat1.subtract(rat1));
        Rational<BigInteger> rat3 = new Rational<>(Z, Z.valueOf(2), Z.valueOf(2));
        Rational<BigInteger> rat4 = new Rational<>(Z, Z.valueOf(13), Z.valueOf(13));
        Assert.assertEquals(rat3, rat4);
    }
}