package cc.r2.core.number;

import org.junit.Assert;
import org.junit.Test;

public class UtilTest {
    @Test
    public void name1() throws Exception {
        Assert.assertEquals(BigInteger.valueOf(3), Util.gcd(BigInteger.valueOf(6), BigInteger.valueOf(9), BigInteger.valueOf(12)));
        Assert.assertEquals(BigInteger.valueOf(36), Util.lcm(BigInteger.valueOf(6), BigInteger.valueOf(9), BigInteger.valueOf(12)));
    }
}