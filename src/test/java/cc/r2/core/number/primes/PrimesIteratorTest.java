package cc.r2.core.number.primes;

import org.junit.Assert;
import org.junit.Test;

/**
 * Created by poslavsky on 21/12/2016.
 */
public class PrimesIteratorTest {

    @Test
    public void test1() throws Exception {
        PrimesIterator it = new PrimesIterator();
        int v = 0;
        while (v < 100) {
            long take = it.take();
            if (take == -1)
                break;
            Assert.assertEquals(SmallPrimes.SmallPrimes9[v++], take);
        }
    }

    @Test
    public void test2() throws Exception {
        PrimesIterator it = new PrimesIterator(Integer.MAX_VALUE - 1000);
        int v = 0;
        while (++v < 1000) {
            long take = it.take();
            if (take == -1)
                break;
            v++;
            Assert.assertTrue(BigPrimes.isPrime(take));
        }
    }

    @Test
    public void test3() throws Exception {
        Assert.assertEquals(7, new PrimesIterator(6).take());
        Assert.assertEquals(7, new PrimesIterator(7).take());
        Assert.assertEquals(4099, new PrimesIterator(4096).take());
        Assert.assertEquals(4111, new PrimesIterator(4099 + 1).take());
        Assert.assertEquals(1677721, new PrimesIterator(1677721).take());
        Assert.assertEquals(1677727, new PrimesIterator(1677722).take());
        Assert.assertEquals(1677727, new PrimesIterator(1677722).take());
        Assert.assertEquals(16777213, new PrimesIterator(16777210).take());
        Assert.assertEquals(16777259, new PrimesIterator(16777214).take());
        Assert.assertEquals(17777239, new PrimesIterator(17777214).take());
        Assert.assertEquals(1771277227, new PrimesIterator(1771277214).take());
    }
}