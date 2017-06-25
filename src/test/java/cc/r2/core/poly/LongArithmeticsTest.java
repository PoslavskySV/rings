package cc.r2.core.poly;

import cc.r2.core.test.AbstractTest;
import cc.r2.core.util.RandomDataGenerator;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class LongArithmeticsTest extends AbstractTest {
    @Test
    public void testPerfectPower() throws Exception {
        int nIterations = its(1000, 1000);
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < nIterations; i++) {
            long base = rndd.nextInt(2, 39);
            long exponent = rndd.nextInt(2, 11);
            long n = LongArithmetics.safePow(base, exponent);

            long[] ipp = LongArithmetics.perfectPowerDecomposition(n);
            Assert.assertEquals(n, LongArithmetics.safePow(ipp[0], ipp[1]));
            Assert.assertTrue(ipp[0] + " > " + base, ipp[0] <= base);
            Assert.assertTrue(ipp[1] >= exponent);

            if (n != 8) // Catalan's conjecture
                Assert.assertNull(base + " ^ " + exponent, LongArithmetics.perfectPowerDecomposition(n + 1));
        }
    }
}