package cc.redberry.rings.poly;

import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.util.RandomDataGenerator;
import org.junit.Assert;
import org.junit.Test;

/**
 * @since 1.0
 */
public class MachineArithmeticTest extends AbstractTest {
    @Test
    public void testPerfectPower() throws Exception {
        int nIterations = its(1000, 1000);
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < nIterations; i++) {
            long base = rndd.nextInt(2, 39);
            long exponent = rndd.nextInt(2, 11);
            long n = MachineArithmetic.safePow(base, exponent);

            long[] ipp = MachineArithmetic.perfectPowerDecomposition(n);
            Assert.assertEquals(n, MachineArithmetic.safePow(ipp[0], ipp[1]));
            Assert.assertTrue(ipp[0] + " > " + base, ipp[0] <= base);
            Assert.assertTrue(ipp[1] >= exponent);

            if (n != 8) // Catalan's conjecture
                Assert.assertNull(base + " ^ " + exponent, MachineArithmetic.perfectPowerDecomposition(n + 1));
        }
    }
}