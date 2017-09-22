package cc.redberry.rings.bigint;

import cc.redberry.rings.test.AbstractTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @since 1.0
 */
public class BigIntegerUtilTest extends AbstractTest {
    @Test
    public void testSqrt1() throws Exception {
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < 100; i++) {
            long v = rndd.nextLong(0, Long.MAX_VALUE);
            if (v < 0) v = -v;
            assertExactSqrt(v);
        }
    }

    @Test
    public void testSqrt2() throws Exception {
        assertEquals(BigInteger.valueOf(3), BigIntegerUtil.sqrtFloor(BigInteger.valueOf(10)));
        assertEquals(BigInteger.valueOf(4), BigIntegerUtil.sqrtCeil(BigInteger.valueOf(10)));
    }

    @Test
    public void testSqrt3() throws Exception {
        assertEquals(BigInteger.valueOf(979), BigIntegerUtil.sqrtFloor(BigInteger.valueOf(979 * 979 + 1)));
        assertEquals(BigInteger.valueOf(980), BigIntegerUtil.sqrtCeil(BigInteger.valueOf(979 * 979 + 1)));
    }

    static void assertExactSqrt(long val) {
        BigInteger bval = BigInteger.valueOf(val);
        assertExactSqrt(bval);
    }

    static void assertExactSqrt(BigInteger bval) {
        assertEquals(bval, BigIntegerUtil.sqrtFloor(bval.pow(2)));
        assertEquals(bval, BigIntegerUtil.sqrtCeil(bval.pow(2)));
    }

    @Test
    public void testPerfectPower() throws Exception {
        int nIterations = its(1000, 1000);
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < nIterations; i++) {
            BigInteger base = BigInteger.valueOf(rndd.nextInt(2, 59));
            BigInteger exponent = BigInteger.valueOf(rndd.nextInt(2, 21));
            BigInteger n = BigIntegerUtil.pow(base, exponent);

            BigInteger[] ipp = BigIntegerUtil.perfectPowerDecomposition(n);
            assertNotNull(ipp);
            assertEquals(base + " ^ " + exponent, n, BigIntegerUtil.pow(ipp[0], ipp[1]));
            assertTrue(ipp[0] + " > " + base, ipp[0].compareTo(base) <= 0);
            assertTrue(ipp[1].compareTo(exponent) >= 0);

            if (!n.equals(BigInteger.valueOf(8))) // Catalan's conjecture
                Assert.assertNull(base + " ^ " + exponent, BigIntegerUtil.perfectPowerDecomposition(n.increment()));
        }
    }
}