package cc.r2.core.number;

import cc.r2.core.test.AbstractTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.Test;

import static cc.r2.core.number.BigIntegerArithmetics.sqrtCeil;
import static cc.r2.core.number.BigIntegerArithmetics.sqrtFloor;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class BigIntegerArithmeticsTest extends AbstractTest {
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
        assertEquals(BigInteger.valueOf(3), sqrtFloor(BigInteger.valueOf(10)));
        assertEquals(BigInteger.valueOf(4), sqrtCeil(BigInteger.valueOf(10)));
    }

    @Test
    public void testSqrt3() throws Exception {
        assertEquals(BigInteger.valueOf(979), sqrtFloor(BigInteger.valueOf(979 * 979 + 1)));
        assertEquals(BigInteger.valueOf(980), sqrtCeil(BigInteger.valueOf(979 * 979 + 1)));
    }

    static void assertExactSqrt(long val) {
        BigInteger bval = BigInteger.valueOf(val);
        assertExactSqrt(bval);
    }

    static void assertExactSqrt(BigInteger bval) {
        assertEquals(bval, sqrtFloor(bval.pow(2)));
        assertEquals(bval, sqrtCeil(bval.pow(2)));
    }
}