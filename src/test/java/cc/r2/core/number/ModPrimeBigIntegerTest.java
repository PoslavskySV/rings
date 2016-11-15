package cc.r2.core.number;

import org.apache.commons.math3.random.Well512a;
import org.junit.Assert;
import org.junit.Test;


public class ModPrimeBigIntegerTest {
    @Test
    public void randomTest() throws Exception {
        long[] primes = {11, 113, 739, 10957, 1177921, 226317929};
        int N = 10000;
        Well512a rnd = new Well512a(System.currentTimeMillis());

        for (long prime : primes) {
            ModPrimeBigIntegerField ring = new ModPrimeBigIntegerField(BigInteger.valueOf(prime));
            for (int i = 0; i < N; i++) {
                ModPrimeBigInteger a = new ModPrimeBigInteger(ring, BigInteger.valueOf(rnd.nextLong()));
                ModPrimeBigInteger b = new ModPrimeBigInteger(ring, BigInteger.valueOf(rnd.nextLong()));
                if (rnd.nextBoolean() && rnd.nextBoolean())
                    a = a.negate();
                if (rnd.nextBoolean() && rnd.nextBoolean())
                    b = b.negate();

                if (a.isZero() || b.isZero())
                    continue;

                Assert.assertEquals(a, a.divide(b).multiply(b));
                Assert.assertEquals(b, b.divide(a).multiply(a));
            }
        }
    }

    @Test
    public void gedsf() throws Exception {
        System.out.println(BigInteger.valueOf(-123).mod(BigInteger.valueOf(12)));
        System.out.println(BigInteger.valueOf(123).mod(BigInteger.valueOf(12)));
//        ModPrimeBigIntegerField ring = new ModPrimeBigIntegerField(BigInteger.valueOf(11));
//        System.out.println(new ModPrimeBigInteger(ring, BigInteger.valueOf(-5)));
    }
}