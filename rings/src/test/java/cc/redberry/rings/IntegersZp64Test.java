package cc.redberry.rings;

import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.TimeUnits;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

/**
 *
 */
public class IntegersZp64Test {

    @Test
    @Ignore
    public void test1() {
        for (int k = 0; k < 1000; ++k) {
            IntegersZp64 zp = new IntegersZp64(SmallPrimes.nextPrime(10000 + k));
            long start = System.nanoTime();
            long sum1 = 0;
            for (int i = 1; i < zp.modulus; ++i)
                sum1 += zp.reciprocal(i);
            long el1 = System.nanoTime() - start;

            zp.buildCachedReciprocals();
            start = System.nanoTime();
            long sum2 = 0;
            for (int i = 1; i < zp.modulus; ++i)
                sum2 += zp.reciprocal(i);
            long el2 = System.nanoTime() - start;

            System.out.println("no cache: " + TimeUnits.nanosecondsToString(el1));
            System.out.println("   cache: " + TimeUnits.nanosecondsToString(el2));
            System.out.println();
            Assert.assertEquals(sum1, sum2);
        }
    }
}