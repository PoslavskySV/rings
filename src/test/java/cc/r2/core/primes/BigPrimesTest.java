package cc.r2.core.primes;

import cc.r2.core.bigint.BigInteger;
import cc.r2.core.test.AbstractTest;
import cc.r2.core.test.TimeConsuming;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.List;

import static cc.r2.core.bigint.BigInteger.ZERO;
import static cc.r2.core.primes.BigPrimes.LucasPrimalityTest;
import static org.junit.Assert.*;

public class BigPrimesTest extends AbstractTest {
    @Test
    public void test1() throws Exception {
        assertFalse(BigPrimes.isPrime(ZERO));
        assertEquals(1, BigPrimes.primeFactors(ZERO).size());
    }

    @Test
    @TimeConsuming
    public void test2() throws Exception {
        int lucasK = 20;
        int dups = 0;
        RandomGenerator rnd = getRandom();
        long[] somePrimes = {15722669197L, 72062552321653L, 41543465813L, 10707835013L, 1631650208641L, 490247130077L, 32726648113L};
        for (int i = 0; i < 1000; i++) {
            for (long prime : somePrimes) {
                BigInteger n = BigInteger.valueOf(prime);
                if (!LucasPrimalityTest(n, lucasK, rnd))
                    ++dups;
                assertTrue(BigPrimes.isPrime(n));
            }
        }
        assertTrue(dups < 10);
    }

    @Test
    public void test3() throws Exception {
        for (int i = 0; i < 1000; i++)
            assertEquals(SmallPrimes.nextPrime(i), BigPrimes.nextPrime(BigInteger.valueOf(i)).intValue());
    }

    @Test
    public void randomTest1() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 500); i++) {
            BigInteger n = new BigInteger(1 + rnd.nextInt(60), rnd);

            boolean isPrime = BigPrimes.isPrime(n);
            List<BigInteger> factors = BigPrimes.primeFactors(n);

            if (isPrime)
                assertEquals(1, factors.size());
            assertPrimeFactors(n, factors);
        }
    }

    @Test
    @TimeConsuming
    public void randomTest2() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < 10; i++) {
            BigInteger n = BigInteger.valueOf(1 + rnd.nextInt(30));
            int its = 3;
            n = n.multiply(n);
            while (its-- > 0)
                n = n.multiply(n.nextProbablePrime());

            boolean isPrime = BigPrimes.isPrime(n);
            List<BigInteger> factors = BigPrimes.primeFactors(n);

            if (isPrime)
                assertEquals(1, factors.size());
            assertPrimeFactors(n, factors);
        }
    }

    static void assertPrimeFactors(BigInteger n, List<BigInteger> factors) {
        if (factors.size() == 1) {
            Assert.assertEquals(n, factors.get(0));
            return;
        }

        BigInteger r = BigInteger.ONE;
        for (BigInteger f : factors) {
            assertFalse(f.isOne());
            r = r.multiply(f);
        }
        assertEquals(n, r);
    }
}