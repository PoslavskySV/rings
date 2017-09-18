package cc.r2.core.primes;

import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well512a;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import static cc.r2.core.primes.SmallPrimes.*;
import static org.junit.Assert.*;

public class SmallPrimesTest {
    @Test
    public void test1() throws Exception {
        for (int i : SmallPrimes12)
            assertTrue(SmallPrimes.isPrime(i));
    }

    @Test
    public void test2() throws Exception {
        assertEquals(SmallPrimes12.length, new TIntHashSet(SmallPrimes12).size());
        assertEquals(SmallPrimes10.length, new TIntHashSet(SmallPrimes10).size());
    }

    @Test
    public void test3() throws Exception {
        SieveOfAtkin sieve = SieveOfAtkin.SmallPrimesSieve;
        for (int i = 0; i < sieve.getLimit(); i++) {
            if (sieve.isPrime(i))
                assertTrue(SmallPrimes.isPrime(i));
            else
                assertPrimeFactors(i, SmallPrimes.primeFactors(i));
        }
    }

    @Test
    public void test4() throws Exception {
        assertFalse(SmallPrimes.isPrime(0));
        assertEquals(1, SmallPrimes.primeFactors(0).length);
    }

    @Test
    public void test5() throws Exception {
        SieveOfAtkin sieve = SieveOfAtkin.createSieve(PRIMES12_LAST);
        int p = 0;
        for (int i = 0; i < PRIMES12_LAST; i++)
            if (sieve.isPrime(i))
                assertEquals(i, SmallPrimes12[p++]);
    }

    @Test
    public void randomTest() throws Exception {
        RandomGenerator rnd = new Well512a();
        for (int i = 0; i < 1000; i++) {
            int n = rnd.nextInt();
            if (n < 0) n = -n;

            boolean isPrime = SmallPrimes.isPrime(n);
            int[] factors = SmallPrimes.primeFactors(n);

            if (isPrime)
                assertEquals(1, factors.length);
            assertPrimeFactors(n, factors);
        }
    }

    static void assertPrimeFactors(int n, int[] factors) {
        int r = 1;
        for (int f : factors)
            r = r * f;
        assertEquals(n, r);
    }

    @Ignore
    @Test
    public void prettyPrint() throws Exception {
        System.out.println("static final int[] SmallPrimes12 = {");
        for (int i = 0; ; i++) {
            for (int j = 0; j < 20; ++j) {
                if (i >= SmallPrimes12.length) {
                    System.out.print("\n};");
                    return;
                }
                int prime = SmallPrimes12[i];
                System.out.print(prime + ", ");
                i++;
            }
            System.out.print("\n");
        }
    }

    @Test
    public void test6() throws Exception {
        Assert.assertEquals(Integer.MAX_VALUE, SmallPrimes.nextPrime(Integer.MAX_VALUE));
    }
}