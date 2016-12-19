package cc.r2.core.number;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SieveOfAtkin;
import cc.r2.core.number.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.BitSet;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.ChineseRemainderAlgorithm.CRT;

/**
 * Created by poslavsky on 06/12/2016.
 */
public class ChineseRemainderAlgorithmTest {
    @Test
    public void test1() throws Exception {
        BigInteger[] coprimes = {BigInteger.valueOf(11), BigInteger.valueOf(13)};
        BigInteger[] remainders = {BigInteger.valueOf(2123), BigInteger.valueOf(7213)};
        BigInteger crt = CRT(coprimes, remainders);
        BigInteger sup = Arrays.stream(coprimes).reduce(ONE, BigInteger::multiply);
        System.out.println(crt);

        assertCRT(coprimes, remainders, crt);

        crt = toSymMod(crt, sup);
        System.out.println(crt);
        assertCRT(coprimes, remainders, crt);
    }


    @Test
    public void test2() throws Exception {
        BigInteger[] coprimes = {BigInteger.valueOf(121), BigInteger.valueOf(13)};
        BigInteger[] remainders = {BigInteger.valueOf(2123), BigInteger.valueOf(7213)};
        BigInteger crt = CRT(coprimes, remainders);
        BigInteger sup = Arrays.stream(coprimes).reduce(ONE, BigInteger::multiply);
        System.out.println(crt);

        assertCRT(coprimes, remainders, crt);

        crt = toSymMod(crt, sup);
        System.out.println(crt);
        assertCRT(coprimes, remainders, crt);
    }

    @Test
    public void testRandom() throws Exception {
        RandomGenerator rnd = new Well1024a(System.currentTimeMillis());
        int bound = 10000;
        SieveOfAtkin primes = SieveOfAtkin.createSieve(bound);
        for (int i = 0; i < 100; i++) {
            int n = 1 + rnd.nextInt(100);
            BitSet used = new BitSet(bound);
            BigInteger[] comprimes = new BigInteger[n];
            BigInteger[] remainders = new BigInteger[n];
            for (int j = 0; j < n; j++) {
                int prime;
                while (!primes.isPrime(prime = 1 + rnd.nextInt(bound - 1)) || used.get(prime)) ;
                used.set(prime);
                comprimes[j] = BigInteger.valueOf(prime);
                Assert.assertTrue(SmallPrimes.isPrime(prime));
                remainders[j] = BigInteger.valueOf(rnd.nextInt(prime));
            }

            BigInteger crt = CRT(comprimes, remainders);
            for (int j = 0; j < n; j++) {
                Assert.assertEquals(crt.mod(comprimes[j]), remainders[j]);
            }
        }

    }

    static void assertCRT(BigInteger[] coprimes, BigInteger[] rems, BigInteger crt) {
        for (int i = 0; i < coprimes.length; i++)
            Assert.assertEquals(rems[i].mod(coprimes[i]), crt.mod(coprimes[i]));
    }

    static BigInteger toSymMod(BigInteger b, BigInteger prime) {
        BigInteger t = prime.decrement().divide(BigInteger.TWO);
        if (b.compareTo(t) <= 0)
            return b;
        else return t.subtract(b);
    }
}