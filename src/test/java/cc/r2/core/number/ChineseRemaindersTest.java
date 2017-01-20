package cc.r2.core.number;

import cc.r2.core.number.primes.SieveOfAtkin;
import cc.r2.core.number.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.BitSet;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.number.ChineseRemainders.CRT;
import static cc.r2.core.number.ChineseRemainders.ChineseRemainders;
import static java.lang.Math.floorMod;

/**
 * Created by poslavsky on 06/12/2016.
 */
public class ChineseRemaindersTest {
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
        assertCRT(coprimes, remainders, crt);
        crt = toSymMod(crt, sup);
        assertCRT(coprimes, remainders, crt);
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a(System.currentTimeMillis());
        int bound = 100;
        SieveOfAtkin sieve = SieveOfAtkin.createSieve(bound);
        for (int i = 0; i < 1000; i++) {
            int n = 2 + rnd.nextInt(9);
            BitSet used = new BitSet(bound);
            long[] primes = new long[n];
            long[] remainders = new long[n];
            for (int j = 0; j < n; j++) {
                int prime;
                while (!sieve.isPrime(prime = 1 + rnd.nextInt(bound - 1)) || used.get(prime)) ;
                used.set(prime);
                primes[j] = prime;
                Assert.assertTrue(SmallPrimes.isPrime(prime));
                remainders[j] = rnd.nextInt(prime);
                remainders[j] = rnd.nextBoolean() ? -remainders[j] : remainders[j];
            }

            long crt = ChineseRemainders(primes, remainders);
            for (int j = 0; j < n; j++)
                Assert.assertEquals(floorMod(crt, primes[j]), floorMod(remainders[j], primes[j]));

            crt = ChineseRemainders(Arrays.copyOf(primes, 2), Arrays.copyOf(remainders, 2));
            Assert.assertEquals(crt, ChineseRemainders(primes[0], primes[1], remainders[0], remainders[1]));
        }
    }


    @Test
    public void test3() throws Exception {
        long crt;
        crt = ChineseRemainders.ChineseRemainders(284407855036305L, 47, 1, -15);
        Assert.assertEquals(8532235651089151L, crt);
        crt = ChineseRemainders.ChineseRemainders(284407855036305L, 47, -2, -17);
        Assert.assertEquals(9669867071234368L, crt);
        crt = ChineseRemainders.ChineseRemainders(284407855036305L, 47, 2, 17);
        Assert.assertEquals(3697302115471967L, crt);
    }

    @Test
    public void testRandomXXX() throws Exception {
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

    static BigInteger toSymMod(BigInteger value, BigInteger modulus) {
        if (modulus.compareTo(ZERO) < 0)
            throw new IllegalArgumentException("Negative modulus");
        value = value.mod(modulus);
        BigInteger mHalf = modulus.divide(TWO);
        return value.compareTo(mHalf) <= 0 ? value : value.subtract(modulus);
    }
}