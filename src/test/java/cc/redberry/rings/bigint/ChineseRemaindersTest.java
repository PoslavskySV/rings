package cc.redberry.rings.bigint;

import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.primes.SieveOfAtkin;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.test.AbstractTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.BitSet;

import static cc.redberry.rings.bigint.BigInteger.*;
import static cc.redberry.rings.bigint.ChineseRemainders.ChineseRemainders;
import static java.lang.Math.floorMod;

/**
 * Created by poslavsky on 06/12/2016.
 */
public class ChineseRemaindersTest extends AbstractTest {
    @Test
    public void test1() throws Exception {
        BigInteger[] coprimes = {BigInteger.valueOf(11), BigInteger.valueOf(13)};
        BigInteger[] remainders = {BigInteger.valueOf(2123), BigInteger.valueOf(7213)};
        BigInteger crt = ChineseRemainders.ChineseRemainders(coprimes, remainders);
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
        BigInteger crt = ChineseRemainders.ChineseRemainders(coprimes, remainders);
        BigInteger sup = Arrays.stream(coprimes).reduce(ONE, BigInteger::multiply);
        assertCRT(coprimes, remainders, crt);
        crt = toSymMod(crt, sup);
        assertCRT(coprimes, remainders, crt);
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
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
    public void test4() throws Exception {
        BigInteger crt;
        crt = ChineseRemainders.ChineseRemainders(BigInteger.valueOf(284407855036305L), BigInteger.valueOf(47), ONE, BigInteger.valueOf(-15));
        Assert.assertEquals(BigInteger.valueOf(8532235651089151L), crt);
        crt = ChineseRemainders.ChineseRemainders(BigInteger.valueOf(284407855036305L), BigInteger.valueOf(47), BigInteger.valueOf(-2), BigInteger.valueOf(-17));
        Assert.assertEquals(BigInteger.valueOf(9669867071234368L), crt);
        crt = ChineseRemainders.ChineseRemainders(BigInteger.valueOf(284407855036305L), BigInteger.valueOf(47), BigInteger.valueOf(2), BigInteger.valueOf(17));
        Assert.assertEquals(BigInteger.valueOf(3697302115471967L), crt);
    }

    @Test
    public void testRandom5() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics fast = new DescriptiveStatistics(), gen = new DescriptiveStatistics();
        long nIterations = its(1000, 10000);
        int nSmallIterations = 100;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                fast.clear();
                gen.clear();
            }

            BigInteger base = new BigInteger(rndd.nextInt(3, 48), rnd);
            BigInteger bPrime1 = (base = base.nextProbablePrime());
            BigInteger bPrime2 = (base = base.nextProbablePrime());
            long prime1 = bPrime1.longValueExact();
            long prime2 = bPrime2.longValueExact();

            if (MachineArithmetic.isOverflowMultiply(prime1, prime2)) {
                --i;
                continue;
            }

            if (prime1 < 10 && prime2 < 10) {
                --i;
                continue;
            }


            long start = System.nanoTime();
            ChineseRemainders.ChineseRemaindersMagic magic = ChineseRemainders.createMagic(prime1, prime2);
            fast.addValue(System.nanoTime() - start);
            for (int j = 0; j < nSmallIterations; j++) {
                long remainder1 = rndd.nextLong(1, prime1 - 1);
                long remainder2 = rndd.nextLong(1, prime2 - 1);
                //System.out.println(String.format("%s %s %s %s", prime1, prime2, remainder1, remainder2));

                start = System.nanoTime();
                long cdrt = ChineseRemainders.ChineseRemainders(magic, remainder1, remainder2);
                fast.addValue(System.nanoTime() - start);
                if (bPrime1.isInt() && bPrime2.isInt()) {
                    start = System.nanoTime();
                    long crt1 = ChineseRemainders.ChineseRemainders(prime1, prime2, remainder1, remainder2);
                    gen.addValue(System.nanoTime() - start);
                    Assert.assertEquals(crt1, cdrt);
                }
                Assert.assertEquals(
                        ChineseRemainders.ChineseRemainders(bPrime1, bPrime2, BigInteger.valueOf(remainder1), BigInteger.valueOf(remainder2)).longValue(),
                        cdrt);
            }
        }
        System.out.println("Fast: " + fast.getPercentile(50));
        System.out.println("Gene: " + gen.getPercentile(50));
    }

    @Test
    public void test5() throws Exception {
        long prime1 = 30223, prime2 = 30241, remainder1 = 21175, remainder2 = 29739;
        ChineseRemainders.ChineseRemaindersMagic magic = ChineseRemainders.createMagic(prime1, prime2);
        Assert.assertEquals(
                ChineseRemainders.ChineseRemainders(prime1, prime2, remainder1, remainder2),
                ChineseRemainders.ChineseRemainders(magic, remainder1, remainder2));
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = getRandom();
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

            BigInteger crt = ChineseRemainders.ChineseRemainders(comprimes, remainders);
            for (int j = 0; j < n; j++) {
                Assert.assertEquals(crt.mod(comprimes[j]), remainders[j]);
            }
        }
    }

    @Test
    public void test6() throws Exception {
        UnivariatePolynomialZp64
                prime1 = UnivariatePolynomialZ64.create(1, 2, 3, 4).modulus(19),
                prime2 = UnivariatePolynomialZ64.create(1, 2, 1, 1, 1, 1).modulus(19);

//        System.out.println(Factorization.factor(prime1));
//        System.out.println(Factorization.factor(prime2));
//        if (true) return;

        UnivariatePolynomialZp64
                remainder1 = UnivariatePolynomialZ64.create(1, 2, 3).modulus(19),
                remainder2 = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5).modulus(19);


        UnivariateRing<UnivariatePolynomialZp64> domain = new UnivariateRing<>(prime1.createOne());
        UnivariatePolynomialZp64 crt = ChineseRemainders.ChineseRemainders(domain, prime1, prime2, remainder1, remainder2);
        assertCRT(domain,
                new UnivariatePolynomialZp64[]{prime1, prime2},
                new UnivariatePolynomialZp64[]{remainder1, remainder2},
                crt);
    }

    static void assertCRT(BigInteger[] coprimes, BigInteger[] rems, BigInteger crt) {
        for (int i = 0; i < coprimes.length; i++)
            Assert.assertEquals(rems[i].mod(coprimes[i]), crt.mod(coprimes[i]));
    }

    static <E> void assertCRT(Ring<E> ring, E[] coprimes, E[] rems, E crt) {
        for (int i = 0; i < coprimes.length; i++)
            Assert.assertEquals("" + i, rems[i], ring.remainder(crt, coprimes[i]));
    }

    static BigInteger toSymMod(BigInteger value, BigInteger modulus) {
        if (modulus.compareTo(ZERO) < 0)
            throw new IllegalArgumentException("Negative modulus");
        value = value.mod(modulus);
        BigInteger mHalf = modulus.divide(TWO);
        return value.compareTo(mHalf) <= 0 ? value : value.subtract(modulus);
    }
}