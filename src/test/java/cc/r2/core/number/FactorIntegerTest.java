package cc.r2.core.number;

import org.apache.commons.math3.random.Well1024a;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.FactorInteger.PollardRho;
import static cc.r2.core.number.FactorInteger.TrialDivision;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 11/11/2016.
 */
public class FactorIntegerTest {

    @Test
    public void asdasd() throws Exception {
//
//        System.out.println(213124311234L % 9L);
//        BigInteger prime = BigInteger.ONE;
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//        System.out.println(prime = prime.nextProbablePrime());
//
//        System.out.println(prime = prime.nextProbablePrime());


//        System.out.println(PollardP1(BigInteger.valueOf(1003 * 1039), TWO, BigInteger.valueOf(10000)));
        BigInteger integer = new BigInteger("1844674073797").multiply(new BigInteger("1844674073173"));
        System.out.println(integer);
        System.out.println(integer.toString().length());
//        System.out.println(TrialDivision(integer));
        BigInteger x = PollardRho(integer, 7);
//        BigInteger x = PollardP1x(integer, ONE, TWO, TWO, BigInteger.valueOf(5000));
        System.out.println(x);
        System.out.println(Arrays.toString(integer.divideAndRemainder(x)));

//        System.out.println(pMinusOneFactor(BigInteger.valueOf(213119)));
    }

    @Test
    public void testTrialDivision1() throws Exception {
        BigInteger integer = BigInteger.valueOf(2 * 2 * 3 * 17 * 19 * 19 * 371);
        assertFactorization(integer, TrialDivision(integer));
    }

    @Ignore
    @Test
    public void testTrialDivision2() throws Exception {
        BigInteger integer = new BigInteger("32148653256391487569");
        List<BigInteger> factors = TrialDivision(integer);
        assertFactorization(integer, factors);
    }

    @Test
    public void testTrialDivision3() throws Exception {
        BigInteger integer = BigInteger.valueOf(2036166157L * 3001L * 3003L);

//        System.out.println(BigInteger.valueOf(2471951).pow(2));
        System.out.println(integer);
        List<BigInteger> factors = TrialDivision(integer);
        System.out.println(factors);
        assertFactorization(integer, factors);
    }

    @Test
    public void testTrialDivisionRandomComposite1() throws Exception {
        Well1024a rnd = new Well1024a(System.currentTimeMillis());
        for (int i = 0; i < 1000; i++) {
            BigInteger num = BigInteger.ONE;
            for (int j = 0; j < 10; j++)
                num = num.multiply(BigInteger.valueOf(rnd.nextInt(1024)));

            List<BigInteger> factors = TrialDivision(num);
            assertFactorization(num, factors);
        }
    }

    @Test
    public void testTrialDivisionRandomPrime() throws Exception {
        BigInteger[] primes = eulerPrimes();
        for (BigInteger num : primes)
            assertFactorization(num, TrialDivision(num));
    }

    static BigInteger[] eulerPrimes() {
        BigInteger[] primes = new BigInteger[40];
        BigInteger bi41 = BigInteger.valueOf(41);
        for (int i = 0; i < 40; i++) {
            BigInteger bi = BigInteger.valueOf(i);
            primes[i] = bi.multiply(bi).add(bi).add(bi41);
        }
        return primes;
    }

    static void assertFactorization(BigInteger integer, List<BigInteger> factors) {
        if (integer.isZero())
            assertTrue(factors.isEmpty());
        else {
            BigInteger actual = ONE;
            for (BigInteger factor : factors) {
                assertTrue(factor.isPrime());
                actual = actual.multiply(factor);
            }
            assertEquals(integer, actual);
        }
    }
}