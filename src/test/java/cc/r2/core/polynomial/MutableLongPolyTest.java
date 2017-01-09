package cc.r2.core.polynomial;

import cc.r2.core.number.primes.SieveOfAtkin;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

public class MutableLongPolyTest {
    @Test
    public void testMultiply1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        testMultiplyRandom(rnd, 10000, 0, 100, 0, 100, 3);
    }

    @Test
    public void testMultiply2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        testMultiplyRandom(rnd, 1000, 100, 500, 100, 500, 3);
    }

    @Test
    public void testMultiply3() throws Exception {
        RandomGenerator rnd = new Well1024a();
        testMultiplyRandom(rnd, 1000, 1000, 5000, 1000, 5000, 3);
    }

    @Test
    public void testMultiplyMod4() throws Exception {
        RandomGenerator rnd = new Well1024a();
        testMultiplyModRandom(rnd, 1000, 0, 100, 0, 100, Integer.MAX_VALUE, 2);
        testMultiplyModRandom(rnd, 1000, 0, 100, 0, 100, Integer.MAX_VALUE, 3);
        testMultiplyModRandom(rnd, 1000, 0, 100, 0, 100, Integer.MAX_VALUE, 5);
    }

    @Test
    public void testMultiplyMod5() throws Exception {
        RandomGenerator rnd = new Well1024a();
        SieveOfAtkin sieve = SieveOfAtkin.createSieve(10000);
        for (int i = 0; i < 5; i++) {
            int modulus;
            do {modulus = sieve.randomPrime(rnd);} while (modulus <= 3);
            testMultiplyModRandom(rnd, 1000, 0, 100, 0, 100, Integer.MAX_VALUE, modulus);
            testMultiplyModRandom(rnd, 100, 100, 500, 100, 500, Integer.MAX_VALUE, modulus);
            testMultiplyModRandom(rnd, 10, 1000, 5000, 1000, 5000, Integer.MAX_VALUE, modulus);
        }
    }

    private static void testMultiplyRandom(RandomGenerator rnd,
                                           int nIterations,
                                           int minADegree, int maxADegree,
                                           int minBDegree, int maxBDegree,
                                           int bound) {
        testMultiplyModRandom(rnd, nIterations, minADegree, maxADegree, minBDegree, maxBDegree, bound, -1);
    }

    private static void testMultiplyModRandom(RandomGenerator rnd,
                                              int nIterations,
                                              int minADegree, int maxADegree,
                                              int minBDegree, int maxBDegree,
                                              int bound, long modulus) {
        DescriptiveStatistics
                classicalTimes = new DescriptiveStatistics(),
                karatsubaTimes = new DescriptiveStatistics();
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                classicalTimes.clear();
                karatsubaTimes.clear();
            }
            int aDegree = minADegree + rnd.nextInt(maxADegree - minADegree);
            int bDegree = minBDegree + rnd.nextInt(maxBDegree - minBDegree);
            long[] a = RandomPolynomials.randomArray(aDegree, bound, rnd);
            long[] b = RandomPolynomials.randomArray(bDegree, bound, rnd);

            long[] classical;
            long start = System.nanoTime();
            if (modulus != -1)
                classical = MutableLongPoly.multiplyModClassical(a, 0, a.length, b, 0, b.length, modulus);
            else
                classical = MutableLongPoly.multiplyClassical(a, 0, a.length, b, 0, b.length);
            classicalTimes.addValue(System.nanoTime() - start);

            long[] karatsuba;
            start = System.nanoTime();
            if (modulus != -1)
                karatsuba = MutableLongPoly.multiplyModKaratsuba(a, 0, a.length, b, 0, b.length, modulus);
            else
                karatsuba = MutableLongPoly.multiplyKaratsuba(a, 0, a.length, b, 0, b.length);
            karatsubaTimes.addValue(System.nanoTime() - start);

            MutableLongPoly polyA = MutableLongPoly.create(a.clone()), polyB = MutableLongPoly.create(b);
            long[] asis = modulus != -1 ? polyA.multiply(polyB, modulus).data : polyA.multiply(polyB).data;

            Assert.assertArrayEquals(classical, karatsuba);
            Assert.assertArrayEquals(classical, asis);
        }

        System.out.println("========== Classical ===========");
        System.out.println(classicalTimes);
        System.out.println("========== Karatsuba ===========");
        System.out.println(karatsubaTimes);
    }
}