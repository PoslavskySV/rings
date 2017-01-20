package cc.r2.core.polynomial;

import cc.r2.core.number.primes.SieveOfAtkin;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

public class MutablePolynomialTest {
    @Test
    public void testShift1() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5, 6);
        MutablePolynomial expected = MutablePolynomial.create(3, 4, 5, 6);
        Assert.assertEquals(expected, a.clone().shiftLeft(2));
    }

    @Test
    public void testShift2() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 5, 6);
        MutablePolynomial expected = MutablePolynomial.create(0, 0, 1, 2, 3, 4, 5, 6);
        Assert.assertEquals(expected, a.clone().shiftRight(2));
    }


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

    @Test
    public void testSquare1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        testSquareRandom(rnd, 15009, 300, 301, 3);
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
            long[] a = RandomPolynomials.randomLongArray(aDegree, bound, rnd);
            long[] b = RandomPolynomials.randomLongArray(bDegree, bound, rnd);

            long[] classical;
            long start = System.nanoTime();
            if (modulus != -1)
                classical = MutablePolynomial.multiplyModClassical(a, 0, a.length, b, 0, b.length, modulus);
            else
                classical = MutablePolynomial.multiplyClassical(a, 0, a.length, b, 0, b.length);
            classicalTimes.addValue(System.nanoTime() - start);

            long[] karatsuba;
            start = System.nanoTime();
            if (modulus != -1)
                karatsuba = MutablePolynomial.multiplyModKaratsuba(a, 0, a.length, b, 0, b.length, modulus);
            else
                karatsuba = MutablePolynomial.multiplyKaratsuba(a, 0, a.length, b, 0, b.length);
            karatsubaTimes.addValue(System.nanoTime() - start);

            MutablePolynomial polyA = MutablePolynomial.create(a.clone()), polyB = MutablePolynomial.create(b);
            long[] asis = modulus != -1 ? polyA.multiply(polyB, modulus).data : polyA.multiply(polyB).data;

            Assert.assertArrayEquals(classical, karatsuba);
            if(!Arrays.equals(classical,asis)) {
                System.out.println(polyA);
                System.out.println(polyB);
            }
            Assert.assertArrayEquals(classical, asis);
            if (modulus != -1) {
                for (int j = 0; j < classical.length; j++) {
                    Assert.assertTrue(asis[j] < modulus && asis[j] >= 0);
                    Assert.assertTrue(classical[j] < modulus && classical[j] >= 0);
                    Assert.assertTrue(karatsuba[j] < modulus && karatsuba[j] >= 0);
                }
            }
        }

        System.out.println("========== Classical ===========");
        System.out.println(classicalTimes);
        System.out.println("========== Karatsuba ===========");
        System.out.println(karatsubaTimes);
    }


    private static void testSquareRandom(RandomGenerator rnd,
                                         int nIterations,
                                         int minADegree, int maxADegree,
                                         int bound) {
        testSquareModRandom(rnd, nIterations, minADegree, maxADegree, bound, -1);
    }

    private static void testSquareModRandom(RandomGenerator rnd,
                                            int nIterations,
                                            int minADegree, int maxADegree,
                                            int bound, long modulus) {
        DescriptiveStatistics
                classicalTimes = new DescriptiveStatistics(),
                karatsubaTimes = new DescriptiveStatistics(),
                multClassicalTimes = new DescriptiveStatistics(),
                multKaratsubaTimes = new DescriptiveStatistics();
        ;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                classicalTimes.clear();
                karatsubaTimes.clear();
                multClassicalTimes.clear();
                multKaratsubaTimes.clear();
            }
            int aDegree = minADegree + rnd.nextInt(maxADegree - minADegree);
            long[] a = RandomPolynomials.randomLongArray(aDegree, bound, rnd);

            for (int j = 0; j < a.length; j++)
                if (a[j] == 0) a[j] = 1;

            long[] classical;
            long start = System.nanoTime();
            if (modulus != -1)
                classical = MutablePolynomial.squareModClassical(a, 0, a.length, modulus);
            else
                classical = MutablePolynomial.squareClassical(a, 0, a.length);
            classicalTimes.addValue(System.nanoTime() - start);

            long[] karatsuba;
            start = System.nanoTime();
            if (modulus != -1)
                karatsuba = MutablePolynomial.squareModKaratsuba(a, 0, a.length, modulus);
            else
                karatsuba = MutablePolynomial.squareKaratsuba(a, 0, a.length);
            karatsubaTimes.addValue(System.nanoTime() - start);


            long[] multClassical;
            start = System.nanoTime();
            if (modulus != -1)
                multClassical = MutablePolynomial.multiplyModClassical(a, 0, a.length, a, 0, a.length, modulus);
            else
                multClassical = MutablePolynomial.multiplyClassical(a, 0, a.length, a, 0, a.length);
            multClassicalTimes.addValue(System.nanoTime() - start);


            long[] multKaratsuba;
            start = System.nanoTime();
            if (modulus != -1)
                multKaratsuba = MutablePolynomial.multiplyModKaratsuba(a, 0, a.length, a, 0, a.length, modulus);
            else
                multKaratsuba = MutablePolynomial.multiplyKaratsuba(a, 0, a.length, a, 0, a.length);
            multKaratsubaTimes.addValue(System.nanoTime() - start);

            Assert.assertArrayEquals(classical, karatsuba);
            Assert.assertArrayEquals(classical, multClassical);
            Assert.assertArrayEquals(classical, multKaratsuba);
        }

        System.out.println("========== Classical square ===========");
        System.out.println(classicalTimes);
        System.out.println("========== Karatsuba square ===========");
        System.out.println(karatsubaTimes);
        System.out.println("========== Classical multiplication ===========");
        System.out.println(multClassicalTimes);
        System.out.println("========== Karatsuba multiplication ===========");
        System.out.println(multKaratsubaTimes);
    }

    @Test
    public void testDerivative1() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 4, 0);
        Assert.assertEquals(MutablePolynomial.create(2, 6, 12), MutablePolynomial.derivative(a));
    }
}