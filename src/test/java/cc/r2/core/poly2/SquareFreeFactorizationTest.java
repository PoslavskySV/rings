package cc.r2.core.poly2;

import cc.r2.core.poly2.FactorizationTestUtil.RandomSource;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import static cc.r2.core.poly2.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.poly2.PolynomialArithmetics.polyPow;
import static cc.r2.core.poly2.RandomPolynomials.randomPoly;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreeFactorization;
import static cc.r2.core.poly2.SquareFreeFactorization.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 20/01/2017.
 */
public class SquareFreeFactorizationTest {
    @Test
    public void test1() throws Exception {
        MutablePolynomialZ poly = polyPow(MutablePolynomialZ.create(1, 3).multiply(2), 3, false).multiply(polyPow(MutablePolynomialZ.create(-3, -5, 7), 2, false));
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomialZ.create(1, 3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomialZ.create(3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomialZ.create(33);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomialZ.create(22, 22).multiply(MutablePolynomialZ.create(12, 12, 12)).multiply(12);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics yun = new DescriptiveStatistics(), musser = new DescriptiveStatistics();

        long overflows = 0;
        for (int i = 0; i < 1000; i++) {
            if (i == 100) {
                yun.clear();
                musser.clear();
            }
            int nbase = rndd.nextInt(1, 5);
            MutablePolynomialZ poly;
            try {
                poly = MutablePolynomialZ.create(rndd.nextLong(1, 10));
                for (int j = 0; j < nbase; j++) {
                    MutablePolynomialZ factor = randomPoly(rndd.nextInt(1, 3), 10, rnd);
                    int exponent = rndd.nextInt(1, 5);
                    poly = poly.multiply(polyPow(factor, exponent, true));
                }
            } catch (ArithmeticException e) {
                --i;
                continue;
            }
            try {
                long start = System.nanoTime();
                FactorDecomposition<MutablePolynomialZ> yunFactorization = SquareFreeFactorizationYun(poly);
                yun.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                FactorDecomposition<MutablePolynomialZ> musserFactorization = SquareFreeFactorizationMusser(poly);
                musser.addValue(System.nanoTime() - start);

                assertEquals(yunFactorization.factors.size(), musserFactorization.factors.size());
                assertFactorization(poly, musserFactorization);
                assertFactorization(poly, yunFactorization);
            } catch (ArithmeticException exc) {
                if (!exc.getMessage().contains("overflow"))
                    throw exc;
                ++overflows;
            }
        }

        System.out.println("Overflows: " + overflows);
        System.out.println("Timings    Yun: " + yun.getMean());
        System.out.println("Timings Musser: " + musser.getMean());
    }

    @Test
    public void test3() throws Exception {
        MutablePolynomialZ poly = MutablePolynomialZ.create(0, 0, -1458, 6561, -6561);
        FactorDecomposition<MutablePolynomialZ> factorization = SquareFreeFactorizationYun(poly);
        assertFactorization(poly, factorization);
    }


    @Test
    public void test4() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(2, 1, 0, 1, 0, 2);
        poly = poly.multiply(MutablePolynomialMod.create(2, 1, 0, 1, 0, 2));
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test5() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(2, 1, 0, 0, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void testRandom6() throws Exception {
        final class test {
            long overflows = 0;
            final DescriptiveStatistics arithmetics = new DescriptiveStatistics();
            final DescriptiveStatistics timings = new DescriptiveStatistics();

            void run(int bound, int maxDegree, int maxNBase, int maxExponent, int nIterations) {
                RandomGenerator rnd = new Well1024a();
                RandomDataGenerator rndd = new RandomDataGenerator(rnd);
                long[] primes = {2, 3, 5, 7, 11, 13, 101};
                long start;
                for (int i = 0; i < nIterations; i++) {
                    if (i == 100)
                        timings.clear();
                    int nbase = rndd.nextInt(1, maxNBase);
                    for (long modulus : primes) {
                        MutablePolynomialMod poly;
                        start = System.nanoTime();
                        poly = MutablePolynomialMod.create(modulus, rndd.nextLong(1, 1000));
                        for (int j = 0; j < nbase; j++) {
                            MutablePolynomialMod f = randomPoly(rndd.nextInt(1, maxDegree), bound, rnd).modulus(modulus);
                            poly = poly.multiply(PolynomialArithmetics.polyPow(f, rndd.nextInt(1, maxExponent), true));
                        }
                        arithmetics.addValue(System.nanoTime() - start);
                        try {
                            start = System.nanoTime();
                            FactorDecomposition<MutablePolynomialMod> factorization = SquareFreeFactorization(poly);
                            timings.addValue(System.nanoTime() - start);
                            assertFactorization(poly, factorization);
                        } catch (ArithmeticException exc) {
                            if (!exc.getMessage().contains("overflow"))
                                throw exc;
                            ++overflows;
                        }

                    }
                }
            }
        }
        test smallPolys = new test();
        smallPolys.run(10, 3, 5, 5, 1000);
        System.out.println("Overflows: " + smallPolys.overflows);
        System.out.println("Timings:\n" + smallPolys.timings.getPercentile(50));
        System.out.println("Arithmetics:\n" + smallPolys.arithmetics.getPercentile(50));

        System.out.println("\n ==== \n");
        test largePolys = new test();
        largePolys.run(1000, 15, 10, 20, 10);
        System.out.println("Overflows: " + largePolys.overflows);
        System.out.println("Timings:\n" + largePolys.timings.getPercentile(50));
        System.out.println("Arithmetics:\n" + largePolys.arithmetics.getPercentile(50));
    }

    @Test
    public void test6a() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(5, 0, 0, 1, 3, 4, 3, 3, 2, 3);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6b() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(3, 0, 0, 0, 2);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6c() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(2, 0, 0, 0, 1, 1);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6d() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6e() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(5, 2, 3, 2, 1, 3, 3, 3);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6f() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(3, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test6g() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialMod.create(101, 0, 0, 0, 8, 20, 67, 55);
        assertFactorization(poly, SquareFreeFactorization(poly));
    }

    @Test
    public void test7() throws Exception {
        assertEquals(1, SquareFreeFactorization(MutablePolynomialZ.create(0, 0, 0, 0, 1)).factors.size());
        assertEquals(1, SquareFreeFactorization(MutablePolynomialMod.create(31, 0, 0, 0, 0, 1)).factors.size());
    }

    @Test
    public void test8() throws Exception {
        RandomSource rnd = new RandomSource(new Well1024a(), 10, 20, false);
        long modulus = 3;
        for (int i = 0; i < 1000; i++) {
            MutablePolynomialMod poly = rnd.take(modulus);
            assertTrue(isSquareFree(SquareFreePart(poly)));
        }
    }
}