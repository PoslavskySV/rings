package cc.r2.core.polynomial;

import cc.r2.core.polynomial.FactorizationTestUtil.RandomSource;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import static cc.r2.core.polynomial.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyPow;
import static cc.r2.core.polynomial.RandomPolynomials.randomPoly;
import static cc.r2.core.polynomial.SquareFreeFactorization.SquareFreeFactorization;
import static cc.r2.core.polynomial.SquareFreeFactorization.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 20/01/2017.
 */
public class SquareFreeFactorizationTest {
    @Test
    public void test1() throws Exception {
        MutablePolynomial poly = polyPow(MutablePolynomial.create(1, 3).multiply(2), 3, false).multiply(polyPow(MutablePolynomial.create(-3, -5, 7), 2, false));
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomial.create(1, 3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomial.create(3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomial.create(33);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutablePolynomial.create(22, 22).multiply(MutablePolynomial.create(12, 12, 12)).multiply(12);
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
            MutablePolynomial poly;
            try {
                poly = MutablePolynomial.create(rndd.nextLong(1, 10));
                for (int j = 0; j < nbase; j++) {
                    MutablePolynomial factor = randomPoly(rndd.nextInt(1, 3), 10, rnd);
                    int exponent = rndd.nextInt(1, 5);
                    poly = poly.multiply(polyPow(factor, exponent, true));
                }
            } catch (ArithmeticException e) {
                --i;
                continue;
            }
            try {
                long start = System.nanoTime();
                FactorDecomposition yunFactorization = SquareFreeFactorizationYun(poly);
                yun.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                FactorDecomposition musserFactorization = SquareFreeFactorizationMusser(poly);
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
        System.out.println("Timings Yun:\n" + yun);
        System.out.println("\nTimings Musser:\n" + musser);
    }

    @Test
    public void test3() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, -1458, 6561, -6561);
        FactorDecomposition factorization = SquareFreeFactorizationYun(poly);
        assertFactorization(poly, factorization);
    }


    @Test
    public void test4() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 0, 1, 0, 2);
        poly = poly.multiply(MutablePolynomial.create(1, 0, 1, 0, 2), 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
    }

    @Test
    public void test5() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 0, 0, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
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
                        MutablePolynomial poly;
                        start = System.nanoTime();
                        poly = MutablePolynomial.create(rndd.nextLong(1, 1000)).modulus(modulus);
                        for (int j = 0; j < nbase; j++) {
                            MutablePolynomial f = randomPoly(rndd.nextInt(1, maxDegree), bound, rnd);
                            f = f.modulus(modulus);
                            poly = poly.multiply(PolynomialArithmetics.polyPowMod(f, rndd.nextInt(1, maxExponent), modulus, true), modulus);
                        }
                        arithmetics.addValue(System.nanoTime() - start);
                        try {
                            start = System.nanoTime();
                            FactorDecomposition factorization = SquareFreeFactorization(poly, modulus);
                            timings.addValue(System.nanoTime() - start);
                            assertFactorization(poly, factorization, modulus);
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
        System.out.println("Timings:\n" + smallPolys.timings);
        System.out.println("Arithmetics:\n" + smallPolys.arithmetics);

        test largePolys = new test();
        largePolys.run(1000, 15, 10, 20, 10);
        System.out.println("Overflows: " + largePolys.overflows);
        System.out.println("Timings:\n" + largePolys.timings);
        System.out.println("Arithmetics:\n" + largePolys.arithmetics);
    }

    @Test
    public void test6a() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 1, 3, 4, 3, 3, 2, 3);
        assertFactorization(poly, SquareFreeFactorization(poly, 5), 5);
    }

    @Test
    public void test6b() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 0, 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test6c() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 0, 1, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
    }

    @Test
    public void test6d() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test6e() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(2, 3, 2, 1, 3, 3, 3);
        assertFactorization(poly, SquareFreeFactorization(poly, 5), 5);
    }

    @Test
    public void test6f() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test6g() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 0, 8, 20, 67, 55);
        assertFactorization(poly, SquareFreeFactorization(poly, 101), 101);
    }

    @Test
    public void test7() throws Exception {
        assertEquals(1, SquareFreeFactorization(MutablePolynomial.create(0, 0, 0, 0, 1)).factors.size());
        assertEquals(1, SquareFreeFactorization(MutablePolynomial.create(0, 0, 0, 0, 1), 31).factors.size());
    }

    @Test
    public void test8() throws Exception {
        RandomSource rnd = new RandomSource(new Well1024a(), 10, 20, false);
        long modulus = 3;
        for (int i = 0; i < 1000; i++) {
            MutablePolynomial poly = rnd.take(modulus);
            assertTrue(isSquareFree(SquareFreePart(poly, modulus), modulus));
        }
    }
}