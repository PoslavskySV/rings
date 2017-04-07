package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.univar.Factorization.factor;
import static cc.r2.core.poly.univar.FactorizationTestUtil.assertFactorization;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class FactorizationTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        assertTrue(Factorization.factor(lMutablePolynomialZ.create(3, 7).modulus(19)).get(0).isMonic());
    }

    @Test
    public void test2() throws Exception {
        BigInteger modulus = BigInteger.LONG_MAX_VALUE;
        modulus = modulus.multiply(modulus).increment().nextProbablePrime();
        bMutablePolynomialZp poly = bMutablePolynomialZ.create(
                BigInteger.valueOf(Long.MAX_VALUE),
                BigInteger.valueOf(Long.MAX_VALUE - 1),
                BigInteger.valueOf(Long.MAX_VALUE - 2)).modulus(modulus);
        for (int i = 0; i < 5; i++)
            poly = poly.square().add(poly.derivative()).increment();
        bFactorDecomposition<bMutablePolynomialZp> fct = factor(poly);
        Assert.assertEquals(7, fct.size());
        assertFactorization(poly, fct);
    }

    @Test
    public void test3() throws Exception {
        long modulus = 13;
        lMutablePolynomialZp poly = lMutablePolynomialZ.create(5, 8, 1, 5, 7, 0, 0, 1, 5, 7, 0, 9, 3, 2).modulus(modulus);
        assertFactorization(poly, Factorization.factor(poly));
    }

    @Test
    public void test4_randomZp() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nIterations = (int) its(1000, 3000);
        for (int n = 0; n < nIterations; n++) {
            int nFactors = rndd.nextInt(4, 8);
            long modulus = getModulusRandom(rndd.nextInt(5, 31));
            lMutablePolynomialZp poly = lMutablePolynomialZp.constant(modulus, rndd.nextLong(1, modulus - 1));
            int expectedNFactors = 0;
            for (int i = 0; i < nFactors; i++) {
                lMutablePolynomialZp m = RandomPolynomials.randomMonicPoly(rndd.nextInt(1, 5), modulus, rnd);
                if (!m.isConstant()) ++expectedNFactors;
                if (m.isZero()) continue;
                poly = poly.multiply(m);
            }

            lFactorDecomposition<lMutablePolynomialZp> lFactors = Factorization.factor(poly);
            assertTrue(lFactors.size() >= expectedNFactors);
            assertFactorization(poly, lFactors);

            if (n % 100 == 0) {
                bFactorDecomposition<bMutablePolynomialZp> bFactors = factor(poly.toBigPoly());
                bFactorDecomposition<bMutablePolynomialZp> converted = lFactorDecomposition.convert(lFactors);
                converted.canonicalForm();
                bFactors.canonicalForm();
                Assert.assertEquals(converted, bFactors);
            }
        }
    }

    @Test
    public void test5_randomZ() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nIterations = (int) its(100, 1000);
        int maxDegree = 30;
        for (int n = 0; n < nIterations; n++) {
            bMutablePolynomialZ poly = bMutablePolynomialZ.one();
            int expectedNFactors = 0;
            while (true) {
                bMutablePolynomialZ m = RandomPolynomials.randomPoly(rndd.nextInt(1, 15), BigInteger.LONG_MAX_VALUE, rnd);
                if (m.isZero()) continue;
                if (!m.isConstant()) ++expectedNFactors;
                poly = poly.multiply(m);
                if (poly.degree() >= maxDegree)
                    break;
            }

            FactorDecomposition<bMutablePolynomialZ> lFactors = Factorization.factor(poly);
            assertTrue(lFactors.size() >= expectedNFactors);
            assertFactorization(poly, lFactors);
        }
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test6_referenceZ() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
        lMutablePolynomialZ poly = a.multiply(b).primitivePart();

        lMutablePolynomialZp polyMod = poly.modulus(SmallPrimes.nextPrime(2131243213));


        DescriptiveStatistics timingZ = new DescriptiveStatistics(), timingZp = new DescriptiveStatistics();
        int nIterations = its(100, 10000);
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000) {timingZ.clear(); timingZp.clear();}

            long start;
            start = System.nanoTime();
            lFactorDecomposition<lMutablePolynomialZp> factorsZp = Factorization.factor(polyMod);
            long timeZp = System.nanoTime() - start;
            timingZp.addValue(timeZp);
            assertFactorization(polyMod, factorsZp);
            assertTrue(factorsZp.size() >= 2);
            assertEquals(5, factorsZp.size());

            start = System.nanoTime();
            lFactorDecomposition<lMutablePolynomialZ> factorsZ = Factorization.factor(poly);
            long timeZ = System.nanoTime() - start;
            timingZ.addValue(timeZ);
            assertFactorization(poly, factorsZ);
            assertTrue(factorsZ.size() >= 2);
            assertEquals(2, factorsZ.size());

            if (i > 9990) {
                System.out.println("Zp: " + TimeUnits.nanosecondsToString(timeZp));
                System.out.println("Z : " + TimeUnits.nanosecondsToString(timeZ));
            }
        }
        System.out.println(TimeUnits.statisticsNanotimeFull(timingZ));
        System.out.println(TimeUnits.statisticsNanotimeFull(timingZp));

//        -XX:+AggressiveOpts
//        DescriptiveStatistics:
//        n: 9000
//        min: 324us
//        max: 16ms
//        mean: 500us
//        std dev: 325us
//        median: 440us
//        skewness: 22ns
//        kurtosis: 855ns
//
//        DescriptiveStatistics:
//        n: 9000
//        min: 1373us
//        max: 39ms
//        mean: 1907us
//        std dev: 895us
//        median: 1719us
//        skewness: 20ns
//        kurtosis: 735ns
//
    }


    @Test
    @Benchmark(runAnyway = true)
    public void test7_referenceZb() throws Exception {
        bMutablePolynomialZ a = bMutablePolynomialZ.create(1, 2, 3, 5, 3, 2, 1),
                b = bMutablePolynomialZ.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = bMutablePolynomialZ.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = bMutablePolynomialZ.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = bMutablePolynomialZ.create(-11, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = bMutablePolynomialZ.create(-33, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        bMutablePolynomialZ poly = a.clone().multiply(b, c, d, e, f, g, g.clone().increment(), f.clone().increment());


        DescriptiveStatistics timing = new DescriptiveStatistics();
        int nIterations = its(100, 100);
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000)
                timing.clear();
            long start = System.nanoTime();
            bFactorDecomposition<bMutablePolynomialZ> factors = Factorization.factor(poly);
            long time = System.nanoTime() - start;
            timing.addValue(time);
            assertEquals(9, factors.size());
            assertFactorization(poly, factors);
        }
        System.out.println(TimeUnits.statisticsNanotimeFull(timing));

        //    -XX:+AggressivrOpts
        //    DescriptiveStatistics:
        //    n: 100
        //    min: 165ms
        //    max: 1327ms
        //    mean: 326ms
        //    std dev: 139ms
        //    median: 300ms
        //    skewness: 4ns
        //    kurtosis: 26ns
    }
}