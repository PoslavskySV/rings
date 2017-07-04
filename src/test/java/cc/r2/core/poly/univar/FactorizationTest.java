package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.FiniteField;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.FactorDecompositionTest.assertFactorization;
import static cc.r2.core.poly.univar.Factorization.factor;
import static cc.r2.core.poly.univar.Factorization.factorInFiniteField;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class FactorizationTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        assertTrue(Factorization.factorInFiniteField(lUnivariatePolynomialZ.create(3, 7).modulus(19)).get(0).isMonic());
    }

    @Test
    public void test2() throws Exception {
        BigInteger modulus = BigInteger.LONG_MAX_VALUE;
        modulus = modulus.multiply(modulus).increment().nextProbablePrime();
        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(new IntegersModulo(modulus),
                BigInteger.valueOf(Long.MAX_VALUE),
                BigInteger.valueOf(Long.MAX_VALUE - 1),
                BigInteger.valueOf(Long.MAX_VALUE - 2));
        for (int i = 0; i < 5; i++)
            poly = poly.square().add(poly.derivative()).increment();
        FactorDecomposition<UnivariatePolynomial<BigInteger>> fct = factorInFiniteField(poly);
        Assert.assertEquals(7, fct.size());
        assertFactorization(poly, fct);
    }

    @Test
    public void test3() throws Exception {
        long modulus = 13;
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(5, 8, 1, 5, 7, 0, 0, 1, 5, 7, 0, 9, 3, 2).modulus(modulus);
        assertFactorization(poly, Factorization.factorInFiniteField(poly));
    }

    @Test
    public void test4_randomZp() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nIterations = its(1000, 3000);
        for (int n = 0; n < nIterations; n++) {
            int nFactors = rndd.nextInt(4, 8);
            long modulus = getModulusRandom(rndd.nextInt(2, 31));
            lUnivariatePolynomialZp poly = lUnivariatePolynomialZp.constant(modulus, rndd.nextLong(1, modulus - 1));
            int expectedNFactors = 0;
            for (int i = 0; i < nFactors; i++) {
                lUnivariatePolynomialZp m = RandomPolynomials.randomMonicPoly(rndd.nextInt(1, 5), modulus, rnd);
                if (m.isZero()) continue;
                if (m.isMonomial()) continue;
                if (!m.isConstant()) ++expectedNFactors;
                poly = poly.multiply(m);
            }

            try {
                FactorDecomposition<lUnivariatePolynomialZp> lFactors = Factorization.factorInFiniteField(poly);
                assertTrue(lFactors.sumExponents() >= expectedNFactors);
                assertFactorization(poly, lFactors);

                if (n % 100 == 0) {
                    FactorDecomposition<UnivariatePolynomial<BigInteger>> bFactors = factorInFiniteField(poly.toBigPoly());
                    FactorDecomposition<UnivariatePolynomial<BigInteger>> converted = Factorization.convertFactorizationToBigIntegers(lFactors);
                    converted.canonicalForm();
                    bFactors.canonicalForm();
                    Assert.assertEquals(converted, bFactors);
                }
            } catch (Throwable e) {
                System.out.println(expectedNFactors);
                System.out.println(modulus);
                System.out.println(poly.toStringForCopy());
                throw e;
            }
        }
    }

    @Test
    public void test4_randomZp_a() throws Exception {
        long modulus = 59;
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(46, 16, 1, 54, 16, 57, 22, 15, 31, 21).modulus(modulus);
        FactorDecomposition<lUnivariatePolynomialZp> fct = factorInFiniteField(poly);
        assertEquals(5, fct.size());
        assertEquals(6, fct.sumExponents());
    }

    @Test
    public void test5_randomZ() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nIterations = (int) its(100, 1000);
        int maxDegree = 30;
        for (int n = 0; n < nIterations; n++) {
            UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(1);
            int expectedNFactors = 0;
            while (true) {
                UnivariatePolynomial<BigInteger> m = RandomPolynomials.randomPoly(rndd.nextInt(1, 15), BigInteger.LONG_MAX_VALUE, rnd);
                if (m.isZero()) continue;
                if (!m.isConstant()) ++expectedNFactors;
                poly = poly.multiply(m);
                if (poly.degree() >= maxDegree)
                    break;
            }

            FactorDecomposition<UnivariatePolynomial<BigInteger>> lFactors = Factorization.factorInZ(poly);
            assertTrue(lFactors.size() >= expectedNFactors);
            assertFactorization(poly, lFactors);
        }
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test6_referenceZ() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
        lUnivariatePolynomialZ poly = a.multiply(b).primitivePart();

        lUnivariatePolynomialZp polyMod = poly.modulus(SmallPrimes.nextPrime(2131243213));


        DescriptiveStatistics timingZ = new DescriptiveStatistics(), timingZp = new DescriptiveStatistics();
        int nIterations = its(100, 10000);
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000) {timingZ.clear(); timingZp.clear();}

            long start;
            start = System.nanoTime();
            FactorDecomposition<lUnivariatePolynomialZp> factorsZp = Factorization.factorInFiniteField(polyMod);
            long timeZp = System.nanoTime() - start;
            timingZp.addValue(timeZp);
            assertFactorization(polyMod, factorsZp);
            assertTrue(factorsZp.size() >= 2);
            assertEquals(5, factorsZp.size());

            start = System.nanoTime();
            FactorDecomposition<lUnivariatePolynomialZ> factorsZ = Factorization.factorInZ(poly);
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
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(1, 2, 3, 5, 3, 2, 1),
                b = UnivariatePolynomial.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = UnivariatePolynomial.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = UnivariatePolynomial.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = UnivariatePolynomial.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = UnivariatePolynomial.create(-11, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = UnivariatePolynomial.create(-33, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        UnivariatePolynomial<BigInteger> poly = a.clone().multiply(b, c, d, e, f, g, g.clone().increment(), f.clone().increment());


        DescriptiveStatistics timing = new DescriptiveStatistics();
        int nIterations = its(100, 100);
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000)
                timing.clear();
            long start = System.nanoTime();
            FactorDecomposition<UnivariatePolynomial<BigInteger>> factors = Factorization.factorInZ(poly);
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

    @Test
    public void test8() throws Exception {
        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(0, -13284, -15390, -33552, -55998, 2151, 381296, 1573628, 1135112, 688800, 358176, 1119300);
        assertEquals(2, factor(poly).size());

        poly = UnivariatePolynomial.create(0, -13284, -15390, -33552, -55998, 2151, 381296, 1573628, 1135112, -688800, 358176, 1119300);
        assertEquals(4, factor(poly).size());
    }

    @Test
    public void testFiniteField1() throws Exception {
        lUnivariatePolynomialZp irreducible = IrreduciblePolynomials.randomIrreduciblePolynomial(2, 10, getRandom());
        FiniteField<lUnivariatePolynomialZp> domain = new FiniteField<>(irreducible);
        lUnivariatePolynomialZp
                c0 = domain.valueOf(lUnivariatePolynomialZ.create(1, 2, 3, 4, 5).modulus(irreducible.domain)),
                c1 = domain.valueOf(lUnivariatePolynomialZ.create(1, -2, 3, -4, -5).modulus(irreducible.domain)),
                c2 = domain.valueOf(lUnivariatePolynomialZ.create(11, 12, 13, 14, 15).modulus(irreducible.domain)),
                c3 = domain.add(c0, c1),
                c4 = domain.subtract(c1, c2),
                c5 = domain.multiply(c0, c1);

        UnivariatePolynomial<lUnivariatePolynomialZp>
                poly1 = UnivariatePolynomial.create(domain, c0, c1, c2, c3, c4, c5),
                poly2 = UnivariatePolynomial.create(domain, c5, c4, c3, c2, c1, c0),
                poly = poly1.clone().multiply(poly2).multiply(poly1.clone().add(poly2));

        FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>> factors = Factorization.factorInFiniteField(poly);
        assertFactorization(poly, factors);
    }

    @Test
    public void testFiniteFieldRandom1() throws Exception {
        testFiniteFieldRandom(2, 30, 4, 8, 4, 8, 100);
    }

    @Test
    public void testFiniteFieldRandom2_GF2p() throws Exception {
        testFiniteFieldRandom(2, 2, 4, 8, 4, 8, its(100, 1000));
    }

    private static void testFiniteFieldRandom(int minModulusBits, int maxModulusBits,
                                              int minDegree, int maxDegree,
                                              int minFactors, int maxFactors,
                                              int nIterations) {
        System.out.println("> Factorization in Galois fields");
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics statistics = new DescriptiveStatistics();
        int prevPercent = -1, currPercent;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                statistics.clear();

            if ((currPercent = (int) (100. * n / nIterations)) != prevPercent) {
                prevPercent = currPercent;
                System.out.print(">");
                System.out.flush();
            }

            long modulus = getModulusRandom(rndd.nextInt(minModulusBits, maxModulusBits));
            lUnivariatePolynomialZp irreducible = IrreduciblePolynomials.randomIrreduciblePolynomial(modulus, rndd.nextInt(minDegree, maxDegree), rnd);
            FiniteField<lUnivariatePolynomialZp> field = new FiniteField<>(irreducible);
            int nFactors = rndd.nextInt(minFactors, maxFactors);

            UnivariatePolynomial<lUnivariatePolynomialZp> poly = UnivariatePolynomial.one(field);
            int expectedNFactors = 0;
            for (int i = 0; i < nFactors; i++) {
                UnivariatePolynomial<lUnivariatePolynomialZp> m = RandomPolynomials.randomPoly(rndd.nextInt(1, 5), field, rnd);
                if (m.isZero()) continue;
                if (m.isMonomial()) continue;
                if (!m.isConstant()) ++expectedNFactors;
                poly = poly.multiply(m);
            }

            FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>> lFactors = null;
            try {
                long start = System.nanoTime();
                lFactors = Factorization.factorInFiniteField(poly);
                statistics.addValue(System.nanoTime() - start);
                assertTrue(lFactors.sumExponents() >= expectedNFactors);
                assertFactorization(poly, lFactors);
            } catch (Throwable e) {
                System.out.println(expectedNFactors);
                System.out.println(modulus);
                System.out.println(irreducible);
                System.out.println(poly.toStringForCopy());
                System.out.println(lFactors);
                throw e;
            }
        }

        System.out.println();
        System.out.println("============ Timing ============ ");
        System.out.println(TimeUnits.statisticsNanotimeFull(statistics));
        System.out.println();
    }

    @Test
    public void testFiniteField2() throws Exception {
        lUnivariatePolynomialZp irreducible = lUnivariatePolynomialZ.create(1, 1, 0, 1).modulus(2);
        FiniteField<lUnivariatePolynomialZp> domain = new FiniteField<>(irreducible);
        UnivariatePolynomial<lUnivariatePolynomialZp> poly =
                UnivariatePolynomial.create(domain,
                        lUnivariatePolynomialZ.create(0, 0, 1).modulus(2),
                        lUnivariatePolynomialZ.create(0, 1, 1).modulus(2),
                        lUnivariatePolynomialZ.create(0, 0, 1).modulus(2),
                        lUnivariatePolynomialZ.create(1).modulus(2),
                        lUnivariatePolynomialZ.create(0, 0, 1).modulus(2),
                        lUnivariatePolynomialZ.create(0).modulus(2),
                        lUnivariatePolynomialZ.create(1, 0, 1).modulus(2),
                        lUnivariatePolynomialZ.create(0, 0, 1).modulus(2));

        FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>> factors = Factorization.factorInFiniteField(poly);
        assertFactorization(poly, factors);
    }
}