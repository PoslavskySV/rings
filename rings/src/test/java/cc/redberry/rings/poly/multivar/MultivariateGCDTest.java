package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.io.IStringifier;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.RandomDataGenerator;
import cc.redberry.rings.util.RandomUtil;
import cc.redberry.rings.util.TimeUnits;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.BitSet;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.zip.GZIPInputStream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.PolynomialMethods.polyPow;
import static cc.redberry.rings.poly.multivar.AMultivariatePolynomial.renameVariables;
import static cc.redberry.rings.poly.multivar.MultivariateDivision.divideExact;
import static cc.redberry.rings.poly.multivar.MultivariateDivision.dividesQ;
import static cc.redberry.rings.poly.multivar.MultivariateGCD.*;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;
import static cc.redberry.rings.poly.multivar.RandomMultivariatePolynomials.randomPolynomial;
import static cc.redberry.rings.util.TimeUnits.nanosecondsToString;
import static org.junit.Assert.*;

/**
 * @since 1.0
 */
public class MultivariateGCDTest extends AMultivariateTest {
    private static void assertBrownGCD(MultivariatePolynomial<BigInteger> gcd,
                                       MultivariatePolynomial<BigInteger> a,
                                       MultivariatePolynomial<BigInteger> b) {
        MultivariatePolynomial<BigInteger> actualGCD = BrownGCD(a, b);
        MultivariatePolynomialZp64 lActualGCD = BrownGCD(asOverZp64(a), asOverZp64(b));
        Assert.assertTrue(dividesQ(actualGCD, gcd));
        Assert.assertEquals(asOverZp64(actualGCD).monic(), lActualGCD.monic());
    }

    @Test // =====> testBrown1   elapsed 578us
    public void testBrown1() throws Exception {
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(1321323));
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c", domain, MonomialOrder.DEFAULT),
                b = parse("a^2+2*b^2 + 2*c", domain, MonomialOrder.DEFAULT),
                gcd = parse("c*a+b+a+ c*a^3", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown2   elapsed 444us
    public void testBrown2() throws Exception {
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(1321323));
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c", domain, MonomialOrder.DEFAULT),
                b = parse("a^2+2*b^2 + 2*c", domain, MonomialOrder.DEFAULT),
                gcd = parse("c*a*b+b*b+a*b+ c*a^3*b", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown3   elapsed 836us
    public void testBrown3() throws Exception {
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(659));
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, MonomialOrder.DEFAULT),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, MonomialOrder.DEFAULT),
                gcd = parse("7*b+655*a*b^3*c^2+6*a^2*b^3*c^2", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown3a   elapsed 1895us
    public void testBrown3a() throws Exception {
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, MonomialOrder.DEFAULT),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, MonomialOrder.DEFAULT),
                gcd = parse("7*b^6+655*a*b^3*c^6+6*a^2*b^3*c^4", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown4   elapsed 3ms
    public void testBrown4() throws Exception {
        IntegersZp domain = new IntegersZp(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*b^5*c+a*b^3+a^2*b^2+a^2*b^2*c+a^3*b*c^3", domain, MonomialOrder.DEFAULT, vars),
                b = parse("9*a*b^2*c^6+a*b^4*c^6+a^2*b^2*c^3+a^5*b^2+a^5*b^6*c^4+a^6*c+a^6*b^2*c", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("653*b^3*c^4+b^4+b^5*c^3+a^2*b*c^2+a^4*b^2*c^4", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown5   elapsed 9ms
    public void testBrown5() throws Exception {
        RandomGenerator rnd = PrivateRandom.getRandom();
        rnd.setSeed(28);
        IntegersZp domain = new IntegersZp(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("561*a^2*c^2+a^2*b^2*c+a^3*b+a^4*b^2+a^4*b^5*c^3+a^5*b", domain, MonomialOrder.DEFAULT, vars),
                b = parse("561*a*c^3+a*b^4*c^5+a^2*c^2+a^2*b^6*c^3+a^3*b^6*c^5+a^5*b^5*c^3+a^5*b^5*c^6", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("4*c^2+b^4+a^2*b^4*c+a^3*b^2*c+a^3*b^6+a^5*b^2*c^6+a^6*b^5", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown6   elapsed 5ms
    public void testBrown6() throws Exception {
        PrivateRandom.getRandom().setSeed(1564);
        IntegersZp domain = new IntegersZp(937);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("931*a^3*b^4*c+a^4+a^4*b^6*c^2+a^5*b*c^3+a^6*b*c^2", domain, MonomialOrder.DEFAULT, vars),
                b = parse("932*b*c+a*b^6*c^2+a^3*b*c^2+a^3*b^3*c^5+a^3*b^5*c+a^5*b^5*c^3+a^6*b^2*c^6+a^6*b^4*c^5+a^6*b^6", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("935*c^2+c^4+a^3*b*c^5+a^3*b^2*c^3+a^4*b^3", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown7   elapsed 3ms
    public void testBrown7() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersZp domain = new IntegersZp(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("563*b^2*c+6*b^4*c^3+4*a*b^4*c^3+563*a*b^4*c^4+560*a^2*b^5*c^2+9*a^3*b^4*c+5*a^4*c^2+7*a^4*b^3*c^5+4*a^5*b^4*c^5+6*a^5*b^5", domain, MonomialOrder.DEFAULT, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("4+8*b^2*c^3+4*b^3+8*b^3*c+7*a*c+a*b*c+7*a^2*b^2+2*a^2*b^2*c^2+5*a^3*c^2+5*a^3*c^3", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown_random1   elapsed 5s
    public void testBrown_random1() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 3, nVarsMax = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test // =====> testBrown_random2   elapsed 257ms
    public void testBrown_random2() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int nVarsMin = 3, nVarsMax = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);

        GCDSampleData sampleData = coprimeData(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd));
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test // =====> testBrown_sparse_variables_random3   elapsed 6s
    public void testBrown_sparse_variables_random3() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<MonomialZp64, MultivariatePolynomialZp64>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test // =====> testBrown8   elapsed 635us
    public void testBrown8() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersZp domain = new IntegersZp(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, MonomialOrder.DEFAULT, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("a^2 + c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown9   elapsed 125us
    public void testBrown9() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersZp domain = new IntegersZp(569);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, MonomialOrder.DEFAULT, vars),
                b = parse("b^2 + d^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("d^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test // =====> testBrown10   elapsed 552us
    public void testBrown10() throws Exception {
        IntegersZp domain = new IntegersZp(5642359);
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger>
                a = parse("1199884 + 4783764*b + b^2 + 3215597*a*b + 2196297*a*b^2 + 4781733*a^4 + a^4*b + 2196297*a^5*b", domain, MonomialOrder.DEFAULT, vars),
                b = parse("4645946 + 3921107*b + b^2 + 5605437*a*b + 2196297*a*b^2 + 4781733*a^3 + a^3*b + 2196297*a^4*b", domain, MonomialOrder.DEFAULT, vars);

        MultivariatePolynomial<BigInteger> gcdActual = BrownGCD(a, b);
        gcdActual = gcdActual.monic().multiply(domain.valueOf(4781733));
        MultivariatePolynomial<BigInteger> expected = parse("1574588 + 4559668*b + 4781733*a*b", domain, MonomialOrder.DEFAULT, vars);
        assertEquals(expected, gcdActual);
    }


    private static void assertZippelGCD(MultivariatePolynomial<BigInteger> gcd,
                                        MultivariatePolynomial<BigInteger> a,
                                        MultivariatePolynomial<BigInteger> b) {
        MultivariatePolynomial<BigInteger> actualGCD = ZippelGCD(a, b);
        MultivariatePolynomialZp64 lActualGCD = ZippelGCD(asOverZp64(a), asOverZp64(b));
        Assert.assertTrue(dividesQ(actualGCD, gcd));
        Assert.assertEquals(asOverZp64(actualGCD).monic(), lActualGCD.monic());
    }

    @Test // =====> testZippel1   elapsed 11ms
    public void testZippel1() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^2 + b^2 + a*c^2", domain, MonomialOrder.DEFAULT, vars),
                b = parse("a^2 + 2*b^2 + b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("a^2 + c*a + b + a*c*b + a*c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        RandomGenerator rnd = getRandom();

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(101);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int i = 0; i < 100; i++) {
            SparseInterpolation<BigInteger> sparseInterpolation
                    = createInterpolation(variable, a, b, skeleton, 1, rnd);
            BigInteger point = domain.randomElement(rnd);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }

        MultivariatePolynomialZp64 la = asOverZp64(a), lb = asOverZp64(b),
                lskeleton = asOverZp64(skeleton), lgcd = asOverZp64(gcd);
        for (int i = 0; i < 100; i++) {
            lSparseInterpolation sparseInterpolation
                    = createInterpolation(variable, la, lb, lskeleton, 1, rnd);
            long point = domain.asMachineRing().randomElement(rnd);
            assertEquals(lgcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test // =====> testZippel2   elapsed 1647us
    public void testZippel2() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^2 + b^2 + a*c^2", domain, MonomialOrder.DEFAULT, vars),
                b = parse("a^2 + 2*b^2 + b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("a^2 + c*a + b + a*c*b + a*c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertZippelGCD(gcd, a, b);
    }

    @Test // =====> testZippel5   elapsed 1139us
    public void testZippel5() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(31579447));
        System.out.println(domain);
        MultivariatePolynomial<BigInteger>
                a = parse("b^2 + a*c + a*b*c^2 + a*b^2 + a^5", domain, MonomialOrder.DEFAULT, vars),
                b = parse("1 + a*b", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("1 + a*b*c^2 + a^2*c^2 + a^2*b*c^2 + a^2*b^2*c^2 + a^7", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int expEvals : Arrays.asList(1, Integer.MAX_VALUE)) {
            SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, expEvals, rnd);
            BigInteger point = domain.valueOf(1324);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test // =====> testZippel6   elapsed 920us
    public void testZippel6() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(31579447));
        MultivariatePolynomial<BigInteger>
                a = parse("1+3*b*c^2+7*b^2*c^2+4*a+7*a^3*c^2+a^6", domain, MonomialOrder.DEFAULT, vars),
                b = parse("b^3*c^2+a*b^2*c+9*a^2*b*c+9*a^3*b*c^2+a^7", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("b^3*c^2+2*a*c+5*a*b+2*a*b^2*c+4*a*b^3*c^2+5*a^3*b+a^7", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int expEvals : Arrays.asList(1, Integer.MAX_VALUE)) {
            SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, expEvals, rnd);
            BigInteger point = domain.valueOf(1324);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test // =====> testZippel_monic_random1   elapsed 948ms
    public void testZippel_monic_random1() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 10, minDegree = 3, maxDegree = 5, minSize = 5, maxSize = 10;

        int nIterations = its(100, 1500);
        GCDSampleData<MonomialZp64, MultivariatePolynomialZp64> sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n % 100 == 0) System.out.println(n);

            GCDSample<MonomialZp64, MultivariatePolynomialZp64> data = sampleData.nextSample(true, true);

            int variable = data.a.nVariables - 1;
            long seed;
            do {seed = data.a.ring.randomElement(rnd);} while (seed == 0);
            MultivariatePolynomialZp64 skeleton = data.gcd.evaluate(variable, seed);

            for (int i = 0; i < 10; i++) {
                int rndSeed = i ^ n;
                rnd.setSeed(rndSeed);
                SparseInterpolation<BigInteger> sparseInterpolation
                        = createInterpolation(variable, data.a.toBigPoly(), data.b.toBigPoly(), skeleton.toBigPoly(), 1 + (1 << 30) * (i % 2), rnd);

                lSparseInterpolation lSparseInterpolation
                        = createInterpolation(variable, data.a, data.b, skeleton, 1, rnd);
                long point = data.a.ring.randomElement(rnd);
                try {
                    MultivariatePolynomialZp64 expected = data.gcd.evaluate(variable, point).monic();

                    MultivariatePolynomialZp64 lActual = lSparseInterpolation.evaluate(point).monic();
                    assertEquals(expected, lActual);

                    MultivariatePolynomial<BigInteger> bActual = sparseInterpolation.evaluate(BigInteger.valueOf(point)).monic();
                    assertEquals(expected.toBigPoly(), bActual);
                } catch (Throwable e) {
                    System.out.println("rnd seed : " + rndSeed);
                    System.out.println("seed point : " + seed);
                    System.out.println("point : " + point);
                    System.out.println("a: " + data.a);
                    System.out.println("b: " + data.b);
                    System.out.println("gcd : " + data.gcd);
                    throw e;
                }
            }
        }
    }

    @Test // =====> testZippel7   elapsed 1207us
    public void testZippel7() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c", "d"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(31579447));
        MultivariatePolynomial<BigInteger>
                a = parse("29322275*b+5*b*c+6*a^2*b^3*c^2+29322274*a^2*b^3*c^2*d^3+5*a^3*b*c^2*d^2+a^11", domain, MonomialOrder.DEFAULT, vars),
                b = parse("7*a^3*b^3*c^3*d^4+9*a^3*b^4*c+29322274*a^3*b^4*c^5*d^2+29322277*a^5*b*c*d+a^5*b^4+a^15", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("4*d^2+8*c^2*d+4*b*c+6*b^3*c^2*d^2+2*a^3*b^2*c^3*d+a^10", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int expEvals : Arrays.asList(1, Integer.MAX_VALUE)) {
            SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, expEvals, rnd);
            BigInteger point = domain.valueOf(1324);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test // =====> testZippel_monic_random2   elapsed 3s
    public void testZippel_monic_random2() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 3, maxDegree = 5, minSize = 5, maxSize = 10;
        int nIterations = its(100, 1500);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        sampleData.monic = true;
        sampleData.primitive = true;

        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Zippel (monic)", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Zippel (monic)", MultivariateGCD::ZippelGCD));

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Zippel (monic)", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test(timeout = 10000) // =====> testZippel9   elapsed 5ms
    public void testZippel9() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(26478253);
        PrivateRandom.getRandom().setSeed(0);
        MultivariatePolynomial<BigInteger>
                a = parse("26478246*a*c^2+7*a*b+26478250*a*b*c^2+26478249*a*b^3*c^2+26478248*a^2*c^2+8*a^3*b*c^2+a^7", domain, MonomialOrder.DEFAULT, vars),
                b = parse("4*b^3*c^2+7*a+5*a*b+8*a*b^2+6*a^3*b^2*c+a^7", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("26478248*a*b^2*c^2+3*a*b^3*c^2+2*a^2*b^3*c^2+5*a^3*c+a^8", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertNotNull(ZippelGCD(a, b));
        assertNotNull(ZippelGCD(asOverZp64(a), asOverZp64(b)));
    }

    @Test // =====> testZippel8   elapsed 927us
    public void testZippel8() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(31579447);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a+29923129*a*b*c^2+3*a*b^2+29923132*a^2*b*c^2+7*a^3*c", domain, MonomialOrder.DEFAULT, vars),
                b = parse("4*c^2+29923126*a*c^2+5*a*b+6*a^2*b^2*c^3+29923128*a^3*c^3", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("29923132+8*b*c^3+29923132*b^2*c^3+8*a*b*c+7*a^3*b^3*c", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int expEvals : Arrays.asList(1, Integer.MAX_VALUE)) {
            SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, expEvals, rnd);
            BigInteger point = domain.valueOf(1324);
            assertEquals(gcd.evaluate(variable, point).monic(), sparseInterpolation.evaluate(point).monic());
        }
    }

    @Test// =====> testZippel_nonmonic_random2   elapsed 6s
    public void testZippel_nonmonic_random2() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 3, maxDegree = 5, minSize = 5, maxSize = 10;

        int nIterations = its(100, 1500);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        sampleData.monic = false;
        sampleData.primitive = false;

        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Zippel (non monic)", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Zippel (non monic)", MultivariateGCD::ZippelGCD));

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Zippel (non monic)", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test // =====> testZippel_nonmonic_random3   elapsed 26s
    public void testZippel_nonmonic_random3() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        RandomGenerator rnd = getRandom();
        int nVarsMin = 7, nVarsMax = 10, minDegree = 7, maxDegree = 10, minSize = 7, maxSize = 10;
        int nIterations = its(500, 1500);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        sampleData.monic = false;
        sampleData.primitive = false;

        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Zippel (non monic)", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Zippel (non monic)", MultivariateGCD::ZippelGCD));
    }

    @Test // =====> testZippel3   elapsed 4ms
    public void testZippel3() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("5*a^2*c^2+5*a^2*b^2*c^2+5*a^2*b^4*c^3+9*a^2*b^5*c^5+25709547*a^3*b^6*c^6+8*a^4*b*c^3+a^4*b^3*c+5*a^4*b^3*c^6", domain, MonomialOrder.DEFAULT, vars),
                b = parse("3*a*b^2*c^2+2*a^2*b^4+25709540*a^4*b*c^6+7*a^5*c^2+8*a^6*b*c^3", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("a + 5*b^2*c^6+2*a^4*b^4*c^5+25709543*a^5*b^2*c^5+9*a^6*c+25709540*a^6*c^3", domain, MonomialOrder.DEFAULT, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        assertZippelGCD(gcd, a, b);
    }

    @Test // =====> testZippel10   elapsed 995us
    public void testZippel10() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a + b + c", domain, MonomialOrder.DEFAULT, vars),
                b = parse("a - b + c", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("a^3*b+2*a^3+c*a^3+12*b^2+24*b+12*b*c", domain, MonomialOrder.DEFAULT, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        assertZippelGCD(gcd, a, b);
    }

    @Ignore
    @Test
    public void testZippel4_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b", "c", "d", "e"};
        IntegersZp domain = new IntegersZp(SmallPrimes.nextPrime(100000));
        MultivariatePolynomial<BigInteger>
                a = parse("2147483167*a^4*b^60*c^57*d^26*e+44*a^8*b^39*c^67*d^22*e^17+38*a^32*b^6*c^13*d^10*e^3+357*a^36*b^34*c^60*d^2*e^59+563*a^42*b^41*c^45*d^52*e^14+257*a^44*b^68*c^43*d^2*e^73+613*a^48*b^55*c^22*d^32*e^19+2147483093*a^52*b^26*c^4*d^72*e^32+19*a^52*b^40*c^26*d^45*e^55+639*a^55*b^72*c^55*d^65", domain, MonomialOrder.DEFAULT, vars),
                b = parse("2147483150*b^25*c^18*d^62*e^59+2147482723*a^4*b^5*c^65*d^26*e^7+261*a^15*b^60*c^59*d^63*e^53+394*a^27*b^22*c^34*d^54*e^13+952*a^39*b^48*c^17*d^54*e^16+243*a^60*b^15*c^3*d^51*e^46+40*a^61*b^56*c^39*d^40*e^21+555*a^62*b^20*c^20*d^60*e^47+627*a^67*b^8*c^22*d^67*e^61+447*a^70*b^59*c^71*d^24*e^5", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("35*a*b^36*c^74*d^62*e^51+376*a^2*b^28*c^64*e^53+893*a^6*b^13*c^60*d^44*e^42+23*a^8*b^71*c^40*d^36*e^11+783*a^20*b^28*c^12*d^31*e^68+2147482938*a^31*b^30*c^40*d^65*e^72+2147482960*a^31*b^49*c^38*d^71*e^55+737*a^47*b^15*c^71*d^13*e^72+868*a^53*b^30*c^40*d^29*e^46+898*a^61*b^71*c^13*d^50*e^66", domain, MonomialOrder.DEFAULT, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        System.out.println(a);
        System.out.println(b);

        MultivariatePolynomialZp64
                aL = asOverZp64(a),
                bL = asOverZp64(b);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertEquals(10, ZippelGCD(aL, bL).size());
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            start = System.nanoTime();
            System.out.println(ZippelGCD(aL.clone().increment(), bL));
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            System.out.println();
//            System.out.println(TimeUnits.nanosecondsToString(MultivariateGCD.BROWN));
        }
    }

    @Ignore
    @Test
    public void testZippel5_bivariate_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b"};
        IntegersZp domain = new IntegersZp(100011111111101L);
//        IntegersZp ring = new IntegersZp(100019);
        MultivariatePolynomial<BigInteger>
                a = parse("38*a^32*b^6 + 2147483093*a^52*b^26 + 357*a^36*b^34 + 44*a^8*b^39 + 19*a^52*b^40 + 563*a^42*b^41 + 613*a^48*b^55 + 2147483167*a^4*b^60 + 257*a^44*b^68 + 639*a^55*b^72", domain, MonomialOrder.DEFAULT, vars),
                b = parse("2147482723*a^4*b^5 + 627*a^67*b^8 + 243*a^60*b^15 + 555*a^62*b^20 + 394*a^27*b^22 + 2147483150*b^25 + 952*a^39*b^48 + 40*a^61*b^56 + 447*a^70*b^59 + 261*a^15*b^60", domain, MonomialOrder.DEFAULT, vars),
                gcd = parse("893*a^6*b^13 + 737*a^47*b^15 + 376*a^2*b^28 + 783*a^20*b^28 + 2147482938*a^31*b^30 + 868*a^53*b^30 + 35*a*b^36 + 2147482960*a^31*b^49 + 23*a^8*b^71 + 898*a^61*b^71", domain, MonomialOrder.DEFAULT, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        System.out.println(a);
        System.out.println(b);

        MultivariatePolynomialZp64
                aL = asOverZp64(a),
                bL = asOverZp64(b);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertEquals(10, ZippelGCD(aL, bL).size());
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            start = System.nanoTime();
            System.out.println(ZippelGCD(aL.clone().increment(), bL));
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            System.out.println();
        }
    }

    @Test // =====> testZippel_sparse_variables_random   elapsed 6s
    public void testZippel_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<MonomialZp64, MultivariatePolynomialZp64>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD));
    }

    @Test // =====> testPairedIterator1   elapsed 27ms
    public void testPairedIterator1() throws Exception {
        RandomGenerator rnd = getRandom();
        int nIterations = its(1000, 1000);
        for (int n = 0; n < nIterations; n++) {
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(5, 50, 20, rnd),
                    b = randomPolynomial(5, 50, 20, rnd);

            PairedIterator<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>, Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> it = new PairedIterator<>(a, b);
            new PairedIterator<>(a, b);

            MultivariatePolynomial<BigInteger> acc = a.createZero();
            while (it.hasNext()) {
                it.advance();
                acc.add(it.aTerm);
                acc.add(it.bTerm);
                assertTrue(it.aTerm.coefficient.isZero() || it.bTerm.coefficient.isZero() || 0 == a.ordering.compare(it.aTerm, it.bTerm));
            }
            assertEquals(a.clone().add(b), acc);
        }
    }

    @Test // =====> testPairedIterator2   elapsed 271us
    public void testPairedIterator2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse(" y  +   0  +  x*y   + x*y^2 + x*y^3 + x*y^4  +   0    "),
                b = parse(" 0  +   x  +   0    +   0   +  0    +   0    + x^2*y^2");

        PairedIterator<
                Monomial<BigInteger>, MultivariatePolynomial<BigInteger>,
                Monomial<BigInteger>, MultivariatePolynomial<BigInteger>
                > it = new PairedIterator<>(a, b);

        MultivariatePolynomial<BigInteger> acc = a.createZero();
        while (it.hasNext()) {
            it.advance();
            acc.add(it.aTerm);
            acc.add(it.bTerm);
            assertTrue(it.aTerm.coefficient.isZero() || it.bTerm.coefficient.isZero() || 0 == a.ordering.compare(it.aTerm, it.bTerm));
        }
        assertEquals(a.clone().add(b), acc);
    }

    @Test // =====> testSparseInterpolation1   elapsed 6ms
    public void testSparseInterpolation1() throws Exception {
        IntegersZp domain = new IntegersZp(31574773);
        MultivariatePolynomial<BigInteger>
                a = parse("31574768*a*b^2*c^4+4*a^4*b^3+3*a^5*b+31574764*a^5*b^5*c^5+6*a^6*b^3*c^2", domain, MonomialOrder.DEFAULT),
                b = parse("7*a^2*b^6*c^3+a^5*b^4*c^4+31574764*a^6*c^3+5*a^6*b^2*c^2", domain, MonomialOrder.DEFAULT),
                gcd = parse("9*c^4+31574766*a*b^2+2*a^2*b*c^2+31574768*a^2*b^3*c^6+9*a^3*b^2*c^3", domain, MonomialOrder.DEFAULT);

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), interpolateGCD(asOverZp64(a), asOverZp64(b), asOverZp64(ZippelGCD(a, b)), getRandom()).monic().toBigPoly());
    }

    @Test // =====> testSparseInterpolation2   elapsed 195us
    public void testSparseInterpolation2() throws Exception {
        IntegersZp domain = new IntegersZp(24001871);
        MultivariatePolynomial<BigInteger>
                a = parse("3*b^4*c^2+7*b^4*c^3+4*a*b^5*c^3+6*a^4*b^6+24001865*a^5*b", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                b = parse("5*a*c^4+9*a^4*b^4*c^2+9*a^6*b*c^6", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                gcd = parse("5*a*b^2*c^2+a*b^2*c^4+24001866*a*b^4*c^3 + 1", domain, MonomialOrder.DEFAULT, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(asOverZp64(gcd).monic(), interpolateGCD(asOverZp64(a), asOverZp64(b), asOverZp64(gcd), getRandom()).monic());
    }

    @Test // =====> testSparseInterpolation3   elapsed 253us
    public void testSparseInterpolation3() throws Exception {
        IntegersZp domain = new IntegersZp(17312587);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a^3*c^6+9*a^5*b^2*c^3+7*a^5*b^6*c^5+8*a^5*b^6*c^6+6*a^6*b^6*c", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                b = parse("17312581*c^6+5*a^2*b^2*c^6+3*a^4*b^6*c^4+2*a^5+4*a^5*b^3*c^6", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                gcd = parse("1 + 5*a^5*b*c^2+6*a^5*b^3*c^6+2*a^5*b^4*c^4", domain, MonomialOrder.DEFAULT, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        MultivariatePolynomialZp64 intrp = interpolateGCD(asOverZp64(a), asOverZp64(b), asOverZp64(gcd), getRandom());
        assertEquals(asOverZp64(gcd).monic(), intrp.monic());
    }

    @Test // =====> testSparseInterpolation4   elapsed 1780us
    public void testSparseInterpolation4() throws Exception {
        IntegersZp domain = new IntegersZp(27445993);
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b*c^3+8*a^3*c+8*a^4*b^2*c^4+8*a^4*b^6*c^6", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                b = parse("1 + a*b^6*c^2+6*a^2*b^3*c^3+27445990*a^3*b^6*c^2", domain, MonomialOrder.DEFAULT, "a", "b", "c"),
                gcd = parse("1 + 5*b*c^3+8*b^5*c+4*b^6+5*a*b^3+4*a^6*b^3*c^3", domain, MonomialOrder.DEFAULT, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        MultivariatePolynomialZp64 la = asOverZp64(a), lb = asOverZp64(b);
        MultivariatePolynomialZp64 lgcd = ZippelGCD(la, lb);
        MultivariatePolynomialZp64 intrp = interpolateGCD(la, lb, lgcd, getRandom());
        assertEquals(lgcd.monic(), intrp.monic());
    }

    @Test // =====> testSparseInterpolation_random1   elapsed 8s
    public void testSparseInterpolation_random1() throws Exception {
        int nIterations = its(1000, 2000);
        RandomGenerator rnd = getRandom();

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 5, 5, 15, 5, 15, rnd);
        for (int n = 0; n < nIterations; n++) {
            GCDSample<MonomialZp64, MultivariatePolynomialZp64> gcdTriplet = sampleData.nextSample(false, false);
            MultivariatePolynomialZp64 contentGCD = MultivariateGCD.contentGCD(gcdTriplet.a, gcdTriplet.b, 0, MultivariateGCD::PolynomialGCD);
            MultivariatePolynomialZp64 gcd = null, actual = null;
            try {
                MultivariatePolynomialZp64 la = divideExact(gcdTriplet.a, contentGCD);
                MultivariatePolynomialZp64 lb = divideExact(gcdTriplet.b, contentGCD);
                gcd = ZippelGCD(la, lb);
                if (la.isConstant() || lb.isConstant() || gcd.degree(0) == 0) {
                    --n;
                    continue;
                }
                actual = interpolateGCD(la, lb, gcd, rnd);
                assertEquals(gcd.monic(), actual.monic());
            } catch (Throwable thr) {
                System.out.println(gcdTriplet.domain);
                System.out.println(gcdTriplet.a);
                System.out.println(gcdTriplet.b);
                System.out.println(gcd);
                System.out.println(actual);
                throw thr;
            }
        }
    }

    @Test // =====> testSparseInterpolation5   elapsed 12ms
    public void testSparseInterpolation5() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b*c^3+8*a^3*c+8*a^4*b^2*c^4+8*a^4*b^6*c^6", Rings.Z, MonomialOrder.DEFAULT, "a", "b", "c"),
                b = parse("1 + a*b^6*c^2+6*a^2*b^3*c^3+27445990*a^3*b^6*c^2", Rings.Z, MonomialOrder.DEFAULT, "a", "b", "c"),
                gcd = parse("1 + 5*b*c^3+8*b^5*c+4*b^6+5*a*b^3+4*a^6*b^3*c^3", Rings.Z, MonomialOrder.DEFAULT, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        IntegersZp domain = new IntegersZp(27445993);

        MultivariatePolynomialZp64
                la = asOverZp64(a.setRing(domain)),
                lb = asOverZp64(b.setRing(domain));
        MultivariatePolynomialZp64 skeleton = ZippelGCD(la, lb);

        IntegersZp domain1 = new IntegersZp(BigPrimes.nextPrime(37445993132451L));
        MultivariatePolynomialZp64
                la1 = asOverZp64(a.setRing(domain1)),
                lb1 = asOverZp64(b.setRing(domain1));

        MultivariatePolynomialZp64 gcd1 = ZippelGCD(la1, lb1);

        skeleton = skeleton.setRing(la1.ring);
        MultivariatePolynomialZp64 intrp = interpolateGCD(la1, lb1, skeleton, getRandom());
        assertEquals(gcd1.monic(), intrp.monic());
    }

    @Test // =====> testSparseInterpolation_random2   elapsed 8s
    public void testSparseInterpolation_random2() throws Exception {
        int nIterations = its(500, 1000);
        int badEvaluations = 0;
        RandomGenerator rnd = getRandom();
        GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> sampleData =
                new GCDSampleDataGeneric<>(Rings.Z, 3, 5, 5, 15, 5, 15, rnd);
        for (int n = 0; n < nIterations; n++) {
            GCDSample<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> gcdTriplet = sampleData.nextSample(false, false);
//            MultivariatePolynomial<BigInteger> contentGCD = MultivariateGCD.contentGCD(gcdTriplet.aCoFactor, gcdTriplet.bCoFactor, 0, MultivariateGCD::PolynomialGCD);
//            gcdTriplet = new GCDSample<>(divideExact(gcdTriplet.aCoFactor, contentGCD), divideExact(gcdTriplet.bCoFactor, contentGCD), divideExact(gcdTriplet.gcd, contentGCD));

            MultivariatePolynomialZp64 skeleton = null, gcd = null, actual = null;
            IntegersZp domain = null, domain1 = null;
            long seed = -1;
            try {

                MultivariatePolynomial<BigInteger> contentGCD = MultivariateGCD.contentGCD(gcdTriplet.a, gcdTriplet.b, 0, MultivariateGCD::PolynomialGCD);
                MultivariatePolynomial<BigInteger> a = divideExact(gcdTriplet.a, contentGCD);
                MultivariatePolynomial<BigInteger> b = divideExact(gcdTriplet.b, contentGCD);

                domain = new IntegersZp(getModulusRandom(20));
                MultivariatePolynomialZp64
                        la = asOverZp64(a.setRing(domain)),
                        lb = asOverZp64(b.setRing(domain));


                skeleton = ZippelGCD(la, lb);
                if (la.isConstant() || lb.isConstant() || skeleton.degree(0) == 0) {
                    --n;
                    continue;
                }

                domain1 = new IntegersZp(getModulusRandom(20));
                MultivariatePolynomialZp64
                        la1 = asOverZp64(a.setRing(domain1)),
                        lb1 = asOverZp64(b.setRing(domain1));

                gcd = ZippelGCD(la1, lb1);
                if (!gcd.sameSkeletonQ(skeleton)) {
                    --n;
                    continue;
                }

                rnd.setSeed(seed = rnd.nextLong());
                actual = interpolateGCD(la1, lb1, skeleton.setRing(la1.ring), rnd);
                if (actual == null) {
                    ++badEvaluations;
                    // bad evaluation point => try over
                    actual = interpolateGCD(la1, lb1, skeleton.setRing(la1.ring), rnd);
                }
                assertEquals(gcd.monic(), actual.monic());
            } catch (Throwable thr) {
                System.out.println(seed);
                System.out.println(gcdTriplet.domain);
                System.out.println(gcdTriplet.a);
                System.out.println(gcdTriplet.b);
                System.out.println(gcdTriplet.gcd);
                System.out.println(domain);
                System.out.println(domain1);
                System.out.println(skeleton);
                System.out.println(actual);
                throw thr;
            }
        }
        System.out.println("Bad evaluations: " + badEvaluations);
    }

    @Ignore
    @Test
    public void testSparseInterpolation6a() throws Exception {
        for (int i = 0; i < 10000; i++) {
            RandomGenerator rnd = PrivateRandom.getRandom();
            rnd.setSeed(i);
            IntegersZp domain = new IntegersZp(1049);
            MultivariatePolynomial<BigInteger>
                    a = parse("15*a*b^5*c^5*d^3+27*a^2*b^10*c^4*d+35*a^3*b^7*c^5*d^4+20*a^3*b^7*c^5*d^5+40*a^4*b^11*c^11*d^10+63*a^4*b^12*c^4*d^2+36*a^4*b^12*c^4*d^3+72*a^5*b^16*c^10*d^8+243*a^6*b^12*c^9*d^4+15*a^6*b^15*c*d^3+15*a^6*b^15*c^11*d^10+45*a^7*b^9*c^11*d^7+5*a^7*b^14*c^7*d^4+12*a^8*b*c*d^11+35*a^8*b^14*c^4*d^8+567*a^8*b^14*c^9*d^5+324*a^8*b^14*c^9*d^6+81*a^8*b^14*c^10*d^5+231*a^8*b^15*c^3*d^10+35*a^8*b^17*c*d^4+20*a^8*b^17*c*d^5+35*a^8*b^17*c^11*d^11+20*a^8*b^17*c^11*d^12+9*a^8*b^19*c^6*d^2+15*a^9*c^10*d^11+415*a^9*b^6*c^4*d^11+390*a^9*b^8*c^8*d^5+648*a^9*b^18*c^15*d^11+63*a^9*b^19*c^3*d^6+40*a^9*b^21*c^7*d^10+40*a^9*b^21*c^17*d^17+24*a^10*c^3*d^11+28*a^10*b^3*c*d^12+16*a^10*b^3*c*d^13+747*a^10*b^11*c^3*d^9+20*a^10*b^12*c^4*d^5+702*a^10*b^13*c^7*d^3+405*a^10*b^14*c^5*d^5+539*a^10*b^17*c^3*d^11+308*a^10*b^17*c^3*d^12+35*a^11*b^2*c^10*d^12+20*a^11*b^2*c^10*d^13+32*a^11*b^7*c^7*d^18+36*a^11*b^17*c^3*d^3+729*a^11*b^19*c^4*d^3+616*a^11*b^21*c^9*d^17+56*a^12*b^2*c^3*d^12+32*a^12*b^2*c^3*d^13+40*a^12*b^6*c^16*d^18+729*a^12*b^16*c^15*d^8+45*a^12*b^19*c^7*d^7+45*a^12*b^19*c^17*d^14+81*a^12*b^21*c^11*d^5+5*a^12*b^24*c^3*d^4+5*a^12*b^24*c^13*d^11+18*a^13*b^3*c^5*d^13+64*a^13*b^6*c^9*d^18+567*a^13*b^21*c^8*d^9+35*a^13*b^24*d^8+35*a^13*b^24*c^10*d^15+36*a^14*b^5*c^7*d^15+4*a^14*b^10*c^3*d^12+429*a^14*b^13*c^8*d^12+24*a^14*b^15*c^12*d^6+415*a^14*b^16*d^11+415*a^14*b^16*c^10*d^18+390*a^14*b^18*c^4*d^5+390*a^14*b^18*c^14*d^12+693*a^14*b^19*c^9*d^14+77*a^14*b^24*c^5*d^11+45*a^15*b^4*c^16*d^15+42*a^15*b^5*c^5*d^14+24*a^15*b^5*c^5*d^15+5*a^15*b^9*c^12*d^12+28*a^15*b^10*d^16+324*a^15*b^19*c^8*d^6+267*a^15*b^21*c^9*d^6+20*a^15*b^22*d^5+20*a^15*b^22*c^10*d^12+405*a^15*b^24*c*d^5+539*a^15*b^24*c^2*d^15+405*a^15*b^24*c^11*d^12+332*a^16*b^2*d^19+312*a^16*b^4*c^4*d^13+72*a^16*b^4*c^9*d^15+8*a^16*b^9*c^5*d^12+35*a^16*b^9*c^9*d^16+48*a^16*b^9*c^11*d^20+97*a^16*b^16*c^2*d^18+761*a^16*b^18*c^6*d^12+415*a^17*b*c^9*d^19+390*a^17*b^3*c^13*d^13+16*a^17*b^8*d^13+56*a^17*b^9*c^2*d^16+324*a^17*b^10*c*d^13+308*a^17*b^22*c^2*d^12+992*a^17*b^24*c^3*d^12+664*a^18*b*c^2*d^19+624*a^18*b^3*c^6*d^13+20*a^18*b^7*c^9*d^13+405*a^18*b^9*c^10*d^13+32*a^19*b^7*c^2*d^13+54*a^19*b^7*c^11*d^17+648*a^19*b^9*c^3*d^13+6*a^19*b^12*c^7*d^14+42*a^20*b^12*c^4*d^18+498*a^21*b^4*c^4*d^21+468*a^21*b^6*c^8*d^15+24*a^22*b^10*c^4*d^15+486*a^22*b^12*c^5*d^15", domain, MonomialOrder.DEFAULT),
                    b = parse("12*c^6*d^2+234*b^8*c^8*d^6+28*a^2*b^2*c^6*d^3+16*a^2*b^2*c^6*d^4+3*a^2*b^7*c^5*d+546*a^2*b^10*c^8*d^7+312*a^2*b^10*c^8*d^8+32*a^3*b^6*c^12*d^9+624*a^3*b^14*c^14*d^13+7*a^4*b^9*c^5*d^2+4*a^4*b^9*c^5*d^3+8*a^5*b^13*c^11*d^8+12*a^6*b*c^6+36*a^6*b^4*c^12*d^6+4*a^6*b^9*c^8*d^3+702*a^6*b^12*c^14*d^10+78*a^6*b^17*c^10*d^7+28*a^7*b^9*c^5*d^7+546*a^7*b^17*c^7*d^11+3*a^8*c*d^3+332*a^8*b*c^5*d^10+28*a^8*b^3*c^6*d+16*a^8*b^3*c^6*d^2+312*a^8*b^3*c^9*d^4+180*a^8*b^9*c^7*d^14+9*a^8*b^11*c^11*d^5+839*a^8*b^11*c^11*d^8+a^8*b^16*c^7*d^2+16*a^9*b^7*c^5*d^4+32*a^9*b^7*c^12*d^7+324*a^9*b^9*c^6*d^4+312*a^9*b^15*c^7*d^8+7*a^9*b^16*c^4*d^6+24*a^9*b^17*c^8*d^8+7*a^10*b^2*c*d^4+4*a^10*b^2*c*d^5+83*a^10*b^8*c^4*d^9+78*a^10*b^10*c^8*d^3+8*a^11*b^6*c^7*d^10+4*a^11*b^14*c^4*d^3+81*a^11*b^16*c^5*d^3+36*a^12*b^5*c^12*d^4+4*a^12*b^10*c^8*d+28*a^13*b^10*c^5*d^5+332*a^14*b^2*c^5*d^8+9*a^14*b^4*c^7*d^7+312*a^14*b^4*c^9*d^2+a^14*b^9*c^3*d^4+16*a^15*b^8*c^5*d^2+7*a^15*b^9*d^8+324*a^15*b^10*c^6*d^2+83*a^16*b*d^11+78*a^16*b^3*c^4*d^5+4*a^17*b^7*d^5+81*a^17*b^9*c*d^5", domain, MonomialOrder.DEFAULT),
                    base = parse("3*c+7*a^2*b^2*c*d+4*a^2*b^2*c*d^2+8*a^3*b^6*c^7*d^7+9*a^6*b^4*c^7*d^4+a^6*b^9*c^3*d+7*a^7*b^9*d^5+17492158*a^8*b*d^8+17492153*a^8*b^3*c^4*d^2+4*a^9*b^7*d^2+17492156*a^9*b^9*c*d^2", domain, MonomialOrder.DEFAULT);

            MultivariatePolynomialZp64
                    la = asOverZp64(a),
                    lb = asOverZp64(b),
                    skeleton = asOverZp64(base);
            MultivariatePolynomialZp64 lgcd = ZippelGCD(la, lb);
            MultivariatePolynomialZp64 intrp = null;
            try {
                rnd.setSeed(i);
                intrp = interpolateGCD(la, lb, skeleton, rnd);
                if (intrp != null)
                    assertEquals(lgcd.monic(), intrp.monic());
            } catch (Throwable e) {
                System.out.println(i);
                System.out.println(lgcd);
                System.out.println(intrp);
            }
        }
    }

    @Test // =====> testSparseInterpolation6   elapsed 11ms
    public void testSparseInterpolation6() throws Exception {
        RandomGenerator rnd = PrivateRandom.getRandom();
        rnd.setSeed(743);
        IntegersZp domain = new IntegersZp(1049);
        MultivariatePolynomial<BigInteger>
                a = parse("15*a*b^5*c^5*d^3+27*a^2*b^10*c^4*d+35*a^3*b^7*c^5*d^4+20*a^3*b^7*c^5*d^5+40*a^4*b^11*c^11*d^10+63*a^4*b^12*c^4*d^2+36*a^4*b^12*c^4*d^3+72*a^5*b^16*c^10*d^8+243*a^6*b^12*c^9*d^4+15*a^6*b^15*c*d^3+15*a^6*b^15*c^11*d^10+45*a^7*b^9*c^11*d^7+5*a^7*b^14*c^7*d^4+12*a^8*b*c*d^11+35*a^8*b^14*c^4*d^8+567*a^8*b^14*c^9*d^5+324*a^8*b^14*c^9*d^6+81*a^8*b^14*c^10*d^5+231*a^8*b^15*c^3*d^10+35*a^8*b^17*c*d^4+20*a^8*b^17*c*d^5+35*a^8*b^17*c^11*d^11+20*a^8*b^17*c^11*d^12+9*a^8*b^19*c^6*d^2+15*a^9*c^10*d^11+415*a^9*b^6*c^4*d^11+390*a^9*b^8*c^8*d^5+648*a^9*b^18*c^15*d^11+63*a^9*b^19*c^3*d^6+40*a^9*b^21*c^7*d^10+40*a^9*b^21*c^17*d^17+24*a^10*c^3*d^11+28*a^10*b^3*c*d^12+16*a^10*b^3*c*d^13+747*a^10*b^11*c^3*d^9+20*a^10*b^12*c^4*d^5+702*a^10*b^13*c^7*d^3+405*a^10*b^14*c^5*d^5+539*a^10*b^17*c^3*d^11+308*a^10*b^17*c^3*d^12+35*a^11*b^2*c^10*d^12+20*a^11*b^2*c^10*d^13+32*a^11*b^7*c^7*d^18+36*a^11*b^17*c^3*d^3+729*a^11*b^19*c^4*d^3+616*a^11*b^21*c^9*d^17+56*a^12*b^2*c^3*d^12+32*a^12*b^2*c^3*d^13+40*a^12*b^6*c^16*d^18+729*a^12*b^16*c^15*d^8+45*a^12*b^19*c^7*d^7+45*a^12*b^19*c^17*d^14+81*a^12*b^21*c^11*d^5+5*a^12*b^24*c^3*d^4+5*a^12*b^24*c^13*d^11+18*a^13*b^3*c^5*d^13+64*a^13*b^6*c^9*d^18+567*a^13*b^21*c^8*d^9+35*a^13*b^24*d^8+35*a^13*b^24*c^10*d^15+36*a^14*b^5*c^7*d^15+4*a^14*b^10*c^3*d^12+429*a^14*b^13*c^8*d^12+24*a^14*b^15*c^12*d^6+415*a^14*b^16*d^11+415*a^14*b^16*c^10*d^18+390*a^14*b^18*c^4*d^5+390*a^14*b^18*c^14*d^12+693*a^14*b^19*c^9*d^14+77*a^14*b^24*c^5*d^11+45*a^15*b^4*c^16*d^15+42*a^15*b^5*c^5*d^14+24*a^15*b^5*c^5*d^15+5*a^15*b^9*c^12*d^12+28*a^15*b^10*d^16+324*a^15*b^19*c^8*d^6+267*a^15*b^21*c^9*d^6+20*a^15*b^22*d^5+20*a^15*b^22*c^10*d^12+405*a^15*b^24*c*d^5+539*a^15*b^24*c^2*d^15+405*a^15*b^24*c^11*d^12+332*a^16*b^2*d^19+312*a^16*b^4*c^4*d^13+72*a^16*b^4*c^9*d^15+8*a^16*b^9*c^5*d^12+35*a^16*b^9*c^9*d^16+48*a^16*b^9*c^11*d^20+97*a^16*b^16*c^2*d^18+761*a^16*b^18*c^6*d^12+415*a^17*b*c^9*d^19+390*a^17*b^3*c^13*d^13+16*a^17*b^8*d^13+56*a^17*b^9*c^2*d^16+324*a^17*b^10*c*d^13+308*a^17*b^22*c^2*d^12+992*a^17*b^24*c^3*d^12+664*a^18*b*c^2*d^19+624*a^18*b^3*c^6*d^13+20*a^18*b^7*c^9*d^13+405*a^18*b^9*c^10*d^13+32*a^19*b^7*c^2*d^13+54*a^19*b^7*c^11*d^17+648*a^19*b^9*c^3*d^13+6*a^19*b^12*c^7*d^14+42*a^20*b^12*c^4*d^18+498*a^21*b^4*c^4*d^21+468*a^21*b^6*c^8*d^15+24*a^22*b^10*c^4*d^15+486*a^22*b^12*c^5*d^15", domain, MonomialOrder.DEFAULT),
                b = parse("12*c^6*d^2+234*b^8*c^8*d^6+28*a^2*b^2*c^6*d^3+16*a^2*b^2*c^6*d^4+3*a^2*b^7*c^5*d+546*a^2*b^10*c^8*d^7+312*a^2*b^10*c^8*d^8+32*a^3*b^6*c^12*d^9+624*a^3*b^14*c^14*d^13+7*a^4*b^9*c^5*d^2+4*a^4*b^9*c^5*d^3+8*a^5*b^13*c^11*d^8+12*a^6*b*c^6+36*a^6*b^4*c^12*d^6+4*a^6*b^9*c^8*d^3+702*a^6*b^12*c^14*d^10+78*a^6*b^17*c^10*d^7+28*a^7*b^9*c^5*d^7+546*a^7*b^17*c^7*d^11+3*a^8*c*d^3+332*a^8*b*c^5*d^10+28*a^8*b^3*c^6*d+16*a^8*b^3*c^6*d^2+312*a^8*b^3*c^9*d^4+180*a^8*b^9*c^7*d^14+9*a^8*b^11*c^11*d^5+839*a^8*b^11*c^11*d^8+a^8*b^16*c^7*d^2+16*a^9*b^7*c^5*d^4+32*a^9*b^7*c^12*d^7+324*a^9*b^9*c^6*d^4+312*a^9*b^15*c^7*d^8+7*a^9*b^16*c^4*d^6+24*a^9*b^17*c^8*d^8+7*a^10*b^2*c*d^4+4*a^10*b^2*c*d^5+83*a^10*b^8*c^4*d^9+78*a^10*b^10*c^8*d^3+8*a^11*b^6*c^7*d^10+4*a^11*b^14*c^4*d^3+81*a^11*b^16*c^5*d^3+36*a^12*b^5*c^12*d^4+4*a^12*b^10*c^8*d+28*a^13*b^10*c^5*d^5+332*a^14*b^2*c^5*d^8+9*a^14*b^4*c^7*d^7+312*a^14*b^4*c^9*d^2+a^14*b^9*c^3*d^4+16*a^15*b^8*c^5*d^2+7*a^15*b^9*d^8+324*a^15*b^10*c^6*d^2+83*a^16*b*d^11+78*a^16*b^3*c^4*d^5+4*a^17*b^7*d^5+81*a^17*b^9*c*d^5", domain, MonomialOrder.DEFAULT);

        MultivariatePolynomialZp64
                la = asOverZp64(a),
                lb = asOverZp64(b);

        MultivariatePolynomialZp64 lgcd = ZippelGCD(la, lb);
        assertTrue(dividesQ(la, lgcd));
        assertTrue(dividesQ(lb, lgcd));


        rnd.setSeed(701);
        MultivariatePolynomialZp64 intrp = interpolateGCD(la, lb, lgcd, rnd);
        if (intrp != null)
            assertEquals(lgcd.monic(), intrp.monic());
    }

    @Test // =====> testSparseInterpolation7   elapsed 88ms
    public void testSparseInterpolation7() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("7*b*c^4*d^6+9*a^2*b*c^8*d^7*e^7+7*a^2*b^6*c^8*d^3*e+8*a^3*c^2*d^6*e^4+7*a^3*b^5*c^7*d^6*e^4+a^4*b^3*c^5*d^8*e^5+25732656*a^5*c^4*d^2*e^3+9*a^5*b^2*c^6*d^5*e^4+25732652*a^6*b^3*c*d*e+25732656*a^7*b^3*c^8*d+a^7*b^3*c^8*d^2*e"),
                b = parse("1 + 25732655*a^9*b^8*c^12*d^18*e^13+9*a^11*b^16*c^19*d^11*e^13+4*a^16*b^20*c^17*d^3*e^4+2*a^20*b^10*d^3*e^13+4*a^20*b^11*c^13*d^17*e^9"),
                gcd = parse("1 + 3*a^2*b^17*c^14*d^6*e^14+4*a^3*b^14*c^15*d^10*e^8+25732658*a^5*b^17*c^10*d^9*e^12+8*a^6*b^10*c^4*d^3*e^10+25732659*a^6*b^10*c^7*d^5*e^15+a^7*b^2*c^3*d+6*a^9*b^9*c^10*d^6*e^5+3*a^11*b^15*c^7*d^17*e^15+25732652*a^13*b^3*c^5*d^13*e^11+2*a^13*b^12*d^2*e^16+9*a^15*b^2*c^2*d^5*e^4+2*a^15*b^2*c^14*d^14*e^14+a^15*b^13*c^8*e^12+a^16*b*c^10*d^13*e^10");

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        RandomGenerator rnd = getRandom();

        IntegersZp domain = new IntegersZp(806213L);
        MultivariatePolynomialZp64
                la = asOverZp64(a.setRing(domain)),
                lb = asOverZp64(b.setRing(domain));

        MultivariatePolynomialZp64 skeleton = ZippelGCD(la, lb);

        IntegersZp domain1 = new IntegersZp(755899L);
        MultivariatePolynomialZp64
                la1 = asOverZp64(a.setRing(domain1)),
                lb1 = asOverZp64(b.setRing(domain1));

        MultivariatePolynomialZp64 gcd0 = ZippelGCD(la1, lb1);
        System.out.println(gcd0.sameSkeletonQ(skeleton));

        rnd.setSeed(-7756222446675659124L);
        MultivariatePolynomialZp64 actual = interpolateGCD(la1, lb1, skeleton.setRing(la1.ring), rnd);
        if (actual == null) {
            System.out.println("bad evaluation");
            // bad evaluation point => try over
            actual = interpolateGCD(la1, lb1, skeleton.setRing(la1.ring), rnd);
        }
        assertEquals(gcd0.monic(), actual.monic());
    }

    @Test // =====> testModularGCD1   elapsed 12ms
    public void testModularGCD1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("a + 17*b + 2*c"),
                b = parse("3*a + b - c"),
                gcd = parse("1273465812736485821734523745*a*b - 21475715234*b - c");
        assertEquals(gcd, ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test // =====> testModularGCD2   elapsed 3ms
    public void testModularGCD2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("1234324234*a + 12317*b + 2*c"),
                b = parse("3*a + 143423423423423412314*b - c"),
                gcd = parse("1273465812736485821734523745*a*b - 21475715234*b - 143423423423423412314123123*c");
        assertEquals(gcd, ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test // =====> testModularGCD3   elapsed 10ms
    public void testModularGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("5*b^6*c^15*d^3+4*b^8*c*d^11+17492152*a^2*b^8*c^15*d^10+8*a^2*b^10*d^11+9*a^3*b^2*c^10*d+5*a^4*b*c^5*d^3+6*a^4*b^13*c^3*d^13+17492156*a^8*b^6*c^12*d^4+5*a^9*b^9*d^11+5*a^10*b^6*c^15*d^10"),
                b = parse("b^8*d^3+a^4*b^2*c^7*d+4*a^5*d^2+4*a^5*b^6*c+17492153*a^7*c^8*d^6"),
                gcd = parse("7*a^2*b^7*c^9*d^6+17492158*a^2*b^8*c*d^9+4*a^2*b^9*c^7*d^3+3*a^3*d+7*a^3*b^2*c^2*d^2+4*a^3*b^2*c^2*d^3+17492156*a^3*b^9*c^9*d^3+a^5*b^6*c^9*d^2+17492153*a^6*b^8*c^3*d^3+8*a^9*b^3*c^6*d^8+9*a^9*b^6*c^4*d^5");
        assertEquals(gcd, ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test // =====> testModularGCD4   elapsed 15ms
    public void testModularGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("4*a*b^7*c*d^2+a^2*b^6*c^8*d^4+7*a^3*b^5*c^6*d^4+5*a^4*b^3*c*d^7+6*a^5*b^4*c^7+8*a^7*c^8*d+a^8*b^5*c^3*d^2"),
                b = parse("25987600*b^18*c^17*d^14+25987597*a*b^9*c^9*d^2+2*a^2*b^7*c^12*d^7+4*a^4*b^14*c^11*d^2+6*a^6*b^2*d+2*a^6*b^3*c^16*d^14+5*a^9*b^17*c^16*d^2+a^14*b^18*c^17*d^4"),
                gcd = parse("25987593*a^4*c^4*d^4+9*a^6*c^3*d^10+7*a^6*b^14*c^4*d^7+8*a^7*b^9*c^13*d+7*a^9*b^2*c^13*d^4+2*a^10*b^6*c^9*d^7+2*a^11*b^5*c^7*d^3+2*a^11*b^12*c^13*d^14+7*a^14*b^8*c^14*d^3+6*a^14*b^13*c^4*d^11");
        assertEquals(gcd, ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test // =====> testModularGCD5   elapsed 14ms
    public void testModularGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a*b^5*c^10*d^7*e^16+2*a^4*b^3*c^9*d^6*e^8+5*a^4*b^6*c^16*d^11*e^2+a^4*b^13*d^5*e^6+30844060*a^5*b*c^9*d^8*e^12+4*a^8*b*c^17*d^11*e^3+9*a^8*b^13*c^16*d^17*e^11+a^9*b^2*c^2*d^10*e^14+5*a^9*b^6*c^3*d^7*e^4+7*a^9*b^8*c^3*d^16*e^2+9*a^14*b^5*c^2*d^3*e^16"),
                b = parse("7*b^6*c^18*d^5*e+30844053*a^2*b^8*c^10*d^8*e^6+a^3*b^14*c^4*d^11*e^7+a^4*b^10*c*d^15*e^18+3*a^15*b^9*c^3*e^11+5*a^18*b^13*c^16*d^15*e^15"),
                gcd = parse("9*a^3*b^11*c^7*d^4*e^6+30844059*a^5*b^6*c^15*d^8*e^10+8*a^5*b^10*c^15*d^2*e^9+5*a^10*b^11*c^7*d^9*e^16+2*a^13*b^3*c^13*d^6*e^2+30844060*a^14*b^3*c^6*d^3*e^13+30844055*a^14*b^6*c^4*d^13+30844055*a^14*b^17*c^2*d^8*e^13+2*a^17*b^5*c^7*d*e^11");
        assertTrue(dividesQ(ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test // =====> testModularGCD6   elapsed 903ms
    public void testModularGCD6() throws Exception {
        for (int i = 46; i < 100; i++) {
            PrivateRandom.getRandom().setSeed(46);
            MultivariatePolynomial<BigInteger>
                    a = parse("8*a*c^5*d^10*e^5+31118523*a*b^3*c^5*d^8*e^10+a^2*b^7*c*d*e^12+4*a^2*b^8*c*d^9*e^10+31118524*a^3*b^5*c^14*d^5*e^13+31118529*a^4*b^3*c^12*d^6*e^8+3*a^5*b^4*d^11*e^9+31118526*a^5*b^8*c^6*d^12*e+4*a^7*b^13*c^11*d^3+31118529*a^9*b^12*c^4*d^2*e^11+5*a^11*b^9*c^2*d*e^11+8*a^13*b^13*c^7*d^2*e^8+8*a^14*b^5*c^14*d^6*e^4"),
                    b = parse("31118526*c^3*d^4*e^2+31118530*b^4*c^6*d^5*e^6+5*a*b*c^4*d^4*e^3+31118527*a*b^3*d*e^2+31118525*a^2*b*c^7*d*e^4+5*a^2*b^4*c^8*d^2*e^5+6*a^2*b^6*d^7*e^5+9*a^2*b^7*c^8*d*e^5+4*a^4*b^6*e^7+3*a^5*b^2*c^6*d^4*e^3+31118529*a^7*b*c^2*d^5*e^8+8*a^7*b^3*c^3*d^4*e^5+7*a^8*b*c^2*d^5*e^8+6*a^8*b^3*c^3*d^5*e^3"),
                    gcd = parse("2*c^3*d*e^5+31118524*b^6*c^2*d^3*e^4+31118528*a^2*b^3*c^2*d^3+7*a^3*b*c^3*d^2+5*a^3*b^3*c^4*d^5*e^2+31118527*a^4*c^2*d^3+7*a^4*b*c*d*e^4+9*a^4*b*c^6*d^3*e^4+5*a^5*d^2*e^2+4*a^6*b^2*c^4*e+7*a^6*b^3*c^5*d^4*e^3");
            assertTrue(dividesQ(ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
        }
    }

    @Test // =====> testModularGCD7   elapsed 249us
    public void testModularGCD7() throws Exception {
        PrivateRandom.getRandom().setSeed(46);
        MultivariatePolynomial<BigInteger>
                a = parse("5*b*d^4+2*a^3*b*c*d^4+a^5*b^2*c^2*d^6+6*a^5*b^5*c^3*d^6+4*a^6*b*c^2+8*a^6*b^5*c^5*d^5"),
                b = parse("8*a*b*c^3*d^6+4*a*b^4*c*d+4*a*b^5*c*d^3+3*a^3*b^4*c^2"),
                gcd = parse("5*a^7*b^2*c^13*d^4");
        assertTrue(dividesQ(ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test // =====> testModularGCD8   elapsed 6ms
    public void testModularGCD8() throws Exception {
        PrivateRandom.getRandom().setSeed(48);
        MultivariatePolynomial<BigInteger>
                a = parse("8*a*b^19*c^11+8*a^3*b^4*c^9+7*a^3*b^10*c^12+3*a^5*b^14*c^21+7*a^9*b^21*c+8*a^10*b^8*c^5+a^14*b^21*c^12+15477328*a^21*b^20*c^8"),
                b = parse("15477335*b^8*c^4+7*b^9*c^8+15477332*a^3*b^13*c^4+15477335*a^9*b^13*c^6+15477328*a^12*c^9"),
                gcd = parse("15477332*a^10*b^13*c^5+7*a^14*b^5*c^3+6*a^19*b^12*c^5+2*a^19*b^12*c^13+15477329*a^20*b*c^19+15477332*a^20*b^8*c^12+7*a^21*b^8*c^2");
        assertTrue(dividesQ(ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test // =====> testModularGCD9   elapsed 3ms
    public void testModularGCD9() throws Exception {
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("428678675174997*b - 576309141757314*c - 1799929908190992*b*c + 43581966762456*b^2*c + 2155012404966050*c^2 + 1356161027210220*b*c^2 - 162945010788840*b^2*c^2 - 579102861059775*b^6*c^2 + 667785318790236*b^5*c^3 + 569898197386650*b^6*c^3 - 41635029186864*b^7*c^3", vars),
                b = parse("-c", vars),
                gcd = parse("-2287341106126463750*b^7*c^5 + 1532098182980478300*b^8*c^5 - 946030127950514950*b^4*c^6 + 633666328207734108*b^5*c^6 + 723818682898978700*b^6*c^6 - 2410587259891460925*b^7*c^6 + 67657454221929000*b^8*c^6 + 299366928422627212*b^3*c^7 - 997003974528718653*b^4*c^7 + 27982704417344040*b^5*c^7 + 31963832309981000*b^6*c^7 - 101990974478857750*b^7*c^7 + 13220043258527560*b^3*c^8 - 42182835947642390*b^4*c^8", vars);
        assertTrue(dividesQ(ZippelGCDInZ(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }
//
//    @Test
//    public void testModularGCD9a() throws Exception {
//        String[] vars = {"a", "b", "c"};
//        IntegersZp64 ring = new IntegersZp64(1033);
//        MultivariatePolynomialZp64
//                a = MultivariatePolynomialZp64.parse("788*b^3+694*b^4+812*b^5+539*a*b^2+1023*a*b^3+681*a*b^4+666*a*b^5+441*a^2*b+625*a^2*b^2+223*a^2*b^3+574*a^2*b^4+431*a^2*b^5+881*a^3*b+181*a^3*b^2+337*a^3*b^3+988*a^3*b^4+504*a^3*b^5+584*a^4*b+388*a^4*b^2+948*a^4*b^3+263*a^4*b^4+290*a^5+842*a^5*b+669*a^5*b^2+689*a^5*b^3+260*a^5*b^4+185*a^5*b^5+186*a^5*b^6+577*a^6+682*a^6*b+963*a^6*b^2+821*a^6*b^3+636*a^6*b^4+285*a^6*b^5+771*a^6*b^6+221*a^7*b+847*a^7*b^2+346*a^7*b^3+360*a^7*b^4+778*a^7*b^5+937*a^7*b^6+825*a^8*b^3+959*a^8*b^4+663*a^8*b^5+540*a^8*b^6+471*a^9*b^3+95*a^9*b^4+823*a^9*b^5+125*a^10*b^2+47*a^10*b^3+268*a^10*b^4+745*a^10*b^5+409*a^11*b^2+179*a^11*b^3+315*a^11*b^4+587*a^11*b^5+819*a^12*b^3+695*a^12*b^4", ring, vars),
//                b = MultivariatePolynomialZp64.parse("427*b^2+742*b^3+733*a*b+544*a*b^2+306*a*b^3+258*a^2*b+533*a^2*b^2+248*a^3*b+160*a^3*b^2+133*a^4+599*a^4*b+918*a^4*b^2+774*a^5+566*a^5*b", ring, vars),
//                skeleton = MultivariatePolynomialZp64.parse("943*b^2+976*b^3+899*a*b+212*a*b^2+764*a*b^3+26*a^2*b+274*a^2*b^2+921*a^3*b+189*a^3*b^2+866*a^4+265*a^4*b+955*a^4*b^2+548*a^5+858*a^5*b", ring, vars),
//                content = MultivariatePolynomialZp64.parse("673+b", ring, vars);
//
//        System.out.println(b.sameSkeletonQ(skeleton));
//        System.out.println(b.clone().setAllCoefficientsToUnit());
//
////        a = divideExact(a, content);
//        b = divideExact(b, content);
//
//        skeleton = skeleton.setAllCoefficientsToUnit();
//        content = content.setAllCoefficientsToUnit();
//        System.out.println(skeleton);
//        System.out.println(content);
//        skeleton = divideSkeletonExact(skeleton, content);
//    }

    @Test // =====> testModularGCD_random1   elapsed 12s
    public void testModularGCD_random1() throws Exception {
        int nIterations = its(1000, 3000);
        RandomGenerator rnd = getRandom();
        GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> sampleData =
                new GCDSampleDataGeneric<>(Rings.Z, 3, 5, 5, 15, 5, 15, rnd);

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd", MultivariateGCD::ZippelGCDInZ));
    }

    @Test // =====> testModularGCD_sparse_variables_random   elapsed 19s
    public void testModularGCD_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>>
                sampleData =
                boundCoefficients(
                        fixVariables(new GCDSampleDataGeneric<>(Rings.Z, nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars),
                        BigInteger.valueOf(100));

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("ZippelGCDInZ", MultivariateGCD::ZippelGCDInZ));
    }

    @Ignore
    @Test
    public void testModularGCD_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b", "c", "d", "e"};
        Ring<BigInteger> ring = Rings.Z;
        MultivariatePolynomial<BigInteger>
                a = parse("2147483167*a^4*b^60*c^57*d^26*e+44*a^8*b^39*c^67*d^22*e^17+38*a^32*b^6*c^13*d^10*e^3+357*a^36*b^34*c^60*d^2*e^59+563*a^42*b^41*c^45*d^52*e^14+257*a^44*b^68*c^43*d^2*e^73+613*a^48*b^55*c^22*d^32*e^19+2147483093*a^52*b^26*c^4*d^72*e^32+19*a^52*b^40*c^26*d^45*e^55+639*a^55*b^72*c^55*d^65", ring, MonomialOrder.DEFAULT, vars),
                b = parse("2147483150*b^25*c^18*d^62*e^59+2147482723*a^4*b^5*c^65*d^26*e^7+261*a^15*b^60*c^59*d^63*e^53+394*a^27*b^22*c^34*d^54*e^13+952*a^39*b^48*c^17*d^54*e^16+243*a^60*b^15*c^3*d^51*e^46+40*a^61*b^56*c^39*d^40*e^21+555*a^62*b^20*c^20*d^60*e^47+627*a^67*b^8*c^22*d^67*e^61+447*a^70*b^59*c^71*d^24*e^5", ring, MonomialOrder.DEFAULT, vars),
                gcd = parse("35*a*b^36*c^74*d^62*e^51+376*a^2*b^28*c^64*e^53+893*a^6*b^13*c^60*d^44*e^42+23*a^8*b^71*c^40*d^36*e^11+783*a^20*b^28*c^12*d^31*e^68+2147482938*a^31*b^30*c^40*d^65*e^72+2147482960*a^31*b^49*c^38*d^71*e^55+737*a^47*b^15*c^71*d^13*e^72+868*a^53*b^30*c^40*d^29*e^46+898*a^61*b^71*c^13*d^50*e^66", ring, MonomialOrder.DEFAULT, vars);

        int[] newVariables = {4, 1, 3, 2, 0};
        a = renameVariables(a, newVariables);
        b = renameVariables(b, newVariables);
        gcd = renameVariables(gcd, newVariables);

        newVariables = MultivariateGCD.inversePermutation(new int[]{1, 3, 2, 0, 4});
        a = renameVariables(a, newVariables);
        b = renameVariables(b, newVariables);
        gcd = renameVariables(gcd, newVariables);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        System.out.println(Arrays.toString(a.degrees()));
        System.out.println(Arrays.toString(b.degrees()));
        System.out.println(a);
        System.out.println(b);

        for (int i = 0; i < 1000; i++) {
            System.out.println();
            long start = System.nanoTime();
            System.out.println(ZippelGCDInZ(a.clone().increment(), b));
            System.out.println(nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            assertTrue(dividesQ(ZippelGCDInZ(a, b), gcd));
            System.out.println(nanosecondsToString(System.nanoTime() - start));
        }


        //        7ms
        //        194ms
        //
        //        11ms
        //        161ms
        //
        //        8ms
        //        138ms
    }

    @Test(timeout = 10000) // =====> testGCDInput   elapsed 60ms
    public void testGCDInput() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b", "c", "d", "e"};
        IntegersZp64 ring = Rings.Zp64(1031);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2147483167*a^4*b^60*c^57*d^26*e+44*a^8*b^39*c^67*d^22*e^17+38*a^32*b^6*c^13*d^10*e^3+357*a^36*b^34*c^60*d^2*e^59+563*a^42*b^41*c^45*d^52*e^14+257*a^44*b^68*c^43*d^2*e^73+613*a^48*b^55*c^22*d^32*e^19+2147483093*a^52*b^26*c^4*d^72*e^32+19*a^52*b^40*c^26*d^45*e^55+639*a^55*b^72*c^55*d^65", ring, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("2147483150*b^25*c^18*d^62*e^59+2147482723*a^4*b^5*c^65*d^26*e^7+261*a^15*b^60*c^59*d^63*e^53+394*a^27*b^22*c^34*d^54*e^13+952*a^39*b^48*c^17*d^54*e^16+243*a^60*b^15*c^3*d^51*e^46+40*a^61*b^56*c^39*d^40*e^21+555*a^62*b^20*c^20*d^60*e^47+627*a^67*b^8*c^22*d^67*e^61+447*a^70*b^59*c^71*d^24*e^5", ring, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("35*a*b^36*c^74*d^62*e^51+376*a^2*b^28*c^64*e^53+893*a^6*b^13*c^60*d^44*e^42+23*a^8*b^71*c^40*d^36*e^11+783*a^20*b^28*c^12*d^31*e^68+2147482938*a^31*b^30*c^40*d^65*e^72+2147482960*a^31*b^49*c^38*d^71*e^55+737*a^47*b^15*c^71*d^13*e^72+868*a^53*b^30*c^40*d^29*e^46+898*a^61*b^71*c^13*d^50*e^66", ring, MonomialOrder.DEFAULT, vars);

        int[] newVariables = {4, 1, 3, 2, 0};
        a = renameVariables(a, newVariables);
        b = renameVariables(b, newVariables);
        gcd = renameVariables(gcd, newVariables);

        newVariables = MultivariateGCD.inversePermutation(new int[]{1, 3, 2, 0, 4});
        a = renameVariables(a, newVariables);
        b = renameVariables(b, newVariables);
        gcd = renameVariables(gcd, newVariables);

        a = a.divideOrNull(a.monomialContent());
        b = b.divideOrNull(b.monomialContent());
        gcd = gcd.divideOrNull(gcd.monomialContent());

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

//        System.out.println(Arrays.toString(a.degrees()));
//        System.out.println(Arrays.toString(b.degrees()));
//        System.out.println(Arrays.toString(gcd.degrees()));

        assertTrue(dividesQ(PolynomialGCD(a, b), gcd));
    }

    @Test // =====> testRationals1   elapsed 363ms
    public void testRationals1() throws Exception {
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("-(2/3)*a*b*c - (7/6)*a^3*c^4 + (2/3)*b^3", Rings.Q),
                b = parse("(2/3)*a^2*b*c + (1/6)*a^3*c^4 + (2/13)*c^3", Rings.Q),
                gcd = parse("(12/3)*a^2*b*c^2 - (11/6)*a^3*b*c^4 + (2/11)*c", Rings.Q);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertTrue(dividesQ(BrownGCD(a, b), gcd));
        assertTrue(dividesQ(ZippelGCD(a, b), gcd));
    }

    @Test // =====> testFiniteField1   elapsed 521ms
    public void testFiniteField1() throws Exception {
        for (int i = 0; i < 10; i++) {
            FiniteField<UnivariatePolynomialZp64> field = FiniteField.GF17p5;
            MultivariatePolynomial<UnivariatePolynomialZp64>
                    a = MultivariatePolynomial.zero(3, field, MonomialOrder.DEFAULT)
                    .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 2, 3, 4, 5).modulus(17)), 1, 1, 3))
                    .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 3, 2, 13).modulus(17)), 3, 2, 1))
                    .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 13, 12, 13).modulus(17)), 0, 2, 1)),
                    b = MultivariatePolynomial.zero(3, field, MonomialOrder.DEFAULT)
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5).modulus(17)), 1, 1, 13))
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 1, 2, 13).modulus(17)), 2, 2, 1))
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 113, 112, 13).modulus(17)), 10, 2, 1)),
                    gcd = MultivariatePolynomial.one(3, field, MonomialOrder.DEFAULT)
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5, 12).modulus(17)), 11, 1, 13))
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(11, 2, 1, 1, 2, 13).modulus(17)), 21, 2, 1))
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(17)), 10, 12, 1))
                            .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(17)), 0, 0, 1));

            a = a.clone().add(b).multiply(gcd);
            b = b.clone().subtract(gcd).multiply(gcd);
            long start = System.nanoTime();
            assertTrue(dividesQ(PolynomialGCD(a, b), gcd));
            System.out.println(nanosecondsToString(System.nanoTime() - start));
        }
    }

//
//    static <E> UnivariatePolynomial<MultivariatePolynomial<E>> asMUnivariate(MultivariatePolynomial<E> poly, int variable) {
//        int[] uDegrees = poly.degrees(variable);
//        MultivariatePolynomial<E>[] data = new MultivariatePolynomial[ArraysUtil.max(uDegrees) + 1];
//        for (int i = 0; i < data.length; i++)
//            data[i] = poly.createZero();
//        for (int degree : uDegrees)
//            data[degree] = poly.coefficientOf(variable, degree);
//        return UnivariatePolynomial.create(new MultivariateRing<>(poly), data);
//    }
//
//    static <E> MultivariatePolynomial<E> fromMUnivariate(UnivariatePolynomial<MultivariatePolynomial<E>> poly, int variable) {
//        MultivariatePolynomial<E> zero = poly.ring.getZero();
//        for (int i = 0; i <= poly.degree(); i++) {
//            if (poly.isZeroAt(i))
//                continue;
//
//            MultivariatePolynomial<E> term = poly.get(i);
//            zero.add(
//                    term.multiply(new Monomial<>(zero.nVariables, variable, i, zero.ring.getOne())));
//        }
//        return zero;
//    }
//
//    @Test
//    public void testPRS() throws Exception {
//        MultivariatePolynomial<BigInteger>
//                a = parse("a + 2*b"),
//                b = parse("a*b + 17*a*b^3"),
//                gcd = parse("1 + a*b - 2*a*b^2");
//
//        a = a.multiply(gcd);
//        b = b.multiply(gcd);
//        IntegersZp ring = new IntegersZp(3);
//        a = a.setRing(ring);
//        b = b.setRing(ring);
//
//
////        System.out.println(MultivariateGCD.PolynomialGCD(a, b));
//        int variable = 0;
////        System.out.println(a);
////        System.out.println(asMUnivariate(a, variable));
////        System.out.println(b);
////        System.out.println(asMUnivariate(b, variable));
//
//        UnivariatePolynomial<MultivariatePolynomial<BigInteger>> result = UnivariateGCD.SubresultantRemainders(asMUnivariate(a, variable), asMUnivariate(b, variable)).gcd();
//        System.out.println(result);
//        System.out.println(fromMUnivariate(result, variable));
//
//    }

    @Test // =====> testSmallDomain1   elapsed 44ms
    public void testSmallDomain1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("a + 2*b + c"),
                b = parse("a*b + 17*a*b^13 + a*c^2"),
                gcd = parse("1 + a^2*b^31*c + c^3*a^3*b - 2*a^5*b^2 - b*c^2");

        for (long modulus : new long[]{2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43}) {
            IntegersZp domain = new IntegersZp(modulus);
            MultivariatePolynomial<BigInteger>
                    a1 = a.clone().setRing(domain),
                    b1 = b.clone().setRing(domain),
                    gcd1 = gcd.clone().setRing(domain);
            a1 = a1.multiply(gcd1);
            b1 = b1.multiply(gcd1);

            assertEquals(gcd1.monic(), KaltofenMonaganSparseModularGCDInGF(a1, b1).monic());
        }
    }

    @Test // =====> testSmallDomain_random1   elapsed 19s
    public void testSmallDomain_random1() throws Exception {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int nIterations = its(200, 1000);
        RandomGenerator rnd = getRandom();

        lGCDSampleDataZp sampleData =
                new lGCDSampleDataZp(3, 5, 5, 15, 5, 15, rnd);
        sampleData.minModulusBits = 2;
        sampleData.maxModulusBits = 5;
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::KaltofenMonaganSparseModularGCDInGF),
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::KaltofenMonaganSparseModularGCDInGF));
    }

    @Test // =====> testSmallDomain_sparse_variables_random   elapsed 21s
    public void testSmallDomain_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(100, 1000);
        lGCDSampleDataZp source = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        source.minModulusBits = 2;
        source.maxModulusBits = 5;

        GCDSampleData<MonomialZp64, MultivariatePolynomialZp64>
                sampleData = filterZeros(fixVariables(source, nVars));

        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::KaltofenMonaganSparseModularGCDInGF),
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::KaltofenMonaganSparseModularGCDInGF));
    }

    @Test // =====> testSmallDomain_sparse_variables_random_a   elapsed 209ms
    public void testSmallDomain_sparse_variables_random_a() throws Exception {
        for (int i = 220; i < 250; ++i) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            IntegersZp64 domain = new IntegersZp64(29);
            String[] vars = IStringifier.defaultVars(100);
            MultivariatePolynomialZp64
                    a = MultivariatePolynomialZp64.parse("x4*x6^3*x7^4*x8*x9^3*x13^4*x17^2*x18^3*x21^5*x22*x23^3*x26*x27*x28^5*x30^2*x31^4*x32^3*x34*x37^4*x40^5*x41*x42^5*x44^2*x48^4*x49^5*x50^5*x52^5*x53^4*x56^2*x58*x59^5*x61*x64^3*x70^5*x72^3*x73^5*x75*x76^5*x79^3*x81*x83^2*x84*x85^5*x88*x90*x92^5*x93^5*x95^2*x97^5*x99+23*x4*x6^3*x7^4*x8*x9^3*x13^4*x18^3*x21^5*x22*x23^3*x26*x27*x28^5*x30^2*x31^4*x32^3*x34*x37^4*x40^5*x41*x42^5*x44^2*x48^4*x49^5*x50^6*x52^5*x53^4*x56^2*x58*x59^5*x61*x64^3*x70^5*x72^3*x73^5*x75*x76^5*x79^3*x81^3*x83^2*x84*x85^5*x88*x90*x92^5*x93^5*x95^2*x97^5*x99+17*x4*x6^3*x7^4*x8*x9^3*x13^4*x18^3*x21^5*x22*x23^3*x26*x27*x28^5*x30^2*x31^4*x32^3*x34*x37^4*x40^5*x41*x42^5*x44^2*x48^4*x49^5*x50^7*x52^5*x53^4*x56^2*x58*x59^5*x61*x64^3*x70^5*x72^3*x73^5*x75*x76^5*x79^3*x81^2*x83^2*x84*x85^5*x88*x90*x92^5*x93^5*x95^2*x97^5*x99+3*x4*x6^3*x7^4*x8*x9^3*x13^4*x17^2*x18^3*x21^5*x22*x23^3*x26*x27*x28^5*x30^2*x31^4*x32^3*x34*x37^4*x40^5*x41*x42^5*x44^2*x48^4*x49^5*x50^5*x52^5*x53^4*x56^2*x58*x59^5*x61*x64^3*x70^5*x72^3*x73^5*x75*x76^5*x79^3*x81^3*x83^2*x84*x85^5*x88*x90*x92^5*x93^5*x95^2*x97^5*x99+15*x4*x6^3*x7^4*x8*x9^3*x13^4*x17^2*x18^3*x21^5*x22*x23^3*x26*x27*x28^5*x30^2*x31^4*x32^3*x34*x37^4*x40^5*x41*x42^5*x44^2*x48^4*x49^5*x50^7*x52^5*x53^4*x56^2*x58*x59^5*x61*x64^3*x70^5*x72^3*x73^5*x75*x76^5*x79^3*x81^2*x83^2*x84*x85^5*x88*x90*x92^5*x93^5*x95^2*x97^5*x99", domain, vars),
                    b = MultivariatePolynomialZp64.parse("24*x1^3*x2^3*x3*x5^3*x10^2*x11*x12^4*x15^4*x16^2*x17^4*x19*x20^4*x33^5*x36^4*x38*x39^4*x43^4*x45^3*x46^2*x47^3*x50*x51^4*x55*x57*x60^3*x62^5*x63^4*x65^4*x66^3*x67^2*x68*x69^2*x71*x74*x77*x78^3*x80*x81*x87*x89*x96^5*x100^2+x1^3*x2^3*x3*x5^3*x10^2*x11*x12^4*x15^4*x16^2*x17^2*x19*x20^4*x33^5*x36^4*x38*x39^4*x43^4*x45^3*x46^2*x47^3*x50^2*x51^4*x55*x57*x60^3*x62^5*x63^4*x65^4*x66^3*x67^2*x68*x69^2*x71*x74*x77*x78^3*x80*x81^3*x87*x89*x96^5*x100^2+2*x1^3*x2^3*x3*x5^3*x10^2*x11*x12^4*x15^4*x16^2*x17^2*x19*x20^4*x33^5*x36^4*x38*x39^4*x43^4*x45^3*x46^2*x47^3*x50^3*x51^4*x55*x57*x60^3*x62^5*x63^4*x65^4*x66^3*x67^2*x68*x69^2*x71*x74*x77*x78^3*x80*x81^2*x87*x89*x96^5*x100^2+14*x1^3*x2^3*x3*x5^3*x10^2*x11*x12^4*x15^4*x16^2*x17^4*x19*x20^4*x33^5*x36^4*x38*x39^4*x43^4*x45^3*x46^2*x47^3*x50*x51^4*x55*x57*x60^3*x62^5*x63^4*x65^4*x66^3*x67^2*x68*x69^2*x71*x74*x77*x78^3*x80*x81^3*x87*x89*x96^5*x100^2+12*x1^3*x2^3*x3*x5^3*x10^2*x11*x12^4*x15^4*x16^2*x17^4*x19*x20^4*x33^5*x36^4*x38*x39^4*x43^4*x45^3*x46^2*x47^3*x50^3*x51^4*x55*x57*x60^3*x62^5*x63^4*x65^4*x66^3*x67^2*x68*x69^2*x71*x74*x77*x78^3*x80*x81^2*x87*x89*x96^5*x100^2", domain, vars),
                    gcd = MultivariatePolynomialZp64.parse("23*x17^2*x50+7*x50^2*x81^2+14*x50^3*x81+11*x17^2*x50*x81^2+26*x17^2*x50^3*x81", domain, vars),
                    actual = MultivariatePolynomialZp64.parse("2*x17^2*x50*x81+17*x50^2*x81^3+5*x50^3*x81^2+6*x17^2*x50*x81^3+x17^2*x50^3*x81^2", domain, vars);

            MultivariatePolynomialZp64 lResult = KaltofenMonaganSparseModularGCDInGF(a, b);
            checkConsistency(lResult);
            assertTrue(dividesQ(lResult, gcd));

            Conversions64bit.SWITCH_TO_64bit = false;
            MultivariatePolynomial<BigInteger> result = KaltofenMonaganSparseModularGCDInGF(a.toBigPoly(), b.toBigPoly());
            checkConsistency(result);
            assertTrue(dividesQ(result, gcd.toBigPoly()));

            assertEquals(result.size(), lResult.size());
        }
    }

    @Test // =====> testSmallDomain2   elapsed 9s
    public void testSmallDomain2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};

        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+2*c^3*d^2+2*b^3*c^3*d^3*e+a*c^3*d*e+2*a^2*b^3*c^2*d^2*e^3+a^2*b^3*c^3*e^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c^2*d^3*e^3+a*c^3*d*e^2+2*a^3*e^3+2*a^3*b^3*d*e^3+2*a^3*b^3*c*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*a*b^3*c+a^2*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*b^3*c^3*d^3*e+2*a*b^2*c*d^2*e^3+a*b^3*c^2*d*e^2+a^3*b^2*c^3*d^2", domain, vars),
        }, base = arr[0].createOne().multiply(arr);

        MultivariatePolynomialZp64 a = base;
        MultivariatePolynomialZp64 b = a.derivative(1);

        for (int i = 0; i < its(5, 5); i++) {
            timestamp();
            MultivariatePolynomialZp64 gcd = MultivariateGCD.KaltofenMonaganSparseModularGCDInGF(a, b);
            timeElapsed();

            assertTrue(dividesQ(a, gcd));
            assertTrue(dividesQ(b, gcd));
        }
    }

    @Test // =====> testSmallDomain3   elapsed 506us
    public void testSmallDomain3() throws Exception {
        IntegersZp64 domain = new IntegersZp64(5);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a*c^3*d^2+4*a*b*c^2*d", domain),
                b = MultivariatePolynomialZp64.parse("a*c^4*d+a^3*b^2*c*d^4+3*a^3*b^2*c^2*d^3+2*a^3*b^3*d^3", domain);

        for (int i = 0; i < 10; i++) {
            PrivateRandom.getRandom().setSeed(i);
            timestamp();
            assertTrue(PolynomialGCD(a, b).isMonomial());
            timeElapsed();
        }
    }

    @Test
    public void testSmallDomain7() {
        IntegersZp64 domain = new IntegersZp64(5);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*x^5+3*y^6+2*x^6+3*x*y^6+4*x^5*y^2+3*x^6*y+2*x*y^7+4*x^6*y^2+2*x^7*y+3*x^2*y^7+x^6*y^3+3*x^7*y^2+3*x^9+2*x^2*y^8+2*x^4*y^6+4*x^6*y^4+4*x^7*y^3+2*x^9*y+x*y^10+3*x^4*y^7+3*x^9*y^2+4*x^10*y+x^2*y^10+3*x^4*y^8+x^5*y^7+3*x^6*y^6+x^7*y^5+4*x^9*y^3+2*x^10*y^2+4*x^2*y^11+3*x^5*y^8+3*x^7*y^6+x^9*y^4+2*x^10*y^3+2*x^12*y+3*x^4*y^10+x^5*y^9+4*x^9*y^5+4*x^10*y^4+3*x^11*y^3+2*x^12*y^2+x^4*y^11+2*x^6*y^9+3*x^7*y^8+3*x^9*y^6+2*x^10*y^5+3*x^11*y^4+2*x^12*y^3+x^4*y^12+x^5*y^11+2*x^6*y^10+2*x^7*y^9+3*x^9*y^7+x^10*y^6+3*x^11*y^5+3*x^12*y^4+4*x^5*y^12+3*x^6*y^11+x^7*y^10+2*x^9*y^8+3*x^10*y^7+4*x^11*y^6+4*x^12*y^5+4*x^14*y^3+x^4*y^14+2*x^6*y^12+2*x^7*y^11+x^9*y^9+3*x^11*y^7+2*x^12*y^6+x^14*y^4+2*x^16*y^2+2*x^5*y^14+x^6*y^13+x^7*y^12+2*x^9*y^10+3*x^11*y^8+4*x^12*y^7+3*x^14*y^5+x^15*y^4+3*x^16*y^3+x^6*y^14+2*x^7*y^13+2*x^11*y^9+3*x^12*y^8+2*x^15*y^5+2*x^16*y^4+2*x^17*y^3+3*x^6*y^15+2*x^9*y^12+3*x^10*y^11+3*x^15*y^6+4*x^16*y^5+x^17*y^4+4*x^7*y^15+4*x^10*y^12+x^11*y^11+4*x^12*y^10+3*x^14*y^8+x^15*y^7+3*x^16*y^6+4*x^17*y^5+3*x^9*y^14+3*x^10*y^13+3*x^11*y^12+2*x^12*y^11+4*x^14*y^9+2*x^17*y^6+x^9*y^15+2*x^10*y^14+x^11*y^13+3*x^12*y^12+3*x^14*y^10+4*x^16*y^8+4*x^19*y^5+x^20*y^4+2*x^7*y^18+x^9*y^16+4*x^10*y^15+x^14*y^11+2*x^16*y^9+2*x^17*y^8+3*x^19*y^6+3*x^9*y^17+x^11*y^15+4*x^12*y^14+x^15*y^11+2*x^16*y^10+2*x^20*y^6+x^10*y^17+3*x^11*y^16+x^14*y^13+x^15*y^12+4*x^16*y^11+x^17*y^10+4*x^19*y^8+3*x^20*y^7+x^9*y^19+4*x^10*y^18+4*x^11*y^17+3*x^12*y^16+2*x^14*y^14+4*x^16*y^12+4*x^20*y^8+2*x^14*y^15+2*x^15*y^14+2*x^20*y^9+x^10*y^20+2*x^11*y^19+4*x^12*y^18+x^14*y^16+4*x^15*y^15+4*x^17*y^13+x^19*y^11+3*x^20*y^10+2*x^22*y^8+x^12*y^19+x^14*y^17+3*x^15*y^16+2*x^16*y^15+2*x^17*y^14+4*x^20*y^11+3*x^21*y^10+x^22*y^9+4*x^11*y^21+3*x^12*y^20+3*x^14*y^18+3*x^15*y^17+2*x^16*y^16+x^17*y^15+4*x^19*y^13+x^20*y^12+x^21*y^11+3*x^22*y^10+x^12*y^21+4*x^14*y^19+4*x^15*y^18+x^16*y^17+x^19*y^14+4*x^20*y^13+x^21*y^12+4*x^22*y^11+3*x^14*y^20+x^17*y^17+4*x^19*y^15+3*x^20*y^14+3*x^21*y^13+3*x^14*y^21+4*x^15*y^20+4*x^16*y^19+3*x^17*y^18+2*x^19*y^16+4*x^20*y^15+x^21*y^14+4*x^22*y^13+2*x^14*y^22+x^15*y^21+4*x^16*y^20+2*x^19*y^17+4*x^20*y^16+4*x^21*y^15+x^25*y^11+2*x^14*y^23+3*x^15*y^22+3*x^16*y^21+4*x^17*y^20+4*x^19*y^18+3*x^20*y^17+2*x^21*y^16+2*x^25*y^12+2*x^14*y^24+4*x^15*y^23+3*x^19*y^19+3*x^20*y^18+3*x^24*y^14+4*x^25*y^13+3*x^14*y^25+4*x^15*y^24+4*x^16*y^23+2*x^17*y^22+3*x^19*y^20+x^20*y^19+2*x^21*y^18+3*x^22*y^17+x^24*y^15+4*x^25*y^14+4*x^15*y^25+3*x^16*y^24+2*x^17*y^23+3*x^19*y^21+2*x^20*y^20+x^21*y^19+2*x^24*y^16+2*x^25*y^15+x^16*y^25+x^17*y^24+4*x^19*y^22+4*x^20*y^21+3*x^21*y^20+2*x^22*y^19+2*x^24*y^17+x^25*y^16+x^16*y^26+4*x^17*y^25+4*x^20*y^22+3*x^21*y^21+2*x^22*y^20+2*x^24*y^18+2*x^25*y^17+x^17*y^26+4*x^20*y^23+3*x^21*y^22+3*x^22*y^21+2*x^25*y^18+4*x^17*y^27+3*x^22*y^22+2*x^25*y^19+2*x^22*y^23", domain, "x", "y"),
                b = MultivariatePolynomialZp64.parse("x^5+4*y^6+2*x^5*y^2+x^7+4*x^2*y^6+4*x^8+x^3*y^6+2*x^7*y^2+x^8*y+4*x^3*y^7+3*x^8*y^2+3*x^9*y+2*x^4*y^7+2*x^8*y^3+2*x^9*y^2+x^10*y+3*x^11+3*x^4*y^8+4*x^5*y^7+2*x^6*y^6+3*x^8*y^4+x^9*y^3+2*x^11*y+2*x^3*y^10+3*x^6*y^7+2*x^10*y^3+3*x^11*y^2+2*x^12*y+4*x^4*y^10+3*x^6*y^8+3*x^7*y^7+x^8*y^6+4*x^9*y^5+4*x^10*y^4+4*x^11*y^3+x^12*y^2+x^4*y^11+x^5*y^10+4*x^7*y^8+2*x^9*y^6+x^11*y^4+x^12*y^3+3*x^14*y+3*x^6*y^10+3*x^7*y^9+3*x^10*y^6+4*x^11*y^5+2*x^12*y^4+x^13*y^3+3*x^14*y^2+4*x^15*y+x^6*y^11+4*x^8*y^9+2*x^9*y^8+x^10*y^7+3*x^11*y^6+x^12*y^5+x^13*y^4+3*x^14*y^3+3*x^15*y^2+x^6*y^12+3*x^7*y^11+4*x^8*y^10+3*x^9*y^9+2*x^10*y^8+3*x^11*y^7+3*x^12*y^6+x^13*y^5+2*x^14*y^4+2*x^15*y^3+2*x^7*y^12+x^8*y^11+4*x^9*y^10+x^10*y^9+2*x^11*y^8+4*x^12*y^7+3*x^13*y^6+x^14*y^5+x^15*y^4+4*x^16*y^3+x^6*y^14+4*x^8*y^12+3*x^9*y^11+x^11*y^9+x^13*y^7+3*x^14*y^6+x^16*y^4+4*x^18*y^2+x^7*y^14+2*x^8*y^13+4*x^9*y^12+3*x^10*y^11+2*x^11*y^10+x^13*y^8+x^14*y^7+3*x^15*y^6+3*x^16*y^5+3*x^17*y^4+x^18*y^3+2*x^8*y^14+3*x^9*y^13+2*x^10*y^12+4*x^13*y^9+2*x^14*y^8+2*x^15*y^7+x^17*y^5+4*x^18*y^4+3*x^19*y^3+x^8*y^15+2*x^10*y^13+2*x^11*y^12+4*x^12*y^11+x^15*y^8+4*x^17*y^6+3*x^18*y^5+4*x^19*y^4+x^9*y^15+2*x^12*y^12+2*x^13*y^11+x^14*y^10+4*x^15*y^9+3*x^16*y^8+3*x^17*y^7+x^18*y^6+x^19*y^5+2*x^10*y^15+3*x^11*y^14+4*x^12*y^13+x^13*y^12+3*x^14*y^11+x^15*y^10+4*x^16*y^9+3*x^19*y^6+4*x^10*y^16+x^11*y^15+x^12*y^14+2*x^13*y^13+2*x^14*y^12+2*x^15*y^11+3*x^16*y^10+3*x^18*y^8+4*x^21*y^5+3*x^22*y^4+3*x^9*y^18+4*x^10*y^17+x^11*y^16+2*x^12*y^15+x^16*y^11+4*x^18*y^9+3*x^19*y^8+x^20*y^7+3*x^21*y^6+2*x^10*y^18+3*x^11*y^17+2*x^13*y^15+x^14*y^14+2*x^15*y^13+3*x^17*y^11+4*x^18*y^10+4*x^20*y^8+x^22*y^6+4*x^10*y^19+3*x^12*y^17+x^13*y^16+2*x^15*y^14+x^16*y^13+3*x^17*y^12+3*x^18*y^11+4*x^19*y^10+3*x^20*y^9+4*x^21*y^8+4*x^22*y^7+x^11*y^19+2*x^12*y^18+3*x^13*y^17+2*x^14*y^16+x^15*y^15+2*x^16*y^14+3*x^18*y^12+2*x^22*y^8+3*x^15*y^16+2*x^16*y^15+x^17*y^14+x^22*y^9+3*x^12*y^20+4*x^13*y^19+x^14*y^18+2*x^15*y^17+x^16*y^16+2*x^17*y^15+x^19*y^13+x^21*y^11+4*x^22*y^10+3*x^24*y^8+3*x^25*y^7+4*x^14*y^19+4*x^15*y^18+x^16*y^17+4*x^17*y^16+4*x^18*y^15+3*x^19*y^14+x^20*y^13+2*x^22*y^11+x^23*y^10+4*x^24*y^9+x^25*y^8+3*x^13*y^21+2*x^14*y^20+2*x^15*y^19+3*x^16*y^18+4*x^17*y^17+4*x^18*y^16+4*x^19*y^15+2*x^20*y^14+4*x^21*y^13+3*x^22*y^12+2*x^23*y^11+2*x^24*y^10+2*x^25*y^9+4*x^14*y^21+4*x^15*y^20+4*x^16*y^19+2*x^17*y^18+2*x^18*y^17+x^21*y^14+2*x^22*y^13+2*x^23*y^12+x^24*y^11+2*x^25*y^10+3*x^16*y^20+4*x^19*y^17+2*x^20*y^16+4*x^21*y^15+4*x^22*y^14+x^23*y^13+x^25*y^11+3*x^16*y^21+2*x^17*y^20+3*x^18*y^19+2*x^19*y^18+2*x^21*y^16+2*x^22*y^15+2*x^23*y^14+x^24*y^13+3*x^25*y^12+x^15*y^23+2*x^16*y^22+3*x^17*y^21+3*x^18*y^20+x^20*y^18+2*x^21*y^17+2*x^22*y^16+3*x^23*y^15+x^25*y^13+3*x^27*y^11+x^15*y^24+2*x^16*y^23+4*x^17*y^22+x^18*y^21+x^19*y^20+4*x^21*y^18+4*x^22*y^17+4*x^23*y^16+x^25*y^14+x^27*y^12+2*x^16*y^24+2*x^17*y^23+3*x^20*y^20+3*x^21*y^19+4*x^22*y^18+2*x^25*y^15+3*x^26*y^14+2*x^27*y^13+3*x^16*y^25+2*x^17*y^24+3*x^18*y^23+3*x^19*y^22+4*x^20*y^21+3*x^21*y^20+3*x^22*y^19+4*x^23*y^18+2*x^24*y^17+2*x^25*y^16+x^26*y^15+2*x^27*y^14+2*x^17*y^25+x^18*y^24+3*x^19*y^23+3*x^20*y^22+3*x^21*y^21+x^22*y^20+2*x^23*y^19+4*x^25*y^17+2*x^26*y^16+x^27*y^15+2*x^18*y^25+4*x^19*y^24+3*x^20*y^23+4*x^21*y^22+2*x^22*y^21+x^23*y^20+3*x^24*y^19+4*x^25*y^18+2*x^26*y^17+3*x^27*y^16+2*x^18*y^26+x^19*y^25+2*x^22*y^22+x^23*y^21+3*x^24*y^20+2*x^25*y^19+2*x^26*y^18+x^27*y^17+4*x^19*y^26+2*x^20*y^25+2*x^22*y^23+x^23*y^22+2*x^24*y^21+x^25*y^20+x^27*y^18+x^19*y^27+4*x^20*y^26+2*x^24*y^22+2*x^25*y^21+x^27*y^19+4*x^20*y^27+3*x^24*y^23+2*x^25*y^22+2*x^25*y^23", domain, "x", "y");

        for (int i = 0; i < 50; ++i) {
            PrivateRandom.getRandom().setSeed(i);
            MultivariatePolynomialZp64 gcd = MultivariateGCD.KaltofenMonaganSparseModularGCDInGF(a, b);
            assertTrue(dividesQ(a, gcd));
            assertTrue(dividesQ(b, gcd));
        }
    }

    @Test // =====> testArrayGCD1   elapsed 1530us
    public void testArrayGCD1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64
                gcd = MultivariatePolynomialZp64.parse("c*a + b + a + c^15*a^3 + b*c*a^5 + d^2*c*a", domain, MonomialOrder.DEFAULT),
                arr[] = {
                        MultivariatePolynomialZp64.parse("c*b*a^2 + b^2 + c + b*a^15 + d", domain, MonomialOrder.DEFAULT).multiply(gcd),
                        MultivariatePolynomialZp64.parse("a^12 + 2*b^12 + 2*c + c*a^5 + d*a", domain, MonomialOrder.DEFAULT).multiply(gcd),
                        MultivariatePolynomialZp64.parse("a^2 + 2*b^12 + 2*c + c*a^5 + d*a", domain, MonomialOrder.DEFAULT).multiply(gcd),
                        MultivariatePolynomialZp64.parse("a^12 - 2*b^2 + 2*c^3 + c*a^5 + d*a", domain, MonomialOrder.DEFAULT).multiply(gcd),
                        MultivariatePolynomialZp64.parse("b^12 - 2*b^2 + 2*c^3 + c*a^5 + d*a", domain, MonomialOrder.DEFAULT).multiply(gcd)
                };

        MultivariatePolynomialZp64 aGcd = MultivariateGCD.PolynomialGCD(arr);
        assertEquals(gcd, aGcd);
    }

    @Test // =====> testArrayGCD2   elapsed 15ms
    @SuppressWarnings("unchecked")
    public void testArrayGCD2() throws Exception {
        // very tricky example with recursive finite fields modulo 2
        FiniteField<UnivariatePolynomialZp64> minorDomain = new FiniteField<>(UnivariatePolynomialZ64.create(1, 0, 1, 1).modulus(2));
        Coder<UnivariatePolynomialZp64, ?, ?> yParser = Coder.mkPolynomialCoder(minorDomain, "y");
        Coder<UnivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> xyParser = Coder.mkUnivariateCoder(Rings.UnivariateRing(minorDomain), yParser, "x");

        FiniteField<UnivariatePolynomial<UnivariatePolynomialZp64>> domain = new FiniteField<>(xyParser.parse("(1+y^2)+(y^2)*x+(y+y^2)*x^2+x^3"));
        Coder<MultivariatePolynomial<UnivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> parser = Coder.mkMultivariateCoder(Rings.MultivariateRing(2, domain), xyParser, "a", "b");
        MultivariatePolynomial<UnivariatePolynomialZp64>
                arr[] = new MultivariatePolynomial[]{
                parser.parse("((y)*x^0)*b+((y+y^2)*x^0)*b^2+b^3+((1+y)*x^0)*b^4"),
                parser.parse("((1+y+y^2)*x^0)*b^4+((1+y)*x^0)*b^8"),
                parser.parse("((y+y^2)*x^0)*b^3+((1+y)*x^0)*b^4+((1+y^2)*x^0)*b^5"),
                parser.parse("((y)*x^0)*b^2+((y+y^2)*x^0)*b^3+((y)*x^0)*b^4+((y^2)*x^0)*b^5+((1+y+y^2)*x^0)*b^6"),
                parser.parse("((y^2)*x^0)+b+((y)*x^0)*b^3+((y+y^2)*x^0)*b^4+((1+y)*x^0)*b^6+((1+y^2)*x^0)*b^7"),
                parser.parse("((y+y^2)*x^0)*b^7"),
                parser.parse("((1+y^2)*x^0)*b^5"),
                parser.parse("((1+y+y^2)*x^0)*b^4+((1+y)*x^0)*b^8"),
                parser.parse("((1+y)*x^0)*b^4+((y^2)*x^0)*b^5+((y+y^2)*x^0)*b^6+((y^2)*x^0)*b^7"),
                parser.parse("((y+y^2)*x^0)*b^3+((1+y)*x^0)*b^4+((1+y^2)*x^0)*b^5"),
                parser.parse("b^2+((y+y^2)*x^0)*b^3+((1+y+y^2)*x^0)*b^5+((y^2)*x^0)*b^6+b^7+((y)*x^0)*b^8"),
                parser.parse("((y)*x^0)*b^2+((y+y^2)*x^0)*b^3+((y)*x^0)*b^4+((y^2)*x^0)*b^5+((1+y+y^2)*x^0)*b^6"),
                parser.parse("((y+y^2)*x^0)*b^7"),
                parser.parse("((y^2)*x^0)+b+((y)*x^0)*b^3+((y+y^2)*x^0)*b^4+((1+y)*x^0)*b^6+((1+y^2)*x^0)*b^7"),
                parser.parse("((1+y^2)*x^0)+((y)*x^0)*b+b^3+((1+y)*x^0)*b^4"),
        };
        for (int i = 0; i < 100; i++)
            assertTrue(MultivariateGCD.PolynomialGCD(arr).isConstant());
    }

    @Test // =====> testArrayGCD3   elapsed 4ms
    @SuppressWarnings("unchecked")
    public void testArrayGCD3() throws Exception {
        // very tricky example with recursive finite fields modulo 2
        FiniteField<UnivariatePolynomialZp64> minorDomain = new FiniteField<>(UnivariatePolynomialZ64.create(1, 0, 1, 1).modulus(2));
        Coder<UnivariatePolynomialZp64, ?, ?> yParser = Coder.mkPolynomialCoder(minorDomain, "y");
        Coder<UnivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> xyParser = Coder.mkUnivariateCoder(Rings.UnivariateRing(minorDomain), yParser, "x");

        FiniteField<UnivariatePolynomial<UnivariatePolynomialZp64>> domain = new FiniteField<>(xyParser.parse("(1+y^2)+(y^2)*x+(y+y^2)*x^2+x^3"));
        Coder<MultivariatePolynomial<UnivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> parser = Coder.mkMultivariateCoder(Rings.MultivariateRing(2, domain), xyParser, "a", "b");
        MultivariatePolynomial<UnivariatePolynomialZp64>
                arr[] = new MultivariatePolynomial[]{
                parser.parse("(x)*b+(x+x^2)*b^2+b^3+(1+x)*b^4"),
                parser.parse("(1+x+x^2)*b^4+(1+x)*b^8"),
                parser.parse("(x+x^2)*b^3+(1+x)*b^4+(1+x^2)*b^5"),
                parser.parse("(x)*b^2+(x+x^2)*b^3+(x)*b^4+(x^2)*b^5+(1+x+x^2)*b^6"),
                parser.parse("(x^2)+b+(x)*b^3+(x+x^2)*b^4+(1+x)*b^6+(1+x^2)*b^7"),
                parser.parse("(x+x^2)*b^7"),
                parser.parse("(1+x^2)*b^5"),
                parser.parse("(1+x+x^2)*b^4+(1+x)*b^8"),
                parser.parse("(1+x)*b^4+(x^2)*b^5+(x+x^2)*b^6+(x^2)*b^7"),
                parser.parse("(x+x^2)*b^3+(1+x)*b^4+(1+x^2)*b^5"),
                parser.parse("b^2+(x+x^2)*b^3+(1+x+x^2)*b^5+(x^2)*b^6+b^7+(x)*b^8"),
                parser.parse("(x)*b^2+(x+x^2)*b^3+(x)*b^4+(x^2)*b^5+(1+x+x^2)*b^6"),
                parser.parse("(x+x^2)*b^7"),
                parser.parse("(x^2)+b+(x)*b^3+(x+x^2)*b^4+(1+x)*b^6+(1+x^2)*b^7"),
                parser.parse("(1+x^2)+(x)*b+b^3+(1+x)*b^4"),
        };
        assertTrue(MultivariateGCD.PolynomialGCD(arr).isConstant());
    }

    @Test // =====> testArrayGCD4   elapsed 2ms
    @SuppressWarnings("unchecked")
    public void testArrayGCD4() throws Exception {
        // very tricky example with recursive finite fields modulo 2
        MultivariatePolynomial<MultivariatePolynomial<BigInteger>>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("2").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("4").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("6").setNVariables(4).asOverMultivariateEliminate(2, 3)};

        assertEquals(pp[0].parsePoly("2"), MultivariateGCD.PolynomialGCD(pp));


        pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("2*x*y*z*t").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("4 - 2").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("6*x*y*z*t - 2").setNVariables(4).asOverMultivariateEliminate(2, 3)};

        assertEquals(pp[0].parsePoly("2"), MultivariateGCD.PolynomialGCD(pp));

        pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("2*x*y*z*t - 16").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("4*x - 2").setNVariables(4).asOverMultivariateEliminate(2, 3),
                MultivariatePolynomial.parse("6*x*y*z*t - 2").setNVariables(4).asOverMultivariateEliminate(2, 3)};

        assertEquals(pp[0].parsePoly("2"), MultivariateGCD.PolynomialGCD(pp));
    }

    @Test // =====> testCommonZeroes1   elapsed 1048us
    public void testCommonZeroes1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1 + c*b*a^2+b^2 + c + a^5", domain, MonomialOrder.DEFAULT),
                b = MultivariatePolynomialZp64.parse("a^2+2*b^2 + 2*c + a^5", domain, MonomialOrder.DEFAULT);
        ZeroVariables pZeros = commonPossibleZeroes(a, b, a.nVariables);
        assertTrue(pZeros.pZeros.size() > 0);
        for (BitSet z : pZeros.pZeros) {
            assertFalse(setZeroes(a, z).isZero());
            assertFalse(setZeroes(b, z).isZero());
        }
    }

    @Test // =====> testCommonZeroes2   elapsed 512us
    public void testCommonZeroes2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1 + c*b*a^2+b^2 + c + a^5*d + e*a + b*g - f", domain, MonomialOrder.DEFAULT),
                b = MultivariatePolynomialZp64.parse("a^2+2*b^2 + 2*c + a^5 + a*b*c*d*e*f + g", domain, MonomialOrder.DEFAULT);
        MultivariatePolynomialZp64 ac = a.clone();
        ZeroVariables pZeros = commonPossibleZeroes(a, b, a.nVariables);
        assertTrue(pZeros.pZeros.size() > 0);
        for (BitSet z : pZeros.pZeros) {
            assertFalse(setZeroes(a, z).isZero());
            assertFalse(setZeroes(b, z).isZero());
        }
    }

    private static MultivariatePolynomialZp64 setZeroes(MultivariatePolynomialZp64 poly, BitSet zeroes) {
        TIntArrayList vars = new TIntArrayList();
        for (int i = 0; i < zeroes.size(); i++)
            if (zeroes.get(i))
                vars.add(i);

        return poly.evaluateAtZero(vars.toArray());
    }

    @Test // =====> testEZEvaluations1   elapsed 3s
    public void testEZEvaluations1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("4*b^5*d^3*e^6+2*a^4*c^2*d^3*e^2+9939223*a^4*b^3*c^3*e+7*a^4*b^3*c^5*d^5*e+3*a^4*b^5*c^5*d^5*e^4+2*a^5*c^3*d^6*e^3+9939225*a^5*c^5*d*e^3+3*a^5*b*c^6*d^4*e^3+7*a^5*b^2*c^3*d^4*e^4+9939223*a^5*b^6*c^4*d", domain),
                b = MultivariatePolynomialZp64.parse("4*b^6*c^4+6*a*b*c^2*e+9939226*a*b^5*c*d^3*e^3+8*a^2*b^4*c^5*d^5*e^3+3*a^3*c^6*d^4*e^2+9939223*a^3*b^2*c*d^2*e^2+2*a^4*b^3*c^6*d^2*e+2*a^5*b^3*d^6*e^2+4*a^5*b^5*c^5*d*e^4+9*a^6*b*c^5*d^2+9939221*a^6*b^2*c^4*d^3*e^5", domain);
        EZGCDEvaluations evals = new EZGCDEvaluations(a, b, a.nVariables, getRandom());
        for (int j = 0; j < 10; j++) {
            System.out.println(j);
            evals.nextEvaluation();
            assertEquals(a, evals.reconstruct(evals.aReduced));
            assertEquals(b, evals.reconstruct(evals.bReduced));
        }
    }

    @Test // =====> testEZEvaluations2   elapsed 52ms
    public void testEZEvaluations2() throws Exception {
        RandomGenerator rnd = getRandom();
        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 4, 2, 3, 5, 7, rnd);
        rnd.setSeed(42);
        GCDSample<MonomialZp64, MultivariatePolynomialZp64> sample = sampleData.nextSample(false, true);
        for (MultivariatePolynomialZp64[] pp : new MultivariatePolynomialZp64[][]{{sample.a, sample.b}, {sample.aCoFactor, sample.bCoFactor}}) {
            MultivariatePolynomialZp64
                    a = pp[0],
                    b = pp[1];
            if (a.isConstant() || b.isConstant())
                continue;
            EZGCDEvaluations evals = new EZGCDEvaluations(a, b, a.nVariables, rnd);
            for (int j = 0; j < 10; j++) {
                evals.nextEvaluation();
                assertEquals(a, evals.reconstruct(evals.aReduced));
                assertEquals(b, evals.reconstruct(evals.bReduced));
            }
        }
    }

    @Test // =====> testEZEvaluationsRandom   elapsed 1643ms
    public void testEZEvaluationsRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        int nIterations = its(10, 100);
        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 4, 2, 3, 5, 7, rnd);
        for (int i = 0; i < nIterations; i++) {
            System.out.println(i);
            GCDSample<MonomialZp64, MultivariatePolynomialZp64> sample = sampleData.nextSample(false, true);
            for (MultivariatePolynomialZp64[] pp : new MultivariatePolynomialZp64[][]{{sample.a, sample.b}, {sample.aCoFactor, sample.bCoFactor}}) {
                MultivariatePolynomialZp64
                        a = pp[0],
                        b = pp[1];
                if (a.isConstant() || b.isConstant())
                    continue;
                EZGCDEvaluations evals = new EZGCDEvaluations(a, b, a.nVariables, rnd);
                for (int j = 0; j < 10; j++) {
                    evals.nextEvaluation();
                    assertEquals(a, evals.reconstruct(evals.aReduced));
                    assertEquals(b, evals.reconstruct(evals.bReduced));
                }
            }
        }
    }

    @Test // =====> testEZGCD1   elapsed 5ms
    public void testEZGCD1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("c*b*a^2 + b^2 + c + b*a^15 + d", domain, MonomialOrder.DEFAULT),
                b = MultivariatePolynomialZp64.parse("a^12 + 2*b^12 + 2*c + c*a^5 + d*a", domain, MonomialOrder.DEFAULT),
                gcd = MultivariatePolynomialZp64.parse("c*a + b + a + c^15*a^3 + b*c*a^5 + d^2*c*a", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);

        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD2   elapsed 160ms
    public void testEZGCD2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64

                u = MultivariatePolynomialZp64.parse("c*b*a + b^2 + c + b*a^2 + 1", domain, MonomialOrder.DEFAULT),
                v = MultivariatePolynomialZp64.parse("2 + a^2 + 2*b^2 + 2*c + c*a^2 + a", domain, MonomialOrder.DEFAULT),
                a = u.clone().square().multiply(u).multiply(v),
                b = v.clone().square().multiply(u),
                gcd = MultivariatePolynomialZp64.parse("c*a + b + a + c*a^3 + b*c*a^2 + c*a", domain, MonomialOrder.DEFAULT)
                        .multiply(u).multiply(v).multiply(v);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD3   elapsed 15ms
    public void testEZGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        IntegersZp64 domain = new IntegersZp64(15627449);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("3*b^2*c*d+5*a*d^2+15627444*a*b^2*c+15627440*a^2*b^2*c^2", domain, MonomialOrder.DEFAULT),
                b = MultivariatePolynomialZp64.parse("4*b*c^2*d^2+8*b^2*c*d^2+15627440*b^3*d^2+3*a*c^3+a^2*b", domain, MonomialOrder.DEFAULT),
                gcd = MultivariatePolynomialZp64.parse("15627440*b^2*c^3+7*b^3*c+6*a*b^2*c^3*d^2+7*a^2*b*d+5*a^3*c^2*d^3+6*a^3*b^2*c^2*d", domain, MonomialOrder.DEFAULT);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD4   elapsed 803us
    public void testEZGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        IntegersZp64 domain = new IntegersZp64(15627449);
        String[] vars = {"b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("13312440+25*d", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("902776+10710357*c+6547542*c^2+4925527*b+2965659*b*c+20*b*c^2+2903103*b^2+40*b^2*c+15627404*b^3", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("102871+8266210*d+5121205*d^2+16248*d^3+1722152*c+2574791*c*d+10788581*c*d^2+15247596*c*d^3+8472569*c^2+898837*c^2*d+14099452*c^2*d^2+c^2*d^3", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD5   elapsed 2ms
    public void testEZGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(0);
        IntegersZp64 domain = new IntegersZp64(24254707);
        String[] vars = {"b", "c", "d", "a"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("24254706*b+b^2*c+7*a*c+5*a*b^2+2*a^2*b*c+24254705*a^2*b*c^2+a^2*b^2*c^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("4*b*c+6*b^2*c+4*a^2+7*a^2*b*c^2+3*a^2*b^2*c", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("9*c+24254705*a*c^2+6*a*b^2+4*a^2+3*a^2*c+24254698*a^2*b^2*c", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD6   elapsed 868us
    public void testEZGCD6() throws Exception {
        PrivateRandom.getRandom().setSeed(0);
        IntegersZp64 domain = new IntegersZp64(24254707);
        String[] vars = {"b", "c", "d", "a"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("11106134 + 20915017*a + 948048*a^2 + 18101438*b + 8523620*a*b + 19589342*a^2*b + 13684993*b^2 + 8219998*a*b^2 + 24254698*a^2*b^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("12760358 + 5425698*a + 5129306*a^2 + 14380683*b + 16257092*a*b + 24254680*a^2*b + 4515223*b^2 + 24254644*a*b^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("9740181 + 21604578*a + 9691573*a^2 + 11040951*b + 17951441*a*b + a^2*b", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD7   elapsed 4ms
    public void testEZGCD7() throws Exception {
        PrivateRandom.getRandom().setSeed(6);
        IntegersZp64 domain = new IntegersZp64(28565767);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6*c^2*d+3*a*d^2+9*a*d^3+28565764*a^2*b*d^2+2*a^2*b^2*c*d^2+9*a^3*c*d^2+3*a^3*b^2*c*d^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("c^2*d^2+28565766*a*b^2*d^2+2*a*b^2*c*d+28565763*a^2*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("6*c*d+7*c*d^2+7*a*b^2+8*a^2*b*c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD8   elapsed 282us
    public void testEZGCD8() throws Exception {
        PrivateRandom.getRandom().setSeed(6);
        IntegersZp64 domain = new IntegersZp64(28565767);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1564663 + 63*d", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("15383307 + 22520298*b + 12522045*b^2 + 9677552*c + 7*c^2 + 3221049*d + 5785123*b*d + 28565760*b^2*d", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("21095373 + c", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD9   elapsed 19ms
    public void testEZGCD9() throws Exception {
        PrivateRandom.getRandom().setSeed(27);
        IntegersZp64 domain = new IntegersZp64(8678027);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6*c^2*d+2*b^2*c^2*d^2+a*c^2+a*b^2*c^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("3*b^3*c*d^2+2*a*b^2*d^2+8678026*a^2*c^3*d^3+8678018*a^2*b+5*a^3*b^2*c^3*d", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("8678020*c^2+3*b*c+6*b^2*c^2*d+2*a^2*b^3*c+8678020*a^3*c*d^3", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD10   elapsed 3ms
    public void testEZGCD10() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        IntegersZp64 domain = new IntegersZp64(28996511);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6*b^2*c^2+b^3*c+28996505*a*b*c+28996507*a*b^2+4*a*b^2*c^3", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("6*b^2*c+12*a*b^2*c+3*a^2*b^2*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("9+28996505*c^2+7*a*b*c+8*a*b*c^2+9*a^2*b^2*c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD11   elapsed 20ms
    public void testEZGCD11() throws Exception {
        PrivateRandom.getRandom().setSeed(11);
        IntegersZp64 domain = new IntegersZp64(22687397);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("7*d+7*b^2*c^2+3*a*c^2+2*a^2*d+6*a^2*d^2+9*a^2*c^2*d", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("2*a*b*c^2+22687392*a*b^2*c^2*d^2+9*a*b^2*c^3+22687395*a*b^3*c*d+3*a^2*c*d^3+6*a^3+8*a^3*b^2*c*d^3", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("22687391*d^2+22687393*b*c^2*d^2+5*a+5*a*b+8*a*b^2*c^2*d^2+4*a^2*b^2*c^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD12   elapsed 791us
    public void testEZGCD12() throws Exception {
        PrivateRandom.getRandom().setSeed(11);
        IntegersZp64 domain = new IntegersZp64(22687397);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("11244809 + 22687361*a^2 + 30*b + 11244809*c + 30*b*c", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("21851038 + 15893672*a + 14760820*a^2 + 12564491*a^3 + 11694181*a^4 + 16683489*a^5 + 11244809*b + 9237198*a*b + 4795625*a^2*b + 6123527*a^3*b + 5432952*a^4*b + 22659623*a^5*b + 15*b^2 + 14306815*a*b^2 + 15896207*a^2*b^2 + 9051*a^3*b^2 + 22645277*a^4*b^2 + 11745*a*b^3 + 35100*a^2*b^3 + 2299511*c + 15893672*a*c + 1652686*a^2*c + 6762669*a^3*c + 19488092*a^4*c + 3757801*a^5*c + 11247299*b*c + 9237198*a*b*c + 12416836*a^2*b*c + 19792351*a^3*b*c + 9233182*a^4*b*c + 22686257*a^5*b*c + 15*b^2*c + 14306815*a*b^2*c + 15890231*a^2*b^2*c + 22668450*a^3*b^2*c + 22687373*a^4*b^2*c + 11745*a*b^3*c + 35100*a^2*b^3*c + 9448492*c^2 + 11808378*a*c^2 + 20057596*a^2*c^2 + 20794627*a^3*c^2 + 18785592*a^4*c^2 + 4960279*a^5*c^2 + 21841949*b*c^2 + 2669531*a*b*c^2 + 13194739*a^2*b*c^2 + 10436283*a^3*b*c^2 + 19195830*a^4*b*c^2 + 5885213*a^5*b*c^2 + 8273012*b^2*c^2 + 20876792*a*b^2*c^2 + 19005620*a^2*b^2*c^2 + 6957603*a^3*b^2*c^2 + 17675248*a^4*b^2*c^2 + 22615613*a^5*b^2*c^2 + 7430983*b^3*c^2 + 18240653*a*b^3*c^2 + 4232985*a^2*b^3*c^2 + 147044*a^3*b^3*c^2 + 48*a^4*b^3*c^2 + 5976*b^4*c^2 + 42092*a*b^4*c^2 + 24*a^2*b^4*c^2 + 9149236*c^3 + 9268605*a*c^3 + 4080916*a^2*c^3 + 6163313*a^3*c^3 + 2920*a^4*c^3 + 3127846*a^5*c^3 + 11966525*b*c^3 + 22671967*a*b*c^3 + 22683508*a^2*b*c^3 + 21437536*a^3*b*c^3 + 4916674*a^5*b*c^3 + 1168*b^2*c^3 + 2213725*a^3*b^2*c^3 + 22684357*a^5*b^2*c^3 + 61720*a^3*b^3*c^3 + 12936437*c^4 + 9268605*a*c^4 + 15448872*a^2*c^4 + 21785946*a^3*c^4 + 4125558*a^4*c^4 + 19072832*a^5*c^4 + 19011898*b*c^4 + 22671967*a*b*c^4 + 4208806*a^2*b*c^4 + 15458267*a^3*b*c^4 + 22681557*a^4*b*c^4 + 3920321*a^5*b*c^4 + 41624*b^2*c^4 + 22684477*a^2*b^2*c^4 + 3920321*a^3*b^2*c^4 + 7937375*a^5*b^2*c^4 + 12854049*a^3*b^3*c^4 + 6080*a^5*b^3*c^4 + 3040*a^3*b^4*c^4 + 12125660*a*c^5 + 12572867*a^3*c^5 + 12572867*a*b*c^5 + 22685877*a^3*b*c^5 + 22686637*a*b^2*c^5", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("1", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD13   elapsed 3ms
    public void testEZGCD13() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(13666309);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("4*b^2*c^3+3*a*b^3+2*a*b^3*c+4*a*b^3*c^3+2*a^3*b^2*c", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("5*b*c+5*a*b+13666307*a*b*c+7*a^2*c+6*a^2*c^2+13666307*a^2*b", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("4*b^2+15*a*b^2*c+4*a^2+13666308*a^2*c", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD14   elapsed 544us
    public void testEZGCD14() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(13666309);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("3338057+430735*b+13248829*b^2+11374034*b^3+3423812*a+698808*a*b+7995810*a*b^2+60*a*b^3+8933188*a^2+13666305*a^2*b", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("7990553+3359122*b+846494*b^2+131346*a+12831229*a*b+6789484*a*b^2+12272859*a^2+7995810*a^2*b+90*a^2*b^2+13666303*a^3", domain, MonomialOrder.DEFAULT, vars);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }


    @Test // =====> testEZGCD15   elapsed 405us
    public void testEZGCD15() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(22687397);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("30*b+30*b*c+22687361*a^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("2490*b*c+6800*b*c^2+4310*b*c^3+15*b^2+15*b^2*c+1168*b^2*c^3+41624*b^2*c^4+14028*b^3*c^2+5976*b^4*c^2+40950*a*b+40950*a*b*c+22671967*a*b*c^3+22671967*a*b*c^4+28*a*b^2*c^2+22686637*a*b^2*c^5+11745*a*b^3+11745*a*b^3*c+42092*a*b^4*c^2+22684409*a^2*c+22681057*a^2*c^2+22645773*a^2*c^3+22687379*a^2*b+22673369*a^2*b*c+22681172*a^2*b*c^2+22683508*a^2*b*c^3+83248*a^2*b*c^4+22681421*a^2*b^2*c+28056*a^2*b^2*c^2+22684477*a^2*b^2*c^4+35100*a^2*b^3+35100*a^2*b^3*c+11952*a^2*b^3*c^2+24*a^2*b^4*c^2+22638257*a^3+22687369*a^3*c+18516*a^3*c^3+760*a^3*c^4+56*a^3*b*c^2+22685877*a^3*b*c^5+9051*a^3*b^2+22668450*a^3*b^2*c+147044*a^3*b^3*c^2+61720*a^3*b^3*c^3+3040*a^3*b^4*c^4+7470*a^4*c^2+2920*a^4*c^3+22681557*a^4*b*c^4+22645277*a^4*b^2+22687373*a^4*b^2*c+48*a^4*b^3*c^2+22659623*a^5*b+22686257*a^5*b*c+22615613*a^5*b^2*c^2+22684357*a^5*b^2*c^3+6080*a^5*b^3*c^4", domain, MonomialOrder.DEFAULT, vars);
//                gcd = MultivariatePolynomialZp64.parse("4*b^2+15*a*b^2*c+4*a^2+13666308*a^2*c", ring, MonomialOrder.DEFAULT, vars);
//        a = a.multiply(gcd);
//        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD16   elapsed 6ms
    public void testEZGCD16() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(9607987);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("8+3*a+6*a^2*b*c^2+6*a^2*b^2*c", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("6*b^2+8*b^2*c+9607986*a*b+a*b*c^2+a^2*c+2*a^2*b+7*a^2*b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("6*b*c+6*a*b*c+9607979*a*b^2*c^2+3*a*b^3*c^3+9607982*a^3*b^3", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD17   elapsed 856us
    public void testEZGCD17() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(9607987);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6998457+9068733*c+9042619*c^2+1283895*c^3+3738482*c^4+4888116*c^5+5574926*b+6161435*b*c+3428490*b*c^2+1423636*b*c^3+7718978*b*c^4+9607957*b*c^5+5803370*b^2+3797801*b^2*c+1899022*b^2*c^2+1286548*b^2*c^3+6635895*b^3+1871925*b^3*c+2687295*b^3*c^2+9406093*b^3*c^3+3101467*b^4+8800611*b^4*c+5291877*b^4*c^2+18*b^4*c^3", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("4572338+8988826*c+5700575*c^2+6517458*c^3+7602268*b+8674751*b*c+4253585*b*c^2+9607947*b*c^3+1248945*b^2+3165308*b^2*c+9363637*b^3+9338813*b^3*c+8757270*b^4+24*b^4*c", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("1", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD18   elapsed 6ms
    public void testEZGCD18() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(592346501);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("56*b+47*b*c+37*a*b*c^2+43*a^2*b*c", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("26*b*c^2+54*a*c+4*a^2*b+2*a^2*b*c+42*a^2*b*c^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("13*c+50*a*b+51*a^2*b*c+65*a^2*b^2*c^3+33*a^3*b^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEZGCD19   elapsed 888us
    public void testEZGCD19() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(592346501);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("125427093 + 287292359*a + 259899124*a^2 + 224214583*b + 423992120*a*b + 1221*a^2*b + 522309957*b^2 + 1419*a*b^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("364341910 + 56968290*a + 134477777*a^2 + 264733241*b + 223672725*a*b + 365910146*a^2*b + 448183856*b^2 + 56041492*a*b^2 + 1386*a^2*b^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("8864159 + 332216825*a + 307171438*a^2 + 574396609*a^3 + b", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD_random1   elapsed 14s
    public void testEEZGCD_random1() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 4, minDegree = 2, maxDegree = 3, minSize = 5, maxSize = 7;
        int nIterations = its(1000, 1000);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EZ-GCD", MultivariateGCD::EZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test // =====> testEEZGCD_sparse_variables_random   elapsed 9s
    public void testEEZGCD_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<MonomialZp64, MultivariatePolynomialZp64>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EZ-GCD", MultivariateGCD::EZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test // =====> testEEZGCD1   elapsed 2ms
    public void testEEZGCD1() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(9607987);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6998457+9068733*c+9042619*c^2+1283895*c^3+3738482*c^4+4888116*c^5+5574926*b+6161435*b*c+3428490*b*c^2+1423636*b*c^3+7718978*b*c^4+9607957*b*c^5+5803370*b^2+3797801*b^2*c+1899022*b^2*c^2+1286548*b^2*c^3+6635895*b^3+1871925*b^3*c+2687295*b^3*c^2+9406093*b^3*c^3+3101467*b^4+8800611*b^4*c+5291877*b^4*c^2+18*b^4*c^3", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("4572338+8988826*c+5700575*c^2+6517458*c^3+7602268*b+8674751*b*c+4253585*b*c^2+9607947*b*c^3+1248945*b^2+3165308*b^2*c+9363637*b^3+9338813*b^3*c+8757270*b^4+24*b^4*c", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("1", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD2   elapsed 6ms
    public void testEEZGCD2() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(16604789);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("6*a*c*d^2+2*a*c^2+2*a*b^2*d+3*a*b^2*c^2*d^2", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("7*b^2*c^3*d+6*a^2*b^2*c^3+16604785*a^2*b^3*d+9*a^3*b*c^2+16604781*a^3*b*c^3", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("4*c^2*d+8*a+4*a*c^3*d^2+16604785*a^2*d^3+2*a^2*b^2*c^3*d", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD3   elapsed 608us
    public void testEEZGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(31012727);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("4*c+3*b*c^2+4*a*c+31012723*a*b+2*a^2*c+7*a^2*b*c", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("7*c^2+2*b*c^2+3*a*b*c^2+31012724*a*b^2+a^2*b^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("31012726*b*c+31012724*b*c^2+3*a*b", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD4   elapsed 490us
    public void testEEZGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(31012727);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("31012718 + 3*b", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("12*a + 47454*a^2 + 30989135*b + 12*a*b + 9801*a^2*b + 41292*a*b^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("20675151*a + 31012726*a^2 + b", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD5   elapsed 22ms
    public void testEEZGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(78);
        IntegersZp64 domain = new IntegersZp64(26662411);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("26662407*a^3*b^3*c*d^5*e^6+3*a^3*b^6*c^3*e^2+4*a^4*b^2*c^3*e^3+8*a^4*b^5*c^6*d^4*e^4+7*a^6*b*c^3*d", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("8*c+26662402*b*c*d+8*b^2*c*d*e+b^2*c^2+a*b*e^2+a^2*d*e^2+8*a^2*d^2*e^2+5*a^2*b*c*d^2*e^2+15*a^2*b*c^2*d", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("5*b^2*d^2*e^2+9*b^2*c^2*e+26662408*b^2*c^2*e^2+3*a*c^2*e+3*a*b*d^2+26662408*a*b^2*d*e+26662407*a*b^2*d^2*e^2+a^2*d^2*e^2+26662402*a^2*b^2*c*e^2", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }


    @Test // =====> testEEZGCD6   elapsed 1235us
    public void testEEZGCD6() throws Exception {
        PrivateRandom.getRandom().setSeed(78);
        IntegersZp64 domain = new IntegersZp64(26662411);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("15*c^2*d+27*b^2+26662402*b^2*d+26662402*a*c+26662399*a*c^2*d+26662384*a^2*b*d", domain, MonomialOrder.DEFAULT, vars),
                b = MultivariatePolynomialZp64.parse("26657131*c^7*d^8+46200*b*c^2*d^2+26593651*b*c^3*d^2+20960*b*c^3*d^3+2620*b^2*c^2*d^2+26652907*b^2*c^5*d^7+3168*b^2*c^5*d^8+83160*b^3*d+26634691*b^3*d^2+26538643*b^3*c*d+78984*b^3*c*d^2+26649835*b^3*c*d^3+4716*b^4*d+26660839*b^4*d^2+7640*a*c^2*d^4+3168*a*c^6*d^7+26651827*a*c^7*d^6+4224*a*c^7*d^8+26634691*a*b*c*d+35952*a*b*c^2+41256*a*b*c^2*d+26612875*a*b*c^2*d^2+26631226*a*b*c^3+36672*a*b*c^3*d+55008*a*b*c^3*d^2+26645643*a*b*c^3*d^3+13752*a*b^2*d^3+26657827*a*b^2*d^4+26660839*a*b^2*c*d+4584*a*b^2*c^2+26660324*a*b^2*c^2*d^2+17640*a*b^2*c^2*d^5+26644891*a*b^2*c^5*d^7+35976*a*b^3*d+26621965*a*b^3*c*d+27720*a*b^3*c*d^2+3465*a*b^4*d+13833*a*b^4*d^3+31752*a*b^4*d^4+26651827*a*b^4*d^5+40*a*b^5*c^6*d^6+72*a*b^7*c^4*d^5+26662387*a*b^7*c^4*d^6+26657827*a^2*c*d^3+3465*a^2*c^2*d^2+26656299*a^2*c^2*d^4+5775*a^2*c^3*d^4+46200*a^2*c^4*d^4+26656571*a^2*c^7*d^8+11992*a^2*b*c^2*d^2+26648929*a^2*b*c^3*d^2+9240*a^2*b*c^3*d^3+38200*a^2*b*c^4*d^4+9504*a^2*b*c^5*d^8+26579251*a^2*b^2*d^2+4494*a^2*b^2*d^3+123768*a^2*b^2*c*d^2+26635078*a^2*b^2*c*d^3+26648362*a^2*b^2*c*d^4+1155*a^2*b^2*c^2*d^2+100680*a^2*b^2*c^2*d^3+26639302*a^2*b^2*c^2*d^4+26648299*a^2*b^2*c^2*d^5+114600*a^2*b^2*c^3*d^2+26657695*a^2*b^3*d^2+68760*a^2*b^3*c^2*d^3+26639491*a^2*b^3*c^2*d^4+3300*a^2*b^4*d^4+206280*a^2*b^4*c*d+26593651*a^2*b^4*c*d^2+26662387*a^2*b^5*c^5*d^5+36888*a^2*b^5*c^6*d^4+26662379*a^2*b^5*c^6*d^6+6336*a^2*b^7*c^4*d^5+26658946*a^3*c^2*d^3+1498*a^3*c^2*d^4+4494*a^3*c^3*d^2+26634691*a^3*c^3*d^3+26657791*a^3*c^3*d^4+35952*a^3*c^4*d^2+26625451*a^3*c^4*d^4+26648659*a^3*b*d^4+26639491*a^3*b*c^3*d^3+17325*a^3*b*c^4*d^2+26631851*a^3*b*c^4*d^4+4497*a^3*b^2*c*d^3+26593651*a^3*b^2*c^2*d+35976*a^3*b^2*c^2*d^3+1100*a^3*b^2*c^2*d^5+51975*a^3*b^2*c^3+26621831*a^3*b^2*c^3*d^2+26630659*a^3*b^3*d^5+22470*a^3*b^3*c^2*d^3+159390*a^3*b^4*c*d+26631751*a^3*b^4*c*d^2+2112*a^3*b^5*c^6*d^6+26662339*a^3*b^6*c^4*d^6+1499*a^4*c^3*d^4+11992*a^4*c^4*d^4+26652016*a^4*b*c*d^4+26579251*a^4*b*c^2*d^4+7490*a^4*b*c^4*d^4+26631751*a^4*b^2*c^2*d+26593651*a^4*b^2*c^2*d^4+5775*a^4*b^2*c^3+26644001*a^4*b^2*c^3*d^2+26456131*a^4*b^3*c*d^2+15456*a^4*b^4*c*d+5152*a^5*b^2*c^3*d^2+26570431*a^5*b^3*c*d^2", domain, MonomialOrder.DEFAULT, vars),
                gcd = MultivariatePolynomialZp64.parse("1", domain, MonomialOrder.DEFAULT, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testEEZGCD_random2   elapsed 6s
    public void testEEZGCD_random2() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 2, maxDegree = 5, minSize = 5, maxSize = 17;
        int nIterations = its(100, 200);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test // =====> testEEZGCD_random3   elapsed 25s
    public void testEEZGCD_random3() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 2, maxDegree = 5, minSize = 5, maxSize = 17;
        int nIterations = its(100, 200);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test // =====> testEEZGCD7   elapsed 132ms
    public void testEEZGCD7() throws Exception {
        IntegersZp64 domain = new IntegersZp64(BigPrimes.nextPrime(1321323));
        MultivariatePolynomialZp64

                u = MultivariatePolynomialZp64.parse("c*b*a + b^2 + c + b*a^2 + 1", domain, MonomialOrder.DEFAULT),
                v = MultivariatePolynomialZp64.parse("2 + a^2 + 2*b^2 + 2*c + c*a^2 + a", domain, MonomialOrder.DEFAULT),
                a = u.clone().square().multiply(u).multiply(v),
                b = v.clone().square().multiply(u),
                gcd = MultivariatePolynomialZp64.parse("c*a + b + a + c*a^3 + b*c*a^2 + c*a", domain, MonomialOrder.DEFAULT)
                        .multiply(u).multiply(v).multiply(v);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test // =====> testPolynomialGCD1   elapsed 2ms
    public void testPolynomialGCD1() throws Exception {
        // nUsedVariables == 1
        IntegersZp64 domain = new IntegersZp64(19);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64 a = MultivariatePolynomialZp64.parse("2 + 3*b + b^2", domain, vars);
        MultivariatePolynomialZp64 b = MultivariatePolynomialZp64.parse("2 + b + a^2 + a^2*b", domain, vars);
        for (int i = 0; i < 1000; i++)
            assertTrue(PolynomialGCD(a, b).isConstant());
    }

    @Test // =====> testPolynomialGCD2   elapsed 3ms
    public void testPolynomialGCD2() throws Exception {
        // nUsedVariables == 1
        IntegersZp64 domain = new IntegersZp64(19);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("a^2 + b*a", domain, vars),
                b = MultivariatePolynomialZp64.parse("b + 1", domain, vars),
                gcd = MultivariatePolynomialZp64.parse("2 + b", domain, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        for (int i = 0; i < 1000; i++)
            assertEquals(gcd, PolynomialGCD(a, b));
    }

    @Test(timeout = 1000) // =====> testPolynomialGCD3   elapsed 633us
    public void testPolynomialGCD3() throws Exception {
        // nUsedVariables == 1
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b"};
        // very tricky example
        // for both a = 0 and a = 1, the gcd is (1 + b^2), while the true gcd is only 1 + b
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("a^2*b + a*b + b + 1", domain, vars),
                b = MultivariatePolynomialZp64.parse("b + 1", domain, vars),
                gcd = MultivariatePolynomialZp64.parse("b + 1", domain, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        for (int i = 0; i < 10; i++)
            assertEquals(gcd, PolynomialGCD(a, b));
    }

    @Test(timeout = 10_000) // =====> testPolynomialGCD4   elapsed 871ms
    public void testPolynomialGCD4() throws Exception {
        Ring<BigInteger> ring = Rings.Z;
        MultivariatePolynomial<BigInteger>
                p1 = polyPow(MultivariatePolynomial.parse("1 + 3*a*b + 5*b*c + 7*c*d + 9*d*e + 11*e*f + 13*f*g + 15*g*a", ring), 3),
                p2 = polyPow(MultivariatePolynomial.parse("1 + 3*a*c + 5*b*d + 7*c*e + 9*f*e + 11*g*f + 13*f*a + 15*g*b", ring), 3),
                p3 = polyPow(MultivariatePolynomial.parse("1 + 3*a*d + 5*b*e + 7*c*f + 9*f*g + 11*g*a + 13*f*b + 15*g*c", ring), 3),
                poly = p1.multiply(p2, p3);
        poly.decrement();
        Assert.assertTrue(poly.asUnivariate(1).content().isConstant());
    }

    @Test(timeout = 1000000) // =====> testSmallDomain4   elapsed 17s
    public void testSmallDomain4() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};

        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+2*c^3*d^2+2*b^3*c^3*d^3*e+a*c^3*d*e+2*a^2*b^3*c^2*d^2*e^3+a^2*b^3*c^3*e^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c^2*d^3*e^3+a*c^3*d*e^2+2*a^3*e^3+2*a^3*b^3*d*e^3+2*a^3*b^3*c*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*a*b^3*c+a^2*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*b^3*c^3*d^3*e+2*a*b^2*c*d^2*e^3+a*b^3*c^2*d*e^2+a^3*b^2*c^3*d^2", domain, vars),
        }, base = arr[0].createOne().multiply(arr);

        MultivariatePolynomialZp64 a = MultivariateFactorization.orderByDegrees(base, false, -1).ordered;
        MultivariatePolynomialZp64 b = a.derivative(0);

        for (int i = 1; i < 5; i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            Assert.assertFalse(PolynomialGCD(a, b).isConstant());
            System.out.println(nanosecondsToString(System.nanoTime() - start));
        }
        //12s
        //3s
        //25s
        //18s
    }

    @Test(timeout = 100000) // =====> testSmallDomain5   elapsed 186ms
    public void testSmallDomain5() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"x", "y"};

        MultivariatePolynomialZp64 a = MultivariatePolynomialZp64.parse("y^2+y^3+x*y^2+x*y^3+2*x^3*y^6+2*x^4*y^2+2*x^4*y^3+x^4*y^5+2*x^5*y^2+x^5*y^3+x^5*y^9+2*x^6+2*x^6*y^3+2*x^6*y^5+2*x^6*y^8+2*x^6*y^9+2*x^7+2*x^7*y^5+x^7*y^6+x^7*y^9+2*x^8*y^5+x^8*y^6+2*x^8*y^12+x^9*y^3+x^9*y^8+2*x^9*y^9+2*x^9*y^12+x^10*y^8+x^10*y^9+x^11*y^5+2*x^11*y^6+x^11*y^9+2*x^11*y^11+x^12*y^5+2*x^12*y^6+2*x^12*y^9+2*x^12*y^11+x^12*y^12+x^13*y^12+2*x^14*y^8+x^14*y^9+x^14*y^12+x^14*y^14+x^16*y^11+x^17*y^11+2*x^19*y^14", domain, vars);
        MultivariatePolynomialZp64 b = MultivariatePolynomialZp64.parse("1+y^4+2*y^5+2*y^9+x+x*y^4+2*x*y^5+2*x*y^9+x^3*y^8+x^3*y^12+2*x^4+x^4*y^3+2*x^4*y^4+x^4*y^5+x^4*y^7+x^4*y^9+2*x^5+2*x^5*y^4+2*x^5*y^5+x^5*y^6+2*x^5*y^9+2*x^5*y^11+x^6*y^2+x^6*y^5+x^6*y^6+2*x^6*y^7+x^6*y^9+2*x^6*y^10+x^6*y^11+x^7*y^2+x^7*y^6+2*x^7*y^7+2*x^7*y^8+2*x^7*y^11+2*x^7*y^12+2*x^8*y^3+2*x^8*y^7+2*x^8*y^8+2*x^8*y^12+x^8*y^14+2*x^9*y^5+2*x^9*y^6+x^9*y^10+x^9*y^11+x^9*y^14+x^10*y^10+2*x^10*y^11+x^11*y^7+x^11*y^8+2*x^11*y^11+2*x^11*y^12+2*x^11*y^13+x^12*y^7+x^12*y^8+x^12*y^11+2*x^12*y^13+2*x^12*y^14+2*x^13*y^9+2*x^13*y^14+2*x^14*y^10+2*x^14*y^11+2*x^14*y^14+x^14*y^16+x^15*y^12+x^16*y^13+x^17*y^13+2*x^19*y^16", domain, vars);

        for (int i = 5; i < 100; i++) {
            PrivateRandom.getRandom().setSeed(i);
            assertFalse(PolynomialGCD(a, b).isConstant());
        }
    }

    @Test
    @SuppressWarnings("unchecked") // =====> testNestedRing1   elapsed 2s
    public void testNestedRing1() throws Exception {
        for (Ring<?> inner1 : Arrays.asList(
                Rings.GF(7, 3),
                Rings.UnivariateRing(Rings.Z),
                Rings.UnivariateRingZp64(17))) {
            MultivariateRing<?> inner2 = Rings.MultivariateRing(2, inner1);
            UnivariateRing<?> inner3 = Rings.UnivariateRing(inner2);
            MultivariateRing<?> ring = Rings.MultivariateRing(2, inner3);

            for (int i = 0; i < its(8, 16); ++i) {
                AMultivariatePolynomial
                        a = ring.randomElement(3, 5),
                        b = ring.randomElement(3, 5),
                        c = ring.randomElement(3, 5);

                a.multiply((IPolynomial) c);
                b.multiply((IPolynomial) c);

                AMultivariatePolynomial gcd = MultivariateGCD.PolynomialGCD(a, b);
                Assert.assertTrue(dividesQ(a, gcd));
                Assert.assertTrue(dividesQ(b, gcd));
                Assert.assertTrue(dividesQ(gcd, c));
            }
        }
    }

    @Test // =====> testModularGCDInZ1   elapsed 2ms
    public void testModularGCDInZ1() throws Exception {
        BiFunction<MultivariatePolynomialZp64, MultivariatePolynomialZp64, MultivariatePolynomialZp64>
                modAlg = MultivariateGCD::ZippelGCD;

        MultivariatePolynomial<BigInteger>
                a = parse("x*y*z^2 + y*z*x^6 + 1"),
                b = parse("x*y*z^2 + 12342134234233*y*z*x^6 - 123*z + x"),
                g = parse("x^6*y + 1212423413*y*z^7*x^6 - 1231241234234*z + 164287246876423*y"),
                ag = a.clone().multiply(g),
                bg = b.clone().multiply(g);


        assertEquals(g.primitivePart(), ModularGCDInZ(ag, bg, modAlg).primitivePart());
    }

    @Test // =====> testModularGCDInZRandom1   elapsed 39s
    public void testModularGCDInZRandom1() throws Exception {
        BiFunction<MultivariatePolynomialZp64, MultivariatePolynomialZp64, MultivariatePolynomialZp64>
                modAlg = MultivariateGCD::ZippelGCD;

        BiFunction<MultivariatePolynomial<BigInteger>, MultivariatePolynomial<BigInteger>, MultivariatePolynomial<BigInteger>>
                zAlg = (a, b) -> MultivariateGCD.ModularGCDInZ(a, b, modAlg);

        int nIterations = its(1000, 1000);
        RandomGenerator rnd = getRandom();
        GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> sampleData =
                new GCDSampleDataGeneric<>(Rings.Z, 3, 5, 5, 15, 5, 15, rnd);

        testGCDAlgorithms(sampleData, nIterations, GCDAlgorithm.named("Modular gcd in Z", zAlg));
    }

    @Test // =====> testEEZGCD8   elapsed 31ms
    public void testEEZGCD8() throws Exception {
        IntegersZp64 cfRing = Rings.Zp64(7);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("y^3+y^4+x*y^4+5*x*y^5+4*x^2+6*x^2*y^2+4*x^2*y^3+4*x^2*y^4+3*x^2*y^5+6*x^2*y^6+x^3+x^3*y+5*x^3*y^2+6*x^3*y^3+x^3*y^4+5*x^3*y^5+6*x^3*y^6+5*x^4*y+4*x^4*y^2+4*x^4*y^3+2*x^4*y^4+2*x^4*y^5+4*x^4*y^6+6*x^4*y^7+4*x^4*y^8+2*x^5+2*x^5*y+x^5*y^2+3*x^5*y^3+x^5*y^4+x^5*y^5+3*x^5*y^6+3*x^5*y^8+3*x^5*y^9+4*x^7*y+5*x^7*y^2+2*x^7*y^3+3*x^7*y^4+4*x^7*y^5+3*x^7*y^6+5*x^7*y^7+5*x^7*y^8+3*x^7*y^9+4*x^7*y^10+x^7*y^11+5*x^8+5*x^8*y+4*x^8*y^2+6*x^8*y^3+x^8*y^4+3*x^8*y^5+6*x^8*y^6+3*x^8*y^7+4*x^8*y^8+x^8*y^9+3*x^8*y^10+5*x^8*y^11+6*x^9+6*x^9*y+6*x^9*y^2+4*x^9*y^3+6*x^9*y^4+3*x^9*y^5+5*x^9*y^6+4*x^9*y^7+4*x^9*y^8+3*x^9*y^9+4*x^9*y^11+6*x^9*y^12+5*x^10*y+2*x^10*y^3+4*x^10*y^4+2*x^10*y^5+2*x^10*y^6+x^10*y^7+2*x^10*y^8+3*x^10*y^10+4*x^10*y^11+5*x^11*y+4*x^11*y^2+3*x^11*y^3+6*x^11*y^4+5*x^11*y^6+5*x^11*y^7+x^11*y^8+5*x^11*y^9+6*x^11*y^10+2*x^12*y^6+5*x^12*y^8+3*x^12*y^9+3*x^12*y^10+x^12*y^11", cfRing),
                b = MultivariatePolynomialZp64.parse("2*y^2+2*y^3+3*x*y^2+3*x*y^4+2*x^2*y^3+4*x^2*y^5+4*x^3+4*x^3*y+3*x^3*y^2+2*x^3*y^4+x^3*y^5+6*x^3*y^6+x^4+5*x^4*y+3*x^4*y^2+2*x^4*y^3+4*x^4*y^5+x^4*y^6+x^5+6*x^5*y+2*x^5*y^2+4*x^5*y^3+6*x^5*y^4+6*x^5*y^5+x^5*y^6+3*x^5*y^7+x^5*y^8+6*x^6+6*x^6*y+2*x^6*y^2+x^6*y^3+6*x^6*y^4+5*x^6*y^6+4*x^6*y^7+6*x^6*y^8+5*x^6*y^9+2*x^7+x^7*y+3*x^7*y^2+x^7*y^4+x^7*y^6+4*x^7*y^8+x^7*y^9+2*x^7*y^10+4*x^8+x^8*y+4*x^8*y^3+x^8*y^4+2*x^8*y^5+2*x^8*y^6+6*x^8*y^7+2*x^8*y^10+3*x^8*y^11+3*x^9+x^9*y+x^9*y^2+4*x^9*y^3+2*x^9*y^4+5*x^9*y^5+2*x^9*y^6+3*x^9*y^7+6*x^9*y^9+4*x^9*y^10+4*x^9*y^11+x^10+3*x^10*y+3*x^10*y^2+5*x^10*y^3+4*x^10*y^4+6*x^10*y^5+5*x^10*y^6+3*x^10*y^7+6*x^10*y^8+3*x^10*y^9+3*x^10*y^10+6*x^10*y^12+3*x^11+2*x^11*y+5*x^11*y^2+2*x^11*y^3+2*x^11*y^4+x^11*y^5+5*x^11*y^6+3*x^11*y^7+5*x^11*y^8+4*x^11*y^9+x^11*y^10+3*x^11*y^11+x^12+6*x^12*y+x^12*y^3+5*x^12*y^4+6*x^12*y^5+3*x^12*y^6+6*x^12*y^7+4*x^12*y^8+x^12*y^9+5*x^12*y^10+2*x^13*y^5+x^13*y^6+2*x^13*y^7+3*x^13*y^9+x^13*y^10+4*x^13*y^11+3*x^14*y^5+5*x^14*y^6+5*x^14*y^7+x^14*y^9+2*x^14*y^10", cfRing),
                g = a.clone().subtract(b);
        assertEquals(ZippelGCD(a, b), EEZGCD(a, b));
        a.multiply(g);
        b.multiply(g);
        assertEquals(ZippelGCD(a, b), EEZGCD(a, b));
    }

    @Ignore
    @Benchmark
    @Test
    public void testEEZGCD_random4() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 4, nVarsMax = 4, minDegree = 1, maxDegree = 5, minSize = 550, maxSize = 560;
        int nIterations = its(10, 200);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD)
        );
    }

    @Ignore("issue #20")
    @Test
    public void testSmallDomain6() throws Exception {
        IntegersZp64 ring = Rings.Zp64(3);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1 + a + 5*b + 7*c + d + 11*e + 13*f + g", ring),
                b = MultivariatePolynomialZp64.parse("1 - a + 5*b - 7*c + d - 11*e + 13*f - g", ring),
                g = MultivariatePolynomialZp64.parse("1 + a - 5*b + 7*c - d + 11*e - 13*f + g", ring);

        int exp = 7;
        a = PolynomialMethods.polyPow(a, exp);
        b = PolynomialMethods.polyPow(b, exp);
        g = PolynomialMethods.polyPow(g, exp);

        MultivariatePolynomialZp64 ag = a.clone().multiply(g);
        MultivariatePolynomialZp64 bg = b.clone().multiply(g);

        System.out.println(PolynomialGCD(ag, bg));
    }

    @Ignore("too expensive (memory limit in CI)")
    @Test
    public void testHugePoly1() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> R = Rings.MultivariateRingZp64(3, SmallPrimes.nextPrime(1 << 25));
        for (int j = 1; j < 2; ++j) {
            try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                    new FileInputStream(
                            "src/test/resources/cc/redberry/rings/poly/multivar/HugeGCDinZp.txt.gz"))))) {
                String s1 = in.readLine();
                String s2 = in.readLine();

                MultivariatePolynomialZp64 p = R.parse(s1);
                MultivariatePolynomialZp64 q = R.parse(s2);
                long t1 = System.nanoTime();
                MultivariatePolynomialZp64 gcd = MultivariateGCD.ZippelGCD(p, q);
                long t2 = System.nanoTime();

                System.out.println("Total time : " + nanosecondsToString(t2 - t1));
                System.out.println("GCD size   : " + gcd.size());
                System.out.println();
            }
        }

        //Total time : 3s
        //GCD size   : 263
    }

    @Ignore("too expensive (memory limit in CI)")
    @Test
    public void testHugePoly2() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> R = Rings.MultivariateRing(3, Rings.Z);
        for (int j = 1; j < 4; ++j) {
            try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                    new FileInputStream(
                            "src/test/resources/cc/redberry/rings/poly/multivar/HugeGCDinZ.txt.gz"))))) {
                String s1 = in.readLine();
                String s2 = in.readLine();

                MultivariatePolynomial<BigInteger> p = R.parse(s1);
                MultivariatePolynomial<BigInteger> q = R.parse(s2);
                long t1 = System.nanoTime();
                MultivariatePolynomial<BigInteger> gcd = MultivariateGCD.PolynomialGCD(p, q);
                long t2 = System.nanoTime();

                System.out.println("Total time : " + nanosecondsToString(t2 - t1));
                System.out.println("GCD size   : " + gcd.size());
                System.out.println();
            }
        }

        //Total time : 113s
        //GCD size   : 3380
    }

    @Ignore("too expensive (memory limit in CI)")
    @Test
    public void testHugePoly3() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> R = Rings.MultivariateRing(7, Rings.Z);
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                new FileInputStream(
                        "src/test/resources/cc/redberry/rings/poly/multivar/HugeGCDinZ.LinZip.txt.gz"))))) {
            int iProblem = 0;
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#"))
                    continue;
                String[] split = line.split("\t");

                MultivariatePolynomial<BigInteger>
                        p = R.parse(split[1]),
                        q = R.parse(split[2]),
                        expected = R.parse(split[3]);

                long t1 = System.nanoTime();
                MultivariatePolynomial<BigInteger> gcd = MultivariateGCD.PolynomialGCD(p, q);
                long t2 = System.nanoTime();
                assertTrue(dividesQ(gcd, expected));

                System.out.println("Problem    : " + iProblem);
                System.out.println("Total time : " + nanosecondsToString(t2 - t1));
                System.out.println("GCD size   : " + gcd.size());
                System.out.println();
                ++iProblem;
            }
        }

        //Problem    : 0
        //Total time : 30s
        //GCD size   : 201
        //
        //Problem    : 1
        //Total time : 36s
        //GCD size   : 201
        //
        //Problem    : 2
        //Total time : 39s
        //GCD size   : 201
        //
        //Problem    : 3
        //Total time : 31s
        //GCD size   : 201
        //
        //Problem    : 4
        //Total time : 29s
        //GCD size   : 201
    }

    @Ignore("too expensive (memory limit in CI)")
    @Test
    public void testHugePoly4() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>>
                R = Rings.MultivariateRing(3, Rings.Z),
                modR = Rings.MultivariateRing(3, Rings.Zp(1073741827));
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                new FileInputStream(
                        "src/test/resources/cc/redberry/rings/poly/multivar/HugeHugeGCDinZ.txt.gz"))))) {
            int iProblem = 0;
            String line;
            while ((line = in.readLine()) != null) {
                if (line.startsWith("#"))
                    continue;
                String[] split = line.split("\t");

                MultivariatePolynomial<BigInteger>
                        p = R.parse(split[1]),
                        q = R.parse(split[2]),
                        expected = R.parse(split[3]);


                System.out.println("Problem    : " + iProblem++);
                long t1 = System.nanoTime();
                MultivariatePolynomial<BigInteger> gcd = MultivariateGCD.PolynomialGCD(p, q);
                long t2 = System.nanoTime();
                assertTrue(dividesQ(gcd, expected));

                System.out.println("Total time (Z)   : " + nanosecondsToString(t2 - t1));
                System.out.println("GCD size         : " + gcd.size());
                System.out.println();

                t1 = System.nanoTime();
                MultivariatePolynomial<BigInteger> mgcd = MultivariateGCD.PolynomialGCD(modR.valueOf(p), modR.valueOf(q));
                t2 = System.nanoTime();
                assertTrue(dividesQ(gcd, expected));

                System.out.println("Total time (mod) : " + nanosecondsToString(t2 - t1));
                System.out.println("MGCD size        : " + mgcd.size());
                System.out.println();
            }
        }
        // > 1444s
        // 355s
    }

    @Test(timeout = 100_000L) // =====> testMediumCharacteristic1   elapsed 12s
    public void testMediumCharacteristic1() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(5, Rings.Zp(4099));
        try (BufferedReader in = new BufferedReader(new InputStreamReader(new GZIPInputStream(
                MultivariateGCDTest.class.getClassLoader().getResourceAsStream(
                        "MediumCharacteristicHugePoly.txt.gz"))))) {

            MultivariatePolynomial<BigInteger> a = ring.parse(in.readLine());
            MultivariatePolynomial<BigInteger> b = ring.parse(in.readLine());

            for (boolean switch64 : Arrays.asList(true, false)) {
                Conversions64bit.SWITCH_TO_64bit = switch64;
                long start = System.nanoTime();
                MultivariatePolynomial<BigInteger> gcd = PolynomialGCD(a, b);
                System.out.println(nanosecondsToString(System.nanoTime() - start));
            }
        }
        //2s
        //11s
    }

    @Test
    public void testTrivial1() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(7, Rings.Z);
        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);
        MultivariatePolynomial<BigInteger> a =
                RandomMultivariatePolynomials.randomPolynomial(ring.nVariables(), 15, 20,
                        2500, Rings.Z, MonomialOrder.GREVLEX, r -> Rings.Z.valueOf(1 + r.nextInt(100)), rnd);

        MultivariatePolynomial<BigInteger> b =
                RandomMultivariatePolynomials.randomPolynomial(ring.nVariables(), 15, 20,
                        2, Rings.Z, MonomialOrder.GREVLEX, r -> Rings.Z.valueOf(1 + r.nextInt(100)), rnd);

        for (int i = 0; i < 10; ++i) {
            long start = System.nanoTime();
            PolynomialGCD(a, b).size();
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
        // 825us
        // 1850us
        // 1911us
        // 2ms
        // 2ms
    }

    @Test
    public void testTrivial2() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(7, Rings.Z);
        MultivariatePolynomial<BigInteger> a = ring.parse("x1*x2^3*x3*x4^3*x5^4-1*x2^5*x3*x4^3*x5^4-1*x1^7*x3^4*x4^3*x5*x6^2*x7^4+x1^6*x2^2*x3^4*x4^3*x5*x6^2*x7^4+x1^7*x3^4*x4^3*x5*x6^4*x7^3+x1^5*x3^2*x4^6*x5^2*x6^6*x7+x1^8*x3*x4^2*x5*x6^7*x7^4+x1^8*x2*x3*x4^3*x6^6*x7^4-1*x1^6*x2^2*x3^4*x4^3*x5*x6^4*x7^3-1*x1^4*x2^2*x3^2*x4^6*x5^2*x6^6*x7-1*x1^5*x3^2*x4^6*x5^2*x6^8-1*x1^7*x2^2*x3*x4^2*x5*x6^7*x7^4-1*x1^9*x3*x4^3*x5*x6^6*x7^4+x1^7*x3*x4^6*x6^6*x7^4-1*x1^7*x2^3*x3*x4^3*x6^6*x7^4-1*x1^8*x3*x4^2*x5*x6^9*x7^3-1*x1^8*x2*x3*x4^3*x6^8*x7^3+x1^4*x2^2*x3^2*x4^6*x5^2*x6^8+x1^8*x2^2*x3*x4^3*x5*x6^6*x7^4+x1^7*x2^2*x3*x4^2*x5*x6^9*x7^3+x1^9*x3*x4^3*x5*x6^8*x7^3-1*x1^7*x3*x4^6*x6^8*x7^3+x1^7*x2^3*x3*x4^3*x6^8*x7^3+x1^4*x2^3*x3^5*x5*x6^6*x7^7-1*x1^8*x2^2*x3*x4^3*x5*x6^8*x7^3-1*x1^3*x2^5*x3^5*x5*x6^6*x7^7-1*x1^4*x2^3*x3^5*x5*x6^8*x7^6+x1^8*x3^5*x4^3*x5^2*x6^5*x7^4+x1^3*x2^5*x3^5*x5*x6^8*x7^6+x1^8*x3^5*x4^4*x5*x6^6*x7^4-1*x1^7*x2^2*x3^5*x4^3*x5^2*x6^5*x7^4-1*x1^8*x3^5*x4^3*x5^2*x6^7*x7^3-1*x1^7*x2^2*x3^5*x4^4*x5*x6^6*x7^4-1*x1^8*x3^5*x4^4*x5*x6^8*x7^3+x1^7*x2^2*x3^5*x4^3*x5^2*x6^7*x7^3-1*x1^8*x2^6*x3^2*x4^3*x5*x6^6*x7^4+x1^7*x2^2*x3^5*x4^4*x5*x6^8*x7^3+x1^7*x2^8*x3^2*x4^3*x5*x6^6*x7^4+x1^8*x2^6*x3^2*x4^3*x5*x6^8*x7^3-1*x1^7*x2^8*x3^2*x4^3*x5*x6^8*x7^3-1*x1^11*x4^3*x5^5*x6^10*x7^4+x1^10*x2^2*x4^3*x5^5*x6^10*x7^4+x1^11*x4^3*x5^5*x6^12*x7^3-1*x1^10*x2^2*x4^3*x5^5*x6^12*x7^3");
        MultivariatePolynomial<BigInteger> b = ring.parse("-1*x1^8*x3*x4^2*x5*x6^6*x7^4+x1^7*x2^2*x3*x4^2*x5*x6^6*x7^4+x1^8*x3*x4^2*x5*x6^8*x7^3-1*x1^7*x2^2*x3*x4^2*x5*x6^8*x7^3");
        for (int i = 0; i < 10; ++i) {
            long start = System.nanoTime();
            PolynomialGCD(a, b).size();
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }

        // 273us
        // 306us
        // 882us
        // 457us
        // 351us
        // 276us
    }

    @Test
    public void testAlgExt1() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(UnivariatePolynomial.create(Q, Q.valueOf(-2), Q.valueOf(0), Q.valueOf(1)));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> mRing = Rings.MultivariateRing(3, field);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(mRing, cfCoder, "x", "y", "z");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + s*x^5 - s*y + z) * ( 1 + s*x^2 + 12*x^5)");
        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + s*x^5 - s*y + z) * ( 14 - s*x + 2*x^17)");
        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> gcd = PolynomialGCD(a, b);
        assertTrue(dividesQ(a, gcd));
        assertTrue(dividesQ(b, gcd));
    }

    @Test
    public void testAlgExt2() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = AlgebraicExtension(UnivariatePolynomial.create(13, 0, 1, 2, 3, 4, 15, 19).mapCoefficients(Q, Q::mkNumerator));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> mRing = Rings.MultivariateRing(3, field);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(mRing, cfCoder, "x", "y", "z");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + s*x^5 - s*y^5*z + (1 - s^6)*z^3*x) * ( 1 + s*x^2 + 12*x^5 + y)^5");
        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + s*x^5 - s*y^5*z + (1 - s^6)*z^3*x) * ( 14 - s*x + 2*s*x^17 + y - z)^5");

        for (int i = 1; i < 10; ++i) {
            long start = System.nanoTime();
            MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> gcd = ZippelGCDInNumberFieldViaRationalReconstruction(a, b);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));
            assertTrue(dividesQ(a, gcd));
            assertTrue(dividesQ(b, gcd));
        }

        //Modular: 1028ms
        //Modular: 303ms
        //Modular: 268ms
        //Modular: 206ms
        //Modular: 216ms
        //Modular: 197ms
        //Modular: 182ms
        //Modular: 241ms
        //Modular: 168ms
    }

    @Test
    public void testZippelModularAlgExt1_random() {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < 10; ++i) {
            UnivariatePolynomial<BigInteger> minimalPolyZ =
                    IrreduciblePolynomials.randomIrreduciblePolynomialOverZ(rndd.nextInt(2, 5), rnd);
            minimalPolyZ.setLC(Z.valueOf(rndd.nextInt(1, 50)));
            if (!IrreduciblePolynomials.irreducibleQ(minimalPolyZ)) {
                --i;
                continue;
            }
            UnivariatePolynomial<Rational<BigInteger>> minimalPoly = minimalPolyZ.mapCoefficients(Q, Q::mkNumerator);
            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> numberField = new AlgebraicNumberField<>(minimalPoly);

            GCDSampleDataGeneric<UnivariatePolynomial<Rational<BigInteger>>> source
                    = new GCDSampleDataGeneric<>(numberField, 3, 6, 10, 20, 10, 30, rnd);

            source.rndCoefficients = rand -> numberField.valueOf(RandomUnivariatePolynomials.
                    randomPoly(minimalPoly.degree(), Q, __ -> Q.mk(rand.nextInt(100), 1 + rand.nextInt(10)), rand));
            System.out.println("\n\n");
            System.out.println("Ring: " + numberField);
            testGCDAlgorithms(source, its(10, 30),
                    GCDAlgorithm.named("Zippel with rational reconstruction", MultivariateGCD::ZippelGCDInNumberFieldViaRationalReconstruction),
                    GCDAlgorithm.named("Zippel with Langemyr & McCallum", MultivariateGCD::ZippelGCDInNumberFieldViaLangemyrMcCallum),
                    GCDAlgorithm.named("Modular with rational reconstruction", (a, b) -> ModularGCDInNumberFieldViaRationalReconstruction(a, b, MultivariateGCD::PolynomialGCD)),
                    GCDAlgorithm.named("Modular with Langemyr & McCallum", (a, b) -> ModularGCDInNumberFieldViaLangemyrMcCallum(a, b, MultivariateGCD::PolynomialGCD))
            );
        }
    }

    @Test
    @Benchmark(runAnyway = true)
    public void testZippelModularAlgExt_2() {
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly =
                UnivariatePolynomial.create(31229703,
                        31584466, 9500649, 9480702,
                        23265262, 5568454, 3392530,
                        30401154, 15298203, 25411939,
                        30401154, 15298203, 25411939,
                        31584466, 9500649, 9480702,
                        30401154, 15298203, 25411939,
                        30401154, 15298203, 25411939, 1)
                        .mapCoefficients(Q, Q::valueOfBigInteger);
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = MultivariateRing(5, field);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(uRing, cfCoder, "x", "y", "z", "p", "q");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("-2629984078 - 2747848492*s - 826556509*s^2 - 824821066*s^3 - 2024077758*s^4 - 484455432*s^5 - 295150100*s^6 - 2644900377*s^7 - 1330943630*s^8 - 2210838750*s^9 + (-2539295050 - 2653095078*s - 798054533*s^2 - 796378949*s^3 - 1954282008*s^4 - 467750096*s^5 - 284972456*s^6 - 2553696896*s^7 - 1285048960*s^8 - 2134602853*s^9)*x*p + (2750903023 + 2874186377*s + 864559081*s^2 + 862743871*s^3 + 2117138752*s^4 + 506729319*s^5 + 308720296*s^6 + 2766505047*s^7 + 1392136388*s^8 + 2312486542*s^9)*x^2*y^3*q + (-2902051419 - 3032108714*s - 912062266*s^2 - 910147339*s^3 - 2233465130*s^4 - 534571658*s^5 - 325682847*s^6 - 2918510751*s^7 - 1468627461*s^8 - 2439546109*s^9)*z^3 + (-2025390071 - 2116159152*s - 636543408*s^2 - 635207023*s^3 - 1558772527*s^4 - 373086362*s^5 - 227299470*s^6 - 2036877303*s^7 - 1024979692*s^8 - 1702599819*s^9)*x*y^4 + (-2962510945 - 3095277661*s - 931063590*s^2 - 929108773*s^3 - 2279995709*s^4 - 545708558*s^5 - 332467948*s^6 - 2979313000*s^7 - 1499223872*s^8 - 2490370012*s^9)*q^5*x + (-1934700903 - 2021405804*s - 608041537*s^2 - 606764926*s^3 - 1488976751*s^4 - 356381007*s^5 - 217121886*s^6 - 1945673813*s^7 - 979084946*s^8 - 1626364082*s^9)*z^6*x*p*y^7*q + (-997580106 - 1042287350*s - 313521458*s^2 - 312863069*s^3 - 767753632*s^4 - 183758950*s^5 - 111953411*s^6 - 1003238154*s^7 - 504840626*s^8 - 838594076*s^9)*z^7 + (-1178958351 - 1231794140*s - 370525402*s^2 - 369747315*s^3 - 907345216*s^4 - 217169661*s^5 - 132308571*s^6 - 1185645024*s^7 - 596629988*s^8 - 991065545*s^9)*p^8*q + (-2025390077 - 2116159255*s - 636543412*s^2 - 635207025*s^3 - 1558772532*s^4 - 373086332*s^5 - 227299454*s^6 - 2036877243*s^7 - 1024979650*s^8 - 1702599890*s^9)*p^9*x^3*z^2 + (1118499087 + 1168625309*s + 351524046*s^2 + 350785991*s^3 + 860814787*s^4 + 206032829*s^5 + 125523644*s^6 + 1124842769*s^7 + 566033562*s^8 + 940241738*s^9)*p^5*y^5");
        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("-2297457369 - 2400419408*s - 722049354*s^2 - 720533283*s^3 - 1768159825*s^4 - 423202465*s^5 - 257832326*s^6 - 2310487767*s^7 - 1162663504*s^8 - 1931307309*s^9 + (-1118498918 - 1168625170*s - 351523998*s^2 - 350785953*s^3 - 860814620*s^4 - 206032844*s^5 - 125523543*s^6 - 1124842734*s^7 - 566033484*s^8 - 940241840*s^9)*p + (-634823841 - 663273730*s - 199513585*s^2 - 199094740*s^3 - 488570500*s^4 - 116937491*s^5 - 71243056*s^6 - 638424217*s^7 - 321262269*s^8 - 533650647*s^9)*y^2*z + (-120918804 - 126337886*s - 38002569*s^2 - 37922791*s^3 - 93061011*s^4 - 22273824*s^5 - 13570040*s^6 - 121604571*s^7 - 61192741*s^8 - 101647722*s^9)*x^3*q^2 + (-2629984238 - 2747848535*s - 826556405*s^2 - 824821002*s^3 - 2024077765*s^4 - 484455463*s^5 - 295150103*s^6 - 2644900316*s^7 - 1330943621*s^8 - 2210838616*s^9)*x^2*y*z*p + (60459502 + 63168844*s + 19001371*s^2 + 18961478*s^3 + 46530559*s^4 + 11136967*s^5 + 6785128*s^6 + 60802294*s^7 + 30596476*s^8 + 50823949*s^9)*x^2*y*z*p^5 + (-2599754450 - 2716264170*s - 817055796*s^2 - 815340331*s^3 - 2000812621*s^4 - 478886953*s^5 - 291757547*s^6 - 2614499288*s^7 - 1315645534*s^8 - 2185426692*s^9)*q*z^6*p + (-1964930626 - 2052990342*s - 617542144*s^2 - 616245697*s^3 - 1512242053*s^4 - 361949415*s^5 - 220514397*s^6 - 1976074963*s^7 - 994383244*s^8 - 1651775954*s^9)*z^7 + (-1571944495 - 1642392288*s - 494033670*s^2 - 492996488*s^3 - 1209793576*s^4 - 289559689*s^5 - 176411479*s^6 - 1580860073*s^7 - 795506501*s^8 - 1321420729*s^9)*y^8 + (2811362414 + 2937355383*s + 883560380*s^2 + 881705200*s^3 + 2163669280*s^4 + 517866158*s^5 + 315505366*s^6 + 2827307413*s^7 + 1422732958*s^8 + 2363310334*s^9)*p^9*x + (2116079261 + 2210912639*s + 665045462*s^2 + 663649187*s^3 + 1628568308*s^4 + 389791803*s^5 + 237477108*s^6 + 2128080876*s^7 + 1070874245*s^8 + 1778835730*s^9)*x^5*z^6 + (-2236998012 - 2337250400*s - 703047951*s^2 - 701571937*s^3 - 1721629439*s^4 - 412065534*s^5 - 251047253*s^6 - 2249685310*s^7 - 1132067031*s^8 - 1880483482*s^9)*q^2*y*z*p");
        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> g = coder.parse("-2811362312 - 2937355319*s - 883560337*s^2 - 881705356*s^3 - 2163669316*s^4 - 517866198*s^5 - 315505348*s^6 - 2827307397*s^7 - 1422732833*s^8 - 2363310305*s^9 + (-1027809863 - 1073871832*s - 323021979*s^2 - 322343852*s^3 - 791018874*s^4 - 189327442*s^5 - 115345990*s^6 - 1033639203*s^7 - 520138994*s^8 - 864005935*s^9)*x*z + (755742670 + 789611730*s + 237516324*s^2 + 237017568*s^3 + 581631609*s^4 + 139211399*s^5 + 84813301*s^6 + 760028909*s^7 + 382455125*s^8 + 635298486*s^9)*z^2*y^2 + (-2569524725 - 2684679604*s - 807555125*s^2 - 805859617*s^3 - 1977547175*s^4 - 473318570*s^5 - 288364973*s^6 - 2584098178*s^7 - 1300347246*s^8 - 2160014731*s^9)*x^2*y*z*p^5 + (-90689044 - 94753334*s - 28502018*s^2 - 28442203*s^3 - 69795734*s^4 - 16705297*s^5 - 10177621*s^6 - 91203421*s^7 - 45894548*s^8 - 76235851*s^9)*z^4*y^6*q*p + (-60459337 - 63168954*s - 19001279*s^2 - 18961469*s^3 - 46530495*s^4 - 11136820*s^5 - 6785055*s^6 - 60802219*s^7 - 30596432*s^8 - 50823814*s^9)*q^5*y + (-2206768278 - 2305665984*s - 693547285*s^2 - 692091191*s^3 - 1698364119*s^4 - 406497048*s^5 - 247654596*s^6 - 2219284280*s^7 - 1116768829*s^8 - 1855071634*s^9)*x^6*y^3*z^2*p^2*q + (-2992740680 - 3126862134*s - 940564345*s^2 - 938589464*s^3 - 2303260921*s^4 - 551276902*s^5 - 335860390*s^6 - 3009714187*s^7 - 1514522192*s^8 - 2515781892*s^9)*x^2*y*z*p + (-937120739 - 979118398*s - 294520061*s^2 - 293901791*s^3 - 721223167*s^4 - 172622016*s^5 - 105168421*s^6 - 942435722*s^7 - 474244292*s^8 - 787770068*s^9)*q^2*y*z*p + (-2569524839 - 2684679579*s - 807555086*s^2 - 805859625*s^3 - 1977547366*s^4 - 473318589*s^5 - 288365123*s^6 - 2584098088*s^7 - 1300347180*s^8 - 2160014805*s^9)*x^2*y*z*q + (-302297052 - 315844685*s - 95006420*s^2 - 94806952*s^3 - 232652706*s^4 - 55684503*s^5 - 33925202*s^6 - 304011605*s^7 - 152982033*s^8 - 254119332*s^9)*z^10*y + (362756428 + 379013624*s + 114007873*s^2 + 113768522*s^3 + 279183228*s^4 + 66821400*s^5 + 40710454*s^6 + 364813756*s^7 + 183578437*s^8 + 304943317*s^9)*x^5*z^5*q + (-2327687116 - 2432003792*s - 731549932*s^2 - 730014098*s^3 - 1791425234*s^4 - 428770976*s^5 - 261224810*s^6 - 2340888824*s^7 - 1177961606*s^8 - 1956719264*s^9)*x^12 + (-2962510886 - 3095277709*s - 931063515*s^2 - 929108855*s^3 - 2279995592*s^4 - 545708444*s^5 - 332467892*s^6 - 2979313007*s^7 - 1499223814*s^8 - 2490369949*s^9)*p^5*q^5");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                ag = a.clone().multiply(g),//.multiply(a.clone().decrement()),
                bg = b.clone().multiply(g);//.multiply(b.clone().decrement());

        for (int i = 0; i < 2; ++i) {
            long start;
            MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> mod;

            start = System.nanoTime();
            mod = ZippelGCDInNumberFieldViaRationalReconstruction(ag, bg);
            System.out.println("Rational reconstruction : " + nanosecondsToString(System.nanoTime() - start));
            assertTrue(dividesQ(mod, g));

            start = System.nanoTime();
            mod = ZippelGCDInNumberFieldViaLangemyrMcCallum(ag, bg);
            System.out.println("Langemyr & McCallum     : " + nanosecondsToString(System.nanoTime() - start));
            assertTrue(dividesQ(mod, g));
            System.out.println();
        }
        // Rational reconstruction : 5s
        // Langemyr & McCallum     : 15s
        // Rational reconstruction : 1868ms
        // Langemyr & McCallum     : 16s
    }

    /* =============================================== Test data =============================================== */

    /** sample data for test of GCD */
    public static final class GCDSample<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** sample polynomials */
        public final Poly a, b;
        /** their GCD */
        public final Poly gcd;
        /** a/gcd and b/gcd */
        public final Poly aCoFactor, bCoFactor;
        /** either generic Ring or IntegersZp64 */
        public final Object domain;

        public GCDSample(Poly aCoFactor, Poly bCoFactor, Poly gcd) {
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
            this.gcd = gcd;
            this.a = aCoFactor.clone().multiply(gcd);
            this.b = bCoFactor.clone().multiply(gcd);
            checkConsistency(a, b, gcd, aCoFactor, bCoFactor);
            this.domain = (gcd instanceof MultivariatePolynomialZp64)
                    ? ((MultivariatePolynomialZp64) gcd).ring
                    : ((MultivariatePolynomial) gcd).ring;
        }

        @Override
        public String toString() {
            return new StringBuilder()
                    .append("\naCoFactor = " + aCoFactor)
                    .append("\nbCoFactor = " + bCoFactor)
                    .append("\ngcd       = " + gcd)
                    .append("\nring    = " + domain)
                    .toString();
        }
    }

    /** substitute random values for some variables */
    public static GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>>
    boundCoefficients(final GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> source, BigInteger bound) {
        final IntegersZp domain = new IntegersZp(bound);
        return new GCDSampleData<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>>() {
            @Override
            GCDSample<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> nextSample0(boolean primitive, boolean monic) {
                GCDSample<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> sample = source.nextSample(primitive, monic);
                return new GCDSample<>(
                        setRandomCoefficients(sample.aCoFactor, bound),
                        setRandomCoefficients(sample.bCoFactor, bound),
                        setRandomCoefficients(sample.gcd, bound));
            }

            private MultivariatePolynomial<BigInteger> setRandomCoefficients(MultivariatePolynomial<BigInteger> poly, BigInteger bound) {
                RandomGenerator rnd = PrivateRandom.getRandom();
                MultivariatePolynomial<BigInteger> r = poly.createZero();
                for (Monomial<BigInteger> term : poly)
                    r.add(term.setCoefficient(RandomUtil.randomInt(bound, rnd)));
                return r;
            }

            @Override
            public String toString() {
                return source.toString();
            }
        };
    }

    /** substitute random values for some variables */
    public static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    GCDSampleData<Term, Poly> filterZeros(final GCDSampleData<Term, Poly> data) {
        return new GCDSampleData<Term, Poly>() {
            @Override
            GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic) {
                GCDSample<Term, Poly> sample;
                do { sample = data.nextSample(primitive, monic);}
                while (sample.b.isZero() || sample.a.isZero());
                return sample;
            }

            @Override
            public String toString() {
                return data.toString();
            }
        };
    }

    /** substitute random values for some variables */
    public static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    GCDSampleData<Term, Poly> fixVariables(final GCDSampleData<Term, Poly> data, int nActuallyUsedVariables) {
        return new GCDSampleData<Term, Poly>() {
            @Override
            GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic) {
                GCDSample<Term, Poly> sample = data.nextSample(primitive, monic);

                RandomGenerator rnd = getRandom();
                Poly gcd = sample.gcd;
                do {
                    gcd = gcd.evaluateAtRandom(rnd.nextInt(gcd.nVariables), rnd);
                } while (gcd.nUsedVariables() > nActuallyUsedVariables);

                Poly
                        aCoFactor = sample.aCoFactor,
                        bCoFactor = sample.bCoFactor;
                int[] gcdDegrees = gcd.degrees();
                for (int i = 0; i < gcdDegrees.length; i++) {
                    if (gcdDegrees[i] == 0) {
                        if (rnd.nextBoolean())
                            aCoFactor = aCoFactor.evaluateAtRandom(i, rnd);
                        else
                            bCoFactor = bCoFactor.evaluateAtRandom(i, rnd);
                    }
                }

                return new GCDSample<>(aCoFactor, bCoFactor, gcd);
            }

            @Override
            public String toString() {
                return data.toString() + ", nUsedVariables = " + nActuallyUsedVariables;
            }
        };
    }

    /** probable coprime polynomials */
    public static <Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    GCDSampleData<Term, Poly> coprimeData(final GCDSampleData<Term, Poly> data) {
        return new GCDSampleData<Term, Poly>() {
            @Override
            GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic) {
                GCDSample<Term, Poly> sample = data.nextSample(primitive, monic);
                return new GCDSample<>(sample.aCoFactor, sample.bCoFactor, sample.aCoFactor.createOne());
            }

            @Override
            public String toString() {
                return data.toString();
            }
        };
    }

    /** provides random data for GCD tests */
    public static abstract class GCDSampleData<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {

        final DescriptiveStatistics
                medFactorsSize = new DescriptiveStatistics(),
                medFactorsDegree = new DescriptiveStatistics(),
                medFactorsSparsity = new DescriptiveStatistics(),
                medFactorsSparsity2 = new DescriptiveStatistics(),
                medGCDSize = new DescriptiveStatistics(),
                medGCDSparsity = new DescriptiveStatistics(),
                medGCDSparsity2 = new DescriptiveStatistics(),
                medGCDDegree = new DescriptiveStatistics(),
                medFactorsUsedVariables = new DescriptiveStatistics(),
                medGCDUsedVariables = new DescriptiveStatistics();

        public boolean primitive = false, monic = false;

        final GCDSample<Term, Poly> nextSample() {
            return nextSample(primitive, monic);
        }

        final GCDSample<Term, Poly> nextSample(boolean primitive, boolean monic) {
            GCDSample<Term, Poly> sample;
            do {sample = nextSample0(primitive, monic);} while (sample.gcd.isZero());
            medFactorsDegree.addValue(sample.a.degreeSum());
            medFactorsDegree.addValue(sample.b.degreeSum());
            medGCDDegree.addValue(sample.gcd.degreeSum());

            medFactorsSize.addValue(sample.a.size());
            medFactorsSize.addValue(sample.b.size());
            medGCDSize.addValue(sample.gcd.size());

            medFactorsUsedVariables.addValue(sample.a.nUsedVariables());
            medFactorsUsedVariables.addValue(sample.b.nUsedVariables());
            medGCDUsedVariables.addValue(sample.gcd.nUsedVariables());

            medFactorsSparsity.addValue(sample.a.sparsity());
            medFactorsSparsity.addValue(sample.b.sparsity());
            medGCDSparsity.addValue(sample.gcd.sparsity());

            medFactorsSparsity2.addValue(sample.a.sparsity2());
            medFactorsSparsity2.addValue(sample.b.sparsity2());
            medGCDSparsity2.addValue(sample.gcd.sparsity2());

            return sample;
        }

        abstract GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic);

        static double med(DescriptiveStatistics stat) {
            double p = stat.getPercentile(0.5);
            if (p <= 1.1)
                return stat.getMean();
            return p;
        }

        final void printSamplesStatistics() {
            System.out.println("Median polys size : " + med(medFactorsSize));
            System.out.println("Median gcd   size : " + med(medGCDSize));
            System.out.println("Median polys deg  : " + med(medFactorsDegree));
            System.out.println("Median gcd   deg  : " + med(medGCDDegree));
            System.out.println("Median poly nVars : " + med(medFactorsUsedVariables));
            System.out.println("Median gcd  nVars : " + med(medGCDUsedVariables));
            System.out.println("Median polys sparsity  : " + med(medFactorsSparsity));
            System.out.println("Median gcd sparsity    : " + med(medGCDSparsity));
            System.out.println("Median polys sparsity2 : " + med(medFactorsSparsity2));
            System.out.println("Median gcd sparsity2   : " + med(medGCDSparsity2));

        }
    }

    /** provides random data for GCD tests */
    public static abstract class AGCDSampleData<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
            extends GCDSampleData<Term, Poly> {
        public final int nVarsMin, nVarsMax,
                minDegree, maxDegree,
                minSize, maxSize;
        public final RandomGenerator rnd;
        public final RandomDataGenerator rndd;

        public AGCDSampleData(int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            this.nVarsMin = nVarsMin;
            this.nVarsMax = nVarsMax;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.minSize = minSize;
            this.maxSize = maxSize;
            this.rnd = rnd;
            this.rndd = new RandomDataGenerator(rnd);
        }

        /** number of generated samples */
        public long counter = 0;

        /**
         * Next sample
         *
         * @param primitive whether to make sample primitive
         * @param monic     ensure that sample polynomials are monic in main variable
         */
        @Override
        public final GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic) {
            rnd.setSeed(++counter);
            return nextSample1(primitive, monic);
        }

        protected abstract GCDSample<Term, Poly> nextSample1(boolean primitive, boolean monic);

        private static String range(int from, int to) {
            return "[" + from + ", " + to + "]";
        }

        @Override
        public String toString() {
            return "Sample data: " +
                    " nVariables ∈ " + range(nVarsMin, nVarsMax) +
                    ", deg ∈ " + range(minDegree, maxDegree) +
                    ", size ∈ " + range(minSize, maxSize) +
                    ", monic input = " + monic +
                    ", primitive input = " + primitive;
        }
    }

    /** sample data for machine-sized Zp polynomials */
    public static final class lGCDSampleDataZp extends AGCDSampleData<MonomialZp64, MultivariatePolynomialZp64> {
        public lGCDSampleDataZp(int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            super(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        }

        public int minModulusBits = 24, maxModulusBits = 32;

        @Override
        public GCDSample<MonomialZp64, MultivariatePolynomialZp64> nextSample1(boolean primitive, boolean monic) {
            long modulus = APolynomialTest.getModulusRandom(rndd.nextInt(minModulusBits, maxModulusBits));
            IntegersZp64 domain = new IntegersZp64(modulus);

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            MultivariatePolynomialZp64[] data = new MultivariatePolynomialZp64[3];
            for (int i = 0; i < 3; i++) {
                MultivariatePolynomialZp64
                        p = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), domain, MonomialOrder.DEFAULT, rnd);

                // make coefficients small
                p = p.setRing(66).setRing(domain);

                if (primitive)
                    p = MultivariatePolynomialZp64.asNormalMultivariate(p.asOverUnivariateEliminate(0).primitivePart(), 0);

                if (monic)
                    // make monic in main variable by adding some monomial
                    p.add(p.createMonomial(0, p.degreeSum() + 1));

                if (p.isZero()) {
                    --i;
                    continue;
                }
                data[i] = p;
            }

            return new GCDSample<>(data[0], data[1], data[2]);
        }
    }

    /** generic sample data */
    public static final class GCDSampleDataGeneric<E> extends AGCDSampleData<Monomial<E>, MultivariatePolynomial<E>> {
        final Ring<E> ring;

        public GCDSampleDataGeneric(Ring<E> ring, int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            super(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
            this.ring = ring;
            this.rndCoefficients = ring::randomElement;
        }

        public Function<RandomGenerator, E> rndCoefficients;

        @Override
        public GCDSample<Monomial<E>, MultivariatePolynomial<E>> nextSample1(boolean primitive, boolean monic) {
            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            @SuppressWarnings("unchecked")
            MultivariatePolynomial<E>[] data = new MultivariatePolynomial[3];
            for (int i = 0; i < 3; i++) {
                MultivariatePolynomial<E>
                        p = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), ring, MonomialOrder.DEFAULT, rndCoefficients, rnd);

                if (primitive)
                    p = MultivariatePolynomial.asNormalMultivariate(p.asOverUnivariateEliminate(0).primitivePart(), 0);

                if (monic)
                    // make monic in main variable by adding some monomial
                    p.add(p.createMonomial(0, p.degreeSum() + 1));

                if (p.isZero()) {
                    --i;
                    continue;
                }
                data[i] = p;
            }

            return new GCDSample<>(data[0], data[1], data[2]);
        }
    }

    /** GCD algorithm with name abbreviation */
    public static final class GCDAlgorithm<Poly extends AMultivariatePolynomial> {
        final String algorithmName;
        final BiFunction<Poly, Poly, Poly> algorithm;

        public GCDAlgorithm(String algorithmName, BiFunction<Poly, Poly, Poly> algorithm) {
            this.algorithmName = algorithmName;
            this.algorithm = algorithm;
        }

        @Override
        public String toString() {
            return algorithmName;
        }

        public static <P extends AMultivariatePolynomial> GCDAlgorithm<P>
        named(String name, BiFunction<P, P, P> algorithm) {
            return new GCDAlgorithm<>(name, algorithm);
        }
    }

    /** run specified gcd algorithms on sample input */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void testGCDAlgorithms(GCDSampleData<Term, Poly> sampleData, int nIterations,
                           GCDAlgorithm<Poly>... algorithms) {
        System.out.println("\nRunning gcd tests for " + Arrays.toString(algorithms));
        System.out.println(sampleData);

        DescriptiveStatistics timings[] = new DescriptiveStatistics[algorithms.length];
        for (int i = 0; i < timings.length; i++)
            timings[i] = new DescriptiveStatistics();


        int progressPercent = 0, previousPercent = -1;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.stream(timings).forEach(DescriptiveStatistics::clear);

            progressPercent = (int) (100.0 * n / nIterations);
            if (progressPercent != previousPercent) {
                previousPercent = progressPercent;
                // print progress
                System.out.print(">");
                System.out.flush();
            }

            GCDSample<Term, Poly> sample = sampleData.nextSample();

            Poly result = null;
            GCDAlgorithm<Poly> algorithm = null;
            long seed = -1;
            try {
                for (int i = 0; i < algorithms.length; i++) {
                    PrivateRandom.getRandom().setSeed(seed = PrivateRandom.getRandom().nextLong());
                    algorithm = algorithms[i];
                    long start = System.nanoTime();
                    result = algorithm.algorithm.apply(sample.a, sample.b);
                    timings[i].addValue(System.nanoTime() - start);
                    checkConsistency(result);
                    assertTrue(dividesQ(result, sample.gcd));
                }
            } catch (Throwable e) {
                printError(algorithm, result, sample, seed);
                throw e;
            }
        }
        System.out.println();
        sampleData.printSamplesStatistics();
        System.out.println("---------- timings ----------");
        for (int i = 0; i < timings.length; i++)
            System.out.println(algorithms[i].algorithmName + ": " + TimeUnits.statisticsNanotime(timings[i]));
    }

    /** run specified gcd algorithms on sample input */
    public static void testGCDAlgorithm(GCDSampleData<MonomialZp64, MultivariatePolynomialZp64> sampleData, int nIterations,
                                        GCDAlgorithm<MultivariatePolynomialZp64> lAlgorithm,
                                        GCDAlgorithm<MultivariatePolynomial<BigInteger>> algorithm) {
        System.out.println("\nRunning gcd tests for " + algorithm);
        System.out.println(sampleData);

        DescriptiveStatistics
                timingLong = new DescriptiveStatistics(),
                timingBigInteger = new DescriptiveStatistics();

        int progressPercent, previousPercent = -1;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10) {
                timingLong.clear();
                timingBigInteger.clear();
            }

            progressPercent = (int) (100.0 * n / nIterations);
            if (progressPercent != previousPercent) {
                previousPercent = progressPercent;
                // print progress
                System.out.print(">");
                System.out.flush();
            }

            GCDSample<MonomialZp64, MultivariatePolynomialZp64> sample = sampleData.nextSample();

            MultivariatePolynomial<BigInteger> result = null;
            MultivariatePolynomialZp64 lResult = null;
            long seed = -1;
            try {
                PrivateRandom.getRandom().setSeed(seed = PrivateRandom.getRandom().nextLong());

                MultivariatePolynomial<BigInteger>
                        aBig = sample.a.toBigPoly(),
                        bBig = sample.b.toBigPoly();

                long start = System.nanoTime();
                result = algorithm.algorithm.apply(aBig, bBig);
                timingBigInteger.addValue(System.nanoTime() - start);
                checkConsistency(result);
                assertTrue(dividesQ(result, sample.gcd.toBigPoly()));


                start = System.nanoTime();
                lResult = lAlgorithm.algorithm.apply(sample.a, sample.b);
                timingLong.addValue(System.nanoTime() - start);
                checkConsistency(lResult);
                assertTrue(dividesQ(lResult, sample.gcd));

                assertEquals(result.size(), lResult.size());
            } catch (Throwable e) {
                printError(algorithm, result, sample, seed);
                throw e;
            }
        }
        System.out.println();
        sampleData.printSamplesStatistics();
        System.out.println("---------- timings ----------");
        System.out.println("Machine integers  : " + TimeUnits.statisticsNanotime(timingLong));
        System.out.println("Big integers      : " + TimeUnits.statisticsNanotime(timingBigInteger));
    }

    private static void printError(GCDAlgorithm algorithm, Object result, GCDSample sample, long seed) {
        System.out.println();
        System.out.println("========== Exception ========");
        System.out.println("algorithm  : " + algorithm.algorithmName);
        System.out.println("rnd seed   : " + seed);
        System.out.println("ring     : " + sample.domain);
        System.out.println("a          : " + sample.a);
        System.out.println("b          : " + sample.b);
        System.out.println("real gcd   : " + sample.gcd);
        System.out.println("actual gcd : " + result);
        System.out.println("=============================");
    }

    @SuppressWarnings("unchecked")
    private static <Poly extends AMultivariatePolynomial> void checkConsistency(Poly... polys) {
        Arrays.stream(polys).forEach(MultivariateGCDTest::checkConsistency);
    }

    @SuppressWarnings("unchecked")
    private static <Poly extends AMultivariatePolynomial> void checkConsistency(Poly poly) {
        if (poly instanceof MultivariatePolynomial) {
            checkConsistency((MultivariatePolynomial) poly);
        } else
            checkConsistency((MultivariatePolynomialZp64) poly);
    }

    private static <E> void checkConsistency(MultivariatePolynomial<E> poly) {
        Ring<E> ring = poly.ring;
        for (Monomial<E> e : poly.terms) {
            E value = e.coefficient;
            assertFalse(ring.isZero(value));
            assertTrue(value.equals(ring.valueOf(value)));
            if (ring instanceof IntegersZp) {
                assertTrue(ring.signum(value) > 0);
                assertTrue(((BigInteger) value).compareTo(((IntegersZp) ring).modulus) <= 0);
            }
        }
    }

    private static void checkConsistency(MultivariatePolynomialZp64 poly) {
        IntegersZp64 domain = poly.ring;
        for (MonomialZp64 e : poly.terms) {
            long value = e.coefficient;
            assertFalse(value == 0);
            assertTrue(value == domain.modulus(value));
        }
    }

    static <E> Monomial<E> createMonomial(E cf, int... exps) {
        return new Monomial<>(exps, ArraysUtil.sum(exps), cf);
    }
}