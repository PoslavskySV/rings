package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.Rational;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.*;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.RandomUtil;
import cc.r2.core.util.TimeUnits;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.BitSet;
import java.util.function.BiFunction;

import static cc.r2.core.poly.Integers.Integers;
import static cc.r2.core.poly.Rationals.Rationals;
import static cc.r2.core.poly.multivar.DegreeVector.LEX;
import static cc.r2.core.poly.multivar.MultivariateGCD.*;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.asLongPolyZp;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;
import static org.junit.Assert.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateGCDTest extends AbstractPolynomialTest {


    private static void assertBrownGCD(MultivariatePolynomial<BigInteger> gcd,
                                       MultivariatePolynomial<BigInteger> a,
                                       MultivariatePolynomial<BigInteger> b) {
        MultivariatePolynomial<BigInteger> actualGCD = BrownGCD(a, b);
        lMultivariatePolynomialZp lActualGCD = BrownGCD(asLongPolyZp(a), asLongPolyZp(b));
        Assert.assertTrue(dividesQ(actualGCD, gcd));
        Assert.assertEquals(asLongPolyZp(actualGCD).monic(), lActualGCD.monic());
    }

    @Test
    public void testBrown1() throws Exception {
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(1321323));
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c", domain, LEX),
                b = parse("a^2+2*b^2 + 2*c", domain, LEX),
                gcd = parse("c*a+b+a+ c*a^3", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown2() throws Exception {
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(1321323));
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c", domain, LEX),
                b = parse("a^2+2*b^2 + 2*c", domain, LEX),
                gcd = parse("c*a*b+b*b+a*b+ c*a^3*b", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown3() throws Exception {
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(659));
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b+655*a*b^3*c^2+6*a^2*b^3*c^2", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown3a() throws Exception {
        IntegersModulo domain = new IntegersModulo(17);
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b^6+655*a*b^3*c^6+6*a^2*b^3*c^4", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown4() throws Exception {
        IntegersModulo domain = new IntegersModulo(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*b^5*c+a*b^3+a^2*b^2+a^2*b^2*c+a^3*b*c^3", domain, LEX, vars),
                b = parse("9*a*b^2*c^6+a*b^4*c^6+a^2*b^2*c^3+a^5*b^2+a^5*b^6*c^4+a^6*c+a^6*b^2*c", domain, LEX, vars),
                gcd = parse("653*b^3*c^4+b^4+b^5*c^3+a^2*b*c^2+a^4*b^2*c^4", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown5() throws Exception {
        RandomGenerator rnd = PrivateRandom.getRandom();
        rnd.setSeed(28);
        IntegersModulo domain = new IntegersModulo(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("561*a^2*c^2+a^2*b^2*c+a^3*b+a^4*b^2+a^4*b^5*c^3+a^5*b", domain, LEX, vars),
                b = parse("561*a*c^3+a*b^4*c^5+a^2*c^2+a^2*b^6*c^3+a^3*b^6*c^5+a^5*b^5*c^3+a^5*b^5*c^6", domain, LEX, vars),
                gcd = parse("4*c^2+b^4+a^2*b^4*c+a^3*b^2*c+a^3*b^6+a^5*b^2*c^6+a^6*b^5", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown6() throws Exception {
        PrivateRandom.getRandom().setSeed(1564);
        IntegersModulo domain = new IntegersModulo(937);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("931*a^3*b^4*c+a^4+a^4*b^6*c^2+a^5*b*c^3+a^6*b*c^2", domain, LEX, vars),
                b = parse("932*b*c+a*b^6*c^2+a^3*b*c^2+a^3*b^3*c^5+a^3*b^5*c+a^5*b^5*c^3+a^6*b^2*c^6+a^6*b^4*c^5+a^6*b^6", domain, LEX, vars),
                gcd = parse("935*c^2+c^4+a^3*b*c^5+a^3*b^2*c^3+a^4*b^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown7() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersModulo domain = new IntegersModulo(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("563*b^2*c+6*b^4*c^3+4*a*b^4*c^3+563*a*b^4*c^4+560*a^2*b^5*c^2+9*a^3*b^4*c+5*a^4*c^2+7*a^4*b^3*c^5+4*a^5*b^4*c^5+6*a^5*b^5", domain, LEX, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, LEX, vars),
                gcd = parse("4+8*b^2*c^3+4*b^3+8*b^3*c+7*a*c+a*b*c+7*a^2*b^2+2*a^2*b^2*c^2+5*a^3*c^2+5*a^3*c^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
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

    @Test
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

    @Test
    public void testBrown_sparse_variables_random3() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD));
    }

    @Test
    public void testBrown8() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersModulo domain = new IntegersModulo(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, LEX, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, LEX, vars),
                gcd = parse("a^2 + c^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown9() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        IntegersModulo domain = new IntegersModulo(569);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, LEX, vars),
                b = parse("b^2 + d^2", domain, LEX, vars),
                gcd = parse("d^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertBrownGCD(gcd, a, b);
    }

    @Test
    public void testBrown10() throws Exception {
        IntegersModulo domain = new IntegersModulo(5642359);
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger>
                a = parse("1199884 + 4783764*b + b^2 + 3215597*a*b + 2196297*a*b^2 + 4781733*a^4 + a^4*b + 2196297*a^5*b", domain, LEX, vars),
                b = parse("4645946 + 3921107*b + b^2 + 5605437*a*b + 2196297*a*b^2 + 4781733*a^3 + a^3*b + 2196297*a^4*b", domain, LEX, vars);

        MultivariatePolynomial<BigInteger> gcdActual = BrownGCD(a, b);
        gcdActual = gcdActual.monic().multiply(domain.valueOf(4781733));
        MultivariatePolynomial<BigInteger> expected = parse("1574588 + 4559668*b + 4781733*a*b", domain, LEX, vars);
        assertEquals(expected, gcdActual);
    }


    private static void assertZippelGCD(MultivariatePolynomial<BigInteger> gcd,
                                        MultivariatePolynomial<BigInteger> a,
                                        MultivariatePolynomial<BigInteger> b) {
        MultivariatePolynomial<BigInteger> actualGCD = ZippelGCD(a, b);
        lMultivariatePolynomialZp lActualGCD = ZippelGCD(asLongPolyZp(a), asLongPolyZp(b));
        Assert.assertTrue(dividesQ(actualGCD, gcd));
        Assert.assertEquals(asLongPolyZp(actualGCD).monic(), lActualGCD.monic());
    }

    @Test
    public void testZippel1() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^2 + b^2 + a*c^2", domain, LEX, vars),
                b = parse("a^2 + 2*b^2 + b*c^2", domain, LEX, vars),
                gcd = parse("a^2 + c*a + b + a*c*b + a*c^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        RandomGenerator rnd = getRandom();

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(101);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        for (int i = 0; i < 100; i++) {
            SparseInterpolation<BigInteger> sparseInterpolation
                    = createInterpolation(variable, a, b, skeleton, rnd);
            BigInteger point = domain.randomElement(rnd);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }

        lMultivariatePolynomialZp la = asLongPolyZp(a), lb = asLongPolyZp(b),
                lskeleton = asLongPolyZp(skeleton), lgcd = asLongPolyZp(gcd);
        for (int i = 0; i < 100; i++) {
            lSparseInterpolation sparseInterpolation
                    = createInterpolation(variable, la, lb, lskeleton, rnd);
            long point = domain.asLong().randomElement(rnd);
            assertEquals(lgcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test
    public void testZippel2() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^2 + b^2 + a*c^2", domain, LEX, vars),
                b = parse("a^2 + 2*b^2 + b*c^2", domain, LEX, vars),
                gcd = parse("a^2 + c*a + b + a*c*b + a*c^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertZippelGCD(gcd, a, b);
    }

    @Test
    public void testZippel5() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(31579447));
        System.out.println(domain);
        MultivariatePolynomial<BigInteger>
                a = parse("b^2 + a*c + a*b*c^2 + a*b^2 + a^5", domain, LEX, vars),
                b = parse("1 + a*b", domain, LEX, vars),
                gcd = parse("1 + a*b*c^2 + a^2*c^2 + a^2*b*c^2 + a^2*b^2*c^2 + a^7", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
        BigInteger point = domain.valueOf(1324);
        assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
    }

    @Test
    public void testZippel6() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(31579447));
        MultivariatePolynomial<BigInteger>
                a = parse("1+3*b*c^2+7*b^2*c^2+4*a+7*a^3*c^2+a^6", domain, LEX, vars),
                b = parse("b^3*c^2+a*b^2*c+9*a^2*b*c+9*a^3*b*c^2+a^7", domain, LEX, vars),
                gcd = parse("b^3*c^2+2*a*c+5*a*b+2*a*b^2*c+4*a*b^3*c^2+5*a^3*b+a^7", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
        BigInteger point = domain.valueOf(1324);
        assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
    }

    @Test
    public void testZippel_monic_random1() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 10, minDegree = 3, maxDegree = 5, minSize = 5, maxSize = 10;

        int nIterations = its(100, 1500);
        GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp> sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n % 100 == 0) System.out.println(n);

            GCDSample<lMonomialTerm, lMultivariatePolynomialZp> data = sampleData.nextSample(true, true);

            int variable = data.a.nVariables - 1;
            long seed;
            do {seed = data.a.domain.randomElement(rnd);} while (seed == 0);
            lMultivariatePolynomialZp skeleton = data.gcd.evaluate(variable, seed);

            for (int i = 0; i < 10; i++) {
                int rndSeed = i^n;
                rnd.setSeed(rndSeed);
                SparseInterpolation<BigInteger> sparseInterpolation
                        = createInterpolation(variable, data.a.toBigPoly(), data.b.toBigPoly(), skeleton.toBigPoly(), rnd);

                lSparseInterpolation lSparseInterpolation
                        = createInterpolation(variable, data.a, data.b, skeleton, rnd);
                long point = data.a.domain.randomElement(rnd);
                try {
                    lMultivariatePolynomialZp expected = data.gcd.evaluate(variable, point).monic();

                    lMultivariatePolynomialZp lActual = lSparseInterpolation.evaluate(point).monic();
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

    @Test
    public void testZippel7() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(31579447));
        MultivariatePolynomial<BigInteger>
                a = parse("29322275*b+5*b*c+6*a^2*b^3*c^2+29322274*a^2*b^3*c^2*d^3+5*a^3*b*c^2*d^2+a^11", domain, LEX, vars),
                b = parse("7*a^3*b^3*c^3*d^4+9*a^3*b^4*c+29322274*a^3*b^4*c^5*d^2+29322277*a^5*b*c*d+a^5*b^4+a^15", domain, LEX, vars),
                gcd = parse("4*d^2+8*c^2*d+4*b*c+6*b^3*c^2*d^2+2*a^3*b^2*c^3*d+a^10", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
        BigInteger point = domain.valueOf(1324);
        assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
    }

    @Test
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

    @Test(timeout = 10000)
    public void testZippel9() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(26478253);
        PrivateRandom.getRandom().setSeed(0);
        MultivariatePolynomial<BigInteger>
                a = parse("26478246*a*c^2+7*a*b+26478250*a*b*c^2+26478249*a*b^3*c^2+26478248*a^2*c^2+8*a^3*b*c^2+a^7", domain, LEX, vars),
                b = parse("4*b^3*c^2+7*a+5*a*b+8*a*b^2+6*a^3*b^2*c+a^7", domain, LEX, vars),
                gcd = parse("26478248*a*b^2*c^2+3*a*b^3*c^2+2*a^2*b^3*c^2+5*a^3*c+a^8", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertNotNull(ZippelGCD(a, b));
        assertNotNull(ZippelGCD(asLongPolyZp(a), asLongPolyZp(b)));
    }

    @Test
    public void testZippel8() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(31579447);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a+29923129*a*b*c^2+3*a*b^2+29923132*a^2*b*c^2+7*a^3*c", domain, LEX, vars),
                b = parse("4*c^2+29923126*a*c^2+5*a*b+6*a^2*b^2*c^3+29923128*a^3*c^3", domain, LEX, vars),
                gcd = parse("29923132+8*b*c^3+29923132*b^2*c^3+8*a*b*c+7*a^3*b^3*c", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int variable = a.nVariables - 1;
        BigInteger seed = BigInteger.valueOf(7893482);
        MultivariatePolynomial<BigInteger> skeleton = gcd.evaluate(variable, seed);

        SparseInterpolation<BigInteger> sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
        BigInteger point = domain.valueOf(1324);
        assertEquals(gcd.evaluate(variable, point).monic(), sparseInterpolation.evaluate(point).monic());
    }

    @Test
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

    @Test
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

    @Test
    public void testZippel3() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("5*a^2*c^2+5*a^2*b^2*c^2+5*a^2*b^4*c^3+9*a^2*b^5*c^5+25709547*a^3*b^6*c^6+8*a^4*b*c^3+a^4*b^3*c+5*a^4*b^3*c^6", domain, LEX, vars),
                b = parse("3*a*b^2*c^2+2*a^2*b^4+25709540*a^4*b*c^6+7*a^5*c^2+8*a^6*b*c^3", domain, LEX, vars),
                gcd = parse("a + 5*b^2*c^6+2*a^4*b^4*c^5+25709543*a^5*b^2*c^5+9*a^6*c+25709540*a^6*c^3", domain, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        assertZippelGCD(gcd, a, b);
    }

    @Test
    public void testZippel10() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a + b + c", domain, LEX, vars),
                b = parse("a - b + c", domain, LEX, vars),
                gcd = parse("a^3*b+2*a^3+c*a^3+12*b^2+24*b+12*b*c", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        assertZippelGCD(gcd, a, b);
    }

    @Ignore
    @Test
    public void testZippel4_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b", "c", "d", "e"};
        IntegersModulo domain = new IntegersModulo(SmallPrimes.nextPrime(100000));
        MultivariatePolynomial<BigInteger>
                a = parse("2147483167*a^4*b^60*c^57*d^26*e+44*a^8*b^39*c^67*d^22*e^17+38*a^32*b^6*c^13*d^10*e^3+357*a^36*b^34*c^60*d^2*e^59+563*a^42*b^41*c^45*d^52*e^14+257*a^44*b^68*c^43*d^2*e^73+613*a^48*b^55*c^22*d^32*e^19+2147483093*a^52*b^26*c^4*d^72*e^32+19*a^52*b^40*c^26*d^45*e^55+639*a^55*b^72*c^55*d^65", domain, LEX, vars),
                b = parse("2147483150*b^25*c^18*d^62*e^59+2147482723*a^4*b^5*c^65*d^26*e^7+261*a^15*b^60*c^59*d^63*e^53+394*a^27*b^22*c^34*d^54*e^13+952*a^39*b^48*c^17*d^54*e^16+243*a^60*b^15*c^3*d^51*e^46+40*a^61*b^56*c^39*d^40*e^21+555*a^62*b^20*c^20*d^60*e^47+627*a^67*b^8*c^22*d^67*e^61+447*a^70*b^59*c^71*d^24*e^5", domain, LEX, vars),
                gcd = parse("35*a*b^36*c^74*d^62*e^51+376*a^2*b^28*c^64*e^53+893*a^6*b^13*c^60*d^44*e^42+23*a^8*b^71*c^40*d^36*e^11+783*a^20*b^28*c^12*d^31*e^68+2147482938*a^31*b^30*c^40*d^65*e^72+2147482960*a^31*b^49*c^38*d^71*e^55+737*a^47*b^15*c^71*d^13*e^72+868*a^53*b^30*c^40*d^29*e^46+898*a^61*b^71*c^13*d^50*e^66", domain, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        System.out.println(a);
        System.out.println(b);

        lMultivariatePolynomialZp
                aL = asLongPolyZp(a),
                bL = asLongPolyZp(b);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertEquals(10, ZippelGCD(aL, bL).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            start = System.nanoTime();
            System.out.println(ZippelGCD(aL.clone().increment(), bL));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();
//            System.out.println(TimeUnits.nanosecondsToString(MultivariateGCD.BROWN));
        }
    }

    @Ignore
    @Test
    public void testZippel5_bivariate_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b"};
        IntegersModulo domain = new IntegersModulo(100011111111101L);
//        IntegersModulo domain = new IntegersModulo(100019);
        MultivariatePolynomial<BigInteger>
                a = parse("38*a^32*b^6 + 2147483093*a^52*b^26 + 357*a^36*b^34 + 44*a^8*b^39 + 19*a^52*b^40 + 563*a^42*b^41 + 613*a^48*b^55 + 2147483167*a^4*b^60 + 257*a^44*b^68 + 639*a^55*b^72", domain, LEX, vars),
                b = parse("2147482723*a^4*b^5 + 627*a^67*b^8 + 243*a^60*b^15 + 555*a^62*b^20 + 394*a^27*b^22 + 2147483150*b^25 + 952*a^39*b^48 + 40*a^61*b^56 + 447*a^70*b^59 + 261*a^15*b^60", domain, LEX, vars),
                gcd = parse("893*a^6*b^13 + 737*a^47*b^15 + 376*a^2*b^28 + 783*a^20*b^28 + 2147482938*a^31*b^30 + 868*a^53*b^30 + 35*a*b^36 + 2147482960*a^31*b^49 + 23*a^8*b^71 + 898*a^61*b^71", domain, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        System.out.println(a);
        System.out.println(b);

        lMultivariatePolynomialZp
                aL = asLongPolyZp(a),
                bL = asLongPolyZp(b);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertEquals(10, ZippelGCD(aL, bL).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            start = System.nanoTime();
            System.out.println(ZippelGCD(aL.clone().increment(), bL));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();
        }
    }

    @Test
    public void testZippel_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD));
    }


    @Test
    public void testPairedIterator() throws Exception {
        RandomGenerator rnd = getRandom();
        int nIterations = its(1000, 1000);
        for (int n = 0; n < nIterations; n++) {
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(5, 50, 20, rnd),
                    b = randomPolynomial(5, 50, 20, rnd);

            PairIterator<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>>
                    it = new PairIterator<>(a, b);

            MultivariatePolynomial<BigInteger> acc = a.createZero();
            while (it.hasNext()) {
                it.advance();
                acc.add(it.aTerm);
                acc.add(it.bTerm);
                assertTrue(it.aTerm.coefficient.isZero() || it.bTerm.coefficient.isZero() || 0 == a.ordering.compare(it.aTerm, it.bTerm));
            }
            assertEquals(acc, a.clone().add(b));
        }
    }

    @Test
    public void testSparseInterpolation1() throws Exception {
        IntegersModulo domain = new IntegersModulo(31574773);
        MultivariatePolynomial<BigInteger>
                a = parse("31574768*a*b^2*c^4+4*a^4*b^3+3*a^5*b+31574764*a^5*b^5*c^5+6*a^6*b^3*c^2", domain, LEX),
                b = parse("7*a^2*b^6*c^3+a^5*b^4*c^4+31574764*a^6*c^3+5*a^6*b^2*c^2", domain, LEX),
                gcd = parse("9*c^4+31574766*a*b^2+2*a^2*b*c^2+31574768*a^2*b^3*c^6+9*a^3*b^2*c^3", domain, LEX);

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), interpolateGCD(asLongPolyZp(a), asLongPolyZp(b), asLongPolyZp(ZippelGCD(a, b)), getRandom()).monic().toBigPoly());
    }

    @Test
    public void testSparseInterpolation2() throws Exception {
        IntegersModulo domain = new IntegersModulo(24001871);
        MultivariatePolynomial<BigInteger>
                a = parse("3*b^4*c^2+7*b^4*c^3+4*a*b^5*c^3+6*a^4*b^6+24001865*a^5*b", domain, LEX, "a", "b", "c"),
                b = parse("5*a*c^4+9*a^4*b^4*c^2+9*a^6*b*c^6", domain, LEX, "a", "b", "c"),
                gcd = parse("5*a*b^2*c^2+a*b^2*c^4+24001866*a*b^4*c^3", domain, LEX, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(asLongPolyZp(gcd).monic(), interpolateGCD(asLongPolyZp(a), asLongPolyZp(b), asLongPolyZp(gcd), getRandom()).monic());
    }

    @Test
    public void testSparseInterpolation3() throws Exception {
        IntegersModulo domain = new IntegersModulo(17312587);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a^3*c^6+9*a^5*b^2*c^3+7*a^5*b^6*c^5+8*a^5*b^6*c^6+6*a^6*b^6*c", domain, LEX, "a", "b", "c"),
                b = parse("17312581*c^6+5*a^2*b^2*c^6+3*a^4*b^6*c^4+2*a^5+4*a^5*b^3*c^6", domain, LEX, "a", "b", "c"),
                gcd = parse("5*a^5*b*c^2+6*a^5*b^3*c^6+2*a^5*b^4*c^4", domain, LEX, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);
        lMultivariatePolynomialZp intrp = interpolateGCD(asLongPolyZp(a), asLongPolyZp(b), asLongPolyZp(gcd), getRandom());
        assertEquals(asLongPolyZp(gcd).monic(), intrp.monic());
    }

    @Test
    public void testSparseInterpolation4() throws Exception {
        IntegersModulo domain = new IntegersModulo(27445993);
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b*c^3+8*a^3*c+8*a^4*b^2*c^4+8*a^4*b^6*c^6", domain, LEX, "a", "b", "c"),
                b = parse("a*b^6*c^2+6*a^2*b^3*c^3+27445990*a^3*b^6*c^2", domain, LEX, "a", "b", "c"),
                gcd = parse("5*b*c^3+8*b^5*c+4*b^6+5*a*b^3+4*a^6*b^3*c^3", domain, LEX, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        lMultivariatePolynomialZp la = asLongPolyZp(a), lb = asLongPolyZp(b);
        lMultivariatePolynomialZp lgcd = ZippelGCD(la, lb);
        lMultivariatePolynomialZp intrp = interpolateGCD(la, lb, lgcd, getRandom());
        assertEquals(lgcd.monic(), intrp.monic());
    }

    @Test
    public void testSparseInterpolation_random1() throws Exception {
        int nIterations = its(1000, 5000);
        RandomGenerator rnd = getRandom();

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 5, 5, 15, 5, 15, rnd);
        for (int n = 0; n < nIterations; n++) {
            GCDSample<lMonomialTerm, lMultivariatePolynomialZp> gcdTriplet = sampleData.nextSample(false, false);
            lMultivariatePolynomialZp gcd = null, actual = null;
            try {
                lMultivariatePolynomialZp la = gcdTriplet.a;
                lMultivariatePolynomialZp lb = gcdTriplet.b;
                gcd = ZippelGCD(la, lb);
                if (la.isConstant() || lb.isConstant() || gcd.degree(0) == 0) {
                    --n; continue;
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

    @Test
    public void testSparseInterpolation5() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b*c^3+8*a^3*c+8*a^4*b^2*c^4+8*a^4*b^6*c^6", Integers, LEX, "a", "b", "c"),
                b = parse("a*b^6*c^2+6*a^2*b^3*c^3+27445990*a^3*b^6*c^2", Integers, LEX, "a", "b", "c"),
                gcd = parse("5*b*c^3+8*b^5*c+4*b^6+5*a*b^3+4*a^6*b^3*c^3", Integers, LEX, "a", "b", "c");

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        IntegersModulo domain = new IntegersModulo(27445993);

        lMultivariatePolynomialZp
                la = asLongPolyZp(a.setDomain(domain)),
                lb = asLongPolyZp(b.setDomain(domain));
        lMultivariatePolynomialZp skeleton = ZippelGCD(la, lb);

        IntegersModulo domain1 = new IntegersModulo(BigPrimes.nextPrime(37445993132451L));
        lMultivariatePolynomialZp
                la1 = asLongPolyZp(a.setDomain(domain1)),
                lb1 = asLongPolyZp(b.setDomain(domain1));

        lMultivariatePolynomialZp gcd1 = ZippelGCD(la1, lb1);

        skeleton = skeleton.setDomain(la1.domain);
        lMultivariatePolynomialZp intrp = interpolateGCD(la1, lb1, skeleton, getRandom());
        assertEquals(gcd1.monic(), intrp.monic());
    }

    @Test
    public void testSparseInterpolation_random2() throws Exception {
        int nIterations = its(500, 1000);
        int badEvaluations = 0;
        RandomGenerator rnd = getRandom();
        GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> sampleData =
                new GCDSampleDataGeneric<>(Integers, 3, 5, 5, 15, 5, 15, rnd);
        for (int n = 0; n < nIterations; n++) {
            GCDSample<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> gcdTriplet = sampleData.nextSample(false, false);
            lMultivariatePolynomialZp skeleton = null, gcd = null, actual = null;
            IntegersModulo domain = null, domain1 = null;
            long seed = -1;
            try {

                domain = new IntegersModulo(getModulusRandom(20));
                lMultivariatePolynomialZp
                        la = asLongPolyZp(gcdTriplet.a.setDomain(domain)),
                        lb = asLongPolyZp(gcdTriplet.b.setDomain(domain));

                skeleton = ZippelGCD(la, lb);
                if (la.isConstant() || lb.isConstant() || skeleton.degree(0) == 0) {
                    --n; continue;
                }

                domain1 = new IntegersModulo(getModulusRandom(20));
                lMultivariatePolynomialZp
                        la1 = asLongPolyZp(gcdTriplet.a.setDomain(domain1)),
                        lb1 = asLongPolyZp(gcdTriplet.b.setDomain(domain1));

                gcd = ZippelGCD(la1, lb1);
                if (!gcd.sameSkeleton(skeleton)) {
                    --n; continue;
                }

                rnd.setSeed(seed = rnd.nextLong());
                actual = interpolateGCD(la1, lb1, skeleton.setDomain(la1.domain), rnd);
                if (actual == null) {
                    ++badEvaluations;
                    // bad evaluation point => try over
                    actual = interpolateGCD(la1, lb1, skeleton.setDomain(la1.domain), rnd);
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
            IntegersModulo domain = new IntegersModulo(1049);
            MultivariatePolynomial<BigInteger>
                    a = parse("15*a*b^5*c^5*d^3+27*a^2*b^10*c^4*d+35*a^3*b^7*c^5*d^4+20*a^3*b^7*c^5*d^5+40*a^4*b^11*c^11*d^10+63*a^4*b^12*c^4*d^2+36*a^4*b^12*c^4*d^3+72*a^5*b^16*c^10*d^8+243*a^6*b^12*c^9*d^4+15*a^6*b^15*c*d^3+15*a^6*b^15*c^11*d^10+45*a^7*b^9*c^11*d^7+5*a^7*b^14*c^7*d^4+12*a^8*b*c*d^11+35*a^8*b^14*c^4*d^8+567*a^8*b^14*c^9*d^5+324*a^8*b^14*c^9*d^6+81*a^8*b^14*c^10*d^5+231*a^8*b^15*c^3*d^10+35*a^8*b^17*c*d^4+20*a^8*b^17*c*d^5+35*a^8*b^17*c^11*d^11+20*a^8*b^17*c^11*d^12+9*a^8*b^19*c^6*d^2+15*a^9*c^10*d^11+415*a^9*b^6*c^4*d^11+390*a^9*b^8*c^8*d^5+648*a^9*b^18*c^15*d^11+63*a^9*b^19*c^3*d^6+40*a^9*b^21*c^7*d^10+40*a^9*b^21*c^17*d^17+24*a^10*c^3*d^11+28*a^10*b^3*c*d^12+16*a^10*b^3*c*d^13+747*a^10*b^11*c^3*d^9+20*a^10*b^12*c^4*d^5+702*a^10*b^13*c^7*d^3+405*a^10*b^14*c^5*d^5+539*a^10*b^17*c^3*d^11+308*a^10*b^17*c^3*d^12+35*a^11*b^2*c^10*d^12+20*a^11*b^2*c^10*d^13+32*a^11*b^7*c^7*d^18+36*a^11*b^17*c^3*d^3+729*a^11*b^19*c^4*d^3+616*a^11*b^21*c^9*d^17+56*a^12*b^2*c^3*d^12+32*a^12*b^2*c^3*d^13+40*a^12*b^6*c^16*d^18+729*a^12*b^16*c^15*d^8+45*a^12*b^19*c^7*d^7+45*a^12*b^19*c^17*d^14+81*a^12*b^21*c^11*d^5+5*a^12*b^24*c^3*d^4+5*a^12*b^24*c^13*d^11+18*a^13*b^3*c^5*d^13+64*a^13*b^6*c^9*d^18+567*a^13*b^21*c^8*d^9+35*a^13*b^24*d^8+35*a^13*b^24*c^10*d^15+36*a^14*b^5*c^7*d^15+4*a^14*b^10*c^3*d^12+429*a^14*b^13*c^8*d^12+24*a^14*b^15*c^12*d^6+415*a^14*b^16*d^11+415*a^14*b^16*c^10*d^18+390*a^14*b^18*c^4*d^5+390*a^14*b^18*c^14*d^12+693*a^14*b^19*c^9*d^14+77*a^14*b^24*c^5*d^11+45*a^15*b^4*c^16*d^15+42*a^15*b^5*c^5*d^14+24*a^15*b^5*c^5*d^15+5*a^15*b^9*c^12*d^12+28*a^15*b^10*d^16+324*a^15*b^19*c^8*d^6+267*a^15*b^21*c^9*d^6+20*a^15*b^22*d^5+20*a^15*b^22*c^10*d^12+405*a^15*b^24*c*d^5+539*a^15*b^24*c^2*d^15+405*a^15*b^24*c^11*d^12+332*a^16*b^2*d^19+312*a^16*b^4*c^4*d^13+72*a^16*b^4*c^9*d^15+8*a^16*b^9*c^5*d^12+35*a^16*b^9*c^9*d^16+48*a^16*b^9*c^11*d^20+97*a^16*b^16*c^2*d^18+761*a^16*b^18*c^6*d^12+415*a^17*b*c^9*d^19+390*a^17*b^3*c^13*d^13+16*a^17*b^8*d^13+56*a^17*b^9*c^2*d^16+324*a^17*b^10*c*d^13+308*a^17*b^22*c^2*d^12+992*a^17*b^24*c^3*d^12+664*a^18*b*c^2*d^19+624*a^18*b^3*c^6*d^13+20*a^18*b^7*c^9*d^13+405*a^18*b^9*c^10*d^13+32*a^19*b^7*c^2*d^13+54*a^19*b^7*c^11*d^17+648*a^19*b^9*c^3*d^13+6*a^19*b^12*c^7*d^14+42*a^20*b^12*c^4*d^18+498*a^21*b^4*c^4*d^21+468*a^21*b^6*c^8*d^15+24*a^22*b^10*c^4*d^15+486*a^22*b^12*c^5*d^15", domain, LEX),
                    b = parse("12*c^6*d^2+234*b^8*c^8*d^6+28*a^2*b^2*c^6*d^3+16*a^2*b^2*c^6*d^4+3*a^2*b^7*c^5*d+546*a^2*b^10*c^8*d^7+312*a^2*b^10*c^8*d^8+32*a^3*b^6*c^12*d^9+624*a^3*b^14*c^14*d^13+7*a^4*b^9*c^5*d^2+4*a^4*b^9*c^5*d^3+8*a^5*b^13*c^11*d^8+12*a^6*b*c^6+36*a^6*b^4*c^12*d^6+4*a^6*b^9*c^8*d^3+702*a^6*b^12*c^14*d^10+78*a^6*b^17*c^10*d^7+28*a^7*b^9*c^5*d^7+546*a^7*b^17*c^7*d^11+3*a^8*c*d^3+332*a^8*b*c^5*d^10+28*a^8*b^3*c^6*d+16*a^8*b^3*c^6*d^2+312*a^8*b^3*c^9*d^4+180*a^8*b^9*c^7*d^14+9*a^8*b^11*c^11*d^5+839*a^8*b^11*c^11*d^8+a^8*b^16*c^7*d^2+16*a^9*b^7*c^5*d^4+32*a^9*b^7*c^12*d^7+324*a^9*b^9*c^6*d^4+312*a^9*b^15*c^7*d^8+7*a^9*b^16*c^4*d^6+24*a^9*b^17*c^8*d^8+7*a^10*b^2*c*d^4+4*a^10*b^2*c*d^5+83*a^10*b^8*c^4*d^9+78*a^10*b^10*c^8*d^3+8*a^11*b^6*c^7*d^10+4*a^11*b^14*c^4*d^3+81*a^11*b^16*c^5*d^3+36*a^12*b^5*c^12*d^4+4*a^12*b^10*c^8*d+28*a^13*b^10*c^5*d^5+332*a^14*b^2*c^5*d^8+9*a^14*b^4*c^7*d^7+312*a^14*b^4*c^9*d^2+a^14*b^9*c^3*d^4+16*a^15*b^8*c^5*d^2+7*a^15*b^9*d^8+324*a^15*b^10*c^6*d^2+83*a^16*b*d^11+78*a^16*b^3*c^4*d^5+4*a^17*b^7*d^5+81*a^17*b^9*c*d^5", domain, LEX),
                    base = parse("3*c+7*a^2*b^2*c*d+4*a^2*b^2*c*d^2+8*a^3*b^6*c^7*d^7+9*a^6*b^4*c^7*d^4+a^6*b^9*c^3*d+7*a^7*b^9*d^5+17492158*a^8*b*d^8+17492153*a^8*b^3*c^4*d^2+4*a^9*b^7*d^2+17492156*a^9*b^9*c*d^2", domain, LEX);

            lMultivariatePolynomialZp
                    la = asLongPolyZp(a),
                    lb = asLongPolyZp(b),
                    skeleton = asLongPolyZp(base);
            lMultivariatePolynomialZp lgcd = ZippelGCD(la, lb);
            lMultivariatePolynomialZp intrp = null;
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

    @Test
    public void testSparseInterpolation6() throws Exception {
        RandomGenerator rnd = PrivateRandom.getRandom();
        rnd.setSeed(743);
        IntegersModulo domain = new IntegersModulo(1049);
        MultivariatePolynomial<BigInteger>
                a = parse("15*a*b^5*c^5*d^3+27*a^2*b^10*c^4*d+35*a^3*b^7*c^5*d^4+20*a^3*b^7*c^5*d^5+40*a^4*b^11*c^11*d^10+63*a^4*b^12*c^4*d^2+36*a^4*b^12*c^4*d^3+72*a^5*b^16*c^10*d^8+243*a^6*b^12*c^9*d^4+15*a^6*b^15*c*d^3+15*a^6*b^15*c^11*d^10+45*a^7*b^9*c^11*d^7+5*a^7*b^14*c^7*d^4+12*a^8*b*c*d^11+35*a^8*b^14*c^4*d^8+567*a^8*b^14*c^9*d^5+324*a^8*b^14*c^9*d^6+81*a^8*b^14*c^10*d^5+231*a^8*b^15*c^3*d^10+35*a^8*b^17*c*d^4+20*a^8*b^17*c*d^5+35*a^8*b^17*c^11*d^11+20*a^8*b^17*c^11*d^12+9*a^8*b^19*c^6*d^2+15*a^9*c^10*d^11+415*a^9*b^6*c^4*d^11+390*a^9*b^8*c^8*d^5+648*a^9*b^18*c^15*d^11+63*a^9*b^19*c^3*d^6+40*a^9*b^21*c^7*d^10+40*a^9*b^21*c^17*d^17+24*a^10*c^3*d^11+28*a^10*b^3*c*d^12+16*a^10*b^3*c*d^13+747*a^10*b^11*c^3*d^9+20*a^10*b^12*c^4*d^5+702*a^10*b^13*c^7*d^3+405*a^10*b^14*c^5*d^5+539*a^10*b^17*c^3*d^11+308*a^10*b^17*c^3*d^12+35*a^11*b^2*c^10*d^12+20*a^11*b^2*c^10*d^13+32*a^11*b^7*c^7*d^18+36*a^11*b^17*c^3*d^3+729*a^11*b^19*c^4*d^3+616*a^11*b^21*c^9*d^17+56*a^12*b^2*c^3*d^12+32*a^12*b^2*c^3*d^13+40*a^12*b^6*c^16*d^18+729*a^12*b^16*c^15*d^8+45*a^12*b^19*c^7*d^7+45*a^12*b^19*c^17*d^14+81*a^12*b^21*c^11*d^5+5*a^12*b^24*c^3*d^4+5*a^12*b^24*c^13*d^11+18*a^13*b^3*c^5*d^13+64*a^13*b^6*c^9*d^18+567*a^13*b^21*c^8*d^9+35*a^13*b^24*d^8+35*a^13*b^24*c^10*d^15+36*a^14*b^5*c^7*d^15+4*a^14*b^10*c^3*d^12+429*a^14*b^13*c^8*d^12+24*a^14*b^15*c^12*d^6+415*a^14*b^16*d^11+415*a^14*b^16*c^10*d^18+390*a^14*b^18*c^4*d^5+390*a^14*b^18*c^14*d^12+693*a^14*b^19*c^9*d^14+77*a^14*b^24*c^5*d^11+45*a^15*b^4*c^16*d^15+42*a^15*b^5*c^5*d^14+24*a^15*b^5*c^5*d^15+5*a^15*b^9*c^12*d^12+28*a^15*b^10*d^16+324*a^15*b^19*c^8*d^6+267*a^15*b^21*c^9*d^6+20*a^15*b^22*d^5+20*a^15*b^22*c^10*d^12+405*a^15*b^24*c*d^5+539*a^15*b^24*c^2*d^15+405*a^15*b^24*c^11*d^12+332*a^16*b^2*d^19+312*a^16*b^4*c^4*d^13+72*a^16*b^4*c^9*d^15+8*a^16*b^9*c^5*d^12+35*a^16*b^9*c^9*d^16+48*a^16*b^9*c^11*d^20+97*a^16*b^16*c^2*d^18+761*a^16*b^18*c^6*d^12+415*a^17*b*c^9*d^19+390*a^17*b^3*c^13*d^13+16*a^17*b^8*d^13+56*a^17*b^9*c^2*d^16+324*a^17*b^10*c*d^13+308*a^17*b^22*c^2*d^12+992*a^17*b^24*c^3*d^12+664*a^18*b*c^2*d^19+624*a^18*b^3*c^6*d^13+20*a^18*b^7*c^9*d^13+405*a^18*b^9*c^10*d^13+32*a^19*b^7*c^2*d^13+54*a^19*b^7*c^11*d^17+648*a^19*b^9*c^3*d^13+6*a^19*b^12*c^7*d^14+42*a^20*b^12*c^4*d^18+498*a^21*b^4*c^4*d^21+468*a^21*b^6*c^8*d^15+24*a^22*b^10*c^4*d^15+486*a^22*b^12*c^5*d^15", domain, LEX),
                b = parse("12*c^6*d^2+234*b^8*c^8*d^6+28*a^2*b^2*c^6*d^3+16*a^2*b^2*c^6*d^4+3*a^2*b^7*c^5*d+546*a^2*b^10*c^8*d^7+312*a^2*b^10*c^8*d^8+32*a^3*b^6*c^12*d^9+624*a^3*b^14*c^14*d^13+7*a^4*b^9*c^5*d^2+4*a^4*b^9*c^5*d^3+8*a^5*b^13*c^11*d^8+12*a^6*b*c^6+36*a^6*b^4*c^12*d^6+4*a^6*b^9*c^8*d^3+702*a^6*b^12*c^14*d^10+78*a^6*b^17*c^10*d^7+28*a^7*b^9*c^5*d^7+546*a^7*b^17*c^7*d^11+3*a^8*c*d^3+332*a^8*b*c^5*d^10+28*a^8*b^3*c^6*d+16*a^8*b^3*c^6*d^2+312*a^8*b^3*c^9*d^4+180*a^8*b^9*c^7*d^14+9*a^8*b^11*c^11*d^5+839*a^8*b^11*c^11*d^8+a^8*b^16*c^7*d^2+16*a^9*b^7*c^5*d^4+32*a^9*b^7*c^12*d^7+324*a^9*b^9*c^6*d^4+312*a^9*b^15*c^7*d^8+7*a^9*b^16*c^4*d^6+24*a^9*b^17*c^8*d^8+7*a^10*b^2*c*d^4+4*a^10*b^2*c*d^5+83*a^10*b^8*c^4*d^9+78*a^10*b^10*c^8*d^3+8*a^11*b^6*c^7*d^10+4*a^11*b^14*c^4*d^3+81*a^11*b^16*c^5*d^3+36*a^12*b^5*c^12*d^4+4*a^12*b^10*c^8*d+28*a^13*b^10*c^5*d^5+332*a^14*b^2*c^5*d^8+9*a^14*b^4*c^7*d^7+312*a^14*b^4*c^9*d^2+a^14*b^9*c^3*d^4+16*a^15*b^8*c^5*d^2+7*a^15*b^9*d^8+324*a^15*b^10*c^6*d^2+83*a^16*b*d^11+78*a^16*b^3*c^4*d^5+4*a^17*b^7*d^5+81*a^17*b^9*c*d^5", domain, LEX);

        lMultivariatePolynomialZp
                la = asLongPolyZp(a),
                lb = asLongPolyZp(b);

        lMultivariatePolynomialZp lgcd = ZippelGCD(la, lb);
        assertTrue(dividesQ(la, lgcd));
        assertTrue(dividesQ(lb, lgcd));


        rnd.setSeed(701);
        lMultivariatePolynomialZp intrp = interpolateGCD(la, lb, lgcd, rnd);
        if (intrp != null)
            assertEquals(lgcd.monic(), intrp.monic());
    }

    @Test
    public void testSparseInterpolation7() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("7*b*c^4*d^6+9*a^2*b*c^8*d^7*e^7+7*a^2*b^6*c^8*d^3*e+8*a^3*c^2*d^6*e^4+7*a^3*b^5*c^7*d^6*e^4+a^4*b^3*c^5*d^8*e^5+25732656*a^5*c^4*d^2*e^3+9*a^5*b^2*c^6*d^5*e^4+25732652*a^6*b^3*c*d*e+25732656*a^7*b^3*c^8*d+a^7*b^3*c^8*d^2*e"),
                b = parse("25732655*a^9*b^8*c^12*d^18*e^13+9*a^11*b^16*c^19*d^11*e^13+4*a^16*b^20*c^17*d^3*e^4+2*a^20*b^10*d^3*e^13+4*a^20*b^11*c^13*d^17*e^9"),
                gcd = parse("3*a^2*b^17*c^14*d^6*e^14+4*a^3*b^14*c^15*d^10*e^8+25732658*a^5*b^17*c^10*d^9*e^12+8*a^6*b^10*c^4*d^3*e^10+25732659*a^6*b^10*c^7*d^5*e^15+a^7*b^2*c^3*d+6*a^9*b^9*c^10*d^6*e^5+3*a^11*b^15*c^7*d^17*e^15+25732652*a^13*b^3*c^5*d^13*e^11+2*a^13*b^12*d^2*e^16+9*a^15*b^2*c^2*d^5*e^4+2*a^15*b^2*c^14*d^14*e^14+a^15*b^13*c^8*e^12+a^16*b*c^10*d^13*e^10");

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        RandomGenerator rnd = getRandom();

        IntegersModulo domain = new IntegersModulo(806213L);
        lMultivariatePolynomialZp
                la = asLongPolyZp(a.setDomain(domain)),
                lb = asLongPolyZp(b.setDomain(domain));

        lMultivariatePolynomialZp skeleton = ZippelGCD(la, lb);

        IntegersModulo domain1 = new IntegersModulo(755899L);
        lMultivariatePolynomialZp
                la1 = asLongPolyZp(a.setDomain(domain1)),
                lb1 = asLongPolyZp(b.setDomain(domain1));

        lMultivariatePolynomialZp gcd0 = ZippelGCD(la1, lb1);
        System.out.println(gcd0.sameSkeleton(skeleton));

        rnd.setSeed(-7756222446675659124L);
        lMultivariatePolynomialZp actual = interpolateGCD(la1, lb1, skeleton.setDomain(la1.domain), rnd);
        if (actual == null) {
            System.out.println("bad evaluation");
            // bad evaluation point => try over
            actual = interpolateGCD(la1, lb1, skeleton.setDomain(la1.domain), rnd);
        }
        assertEquals(gcd0.monic(), actual.monic());
    }

    @Test
    public void testModularGCD1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("a + 17*b + 2*c"),
                b = parse("3*a + b - c"),
                gcd = parse("1273465812736485821734523745*a*b - 21475715234*b - c");
        assertEquals(gcd, ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test
    public void testModularGCD2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("1234324234*a + 12317*b + 2*c"),
                b = parse("3*a + 143423423423423412314*b - c"),
                gcd = parse("1273465812736485821734523745*a*b - 21475715234*b - 143423423423423412314123123*c");
        assertEquals(gcd, ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test
    public void testModularGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("5*b^6*c^15*d^3+4*b^8*c*d^11+17492152*a^2*b^8*c^15*d^10+8*a^2*b^10*d^11+9*a^3*b^2*c^10*d+5*a^4*b*c^5*d^3+6*a^4*b^13*c^3*d^13+17492156*a^8*b^6*c^12*d^4+5*a^9*b^9*d^11+5*a^10*b^6*c^15*d^10"),
                b = parse("b^8*d^3+a^4*b^2*c^7*d+4*a^5*d^2+4*a^5*b^6*c+17492153*a^7*c^8*d^6"),
                gcd = parse("7*a^2*b^7*c^9*d^6+17492158*a^2*b^8*c*d^9+4*a^2*b^9*c^7*d^3+3*a^3*d+7*a^3*b^2*c^2*d^2+4*a^3*b^2*c^2*d^3+17492156*a^3*b^9*c^9*d^3+a^5*b^6*c^9*d^2+17492153*a^6*b^8*c^3*d^3+8*a^9*b^3*c^6*d^8+9*a^9*b^6*c^4*d^5");
        assertEquals(gcd, ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test
    public void testModularGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("4*a*b^7*c*d^2+a^2*b^6*c^8*d^4+7*a^3*b^5*c^6*d^4+5*a^4*b^3*c*d^7+6*a^5*b^4*c^7+8*a^7*c^8*d+a^8*b^5*c^3*d^2"),
                b = parse("25987600*b^18*c^17*d^14+25987597*a*b^9*c^9*d^2+2*a^2*b^7*c^12*d^7+4*a^4*b^14*c^11*d^2+6*a^6*b^2*d+2*a^6*b^3*c^16*d^14+5*a^9*b^17*c^16*d^2+a^14*b^18*c^17*d^4"),
                gcd = parse("25987593*a^4*c^4*d^4+9*a^6*c^3*d^10+7*a^6*b^14*c^4*d^7+8*a^7*b^9*c^13*d+7*a^9*b^2*c^13*d^4+2*a^10*b^6*c^9*d^7+2*a^11*b^5*c^7*d^3+2*a^11*b^12*c^13*d^14+7*a^14*b^8*c^14*d^3+6*a^14*b^13*c^4*d^11");
        assertEquals(gcd, ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)));
    }

    @Test
    public void testModularGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(29);
        MultivariatePolynomial<BigInteger>
                a = parse("5*a*b^5*c^10*d^7*e^16+2*a^4*b^3*c^9*d^6*e^8+5*a^4*b^6*c^16*d^11*e^2+a^4*b^13*d^5*e^6+30844060*a^5*b*c^9*d^8*e^12+4*a^8*b*c^17*d^11*e^3+9*a^8*b^13*c^16*d^17*e^11+a^9*b^2*c^2*d^10*e^14+5*a^9*b^6*c^3*d^7*e^4+7*a^9*b^8*c^3*d^16*e^2+9*a^14*b^5*c^2*d^3*e^16"),
                b = parse("7*b^6*c^18*d^5*e+30844053*a^2*b^8*c^10*d^8*e^6+a^3*b^14*c^4*d^11*e^7+a^4*b^10*c*d^15*e^18+3*a^15*b^9*c^3*e^11+5*a^18*b^13*c^16*d^15*e^15"),
                gcd = parse("9*a^3*b^11*c^7*d^4*e^6+30844059*a^5*b^6*c^15*d^8*e^10+8*a^5*b^10*c^15*d^2*e^9+5*a^10*b^11*c^7*d^9*e^16+2*a^13*b^3*c^13*d^6*e^2+30844060*a^14*b^3*c^6*d^3*e^13+30844055*a^14*b^6*c^4*d^13+30844055*a^14*b^17*c^2*d^8*e^13+2*a^17*b^5*c^7*d*e^11");
        assertTrue(dividesQ(ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test
    public void testModularGCD6() throws Exception {
        for (int i = 46; i < 100; i++) {
            PrivateRandom.getRandom().setSeed(46);
            MultivariatePolynomial<BigInteger>
                    a = parse("8*a*c^5*d^10*e^5+31118523*a*b^3*c^5*d^8*e^10+a^2*b^7*c*d*e^12+4*a^2*b^8*c*d^9*e^10+31118524*a^3*b^5*c^14*d^5*e^13+31118529*a^4*b^3*c^12*d^6*e^8+3*a^5*b^4*d^11*e^9+31118526*a^5*b^8*c^6*d^12*e+4*a^7*b^13*c^11*d^3+31118529*a^9*b^12*c^4*d^2*e^11+5*a^11*b^9*c^2*d*e^11+8*a^13*b^13*c^7*d^2*e^8+8*a^14*b^5*c^14*d^6*e^4"),
                    b = parse("31118526*c^3*d^4*e^2+31118530*b^4*c^6*d^5*e^6+5*a*b*c^4*d^4*e^3+31118527*a*b^3*d*e^2+31118525*a^2*b*c^7*d*e^4+5*a^2*b^4*c^8*d^2*e^5+6*a^2*b^6*d^7*e^5+9*a^2*b^7*c^8*d*e^5+4*a^4*b^6*e^7+3*a^5*b^2*c^6*d^4*e^3+31118529*a^7*b*c^2*d^5*e^8+8*a^7*b^3*c^3*d^4*e^5+7*a^8*b*c^2*d^5*e^8+6*a^8*b^3*c^3*d^5*e^3"),
                    gcd = parse("2*c^3*d*e^5+31118524*b^6*c^2*d^3*e^4+31118528*a^2*b^3*c^2*d^3+7*a^3*b*c^3*d^2+5*a^3*b^3*c^4*d^5*e^2+31118527*a^4*c^2*d^3+7*a^4*b*c*d*e^4+9*a^4*b*c^6*d^3*e^4+5*a^5*d^2*e^2+4*a^6*b^2*c^4*e+7*a^6*b^3*c^5*d^4*e^3");
            assertTrue(dividesQ(ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
        }
    }

    @Test
    public void testModularGCD7() throws Exception {
        PrivateRandom.getRandom().setSeed(46);
        MultivariatePolynomial<BigInteger>
                a = parse("5*b*d^4+2*a^3*b*c*d^4+a^5*b^2*c^2*d^6+6*a^5*b^5*c^3*d^6+4*a^6*b*c^2+8*a^6*b^5*c^5*d^5"),
                b = parse("8*a*b*c^3*d^6+4*a*b^4*c*d+4*a*b^5*c*d^3+3*a^3*b^4*c^2"),
                gcd = parse("5*a^7*b^2*c^13*d^4");
        assertTrue(dividesQ(ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test
    public void testModularGCD8() throws Exception {
        PrivateRandom.getRandom().setSeed(48);
        MultivariatePolynomial<BigInteger>
                a = parse("8*a*b^19*c^11+8*a^3*b^4*c^9+7*a^3*b^10*c^12+3*a^5*b^14*c^21+7*a^9*b^21*c+8*a^10*b^8*c^5+a^14*b^21*c^12+15477328*a^21*b^20*c^8"),
                b = parse("15477335*b^8*c^4+7*b^9*c^8+15477332*a^3*b^13*c^4+15477335*a^9*b^13*c^6+15477328*a^12*c^9"),
                gcd = parse("15477332*a^10*b^13*c^5+7*a^14*b^5*c^3+6*a^19*b^12*c^5+2*a^19*b^12*c^13+15477329*a^20*b*c^19+15477332*a^20*b^8*c^12+7*a^21*b^8*c^2");
        assertTrue(dividesQ(ModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd)), gcd));
    }

    @Test
    public void testModularGCD_random1() throws Exception {
        int nIterations = its(1000, 3000);
        RandomGenerator rnd = getRandom();
        GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> sampleData =
                new GCDSampleDataGeneric<>(Integers, 3, 5, 5, 15, 5, 15, rnd);

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd", MultivariateGCD::ModularGCD));
    }

    @Test
    public void testModularGCD_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>>
                sampleData =
                boundCoefficients(
                        fixVariables(new GCDSampleDataGeneric<>(Integers, nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars),
                        BigInteger.valueOf(100));

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("ModularGCD", MultivariateGCD::ModularGCD));
    }

    @Ignore
    @Test
    public void testModularGCD_performance() throws Exception {
        PrivateRandom.getRandom().setSeed(1232);
        String[] vars = {"a", "b", "c", "d", "e"};
        Domain<BigInteger> domain = Integers;
        MultivariatePolynomial<BigInteger>
                a = parse("2147483167*a^4*b^60*c^57*d^26*e+44*a^8*b^39*c^67*d^22*e^17+38*a^32*b^6*c^13*d^10*e^3+357*a^36*b^34*c^60*d^2*e^59+563*a^42*b^41*c^45*d^52*e^14+257*a^44*b^68*c^43*d^2*e^73+613*a^48*b^55*c^22*d^32*e^19+2147483093*a^52*b^26*c^4*d^72*e^32+19*a^52*b^40*c^26*d^45*e^55+639*a^55*b^72*c^55*d^65", domain, LEX, vars),
                b = parse("2147483150*b^25*c^18*d^62*e^59+2147482723*a^4*b^5*c^65*d^26*e^7+261*a^15*b^60*c^59*d^63*e^53+394*a^27*b^22*c^34*d^54*e^13+952*a^39*b^48*c^17*d^54*e^16+243*a^60*b^15*c^3*d^51*e^46+40*a^61*b^56*c^39*d^40*e^21+555*a^62*b^20*c^20*d^60*e^47+627*a^67*b^8*c^22*d^67*e^61+447*a^70*b^59*c^71*d^24*e^5", domain, LEX, vars),
                gcd = parse("35*a*b^36*c^74*d^62*e^51+376*a^2*b^28*c^64*e^53+893*a^6*b^13*c^60*d^44*e^42+23*a^8*b^71*c^40*d^36*e^11+783*a^20*b^28*c^12*d^31*e^68+2147482938*a^31*b^30*c^40*d^65*e^72+2147482960*a^31*b^49*c^38*d^71*e^55+737*a^47*b^15*c^71*d^13*e^72+868*a^53*b^30*c^40*d^29*e^46+898*a^61*b^71*c^13*d^50*e^66", domain, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        System.out.println(a);
        System.out.println(b);

        for (int i = 0; i < 1000; i++) {
            System.out.println();
            long start = System.nanoTime();
            System.out.println(ModularGCD(a.clone().increment(), b));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            assertTrue(dividesQ(ModularGCD(a, b), gcd));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
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

    @Test
    public void testRationals1() throws Exception {
        MultivariatePolynomial<Rational>
                a = parse("-2/3*a*b*c - 7/6*a^3*c^4 + 2/3*b^3", Rationals),
                b = parse("2/3*a^2*b*c + 1/6*a^3*c^4 + 2/13*c^3", Rationals),
                gcd = parse("12/3*a^2*b*c^2 - 11/6*a^3*b*c^4 + 2/11*c", Rationals);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertTrue(dividesQ(BrownGCD(a, b), gcd));
        assertTrue(dividesQ(ZippelGCD(a, b), gcd));
    }

    @Ignore//too long for rationals
    @Test
    public void testRationals2() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomial<Rational>
                a = parse("(-2/9)*b*c*d^5*e^5-1/10*b*c^5*d^3*e", Rationals, LEX, vars),
                b = parse("b*c^2*e^3-b*c^3*d*e^3+11/7*b^3*d^2*e-2/11*a^3*b^2*c*d^2+4/11*a^3*b^3*e", Rationals, LEX, vars),
                gcd = parse("11/10*b^4*c^4*e+3/2*a*b^3*c*d^4*e^5-1/3*a^4*b^4*c^2*d^6*e^5+14/5*a^6*b^5*c^3*d^6*e+13*a^6*b^6*c^3*e", Rationals, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertTrue(dividesQ(PolynomialGCD(a, b), gcd));
    }

    @Test
    public void testFiniteField1() throws Exception {
        for (int i = 0; i < 10; i++) {
            FiniteField<lUnivariatePolynomialZp> field = FiniteField.GF17p5;
            MultivariatePolynomial<lUnivariatePolynomialZp>
                    a = MultivariatePolynomial.zero(3, field, LEX)
                    .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(1, 2, 3, 4, 5).modulus(17)), 1, 1, 3))
                    .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 1, 3, 2, 13).modulus(17)), 3, 2, 1))
                    .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 11, 13, 12, 13).modulus(17)), 0, 2, 1)),
                    b = MultivariatePolynomial.zero(3, field, LEX)
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(1, 1, 3, 4, 5).modulus(17)), 1, 1, 13))
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 1, 1, 2, 13).modulus(17)), 2, 2, 1))
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 11, 113, 112, 13).modulus(17)), 10, 2, 1)),
                    gcd = MultivariatePolynomial.one(3, field, LEX)
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(1, 1, 3, 4, 5, 12).modulus(17)), 11, 1, 13))
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(11, 2, 1, 1, 2, 13).modulus(17)), 21, 2, 1))
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 111, 113, 112, 13, 12).modulus(17)), 10, 12, 1))
                            .add(MonomialTerm.create(field.valueOf(lUnivariatePolynomialZ.create(2, 111, 113, 112, 13, 12).modulus(17)), 0, 0, 1));

            a = a.clone().add(b).multiply(gcd);
            b = b.clone().subtract(gcd).multiply(gcd);
            long start = System.nanoTime();
            assertTrue(dividesQ(PolynomialGCD(a, b), gcd));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
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
//        return UnivariatePolynomial.create(new MultivariatePolynomials<>(poly), data);
//    }
//
//    static <E> MultivariatePolynomial<E> fromMUnivariate(UnivariatePolynomial<MultivariatePolynomial<E>> poly, int variable) {
//        MultivariatePolynomial<E> zero = poly.domain.getZero();
//        for (int i = 0; i <= poly.degree(); i++) {
//            if (poly.isZeroAt(i))
//                continue;
//
//            MultivariatePolynomial<E> term = poly.get(i);
//            zero.add(
//                    term.multiply(new MonomialTerm<>(zero.nVariables, variable, i, zero.domain.getOne())));
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
//        IntegersModulo domain = new IntegersModulo(3);
//        a = a.setDomain(domain);
//        b = b.setDomain(domain);
//
//
////        System.out.println(MultivariateGCD.PolynomialGCD(a, b));
//        int variable = 0;
////        System.out.println(a);
////        System.out.println(asMUnivariate(a, variable));
////        System.out.println(b);
////        System.out.println(asMUnivariate(b, variable));
//
//        UnivariatePolynomial<MultivariatePolynomial<BigInteger>> result = UnivariateGCD.SubresultantEuclid(asMUnivariate(a, variable), asMUnivariate(b, variable)).gcd();
//        System.out.println(result);
//        System.out.println(fromMUnivariate(result, variable));
//
//    }

    @Test
    public void testSmallDomain1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("a + 2*b + c"),
                b = parse("a*b + 17*a*b^13 + a*c^2"),
                gcd = parse("1 + a^2*b^31*c + c^3*a^3*b - 2*a^5*b^2 - b*c^2");

        for (long modulus : new long[]{2, 3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43}) {
            IntegersModulo domain = new IntegersModulo(modulus);
            MultivariatePolynomial<BigInteger>
                    a1 = a.clone().setDomain(domain),
                    b1 = b.clone().setDomain(domain),
                    gcd1 = gcd.clone().setDomain(domain);
            a1 = a1.multiply(gcd1);
            b1 = b1.multiply(gcd1);

            assertEquals(gcd1.monic(), ModularGCDFiniteField(a1, b1).monic());
        }
    }

    @Test
    public void testSmallDomain_random1() throws Exception {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        int nIterations = its(200, 1000);
        RandomGenerator rnd = getRandom();

        lGCDSampleDataZp sampleData =
                new lGCDSampleDataZp(3, 5, 5, 15, 5, 15, rnd);
        sampleData.minModulusBits = 2;
        sampleData.maxModulusBits = 5;
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::ModularGCDFiniteField),
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::ModularGCDFiniteField));
    }

    @Test
    public void testSmallDomain_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(100, 1000);
        lGCDSampleDataZp source = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        source.minModulusBits = 2;
        source.maxModulusBits = 5;

        GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp>
                sampleData = filterZeros(fixVariables(source, nVars));

        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::ModularGCDFiniteField),
                GCDAlgorithm.named("Modular gcd (small cardinality)", MultivariateGCD::ModularGCDFiniteField));
    }

    @Test
    public void testMultipleGCD1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp
                gcd = lMultivariatePolynomialZp.parse("c*a + b + a + c^15*a^3 + b*c*a^5 + d^2*c*a", domain, LEX),
                arr[] = {
                        lMultivariatePolynomialZp.parse("c*b*a^2 + b^2 + c + b*a^15 + d", domain, LEX).multiply(gcd),
                        lMultivariatePolynomialZp.parse("a^12 + 2*b^12 + 2*c + c*a^5 + d*a", domain, LEX).multiply(gcd),
                        lMultivariatePolynomialZp.parse("a^2 + 2*b^12 + 2*c + c*a^5 + d*a", domain, LEX).multiply(gcd),
                        lMultivariatePolynomialZp.parse("a^12 - 2*b^2 + 2*c^3 + c*a^5 + d*a", domain, LEX).multiply(gcd),
                        lMultivariatePolynomialZp.parse("b^12 - 2*b^2 + 2*c^3 + c*a^5 + d*a", domain, LEX).multiply(gcd)
                };

        lMultivariatePolynomialZp aGcd = MultivariateGCD.PolynomialGCD(arr);
        assertEquals(gcd, aGcd);
    }

    @Test
    public void testCommonZeroes1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1 + c*b*a^2+b^2 + c + a^5", domain, LEX),
                b = lMultivariatePolynomialZp.parse("a^2+2*b^2 + 2*c + a^5", domain, LEX);
        ZeroVariables pZeros = commonPossibleZeroes(a, b, a.nVariables);
        assertTrue(pZeros.pZeros.size() > 0);
        for (BitSet z : pZeros.pZeros) {
            assertFalse(setZeroes(a, z).isZero());
            assertFalse(setZeroes(b, z).isZero());
        }
    }

    @Test
    public void testCommonZeroes2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1 + c*b*a^2+b^2 + c + a^5*d + e*a + b*g - f", domain, LEX),
                b = lMultivariatePolynomialZp.parse("a^2+2*b^2 + 2*c + a^5 + a*b*c*d*e*f + g", domain, LEX);
        lMultivariatePolynomialZp ac = a.clone();
        ZeroVariables pZeros = commonPossibleZeroes(a, b, a.nVariables);
        assertTrue(pZeros.pZeros.size() > 0);
        for (BitSet z : pZeros.pZeros) {
            assertFalse(setZeroes(a, z).isZero());
            assertFalse(setZeroes(b, z).isZero());
        }
    }

    private static lMultivariatePolynomialZp setZeroes(lMultivariatePolynomialZp poly, BitSet zeroes) {
        TIntArrayList vars = new TIntArrayList();
        for (int i = 0; i < zeroes.size(); i++)
            if (zeroes.get(i))
                vars.add(i);

        return poly.evaluateAtZero(vars.toArray());
    }

    @Test
    public void testEZEvaluations1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("4*b^5*d^3*e^6+2*a^4*c^2*d^3*e^2+9939223*a^4*b^3*c^3*e+7*a^4*b^3*c^5*d^5*e+3*a^4*b^5*c^5*d^5*e^4+2*a^5*c^3*d^6*e^3+9939225*a^5*c^5*d*e^3+3*a^5*b*c^6*d^4*e^3+7*a^5*b^2*c^3*d^4*e^4+9939223*a^5*b^6*c^4*d", domain),
                b = lMultivariatePolynomialZp.parse("4*b^6*c^4+6*a*b*c^2*e+9939226*a*b^5*c*d^3*e^3+8*a^2*b^4*c^5*d^5*e^3+3*a^3*c^6*d^4*e^2+9939223*a^3*b^2*c*d^2*e^2+2*a^4*b^3*c^6*d^2*e+2*a^5*b^3*d^6*e^2+4*a^5*b^5*c^5*d*e^4+9*a^6*b*c^5*d^2+9939221*a^6*b^2*c^4*d^3*e^5", domain);
        EZGCDEvaluations evals = new EZGCDEvaluations(a, b, a.nVariables, getRandom());
        for (int j = 0; j < 10; j++) {
            System.out.println(j);
            evals.nextEvaluation();
            assertEquals(a, evals.reconstruct(evals.aReduced));
            assertEquals(b, evals.reconstruct(evals.bReduced));
        }
    }

    @Test
    public void testEZEvaluations2() throws Exception {
        RandomGenerator rnd = getRandom();
        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 4, 2, 3, 5, 7, rnd);
        rnd.setSeed(42);
        GCDSample<lMonomialTerm, lMultivariatePolynomialZp> sample = sampleData.nextSample(false, true);
        for (lMultivariatePolynomialZp[] pp : new lMultivariatePolynomialZp[][]{{sample.a, sample.b}, {sample.aCoFactor, sample.bCoFactor}}) {
            lMultivariatePolynomialZp
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

    @Test
    public void testEZEvaluationsRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        int nIterations = its(10, 100);
        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(3, 4, 2, 3, 5, 7, rnd);
        for (int i = 0; i < nIterations; i++) {
            System.out.println(i);
            GCDSample<lMonomialTerm, lMultivariatePolynomialZp> sample = sampleData.nextSample(false, true);
            for (lMultivariatePolynomialZp[] pp : new lMultivariatePolynomialZp[][]{{sample.a, sample.b}, {sample.aCoFactor, sample.bCoFactor}}) {
                lMultivariatePolynomialZp
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

    @Test
    public void testEZGCD1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^2 + b^2 + c + b*a^15 + d", domain, LEX),
                b = lMultivariatePolynomialZp.parse("a^12 + 2*b^12 + 2*c + c*a^5 + d*a", domain, LEX),
                gcd = lMultivariatePolynomialZp.parse("c*a + b + a + c^15*a^3 + b*c*a^5 + d^2*c*a", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);

        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp

                u = lMultivariatePolynomialZp.parse("c*b*a + b^2 + c + b*a^2 + 1", domain, LEX),
                v = lMultivariatePolynomialZp.parse("2 + a^2 + 2*b^2 + 2*c + c*a^2 + a", domain, LEX),
                a = u.clone().square().multiply(u).multiply(v),
                b = v.clone().square().multiply(u),
                gcd = lMultivariatePolynomialZp.parse("c*a + b + a + c*a^3 + b*c*a^2 + c*a", domain, LEX)
                        .multiply(u).multiply(v).multiply(v);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        lIntegersModulo domain = new lIntegersModulo(15627449);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("3*b^2*c*d+5*a*d^2+15627444*a*b^2*c+15627440*a^2*b^2*c^2", domain, LEX),
                b = lMultivariatePolynomialZp.parse("4*b*c^2*d^2+8*b^2*c*d^2+15627440*b^3*d^2+3*a*c^3+a^2*b", domain, LEX),
                gcd = lMultivariatePolynomialZp.parse("15627440*b^2*c^3+7*b^3*c+6*a*b^2*c^3*d^2+7*a^2*b*d+5*a^3*c^2*d^3+6*a^3*b^2*c^2*d", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        lIntegersModulo domain = new lIntegersModulo(15627449);
        String[] vars = {"b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("13312440+25*d", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("902776+10710357*c+6547542*c^2+4925527*b+2965659*b*c+20*b*c^2+2903103*b^2+40*b^2*c+15627404*b^3", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("102871+8266210*d+5121205*d^2+16248*d^3+1722152*c+2574791*c*d+10788581*c*d^2+15247596*c*d^3+8472569*c^2+898837*c^2*d+14099452*c^2*d^2+c^2*d^3", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(0);
        lIntegersModulo domain = new lIntegersModulo(24254707);
        String[] vars = {"b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("24254706*b+b^2*c+7*a*c+5*a*b^2+2*a^2*b*c+24254705*a^2*b*c^2+a^2*b^2*c^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("4*b*c+6*b^2*c+4*a^2+7*a^2*b*c^2+3*a^2*b^2*c", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("9*c+24254705*a*c^2+6*a*b^2+4*a^2+3*a^2*c+24254698*a^2*b^2*c", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD6() throws Exception {
        PrivateRandom.getRandom().setSeed(0);
        lIntegersModulo domain = new lIntegersModulo(24254707);
        String[] vars = {"b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("11106134 + 20915017*a + 948048*a^2 + 18101438*b + 8523620*a*b + 19589342*a^2*b + 13684993*b^2 + 8219998*a*b^2 + 24254698*a^2*b^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("12760358 + 5425698*a + 5129306*a^2 + 14380683*b + 16257092*a*b + 24254680*a^2*b + 4515223*b^2 + 24254644*a*b^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("9740181 + 21604578*a + 9691573*a^2 + 11040951*b + 17951441*a*b + a^2*b", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD7() throws Exception {
        PrivateRandom.getRandom().setSeed(6);
        lIntegersModulo domain = new lIntegersModulo(28565767);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6*c^2*d+3*a*d^2+9*a*d^3+28565764*a^2*b*d^2+2*a^2*b^2*c*d^2+9*a^3*c*d^2+3*a^3*b^2*c*d^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("c^2*d^2+28565766*a*b^2*d^2+2*a*b^2*c*d+28565763*a^2*c^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("6*c*d+7*c*d^2+7*a*b^2+8*a^2*b*c^2", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD8() throws Exception {
        PrivateRandom.getRandom().setSeed(6);
        lIntegersModulo domain = new lIntegersModulo(28565767);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1564663 + 63*d", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("15383307 + 22520298*b + 12522045*b^2 + 9677552*c + 7*c^2 + 3221049*d + 5785123*b*d + 28565760*b^2*d", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("21095373 + c", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD9() throws Exception {
        PrivateRandom.getRandom().setSeed(27);
        lIntegersModulo domain = new lIntegersModulo(8678027);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6*c^2*d+2*b^2*c^2*d^2+a*c^2+a*b^2*c^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("3*b^3*c*d^2+2*a*b^2*d^2+8678026*a^2*c^3*d^3+8678018*a^2*b+5*a^3*b^2*c^3*d", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("8678020*c^2+3*b*c+6*b^2*c^2*d+2*a^2*b^3*c+8678020*a^3*c*d^3", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD10() throws Exception {
        PrivateRandom.getRandom().setSeed(4);
        lIntegersModulo domain = new lIntegersModulo(28996511);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6*b^2*c^2+b^3*c+28996505*a*b*c+28996507*a*b^2+4*a*b^2*c^3", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("6*b^2*c+12*a*b^2*c+3*a^2*b^2*c^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("9+28996505*c^2+7*a*b*c+8*a*b*c^2+9*a^2*b^2*c^2", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD11() throws Exception {
        PrivateRandom.getRandom().setSeed(11);
        lIntegersModulo domain = new lIntegersModulo(22687397);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("7*d+7*b^2*c^2+3*a*c^2+2*a^2*d+6*a^2*d^2+9*a^2*c^2*d", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("2*a*b*c^2+22687392*a*b^2*c^2*d^2+9*a*b^2*c^3+22687395*a*b^3*c*d+3*a^2*c*d^3+6*a^3+8*a^3*b^2*c*d^3", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("22687391*d^2+22687393*b*c^2*d^2+5*a+5*a*b+8*a*b^2*c^2*d^2+4*a^2*b^2*c^2", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD12() throws Exception {
        PrivateRandom.getRandom().setSeed(11);
        lIntegersModulo domain = new lIntegersModulo(22687397);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("11244809 + 22687361*a^2 + 30*b + 11244809*c + 30*b*c", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("21851038 + 15893672*a + 14760820*a^2 + 12564491*a^3 + 11694181*a^4 + 16683489*a^5 + 11244809*b + 9237198*a*b + 4795625*a^2*b + 6123527*a^3*b + 5432952*a^4*b + 22659623*a^5*b + 15*b^2 + 14306815*a*b^2 + 15896207*a^2*b^2 + 9051*a^3*b^2 + 22645277*a^4*b^2 + 11745*a*b^3 + 35100*a^2*b^3 + 2299511*c + 15893672*a*c + 1652686*a^2*c + 6762669*a^3*c + 19488092*a^4*c + 3757801*a^5*c + 11247299*b*c + 9237198*a*b*c + 12416836*a^2*b*c + 19792351*a^3*b*c + 9233182*a^4*b*c + 22686257*a^5*b*c + 15*b^2*c + 14306815*a*b^2*c + 15890231*a^2*b^2*c + 22668450*a^3*b^2*c + 22687373*a^4*b^2*c + 11745*a*b^3*c + 35100*a^2*b^3*c + 9448492*c^2 + 11808378*a*c^2 + 20057596*a^2*c^2 + 20794627*a^3*c^2 + 18785592*a^4*c^2 + 4960279*a^5*c^2 + 21841949*b*c^2 + 2669531*a*b*c^2 + 13194739*a^2*b*c^2 + 10436283*a^3*b*c^2 + 19195830*a^4*b*c^2 + 5885213*a^5*b*c^2 + 8273012*b^2*c^2 + 20876792*a*b^2*c^2 + 19005620*a^2*b^2*c^2 + 6957603*a^3*b^2*c^2 + 17675248*a^4*b^2*c^2 + 22615613*a^5*b^2*c^2 + 7430983*b^3*c^2 + 18240653*a*b^3*c^2 + 4232985*a^2*b^3*c^2 + 147044*a^3*b^3*c^2 + 48*a^4*b^3*c^2 + 5976*b^4*c^2 + 42092*a*b^4*c^2 + 24*a^2*b^4*c^2 + 9149236*c^3 + 9268605*a*c^3 + 4080916*a^2*c^3 + 6163313*a^3*c^3 + 2920*a^4*c^3 + 3127846*a^5*c^3 + 11966525*b*c^3 + 22671967*a*b*c^3 + 22683508*a^2*b*c^3 + 21437536*a^3*b*c^3 + 4916674*a^5*b*c^3 + 1168*b^2*c^3 + 2213725*a^3*b^2*c^3 + 22684357*a^5*b^2*c^3 + 61720*a^3*b^3*c^3 + 12936437*c^4 + 9268605*a*c^4 + 15448872*a^2*c^4 + 21785946*a^3*c^4 + 4125558*a^4*c^4 + 19072832*a^5*c^4 + 19011898*b*c^4 + 22671967*a*b*c^4 + 4208806*a^2*b*c^4 + 15458267*a^3*b*c^4 + 22681557*a^4*b*c^4 + 3920321*a^5*b*c^4 + 41624*b^2*c^4 + 22684477*a^2*b^2*c^4 + 3920321*a^3*b^2*c^4 + 7937375*a^5*b^2*c^4 + 12854049*a^3*b^3*c^4 + 6080*a^5*b^3*c^4 + 3040*a^3*b^4*c^4 + 12125660*a*c^5 + 12572867*a^3*c^5 + 12572867*a*b*c^5 + 22685877*a^3*b*c^5 + 22686637*a*b^2*c^5", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("1", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD13() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(13666309);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("4*b^2*c^3+3*a*b^3+2*a*b^3*c+4*a*b^3*c^3+2*a^3*b^2*c", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("5*b*c+5*a*b+13666307*a*b*c+7*a^2*c+6*a^2*c^2+13666307*a^2*b", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("4*b^2+15*a*b^2*c+4*a^2+13666308*a^2*c", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD14() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(13666309);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("3338057+430735*b+13248829*b^2+11374034*b^3+3423812*a+698808*a*b+7995810*a*b^2+60*a*b^3+8933188*a^2+13666305*a^2*b", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("7990553+3359122*b+846494*b^2+131346*a+12831229*a*b+6789484*a*b^2+12272859*a^2+7995810*a^2*b+90*a^2*b^2+13666303*a^3", domain, LEX, vars);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }


    @Test
    public void testEZGCD15() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(22687397);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("30*b+30*b*c+22687361*a^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("2490*b*c+6800*b*c^2+4310*b*c^3+15*b^2+15*b^2*c+1168*b^2*c^3+41624*b^2*c^4+14028*b^3*c^2+5976*b^4*c^2+40950*a*b+40950*a*b*c+22671967*a*b*c^3+22671967*a*b*c^4+28*a*b^2*c^2+22686637*a*b^2*c^5+11745*a*b^3+11745*a*b^3*c+42092*a*b^4*c^2+22684409*a^2*c+22681057*a^2*c^2+22645773*a^2*c^3+22687379*a^2*b+22673369*a^2*b*c+22681172*a^2*b*c^2+22683508*a^2*b*c^3+83248*a^2*b*c^4+22681421*a^2*b^2*c+28056*a^2*b^2*c^2+22684477*a^2*b^2*c^4+35100*a^2*b^3+35100*a^2*b^3*c+11952*a^2*b^3*c^2+24*a^2*b^4*c^2+22638257*a^3+22687369*a^3*c+18516*a^3*c^3+760*a^3*c^4+56*a^3*b*c^2+22685877*a^3*b*c^5+9051*a^3*b^2+22668450*a^3*b^2*c+147044*a^3*b^3*c^2+61720*a^3*b^3*c^3+3040*a^3*b^4*c^4+7470*a^4*c^2+2920*a^4*c^3+22681557*a^4*b*c^4+22645277*a^4*b^2+22687373*a^4*b^2*c+48*a^4*b^3*c^2+22659623*a^5*b+22686257*a^5*b*c+22615613*a^5*b^2*c^2+22684357*a^5*b^2*c^3+6080*a^5*b^3*c^4", domain, LEX, vars);
//                gcd = lMultivariatePolynomialZp.parse("4*b^2+15*a*b^2*c+4*a^2+13666308*a^2*c", domain, LEX, vars);
//        a = a.multiply(gcd);
//        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD16() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("8+3*a+6*a^2*b*c^2+6*a^2*b^2*c", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("6*b^2+8*b^2*c+9607986*a*b+a*b*c^2+a^2*c+2*a^2*b+7*a^2*b*c^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("6*b*c+6*a*b*c+9607979*a*b^2*c^2+3*a*b^3*c^3+9607982*a^3*b^3", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEZGCD17() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6998457+9068733*c+9042619*c^2+1283895*c^3+3738482*c^4+4888116*c^5+5574926*b+6161435*b*c+3428490*b*c^2+1423636*b*c^3+7718978*b*c^4+9607957*b*c^5+5803370*b^2+3797801*b^2*c+1899022*b^2*c^2+1286548*b^2*c^3+6635895*b^3+1871925*b^3*c+2687295*b^3*c^2+9406093*b^3*c^3+3101467*b^4+8800611*b^4*c+5291877*b^4*c^2+18*b^4*c^3", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("4572338+8988826*c+5700575*c^2+6517458*c^3+7602268*b+8674751*b*c+4253585*b*c^2+9607947*b*c^3+1248945*b^2+3165308*b^2*c+9363637*b^3+9338813*b^3*c+8757270*b^4+24*b^4*c", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("1", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD_random1() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 4, minDegree = 2, maxDegree = 3, minSize = 5, maxSize = 7;
        int nIterations = its(100, 1000);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EZ-GCD", MultivariateGCD::EZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test
    public void testEEZGCD_sparse_variables_random() throws Exception {
        RandomGenerator rnd = getRandom();

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp>
                sampleData = fixVariables(
                new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd), nVars);

        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Brown", MultivariateGCD::BrownGCD),
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EZ-GCD", MultivariateGCD::EZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test
    public void testEEZGCD1() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6998457+9068733*c+9042619*c^2+1283895*c^3+3738482*c^4+4888116*c^5+5574926*b+6161435*b*c+3428490*b*c^2+1423636*b*c^3+7718978*b*c^4+9607957*b*c^5+5803370*b^2+3797801*b^2*c+1899022*b^2*c^2+1286548*b^2*c^3+6635895*b^3+1871925*b^3*c+2687295*b^3*c^2+9406093*b^3*c^3+3101467*b^4+8800611*b^4*c+5291877*b^4*c^2+18*b^4*c^3", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("4572338+8988826*c+5700575*c^2+6517458*c^3+7602268*b+8674751*b*c+4253585*b*c^2+9607947*b*c^3+1248945*b^2+3165308*b^2*c+9363637*b^3+9338813*b^3*c+8757270*b^4+24*b^4*c", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("1", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD2() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(16604789);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("6*a*c*d^2+2*a*c^2+2*a*b^2*d+3*a*b^2*c^2*d^2", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("7*b^2*c^3*d+6*a^2*b^2*c^3+16604785*a^2*b^3*d+9*a^3*b*c^2+16604781*a^3*b*c^3", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("4*c^2*d+8*a+4*a*c^3*d^2+16604785*a^2*d^3+2*a^2*b^2*c^3*d", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD3() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(31012727);
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("4*c+3*b*c^2+4*a*c+31012723*a*b+2*a^2*c+7*a^2*b*c", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("7*c^2+2*b*c^2+3*a*b*c^2+31012724*a*b^2+a^2*b^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("31012726*b*c+31012724*b*c^2+3*a*b", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD4() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        lIntegersModulo domain = new lIntegersModulo(31012727);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("31012718 + 3*b", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("12*a + 47454*a^2 + 30989135*b + 12*a*b + 9801*a^2*b + 41292*a*b^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("20675151*a + 31012726*a^2 + b", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD5() throws Exception {
        PrivateRandom.getRandom().setSeed(78);
        lIntegersModulo domain = new lIntegersModulo(26662411);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("26662407*a^3*b^3*c*d^5*e^6+3*a^3*b^6*c^3*e^2+4*a^4*b^2*c^3*e^3+8*a^4*b^5*c^6*d^4*e^4+7*a^6*b*c^3*d", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("8*c+26662402*b*c*d+8*b^2*c*d*e+b^2*c^2+a*b*e^2+a^2*d*e^2+8*a^2*d^2*e^2+5*a^2*b*c*d^2*e^2+15*a^2*b*c^2*d", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("5*b^2*d^2*e^2+9*b^2*c^2*e+26662408*b^2*c^2*e^2+3*a*c^2*e+3*a*b*d^2+26662408*a*b^2*d*e+26662407*a*b^2*d^2*e^2+a^2*d^2*e^2+26662402*a^2*b^2*c*e^2", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        System.out.println(EEZGCD(a, b));
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }


    @Test
    public void testEEZGCD6() throws Exception {
        PrivateRandom.getRandom().setSeed(78);
        lIntegersModulo domain = new lIntegersModulo(26662411);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("15*c^2*d+27*b^2+26662402*b^2*d+26662402*a*c+26662399*a*c^2*d+26662384*a^2*b*d", domain, LEX, vars),
                b = lMultivariatePolynomialZp.parse("26657131*c^7*d^8+46200*b*c^2*d^2+26593651*b*c^3*d^2+20960*b*c^3*d^3+2620*b^2*c^2*d^2+26652907*b^2*c^5*d^7+3168*b^2*c^5*d^8+83160*b^3*d+26634691*b^3*d^2+26538643*b^3*c*d+78984*b^3*c*d^2+26649835*b^3*c*d^3+4716*b^4*d+26660839*b^4*d^2+7640*a*c^2*d^4+3168*a*c^6*d^7+26651827*a*c^7*d^6+4224*a*c^7*d^8+26634691*a*b*c*d+35952*a*b*c^2+41256*a*b*c^2*d+26612875*a*b*c^2*d^2+26631226*a*b*c^3+36672*a*b*c^3*d+55008*a*b*c^3*d^2+26645643*a*b*c^3*d^3+13752*a*b^2*d^3+26657827*a*b^2*d^4+26660839*a*b^2*c*d+4584*a*b^2*c^2+26660324*a*b^2*c^2*d^2+17640*a*b^2*c^2*d^5+26644891*a*b^2*c^5*d^7+35976*a*b^3*d+26621965*a*b^3*c*d+27720*a*b^3*c*d^2+3465*a*b^4*d+13833*a*b^4*d^3+31752*a*b^4*d^4+26651827*a*b^4*d^5+40*a*b^5*c^6*d^6+72*a*b^7*c^4*d^5+26662387*a*b^7*c^4*d^6+26657827*a^2*c*d^3+3465*a^2*c^2*d^2+26656299*a^2*c^2*d^4+5775*a^2*c^3*d^4+46200*a^2*c^4*d^4+26656571*a^2*c^7*d^8+11992*a^2*b*c^2*d^2+26648929*a^2*b*c^3*d^2+9240*a^2*b*c^3*d^3+38200*a^2*b*c^4*d^4+9504*a^2*b*c^5*d^8+26579251*a^2*b^2*d^2+4494*a^2*b^2*d^3+123768*a^2*b^2*c*d^2+26635078*a^2*b^2*c*d^3+26648362*a^2*b^2*c*d^4+1155*a^2*b^2*c^2*d^2+100680*a^2*b^2*c^2*d^3+26639302*a^2*b^2*c^2*d^4+26648299*a^2*b^2*c^2*d^5+114600*a^2*b^2*c^3*d^2+26657695*a^2*b^3*d^2+68760*a^2*b^3*c^2*d^3+26639491*a^2*b^3*c^2*d^4+3300*a^2*b^4*d^4+206280*a^2*b^4*c*d+26593651*a^2*b^4*c*d^2+26662387*a^2*b^5*c^5*d^5+36888*a^2*b^5*c^6*d^4+26662379*a^2*b^5*c^6*d^6+6336*a^2*b^7*c^4*d^5+26658946*a^3*c^2*d^3+1498*a^3*c^2*d^4+4494*a^3*c^3*d^2+26634691*a^3*c^3*d^3+26657791*a^3*c^3*d^4+35952*a^3*c^4*d^2+26625451*a^3*c^4*d^4+26648659*a^3*b*d^4+26639491*a^3*b*c^3*d^3+17325*a^3*b*c^4*d^2+26631851*a^3*b*c^4*d^4+4497*a^3*b^2*c*d^3+26593651*a^3*b^2*c^2*d+35976*a^3*b^2*c^2*d^3+1100*a^3*b^2*c^2*d^5+51975*a^3*b^2*c^3+26621831*a^3*b^2*c^3*d^2+26630659*a^3*b^3*d^5+22470*a^3*b^3*c^2*d^3+159390*a^3*b^4*c*d+26631751*a^3*b^4*c*d^2+2112*a^3*b^5*c^6*d^6+26662339*a^3*b^6*c^4*d^6+1499*a^4*c^3*d^4+11992*a^4*c^4*d^4+26652016*a^4*b*c*d^4+26579251*a^4*b*c^2*d^4+7490*a^4*b*c^4*d^4+26631751*a^4*b^2*c^2*d+26593651*a^4*b^2*c^2*d^4+5775*a^4*b^2*c^3+26644001*a^4*b^2*c^3*d^2+26456131*a^4*b^3*c*d^2+15456*a^4*b^4*c*d+5152*a^5*b^2*c^3*d^2+26570431*a^5*b^3*c*d^2", domain, LEX, vars),
                gcd = lMultivariatePolynomialZp.parse("1", domain, LEX, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testEEZGCD_random2() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 2, maxDegree = 5, minSize = 5, maxSize = 17;
        int nIterations = its(100, 200);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithms(sampleData, nIterations,
                GCDAlgorithm.named("Zippel", MultivariateGCD::ZippelGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test
    public void testEEZGCD_random3() throws Exception {
        RandomGenerator rnd = getRandom();
        int nVarsMin = 3, nVarsMax = 5, minDegree = 2, maxDegree = 5, minSize = 5, maxSize = 17;
        int nIterations = its(100, 200);

        lGCDSampleDataZp sampleData = new lGCDSampleDataZp(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        testGCDAlgorithm(sampleData, nIterations,
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD),
                GCDAlgorithm.named("EEZ-GCD", MultivariateGCD::EEZGCD));
    }

    @Test
    public void testEEZGCD7() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(BigPrimes.nextPrime(1321323));
        lMultivariatePolynomialZp

                u = lMultivariatePolynomialZp.parse("c*b*a + b^2 + c + b*a^2 + 1", domain, LEX),
                v = lMultivariatePolynomialZp.parse("2 + a^2 + 2*b^2 + 2*c + c*a^2 + a", domain, LEX),
                a = u.clone().square().multiply(u).multiply(v),
                b = v.clone().square().multiply(u),
                gcd = lMultivariatePolynomialZp.parse("c*a + b + a + c*a^3 + b*c*a^2 + c*a", domain, LEX)
                        .multiply(u).multiply(v).multiply(v);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        assertEquals(ZippelGCD(a, b).monic(), EEZGCD(a, b).monic());
    }

    @Test
    public void testPolynomialGCD1() throws Exception {
        // nUsedVariables == 1
        lIntegersModulo domain = new lIntegersModulo(19);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp a = lMultivariatePolynomialZp.parse("2 + 3*b + b^2", domain, vars);
        lMultivariatePolynomialZp b = lMultivariatePolynomialZp.parse("2 + b + a^2 + a^2*b", domain, vars);
        for (int i = 0; i < 1000; i++)
            assertTrue(PolynomialGCD(a, b).isConstant());
    }

    @Test
    public void testPolynomialGCD2() throws Exception {
        // nUsedVariables == 1
        lIntegersModulo domain = new lIntegersModulo(19);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^2 + b*a", domain, vars),
                b = lMultivariatePolynomialZp.parse("b + 1", domain, vars),
                gcd = lMultivariatePolynomialZp.parse("2 + b", domain, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        for (int i = 0; i < 1000; i++)
            assertEquals(gcd, PolynomialGCD(a, b));
    }

    @Test(timeout = 1000)
    public void testPolynomialGCD3() throws Exception {
        // nUsedVariables == 1
        lIntegersModulo domain = new lIntegersModulo(2);
        String[] vars = {"a", "b"};
        // very tricky example
        // for both a = 0 and a = 1, the gcd is (1 + b^2), while the true gcd is only 1 + b
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^2*b + a*b + b + 1", domain, vars),
                b = lMultivariatePolynomialZp.parse("b + 1", domain, vars),
                gcd = lMultivariatePolynomialZp.parse("b + 1", domain, vars);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        for (int i = 0; i < 10; i++)
            assertEquals(gcd, PolynomialGCD(a, b));
    }


    /* =============================================== Test data =============================================== */

    /** sample data for test of GCD */
    public static final class GCDSample<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** sample polynomials */
        public final Poly a, b;
        /** their GCD */
        public final Poly gcd;
        /** a/gcd and b/gcd */
        public final Poly aCoFactor, bCoFactor;
        /** either generic Domain or lIntegersModulo */
        public final Object domain;

        public GCDSample(Poly aCoFactor, Poly bCoFactor, Poly gcd) {
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
            this.gcd = gcd;
            this.a = aCoFactor.clone().multiply(gcd);
            this.b = bCoFactor.clone().multiply(gcd);
            checkConsistency(a, b, gcd, aCoFactor, bCoFactor);
            this.domain = (gcd instanceof lMultivariatePolynomialZp)
                    ? ((lMultivariatePolynomialZp) gcd).domain
                    : ((MultivariatePolynomial) gcd).domain;
        }
    }

    /** substitute random values for some variables */
    public static GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>>
    boundCoefficients(final GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> source, BigInteger bound) {
        final IntegersModulo domain = new IntegersModulo(bound);
        return new GCDSampleData<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>>() {
            @Override
            GCDSample<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> nextSample0(boolean primitive, boolean monic) {
                GCDSample<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> sample = source.nextSample(primitive, monic);
                return new GCDSample<>(
                        setRandomCoefficients(sample.aCoFactor, bound),
                        setRandomCoefficients(sample.bCoFactor, bound),
                        setRandomCoefficients(sample.gcd, bound));
            }

            private MultivariatePolynomial<BigInteger> setRandomCoefficients(MultivariatePolynomial<BigInteger> poly, BigInteger bound) {
                RandomGenerator rnd = PrivateRandom.getRandom();
                MultivariatePolynomial<BigInteger> r = poly.createZero();
                for (MonomialTerm<BigInteger> term : poly)
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
    public static <Term extends DegreeVector<Term>,
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
    public static <Term extends DegreeVector<Term>,
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
    public static <Term extends DegreeVector<Term>,
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
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {

        final DescriptiveStatistics
                medFactorsSize = new DescriptiveStatistics(),
                medFactorsDegree = new DescriptiveStatistics(),
                medGCDSize = new DescriptiveStatistics(),
                medGCDDegree = new DescriptiveStatistics(),
                medFactorsUsedVariables = new DescriptiveStatistics(),
                medGCDUsedVariables = new DescriptiveStatistics();

        public boolean primitive = false, monic = false;

        final GCDSample<Term, Poly> nextSample() {
            return nextSample(primitive, monic);
        }

        final GCDSample<Term, Poly> nextSample(boolean primitive, boolean monic) {
            GCDSample<Term, Poly> sample = nextSample0(primitive, monic);
            medFactorsDegree.addValue(sample.a.degreeSum());
            medFactorsDegree.addValue(sample.b.degreeSum());
            medGCDDegree.addValue(sample.gcd.degreeSum());

            medFactorsSize.addValue(sample.a.size());
            medFactorsSize.addValue(sample.b.size());
            medGCDSize.addValue(sample.gcd.size());

            medFactorsUsedVariables.addValue(sample.a.nUsedVariables());
            medFactorsUsedVariables.addValue(sample.b.nUsedVariables());
            medGCDUsedVariables.addValue(sample.gcd.nUsedVariables());
            return sample;
        }

        abstract GCDSample<Term, Poly> nextSample0(boolean primitive, boolean monic);

        final void printSamplesStatistics() {
            System.out.println("Median polys size : " + medFactorsSize.getPercentile(0.5));
            System.out.println("Median gcd   size : " + medGCDSize.getPercentile(0.5));
            System.out.println("Median polys deg  : " + medFactorsDegree.getPercentile(0.5));
            System.out.println("Median gcd   deg  : " + medGCDDegree.getPercentile(0.5));
            System.out.println("Median poly nVars : " + medFactorsUsedVariables.getPercentile(0.5));
            System.out.println("Median gcd  nVars : " + medGCDUsedVariables.getPercentile(0.5));
        }
    }

    /** provides random data for GCD tests */
    public static abstract class AGCDSampleData<
            Term extends DegreeVector<Term>,
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
    public static final class lGCDSampleDataZp extends AGCDSampleData<lMonomialTerm, lMultivariatePolynomialZp> {
        public lGCDSampleDataZp(int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            super(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        }

        public int minModulusBits = 24, maxModulusBits = 32;

        @Override
        public GCDSample<lMonomialTerm, lMultivariatePolynomialZp> nextSample1(boolean primitive, boolean monic) {
            long modulus = AbstractPolynomialTest.getModulusRandom(rndd.nextInt(minModulusBits, maxModulusBits));
            lIntegersModulo domain = new lIntegersModulo(modulus);

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            lMultivariatePolynomialZp[] data = new lMultivariatePolynomialZp[3];
            for (int i = 0; i < 3; i++) {
                lMultivariatePolynomialZp
                        p = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), domain, LEX, rnd);

                // make coefficients small
                p = p.setDomain(66).setDomain(domain);

                if (primitive)
                    p = lMultivariatePolynomialZp.asNormalMultivariate(p.asOverUnivariateEliminate(0).primitivePart(), 0);

                if (monic)
                    // make monic in main variable by adding some monomial
                    p.add(p.createUnivariateMonomial(0, p.degreeSum() + 1));

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
    public static final class GCDSampleDataGeneric<E> extends AGCDSampleData<MonomialTerm<E>, MultivariatePolynomial<E>> {
        final Domain<E> domain;

        public GCDSampleDataGeneric(Domain<E> domain, int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            super(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
            this.domain = domain;
        }

        @Override
        public GCDSample<MonomialTerm<E>, MultivariatePolynomial<E>> nextSample1(boolean primitive, boolean monic) {
            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            @SuppressWarnings("unchecked")
            MultivariatePolynomial<E>[] data = new MultivariatePolynomial[3];
            for (int i = 0; i < 3; i++) {
                MultivariatePolynomial<E>
                        p = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), domain, LEX, rnd);

                if (primitive)
                    p = MultivariatePolynomial.asNormalMultivariate(p.asOverUnivariateEliminate(0).primitivePart(), 0);

                if (monic)
                    // make monic in main variable by adding some monomial
                    p.add(p.createUnivariateMonomial(0, p.degreeSum() + 1));

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
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static void testGCDAlgorithm(GCDSampleData<lMonomialTerm, lMultivariatePolynomialZp> sampleData, int nIterations,
                                        GCDAlgorithm<lMultivariatePolynomialZp> lAlgorithm,
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

            GCDSample<lMonomialTerm, lMultivariatePolynomialZp> sample = sampleData.nextSample();

            MultivariatePolynomial<BigInteger> result = null;
            lMultivariatePolynomialZp lResult = null;
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
        System.out.println("domain     : " + sample.domain);
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
            checkConsistency((lMultivariatePolynomialZp) poly);
    }

    private static <E> void checkConsistency(MultivariatePolynomial<E> poly) {
        Domain<E> domain = poly.domain;
        for (MonomialTerm<E> e : poly.terms) {
            E value = e.coefficient;
            assertFalse(domain.isZero(value));
            assertTrue(value == domain.valueOf(value));
            if (domain instanceof IntegersModulo) {
                assertTrue(domain.signum(value) > 0);
                assertTrue(((BigInteger) value).compareTo(((IntegersModulo) domain).modulus) <= 0);
            }
        }
    }

    private static void checkConsistency(lMultivariatePolynomialZp poly) {
        lIntegersModulo domain = poly.domain;
        for (lMonomialTerm e : poly.terms) {
            long value = e.coefficient;
            assertFalse(value == 0);
            assertTrue(value == domain.modulus(value));
        }
    }
}