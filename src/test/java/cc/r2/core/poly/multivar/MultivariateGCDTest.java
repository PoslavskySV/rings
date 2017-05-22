package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.DegreeVector.LEX;
import static cc.r2.core.poly.multivar.MultivariateGCD.*;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
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
        lMultivariatePolynomial lActualGCD = BrownGCD(asLongPolyZp(a), asLongPolyZp(b));
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
        rnd.setSeed(123);

        int nVarsMin = 3, nVarsMax = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            PrivateRandom.getRandom().setSeed(n);
            MultivariatePolynomial<BigInteger> gcdActual = null;
            GCDTriplet data = sampleData.nextSample(false, false);
            checkConsistency(data.a, data.b, data.gcd, data.aGCD, data.bGCD);
            try {
                PrivateRandom.getRandom().setSeed(n);
                gcdActual = BrownGCD(data.aGCD, data.bGCD);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(gcdActual, data.gcd));

                PrivateRandom.getRandom().setSeed(n);
                assertTrue(dividesQ(BrownGCD(asLongPolyZp(data.aGCD), asLongPolyZp(data.bGCD)), asLongPolyZp(data.gcd)));
            } catch (Throwable err) {
                System.out.println("seed: " + n);
                System.out.println("modulus: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("aGCD: " + data.aGCD);
                System.out.println("bGCD: " + data.bGCD);
                System.out.println("expected gcd: " + data.gcd);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
    }

    @Test
    public void testBrown_random2() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int nVarsMin = 3, nVarsMax = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            PrivateRandom.getRandom().setSeed(n);
            MultivariatePolynomial<BigInteger> gcdActual = null;
            GCDTriplet data = sampleData.nextSample(false, false);

            checkConsistency(data.a, data.b);
            try {
                PrivateRandom.getRandom().setSeed(n);
                gcdActual = BrownGCD(data.a, data.b);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(data.a, gcdActual));
                assertTrue(dividesQ(data.b, gcdActual));

                PrivateRandom.getRandom().setSeed(n);
                lMultivariatePolynomial lGcdActual = BrownGCD(asLongPolyZp(data.a), asLongPolyZp(data.b));
                assertTrue(dividesQ(asLongPolyZp(data.a), lGcdActual));
                assertTrue(dividesQ(asLongPolyZp(data.b), lGcdActual));
            } catch (Throwable err) {
                System.out.println("seed: " + n);
                System.out.println("modulus: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
    }

    @Test
    public void testBrown_random3() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);

        int nVarsMin = 100, nVarsMax = 100, nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            PrivateRandom.getRandom().setSeed(n);
            MultivariatePolynomial<BigInteger> gcdActual = null;
            GCDTriplet data = sampleData.nextSample(false, false);

            MultivariatePolynomial<BigInteger> aGCD = data.a;
            MultivariatePolynomial<BigInteger> bGCD = data.b;
            MultivariatePolynomial<BigInteger> gcd = data.gcd;

            do {
                gcd = gcd.evaluate(rnd.nextInt(gcd.nVariables), data.domain.randomElement(rnd));
            } while (gcd.nUsedVariables() > nVars);

            int[] gcdDegrees = gcd.degrees();
            for (int i = 0; i < gcdDegrees.length; i++) {
                if (gcdDegrees[i] == 0) {
                    if (rnd.nextBoolean())
                        aGCD = aGCD.evaluate(i, data.domain.randomElement(rnd));
                    else
                        bGCD = bGCD.evaluate(i, data.domain.randomElement(rnd));
                }
            }

            aGCD = aGCD.multiply(gcd);
            bGCD = bGCD.multiply(gcd);

            checkConsistency(data.a, data.b, gcd, aGCD, bGCD);
            try {
                PrivateRandom.getRandom().setSeed(n);
                gcdActual = BrownGCD(aGCD, bGCD);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(gcdActual, gcd));

                PrivateRandom.getRandom().setSeed(n);
                lMultivariatePolynomial lGcdActual = BrownGCD(asLongPolyZp(aGCD), asLongPolyZp(bGCD));
                assertTrue(dividesQ(lGcdActual, asLongPolyZp(gcd)));
            } catch (Throwable err) {
                System.out.println("seed: " + n);
                System.out.println("modulus: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("aGCD: " + aGCD);
                System.out.println("bGCD: " + bGCD);
                System.out.println("expected gcd: " + gcd);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
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
        lMultivariatePolynomial lActualGCD = ZippelGCD(asLongPolyZp(a), asLongPolyZp(b));
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

        lMultivariatePolynomial la = asLongPolyZp(a), lb = asLongPolyZp(b),
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
        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n % 100 == 0) System.out.println(n);

            GCDTriplet data = sampleData.nextSample(true, true);

            int variable = data.a.nVariables - 1;
            BigInteger seed;
            do {seed = data.a.domain.randomElement(rnd);} while (seed.isZero());
            MultivariatePolynomial<BigInteger> skeleton = data.gcd.evaluate(variable, seed);

            for (int i = 0; i < 10; i++) {
                int rndSeed = i^n;
                rnd.setSeed(rndSeed);
                SparseInterpolation<BigInteger> sparseInterpolation
                        = createInterpolation(variable, data.aGCD, data.bGCD, skeleton, rnd);

                lSparseInterpolation lSparseInterpolation
                        = createInterpolation(variable, asLongPolyZp(data.aGCD), asLongPolyZp(data.bGCD), asLongPolyZp(skeleton), rnd);
                BigInteger point = data.a.domain.randomElement(rnd);
                try {
                    MultivariatePolynomial<BigInteger> expected = data.gcd.evaluate(variable, point).monic();
                    MultivariatePolynomial<BigInteger> actual = sparseInterpolation.evaluate(point).monic();
                    assertEquals(expected, actual);

                    lMultivariatePolynomial lActual = lSparseInterpolation.evaluate(point.longValueExact()).monic();
                    assertEquals(asLongPolyZp(expected), lActual);
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

        DescriptiveStatistics zippel = new DescriptiveStatistics(), brown = new DescriptiveStatistics();

        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10) {
                zippel.clear(); brown.clear();
            }
            if (n % 100 == 0) System.out.println(n);

            GCDTriplet data = sampleData.nextSample(false, true);
            MultivariatePolynomial<BigInteger> gcdZippel = null, gcdBrown = null;
            try {
                PrivateRandom.getRandom().setSeed(n);
                long start = System.nanoTime();
                gcdZippel = ZippelGCD(data.aGCD, data.bGCD);
                zippel.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                gcdBrown = BrownGCD(data.aGCD, data.bGCD);
                brown.addValue(System.nanoTime() - start);

                assertTrue(dividesQ(gcdZippel, data.gcd));
                assertTrue(dividesQ(gcdBrown, data.gcd));

                PrivateRandom.getRandom().setSeed(n);
                lMultivariatePolynomial lGcdZippel = ZippelGCD(asLongPolyZp(data.aGCD), asLongPolyZp(data.bGCD));
                assertTrue(dividesQ(lGcdZippel, asLongPolyZp(data.gcd)));
            } catch (Throwable e) {
                System.out.println("rnd seed : " + n);
                System.out.println("domain: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("gcd : " + data.gcd);
                System.out.println("gcdActual : " + gcdZippel);
                throw e;
            }
        }
        System.out.println("Zippel: " + TimeUnits.statisticsNanotime(zippel));
        System.out.println("Brown: " + TimeUnits.statisticsNanotime(brown));
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

        DescriptiveStatistics zippel = new DescriptiveStatistics(), brown = new DescriptiveStatistics();

        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10) {
                zippel.clear(); brown.clear();
            }
            if (n % 10 == 0)
                System.out.println(n);

            GCDTriplet data = sampleData.nextSample(false, false);
            MultivariatePolynomial<BigInteger> gcdZippel = null, gcdBrown = null;
            try {
                PrivateRandom.getRandom().setSeed(n);
                long start = System.nanoTime();
                gcdZippel = ZippelGCD(data.aGCD, data.bGCD);
                zippel.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                gcdBrown = BrownGCD(data.aGCD, data.bGCD);
                brown.addValue(System.nanoTime() - start);

                assertTrue(dividesQ(gcdZippel, data.gcd));
                assertTrue(dividesQ(gcdBrown, data.gcd));

                PrivateRandom.getRandom().setSeed(n);
                lMultivariatePolynomial lGcdZippel = ZippelGCD(asLongPolyZp(data.aGCD), asLongPolyZp(data.bGCD));
                assertTrue(dividesQ(lGcdZippel, asLongPolyZp(data.gcd)));
            } catch (Throwable e) {
                System.out.println("rnd seed : " + n);
                System.out.println("domain: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("gcd : " + data.gcd);
                System.out.println("gcdActual : " + gcdZippel);
                throw e;
            }
        }
        System.out.println("Zippel: " + TimeUnits.statisticsNanotime(zippel));
        System.out.println("Brown: " + TimeUnits.statisticsNanotime(brown));
    }

    @Test
    public void assadasd() throws Exception {
        IntegersModulo domain = new IntegersModulo(21535757L);
        MultivariatePolynomial<BigInteger> a = parse("13659400*b^3*c*d+6829700*a*b^3*c^3*d+3855362*a^2*b^3*c^3*d^3+2974338*a^3*b^3*c^2*d^3", domain, LEX);
        MultivariatePolynomial<BigInteger> b = parse("6107385*b^3*c*d+13821571*a*b^3*c^3*d+8143180*a^2*b^3*c^3*d^3+5678391*a^3*b^3*c^2*d^3", domain, LEX);
        System.out.println(a.monic());
        System.out.println(b.monic());
    }

    @Test
    public void testZippel_nonmonic_random3() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        RandomGenerator rnd = getRandom();
        int nVarsMin = 7, nVarsMax = 10, minDegree = 7, maxDegree = 10, minSize = 7, maxSize = 10;
        int nIterations = its(100, 1500);

        DescriptiveStatistics zippel = new DescriptiveStatistics();

        TripletPort sampleData = new TripletPort(nVarsMin, nVarsMax, minDegree, maxDegree, minSize, maxSize, rnd);
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                zippel.clear();
            if (n % 10 == 0)
                System.out.println(n);

            GCDTriplet data = sampleData.nextSample(false, false);
            MultivariatePolynomial<BigInteger> gcdZippel = null, gcdBrown = null;
            try {
                PrivateRandom.getRandom().setSeed(n);
                long start = System.nanoTime();
                gcdZippel = ZippelGCD(data.aGCD, data.bGCD);
                zippel.addValue(System.nanoTime() - start);
                assertTrue(dividesQ(gcdZippel, data.gcd));

                PrivateRandom.getRandom().setSeed(n);
                lMultivariatePolynomial lGcdZippel = ZippelGCD(asLongPolyZp(data.aGCD), asLongPolyZp(data.bGCD));
                assertTrue(dividesQ(lGcdZippel, asLongPolyZp(data.gcd)));
            } catch (Throwable e) {
                System.out.println("rnd seed : " + n);
                System.out.println("domain: " + data.domain);
                System.out.println("a: " + data.a);
                System.out.println("b: " + data.b);
                System.out.println("gcd : " + data.gcd);
                System.out.println("gcdActual : " + gcdZippel);
                throw e;
            }
        }
        System.out.println("Zippel: " + TimeUnits.statisticsNanotime(zippel));
    }

    @Test
    public void testZippel3() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("5*a^2*c^2+5*a^2*b^2*c^2+5*a^2*b^4*c^3+9*a^2*b^5*c^5+25709547*a^3*b^6*c^6+8*a^4*b*c^3+a^4*b^3*c+5*a^4*b^3*c^6", domain, LEX, vars),
                b = parse("3*a*b^2*c^2+2*a^2*b^4+25709540*a^4*b*c^6+7*a^5*c^2+8*a^6*b*c^3", domain, LEX, vars),
                gcd = parse("a + 5*b^2*c^6+2*a^4*b^4*c^5+25709543*a^5*b^2*c^5+9*a^6*c+25709540*a^6*c^3", domain, LEX, vars);
        RandomGenerator rnd = getRandom();

        int variable = a.nVariables - 1;
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        System.out.println(ZippelGCD(a, b));
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

        lMultivariatePolynomial
                aL = asLongPolyZp(a),
                bL = asLongPolyZp(b);

//        System.out.println(a);
//        System.out.println(b);
//        System.out.println(RandomMultivariatePolynomial.randomPolynomial(5, 50, 10, BigInteger.valueOf(1000), domain, LEX, getRandom()));
//        System.out.println(RandomMultivariatePolynomial.randomPolynomial(5, 50, 10, BigInteger.valueOf(1000), domain, LEX, getRandom()));
//        System.out.println(RandomMultivariatePolynomial.randomPolynomial(5, 50, 10, BigInteger.valueOf(1000), domain, LEX, getRandom()));

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertEquals(10, ZippelGCD(aL, bL).size());
            System.out.println(System.nanoTime() - start);
        }

    }


    private static final class GCDTriplet {
        final MultivariatePolynomial<BigInteger> a, b, gcd, aGCD, bGCD;
        final Domain<BigInteger> domain;

        public GCDTriplet(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, MultivariatePolynomial<BigInteger> gcd) {
            this.a = a;
            this.b = b;
            this.gcd = gcd;
            this.domain = a.domain;
            this.aGCD = a.clone().multiply(gcd);
            this.bGCD = b.clone().multiply(gcd);
        }

        private GCDTriplet(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b, MultivariatePolynomial<BigInteger> gcd, MultivariatePolynomial<BigInteger> aGCD, MultivariatePolynomial<BigInteger> bGCD, Domain<BigInteger> domain) {
            this.a = a;
            this.b = b;
            this.gcd = gcd;
            this.aGCD = aGCD;
            this.bGCD = bGCD;
            this.domain = domain;
        }

        public GCDTriplet evaluate(int variable, BigInteger value) {
            return new GCDTriplet(a.evaluate(variable, value), b.evaluate(variable, value), gcd.evaluate(variable, value),
                    aGCD.evaluate(variable, value), bGCD.evaluate(variable, value), domain);
        }
    }

    private static final class TripletPort {
        final int nVarsMin, nVarsMax,
                minDegree, maxDegree,
                minSize, maxSize;
        final RandomGenerator rnd;
        final RandomDataGenerator rndd;

        public TripletPort(int nVarsMin, int nVarsMax, int minDegree, int maxDegree, int minSize, int maxSize, RandomGenerator rnd) {
            this.nVarsMin = nVarsMin;
            this.nVarsMax = nVarsMax;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
            this.minSize = minSize;
            this.maxSize = maxSize;
            this.rnd = rnd;
            this.rndd = new RandomDataGenerator(rnd);
        }

        long counter = 0;

        public GCDTriplet nextSample(boolean primitive, boolean monic) {
            PrivateRandom.getRandom().setSeed(counter++);
            long modulus = getModulusRandom(25);
            IntegersModulo domain = new IntegersModulo(modulus);
            BigInteger bound = BigInteger.valueOf(10);

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    b = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    gcd = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd);


            if (primitive) {
                a = asNormalMultivariate(a.asOverUnivariate(0).primitivePart(), 0);
                b = asNormalMultivariate(b.asOverUnivariate(0).primitivePart(), 0);
                gcd = asNormalMultivariate(gcd.asOverUnivariate(0).primitivePart(), 0);
            }
            if (monic) {
                a.add(new MonomialTerm<>(a.nVariables, 0, a.degree() + 1, BigInteger.ONE));
                b.add(new MonomialTerm<>(a.nVariables, 0, b.degree() + 1, BigInteger.ONE));
                gcd.add(new MonomialTerm<>(a.nVariables, 0, gcd.degree() + 1, BigInteger.ONE));
            }
            return new GCDTriplet(a, b, gcd);
        }
    }

    @SuppressWarnings("unchecked")
    private static void checkConsistency(MultivariatePolynomial... polys) {
        Arrays.stream(polys).forEach(MultivariateGCDTest::checkConsistency);
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

//    @Test
//    public void trash() throws Exception {
//        String[] vars = {"a", "b", "c"};
//        IntegersModulo domain = new IntegersModulo(BigPrimes.nextPrime(56423421232L));
//        MultivariatePolynomial<BigInteger>
//                a = parse("5*a^1123*c^2+5*a^2*b^2*c^2+5*a^2*b^4*c^3+9*a^2213*b^523*c^5+25709547*a^3*b^6*c^611+8*a^4*b*c^3+a^4*b^3*c+5*a^4*b^3*c^6+a^1500", domain, GREVLEX, vars),
//                b = parse("3*a*b^2*c^2+2*a^2*b^421+25709540*a^4*b*c^6+7*a^5*c^1232+8*a^6*b^876*c^3+a^1500", domain, GREVLEX, vars),
//                gcd = parse("5*b^2*c^6+2*a^412*b^4*c^5+25709543*a^5*b^892*c^512+9*a^6*c+25709540*a^6*c^3+a^1500", domain, GREVLEX, vars);
//        RandomGenerator rnd = getRandom();
//
//        int variable = a.nVariables - 1;
//        a = fromZp(convertZp(a, 0).primitivePart(), domain, 0);
//        b = fromZp(convertZp(b, 0).primitivePart(), domain, 0);
//        gcd = fromZp(convertZp(gcd, 0).primitivePart(), domain, 0);
//
//        gcd = gcd.monic();
//        a = a.clone().monic().multiply(gcd);
//        b = b.clone().monic().multiply(gcd);
//
//
//        System.out.println(a);
//        System.out.println(b);
//        System.out.println(domain.modulus);
//        System.out.println(Zippel(a, b));
//    }
}