package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar.MultivariateGCD.*;
import cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.Map;

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
    @Test
    public void testBrown1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c"),
                b = parse("a^2+2*b^2 + 2*c"),
                gcd = parse("c*a+b+a+ c*a^3");
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(1321323));
        a = a.multiply(gcd).setDomain(domain);
        b = b.multiply(gcd).setDomain(domain);
        Assert.assertEquals(gcd, MultivariateGCD.BrownGCD(a, b));
    }

    @Test
    public void testBrown2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c"),
                b = parse("a^2+2*b^2 + 2*c"),
                gcd = parse("c*a*b+b*b+a*b+ c*a^3*b");
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(1321323));
        a = a.multiply(gcd).setDomain(domain);
        b = b.multiply(gcd).setDomain(domain);
        Assert.assertEquals(gcd, MultivariateGCD.BrownGCD(a, b));
    }


    @Test
    public void testBrown3() throws Exception {
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(659));
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b+655*a*b^3*c^2+6*a^2*b^3*c^2", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        Assert.assertEquals(gcd.monic(), MultivariateGCD.BrownGCD(a, b).monic());
    }


    @Test
    public void testBrown3a() throws Exception {
        ModularDomain domain = new ModularDomain(17);
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b^6+655*a*b^3*c^6+6*a^2*b^3*c^4", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        Assert.assertEquals(gcd.monic(), MultivariateGCD.BrownGCD(a, b).monic());
    }

    @Test
    public void testBrown4() throws Exception {
        ModularDomain domain = new ModularDomain(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*b^5*c+a*b^3+a^2*b^2+a^2*b^2*c+a^3*b*c^3", domain, LEX, vars),
                b = parse("9*a*b^2*c^6+a*b^4*c^6+a^2*b^2*c^3+a^5*b^2+a^5*b^6*c^4+a^6*c+a^6*b^2*c", domain, LEX, vars),
                gcd = parse("653*b^3*c^4+b^4+b^5*c^3+a^2*b*c^2+a^4*b^2*c^4", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
        Assert.assertEquals(gcd.monic(), gcdActual.monic());
    }

    @Test
    public void testBrown5() throws Exception {
        RandomGenerator rnd = PrivateRandom.getRandom();
        rnd.setSeed(28);
        ModularDomain domain = new ModularDomain(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("561*a^2*c^2+a^2*b^2*c+a^3*b+a^4*b^2+a^4*b^5*c^3+a^5*b", domain, LEX, vars),
                b = parse("561*a*c^3+a*b^4*c^5+a^2*c^2+a^2*b^6*c^3+a^3*b^6*c^5+a^5*b^5*c^3+a^5*b^5*c^6", domain, LEX, vars),
                gcd = parse("4*c^2+b^4+a^2*b^4*c+a^3*b^2*c+a^3*b^6+a^5*b^2*c^6+a^6*b^5", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void testBrown6() throws Exception {
        PrivateRandom.getRandom().setSeed(1564);
        ModularDomain domain = new ModularDomain(937);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("931*a^3*b^4*c+a^4+a^4*b^6*c^2+a^5*b*c^3+a^6*b*c^2", domain, LEX, vars),
                b = parse("932*b*c+a*b^6*c^2+a^3*b*c^2+a^3*b^3*c^5+a^3*b^5*c+a^5*b^5*c^3+a^6*b^2*c^6+a^6*b^4*c^5+a^6*b^6", domain, LEX, vars),
                gcd = parse("935*c^2+c^4+a^3*b*c^5+a^3*b^2*c^3+a^4*b^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void testBrown7() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        ModularDomain domain = new ModularDomain(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("563*b^2*c+6*b^4*c^3+4*a*b^4*c^3+563*a*b^4*c^4+560*a^2*b^5*c^2+9*a^3*b^4*c+5*a^4*c^2+7*a^4*b^3*c^5+4*a^5*b^4*c^5+6*a^5*b^5", domain, LEX, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, LEX, vars),
                gcd = parse("4+8*b^2*c^3+4*b^3+8*b^3*c+7*a*c+a*b*c+7*a^2*b^2+2*a^2*b^2*c^2+5*a^3*c^2+5*a^3*c^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
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
                gcdActual = MultivariateGCD.BrownGCD(data.aGCD, data.bGCD);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(gcdActual, data.gcd));
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
                gcdActual = MultivariateGCD.BrownGCD(data.a, data.b);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(data.a, gcdActual));
                assertTrue(dividesQ(data.b, gcdActual));
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
            } while (gcd.usedVariables() > nVars);

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
                gcdActual = MultivariateGCD.BrownGCD(aGCD, bGCD);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(gcdActual, gcd));
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
        ModularDomain domain = new ModularDomain(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, LEX, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, LEX, vars),
                gcd = parse("a^2 + c^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }


    @Test
    public void testBrown9() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        ModularDomain domain = new ModularDomain(569);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + c^4", domain, LEX, vars),
                b = parse("b^2 + d^2", domain, LEX, vars),
                gcd = parse("d^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void testBrown10() throws Exception {
        ModularDomain domain = new ModularDomain(5642359);
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger>
                a = parse("1199884 + 4783764*b + b^2 + 3215597*a*b + 2196297*a*b^2 + 4781733*a^4 + a^4*b + 2196297*a^5*b", domain, LEX, vars),
                b = parse("4645946 + 3921107*b + b^2 + 5605437*a*b + 2196297*a*b^2 + 4781733*a^3 + a^3*b + 2196297*a^4*b", domain, LEX, vars);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.BrownGCD(a, b);
        gcdActual = gcdActual.monic().multiply(domain.valueOf(4781733));
        MultivariatePolynomial<BigInteger> expected = parse("1574588 + 4559668*b + 4781733*a*b", domain, LEX, vars);
        assertEquals(expected, gcdActual);
    }

    @Test
    public void testZippel1() throws Exception {
        String[] vars = {"a", "b", "c"};
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(5642342L));
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
            SparseInterpolation sparseInterpolation
                    = createInterpolation(variable, a, b, skeleton, rnd);
            BigInteger point = domain.randomElement(rnd);
            assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
        }
    }

    @Test
    public void testZippel2() throws Exception {
        String[] vars = {"a", "b", "c"};
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^2 + b^2 + a*c^2", domain, LEX, vars),
                b = parse("a^2 + 2*b^2 + b*c^2", domain, LEX, vars),
                gcd = parse("a^2 + c*a + b + a*c*b + a*c^2", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = ZippelGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void testZippel5() throws Exception {
        String[] vars = {"a", "b", "c"};
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(31579447));
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


        SparseInterpolation sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
        BigInteger point = domain.valueOf(1324);
        assertEquals(gcd.evaluate(variable, point), sparseInterpolation.evaluate(point));
    }

    @Test
    public void testZippel6() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(31579447));
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

        SparseInterpolation sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
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
                SparseInterpolation sparseInterpolation
                        = createInterpolation(variable, data.aGCD, data.bGCD, skeleton, rnd);
                BigInteger point = data.a.domain.randomElement(rnd);
                try {
                    MultivariatePolynomial<BigInteger> expected = data.gcd.evaluate(variable, point);
                    MultivariatePolynomial<BigInteger> actual = sparseInterpolation.evaluate(point);
                    expected = expected.monic();
                    actual = actual.monic();
                    assertEquals(expected, actual);
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
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(31579447));
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

        SparseInterpolation sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
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
            System.out.println(n);

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
    public void testZippel8() throws Exception {
        MultivariateGCD.ALWAYS_LINZIP = true;
        String[] vars = {"a", "b", "c"};
        ModularDomain domain = new ModularDomain(31579447);
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

        SparseInterpolation sparseInterpolation = createInterpolation(variable, a, b, skeleton, rnd);
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
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(5642342L));
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

    @Test
    public void testZippel4() throws Exception {
        PrivateRandom.getRandom().setSeed(12323);
        String[] vars = {"a", "b", "c", "d"};
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(5642342L));
        MultivariatePolynomial<BigInteger>
                a = parse("a^4 + b + c + d", domain, LEX, vars),
                b = parse("a^3 + c + d", domain, LEX, vars),
                gcd = parse("a * b*c*d + c + d", domain, LEX, vars);

        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);


        System.out.println(ZippelGCD(a, b));
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
            ModularDomain domain = new ModularDomain(modulus);
            BigInteger bound = BigInteger.valueOf(10);

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    b = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    gcd = randomPolynomial(nVariables, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd);


            if (primitive) {
                a = fromZp(convertZp(a, 0).primitivePart(), domain, 0);
                b = fromZp(convertZp(b, 0).primitivePart(), domain, 0);
                gcd = fromZp(convertZp(gcd, 0).primitivePart(), domain, 0);
            }
            if (monic) {
                a.add(new DegreeVector(a.nVariables, 0, a.degree() + 1), BigInteger.ONE);
                b.add(new DegreeVector(a.nVariables, 0, b.degree() + 1), BigInteger.ONE);
                gcd.add(new DegreeVector(a.nVariables, 0, gcd.degree() + 1), BigInteger.ONE);
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
        for (Map.Entry<DegreeVector, E> e : poly.data.entrySet()) {
            E value = e.getValue();
            assertFalse(domain.isZero(value));
            assertTrue(value == domain.valueOf(value));
            if (domain instanceof ModularDomain) {
                assertTrue(domain.signum(value) > 0);
                assertTrue(((BigInteger) value).compareTo(((ModularDomain) domain).modulus) <= 0);
            }
        }
    }

//    @Test
//    public void trash() throws Exception {
//        String[] vars = {"a", "b", "c"};
//        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(56423421232L));
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