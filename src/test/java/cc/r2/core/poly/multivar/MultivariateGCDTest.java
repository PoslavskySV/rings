package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.Map;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.LEX;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateGCDTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c"),
                b = parse("a^2+2*b^2 + 2*c"),
                gcd = parse("c*a+b+a+ c*a^3");
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(1321323));
        a = a.multiply(gcd).setDomain(domain);
        b = b.multiply(gcd).setDomain(domain);
        Assert.assertEquals(gcd, MultivariateGCD.denseModularGCD(a, b));
    }

    @Test
    public void test2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c"),
                b = parse("a^2+2*b^2 + 2*c"),
                gcd = parse("c*a*b+b*b+a*b+ c*a^3*b");
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(1321323));
        a = a.multiply(gcd).setDomain(domain);
        b = b.multiply(gcd).setDomain(domain);
        Assert.assertEquals(gcd, MultivariateGCD.denseModularGCD(a, b));
    }


    @Test
    public void test3() throws Exception {
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(659));
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b+655*a*b^3*c^2+6*a^2*b^3*c^2", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        Assert.assertEquals(gcd.monic(), MultivariateGCD.denseModularGCD(a, b).monic());
    }


    @Test
    public void test3a() throws Exception {
        ModularDomain domain = new ModularDomain(17);
        MultivariatePolynomial<BigInteger>
                a = parse("656*c^2+7*b^3*c+656*a*b^2+2*a^3*c+5*a^3*b^3", domain, LEX),
                b = parse("654+654*a*b^2*c^2+a*b^3*c^2+652*a^2*b*c^2+656*a^2*b^2*c", domain, LEX),
                gcd = parse("7*b^6+655*a*b^3*c^6+6*a^2*b^3*c^4", domain, LEX);
        a = a.multiply(gcd);
        b = b.multiply(gcd);
        Assert.assertEquals(gcd.monic(), MultivariateGCD.denseModularGCD(a, b).monic());
    }

    @Test
    public void test4() throws Exception {
        ModularDomain domain = new ModularDomain(653);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*b^5*c+a*b^3+a^2*b^2+a^2*b^2*c+a^3*b*c^3", domain, LEX, vars),
                b = parse("9*a*b^2*c^6+a*b^4*c^6+a^2*b^2*c^3+a^5*b^2+a^5*b^6*c^4+a^6*c+a^6*b^2*c", domain, LEX, vars),
                gcd = parse("653*b^3*c^4+b^4+b^5*c^3+a^2*b*c^2+a^4*b^2*c^4", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
        Assert.assertEquals(gcd.monic(), gcdActual.monic());
    }

    @Test
    public void test5() throws Exception {
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

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void test6() throws Exception {
        PrivateRandom.getRandom().setSeed(1564);
        ModularDomain domain = new ModularDomain(937);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("931*a^3*b^4*c+a^4+a^4*b^6*c^2+a^5*b*c^3+a^6*b*c^2", domain, LEX, vars),
                b = parse("932*b*c+a*b^6*c^2+a^3*b*c^2+a^3*b^3*c^5+a^3*b^5*c+a^5*b^5*c^3+a^6*b^2*c^6+a^6*b^4*c^5+a^6*b^6", domain, LEX, vars),
                gcd = parse("935*c^2+c^4+a^3*b*c^5+a^3*b^2*c^3+a^4*b^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void test7() throws Exception {
        PrivateRandom.getRandom().setSeed(2369);
        ModularDomain domain = new ModularDomain(569);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                a = parse("563*b^2*c+6*b^4*c^3+4*a*b^4*c^3+563*a*b^4*c^4+560*a^2*b^5*c^2+9*a^3*b^4*c+5*a^4*c^2+7*a^4*b^3*c^5+4*a^5*b^4*c^5+6*a^5*b^5", domain, LEX, vars),
                b = parse("4*b^2*c+5*b^2*c^3+5*b^3*c+3*a^2*b+3*a^2*b*c^2+565*a^3*b*c^2", domain, LEX, vars),
                gcd = parse("4+8*b^2*c^3+4*b^3+8*b^3*c+7*a*c+a*b*c+7*a^2*b^2+2*a^2*b^2*c^2+5*a^3*c^2+5*a^3*c^3", domain, LEX, vars);
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a, b);
        assertTrue(dividesQ(gcdActual, gcd));
    }

    @Test
    public void testDenseGCD1_random() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);
        RandomDataGenerator rndd = getRandomData();

        int nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        for (int n = 0; n < nIterations; n++) {
            PrivateRandom.getRandom().setSeed(n);
            long modulus = getModulusRandom(10);
            ModularDomain domain = new ModularDomain(modulus);
            BigInteger bound = BigInteger.valueOf(10);
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    b = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    gcd = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    aGCD = a.clone().multiply(gcd),
                    bGCD = b.clone().multiply(gcd),
                    gcdActual = null;

            checkConsistency(a, b, gcd, aGCD, bGCD);
            try {
                PrivateRandom.getRandom().setSeed(n);
                gcdActual = MultivariateGCD.denseModularGCD(aGCD, bGCD);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(gcdActual, gcd));
            } catch (Throwable err) {
                System.out.println("seed: " + n);
                System.out.println("modulus: " + modulus);
                System.out.println("a: " + a);
                System.out.println("b: " + b);
                System.out.println("aGCD: " + aGCD);
                System.out.println("bGCD: " + bGCD);
                System.out.println("expected gcd: " + gcd);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
    }

    @Test
    public void testDenseGCD2_random() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(123);
        RandomDataGenerator rndd = getRandomData();

        int nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(1000, 10000);
        for (int n = 0; n < nIterations; n++) {
            PrivateRandom.getRandom().setSeed(n);
            long modulus = getModulusRandom(10);
            ModularDomain domain = new ModularDomain(modulus);
            BigInteger bound = BigInteger.valueOf(10);
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    b = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    gcdActual = null;

            checkConsistency(a, b);
            try {
                PrivateRandom.getRandom().setSeed(n);
                gcdActual = MultivariateGCD.denseModularGCD(a, b);
                checkConsistency(gcdActual);
                assertTrue(dividesQ(a, gcdActual));
                assertTrue(dividesQ(b, gcdActual));
            } catch (Throwable err) {
                System.out.println("seed: " + n);
                System.out.println("modulus: " + modulus);
                System.out.println("a: " + a);
                System.out.println("b: " + b);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
    }

    @SuppressWarnings("unchecked")
    private static void checkConsistency(MultivariatePolynomial... polys) {
        Arrays.stream(polys).forEach(MultivariateGCDTest::checkConsistency);
    }

    private static <E> void checkConsistency(MultivariatePolynomial<E> poly) {
        Domain<E> domain = poly.domain;
        for (Map.Entry<MultivariatePolynomial.DegreeVector, E> e : poly.data.entrySet()) {
            E value = e.getValue();
            assertFalse(domain.isZero(value));
            assertTrue(value == domain.valueOf(value));
            if (domain instanceof ModularDomain) {
                assertTrue(domain.signum(value) > 0);
                assertTrue(((BigInteger) value).compareTo(((ModularDomain) domain).modulus) <= 0);
            }
        }
    }
}