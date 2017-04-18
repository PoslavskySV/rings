package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.generics.ModularDomain;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.multivar.DivisionWithRemainderMultivariate.dividesQ;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.LEX;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;

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
        for (int i = 0; i < 1000; i++) {

            ModularDomain domain = new ModularDomain(653);
            String[] vars = {"a", "b", "c"};
            MultivariatePolynomial<BigInteger>
                    a = parse("6*b^5*c+a*b^3+a^2*b^2+a^2*b^2*c+a^3*b*c^3", domain, LEX, vars),
                    b = parse("9*a*b^2*c^6+a*b^4*c^6+a^2*b^2*c^3+a^5*b^2+a^5*b^6*c^4+a^6*c+a^6*b^2*c", domain, LEX, vars),
                    gcd = parse("653*b^3*c^4+b^4+b^5*c^3+a^2*b*c^2+a^4*b^2*c^4", domain, LEX, vars);
            a = a.clone().multiply(gcd);
            b = b.clone().multiply(gcd);

            MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a, b);
            Assert.assertTrue(dividesQ(gcdActual, gcd));
//            Assert.assertEquals(gcd.monic(), gcdActual.monic());

        }
    }


    @Test
    public void test_PGCD_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nVars = 3;
        int minDegree = 3, maxDegree = 5;
        int minSize = 5, maxSize = 10;

        int nIterations = its(100, 1000);
        for (int n = 0; n < nIterations; n++) {
            long modulus = getModulusRandom(10);
            ModularDomain domain = new ModularDomain(modulus);
            BigInteger bound = BigInteger.valueOf(10);
            MultivariatePolynomial<BigInteger>
                    a = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    b = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd),
                    gcd = randomPolynomial(nVars, rndd.nextInt(minDegree, maxDegree), rndd.nextInt(minSize, maxSize), bound, domain, LEX, rnd);

            MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.denseModularGCD(a.clone().multiply(gcd), b.clone().multiply(gcd));
            try {
                Assert.assertTrue(dividesQ(gcdActual, gcd));
            } catch (AssertionError err) {
                System.out.println(modulus);
                System.out.println(a);
                System.out.println(b);
                System.out.println("expected gcd: " + gcd);
                System.out.println("actual gcd  : " + gcdActual);
                throw err;
            }
        }
    }
}