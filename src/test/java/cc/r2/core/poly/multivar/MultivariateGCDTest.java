package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.generics.ModularDomain;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Comparator;

import static cc.r2.core.poly.multivar.DivisionWithRemainderMultivariate.dividesQ;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.GREVLEX;
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
        Assert.assertEquals(gcd, MultivariateGCD.modularZpGCD(a, b));
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
        Assert.assertEquals(gcd, MultivariateGCD.modularZpGCD(a, b));
    }

    @Test
    public void test3() throws Exception {
        long modulus = 881;
        String[] vars = {"a", "b"};
        ModularDomain domain = new ModularDomain(modulus);
        Comparator<MultivariatePolynomial.DegreeVector> ordering = LEX;
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b^3+4*a^2*b^3+6*a^3*b", domain, ordering, vars),
                b = parse("b^2+9*a*b^2-8*a^2*b", domain, ordering, vars),
                gcd = parse("6*a+18*a^2+7*a^3+7*a^3*b^2", domain, ordering, vars);

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        System.out.println(MultivariateGCD.modularZpGCD(a, b));
    }

    @Test
    public void test_PGCD_random() throws Exception {
        RandomGenerator rnd = getRandom();
//        rnd.setSeed(123);
        long modulus = getModulusRandom(10);
        ModularDomain domain = new ModularDomain(modulus);
        BigInteger bound = BigInteger.valueOf(10);
        int nVars = 2;
        MultivariatePolynomial<BigInteger>
                a = randomPolynomial(nVars, 3, 5, bound, domain, LEX, rnd),
                b = randomPolynomial(nVars, 3, 5, bound, domain, LEX, rnd),
                gcd = randomPolynomial(nVars, 3, 5, bound, domain, LEX, rnd);

        System.out.println(a);
        System.out.println(b);

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        System.out.println(modulus);
        System.out.println(a);
        System.out.println(b);
        System.out.println(gcd);

        MultivariatePolynomial<BigInteger> gcdActual = MultivariateGCD.modularZpGCD(a, b);
        System.out.println(gcdActual);
        Assert.assertTrue(dividesQ(gcdActual, gcd));

    }
}