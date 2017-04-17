package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.poly.generics.ModularDomain;
import org.junit.Test;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateGCDTest {
    @Test
    public void test1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = parse("c*b*a^2+b^2 + c"),
                b = parse("a^2+2*b^2 + 2*c"),
                gcd = parse("c*a+b+a+ c*a^3");
        ModularDomain domain = new ModularDomain(BigPrimes.nextPrime(1321323));
        a = a.multiply(gcd).setDomain(domain);
        b = b.multiply(gcd).setDomain(domain);

        System.out.println(a);
        System.out.println(b);
        System.out.println(MultivariateGCD.modularZpGCD(a, b));
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

        System.out.println(MultivariateGCD.modularZpGCD(a, b));
    }
}