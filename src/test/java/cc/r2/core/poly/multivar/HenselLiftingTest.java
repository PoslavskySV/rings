package cc.r2.core.poly.multivar;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.lIntegersModulo;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import static cc.r2.core.poly.multivar.HenselLifting.modImage;
import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
@Ignore
public class HenselLiftingTest {
    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^2*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("1 + x - x^2 + x^2*y", domain, vars),
                b = parse("2 + x + x^2 + 2*y*x^2 + x^2*y^5 + y", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 + a^15*b^2*c^3 + a^15*c + 7*a^15- 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5  + a^5*b*c - a^5*b + 2*a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();


        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }
}