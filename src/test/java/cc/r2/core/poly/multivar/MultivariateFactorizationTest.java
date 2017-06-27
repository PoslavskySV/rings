package cc.r2.core.poly.multivar;

import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.util.TimeUnits;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import static cc.r2.core.poly.multivar.MultivariateFactorization.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
@Ignore
public class MultivariateFactorizationTest {
    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(6711111));
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("b*a + b^2 + 1", domain),
                b = lMultivariatePolynomialZp.parse("a*b + 66*b + 1", domain),
                c = lMultivariatePolynomialZp.parse("b*a^2 + a^2 + b", domain),
                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);

        System.out.println(base);
        System.out.println(MultivariateFactorization.biVariateFactorSquareFree(base));
        for (int i = 0; i < 1000; i++) {
            UNI = LIFT = RECOMB = 0;
            long start = System.nanoTime();
            Assert.assertEquals(4, MultivariateFactorization.biVariateFactorSquareFree(base).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println(TimeUnits.nanosecondsToString(UNI));
            System.out.println(TimeUnits.nanosecondsToString(MultivariateFactorization.LIFT));
            System.out.println(TimeUnits.nanosecondsToString(MultivariateFactorization.RECOMB));
            System.out.println();
        }

    }
}