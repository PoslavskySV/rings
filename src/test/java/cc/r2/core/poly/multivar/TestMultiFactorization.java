package cc.r2.core.poly.multivar;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.FactorDecompositionTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.util.TimeUnits;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class TestMultiFactorization extends AbstractPolynomialTest {

    @Test
    public void testsMultivariateFactorization1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(1000, 1000); i++) {
            PrivateRandom.getRandom().setSeed(i);
//            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization1a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(42);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testsMultivariateFactorization1b() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(328);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testsMultivariateFactorization2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 1000); i++) {
            PrivateRandom.getRandom().setSeed(i);
//            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 1000); i++) {
            PrivateRandom.getRandom().setSeed(i);
//            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(1000, 1000); i++) {
//            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
//            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization4a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(31);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testsMultivariateFactorization4b() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(88);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
//        System.out.println(factors);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testsMultivariateFactorization4c() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(531);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
//        System.out.println(factors);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testsMultivariateFactorization5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c^2*b^2 + 11*c + b^3 + 1", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b*c^3 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(1000, 1000); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }
}
