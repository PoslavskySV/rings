package cc.r2.core.poly.multivar;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.FactorDecompositionTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.test.Benchmark;
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
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
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
        for (int i = 0; i < its(100, 1000); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
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
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization6() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(27239);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+10435*c^2+21950*a*b^2*c^2+17887*a*b^3*c+4648*a^2*c+862*a^2*b*c", domain),
                b = lMultivariatePolynomialZp.parse("1+21170*b^2*c+7162*b^3+18183*a^2*b^2*c+16794*a^3*b+3096*a^3*b^3*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization7() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(63185123);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("7962775*c^3+54287330*b^3+48565396*a^2+26248711*a^3*b^3+10971203*a^3*b^3*c", domain),
                b = lMultivariatePolynomialZp.parse("1+48442198*b^2+36965231*b^3+35212338*a*b^2*c^3+62918195*a^2*b^2*c+47759030*a^3*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization8() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(829657);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+81633*a+270565*a^2*b*c+799187*a^2*b^2+159093*a^3*b+562717*a^3*b^3*c", domain),
                b = lMultivariatePolynomialZp.parse("1+73615*a*b^2+92694*a*b^2*c^3+582676*a^3*b*c^3+144867*a^3*b^2*c^2+132332*a^3*b^2*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testsMultivariateFactorization9() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+1179031*b^3+548360*a*b*c^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b", domain),
                b = lMultivariatePolynomialZp.parse("439433*b*c+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Benchmark
    @Test
    public void testsMultivariateFactorization10() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+1179031*b^3+548360*a*b*c^2*e^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b*d^3", domain, vars),
                b = lMultivariatePolynomialZp.parse("439433*b*c*d*e+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b", domain, vars),
                c = lMultivariatePolynomialZp.parse("439433*d*c+197065*d*e+264505*a*c^3*d+1075508*a*d^3*e^3+1338483*a^15*e + b^2", domain, vars),
                d = lMultivariatePolynomialZp.parse("433*d^2*c+165*d*e+265*a^4*c^3*d+107*a*d^3*b+1338*a^15*e +a^15*b + b^2", domain, vars),
                e = lMultivariatePolynomialZp.parse("433*d^2*e+165*d*e+265*b^4*c^3*d+107*c*d^3*a+1338*a^15*e +a^5*e + c^2", domain, vars),
                base = a.clone().multiply(b, c, d, e);
        System.out.println(base);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 10); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorSquareFree(base, false);
            Assert.assertEquals(5, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

//    Domain: Z/829657
//    Polynomial: 1+81633*a+73615*a*b^2+92694*a*b^2*c^3+270565*a^2*b*c+177174*a^2*b^2+417462*a^2*b^2*c^3+159093*a^3*b+582676*a^3*b*c^3+144867*a^3*b^2*c^2+132332*a^3*b^2*c^3+629593*a^3*b^3*c+50657*a^3*b^3*c^4+343478*a^3*b^4+595905*a^3*b^4*c^3+524441*a^4*b*c^3+826590*a^4*b^2*c^2+524016*a^4*b^2*c^3+192983*a^4*b^3+643024*a^4*b^3*c^3+467602*a^4*b^5*c+783665*a^4*b^5*c^4+308800*a^5*b^2*c^4+146627*a^5*b^3*c^3+559745*a^5*b^3*c^4+507407*a^5*b^4*c^2+806637*a^5*b^4*c^3+436944*a^6*b^2*c^3+283828*a^6*b^3*c^2+548501*a^6*b^3*c^3+414635*a^6*b^4*c^4+345447*a^6*b^5*c^3+431666*a^6*b^5*c^4
//    Expected factorization: [1+81633*a+270565*a^2*b*c+799187*a^2*b^2+159093*a^3*b+562717*a^3*b^3*c, 1+73615*a*b^2+92694*a*b^2*c^3+582676*a^3*b*c^3+144867*a^3*b^2*c^2+132332*a^3*b^2*c^3]

}
