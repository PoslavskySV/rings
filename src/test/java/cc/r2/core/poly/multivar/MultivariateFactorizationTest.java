package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.factorization.FactorizationTestData;
import cc.r2.core.poly.factorization.FactorizationTestData.FactorizationAlgorithm;
import cc.r2.core.poly.factorization.FactorizationTestData.SampleDecomposition;
import cc.r2.core.poly.factorization.FactorizationTestData.SampleDecompositionSource;
import cc.r2.core.poly.univar.DivisionWithRemainder;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.ArraysUtil;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.AMultivariatePolynomial.renameVariables;
import static cc.r2.core.poly.multivar.MultivariateFactorization.GCDFreeBasis;
import static cc.r2.core.poly.multivar.MultivariateFactorization.bivariateDenseFactorSquareFree;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;
import static cc.r2.core.util.ArraysUtil.negate;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateFactorizationTest extends AbstractPolynomialTest {
    @Ignore
    @Test
    public void testBivariate1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(67);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("b*a^5 + a^4*b^2 + 11 + b^3", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b + 17*b^2 + 1", domain),
                c = lMultivariatePolynomialZp.parse("b^3*a^4 + a^4 + b", domain),
                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);

        System.out.println(base);
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            Assert.assertEquals(4, bivariateDenseFactorSquareFree(base).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testBivariate2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(62653);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp[] factors = {
                lMultivariatePolynomialZp.parse("17096+6578*a*b^2+54905*a^3", domain, vars),
                lMultivariatePolynomialZp.parse("43370+32368*a^2*b^2+45712*a^2*b^4+52302*a^4+23776*a^4*b^2", domain, vars)
        };
        lMultivariatePolynomialZp poly = factors[0].createOne().multiply(factors);
        FactorDecomposition<lMultivariatePolynomialZp> factorization = bivariateDenseFactorSquareFree(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }

    @Test
    public void testBivariateRandom3() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 5);
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(100, 1000),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"));
    }

    @Test
    public void testBivariateRandom4() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                5, 10,
                3, 6);
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(10, 100),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"));
    }

    @Test
    public void testBivariate5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1336151);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp[] factors = {
                lMultivariatePolynomialZp.parse("319792+402081*a^3", domain, vars),
                lMultivariatePolynomialZp.parse("685686+694157*a", domain, vars),
                lMultivariatePolynomialZp.parse("616781+1057293*b^2+158725*a+730076*a*b^2", domain, vars)
        };
        lMultivariatePolynomialZp poly = factors[0].createOne().multiply(factors);
        FactorDecomposition<lMultivariatePolynomialZp> factorization = bivariateDenseFactorSquareFree(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }

    @Test
    public void testBivariate6() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(57352861);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp[] factors = {
                lMultivariatePolynomialZp.parse("15042434+15817122*b", domain, vars),
                lMultivariatePolynomialZp.parse("39330400+51579304*a^2", domain, vars)
        };
        lMultivariatePolynomialZp poly = factors[0].createOne().multiply(factors);
        FactorDecomposition<lMultivariatePolynomialZp> factorization = bivariateDenseFactorSquareFree(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }


    @Test
    public void testBivaraiteSmallDomain7() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(2);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1 + b*a^5 + a^4*b^2 + b^3", domain, vars),
                b = lMultivariatePolynomialZp.parse("a + a^2 + b + a^6*b + 66*b + 17*b^2 + 1", domain, vars),
                c = lMultivariatePolynomialZp.parse("b^3*a^4 + a^4 + b", domain, vars),
                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain, vars),
                base = a.clone().multiply(b, c, d);
        System.out.println(base);
        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> decomposition = bivariateDenseFactorSquareFree(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(5, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }

        //49ms 35ms 30ms 28ms 38ms 38ms 28ms 28ms 33ms 33ms 33ms 34ms
    }

    @Test
    public void testBivaraiteSmallDomain5Random8() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 5);
        source.minModulusBits = 2;
        source.maxModulusBits = 3;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(50, 500),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"));
    }

    @Test
    public void testBivaraiteSmallDomain5Random9() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                5, 10,
                3, 6);
        source.minModulusBits = 2;
        source.maxModulusBits = 3;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(10, 100),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"),
                FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFree, "Bivariate dense factorization"));
    }

    @Ignore
    @Benchmark
    @Test
    public void testBivariateBenchmarkSingular() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(2);
        lMultivariatePolynomialZp poly = lMultivariatePolynomialZp.parse("x^4120 + x^4118*y^2 + x^3708*y^400 + x^3706*y^402 + x^2781*y^1300 + x^2779*y^1302 + x^1339*y^2700 + x^927*y^3100 + y^4000 + x^7172*y^4167 + x^8349*y^4432 + x^8347*y^4434 + x^6760*y^4567 + x^5833*y^5467 + x^5568*y^7132 + x^11401*y^8599", domain);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            System.out.println(bivariateDenseFactorSquareFree(poly));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start)); ;
        }
//        System.out.println(poly);
    }

    static <Poly extends IGeneralPolynomial<Poly>>
    Poly multiply(Poly... factors) {
        return factors[0].createOne().multiply(factors);
    }

    @Test
    public void testGCDFreeBasis1() throws Exception {
        long modulus = 17;
        lUnivariatePolynomialZp
                a = lUnivariatePolynomialZ.create(1, 2, 3).modulus(modulus),
                b = lUnivariatePolynomialZ.create(3, 2, 1, 2).modulus(modulus),
                c = lUnivariatePolynomialZ.create(1, 0, 0, 2, 3).modulus(modulus),
                d = lUnivariatePolynomialZ.create(1, 11, 0, 12, 4).modulus(modulus);

        FactorDecomposition<lUnivariatePolynomialZp>
                d1 = FactorDecomposition.create(Arrays.asList(multiply(a, b, b, b), multiply(a, c, c), multiply(a, a, d))),
                d2 = FactorDecomposition.create(Arrays.asList(multiply(a, c, b), multiply(b, d, c), multiply(d, c, d))),
                d3 = FactorDecomposition.create(Arrays.asList(multiply(c, c, d), multiply(b, b, c), multiply(a, b, c, d)));

        @SuppressWarnings("unchecked")
        FactorDecomposition<lUnivariatePolynomialZp>[] decomps = new FactorDecomposition[]{
                d1.clone(), d2.clone(), d3.clone()
        };

        GCDFreeBasis(decomps);

        Assert.assertEquals(d1.toPolynomial(), decomps[0].toPolynomial());
        Assert.assertEquals(d2.toPolynomial(), decomps[1].toPolynomial());
        Assert.assertEquals(d3.toPolynomial(), decomps[2].toPolynomial());

        System.out.println(d1.toPolynomial().equals(decomps[0].toPolynomial()));
        System.out.println(Arrays.toString(DivisionWithRemainder.divideAndRemainder(d1.toPolynomial(), decomps[0].toPolynomial(), true)));
        for (FactorDecomposition<lUnivariatePolynomialZp> decomp : decomps)
            System.out.println(decomp.size() + " => " + decomp);
    }

    @Test
    public void testMultivariateFactorization1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization1a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(42);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization1b() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(328);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization4a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(31);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization4b() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(88);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization4c() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = lMultivariatePolynomialZp.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(531);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("2*a^5*c^2*b^2 + 11*c + b^3 + 1", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b*c^3 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization6() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(27239);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+10435*c^2+21950*a*b^2*c^2+17887*a*b^3*c+4648*a^2*c+862*a^2*b*c", domain),
                b = lMultivariatePolynomialZp.parse("1+21170*b^2*c+7162*b^3+18183*a^2*b^2*c+16794*a^3*b+3096*a^3*b^3*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization7() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(63185123);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("7962775*c^3+54287330*b^3+48565396*a^2+26248711*a^3*b^3+10971203*a^3*b^3*c", domain),
                b = lMultivariatePolynomialZp.parse("1+48442198*b^2+36965231*b^3+35212338*a*b^2*c^3+62918195*a^2*b^2*c+47759030*a^3*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization8() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(829657);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+81633*a+270565*a^2*b*c+799187*a^2*b^2+159093*a^3*b+562717*a^3*b^3*c", domain),
                b = lMultivariatePolynomialZp.parse("1+73615*a*b^2+92694*a*b^2*c^3+582676*a^3*b*c^3+144867*a^3*b^2*c^2+132332*a^3*b^2*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization9() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+1179031*b^3+548360*a*b*c^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b", domain),
                b = lMultivariatePolynomialZp.parse("439433*b*c+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Benchmark
    @Test
    public void testMultivariateFactorization10() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                factors[] =
                {
                        lMultivariatePolynomialZp.parse("1+1179031*b^3+548360*a*b*c^2*e^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b*d^3+a^15", domain, vars),
                        lMultivariatePolynomialZp.parse("439433*b*c*d*e+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b+a^15", domain, vars),
                        lMultivariatePolynomialZp.parse("439433*d*c+197065*d*e+264505*a*c^3*d+1075508*a*d^3*e^3+1338483*a^15*e +a^15 + b^2", domain, vars),
                        lMultivariatePolynomialZp.parse("433*d^2*c+165*d*e+265*a^4*c^3*d+107*a*d^3*b+1338*a^15*e +a^15*b +a^15 + b^2", domain, vars),
                        lMultivariatePolynomialZp.parse("433*d^2*e+165*d*e+265*b^4*c^3*d+107*c*d^3*a+1338*a^15*e +a^15*e + c^2 + a^15*b + a^15 + a^15*d + a^15*e", domain, vars),
                },
                base = factors[0].createOne().multiply(factors);
        System.out.println(FactorDecomposition.create(Arrays.asList(factors)));
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorPrimitive(base, false);
            Assert.assertEquals(5, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
        //1391 ms
        //1529 ms
        //1524 ms
        //1351 ms
    }

    @Test
    public void testMultivariateFactorization10a() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                factors[] =
                {
                        lMultivariatePolynomialZp.parse("1+1179031*b^3+548360*a*b*c^2*e^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b*d^3+a^15", domain, vars),
                        lMultivariatePolynomialZp.parse("439433*b*c*d*e+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b+a^15", domain, vars),
                        lMultivariatePolynomialZp.parse("439433*d*c+197065*d*e+264505*a*c^3*d+1075508*a*d^3*e^3+1338483*a^15*e +a^15 + b^2", domain, vars),
                        lMultivariatePolynomialZp.parse("433*d^2*c+165*d*e+265*a^4*c^3*d+107*a*d^3*b+1338*a^15*e +a^15*b +a^15 + b^2", domain, vars),
                        lMultivariatePolynomialZp.parse("433*d^2*e+165*d*e+265*b^4*c^3*d+107*c*d^3*a+1338*a^15*e +a^15*e + c^2 + a^15*b + a^15 + a^15*d + a^15*e", domain, vars),
                };

        lMultivariatePolynomialZp base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(5, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization11() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("a + b + c + d", domain, vars),
                b = lMultivariatePolynomialZp.parse("b + c + d + e", domain, vars),
                c = lMultivariatePolynomialZp.parse("c + d + e + a", domain, vars),
                d = lMultivariatePolynomialZp.parse("2*d + 2*e + 2*a + 2*b", domain, vars),
                e = lMultivariatePolynomialZp.parse("e + a + b + 2*c", domain, vars),
                base = a.clone().multiply(b, c, d, e);
        FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorToPrimitive(base);
        Assert.assertEquals(5, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization12() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(74017943);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("18370804*b^3*c^3+40724543*a+25831118*a*b+28120978*a*b*c^3+49314822*a*b^3*c^3", domain),
                b = lMultivariatePolynomialZp.parse("58629076*b*c^3+37897966*a*b^3*c^2+55834047*a^2*c^3+18265939*a^3*c^3+43535405*a^3*b", domain);

        a = AMultivariatePolynomial.renameVariables(a, new int[]{1, 2, 0});
        b = AMultivariatePolynomial.renameVariables(b, new int[]{1, 2, 0});
        lMultivariatePolynomialZp base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);

        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization13() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(192149);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+81770*b+19081*b*c+23953*a*b^3*c+7807*a^3*b+14026*a^3*b^2*c^3", domain),
                b = lMultivariatePolynomialZp.parse("1+105163*a^2+81015*a^2*c+166076*a^3*c^3+106464*a^3*b^2*c^2+43621*a^3*b^3", domain);
        lMultivariatePolynomialZp base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization14() throws Exception {

        lIntegersModulo domain = new lIntegersModulo(386039);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1+377446*b*c+302126*a*b^2+97219*a*b^2*c^2+92497*a*b^2*c^3+84001*a^3*b^2", domain),
                b = lMultivariatePolynomialZp.parse("1+248663*a*b*c^2+10589*a*b^3*c^3+62097*a^2*c+81842*a^2*b^2*c+51504*a^3", domain);
        a = AMultivariatePolynomial.renameVariables(a, new int[]{0, 2, 1});
        b = AMultivariatePolynomial.renameVariables(b, new int[]{0, 2, 1});
        lMultivariatePolynomialZp base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        FactorDecomposition<lMultivariatePolynomialZp> factors = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization15() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(34957081);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                factors[] =
                {
                        lMultivariatePolynomialZp.parse("2999550*b*c^3*d^2+14700809*a*c^2+13282494*a*b^3*c*d^3+30075047*a^2*c^2*d+2736476*a^3*d^2", domain, vars),
                        lMultivariatePolynomialZp.parse("1+7381919*c^2+33667094*b*c^2*d^2+26355114*b^3*c*d^3+30536438*b^3*c^2*d^3+10734561*a*b*c", domain, vars),
                        lMultivariatePolynomialZp.parse("1+24559556*b^2*c^2*d^3+13753085*a*b*c^2*d+22081133*a*b*c^3*d^3+30781594*a^3*b*c^2+27334226*a^3*b^3*d^3", domain, vars),
                };

        for (int i = 0; i < factors.length; i++) {
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 0, 1);
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 1, 3);
        }

        lMultivariatePolynomialZp base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(3, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization16() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(316797977);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                factors[] =
                {
                        lMultivariatePolynomialZp.parse("157129769*b*c*d^3+234842760*b*c^3*d^3+105538252*a+55274980*a*b^2*c^2+89647854*a^3*b^3*c^2*d^2", domain, vars),
                        lMultivariatePolynomialZp.parse("241626121*d^2+47627151*b^2*c*d+150262012*a^2*b*c^2+299159387*a^2*b^3*c^2+53788517*a^3*b*c^3*d^3", domain, vars),
                        lMultivariatePolynomialZp.parse("1+76011411*b*d^3+189305430*b*c^3*d+218732499*a^2*b^2*c*d^2+125990992*a^2*b^3*c^2*d+36953173*a^3*b*c^2*d", domain, vars),
                        lMultivariatePolynomialZp.parse("1+299415864*a*b^2*c^3*d^3+154985851*a^2*c*d+157246866*a^2*b^2*c^3*d^3+32838497*a^3*b^3*c*d+41239905*a^3*b^3*c*d^2", domain, vars),
                };

        lMultivariatePolynomialZp base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorPrimitive(base, false);
        Assert.assertEquals(4, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization17() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(125617);
        String[] vars = {"a", "b", "c", "d"};
        lMultivariatePolynomialZp
                factors[] =
                {
                        lMultivariatePolynomialZp.parse("1+71543*b*c+89032*b*c*d+101233*b^2*c+69912*b^2*c^2*d^2+122146*a*c^2", domain, vars),
                        lMultivariatePolynomialZp.parse("1+62395*b*c^2*d^2+111331*a*b*c*d^2+13129*a^3*b^2*c^2*d^2+54277*a^3*b^3*c^2+36488*a^3*b^3*c^2*d^3", domain, vars),
                        lMultivariatePolynomialZp.parse("1+64768*a*b^2*d+66817*a*b^3*c*d+19563*a^2+13861*a^3*b*c^3+76958*a^3*b^3*d", domain, vars),
                };

        for (int i = 0; i < factors.length; i++)
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 0, 1);

        lMultivariatePolynomialZp base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        for (int i = 0; i < its(10, 10); i++) {
            //long start = System.nanoTime();
            FactorDecomposition<lMultivariatePolynomialZp> decomposition = MultivariateFactorization.factorPrimitive(base, false);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }
    
    @Test
    public void testMultivariateFactorizationRandom1() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 2,
                3, 3,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(100, 100),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"));
    }

    @Test
    public void testMultivariateFactorizationRandom2() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 4,
                3, 4,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(100, 100),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"));
    }

    @Test
    public void testMultivariateFactorizationRandom3() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 4,
                5, 6,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(100, 100),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"),
                FactorizationTestData.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitive, "Multivariate factorization in finite fields"));
    }

    /* ==================================== Test data =============================================== */

    public static void
    testFactorizationAlgorithm(SampleDecompositionSource<lMultivariatePolynomialZp> source,
                               int nIterations,
                               FactorizationAlgorithm<lMultivariatePolynomialZp> lAlgorithm,
                               FactorizationAlgorithm<MultivariatePolynomial<BigInteger>> bAlgorithm) {
        System.out.println("Testing factorization algorithm " + lAlgorithm.name);
        System.out.println("Input source: " + source);

        DescriptiveStatistics
                lTiming = new DescriptiveStatistics(),
                bTiming = new DescriptiveStatistics();

        int prevProgress = -1, currProgress;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.asList(lTiming, bTiming).forEach(DescriptiveStatistics::clear);

            if ((currProgress = (int) (100.0 * n / nIterations)) != prevProgress) {
                prevProgress = currProgress;
                System.out.print(">");
                System.out.flush();
            }
            SampleDecomposition<lMultivariatePolynomialZp> sample = source.next();
            FactorDecomposition<lMultivariatePolynomialZp> lDecomposition = null;
            FactorDecomposition<MultivariatePolynomial<BigInteger>> bDecomposition = null;
            try {
                long start;
                if (lAlgorithm != null) {
                    start = System.nanoTime();
                    lDecomposition = lAlgorithm.algorithm.apply(sample.poly);
                    lTiming.addValue(System.nanoTime() - start);
                    sample.assertFactorization(lDecomposition);
                }

                if (bAlgorithm != null) {
                    SampleDecomposition<MultivariatePolynomial<BigInteger>>
                            bSample = toBigPoly(sample);
                    start = System.nanoTime();
                    bDecomposition = bAlgorithm.algorithm.apply(bSample.poly);
                    bTiming.addValue(System.nanoTime() - start);
                    bSample.assertFactorization(bDecomposition);
                }
            } catch (Throwable throwable) {
                System.out.println("============ Error ============");
                System.out.println("Domain: " + sample.poly.domain);
                System.out.println("Polynomial: " + sample.poly);
                System.out.println("Expected factorization: " + Arrays.toString(sample.factors));
                System.out.println("Actual decomposition (longs): " + lDecomposition);
                System.out.println("Actual decomposition (BigInts): " + bDecomposition);
                throw throwable;
            }
        }

        System.out.println(source.statisticsToString());

        System.out.println("\n============ Timings ============");
        System.out.println("Machine integers: " + TimeUnits.statisticsNanotime(lTiming));
        System.out.println("Big integers    : " + TimeUnits.statisticsNanotime(bTiming));
    }

    @SuppressWarnings("unchecked")
    static SampleDecomposition<MultivariatePolynomial<BigInteger>>
    toBigPoly(SampleDecomposition<lMultivariatePolynomialZp> decomposition) {
        return new SampleDecomposition<>(
                Arrays.stream(decomposition.factors)
                        .map(lMultivariatePolynomialZp::toBigPoly)
                        .toArray(MultivariatePolynomial[]::new));
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SampleDecompositionSource<Poly> orderVarsByDegree(SampleDecompositionSource<Poly> source) {
        return new SampleDecompositionSource<Poly>() {
            @Override
            public SampleDecomposition<Poly> next0() {
                SampleDecomposition<Poly> sample = source.next0();
                int[] degrees = sample.poly.degrees();
                int[] variables = ArraysUtil.sequence(0, degrees.length);
                ArraysUtil.quickSort(negate(degrees), variables);
                return new SampleDecomposition<>(
                        Arrays.stream(sample.factors)
                                .map(p -> renameVariables(p, variables))
                                .toArray(sample.poly::arrayNewInstance));
            }

            @Override
            public String toString() {
                return super.toString() + " (vars ordered)";
            }
        };
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SampleDecompositionSource<Poly> filterMonomialContent(SampleDecompositionSource<Poly> source) {
        return new SampleDecompositionSource<Poly>() {
            @Override
            public SampleDecomposition<Poly> next0() {
                SampleDecomposition<Poly> sample = source.next0();
                for (Poly factor : sample.factors) {
                    if (!factor.monomialContent().isZeroVector())
                        factor.increment();
                    assert factor.monomialContent().isZeroVector();
                }
                return new SampleDecomposition<>(sample.factors);
            }

            @Override
            public String toString() {
                return super.toString() + " (content filtered)";
            }
        };
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SampleDecompositionSource<Poly> filterNonSquareFree(SampleDecompositionSource<Poly> source) {
        return new SampleDecompositionSource<Poly>() {
            @Override
            public SampleDecomposition<Poly> next0() {
                SampleDecomposition<Poly> sample = source.next0();
                if (MultivariateSquareFreeFactorization.isSquareFree(sample.poly))
                    return sample;
                else return next0();
            }

            @Override
            public String toString() {
                return super.toString() + " (square-free)";
            }
        };
    }

    public static final class lSampleDecompositionSource
            extends SampleDecompositionSource<lMultivariatePolynomialZp> {
        final int
                nFactorsMin, nFactorsMax,
                nVarsMin, nVarsMax,
                minSize, maxSize,
                minDegree, maxDegree;

        int minModulusBits = 10, maxModulusBits = 32;

        public lSampleDecompositionSource(int nFactorsMin, int nFactorsMax,
                                          int nVarsMin, int nVarsMax,
                                          int minSize, int maxSize,
                                          int minDegree, int maxDegree) {
            this.nFactorsMin = nFactorsMin;
            this.nFactorsMax = nFactorsMax;
            this.nVarsMin = nVarsMin;
            this.nVarsMax = nVarsMax;
            this.minSize = minSize;
            this.maxSize = maxSize;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
        }

        final RandomGenerator rnd = getRandom();
        final RandomDataGenerator rndd = getRandomData();

        @Override
        public SampleDecomposition<lMultivariatePolynomialZp> next0() {
            lIntegersModulo domain = new lIntegersModulo(
                    getModulusRandom(rndd.nextInt(minModulusBits, maxModulusBits)));

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            int nFactors = rndd.nextInt(nFactorsMin, nFactorsMax);

            lMultivariatePolynomialZp[] factors = new lMultivariatePolynomialZp[nFactors];
            while (nFactors > 0) {
                lMultivariatePolynomialZp sample =
                        randomPolynomial(nVariables,
                                rndd.nextInt(minDegree, maxDegree),
                                rndd.nextInt(minSize, maxSize),
                                domain, DegreeVector.LEX, rnd);

                if (sample.isConstant() || sample.isMonomial())
                    continue;
                factors[--nFactors] = sample;
            }
            return new SampleDecomposition<>(factors);
        }

        private static String range(int from, int to) {
            return "[" + from + ", " + to + "]";
        }

        @Override
        public String toString() {
            return "Sample data: " +
                    " #factors ∈ " + range(nFactorsMin, nFactorsMax) +
                    " #variables ∈ " + range(nVarsMin, nVarsMax) +
                    ", deg ∈ " + range(minDegree, maxDegree) +
                    ", size ∈ " + range(minSize, maxSize);
        }
    }
}