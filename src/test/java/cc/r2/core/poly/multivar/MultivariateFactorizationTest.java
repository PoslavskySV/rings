package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.factorization.FactorizationTestData.FactorizationAlgorithm;
import cc.r2.core.poly.factorization.FactorizationTestData.SampleDecomposition;
import cc.r2.core.poly.factorization.FactorizationTestData.SampleDecompositionSource;
import cc.r2.core.poly.univar.DivisionWithRemainder;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.multivar.MultivariateFactorization.GCDFreeBasis;
import static cc.r2.core.poly.multivar.MultivariateFactorization.bivariateDenseFactorSquareFree;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;

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

        System.out.println(d1.toPolynomial().equals(decomps[0].toPolynomial()));
        System.out.println(Arrays.toString(DivisionWithRemainder.divideAndRemainder(d1.toPolynomial(), decomps[0].toPolynomial(), true)));
        for (FactorDecomposition<lUnivariatePolynomialZp> decomp : decomps)
            System.out.println(decomp.size() + " => " + decomp);
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
                long start = System.nanoTime();
                lDecomposition = lAlgorithm.algorithm.apply(sample.poly);
                lTiming.addValue(System.nanoTime() - start);
                sample.assertFactorization(lDecomposition);


                SampleDecomposition<MultivariatePolynomial<BigInteger>>
                        bSample = toBigPoly(sample);
                start = System.nanoTime();
                bDecomposition = bAlgorithm.algorithm.apply(bSample.poly);
                bTiming.addValue(System.nanoTime() - start);
                bSample.assertFactorization(bDecomposition);

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

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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