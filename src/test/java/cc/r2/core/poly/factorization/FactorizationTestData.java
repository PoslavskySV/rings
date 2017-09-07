package cc.r2.core.poly.factorization;

import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.IPolynomial;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;

import java.util.Arrays;
import java.util.function.Function;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class FactorizationTestData {

    /** Sample factorization */
    public static class SampleDecomposition<Poly extends IPolynomial<Poly>> {
        public final Poly[] factors;
        public final Poly poly;

        public SampleDecomposition(Poly[] factors) {
            this.factors = factors;
            this.poly = factors[0].createOne().multiply(factors);
            assert Arrays.stream(factors).noneMatch(IPolynomial::isConstant);
            assert Arrays.stream(factors).noneMatch(IPolynomial::isMonomial);
        }

        public void assertFactorization(FactorDecomposition<Poly> factorization) {
            Assert.assertEquals(poly, factorization.toPolynomial());
            Assert.assertTrue(factors.length <= factorization.sumExponents());
        }
    }

    /** Source of the data */
    public static abstract class SampleDecompositionSource<Poly extends IPolynomial<Poly>> {
        public final DescriptiveStatistics
                nFactorsStats = new DescriptiveStatistics(),
                factorsSizeStats = new DescriptiveStatistics(),
                polySizeStats = new DescriptiveStatistics(),
                factorsDegreeStats = new DescriptiveStatistics(),
                polyDegreeStats = new DescriptiveStatistics();

        public final SampleDecomposition<Poly> next() {
            SampleDecomposition<Poly> sample = next0();

            nFactorsStats.addValue(sample.factors.length);
            polySizeStats.addValue(sample.poly.size());
            polyDegreeStats.addValue(sample.poly.degree());
            Arrays.stream(sample.factors).forEach(f -> {
                factorsSizeStats.addValue(f.size());
                factorsDegreeStats.addValue(f.degree());
            });

            return sample;
        }

        public abstract SampleDecomposition<Poly> next0();

        private static String medStr(DescriptiveStatistics stats) {
            return stats.getPercentile(0.5) + " (" + stats.getStandardDeviation() + ")";
        }

        public String statisticsToString() {
            return new StringBuilder()
                    .append("\n========== Samples statistics ==========")
                    .append("\nMed. #factors       : " + medStr(nFactorsStats))
                    .append("\nMed. poly size      : " + medStr(polySizeStats))
                    .append("\nMed. poly degree    : " + medStr(polyDegreeStats))
                    .append("\nMed. factors size   : " + medStr(factorsSizeStats))
                    .append("\nMed. factors degree : " + medStr(factorsDegreeStats))
                    .toString();
        }
    }

    public static class FactorizationAlgorithm<Poly extends IPolynomial<Poly>> {
        public final Function<Poly, FactorDecomposition<Poly>> algorithm;
        public final String name;

        public FactorizationAlgorithm(Function<Poly, FactorDecomposition<Poly>> algorithm, String name) {
            this.algorithm = algorithm;
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        public static <Poly extends IPolynomial<Poly>>
        FactorizationAlgorithm<Poly> named(Function<Poly, FactorDecomposition<Poly>> algorithm, String name) {
            return new FactorizationAlgorithm<>(algorithm, name);
        }
    }

    private static boolean toStringDefined(Object o) {
        String def = o.getClass().getName() + "@" + Integer.toHexString(o.hashCode());
        return !o.toString().equals(def);
    }

    public static <Poly extends IPolynomial<Poly>>
    void testFactorizationAlgorithms(SampleDecompositionSource<Poly> source,
                                     int nIterations,
                                     FactorizationAlgorithm<Poly>... algorithms) {
        System.out.println("Testing factorization algorithms " + Arrays.toString(algorithms));
        if (toStringDefined(source))
            System.out.println("Input source: " + source);

        DescriptiveStatistics[] timings = new DescriptiveStatistics[algorithms.length];
        for (int i = 0; i < timings.length; i++)
            timings[i] = new DescriptiveStatistics();

        int prevProgress = -1, currProgress;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.stream(timings).forEach(DescriptiveStatistics::clear);

            if ((currProgress = (int) (100.0 * n / nIterations)) != prevProgress) {
                prevProgress = currProgress;
                System.out.print(">");
                System.out.flush();
            }
            SampleDecomposition<Poly> sample = source.next();
            for (int i = 0; i < algorithms.length; i++) {
                FactorDecomposition<Poly> decomposition = null;
                try {
                    long start = System.nanoTime();
                    decomposition = algorithms[i].algorithm.apply(sample.poly);
                    timings[i].addValue(System.nanoTime() - start);
                    sample.assertFactorization(decomposition);
                } catch (Throwable throwable) {
                    System.out.println("============ Error ============");
                    System.out.println("Algorithm: " + algorithms[i]);
                    System.out.println("Polynomial: " + sample.poly);
                    System.out.println("Expected factorization: " + Arrays.toString(sample.factors));
                    System.out.println("Actual factorization  : " + decomposition);
                    throw throwable;
                }
            }
        }

        System.out.println(source.statisticsToString());

        System.out.println("\n============ Timings ============");
        for (int i = 0; i < algorithms.length; i++)
            System.out.println(algorithms[i].name + ": " + TimeUnits.statisticsNanotime(timings[i]));
        System.out.println();
    }
}
