package cc.r2.core.poly2;

import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly2.FactorizationTestUtil.PolynomialSource;
import cc.r2.core.poly2.FactorizationTestUtil.RandomSource;
import cc.r2.core.poly2.FactorizationTestUtil.ShoupSource;
import cc.r2.core.polynomial.FactorizationTestDataTest;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.zip.GZIPInputStream;

import static cc.r2.core.poly2.DistinctDegreeFactorization.*;
import static cc.r2.core.poly2.FactorizationTestUtil.assertDistinctDegreeFactorization;
import static cc.r2.core.poly2.FactorizationTestUtil.assertFactorization;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 19/02/2017.
 */
public class DistinctDegreeFactorizationTest {
    static final int bigModulus = 5659;
    static final MutablePolynomialMod bigPoly = MutablePolynomialZ.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1).modulus(bigModulus);

    public void test1() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialZ.create(0, -1, -1, -1, 0, 1, -1, 1, 1).modulus(3);
        FactorDecomposition<MutablePolynomialMod> f = DistinctDegreeFactorizationPlain(poly);
        assertDistinctDegreeFactorization(poly, f);
    }

    @Test
    public void test2() throws Exception {
        long modulus = 11;
        MutablePolynomialMod poly = MutablePolynomialZ
                .create(32, 1, 4)
                .modulus(modulus)
                .multiply(MutablePolynomialZ
                        .create(3, 3, 1)
                        .modulus(modulus));
        FactorDecomposition<MutablePolynomialMod> plain = DistinctDegreeFactorizationPlain(poly);
        FactorDecomposition<MutablePolynomialMod> shoup = DistinctDegreeFactorizationShoup(poly);
        assertDistinctDegreeFactorization(poly, plain);
        assertDistinctDegreeFactorization(poly, shoup);
        assertEquals(plain, shoup);
    }

    @Test
    public void test3() throws Exception {
        long modulus = 11;
        MutablePolynomialMod poly = MutablePolynomialZ
                .create(1, 1)
                .modulus(modulus)
                .multiply(MutablePolynomialZ
                        .create(4, 9, 1)
                        .modulus(modulus));
        System.out.println(poly);
        FactorDecomposition<MutablePolynomialMod> plain = DistinctDegreeFactorizationPlain(poly);
        FactorDecomposition<MutablePolynomialMod> shoup = DistinctDegreeFactorizationShoup(poly);
        System.out.println(plain);
        System.out.println(shoup);
        assertDistinctDegreeFactorization(poly, plain);
        assertDistinctDegreeFactorization(poly, shoup);
        assertEquals(plain, shoup);

    }

    @Ignore
    @Test
    public void test2a() throws Exception {
        for (int i = 0; i < 1000; i++) {
            System.out.println("----");
            long start;
            start = System.nanoTime();
            FactorDecomposition<MutablePolynomialMod> f1 = DistinctDegreeFactorizationPlain(bigPoly);
            System.out.println(System.nanoTime() - start);

            start = System.nanoTime();
            FactorDecomposition<MutablePolynomialMod> f2 = DistinctDegreeFactorizationShoup(bigPoly);
            System.out.println(System.nanoTime() - start);

            assertEquals(f1, f2);
        }
    }

    @Test
    public void test2_benchmark_random_small_polys() throws Exception {
        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(10), SmallPrimes.nextPrime(100), SmallPrimes.nextPrime(200)})
//                .setPrimes(new long[]{BigPrimes.nextPrime(1L << 4)})
                .setSource(new RandomSource(new Well1024a(), 5, 15, true))
                .setnIterations(1000)
                .setPrintProgress(false)
                .setName("Small random polynomials")
                .run();
    }

    @Test
    public void test2_benchmark_random_small_polys_large_modulus() throws Exception {
        new Benchmark()
                .setSource(new RandomSource(new Well1024a(), 150, 200, true))
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(10)})
                .setnIterations(1000)
                .setPrintProgress(false)
                .setName("Small random polynomials")
                .run();
    }

    @Test
    public void test2_benchmark_random_medium1_polys() throws Exception {
        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(10), SmallPrimes.nextPrime(100), SmallPrimes.nextPrime(200)})
                .setSource(new RandomSource(new Well1024a(), 15, 30, true))
                .setnIterations(10000)
                .setPrintProgress(false)
                .setName("Medium 15-30 random polynomials")
                .run();
    }

    @Test
    public void test2_benchmark_random_medium2_polys() throws Exception {
        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(10), SmallPrimes.nextPrime(100), SmallPrimes.nextPrime(200)})
                .setSource(new RandomSource(new Well1024a(), 30, 60, true))
                .setnIterations(10000)
                .setPrintProgress(false)
                .setName("Medium 30-60 random polynomials")
                .run();
    }

    @Test
    public void test2_benchmark_random_medium3_polys() throws Exception {
        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(10), SmallPrimes.nextPrime(100), SmallPrimes.nextPrime(200)})
                .setSource(new RandomSource(new Well1024a(), 60, 120, true))
                .setnIterations(10000)
                .setPrintProgress(false)
                .setName("Medium 60-120 random polynomials")
                .run();
    }

    @Test
    public void test2_shoup_source() throws Exception {
        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(1234)})
                .setSource(new ShoupSource(new Well1024a(), 10, 20))
                .setnIterations(100)
                .setPrintProgress(true)
                .setName("Shoup's benchmark")
                .run();

    }

    @Test
    public void test2sa() throws Exception {

        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(1234)})
                .setSource(new RandomSource(new Well1024a(), 50, 100, true))
                .setnIterations(100)
                .setPrintProgress(true)
                .run();

        new Benchmark()
                .setAlgorithms(EnumSet.allOf(DDFAlgorithm.class))
                .setPrimes(new long[]{SmallPrimes.nextPrime(16234)})
                .setSource(new FactorizationTestUtil.FactorableSource(new Well1024a(), 22, 25, true))
                .setnIterations(100)
                .setPrintProgress(true)
                .run();
    }


    enum DDFAlgorithm {
        Plain {
            FactorizationTestUtil.WithTiming<FactorDecomposition<MutablePolynomialMod>> factor(MutablePolynomialMod poly) {
                long start = System.nanoTime();
                FactorDecomposition<MutablePolynomialMod> r = DistinctDegreeFactorizationPlain(poly);
                long timing = System.nanoTime() - start;
                return new FactorizationTestUtil.WithTiming<>(r, timing);
            }
        },
        XPowers {
            FactorizationTestUtil.WithTiming<FactorDecomposition<MutablePolynomialMod>> factor(MutablePolynomialMod poly) {
                long start = System.nanoTime();
                FactorDecomposition<MutablePolynomialMod> r = DistinctDegreeFactorizationPrecomputedExponents(poly);
                long timing = System.nanoTime() - start;
                return new FactorizationTestUtil.WithTiming<>(r, timing);
            }
        },
        Shoup {
            FactorizationTestUtil.WithTiming<FactorDecomposition<MutablePolynomialMod>> factor(MutablePolynomialMod poly) {
                long start = System.nanoTime();
                FactorDecomposition<MutablePolynomialMod> r = DistinctDegreeFactorizationShoup(poly);
                long timing = System.nanoTime() - start;
                return new FactorizationTestUtil.WithTiming<>(r, timing);
            }
        };
//        ,
//        Shoup2 {
//            WithTiming<FactorDecomposition<MutablePolynomialMod>> factor(MutablePolynomialMod poly, long modulus) {
//                long start = System.nanoTime();
//                FactorDecomposition<MutablePolynomialMod> r = ShoupDDFFactorization2(poly, modulus);
//                long timing = System.nanoTime() - start;
//                return new WithTiming<>(r, timing);
//            }
//        };

        abstract FactorizationTestUtil.WithTiming<FactorDecomposition<MutablePolynomialMod>> factor(MutablePolynomialMod poly);
    }

    final class Benchmark {
        Set<DDFAlgorithm> algorithms;
        long[] primes;
        PolynomialSource source;
        int nIterations;
        boolean printProgress = false;
        String name;

        public Benchmark setAlgorithms(Set<DDFAlgorithm> algorithms) {
            this.algorithms = algorithms;
            return this;
        }

        public Benchmark setPrimes(long[] primes) {
            this.primes = primes;
            return this;
        }

        public Benchmark setSource(PolynomialSource source) {
            this.source = source;
            return this;
        }

        public Benchmark setnIterations(int nIterations) {
            this.nIterations = nIterations;
            return this;
        }

        public Benchmark setPrintProgress(boolean printProgress) {
            this.printProgress = printProgress;
            return this;
        }

        public Benchmark setName(String name) {
            this.name = name;
            return this;
        }

        void run() throws IOException {
            EnumMap<DDFAlgorithm, DescriptiveStatistics> performance = new EnumMap<>(DDFAlgorithm.class);
            algorithms.forEach(al -> performance.put(al, new DescriptiveStatistics()));
            DescriptiveStatistics degrees = new DescriptiveStatistics(), factors = new DescriptiveStatistics();

            long overflows = 0, iterations = 0;
            int prevPercent = 0;
            out:
            for (int i = 0; i < nIterations; i++) {
                if (printProgress) {
                    int currPercent = (100 * i) / nIterations;
                    if (currPercent != prevPercent) {
                        prevPercent = currPercent;
                        System.out.println("->" + currPercent + "%");
                    }
                }
                if (i == nIterations / 10) {
                    performance.values().forEach(DescriptiveStatistics::clear);
                    degrees.clear();
                    factors.clear();
                }

                for (long modulus : primes) {
                    MutablePolynomialMod poly = source.take(modulus);
                    degrees.addValue(poly.degree);
                    try {
                        Map<DDFAlgorithm, FactorizationTestUtil.WithTiming<FactorDecomposition<MutablePolynomialMod>>> results = algorithms.stream().collect(Collectors.toMap(al -> al, al -> al.factor(poly)));
                        DDFAlgorithm baseAlgorithm = algorithms.stream().findFirst().get();
                        FactorDecomposition<MutablePolynomialMod> baseDecomposition = results.get(baseAlgorithm).val;
                        factors.addValue(baseDecomposition.factors.size());
                        results.forEach((alg, res) -> {
                            performance.get(alg).addValue(res.nanoSeconds);

                            FactorDecomposition<MutablePolynomialMod> c = res.val;
                            try {
                                assertDistinctDegreeFactorization(poly, c);
                                assertEquals(baseDecomposition.factors.size(), c.factors.size());
                                assertEquals(baseDecomposition, c);
                            } catch (AssertionError err) {
                                System.out.println("Modulus: " + modulus);
                                System.out.println("Poly: " + poly);
                                System.out.println(baseAlgorithm + ":  " + baseDecomposition.factors.size());
                                System.out.println(alg + ":  " + res.val.factors.size());

                                System.out.println(baseAlgorithm + " : " + baseDecomposition);
                                System.out.println(alg + " : " + res.val);
                                throw err;
                            }
                        });
                    } catch (ArithmeticException exc) {
                        if (!exc.getMessage().contains("overflow"))
                            throw exc;
                        ++overflows;
                    }
                    ++iterations;
                }
            }
            if (printProgress) System.out.println();

            List<String> headers = algorithms.stream().map(DDFAlgorithm::toString).collect(Collectors.toList());
            List<String> text = performance.values().stream().map(DescriptiveStatistics::toString).collect(Collectors.toList());
            headers.add(0, "Factorization");
            text.add(0, factors.toString());
            headers.add(0, "Degrees");
            text.add(0, degrees.toString());

            StringBuilder sb = new StringBuilder();
            String sep = printColumns(headers, text, sb);
            if (name != null) {
                System.out.println(constString('=', sep.length()));
                System.out.println(constString(' ', sep.length() / 2 - name.length() / 2) + name);
            }
            System.out.println(sep);
            System.out.println("Total iterations: " + iterations);
            System.out.println("Modulus: " + Arrays.toString(primes));
            System.out.println("Overflows: " + overflows);
            System.out.println(sb.toString());
        }
    }

    static String printColumns(List<String> headers, List<String> text, StringBuilder sb) {
        List<List<String>> lines = text.stream().map(x -> Arrays.asList(x.split("\n"))).collect(Collectors.toList());
        String fmt = "%-30.30s";
        String fmtString = "";
        for (int i = 0; i < headers.size(); i++)
            fmtString = fmtString + fmt + "     ";
        fmtString += "%n";


        String header = String.format(fmtString, headers.toArray());
        int sepLength = 0;
        for (int i = 0; i < header.length(); i++)
            ++sepLength;
        String separator = constString('-', sepLength);
        sb.append(separator).append("\n");
        sb.append(header);
        sb.append(separator).append("\n");
        for (int i = 0; i < lines.get(0).size(); i++) {
            String[] d = new String[headers.size()];
            for (int j = 0; j < d.length; j++)
                d[j] = lines.get(j).get(i);

            sb.append(String.format(fmtString, d));
        }
        sb.append(separator).append("\n");
        return separator;
    }

    static String constString(char ch, int len) {
        char[] data = new char[len];
        Arrays.fill(data, ch);
        return new String(data);
    }


    @Test
    public void test35a() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialZ
                .create(19, 20, 13, 10, 26, 19, 6, 29, 13, 20, 10, 12, 20, 3, 21, 16, 25, 10, 26, 22, 25, 2, 23, 29, 21, 14, 8, 26, 16, 7, 7, 1)
                .modulus(31);
        FactorDecomposition<MutablePolynomialMod> factorization = DistinctDegreeFactorizationPlain(poly);
        assertEquals(5, factorization.factors.size());
        assertDistinctDegreeFactorization(poly, factorization);
    }

    @Test
    public void test37a() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialZ
                .create(9, 7, 2, 3, 10, 1, 1)
                .modulus(11);
        FactorDecomposition<MutablePolynomialMod> factorization = DistinctDegreeFactorizationPlain(poly);
        assertDistinctDegreeFactorization(poly, factorization);
    }

    @Test
    public void test38() throws Exception {
        MutablePolynomialMod poly = MutablePolynomialZ
                .create(172, 85, 84, 151, 122, 53, 107, 117, 82, 152, 133, 151, 178, 1)
                .modulus(181);
        FactorDecomposition<MutablePolynomialMod> fct = DistinctDegreeFactorizationComplete(poly);
        assertFactorization(poly, fct);
    }


    static void assertMMATest(int nLines, String file) throws Exception {
        InputStream resource = new GZIPInputStream(FactorizationTestDataTest.class.getClassLoader()
                .getResourceAsStream(file));
        int nEntries = 0;
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (FactorizationTestData fct : FactorizationTestData.allMod(resource)) {
            if (nEntries == 100)
                stats.clear();
            FactorDecomposition<MutablePolynomialMod> expected = fct.factorization;
            long start = System.nanoTime();
//            System.out.println(Arrays.toString(fct.poly.data));
            FactorDecomposition<MutablePolynomialMod> actual = DistinctDegreeFactorizationComplete(fct.poly);
            stats.addValue(System.nanoTime() - start);
            Assert.assertEquals("Modulus: " + fct.modulus, expected.canonicalForm(), actual.canonicalForm());
            ++nEntries;
        }
        assertTrue(nEntries > 0);
        if (nLines != -1)
            assertEquals(nLines, nEntries);
        System.out.println(stats);
    }

    @Ignore
    @Test
    public void testDDFLarge() throws Exception {
        for (int i = 0; i < 100; i++) {
            long modulus = 5659;
            MutablePolynomialMod bigPoly = MutablePolynomialZ.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1).modulus(modulus);
            System.out.println(bigPoly);
            long start = System.nanoTime();
            FactorDecomposition<MutablePolynomialMod> factorization = DistinctDegreeFactorizationComplete(bigPoly);
            System.out.println(factorization);
            System.out.println(System.nanoTime() - start);
            assertFactorization(bigPoly, factorization);
            assertEquals(2, factorization.factors.size());
        }
    }


    @Ignore
    @Test
    public void test40() throws Exception {
        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
    }

    @Ignore
    @Test
    public void test41() throws Exception {
        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
        assertMMATest(1000, "cc/r2/core/polynomial/DistinctDegreeFactorizationLarge.gz");
        assertMMATest(100, "cc/r2/core/polynomial/DistinctDegreeFactorizationHuge.gz");
    }
}