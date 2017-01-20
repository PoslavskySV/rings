package cc.r2.core.polynomial;

import cc.r2.core.number.primes.SieveOfAtkin;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static cc.r2.core.polynomial.DistinctDegreeFactorization.DistinctDegreeFactorizationComplete;
import static cc.r2.core.polynomial.DistinctDegreeFactorization.DistinctDegreeFactorizationPlain;
import static cc.r2.core.polynomial.FactorizationTestUtil.assertDistinctDegreeFactorization;
import static cc.r2.core.polynomial.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.polynomial.RandomPolynomials.randomPoly;
import static cc.r2.core.polynomial.SquareFreeFactorization.SquareFreeFactorization;
import static cc.r2.core.polynomial.SquareFreeFactorization.isSquareFree;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by poslavsky on 20/01/2017.
 */
public class DistinctDegreeFactorizationTest {

    //////

    @Test
    public void test34() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, -1, -1, -1, 0, 1, -1, 1, 1);
        Factorization f = DistinctDegreeFactorizationPlain(poly, 3);
        assertDistinctDegreeFactorization(poly, f, 3);
    }

    @Test
    public void test35() throws Exception {
        final class test {
            long overflows = 0;
            long nonSquareFree = 0;
            long iterations = 0;
            final DescriptiveStatistics arithmetics = new DescriptiveStatistics();
            final DescriptiveStatistics timings = new DescriptiveStatistics();
            final DescriptiveStatistics nFactors = new DescriptiveStatistics();
            long[] primes = {2, 3, 5, 7, 11, 13, 101};


            void run(int bound, int minNBase, int maxNBase, int nIterations, boolean ensureSquareFree) {

                RandomGenerator rnd = new Well1024a();
                RandomDataGenerator rndd = new RandomDataGenerator(rnd);

                long start;
                out:
                for (int i = 0; i < nIterations; i++) {
                    if (i == 100)
                        timings.clear();
                    int nFactors = rndd.nextInt(minNBase, maxNBase);
                    for (long modulus : primes) {
                        MutablePolynomial poly;
                        start = System.nanoTime();
                        poly = MutablePolynomial.create(rndd.nextLong(1, 1000)).modulus(modulus);
                        List<MutablePolynomial> factors = new ArrayList<>();
                        int nBases = rndd.nextInt(minNBase, maxNBase);
                        for (int j = 1; j <= nBases; ++j) {
                            MutablePolynomial f;
                            do { f = randomPoly(j, bound, rnd).modulus(modulus);} while (f.degree != j);
                            factors.add(f);
                            poly = poly.multiply(f, modulus);
                        }
                        arithmetics.addValue(System.nanoTime() - start);
                        if (poly.isConstant()) {
                            --i;
                            continue out;
                        }
                        boolean squareFree = false;
                        List<MutablePolynomial> polynomials = new ArrayList<>();
                        if (!(squareFree = isSquareFree(poly, modulus))) {
                            ++nonSquareFree;
                            if (ensureSquareFree) {
                                Factorization factorization = SquareFreeFactorization(poly, modulus);
                                polynomials.addAll(Arrays.asList(factorization.factors));
                            }
                        }
                        if (polynomials.isEmpty()) polynomials.add(poly);

                        for (MutablePolynomial toTest : polynomials) {
                            try {
                                toTest = toTest.multiply(rndd.nextInt(1, 1000), modulus);
                                start = System.nanoTime();
                                Factorization factorization = DistinctDegreeFactorizationPlain(toTest, modulus);
                                timings.addValue(System.nanoTime() - start);
                                this.nFactors.addValue(factorization.factors.length);
                                assertDistinctDegreeFactorization(toTest, factorization, modulus);
                            } catch (ArithmeticException exc) {
                                if (!exc.getMessage().contains("overflow"))
                                    throw exc;
                                ++overflows;
                            }
                        }
                        ++iterations;
                    }
                }
            }

            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();
                sb.append("Total iterations: " + iterations).append("\n");
                sb.append("Overflows: " + overflows).append("\n");
                sb.append("Not square-free: " + nonSquareFree).append("\n");
                sb.append("\nNumber of factors:\n" + nFactors).append("\n");
                sb.append("\nTimings:\n" + timings).append("\n");
                sb.append("\nArithmetics:\n" + arithmetics).append("\n");
                return sb.toString();
            }
        }

        test smallPolys = new test();
        smallPolys.run(10, 1, 5, 10000, false);
        System.out.println(smallPolys);

        test largePolys = new test();
        largePolys.primes = new long[]{29, 31, 89, 101, 107, 139, 223};
        largePolys.run(1000, 8, 10, 20, true);
        System.out.println(largePolys);
    }

    @Test
    public void test36() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        class test {
            final DescriptiveStatistics timings = new DescriptiveStatistics();
            DescriptiveStatistics nFactors = new DescriptiveStatistics();

            void run(int maxDegree, int nIterations, int modBound) {
                SieveOfAtkin sieve = SieveOfAtkin.createSieve(modBound);
                for (int i = 0; i < nIterations; i++) {
                    if (i == 101)
                        timings.clear();
                    int modulus = sieve.randomPrime(rnd);
                    MutablePolynomial poly = randomPoly(rndd.nextInt(0, maxDegree), 10000, rnd).modulus(modulus);
                    long start = System.nanoTime();
                    Factorization factorization = DistinctDegreeFactorizationPlain(poly, modulus);
                    timings.addValue(System.nanoTime() - start);
                    nFactors.addValue(factorization.factors.length);
                    assertDistinctDegreeFactorization(poly, factorization, modulus);
                }
            }

            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();
                sb.append("\nNumber of factors:\n" + nFactors).append("\n");
                sb.append("\nTimings:\n" + timings).append("\n");
                return sb.toString();
            }
        }
        test small = new test();
        small.run(10, 1000, 1000);
        System.out.println(small);
        test big = new test();
        big.run(100, 100, 20);
        System.out.println(big);
    }

    @Test
    public void test35a() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(19, 20, 13, 10, 26, 19, 6, 29, 13, 20, 10, 12, 20, 3, 21, 16, 25, 10, 26, 22, 25, 2, 23, 29, 21, 14, 8, 26, 16, 7, 7, 1);
        Factorization factorization = DistinctDegreeFactorizationPlain(poly, 31);
        assertEquals(5, factorization.factors.length);
        assertDistinctDegreeFactorization(poly, factorization, 31);
    }

    @Test
    public void test37a() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(9, 7, 2, 3, 10, 1, 1);
        Factorization factorization = DistinctDegreeFactorizationPlain(poly, 11);
        assertDistinctDegreeFactorization(poly, factorization, 11);
    }

    @Test
    public void test38() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(172, 85, 84, 151, 122, 53, 107, 117, 82, 152, 133, 151, 178, 1);
        Factorization fct = DistinctDegreeFactorizationComplete(poly, 181);
        assertFactorization(poly, fct, 181);
    }


    static void assertMMATest(int nLines, String file) throws Exception {
        InputStream resource = new GZIPInputStream(FactorizationTestDataTest.class.getClassLoader()
                .getResourceAsStream(file));
        int nEntries = 0;
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (FactorizationTestData fct : FactorizationTestData.allMod(resource)) {
            if (nEntries == 100)
                stats.clear();
            Factorization expected = fct.factorization;
            long start = System.nanoTime();
            System.out.println(Arrays.toString(fct.poly.data));
            Factorization actual = DistinctDegreeFactorizationComplete(fct.poly, fct.modulus);
            stats.addValue(System.nanoTime() - start);
            Assert.assertEquals("Modulus: " + fct.modulus, expected.canonical(), actual.canonical());
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
            MutablePolynomial bigPoly = MutablePolynomial.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1);
            long modulus = 5659;
            long start = System.nanoTime();
            Factorization factorization = DistinctDegreeFactorizationComplete(bigPoly, modulus);
            System.out.println(factorization);
            System.out.println(System.nanoTime() - start);
            assertFactorization(bigPoly, factorization, modulus);
            assertEquals(2, factorization.factors.length);
        }
    }


    //
//    @Ignore
//    @Test
//    public void test40() throws Exception {
//        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
//    }
//
//    @Ignore
//    @Test
//    public void test41() throws Exception {
//        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
//        assertMMATest(1000, "cc/r2/core/polynomial/DistinctDegreeFactorizationLarge.gz");
//        assertMMATest(100, "cc/r2/core/polynomial/DistinctDegreeFactorizationHuge.gz");
//    }
//
//
}