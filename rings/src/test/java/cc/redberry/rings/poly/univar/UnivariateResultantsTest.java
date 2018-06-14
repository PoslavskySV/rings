package cc.redberry.rings.poly.univar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.util.RandomDataGenerator;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.function.BiFunction;

import static cc.redberry.rings.Rings.Z;
import static cc.redberry.rings.Rings.Zp;
import static cc.redberry.rings.poly.univar.RandomUnivariatePolynomials.randomPoly;
import static cc.redberry.rings.poly.univar.UnivariatePolynomial.asOverZ64;
import static cc.redberry.rings.poly.univar.UnivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.univar.UnivariateResultants.*;
import static org.junit.Assert.assertEquals;

/**
 *
 */
public class UnivariateResultantsTest extends APolynomialTest {
    @Test
    public void test1() throws Exception {
        UnivariatePolynomial<BigInteger>
                a = UnivariatePolynomialZ64.create(1, 2, 13).toBigPoly(),
                b = UnivariatePolynomialZ64.create(3, 2, 17).toBigPoly();

        BigInteger expected = SingularResultant(0, a, b);
        assertEquals(expected, PseudoPRS(a, b).resultant());
        assertEquals(expected, PrimitivePRS(a, b).resultant());
        assertEquals(expected, ReducedPRS(a, b).resultant());
        assertEquals(expected, SubresultantPRS(a, b).resultant());
    }

    @Test
    public void test2() throws Exception {
        UnivariatePolynomial<BigInteger>
                a = UnivariatePolynomialZ64.create(-45, 18, 72, -27, -27, 0, 9).toBigPoly(),
                b = UnivariatePolynomialZ64.create(21, -9, -4, 0, 3).toBigPoly();

        BigInteger expected = SingularResultant(0, a, b);
        assertEquals(expected, PseudoPRS(a, b).resultant());
        assertEquals(expected, PrimitivePRS(a, b).resultant());
        assertEquals(expected, ReducedPRS(a, b).resultant());
        assertEquals(expected, SubresultantPRS(a, b).resultant());

        int mod = 17;
        a = a.setRing(Zp(mod));
        b = b.setRing(Zp(mod));
        expected = SingularResultant(mod, a, b);
        assertEquals(expected, PseudoPRS(a, b).resultant());
        assertEquals(expected, PrimitivePRS(a, b).resultant());
        assertEquals(expected, ReducedPRS(a, b).resultant());
        assertEquals(expected, SubresultantPRS(a, b).resultant());
        assertEquals(expected.longValueExact(), ClassicalPRS(asOverZp64(a), asOverZp64(b)).resultant());
    }

    @Test
    @RequiresSingular
    public void test3_random() throws Exception {
        testAlgorithmsRandom(10, getRandom(), 15, 0, true, UnivariateResultants::PseudoPRS);
    }

    @Test
    @RequiresSingular
    public void test4_random() throws Exception {
        for (int characteristic : Arrays.asList(0, 17)) {
            System.out.println("Characteristic " + characteristic + ":");
            testAlgorithmsRandom(
                    its(100, 100), getRandom(),
                    100, characteristic, true,
                    UnivariateResultants::PrimitivePRS,
                    UnivariateResultants::ReducedPRS,
                    UnivariateResultants::SubresultantPRS);
        }
    }

    static void testAlgorithmsRandom(int nIterations,
                                     RandomGenerator rnd,
                                     int maxDegree,
                                     int characteristic,
                                     boolean testZero,
                                     BiFunction<UnivariatePolynomial<BigInteger>,
                                             UnivariatePolynomial<BigInteger>,
                                             PolynomialRemainderSequence<BigInteger>>... algorithms) throws Exception {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics[] stats = new DescriptiveStatistics[algorithms.length];
        for (int i = 0; i < stats.length; ++i)
            stats[i] = new DescriptiveStatistics();

        for (int i = 0; i < nIterations; ++i) {
            UnivariatePolynomial<BigInteger>
                    a = randomPoly(rndd.nextInt(0, maxDegree), rnd).toBigPoly(),
                    b = randomPoly(rndd.nextInt(0, maxDegree), rnd).toBigPoly();

            if (characteristic == 0) {
                a = a.setRing(Zp(10)).setRing(Z);
                b = b.setRing(Zp(10)).setRing(Z);
            } else {
                a = a.setRing(Zp(characteristic));
                b = b.setRing(Zp(characteristic));
            }

            if (a.isZero() || b.isZero()) {
                --i;
                continue;
            }

            BigInteger expected = a.ring.valueOf(SingularResultant(characteristic, a, b));
            int c = 0;
            for (BiFunction<UnivariatePolynomial<BigInteger>, UnivariatePolynomial<BigInteger>, PolynomialRemainderSequence<BigInteger>> algorithm : algorithms) {
                long start = System.nanoTime();
                try {
                    assertEquals(expected, algorithm.apply(a, b).resultant());
                } catch (Throwable t) {
                    System.out.println(asOverZ64(a).toStringForCopy());
                    System.out.println(asOverZ64(b).toStringForCopy());
                    throw t;
                }
                stats[c].addValue(System.nanoTime() - start);
                ++c;
            }

            if (characteristic != 0)
                assertEquals(expected.longValueExact(), ClassicalPRS(asOverZp64(a), asOverZp64(b)).resultant());

            if (testZero) {
                UnivariatePolynomial<BigInteger> g;
                do {
                    g = randomPoly(rndd.nextInt(1, maxDegree), rnd).toBigPoly();
                    if (characteristic == 0)
                        g = g.setRing(Zp(10)).setRing(Z);
                    else
                        g = g.setRing(Zp(characteristic));
                } while (g.degree == 0);
                a.multiply(g);
                b.multiply(g);
                expected = Z.getZero();
                c = 0;
                for (BiFunction<UnivariatePolynomial<BigInteger>, UnivariatePolynomial<BigInteger>, PolynomialRemainderSequence<BigInteger>> algorithm : algorithms) {
                    long start = System.nanoTime();
                    assertEquals(expected, algorithm.apply(a, b).resultant());
                    stats[c].addValue(System.nanoTime() - start);
                }
            }
        }

        for (int i = 0; i < stats.length; ++i)
            System.out.println("Algorithm #" + i + ": " + TimeUnits.statisticsNanotime(stats[i]));
    }

    @Test
    public void test5() throws Exception {
        UnivariatePolynomial<BigInteger>
                a = UnivariatePolynomialZ64.create(9, 0, 0, 9, 1, 6).toBigPoly(),
                b = UnivariatePolynomialZ64.create(5, 3, 0, 8, 8, 3, 0, 1, 5, 7, 2, 3, 2, 2, 5, 9).toBigPoly();
        assertEquals(Z.parse("-34729627388605877211"), PseudoPRS(a, b).resultant());
    }

    static BigInteger SingularResultant(int characteristic,
                                        UnivariatePolynomial<BigInteger> a,
                                        UnivariatePolynomial<BigInteger> b)
            throws IOException, InterruptedException {
        if (!isSingularAvailable())
            throw new RuntimeException("no Singular");

        // prepare tmp file
        File tmpSource = File.createTempFile("singular_" + Integer.toHexString(a.hashCode() + b.hashCode()), null);
        File tmpOutput = File.createTempFile("singular_" + Integer.toHexString(a.hashCode() + b.hashCode()), "output");
        tmpSource.deleteOnExit();
        tmpOutput.deleteOnExit();

        try {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(tmpSource))) {

                writer.write(String.format("ring r = %s,(x),dp; \n", characteristic));
                writer.write(String.format("poly a = %s;\n", a));
                writer.write(String.format("poly b = %s;\n", b));
                writer.write("poly answ = resultant(a, b, x);\n");
                writer.write("print(\"OUTPUTSTARTSHERE\");\n");
                writer.write("print(answ);");
                writer.write("exit;");
            }

            //Files.readAllLines(tmpSource.toPath()).forEach(System.out::println);

            Process process = new ProcessBuilder(getSingularPath(), "-q")
                    .redirectOutput(tmpOutput)
                    .redirectErrorStream(true)
                    .start();

            process.getOutputStream().write(String.format("< \"%s\";", tmpSource.getAbsolutePath()).getBytes());
            process.getOutputStream().flush();
            process.getOutputStream().close();

            process.waitFor();
            String singularOut = Files.readAllLines(tmpOutput.toPath()).stream().reduce((l, r) -> l + r).orElse("");
            if (!singularOut.contains("OUTPUTSTARTSHERE"))
                //<- there was a error in Singular
                return null;

            return new BigInteger(singularOut.split("OUTPUTSTARTSHERE")[1]);
        } finally {
            tmpSource.delete();
            tmpOutput.delete();
        }
    }
}