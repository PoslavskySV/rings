package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Ignore;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.multivar.MultivariateResultants.*;
import static cc.redberry.rings.util.TimeUnits.nanosecondsToString;
import static cc.redberry.rings.util.TimeUnits.statisticsNanotime;
import static org.junit.Assert.assertEquals;

/**
 *
 */
public class MultivariateResultantsTest extends AMultivariateTest {
    @Test
    public void test1() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        MultivariatePolynomial<BigInteger> a = coder.parse("(2*x + y + z)^3 + 1");
        MultivariatePolynomial<BigInteger> b = coder.parse("(x - 3*y - z)^13 + 1");
        MultivariatePolynomial<BigInteger> expected = coder.parse("343 - 124509*y^3 + 2470629*y^6 - 40353607*y^9 - 160083*y^2*z + 6353046*y^5*z - 155649627*y^8*z - 68607*y*z^2 + 6806835*y^4*z^2 - 266827932*y^7*z^2 - 9801*z^3 + 3889620*y^3*z^3 - 266827932*y^6*z^3 + 1250235*y^2*z^4 - 171532242*y^5*z^4 + 214326*y*z^5 - 73513818*y^4*z^5 + 15309*z^6 - 21003948*y^3*z^6 - 3857868*y^2*z^7 - 413343*y*z^8 - 19683*z^9");
        assertEquals(expected, PlainResultant(a, b, 0));
    }

    @Test
    @RequiresSingular
    public void testBrown1_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        IntegersZp64 ring = Zp64(SmallPrimes.nextPrime(1 << 20));
        DescriptiveStatistics
                brown = new DescriptiveStatistics(),
                singular = new DescriptiveStatistics();
        int nIterations = its(100, 100);
        int nVars = 4, degree = 7;
        for (int i = 0; i < nIterations; ++i) {
            MultivariatePolynomialZp64
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, MonomialOrder.LEX, rnd),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, MonomialOrder.LEX, rnd);
            int variable = rnd.nextInt(nVars);

            long start = System.nanoTime();
            MultivariatePolynomialZp64 actual = BrownResultant(a, b, variable);
            brown.addValue(System.nanoTime() - start);

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> expected = SingularResultant(a, b, variable);
            singular.addValue(expected.nanotime);

            assertEquals(expected.resultant, actual);
        }

        System.out.println("Brown: " + statisticsNanotime(brown));
        System.out.println("Singular: " + statisticsNanotime(singular));
    }

    @Test
    public void testBrown2() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(5, Zp64(1048583));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4", "x5");

        MultivariatePolynomialZp64 a = coder.parse("637881*x1^3*x2^2*x3^4*x4*x5^2");
        MultivariatePolynomialZp64 b = coder.parse("220444*x1*x2*x3*x5^2+188423*x3^2*x4*x5^3+1021013*x2*x4*x5^5+27846*x1^3*x2^2*x3*x5^4+932116*x1^4*x2^3*x3*x5^3+1029973*x2^3*x4^6*x5^2+405192*x1^6*x2*x3^2*x4^3+982184*x1^2*x3^5*x4^5*x5+867884*x1^4*x2^6*x5^5+607543*x1*x2^4*x3^2*x4^5*x5^3+124055*x2^2*x3^5*x4^6*x5^3+211592*x1^2*x2^3*x3^2*x4^6*x5^3+273881*x1*x2^6*x3^2*x4^4*x5^5+534279*x1^4*x2^3*x3^4*x4^2*x5^6+301771*x1^5*x2^4*x3^5*x4*x5^5+482564*x1^3*x2*x3^6*x4^6*x5^5+422454*x1^6*x2^3*x3^5*x4^2*x5^6");
        MultivariatePolynomialZp64 expected = coder.parse("11341*x1^18*x3^28*x4^8*x5^18+939807*x1^20*x3^31*x4^12*x5^16+83676*x1^22*x3^34*x4^16*x5^14");
        assertEquals(expected, BrownResultant(a, b, 1));
    }

    @Test
    public void testBrown3() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(5, Zp64(1048583));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4", "x5");

        MultivariatePolynomialZp64 a = coder.parse("303182*x2^2*x3^4*x4^3*x5^2+740361*x1^4*x3^3*x4^2*x5^4");
        MultivariatePolynomialZp64 b = coder.parse("446331*x1^4*x3^4*x4^2+972452*x2^5*x3^2*x4^3*x5^5");
        MultivariatePolynomialZp64 expected = coder.parse("822600*x2^8*x3^32*x4^20*x5^8+329308*x2^11*x3^29*x4^20*x5^15+774779*x2^14*x3^26*x4^20*x5^22+72726*x2^17*x3^23*x4^20*x5^29+1007096*x2^20*x3^20*x4^20*x5^36");
        assertEquals(expected, BrownResultant(a, b, 0));
    }

    @Ignore
    @Test
    public void test4() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(5, Zp64(1048583));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z", "p", "q");

        MultivariatePolynomialZp64 a = coder.parse("(2*x + y + z)^3 + p");
        MultivariatePolynomialZp64 b = coder.parse("(x - 3*y - z)^13 + q");
        MultivariatePolynomialZp64 expected = coder.parse("343 - 124509*y^3 + 2470629*y^6 - 40353607*y^9 - 160083*y^2*z + 6353046*y^5*z - 155649627*y^8*z - 68607*y*z^2 + 6806835*y^4*z^2 - 266827932*y^7*z^2 - 9801*z^3 + 3889620*y^3*z^3 - 266827932*y^6*z^3 + 1250235*y^2*z^4 - 171532242*y^5*z^4 + 214326*y*z^5 - 73513818*y^4*z^5 + 15309*z^6 - 21003948*y^3*z^6 - 3857868*y^2*z^7 - 413343*y*z^8 - 19683*z^9");
        for (int i = 0; i < 1000; ++i) {
            long start;

            start = System.nanoTime();
            MultivariatePolynomialZp64 br = BrownResultant(a, b, 0);
            System.out.println("Brown: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            MultivariatePolynomialZp64 zp = ZippelResultant(a, b, 0);
            System.out.println("Zippel: " + nanosecondsToString(System.nanoTime() - start));
            assertEquals(br, zp);
            System.out.println();
        }
    }

    @Test
    public void testSparseInterpolation1() throws IOException, InterruptedException {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(3, Zp64(1048583));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");
        MultivariatePolynomialZp64 a = coder.parse("(2*x + y + z)^3 + 1");
        MultivariatePolynomialZp64 b = coder.parse("(x - 3*y - z)^13 + 1");
        MultivariatePolynomialZp64 expected = SingularResultant(a, b, 0).resultant;

        MultivariatePolynomialZp64
                aEval = a.eliminate(2, 1),
                bEval = b.eliminate(2, 1),
                rEval = expected.dropVariable(0).evaluate(1, 1);


        MultivariateGCD.lSparseInterpolation interpolation =
                createInterpolation(1,
                        a.asUnivariateEliminate(0),
                        b.asUnivariateEliminate(0),
                        rEval, expected.degree(2), getRandom());

        System.out.println(interpolation.evaluate(2).insertVariable(0));
        System.out.println(expected.eliminate(2, 2));
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SingularResult<Term, Poly> SingularResultant(Poly a, Poly b, int variable)
            throws IOException, InterruptedException {
        if (!isSingularAvailable())
            throw new RuntimeException("no Singular");

        Poly factory = a;
        IPolynomialRing<Poly> ring = Rings.PolynomialRing(factory);
        int nVars = factory.nVariables;
        String[] vars = new String[nVars];
        for (int i = 1; i <= nVars; i++)
            vars[i - 1] = "x" + i;

        // prepare tmp file
        File tmpSource = File.createTempFile("singular_" + Integer.toHexString(a.hashCode() + b.hashCode()), null);
        File tmpOutput = File.createTempFile("singular_" + Integer.toHexString(a.hashCode() + b.hashCode()), "output");
        tmpSource.deleteOnExit();
        tmpOutput.deleteOnExit();

        try {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(tmpSource))) {

                writer.write("system(\"--ticks-per-sec\",1000);");
                writer.newLine();

                writer.write(String.format("ring r = %s, (%s), %s;",
                        factory.coefficientRingCardinality() == null ? "0" : factory.coefficientRingCardinality(),
                        Arrays.stream(vars).collect(Collectors.joining(",")),
                        "dp"));
                writer.newLine();

                writer.write(String.format("poly poly1 = %s;\n", a.toString(vars)));
                writer.write(String.format("poly poly2 = %s;\n", b.toString(vars)));

                writer.write("int t = timer;");
                writer.newLine();
                writer.write("poly answ = resultant(poly1, poly2, x" + (variable + 1) + ");\n");
                writer.newLine();
                writer.write("int elapsed = timer-t;");
                writer.newLine();

                writer.write("print(\"OUTPUTSTARTSHERE\");\n");
                writer.newLine();
                writer.write("print(answ);");
                writer.newLine();

                writer.write("print(\"TIMESEPARATOR\");");
                writer.newLine();
                writer.write("print(elapsed);");
                writer.newLine();

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

            singularOut = singularOut.split("OUTPUTSTARTSHERE")[1];
            String[] answ = singularOut.split("TIMESEPARATOR");

            return new SingularResult<>(
                    ring.parse(answ[0].replace("\n", "").trim(), vars),
                    1000 * 1000 * Long.valueOf(answ[1].replace("\n", "")));
        } finally {
            tmpSource.delete();
            tmpOutput.delete();
        }
    }

    static final class SingularResult<
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        final Poly resultant;
        final long nanotime;

        public SingularResult(Poly resultant, long nanotime) {
            this.resultant = resultant;
            this.nanotime = nanotime;
        }

        @Override
        public String toString() {
            return nanosecondsToString(nanotime) + " " + resultant;
        }
    }
}