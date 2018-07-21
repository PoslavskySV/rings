package cc.redberry.rings.poly.univar;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.AlgebraicNumberField;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
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
import java.util.function.BiFunction;
import java.util.function.Supplier;

import static cc.redberry.rings.Rings.*;
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
        testAlgorithmsRandom(10, getRandom(), 15, 0, true, (p, q) -> PseudoPRS(p, q).resultant());
    }

    @Test
    @RequiresSingular
    public void test4_random() throws Exception {
        System.out.println("Characteristic 0:");
        testAlgorithmsRandom(
                its(100, 100), getRandom(),
                100, 0, true,
                UnivariateResultants::ModularResultant,
                (p, q) -> PrimitivePRS(p, q).resultant(),
                (p, q) -> ReducedPRS(p, q).resultant(),
                (p, q) -> SubresultantPRS(p, q).resultant());

        System.out.println("Characteristic 17:");
        testAlgorithmsRandom(
                its(100, 100), getRandom(),
                100, 17, true,
                (p, q) -> ClassicalPRS(p, q).resultant(),
                (p, q) -> PrimitivePRS(p, q).resultant(),
                (p, q) -> ReducedPRS(p, q).resultant(),
                (p, q) -> SubresultantPRS(p, q).resultant());
    }

    static void testAlgorithmsRandom(int nIterations,
                                     RandomGenerator rnd,
                                     int maxDegree,
                                     int characteristic,
                                     boolean testZero,
                                     BiFunction<UnivariatePolynomial<BigInteger>,
                                             UnivariatePolynomial<BigInteger>,
                                             BigInteger>... algorithms) throws Exception {
        testAlgorithmsRandom(nIterations, rnd, maxDegree, characteristic, testZero, BigInteger.valueOf(5), algorithms);
    }

    static void testAlgorithmsRandom(int nIterations,
                                     RandomGenerator rnd,
                                     int maxDegree,
                                     int characteristic,
                                     boolean testZero,
                                     BigInteger cfBound,
                                     BiFunction<UnivariatePolynomial<BigInteger>,
                                             UnivariatePolynomial<BigInteger>,
                                             BigInteger>... algorithms) throws Exception {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics[] stats = new DescriptiveStatistics[algorithms.length];
        for (int i = 0; i < stats.length; ++i)
            stats[i] = new DescriptiveStatistics();

        for (int i = 0; i < nIterations; ++i) {
            UnivariatePolynomial<BigInteger>
                    a = randomPoly(rndd.nextInt(0, maxDegree), rnd).toBigPoly(),
                    b = randomPoly(rndd.nextInt(0, maxDegree), rnd).toBigPoly();

            if (characteristic == 0) {
                a = a.setRing(Zp(cfBound)).setRing(Z);
                b = b.setRing(Zp(cfBound)).setRing(Z);
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
            for (BiFunction<UnivariatePolynomial<BigInteger>, UnivariatePolynomial<BigInteger>, BigInteger> algorithm : algorithms) {
                long start = System.nanoTime();
                try {
                    assertEquals(expected, algorithm.apply(a, b));
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
                        g = g.setRing(Zp(cfBound)).setRing(Z);
                    else
                        g = g.setRing(Zp(characteristic));
                } while (g.degree == 0);
                a.multiply(g);
                b.multiply(g);
                expected = Z.getZero();
                c = 0;
                for (BiFunction<UnivariatePolynomial<BigInteger>, UnivariatePolynomial<BigInteger>, BigInteger> algorithm : algorithms) {
                    long start = System.nanoTime();
                    assertEquals(expected, algorithm.apply(a, b));
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

    @Test
    public void test6() throws Exception {
        UnivariatePolynomial<BigInteger>
                a = UnivariatePolynomialZ64.create(9, 0, 123, 0, 9, -1231231431, 6, 5).toBigPoly().square().square(),
                b = UnivariatePolynomialZ64.create(5, 3, 0, 8, 12123434, 8, 3, 0, 2341, 5, 7, 2, 3, 2, 2, 5, 1234329).toBigPoly().square();

        a.square();
        b.square();
        for (int i = 0; i < 2; ++i) {
            long start;

            start = System.nanoTime();
            BigInteger direct = SubresultantPRS(a, b).resultant();
            System.out.println("Direct: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            BigInteger modular = ModularResultant(a, b);
            System.out.println("Modular: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            assertEquals(direct, modular);
            System.out.println();
        }
    }

    @Test
    public void test7() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>>
                a = UnivariatePolynomial.parse("1/22 + 13*x^6 - x^10/91", Q, "x"),
                b = UnivariatePolynomial.parse("22 - 13*x^3/24 - 23*x^11/91", Q, "x");
        assertEquals(Q.parse("425247151599448037618564488318977897156050758819345410383614772923/512834526595668226887325363840230565837440811008"), Resultant(a, b));
    }

    @Test
    public void test8() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        MultivariatePolynomial<BigInteger>
                a = coder.parse("(2*x + y + z)^11 - 1"),
                b = coder.parse("(x - 3*y - z)^3 + 1");

        UnivariatePolynomial<MultivariatePolynomial<BigInteger>>
                ua = a.asUnivariateEliminate(0),
                ub = b.asUnivariateEliminate(0);
        for (int i = 0; i < its(2, 2); ++i) {
            long start;

            start = System.nanoTime();
            MultivariatePolynomial<BigInteger> sr = SubresultantPRS(ua, ub).resultant();
            System.out.println("Subresultant PRS: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            MultivariatePolynomial<BigInteger> rr = ReducedPRS(ua, ub).resultant();
            System.out.println("Reduced      PRS: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();

            assertEquals(sr, rr);
        }
    }

    @Test
    public void test9() throws Exception {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> cfRing = new AlgebraicNumberField<>(Coder.mkUnivariateCoder(UnivariateRing(Q), "x").parse("-2 + x^2"));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(cfRing, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> ring = UnivariateRing(cfRing);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(ring, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                a = coder.parse("s + x^2*s/3"),
                b = coder.parse("s + x*s/3 + x^3");
        assertEquals(SubresultantPRS(a, b).resultant(), ModularResultantInNumberField(a, b));

        a = a.clone().multiply(b.clone().increment()).increment();
        b = b.clone().multiply(a.clone().increment()).decrement();
        assertEquals(SubresultantPRS(a, b).resultant(), ModularResultantInNumberField(a, b));
    }

    @Test
    public void test10() throws Exception {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> cfRing
                = new AlgebraicNumberField<>(Coder.mkUnivariateCoder(UnivariateRing(Q), "x").parse("-2 + 17*x^2"));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(cfRing, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> ring = UnivariateRing(cfRing);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(ring, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                a = coder.parse("s + x^2*s/3"),
                b = coder.parse("s + x*s/3 + x^3");
        assertEquals(SubresultantPRS(a, b).resultant(), ModularResultantInNumberField(a, b));

        a = a.clone().multiply(b.clone().increment()).increment();
        b = b.clone().multiply(a.clone().increment()).decrement();
        assertEquals(SubresultantPRS(a, b).resultant(), ModularResultantInNumberField(a, b));
    }

    @Test
    public void test11_random_in_number_field() {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                subresultant = new DescriptiveStatistics(),
                modular = new DescriptiveStatistics();

        for (int i = 0; i < 10; ++i) {
            UnivariatePolynomial<Rational<BigInteger>> minimalPoly = IrreduciblePolynomials
                    .randomIrreduciblePolynomialOverZ(rndd.nextInt(2, 5), rnd)
                    .mapCoefficients(Q, Q::mkNumerator);

            minimalPoly.setLC(Q.mkNumerator(rndd.nextInt(1, 10)));
            if (!IrreduciblePolynomials.irreducibleQ(minimalPoly)) {
                --i;
                continue;
            }

            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = new AlgebraicNumberField<>(minimalPoly);
            System.out.println("Field        : " + field);
            Supplier<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> rndPoly = () ->
                    RandomUnivariatePolynomials.randomPoly(
                            rndd.nextInt(2, 8),
                            field,
                            __ -> RandomUnivariatePolynomials.randomPoly(minimalPoly.degree, Q, ___ -> Q.mk(rnd.nextInt(10), 1 + rnd.nextInt(10)), rnd),
                            rnd);

            for (int j = 0; j < 5; ++j) {
                UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                        a = rndPoly.get(),
                        b = rndPoly.get();

                long start, elapsed;

                start = System.nanoTime();
                UnivariatePolynomial<Rational<BigInteger>> expected = SubresultantPRS(a, b).resultant();
                elapsed = System.nanoTime() - start;
                subresultant.addValue(elapsed);
                System.out.println("Subresultant : " + TimeUnits.nanosecondsToString(elapsed));

                start = System.nanoTime();
                UnivariatePolynomial<Rational<BigInteger>> actual = ModularResultantInNumberField(a, b);
                elapsed = System.nanoTime() - start;
                modular.addValue(elapsed);
                assertEquals(expected, actual);
                System.out.println("Modular      : " + TimeUnits.nanosecondsToString(elapsed));
                System.out.println();

            }
        }
        System.out.println("Subresultant : " + TimeUnits.statisticsNanotime(subresultant));
        System.out.println("Modular      : " + TimeUnits.statisticsNanotime(modular));
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