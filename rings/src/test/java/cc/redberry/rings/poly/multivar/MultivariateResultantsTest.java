package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.AlgebraicNumberField;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.RandomDataGenerator;
import cc.redberry.rings.util.TimeUnits;
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
import java.util.function.Supplier;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
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
        MultivariatePolynomial<BigInteger> b = coder.parse("(x - 3*y - z)^3 + 1");
        MultivariatePolynomial<BigInteger> expected = coder.parse("343 - 124509*y^3 + 2470629*y^6 - 40353607*y^9 - 160083*y^2*z + 6353046*y^5*z - 155649627*y^8*z - 68607*y*z^2 + 6806835*y^4*z^2 - 266827932*y^7*z^2 - 9801*z^3 + 3889620*y^3*z^3 - 266827932*y^6*z^3 + 1250235*y^2*z^4 - 171532242*y^5*z^4 + 214326*y*z^5 - 73513818*y^4*z^5 + 15309*z^6 - 21003948*y^3*z^6 - 3857868*y^2*z^7 - 413343*y*z^8 - 19683*z^9");
        assertEquals(expected, ClassicalResultant(a, b, 0));
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
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd);
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

    @Ignore
    @Test
    public void test5_bivariate_performance() throws Exception {
        IntegersZp64 cfRing = Zp64(1048583);
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(2, cfRing);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y");

        MultivariatePolynomialZp64 a = coder.parse("(2*x + y + 2)^47 + 12 + (x - 1)^46 + (y + 2)^43");
        MultivariatePolynomialZp64 b = coder.parse("(x - 3*y - 3)^43 + 13 + x^45*y^45");
        for (int i = 0; i < 100; ++i) {
            long start;

            start = System.nanoTime();
            MultivariatePolynomialZp64 br = BrownResultant(a, b, 0);
            System.out.println("Brown: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            MultivariatePolynomialZp64 zp = ZippelResultant(a, b, 0);
            System.out.println("Zippel: " + nanosecondsToString(System.nanoTime() - start));
            assertEquals(br, zp);

            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomialZp64>
                    aUni = a.asUnivariate(0).mapCoefficients(UnivariateRingZp64(cfRing), p -> p.asUnivariate()),
                    bUni = b.asUnivariate(0).mapCoefficients(UnivariateRingZp64(cfRing), p -> p.asUnivariate());
            MultivariatePolynomialZp64 uniRes = UnivariateResultants.Resultant(aUni, bUni).asMultivariate(ring.ordering()).insertVariable(0);
            assertEquals(br, uniRes);
            System.out.println("Univar: " + nanosecondsToString(System.nanoTime() - start));
            System.out.println(br.size() + "   " + br.degree());

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
                rEval = expected.dropVariable(0).evaluate(1, 1);


        SparseInterpolationZp64 interpolation =
                createInterpolation(1,
                        a.asUnivariateEliminate(0),
                        b.asUnivariateEliminate(0),
                        rEval, expected.degree(2), getRandom());

        assertEquals(expected.evaluate(2, 2),
                interpolation.evaluate(2).insertVariable(0));
    }

    @Test
    public void testZippel1_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        IntegersZp64 ring = Zp64(SmallPrimes.nextPrime(1 << 20));
        DescriptiveStatistics
                zippel = new DescriptiveStatistics(),
                brown = new DescriptiveStatistics();
        int nIterations = its(100, 100);
        int nVars = 4, degree = 7;
        for (int i = 0; i < nIterations; ++i) {
            MultivariatePolynomialZp64
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd);
            int variable = rnd.nextInt(nVars);

            long start;
            start = System.nanoTime();
            MultivariatePolynomialZp64 actual = ZippelResultant(a, b, variable);
            zippel.addValue(System.nanoTime() - start);


            start = System.nanoTime();
            MultivariatePolynomialZp64 expected = BrownResultant(a, b, variable);
            brown.addValue(System.nanoTime() - start);

            assertEquals(expected, actual);
        }

        System.out.println("Zippel's " + zippel);
        System.out.println("Brown's " + brown);
    }

    @Test
    public void testZippel1() {
        IntegersZp64 cfRing = Zp64(SmallPrimes.nextPrime(1 << 20));
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(4, cfRing);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4");

        MultivariatePolynomialZp64
                a = coder.parse("1064*x2^4*x3^2+813059*x2^5*x4+982139*x1*x2^2+772988*x1^2+966476*x1^2*x2*x3^2*x4+207078*x1^3*x3*x4^2+1037301*x1^3*x2^3+1041799*x1^5+37593*x1^6"),
                b = coder.parse("698234*x3^4*x4^2+847486*x2*x3^2*x4+499885*x2^4*x3+748942*x2^6+537106*x1*x3^5+703951*x1*x2^3*x3+131882*x1*x2^5+87996*x1^2*x3^2*x4+962154*x1^2*x2^3*x3+55977*x1^3*x2*x3^2+16786*x1^4*x3+694160*x1^4*x2^2+554118*x1^5*x3+966084*x1^5*x2");
        assertEquals(BrownResultant(a, b, 0), ZippelResultant(a, b, 0));
    }

    @Test
    public void testZippel2() {
        IntegersZp64 cfRing = Zp64(SmallPrimes.nextPrime(1 << 20));
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(4, cfRing);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4");

        MultivariatePolynomialZp64
                a = coder.parse("512502*x3^5+59504*x3^6+244225*x1*x2^5+530513*x1^2*x2*x3+435564*x1^2*x2^3*x3+546487*x1^3*x2^2*x3+317466*x1^4*x2^2+111459*x1^5*x2+146247*x1^6"),
                b = coder.parse("681651*x1^2+400314*x1^2*x2^4+481090*x1^4*x2*x3");
        assertEquals(BrownResultant(a, b, 0), ZippelResultant(a, b, 0));
    }

    @Test
    public void testZippel3_random_sparse_vars() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        IntegersZp64 ring = Zp64(SmallPrimes.nextPrime(1 << 20));
        DescriptiveStatistics
                zippel = new DescriptiveStatistics(),
                brown = new DescriptiveStatistics();
        int nIterations = its(50, 100);
        int nVars = 100, degree = 7;
        for (int i = 0; i < nIterations; ++i) {
            MultivariatePolynomialZp64
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd);
            while (a.nUsedVariables() >= 5)
                a = a.evaluate(rnd.nextInt(nVars), 1);
            while (b.nUsedVariables() >= 5)
                b = b.evaluate(rnd.nextInt(nVars), 1);
            int variable;
            do { variable = rnd.nextInt(nVars);} while (a.degree(variable) == 0 && b.degree(variable) == 0);

            long start;
            start = System.nanoTime();
            MultivariatePolynomialZp64 actual = ZippelResultant(a, b, variable);
            zippel.addValue(System.nanoTime() - start);
            System.out.println("Zippel: " + nanosecondsToString(System.nanoTime() - start));


            start = System.nanoTime();
            MultivariatePolynomialZp64 expected = BrownResultant(a, b, variable);
            brown.addValue(System.nanoTime() - start);
            System.out.println("Brown: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(expected, actual);
            System.out.println();
        }

        System.out.println("Zippel's " + zippel);
        System.out.println("Brown's " + brown);
    }

    @Test
    public void testZippel4_generic() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        IntegersZp64 ring = Zp64(SmallPrimes.nextPrime(1 << 20));
        int nIterations = its(25, 50);
        int nVars = 4, degree = 5;
        for (int i = 0; i < nIterations; ++i) {
            MultivariatePolynomialZp64
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), ring, LEX, rnd);
            int variable = rnd.nextInt(nVars);

            MultivariatePolynomialZp64 expected = ZippelResultant(a, b, variable);

            MultivariatePolynomial<BigInteger> g;
            g = BrownResultant(a.toBigPoly(), b.toBigPoly(), variable);
            assertEquals(expected, MultivariatePolynomial.asOverZp64(g));
            g = ZippelResultant(a.toBigPoly(), b.toBigPoly(), variable);
            assertEquals(expected, MultivariatePolynomial.asOverZp64(g));
        }
    }

    @Test
    public void testZippel4_generic_a() throws Exception {
        IntegersZp64 ring = Zp64(SmallPrimes.nextPrime(1 << 20));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("642521*x1*x3+544253*x2", ring, "x1", "x2", "x3", "x4"),
                b = MultivariatePolynomialZp64.parse("653650*x1^4+667997*x2^4", ring, "x1", "x2", "x3", "x4");

        int variable = 0;
        MultivariatePolynomialZp64 expected = MultivariatePolynomialZp64.parse("90619*x2^4 + 797965*x2^4*x3^4", ring, "x1", "x2", "x3", "x4");
        assertEquals(expected, BrownResultant(a, b, variable));
        assertEquals(expected, ZippelResultant(a, b, variable));
    }

    @Ignore
    @Test
    public void testZippel5_dense() throws Exception {
        IntegersZp64 cfRing = Zp64(1048583);
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(4, cfRing);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4");

        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);

        MultivariatePolynomialZp64
                a = ring.getOne(), b = ring.getOne();
        for (int i = 0; i < 1; ++i) {
            a.multiply(RandomMultivariatePolynomials.randomPolynomial(ring.nVariables(), 5, 100 * 5 * 5 * 5 * 5 * 5, cfRing, LEX, rnd));
            b.multiply(RandomMultivariatePolynomials.randomPolynomial(ring.nVariables(), 5, 100 * 5 * 5 * 5 * 5 * 5, cfRing, LEX, rnd));
        }

        System.out.println(a.sparsity());
        System.out.println(b.sparsity());

        for (int i = 0; i < 1000; ++i) {
            long start;
            start = System.nanoTime();
            MultivariatePolynomialZp64 actual = ZippelResultant(a, b, 0);
            System.out.println("Zippel: " + nanosecondsToString(System.nanoTime() - start));


            start = System.nanoTime();
            MultivariatePolynomialZp64 expected = BrownResultant(a, b, 0);
            System.out.println("Brown: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(expected, actual);

            start = System.nanoTime();
            MultivariatePolynomialZp64 classic = ClassicalResultant(a, b, 0);
            System.out.println("Classical: " + nanosecondsToString(System.nanoTime() - start));

            System.out.println();
        }
    }


    @Test
    public void testModular1() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        MultivariatePolynomial<BigInteger> a = coder.parse("(2*x + y + z)^13 + 1");
        MultivariatePolynomial<BigInteger> b = coder.parse("(x - 3*y - z)^3 + 1");
        assertEquals(ClassicalResultant(a, b, 0), ModularResultantInZ(a, b, 0));
    }

    @Test
    @RequiresSingular
    public void testModular2_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(25, 100);
        int nVars = 4, degree = 7;
        DescriptiveStatistics
                modular = new DescriptiveStatistics(),
                singular = new DescriptiveStatistics();
        int cfBound = 10_000;
        for (int i = 0; i < nIterations; ++i) {
            MultivariatePolynomial<BigInteger>
                    a = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), Zp64(cfBound), LEX, rnd).toBigPoly().setRing(Z),
                    b = RandomMultivariatePolynomials.randomSharpPolynomial(nVars, degree, rndd.nextInt(1, 20), Zp64(cfBound), LEX, rnd).toBigPoly().setRing(Z);
            int variable = rnd.nextInt(nVars);

            try {
                SingularResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> si = SingularResultant(a, b, variable);
                MultivariatePolynomial<BigInteger> expected = null;
                if (si != null) {
                    expected = si.resultant;
                    singular.addValue(si.nanotime);
                    System.out.println("si: " + nanosecondsToString(si.nanotime));
                }

                long start;
                start = System.nanoTime();
                MultivariatePolynomial<BigInteger> actual = ModularResultantInZ(a, b, variable);
                modular.addValue(System.nanoTime() - start);
                System.out.println("ri: " + nanosecondsToString(System.nanoTime() - start));
                System.out.println();
                if (expected != null)
                    assertEquals(expected, actual);
            } catch (Throwable e) {
                System.out.println(a);
                System.out.println(b);
                System.out.println(variable);
                throw new RuntimeException(e);
            }
        }

        System.out.println("Modular's " + modular);
        System.out.println("Singulars's " + singular);
    }

    @Test
    public void testModular3() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(4, Z, LEX);
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4");

        MultivariatePolynomial<BigInteger> a = coder.parse("94*x2*x4^2+56*x2*x3*x4^2+27*x2^3*x4+27*x1*x2+54*x1*x2*x3^2+60*x1^2*x3*x4+36*x1^2*x2^2");
        MultivariatePolynomial<BigInteger> b = coder.parse("49*x3^3+43*x1*x4+63*x1*x3*x4+14*x1*x2+9*x1*x2*x3^2+44*x1*x2^2+3*x1^2*x3^2+77*x1^3+56*x1^3*x4+83*x1^3*x3+95*x1^3*x2+94*x1^4");

        int variable = 1;
        MultivariatePolynomial<BigInteger> ord = ClassicalResultant(a, b, variable);
        MultivariatePolynomial<BigInteger> mod = ModularResultantInZ(a, b, variable);
        assertEquals(ord, mod);
    }

    @Test
    @RequiresSingular
    public void testModular4() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(4, Z, LEX);
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x1", "x2", "x3", "x4");

        MultivariatePolynomial<BigInteger> a = coder.parse("9789*x1^3*x3*x4+5618*x1^4*x4+5574*x1^4*x2*x3+9813*x1^6");
        MultivariatePolynomial<BigInteger> b = coder.parse("3425*x2^2*x4+2575*x2^3*x3^2+1157*x2^4*x3+6413*x1*x4^4+6401*x1*x2^2*x3^2+296*x1^2+42*x1^2*x2*x4^3+3941*x1^2*x2*x3+3918*x1^3*x3^3+5245*x1^3*x2*x3^2+1114*x1^4*x3^2+8801*x1^4*x2+773*x1^4*x2^2+2744*x1^5*x3+7473*x1^5*x2+7395*x1^6");

        int variable = 0;
        MultivariatePolynomial<BigInteger> sin = SingularResultant(a, b, variable).resultant;
        MultivariatePolynomial<BigInteger> mod = ModularResultantInZ(a, b, variable);
        assertEquals(sin, mod);
    }

    @Test
    public void testInQ1() throws Exception {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(3, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        MultivariatePolynomial<Rational<BigInteger>> a = coder.parse("(2*x/3 + y/4 + z)^3 + 1/5");
        MultivariatePolynomial<Rational<BigInteger>> b = coder.parse("(x - 3*y/7 - z)^3 + 1/2");
        MultivariatePolynomial<Rational<BigInteger>> expected = coder.parse("-343/2460375 - (24545*y^3)/197568 - (253125*y^6)/68841472 - (38443359375*y^9)/10578455953408 - (24545*y^2*z)/21168 - (84375*y^5*z)/1229312 - (38443359375*y^8*z)/377801998336 - (24545*y*z^2)/6804 - (46875*y^4*z^2)/87808 - (4271484375*y^7*z^2)/3373232128 - (24545*z^3)/6561 - (15625*y^3*z^3)/7056 - (158203125*y^6*z^3)/17210368 - (15625*y^2*z^4)/3024 - (52734375*y^5*z^4)/1229312 - (3125*y*z^5)/486 - (5859375*y^4*z^5)/43904 - (21875*z^6)/6561 - (1953125*y^3*z^6)/7056 - (1953125*y^2*z^7)/5292 - (1953125*y*z^8)/6804 - (1953125*z^9)/19683");
        assertEquals(expected, ResultantInQ(a, b, 0));
    }

    @Test
    public void testInNumberField1() throws Exception {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> cfRing
                = new AlgebraicNumberField<>(Coder.mkUnivariateCoder(UnivariateRing(Q), "x").parse("-2 + 17*x^2"));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(cfRing, "s");

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> ring = MultivariateRing(3, cfRing);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, cfCoder, "x", "y", "z");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                a = coder.parse("s*y*z^2 + x^2*z*s/3 - 2"),
                b = coder.parse("s + x*s/3 + z^2*y*x^3 + y + z");
        int variable = 0;
        assertEquals(ClassicalResultant(a, b, variable), ModularResultantInNumberField(a, b, variable));

        a = a.clone().multiply(b.clone().increment()).increment();
        b = b.clone().multiply(a.clone().increment()).decrement();
        assertEquals(ClassicalResultant(a, b, variable), ModularResultantInNumberField(a, b, variable));
    }

    @Test
    public void testInNumberField2() throws Exception {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> cfRing
                = new AlgebraicNumberField<>(Coder.mkUnivariateCoder(UnivariateRing(Q), "x").parse("16471497/10+8378273/5*x+x^2"));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(cfRing, "s");

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> ring = MultivariateRing(4, cfRing);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, cfCoder, "x1", "x2", "x3", "x4");

        MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                a = coder.parse("(-65885963/50-33513087/25*s)*x3+(-13193646697/3600-745664197/200*s)*x1"),
                b = coder.parse("(-8235746/5-75404437/45*s)*x3+(-91488124937/8400-139606973399/12600*s)*x1");
        int variable = 0;
        assertEquals(ClassicalResultant(a, b, variable), ModularResultantInNumberField(a, b, variable));

        a = a.clone().multiply(b.clone().increment()).increment();
        b = b.clone().multiply(a.clone().increment()).decrement();
        assertEquals(ClassicalResultant(a, b, variable), ModularResultantInNumberField(a, b, variable));
    }

    @Test
    public void testInNumberField_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                subresultant = new DescriptiveStatistics(),
                modular = new DescriptiveStatistics();

        int nVars = 4, degree = 7;
        for (int i = 0; i < its(30, 100); ++i) {
            UnivariatePolynomial<Rational<BigInteger>> minimalPoly = IrreduciblePolynomials
                    .randomIrreduciblePolynomialOverZ(rndd.nextInt(2, 5), rnd)
                    .mapCoefficients(Q, Q::mkNumerator);

            minimalPoly.setLC(Q.mkNumerator(rndd.nextInt(1, 10)));
            if (!IrreduciblePolynomials.irreducibleQ(minimalPoly)) {
                --i;
                continue;
            }

            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = new AlgebraicNumberField<>(minimalPoly);
            System.out.println(i);
            System.out.println("Field        : " + field);
            Supplier<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> rndPoly = () ->
                    RandomMultivariatePolynomials.randomSharpPolynomial(
                            nVars,
                            rndd.nextInt(2, 4),
                            rndd.nextInt(degree / 2, degree),
                            field,
                            GREVLEX,
                            __ -> RandomUnivariatePolynomials.randomPoly(minimalPoly.degree(), Q, ___ -> Q.mk(rnd.nextInt(10), 1 + rnd.nextInt(10)), rnd),
                            rnd);

            Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");
            Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(MultivariateRing(nVars, field), cfCoder, "x1", "x2", "x3", "x4");

            for (int variable = 0; variable < nVars; ++variable) {
                MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                        a = rndPoly.get(),
                        b = rndPoly.get();

//                System.out.println(coder.stringify(a));
//                System.out.println(coder.stringify(b));
//                System.out.println(variable);
//                System.out.println();

                long start, elapsed;

                start = System.nanoTime();
                MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> expected = ClassicalResultant(a, b, variable);
                elapsed = System.nanoTime() - start;
                subresultant.addValue(elapsed);
                System.out.println("Subresultant : " + TimeUnits.nanosecondsToString(elapsed));

                start = System.nanoTime();
                MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> actual = ModularResultantInNumberField(a, b, variable);
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

    @Test
    public void test9() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> cfRing = MultivariateRing(2, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkMultivariateCoder(cfRing, "c", "d");

        MultivariateRing<MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> ring = MultivariateRing(2, cfRing);
        Coder<MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, cfCoder, "x", "y");

        MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> a = coder.parse("2 + x^2");
        MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> b = coder.parse("3 + y^2");
        assertEquals(coder.parse("(2 + x^2)^2"), Resultant(a, b, 1));
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