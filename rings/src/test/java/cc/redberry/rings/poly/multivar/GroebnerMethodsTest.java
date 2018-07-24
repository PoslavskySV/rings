package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Ignore;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.multivar.GroebnerBases.GroebnerBasis;
import static cc.redberry.rings.poly.multivar.GroebnerBases.optimalOrder;
import static cc.redberry.rings.poly.multivar.GroebnerBasesTest.SingularGB;
import static cc.redberry.rings.poly.multivar.GroebnerMethods.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.util.TimeUnits.nanosecondsToString;
import static org.junit.Assert.*;

/**
 *
 */
public class GroebnerMethodsTest extends AMultivariateTest {

    @Test
    public void test1() {
        List<DegreeVector> dvs;

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}),
                new DegreeVector(new int[]{3, 4, 5}),
                new DegreeVector(new int[]{9, 8, 7}),
                new DegreeVector(new int[]{1, -1, -1}));
        assertTrue(algebraicallyDependentMonomialsQ(dvs));

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}));
        assertFalse(algebraicallyDependentMonomialsQ(dvs));

        dvs = Arrays.asList(
                new DegreeVector(new int[]{1, 2, 3}),
                new DegreeVector(new int[]{2, 1, 2}),
                new DegreeVector(new int[]{3, 1, 2}),
                new DegreeVector(new int[]{4, 1, 2}),
                new DegreeVector(new int[]{5, 1, 2}));
        assertTrue(algebraicallyDependentMonomialsQ(dvs));
    }

    @Test
    public void test2() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(3, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        assertFalse(probablyAlgebraicallyDependentQ(Arrays.asList(coder.parse("x"), coder.parse("y"))));
        assertFalse(probablyAlgebraicallyDependentQ(Arrays.asList(coder.parse("x^2 + x"), coder.parse("x^2 + y"))));
    }

    @Ignore
    @Test
    public void test3() {
        List<MultivariatePolynomial<Rational<BigInteger>>> katsura = GroebnerBasesData.katsura(8);
        eliminate(katsura.stream().map(p -> p.mapCoefficients(Zp64(17), r -> r.numerator().longValueExact())).collect(Collectors.toList()), 0)
                .stream().sorted(Comparator.comparingInt(AMultivariatePolynomial::nUsedVariables))
                .forEach(System.out::println);
    }

    @Test
    public void test4() {
        MultivariateRing<MultivariatePolynomialZp64> ring = Rings.MultivariateRingZp64(3, Zp64(32003));
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "z");

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(
                coder.parse("x^2"),
                coder.parse("x*y"),
                coder.parse("y^5")
        );

        assertTrue(Ideal.create(eliminate(ideal, 0)).isPrincipal());
    }

    @Test
    public void test5() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(5, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "t", "s", "z");

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(
                coder.parse("x - t"),
                coder.parse("y - t^2"),
                coder.parse("z - t^3"),
                coder.parse("s - x + y^3")
        );

        Ideal eliminate = Ideal.create(eliminate(ideal, 2, 3), GREVLEX);
        Ideal el1 = Ideal.create(eliminate(eliminate(ideal, 2), 3), GREVLEX);
        Ideal el2 = Ideal.create(eliminate(eliminate(ideal, 3), 2), GREVLEX);

        assertEquals(eliminate, el1);
        assertEquals(eliminate, el2);
    }

    @Test
    public void test6() {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(3, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "z");

        List<MultivariatePolynomial<Rational<BigInteger>>> sys;

        sys = Arrays.asList(
                coder.parse("x - y"),
                coder.parse("x + y"),
                coder.parse("x^2 + y^2"));
        assertNonEmptyAnnihilators(sys, algebraicRelations(sys));

        assertTrue(algebraicRelations(Arrays.asList(
                coder.parse("x * z"),
                coder.parse("y * z"))).isEmpty());

        sys = Arrays.asList(
                coder.parse("x * y + z^2"),
                coder.parse("z^2 + y^2"),
                coder.parse("x^2*y^2 -2 * x * y^3 + y^4"));
        assertNonEmptyAnnihilators(sys, algebraicRelations(sys));
    }

    @Test
    @RequiresSingular
    public void test7() throws IOException, InterruptedException {
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(5, Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "t", "s", "z");

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(
                coder.parse("x - t"),
                coder.parse("y - t^2"),
                coder.parse("z - t^3"),
                coder.parse("s - x + y^3")
        );

        assertEquals(
                SingularGB(ideal, GREVLEX, 2, 3),
                GroebnerBasis(eliminate(ideal, 2, 3), GREVLEX));
    }

    @Ignore
    @Test
    @RequiresSingular
    public void test8_eliminate_singular_random() throws IOException, InterruptedException {
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(100, 1);
        int nEls = 4;
        DescriptiveStatistics
                tF4 = new DescriptiveStatistics(),
                tSingular = new DescriptiveStatistics();
        for (int i = 0; i < nIterations; i++) {
            System.out.println();
            List<MultivariatePolynomialZp64> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++)
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(4, 3, 4, domain, GREVLEX, rnd));

            System.out.println(ideal);

            int eliminateVar = rnd.nextInt(4);
            System.out.println("var: " + eliminateVar);

            GroebnerBasesTest.SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX, eliminateVar);
            tSingular.addValue(singular.nanoseconds);
            System.out.println("Empty:" + singular.list.isEmpty());
            System.out.println("Singular : " + nanosecondsToString(singular.nanoseconds));

            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> eliminate = eliminate(ideal, eliminateVar);
            System.out.println("Rings (eliminate): " + nanosecondsToString(System.nanoTime() - start));
            List<MultivariatePolynomialZp64> actual = GroebnerBasis(eliminate, GREVLEX);
            long time = System.nanoTime() - start;
            tF4.addValue(time);
            System.out.println("Rings       : " + nanosecondsToString(time));


            List<MultivariatePolynomialZp64> expected = singular;
            if (!actual.equals(expected)) {
                System.out.println(ideal);
                System.out.println(actual.size());
            }
            assertEquals(expected, actual);
        }

        System.out.println("F4       : " + TimeUnits.statisticsNanotime(tF4));
        System.out.println("Singular : " + TimeUnits.statisticsNanotime(tSingular));
    }

    @Test
    @RequiresSingular
    public void test8() throws IOException, InterruptedException {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(4, 17);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x1", "x2", "x3", "x4");

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(
                coder.parse("8*x1^3*x2^2+6*x1^3*x4^3"),
                coder.parse("16*x2^2*x4+11*x1*x2^2*x4^2+15*x1*x3^3*x4+5*x1*x3^2*x4^3"),
                coder.parse("8*x2^2*x3^2+2*x2^2*x4^3+2*x1^3*x2^2+12*x1*x2^3*x3^3*x4^2"),
                coder.parse("2*x1^3*x2*x3+12*x2^3*x3*x4^2+11*x1*x2^3*x4^3+12*x1*x2^2*x3^2*x4^2"));

        for (int eliminateVar = 0; eliminateVar < 4; eliminateVar++) {
            GroebnerBasesTest.SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX, eliminateVar);
            System.out.println("Singular : " + nanosecondsToString(singular.nanoseconds));

            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> eliminate = eliminate(ideal, eliminateVar);
            System.out.println("Rings (eliminate): " + nanosecondsToString(System.nanoTime() - start));
            List<MultivariatePolynomialZp64> actual = GroebnerBasis(eliminate, GREVLEX);
            long time = System.nanoTime() - start;
            System.out.println("Rings       : " + nanosecondsToString(time));


            List<MultivariatePolynomialZp64> expected = singular;
            if (!actual.equals(expected)) {
                System.out.println(ideal);
                System.out.println(actual.size());
            }
            assertEquals(expected, actual);
        }
    }

    @Test
    public void testNullstellensatzCertificate1() {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(3, 17);
        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkPolynomialCoder(ring, "x", "y", "z");

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(
                coder.parse("x^2 + 1"),
                coder.parse("z + x * y + 2"),
                coder.parse("z^2 - x * y + 2"),
                coder.parse("z - x * y + x + y"));

        assertNullstellensatzCertificate(ideal, NullstellensatzCertificate(ideal, true));
        assertNullstellensatzCertificate(ideal, NullstellensatzCertificate(ideal, false));
    }

    @Test
    public void testNullstellensatzCertificate2_random() {
        RandomGenerator rnd = getRandom();

        testNullstellensatzCertificateRandom(MultivariateRingZp64(3, 17).getOne(), 3, 100, rnd);
        testNullstellensatzCertificateRandom(MultivariateRingZp64(4, 17).getOne(), 2, 10, rnd);

        // this is too long for now due to expression swell in coefficients; fixme with Bareiss
        // testNullstellensatzCertificateRandom(MultivariateRing(3, Q).getOne(), 3, 10, rnd);
        // testNullstellensatzCertificateRandom(MultivariateRing(4, Q).getOne(), 2, 10, rnd);
    }


    @Test
    public void testLeinartDecomposition1() {
        MultivariateRing<MultivariatePolynomialZp64> mRing = MultivariateRingZp64(3, 17);
        Coder<MultivariatePolynomialZp64, ?, ?> mCoder = Coder.mkPolynomialCoder(mRing, "x", "y", "z");
        Rationals<MultivariatePolynomialZp64> fRing = Frac(mRing);
        Coder<Rational<MultivariatePolynomialZp64>, ?, ?> fCoder = Coder.mkRationalsCoder(fRing, mCoder);

        Rational<MultivariatePolynomialZp64> f = fCoder.parse("(x + y) / (x^2 + y^2) / (x^3 - x * y - 1) / (x - y)");
        List<Rational<MultivariatePolynomialZp64>> decomposition = LeinartDecomposition(f);

        assertTrue(f.subtract(decomposition.stream().reduce(fRing.getZero(), fRing::add)).isZero());
    }

    @SuppressWarnings("unchecked")
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly toZ(Poly p) {
        return p.isOverZ()
                ? (Poly) ((MultivariatePolynomial<Rational<BigInteger>>) p).mapCoefficients(Q, c -> Q.mkNumerator(c.numerator().mod(BigInteger.FIVE)))
                : p;
    }

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void testNullstellensatzCertificateRandom(Poly factory, int varDegB, int nIterations, RandomGenerator rnd) {
        DescriptiveStatistics
                boundDeg = new DescriptiveStatistics(),
                boundVar = new DescriptiveStatistics(),
                certBound = new DescriptiveStatistics();
        for (int iter = 0; iter < nIterations; ++iter) {
            List<Poly> ideal = new ArrayList<>();
            for (int i = 0; i < factory.nVariables + 2; ++i) {
                Poly p = toZ(RandomMultivariatePolynomials.randomPolynomial(factory, varDegB, 1 + rnd.nextInt(5), rnd));
                ideal.add(p);
            }
            Ideal<Term, Poly> gb = Ideal.create(ideal, optimalOrder(ideal));
            while (!gb.isTrivial()) {
                Poly p = toZ(RandomMultivariatePolynomials.randomPolynomial(factory, varDegB, 1 + rnd.nextInt(5), rnd));
                ideal.add(p);
                gb = gb.union(p);
            }

            long start;

            start = System.nanoTime();
            List<Poly> certificateA = NullstellensatzCertificate(ideal, true);
            assertNullstellensatzCertificate(ideal, certificateA);
            boundDeg.addValue(System.nanoTime() - start);

            certificateA.forEach(p -> certBound.addValue(p.degree()));

            start = System.nanoTime();
            List<Poly> certificateB = NullstellensatzCertificate(ideal, false);
            assertNullstellensatzCertificate(ideal, certificateB);
            boundVar.addValue(System.nanoTime() - start);

        }
        System.out.println("Bound total degree: " + TimeUnits.statisticsNanotime(boundDeg));
        System.out.println("Bound vars degree:  " + TimeUnits.statisticsNanotime(boundVar));
        System.out.println("Certificate degrees:  " + certBound.getMin() + " <= " + certBound.getMean() + " <= " + certBound.getMax());
        System.out.println();
    }


    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void assertAnnihilators(List<Poly> initialPolys, List<Poly> algebraicRelations) {
        assertTrue(algebraicRelations.stream().map(p -> p.composition(initialPolys)).allMatch(AMultivariatePolynomial::isZero));
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void assertNonEmptyAnnihilators(List<Poly> initialPolys, List<Poly> algebraicRelations) {
        assertFalse(algebraicRelations.isEmpty());
        assertAnnihilators(initialPolys, algebraicRelations);
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void assertNullstellensatzCertificate(List<Poly> initialPolys, List<Poly> certificate) {
        assertTrue(IntStream.range(0, initialPolys.size())
                .mapToObj(i -> initialPolys.get(i).clone().multiply(certificate.get(i)))
                .reduce(initialPolys.get(0).createZero(), (a, b) -> a.add(b)).isOne());
    }
}