package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.multivar.GroebnerBases.*;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.test.TimeConsuming;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.ListWrapper;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;
import java.nio.file.Files;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.multivar.GroebnerBases.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;
import static cc.redberry.rings.primes.SmallPrimes.nextPrime;
import static cc.redberry.rings.util.TimeConstrained.timeConstrained;
import static cc.redberry.rings.util.TimeUnits.nanosecondsToString;
import static cc.redberry.rings.util.TimeUnits.statisticsNanotime;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class GroebnerBasesTest extends AMultivariateTest {
    @Before
    public void beforeMethod() throws Exception {
        if (getClass().getMethod(name.getMethodName()).isAnnotationPresent(RequiresSingular.class))
            Assume.assumeTrue(isSingularAvailable());
        super.beforeMethod();
    }

    private static String getSingularPath() {
        return System.getProperty("singularPath", "/Applications/Singular.app/Contents/bin/Singular");
    }

    private static boolean isSingularAvailable() {
        return new File(getSingularPath()).exists();
    }

    @Test
    public void test1() throws Exception {
        String[] vars = {"u0", "u1", "u2", "u3"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("u3*u3 + u2*u2 + u1*u1 + u0*u0 + u1*u1 + u2*u2 + u3*u3 - u0", Q, GREVLEX, vars),
                b = parse("u3*0 + u2*0 + u1*u3 + u0*u2 + u1*u1 + u2*u0 + u3*u1 - u2", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> expected = new ArrayList<>(Arrays.asList(
                parse("u1^2 - u2 + 2*u0*u2 + 2*u1*u3", Q, GREVLEX, vars),
                parse("-u0 + u0^2 + 2*u2 - 4*u0*u2 + 2*u2^2 - 4*u1*u3 + 2*u3^2", Q, GREVLEX, vars)
        ));

        canonicalize(expected);
        List<MultivariatePolynomial<Rational<BigInteger>>>
                buchberger = BuchbergerGB(Arrays.asList(a, b), GREVLEX),
                f4 = F4GB(Arrays.asList(a, b), GREVLEX);
        assertEquals(expected, buchberger);
        assertEquals(expected, f4);
    }

    @Test
    public void test2() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("x1^2*x2^2*x3 + x5^4 - 1 + x2^3*x4 + x3^5 + x4", Q, GREVLEX, vars),
                b = parse("x1*x2*x3^2 + x2*x4*x1^2*x5 + x3*x2^3", Q, GREVLEX, vars),
                c = parse("x1^3*x2^3*x3^3*x4 - x2*x4^2*x1^4*x5 + x3^3*x2 - x4 - 1", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Stream.of(a, b, c)
                .map(p -> p.homogenize(vars.length)).collect(Collectors.toList());
        Assert.assertTrue(gens.stream().allMatch(MultivariatePolynomial::isHomogeneous));

        List<MultivariatePolynomial<Rational<BigInteger>>> bHom = BuchbergerGB(gens, GREVLEX,
                () -> new GradedSyzygyTreeSet<>(new TreeMap<>(), defaultSelectionStrategy(GREVLEX), SyzygyPair::degree),
                null);
        Assert.assertTrue(bHom.stream().allMatch(MultivariatePolynomial::isHomogeneous));

        List<MultivariatePolynomial<Rational<BigInteger>>> bNonHom = BuchbergerGB(gens, GREVLEX);
        Assert.assertTrue(bNonHom.stream().allMatch(MultivariatePolynomial::isHomogeneous));

        assertEquals(bNonHom, bHom);

        List<MultivariatePolynomial<Rational<BigInteger>>> f4 = F4GB(gens, GREVLEX);
        Assert.assertTrue(bNonHom.stream().allMatch(MultivariatePolynomial::isHomogeneous));

        assertEquals(bNonHom, f4);
    }

    @Test
    @RequiresSingular
    public void test2a() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("x1^2*x2^2*x3 + x5^4 - 1 + x2^3*x4 + x3^5 + x4", Q, GREVLEX, vars),
                b = parse("x1*x2*x3^2 + x2*x4*x1^2*x5 + x3*x2^3", Q, GREVLEX, vars),
                c = parse("x1^3*x2^3*x3^3*x4 - x2*x4^2*x1^4*x5 + x3^3*x2 - x4 - 1", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Stream.of(a, b, c)
                .map(p -> p.homogenize(vars.length)).collect(Collectors.toList());

        for (int i = 0; i < 2; i++) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> hm = BuchbergerGB(gens, GREVLEX,
                    () -> new GradedSyzygyTreeSet<>(new TreeMap<>(), defaultSelectionStrategy(GREVLEX), SyzygyPair::degree),
                    null);
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(hm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> nhm = BuchbergerGB(gens, GREVLEX);
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(nhm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            List<MultivariatePolynomial<Rational<BigInteger>>> ref = SingularGB(gens, GREVLEX);
            Assert.assertEquals(hm, nhm);
            Assert.assertEquals(ref, nhm);
        }

//        788ms
//        482ms
//        389ms
//        419ms
    }

    @Test
    public void test3() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp64 ring = new IntegersZp64(17);
        MultivariatePolynomialZp64
                f1 = MultivariatePolynomialZp64.parse("x^2*y^2 + x*y + 5*z^3*y^2", ring, GREVLEX, vars),
                f2 = MultivariatePolynomialZp64.parse("x*y^4 - y^2 - 5*z^3*x^2", ring, GREVLEX, vars);

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2);
        assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
        List<MultivariatePolynomial<BigInteger>> idealE = ideal.stream().map(MultivariatePolynomialZp64::toBigPoly).collect(Collectors.toList());
        assertEquals(BuchbergerGB(idealE, GREVLEX), F4GB(idealE, GREVLEX));
    }

    @Test
    public void test4() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2 + x*y - 1", Q, GREVLEX, vars),
                f2 = parse("x^2 - z^2", Q, GREVLEX, vars),
                f3 = parse("x*y + 1", Q, GREVLEX, vars);

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(f1, f2, f3);
        assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
    }

    @Test
    public void test5() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp64 ring = new IntegersZp64(17);
        MultivariatePolynomialZp64
                f1 = MultivariatePolynomialZp64.parse("x^2*y^2 + x*y + 5*z^13*y^2", ring, GREVLEX, vars),
                f2 = MultivariatePolynomialZp64.parse("x*y^24 - y^2*z^6 - 5*z^3*x^2 - 1", ring, GREVLEX, vars);

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2);
        assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
    }

    @Test
    @RequiresSingular
    public void test6_katsura() throws Exception {
        IntegersZp64 ring = new IntegersZp64(17);
        for (int i = 3; i <= its(8, 9); i++) {
            System.out.println(String.format("=> Katsura%s:", i));
            List<MultivariatePolynomialZp64> ideal =
                    GroebnerBasisData.katsura(i)
                            .stream()
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator())))
                            .map(p -> p.setOrdering(GREVLEX))
                            .collect(Collectors.toList());

            long start;

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = singular;
            System.out.println("   Singular  : " + nanosecondsToString(singular.nanoseconds));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + nanosecondsToString(f4));
            assertEquals(expected, actualF4);

            if (i <= 8) {
                List<MultivariatePolynomialZp64> actualE =
                        F4GB(ideal.stream().map(MultivariatePolynomialZp64::toBigPoly).collect(Collectors.toList()), GREVLEX)
                                .stream().map(MultivariatePolynomial::asOverZp64).collect(Collectors.toList());
                assertEquals(expected, actualE);
            }
        }
    }

    @Test
    @RequiresSingular
    public void test6_cyclic() throws Exception {
        IntegersZp64 ring = new IntegersZp64(17);
        for (int i = 5; i <= its(7, 7); i++) {
            System.out.println(String.format("=> Cyclic%s:", i));
            List<MultivariatePolynomialZp64> ideal =
                    GroebnerBasisData.cyclic(i)
                            .stream()
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator())))
                            .map(p -> p.setOrdering(GREVLEX))
                            .collect(Collectors.toList());

            long start;

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = singular;
            System.out.println("   Singular  : " + nanosecondsToString(singular.nanoseconds));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + nanosecondsToString(f4));
            assertEquals(expected, actualF4);

            if (i <= 7) {
                List<MultivariatePolynomialZp64> actualE =
                        F4GB(ideal.stream().map(MultivariatePolynomialZp64::toBigPoly).collect(Collectors.toList()), GREVLEX)
                                .stream().map(MultivariatePolynomial::asOverZp64).collect(Collectors.toList());
                assertEquals(expected, actualE);
            }
        }
    }

    @Test
    @RequiresSingular
    public void test7() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("11*x^3+3*x^4*y*z^2+12*x^3*y^4+7*x^2*y^6*z+16*x^3*y*z^6+13*x^4*y*z^5+7*x^5*y*z^5+9*x^4*y^3*z^5+7*x^6*y^4*z^2+10*x^4*y^4*z^5", domain, GREVLEX, vars)),
                f2 = asOverZp64(parse("5*y+14*x*y*z^6+5*x^2*y^2*z^4+16*x*y^4*z^3+9*x^5*y^3*z+7*x^2*y^5*z^3+15*x^5*y^6+16*x^5*y*z^6+9*x^2*y^5*z^6+10*x^5*y^5*z^3", domain, GREVLEX, vars)),
                f3 = asOverZp64(parse("8*y^5+14*y^3*z^3+2*x*y^2*z^4+13*x*y^3*z^3+2*x^2*y^4*z^2+16*x^3*y^3*z^4+8*x^3*y^3*z^5+11*x^6*y^4*z^2+8*x^5*y^5*z^3", domain, GREVLEX, vars));


        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2, f3);
        setMonomialOrder(ideal, GREVLEX);

        for (int i = 0; i < its(2, 2); ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX);
            System.out.println("Singular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> f4 = GroebnerBases.F4GB(ideal, GREVLEX);
            System.out.println("F4: " + nanosecondsToString(System.nanoTime() - start));
            assertEquals(expected, f4);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> buchberger = BuchbergerGB(ideal, GREVLEX);
            System.out.println("Buchberger: " + nanosecondsToString(System.nanoTime() - start));
            assertEquals(expected, buchberger);
            System.out.println();
        }

//        Singular: 955ms
//        F4: 3s
//        Buchberger: 4s
//
//        Singular: 633ms
//        F4: 2s
//        Buchberger: 4s
//
//        Singular: 638ms
    }

    @Test
    @RequiresSingular
    public void test8() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("5*x*y^3+10*x^4*z^4+4*x^4*y*z^5+x^3*y^5*z^2+9*x^4*y^5*z^3", domain, GREVLEX, vars)),
                f2 = asOverZp64(parse("9*x*y^2*z^4+12*x^4*y^2*z+7*x^4*y^3*z^2+10*x^2*y^5*z^3+10*x^4*y^3*z^4", domain, GREVLEX, vars)),
                f3 = asOverZp64(parse("6*y+7*x^2*z^2+5*x^5*y+2*x^2*y^2*z^3+12*x^2*y*z^5", domain, GREVLEX, vars));

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2, f3);
        setMonomialOrder(ideal, GREVLEX);

        SingularResult<MonomialZp64, MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX);
        List<MultivariatePolynomialZp64> buchberger = BuchbergerGB(ideal, GREVLEX);
        assertEquals(expected, buchberger);
    }

    @Test
    @RequiresSingular
    public void test9() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("11*x^3*z^5+13*x^5*y^5*z^5", domain, GREVLEX, vars)),
                f2 = asOverZp64(parse("14*x^2*y^3*z^3+10*x^4*y^3*z+12*x^4*y*z^4", domain, GREVLEX, vars)),
                f3 = asOverZp64(parse("9*x*z^3+8*y^5*z+14*x*y^3*z^4", domain, GREVLEX, vars));

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2, f3);
        setMonomialOrder(ideal, GREVLEX);

        SingularResult<MonomialZp64, MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX);
        List<MultivariatePolynomialZp64> buchberger = BuchbergerGB(ideal, GREVLEX);
        assertEquals(expected, buchberger);
    }

    @Test
    @RequiresSingular
    public void test10_random() throws Exception {
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(100, 300);
        int nEls = 3;
        for (int i = 0; i < nIterations; i++) {
            List<MultivariatePolynomialZp64> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++)
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(3, 4, 3, domain, GREVLEX, rnd));
            List<MultivariatePolynomialZp64> actual = F4GB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX);
            if (!actual.equals(expected)) {
                System.out.println(ideal);
                System.out.println(actual.size());
            }
            assertEquals(expected, actual);
        }
    }

    @Test
    public void test11() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^17*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x^5*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), GREVLEX, vars);

        for (int i = 0; i < its(1, 1); i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> buch = GroebnerBases.BuchbergerGB(ideal, GREVLEX);
            System.out.println("Buchberger: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> f4 = GroebnerBases.F4GB(ideal, GREVLEX);
            System.out.println("F4: " + nanosecondsToString(System.nanoTime() - start));

            Assert.assertEquals(buch, f4);
        }
    }

    @Test
    public void test12() throws Exception {
        String[] vars = {"t", "x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2*t + x*y*z - y*z^4 + x", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^17*y^4 - z^4*y^2 + x*z*y^4 - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x^15*y^4*z - x*y^2 + x*z + t^3", Rings.Zp(17), GREVLEX, vars);
        List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));

        for (int i = 0; i < its(1, 1); i++) {
            long start;
            Comparator<SyzygyPair> strategy = GroebnerBases.normalSelectionStrategy(GREVLEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> norm = BuchbergerGB(ideal, GREVLEX, strategy);
            System.out.println("Normal strategy: " + nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBases.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> sugar = BuchbergerGB(ideal, GREVLEX, strategy);
            System.out.println("Sugar strategy:  " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(norm, sugar);
        }
    }

    @Test
    public void test13() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y*z - y*z^2 + x", Rings.Zp(17), vars),
                f2 = parse("x^3*y^4 - z*y^2 + x*z*y - 1", Rings.Zp(17), vars),
                f3 = parse("x*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), vars);

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));

        for (int i = 0; i < its(1, 1); i++) {
            long start;
            Comparator<SyzygyPair> strategy = GroebnerBases.normalSelectionStrategy(LEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> norm = BuchbergerGB(ideal, LEX, strategy);
            System.out.println("Normal strategy: " + nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBases.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> sugar = BuchbergerGB(ideal, LEX, strategy);
            System.out.println("Sugar strategy:  " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(norm, sugar);
        }
    }

    @Test
    @RequiresSingular
    public void test14_random() throws Exception {
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(1, 1);
        int nEls = 4;
        DescriptiveStatistics
                tF4 = new DescriptiveStatistics(),
                tSingular = new DescriptiveStatistics();
        for (int i = 0; i < nIterations; i++) {
            System.out.println();
            List<MultivariatePolynomialZp64> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++)
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(4, 4, 4, domain, GREVLEX, rnd));

            long start = System.nanoTime();
            List<MultivariatePolynomialZp64> actual = F4GB(ideal, GREVLEX);
            long time = System.nanoTime() - start;
            tF4.addValue(time);
            System.out.println("F4       : " + nanosecondsToString(time));

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            tSingular.addValue(singular.nanoseconds);
            System.out.println("Singular : " + nanosecondsToString(singular.nanoseconds));

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
    public void test15_random() throws Exception {
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(50, 50);
        int nEls = 3;
        DescriptiveStatistics
                tF4 = new DescriptiveStatistics(),
                tBuch = new DescriptiveStatistics(),
                tSingular = new DescriptiveStatistics();
        for (int i = 0; i < nIterations; i++) {
            List<MultivariatePolynomialZp64> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++)
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(4, 3, 4, domain, GREVLEX, rnd));

            //System.out.println();
            //System.out.println(ideal);
            long start = System.nanoTime();
            List<MultivariatePolynomialZp64> actual = F4GB(ideal, GREVLEX);
            long time = System.nanoTime() - start;
            tF4.addValue(time);
            //System.out.println("F4       : " + TimeUnits.nanosecondsToString(time));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> buch = BuchbergerGB(ideal, GREVLEX);
            time = System.nanoTime() - start;
            tBuch.addValue(time);
            //System.out.println("Buch     : " + TimeUnits.nanosecondsToString(time));

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            tSingular.addValue(singular.nanoseconds);
            //System.out.println("Singular : " + TimeUnits.nanosecondsToString(singular.nanoseconds));

            List<MultivariatePolynomialZp64> expected = singular;
            if (!actual.equals(expected)) {
                System.out.println("===>> ERROR");
                System.out.println(ideal);
                System.out.println(actual.size());
            }
            assertEquals(expected, actual);
            assertEquals(expected, buch);
        }

        System.out.println("F4         : " + TimeUnits.statisticsNanotime(tF4));
        System.out.println("Buchberger : " + TimeUnits.statisticsNanotime(tBuch));
        System.out.println("Singular   : " + TimeUnits.statisticsNanotime(tSingular));
    }

    @Test
    @RequiresSingular
    public void test15a_random() throws Exception {
        RandomGenerator rnd = getRandom();
        int nIterations = its(50, 100);
        int nEls = 3;
        DescriptiveStatistics
                tF4 = new DescriptiveStatistics(),
                tBuch = new DescriptiveStatistics(),
                tSingular = new DescriptiveStatistics();
        for (int i = 0; i < nIterations; i++) {
            List<MultivariatePolynomial<BigInteger>> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++) {
                MultivariatePolynomial<BigInteger> p = RandomMultivariatePolynomials.randomPolynomial(3, 3, 4, Z, GREVLEX, rnd);
                // bound coefficients
                p = p.setRing(Zp(100)).setRing(Z);
                ideal.add(p);
            }

            //System.out.println();
            //System.out.println(ideal);
            long start, time;
            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> f4 = F4GB(ideal, GREVLEX);
            time = System.nanoTime() - start;
            tF4.addValue(time);
            System.out.println("F4       : " + nanosecondsToString(time));

            //start = System.nanoTime();
            //List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(ideal, GREVLEX);
            //time = System.nanoTime() - start;
            //tBuch.addValue(time);
            //System.out.println("Buch     : " + TimeUnits.nanosecondsToString(time));

            SingularResult<?, MultivariatePolynomial<BigInteger>> singular = SingularGB(ideal, GREVLEX);
            tSingular.addValue(singular.nanoseconds);
            System.out.println("Singular : " + nanosecondsToString(singular.nanoseconds));
            System.out.println();

            List<MultivariatePolynomial<BigInteger>> expected = singular;
            //assertEquals(expected, buch);
            assertEquals(expected, f4);
        }

        System.out.println("F4         : " + TimeUnits.statisticsNanotime(tF4));
        //System.out.println("Buchberger : " + TimeUnits.statisticsNanotime(tBuch));
        System.out.println("Singular   : " + TimeUnits.statisticsNanotime(tSingular));
    }

    @Test
    public void test16() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("x1^2*x2^2*x3 + x5^4 - (1/2) + x2^3*x4 + x3^5 + x4", Q, GREVLEX, vars),
                b = parse("(1/123)*x1*x2*x3^2 + x2*x4*x1^2*x5 + x3*x2^3", Q, GREVLEX, vars),
                c = parse("x1^3*x2^3*x3^3*x4 - (1/17)*x2*x4^2*x1^4*x5 + x3^3*x2 - x4 - 1", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Arrays.asList(a, b, c);

        List<MultivariatePolynomial<Rational<BigInteger>>> buch = BuchbergerGB(gens, GREVLEX);
        assertTrue(buch.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral)));
        canonicalize(buch);

        List<MultivariatePolynomial<Rational<BigInteger>>> expected = BuchbergerGB(gens, GREVLEX,
                () -> new SyzygyTreeSet<>(new TreeSet<>(normalSelectionStrategy(GREVLEX))), null);

        assertEquals(expected, buch);

        List<MultivariatePolynomial<Rational<BigInteger>>> f4 = F4GB(gens, GREVLEX);
        assertEquals(expected, f4);
    }

    @Test
    public void test16a() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("x1^2*x2^2*x3 + x5^4 - (1/2) + x2^3*x4 + x3^5 + x4", Q, GREVLEX, vars),
                b = parse("(1/123)*x1*x2*x3^2 + x2*x4*x1^2*x5 + x3*x2^3", Q, GREVLEX, vars),
                c = parse("x1^3*x2^3*x3^3*x4 - (1/17)*x2*x4^2*x1^4*x5 + x3^3*x2 - x4 - 1", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Arrays.asList(c, b);

        List<MultivariatePolynomial<Rational<BigInteger>>> actual = BuchbergerGB(gens, GREVLEX);

        assertTrue(actual.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral)));
        canonicalize(actual);
        List<MultivariatePolynomial<Rational<BigInteger>>> expected = BuchbergerGB(gens, GREVLEX,
                () -> new SyzygyTreeSet<>(new TreeSet<>(normalSelectionStrategy(GREVLEX))), null);

        assertEquals(expected, actual);
        List<MultivariatePolynomial<Rational<BigInteger>>> f4 = F4GB(gens, GREVLEX);
        assertEquals(expected, f4);
    }

    @Test
    public void test17() throws Exception {
        for (int i = 0; i < 1; ++i) {
            System.out.println(i);
//            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            String[] vars = {"x1", "x2", "x3", "x4"};
            MultivariatePolynomial<BigInteger>
                    a = parse("6*x2*x4^3 + 11*x1*x3^3*x4 + 15*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
                    b = parse("11*x1^3 + 13*x3^2*x4^2 + x1^3*x2^3*x3 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
                    c = parse("7*x1^2*x2*x4 + 4*x1^3*x3 + 12*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
            List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

            List<MultivariatePolynomial<Rational<MultivariatePolynomial<BigInteger>>>> funcGens = gens.stream()
                    .map(p -> MultivariateConversions.split(p, 2, 3))
                    .map(p -> p.mapCoefficients(Frac(p.ring), cf -> new Rational<>(p.ring, cf)))
                    .collect(Collectors.toList());

            List<MultivariatePolynomial<Rational<MultivariatePolynomial<BigInteger>>>> funcGB = BuchbergerGB(funcGens, GREVLEX);
            assertTrue(funcGB.get(0).isConstant());
            System.out.println(nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @RequiresSingular
    public void test18() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2*y^2 + x*y - z", Q, GREVLEX, vars),
                f2 = parse("x^17*y - y^2 + x*z*y - 1", Q, GREVLEX, vars);

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(f1, f2);
        List<MultivariatePolynomial<Rational<BigInteger>>> sing = SingularGB(ideal, GREVLEX);
        List<MultivariatePolynomial<Rational<BigInteger>>> buch = canonicalize(GroebnerBases.BuchbergerGB(ideal, GREVLEX));
        Assert.assertEquals(sing, buch);
    }

    @Test
    @RequiresSingular
    public void test19() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^4*y^2 + x*y - z^2", Q, GREVLEX, vars),
                f2 = parse("x^17*y - x*y^2 + x*z*y - 1", Q, GREVLEX, vars),
                f3 = parse("x*y^17 - x*y^2 + x*z*y - 1", Q, GREVLEX, vars);

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(f1, f2, f3);
        List<MultivariatePolynomial<Rational<BigInteger>>> sing = SingularGB(ideal, GREVLEX);
        List<MultivariatePolynomial<Rational<BigInteger>>> buch = canonicalize(GroebnerBases.BuchbergerGB(ideal, GREVLEX));
        Assert.assertEquals(sing, buch);
    }

    @Test
    @Ignore // fixme: f4 meets expression swell, modular is too long
    public void test20() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("(-807152359/1175978805)*z+(1708357903/571090061)*x^2*y^2*z^2+(39119166838761599/323038390588954371)*x^3*y^2*z^3", Q, GREVLEX, vars),
                b = parse("(-960519798/1555504243)*x*z^3-(1278846706/1239147733)*x*y^3-(62586766/904196831)*x*y*z^3-(792306301/1609075855)*x^3*y^3*z^3", Q, GREVLEX, vars),
                c = parse("(9306287/4567935)*x^2-(1422841761/1340607578)*x*y*z^3-(115093936/778347949)*x^3*z^3-(44182447/32319755)*x^2*y^3*z^3", Q, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = toIntegral(Arrays.asList(a, b, c));

        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(gens, GREVLEX);
        // failed
        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> f4 = F4GB(gens, GREVLEX);
        // failed
        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX);
    }

    @Test(timeout = 10_000)
    public void test21() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("x^2 - y*z*x^2 + 2", Z, LEX, vars),
                b = parse("x^2*y*z^3 - y*z^2 + 2*x^5", Z, LEX, vars),
                c = parse("x*y*z^3 - y*z^12 + 2*x*y*z^5", Z, LEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        List<MultivariatePolynomial<BigInteger>> gb = GroebnerBasis(gens, LEX);
        assertEquals(3, gb.size());
    }

    @Test(timeout = 10_000)
    public void test22() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("x^12 - y*z*x^2 + 2", Z, LEX, vars),
                b = parse("x^2*y*z^3 - y*z^2 + 2*x^5", Z, LEX, vars),
                c = parse("x*y*z^3 - y*z^12 + 2*x*y*z^5", Z, LEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);
        assertEquals(UnivariatePolynomial.parse("1+3*x+6*x^2+10*x^3+15*x^4+21*x^5+27*x^6+33*x^7+11*x^8", Q),
                HilbertSeriesOfLeadingTermsSet(GroebnerBasisWithOptimizedGradedOrder(gens)).numerator);
    }

    @Test
    public void test23() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("-y^2 - 3*y + z^2 + 3", Z, LEX, vars),
                b = parse("x*z + y*z + z^2", Z, LEX, vars),
                c = parse("-3*x*y + 2*y*z + 6*z^2", Z, LEX, vars),
                d = parse("-2*y*z + z^2 + 2*z + 1", Z, LEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c, d);
        System.out.println(GroebnerBasis(mod(gens, 11777), LEX));
    }

    @Test
    public void testHilbertSeries1() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial
                a = parse("x^3*y^2*z^3", Q, GREVLEX, vars),
                b = parse("x^4*y*z^3", Q, GREVLEX, vars),
                c = parse("x^2*y^4*z^3", Q, GREVLEX, vars);
        List<MultivariatePolynomial> ideal = new ArrayList<>(Arrays.asList(a, b, c));
        HilbertSeries hps = HilbertSeriesOfLeadingTermsSet(ideal);
        assertEquals(2, hps.dimension());
        assertEquals(UnivariatePolynomial.parse("-5 + 6*x", Q), hps.hilbertPolynomial());
    }

    @Test
    public void testHilbertSeries2() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial
                a = parse("x*z", Q, GREVLEX, vars),
                b = parse("y*z", Q, GREVLEX, vars);
        List<MultivariatePolynomial> ideal = new ArrayList<>(Arrays.asList(a, b));
        HilbertSeries hps = HilbertSeriesOfLeadingTermsSet(ideal);
        assertEquals(2, hps.dimension());
        assertEquals(UnivariatePolynomial.parse("x + 2", Q), hps.hilbertPolynomial());
    }

    @Test
    public void testHilbertSeries3() throws Exception {
        String[] vars = {"x", "y", "z", "w", "u"};
        MultivariatePolynomial
                a = parse("x^13*y^2*z^3*w^7", Q, GREVLEX, vars),
                b = parse("x^4*y*z^3*w^6*u^8", Q, GREVLEX, vars),
                c = parse("x^2*y^4*z^3", Q, GREVLEX, vars),
                d = parse("x*y^14*z^3*u^6*w^6", Q, GREVLEX, vars);
        List<MultivariatePolynomial> ideal = new ArrayList<>(Arrays.asList(a, b, c, d));
        HilbertSeries hps = HilbertSeriesOfLeadingTermsSet(ideal);
        assertEquals(4, hps.dimension());
        assertEquals(UnivariatePolynomial.parse("2134+(-1702/3)*x+(57/2)*x^2+(5/6)*x^3", Q), hps.hilbertPolynomial());
    }

    @Test
    public void testHilbertSeries4() throws Exception {
        HilbertSeries hps = HilbertSeriesOfLeadingTermsSet(Collections.singletonList(MultivariatePolynomial.parse("1")));
        assertEquals(0, hps.dimension());
        assertEquals(0, hps.degree());
    }

    @Test
    public void testHilbertSeries5() throws Exception {
        HilbertSeries hps = HilbertSeriesOfLeadingTermsSet(Collections.singletonList(MultivariatePolynomial.parse("0*x*y*z")));
        assertEquals(3, hps.dimension());
        assertEquals(1, hps.degree());
    }

    @Test
    public void testHilbertGB5() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2 + 11*x1*x4 + 15*x2^3*x3*x4 + 13*x1*x4", Z, GREVLEX, vars),
                b = parse("11*x1 + 13*x3^2 + x1*x2^3*x3", Z, GREVLEX, vars),
                c = parse("12 + 7*x4 + 4*x1*x3 + 12*x1^2*x2^2", Z, GREVLEX, vars);

        List<MultivariatePolynomialZp64> gens = homogenize(mod(Arrays.asList(a, b, c), 17));

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            GBResult<MonomialZp64, MultivariatePolynomialZp64> expected = BuchbergerGB(gens, LEX);
            System.out.println("Buchberger : " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            GBResult<MonomialZp64, MultivariatePolynomialZp64> grevlex = F4GB(gens, GREVLEX);
            HilbertSeries hps = HilbertSeries(leadTermsSet(grevlex));
            GBResult<MonomialZp64, MultivariatePolynomialZp64> actual = HilbertGB(gens, LEX, hps);
            System.out.println("Hilbert    : " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(expected, actual);
            System.out.println();
        }
    }

    @Test
    public void testHilbertGB6() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2 + 11*x1*x4 + 15*x2^3*x3*x4 + 13*x1*x4", Z, GREVLEX, vars),
                b = parse("11*x1 + 13*x3^2 + x1*x2^3*x3", Z, GREVLEX, vars),
                c = parse("12 + 7*x4 + 4*x1*x3 + 12*x1^2*x2^2", Z, GREVLEX, vars);

        List<MultivariatePolynomialZp64> gens = mod(Arrays.asList(a, b, c), 17);

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> expected = SingularGB(gens, LEX);
            System.out.println("Singular   : " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            GBResult<MonomialZp64, MultivariatePolynomialZp64> actual = HilbertGB(gens, LEX);
            System.out.println("Hilbert    : " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(expected, actual);
            System.out.println();
        }
    }

    @Test
    public void testHilbertGB7() throws Exception {
        List<MultivariatePolynomialZp64> ideal =
                mod(GroebnerBasisData.katsura(8).stream().map(p -> p.mapCoefficients(Z, Rational::numerator)).collect(Collectors.toList()), 17);

        MultivariatePolynomialZp64 factory = ideal.get(0);

        ideal.get(1).add(factory.parsePoly("x1^5"));
//        ideal.get(2).add(factory.parsePoly("x2^5"));
//        ideal.get(3).add(factory.parsePoly("x3^5"));

        int[] degs = ideal.stream().map(AMultivariatePolynomial::degrees).reduce(new int[factory.nVariables], ArraysUtil::sum);
        System.out.println(Arrays.toString(degs));
        int[] perm = ArraysUtil.sequence(0, degs.length);
        ArraysUtil.quickSort(degs, perm);
//        ArraysUtil.reverse(perm);

        List<MultivariatePolynomialZp64> renamed = ideal.stream().map(p -> AMultivariatePolynomial.renameVariables(p, perm)).collect(Collectors.toList());
        System.out.println(Arrays.toString(renamed.stream().map(AMultivariatePolynomial::degrees).reduce(new int[factory.nVariables], ArraysUtil::sum)));

        long start;
        for (int i = 0; i < 20; ++i) {
            start = System.nanoTime();
            GBResult<MonomialZp64, MultivariatePolynomialZp64> g2 = F4GB(renamed, GREVLEX);
            System.out.println("Renamed: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            GBResult<MonomialZp64, MultivariatePolynomialZp64> g1 = F4GB(ideal, GREVLEX);
            System.out.println("Default: " + nanosecondsToString(System.nanoTime() - start));


            System.out.println();
        }
    }

    @Test
    public void testModularGB1() throws Exception {
        System.out.println(parse("0").size());
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<BigInteger>
                b = parse("x1*x2*x3^2 + 123*x2*x4*x1^2*x5 + 123*x3*x2^3", Z, GREVLEX, vars),
                c = parse("17*x1^3*x2^3*x3^3*x4 - x2*x4^2*x1^4*x5 + 17*x3^3*x2 - 17*x4 - 17", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(c, b);

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(gens, GREVLEX);
            System.out.println("Buchberger: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(buch, mod);
        }
    }

    @Test
    public void testModularGB2() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z", Z, GREVLEX, vars),
                b = parse("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4", Z, GREVLEX, vars),
                c = parse("8*x^3 + 12*y^3 + x*z^2 + 3", Z, GREVLEX, vars),
                d = parse("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c, d);


        List<MultivariatePolynomial<BigInteger>> expected = Arrays.asList(
                parse("x", Z, GREVLEX, vars),
                parse("z^2", Z, GREVLEX, vars),
                parse("1 + 4*y^3", Z, GREVLEX, vars)
        );

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));
            assertEquals(expected, mod);

            // Buchberger is too long here

//            start = System.nanoTime();
//            List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(gens, GREVLEX);
//            System.out.println("Buchberger: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
//
//            assertEquals(buch, mod);
        }
    }

    @Test
    public void testModularGB3() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2*x4^3 + 11*x1*x3^3*x4 + 15*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
                b = parse("11*x1^3 + 13*x3^2*x4^2 + x1^3*x2^3*x3 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
                c = parse("7*x1^2*x2*x4 + 4*x1^3*x3 + 12*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(gens, GREVLEX);
            System.out.println("Buchberger: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(buch, mod);
        }
    }

    @Test
    public void testModularGB4() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z", Z, GREVLEX, vars),
                b = parse("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4", Z, GREVLEX, vars),
                c = parse("8*x^3 + 12*y^3 + x*z^2 + 3", Z, GREVLEX, vars),
                d = parse("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c, d);

        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(1, 2); ++i) {
            rnd.setSeed(123 + i);
            List<MultivariatePolynomial<BigInteger>> shuffled = shuffleGB(gens, rnd, 2, 3);

            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> mod = ModularGB(shuffled, GREVLEX, null, false);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomial<BigInteger>> sparse = ModularGB(shuffled, GREVLEX, null, true);
            System.out.println("Sparse: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(mod, sparse);
        }
    }

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    boolean isReducible(Poly poly, List<Poly> reducers) {
        return reducers.stream().anyMatch(r ->
                poly.terms.keySet().stream().anyMatch(k -> k.dvDivisibleBy(r.lt())));
    }

    @Test
    public void testModularGB5_random() throws Exception {
        RandomGenerator rnd = getRandom();

        long constraint = its(50_000, 30_000);
        boolean print = true;
        for (int nVariables = 3; nVariables <= 3/*4*/; ++nVariables) {
            DescriptiveStatistics
                    tModInitialGB = new DescriptiveStatistics(),
                    tModShuffledGB = new DescriptiveStatistics(),
                    tBuchInitialGB = new DescriptiveStatistics(),
                    tBuchShuffledGB = new DescriptiveStatistics(),
                    tF4InitialGB = new DescriptiveStatistics(),
                    tF4ShuffledGB = new DescriptiveStatistics();

            int nIdealSize = 3, degree = 3;
            for (int i = 0; i < its(5, 10); ++i) {
                List<MultivariatePolynomial<BigInteger>> ideal = new ArrayList<>();
                for (int j = 0; j < nIdealSize; j++) {
                    MultivariatePolynomial<BigInteger> p;
                    int c = 0;
                    do {
                        p = RandomMultivariatePolynomials.randomPolynomial(nVariables, degree, 4, Z, GREVLEX, rnd);
                        // bound coefficients
                        p = p.setRing(Zp(5)).setRing(Z);
                        ++c;
                    } while (c < 10 && (p.isConstant() || isReducible(p, ideal)));

                    ideal.add(p);
                    for (int k = 0; k < 3; ++k) {
                        removeRedundant(ideal);
                        ideal.replaceAll(t -> t.setRing(Zp(5)).setRing(Z));
                        canonicalize(ideal);
                    }
                }

                long start;

                if (print) System.out.println("Initial ideal         : " + ideal);
                start = System.nanoTime();
                List<MultivariatePolynomial<BigInteger>> buchInitialGB = timeConstrained(() -> BuchbergerGB(ideal, GREVLEX), constraint);
                long tBuchIni = System.nanoTime() - start;
                tBuchInitialGB.addValue(tBuchIni);
                if (print) System.out.println("Initial Buchberger    : " + nanosecondsToString(tBuchIni));

                start = System.nanoTime();
                List<MultivariatePolynomial<BigInteger>> f4InitialGB = timeConstrained(() -> F4GB(ideal, GREVLEX), constraint);
                long tF4Ini = System.nanoTime() - start;
                tF4InitialGB.addValue(tF4Ini);
                if (print) System.out.println("Initial F4            : " + nanosecondsToString(tF4Ini));

                start = System.nanoTime();
                List<MultivariatePolynomial<BigInteger>> modInitialGB = timeConstrained(() -> ModularGB(ideal, GREVLEX), constraint);
                long tModIni = System.nanoTime() - start;
                tModInitialGB.addValue(tModIni);
                if (print) System.out.println("Initial modular       : " + nanosecondsToString(tModIni));

                if (buchInitialGB != null && modInitialGB != null)
                    assertEquals(buchInitialGB, modInitialGB);
                if (f4InitialGB != null && modInitialGB != null)
                    assertEquals(f4InitialGB, modInitialGB);

                long seed = rnd.nextLong();
                rnd.setSeed(seed);
                List<MultivariatePolynomial<BigInteger>> shuffled = shuffleGB(ideal, rnd, 2, 2);

                if (print) System.out.println("Shuffled seed         : " + seed);
                //if (print) System.out.println("Shuffled ideal      : " + shuffled);
                start = System.nanoTime();
                GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> buchShuffledGB = timeConstrained(() -> BuchbergerGB(shuffled, GREVLEX), constraint);
                long tBuchShu = System.nanoTime() - start;
                tBuchShuffledGB.addValue(tBuchShu);
                if (print) System.out.println("Shuffled Buchberger   : " + nanosecondsToString(tBuchShu));

                start = System.nanoTime();
                GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> f4ShuffledGB = timeConstrained(() -> F4GB(shuffled, GREVLEX), constraint);
                long tF4Shu = System.nanoTime() - start;
                tF4ShuffledGB.addValue(tF4Shu);
                if (print) System.out.println("Shuffled F4           : " + nanosecondsToString(tF4Shu));

                start = System.nanoTime();
                List<MultivariatePolynomial<BigInteger>> modShuffledGB = timeConstrained(() -> ModularGB(shuffled, GREVLEX), constraint);
                long tModShu = System.nanoTime() - start;
                tModShuffledGB.addValue(tModShu);
                if (print) System.out.println("Shuffled modular    : " + nanosecondsToString(tModShu));

                if (buchShuffledGB != null && modShuffledGB != null)
                    assertEquals(buchShuffledGB, modShuffledGB);
                if (f4ShuffledGB != null && modShuffledGB != null)
                    assertEquals(f4ShuffledGB, modShuffledGB);

                if (print && buchShuffledGB != null) {
                    System.out.println("Buchberger #processed : " + buchShuffledGB.nProcessedPolynomials);
                    System.out.println("Buchberger #redundant : " + buchShuffledGB.nZeroReductions);
                    System.out.println("Buchberger GB size    : " + buchShuffledGB.size());
                }
                if (print && f4ShuffledGB != null) {
                    System.out.println("F4 #processed         : " + f4ShuffledGB.nProcessedPolynomials);
                    System.out.println("F4 #redundant         : " + f4ShuffledGB.nZeroReductions);
                    System.out.println("F4 GB size            : " + f4ShuffledGB.size());
                }
                if (print) System.out.println("\n\n");
            }

            System.out.println(" =========================== " + MultivariateRing(nVariables, Z) + " ================================= ");
            System.out.println("Ring                : " + MultivariateRing(nVariables, Z));
            System.out.println("Initial Buchberger  : " + statisticsNanotime(tBuchInitialGB));
            System.out.println("Initial F4          : " + statisticsNanotime(tF4InitialGB));
            System.out.println("Initial modular     : " + statisticsNanotime(tModInitialGB));
            System.out.println("Shuffled Buchberger : " + statisticsNanotime(tBuchShuffledGB));
            System.out.println("Shuffled F4         : " + statisticsNanotime(tF4ShuffledGB));
            System.out.println("Shuffled modular    : " + statisticsNanotime(tModShuffledGB));
            System.out.println("\n\n");
        }
    }

    @Test
    public void testModularGB6() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("x*y^2*z^3", Z, GREVLEX, vars),
                b = parse("1+x^3*y*z^2", Z, GREVLEX, vars),
                c = parse("x*z+3*x^2*z+y^3*z", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(gens, GREVLEX);
        List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX, null, false);
        assertEquals(buch, mod);
    }

    @Test
    public void testModularGB7() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                g1 = parse("2*y*z+y^2+16*y^2*z^5+16*y^3*z^4+4*y^4*z^3+12*y^2*z^6+12*y^3*z^5+3*y^4*z^4+8*y^2*z^7-24*y^3*z^6-46*y^4*z^5-24*y^5*z^4+8*x^2*y^3*z^4-4*y^6*z^3+8*x^2*y^4*z^3+2*x^2*y^5*z^2-24*y^3*z^7-36*y^4*z^6-18*y^5*z^5+16*x^2*y^3*z^5-3*y^6*z^4+16*x^2*y^4*z^4+4*x^2*y^5*z^3-16*y^3*z^8-24*y^4*z^7-12*y^5*z^6+4*x^2*y^3*z^6-2*y^6*z^5-12*x^2*y^4*z^5-23*x^2*y^5*z^4-12*x^2*y^6*z^3-2*x^2*y^7*z^2+64*y^3*z^9+96*y^4*z^8+48*y^5*z^7+8*y^6*z^6-32*x^2*y^4*z^6-48*x^2*y^5*z^5-24*x^2*y^6*z^4-4*x^2*y^7*z^3+48*y^3*z^10+72*y^4*z^9+36*y^5*z^8+6*y^6*z^7-8*x^2*y^4*z^7-12*x^2*y^5*z^6-6*x^2*y^6*z^5-x^2*y^7*z^4+32*y^3*z^11+48*y^4*z^10+24*y^5*z^9+4*y^6*z^8+64*x^2*y^4*z^8+96*x^2*y^5*z^7+48*x^2*y^6*z^6+8*x^2*y^7*z^5+88*x^2*y^4*z^9+132*x^2*y^5*z^8+32*x^3*y^4*z^8+66*x^2*y^6*z^7+48*x^3*y^5*z^7+11*x^2*y^7*z^6+24*x^3*y^6*z^6+4*x^3*y^7*z^5+32*x^2*y^4*z^10+48*x^2*y^5*z^9+24*x^3*y^4*z^9+24*x^2*y^6*z^8+36*x^3*y^5*z^8+4*x^2*y^7*z^7+18*x^3*y^6*z^7+16*x^4*y^5*z^7+3*x^3*y^7*z^6+24*x^4*y^6*z^6+12*x^4*y^7*z^5+2*x^4*y^8*z^4+16*x^3*y^4*z^10+24*x^3*y^5*z^9+12*x^3*y^6*z^8+32*x^4*y^5*z^8+2*x^3*y^7*z^7+48*x^4*y^6*z^7+16*x^5*y^5*z^7+24*x^4*y^7*z^6+24*x^5*y^6*z^6+4*x^4*y^8*z^5+12*x^5*y^7*z^5+2*x^5*y^8*z^4+8*x^4*y^5*z^9+12*x^4*y^6*z^8+32*x^5*y^5*z^8+6*x^4*y^7*z^7+48*x^5*y^6*z^7+x^4*y^8*z^6+24*x^5*y^7*z^6+4*x^5*y^8*z^5+8*x^5*y^5*z^9+12*x^5*y^6*z^8+6*x^5*y^7*z^7+x^5*y^8*z^6", Z, GREVLEX, vars),
                g2 = parse("2*y*z+y^2+6*z^3+3*z^4+2*z^5+3*x^2*y*z^2+4*x^2*y*z^3+x^3*y*z^2+x^2*y*z^4-64*y^2*z^11-64*y^3*z^10-16*y^4*z^9-48*y^2*z^12-48*y^3*z^11-12*y^4*z^10-32*y^2*z^13-32*y^3*z^12-8*y^4*z^11-96*x^2*y^3*z^10-96*x^2*y^4*z^9-24*x^2*y^5*z^8-112*x^2*y^3*z^11-112*x^2*y^4*z^10-64*x^3*y^3*z^10-28*x^2*y^5*z^9-64*x^3*y^4*z^9-16*x^3*y^5*z^8-48*x^2*y^3*z^12-48*x^2*y^4*z^11-48*x^3*y^3*z^11-12*x^2*y^5*z^10-48*x^3*y^4*z^10-12*x^3*y^5*z^9-48*x^4*y^4*z^9-48*x^4*y^5*z^8-12*x^4*y^6*z^7-32*x^3*y^3*z^12-32*x^3*y^4*z^11-8*x^3*y^5*z^10-76*x^4*y^4*z^10-76*x^4*y^5*z^9-64*x^5*y^4*z^9-19*x^4*y^6*z^8-64*x^5*y^5*z^8-16*x^5*y^6*z^7-24*x^4*y^4*z^11-24*x^4*y^5*z^10-88*x^5*y^4*z^10-6*x^4*y^6*z^9-88*x^5*y^5*z^9-16*x^6*y^4*z^9-22*x^5*y^6*z^8-24*x^6*y^5*z^8-12*x^6*y^6*z^7-2*x^6*y^7*z^6-32*x^5*y^4*z^11-32*x^5*y^5*z^10-12*x^6*y^4*z^10-8*x^5*y^6*z^9-28*x^6*y^5*z^9-19*x^6*y^6*z^8-16*x^7*y^5*z^8-4*x^6*y^7*z^7-16*x^7*y^6*z^7-4*x^7*y^7*z^6-8*x^6*y^4*z^11-12*x^6*y^5*z^10-6*x^6*y^6*z^9-32*x^7*y^5*z^9-x^6*y^7*z^8-32*x^7*y^6*z^8-8*x^8*y^5*z^8-8*x^7*y^7*z^7-8*x^8*y^6*z^7-2*x^8*y^7*z^6-8*x^7*y^5*z^10-8*x^7*y^6*z^9-16*x^8*y^5*z^9-2*x^7*y^7*z^8-16*x^8*y^6*z^8-4*x^8*y^7*z^7-4*x^8*y^5*z^10-4*x^8*y^6*z^9-x^8*y^7*z^8", Z, GREVLEX, vars),
                g3 = parse("4*z^3+3*z^4+2*z^5+2*x^2*y*z^2+4*x^2*y*z^3+16*y^2*z^5+16*y^3*z^4+x^2*y*z^4+4*y^4*z^3+12*y^2*z^6+12*y^3*z^5+3*y^4*z^4+8*y^2*z^7+8*y^3*z^6+2*y^4*z^5+8*x^2*y^3*z^4+8*x^2*y^4*z^3+2*x^2*y^5*z^2+64*y^2*z^8+64*y^3*z^7+16*y^4*z^6+16*x^2*y^3*z^5+16*x^2*y^4*z^4+4*x^2*y^5*z^3+96*y^2*z^9+96*y^3*z^8+24*y^4*z^7+4*x^2*y^3*z^6+4*x^2*y^4*z^5+x^2*y^5*z^4+100*y^2*z^10+100*y^3*z^9+25*y^4*z^8+64*x^2*y^3*z^7+64*x^2*y^4*z^6+16*x^2*y^5*z^5+48*y^2*z^11+48*y^3*z^10+12*y^4*z^9+176*x^2*y^3*z^8+176*x^2*y^4*z^7+44*x^2*y^5*z^6+16*y^2*z^12+16*y^3*z^11+4*y^4*z^10+160*x^2*y^3*z^9+160*x^2*y^4*z^8+40*x^2*y^5*z^7+16*x^4*y^4*z^6+16*x^4*y^5*z^5+4*x^4*y^6*z^4+88*x^2*y^3*z^10+88*x^2*y^4*z^9+22*x^2*y^5*z^8+64*x^4*y^4*z^7+64*x^4*y^5*z^6+16*x^4*y^6*z^5+16*x^2*y^3*z^11+16*x^2*y^4*z^10+4*x^2*y^5*z^9+80*x^4*y^4*z^8+80*x^4*y^5*z^7+20*x^4*y^6*z^6+32*x^4*y^4*z^9+32*x^4*y^5*z^8+8*x^4*y^6*z^7+4*x^4*y^4*z^10+4*x^4*y^5*z^9+x^4*y^6*z^8", Z, GREVLEX, vars),
                g4 = parse("2*y*z+y^2-16*z^6-24*z^7-25*z^8-16*x^2*y*z^5-12*z^9-44*x^2*y*z^6-4*z^10-40*x^2*y*z^7-4*x^4*y^2*z^4-22*x^2*y*z^8-16*x^4*y^2*z^5+64*y^3*z^9-4*x^2*y*z^9+96*y^4*z^8+48*y^5*z^7+8*y^6*z^6-20*x^4*y^2*z^6-128*y^2*z^11-80*y^3*z^10+40*y^4*z^9+36*y^5*z^8+6*y^6*z^7-8*x^4*y^2*z^7-192*y^2*z^12-160*y^3*z^11+24*y^5*z^9+4*y^6*z^8+64*x^2*y^4*z^8-x^4*y^2*z^8+96*x^2*y^5*z^7+48*x^2*y^6*z^6+8*x^2*y^7*z^5-200*y^2*z^13-200*y^3*z^12-50*y^4*z^11-192*x^2*y^3*z^10-104*x^2*y^4*z^9+84*x^2*y^5*z^8+32*x^3*y^4*z^8+66*x^2*y^6*z^7+48*x^3*y^5*z^7+11*x^2*y^7*z^6+24*x^3*y^6*z^6+4*x^3*y^7*z^5-96*y^2*z^14-96*y^3*z^13-24*y^4*z^12-448*x^2*y^3*z^11-416*x^2*y^4*z^10-64*x^3*y^3*z^10-64*x^2*y^5*z^9-40*x^3*y^4*z^9+24*x^2*y^6*z^8+20*x^3*y^5*z^8+4*x^2*y^7*z^7+18*x^3*y^6*z^7+16*x^4*y^5*z^7+3*x^3*y^7*z^6+24*x^4*y^6*z^6+12*x^4*y^7*z^5+2*x^4*y^8*z^4-32*y^2*z^15-32*y^3*z^14-8*y^4*z^13-420*x^2*y^3*z^12-420*x^2*y^4*z^11-96*x^3*y^3*z^11-105*x^2*y^5*z^10-80*x^3*y^4*z^10-96*x^4*y^4*z^9+12*x^3*y^6*z^8-64*x^4*y^5*z^8+2*x^3*y^7*z^7+24*x^4*y^6*z^7+16*x^5*y^5*z^7+24*x^4*y^7*z^6+24*x^5*y^6*z^6+4*x^4*y^8*z^5+12*x^5*y^7*z^5+2*x^5*y^8*z^4-224*x^2*y^3*z^13-224*x^2*y^4*z^12-100*x^3*y^3*z^12-56*x^2*y^5*z^11-100*x^3*y^4*z^11-25*x^3*y^5*z^10-304*x^4*y^4*z^10-296*x^4*y^5*z^9-64*x^5*y^4*z^9-64*x^4*y^6*z^8-32*x^5*y^5*z^8+6*x^4*y^7*z^7+32*x^5*y^6*z^7+x^4*y^8*z^6+24*x^5*y^7*z^6+4*x^5*y^8*z^5-48*x^2*y^3*z^14-48*x^2*y^4*z^13-48*x^3*y^3*z^13-12*x^2*y^5*z^12-48*x^3*y^4*z^12-12*x^3*y^5*z^11-320*x^4*y^4*z^11-320*x^4*y^5*z^10-176*x^5*y^4*z^10-80*x^4*y^6*z^9-168*x^5*y^5*z^9-32*x^5*y^6*z^8-16*x^6*y^5*z^8+6*x^5*y^7*z^7-16*x^6*y^6*z^7+x^5*y^8*z^6-4*x^6*y^7*z^6-16*x^3*y^3*z^14-16*x^3*y^4*z^13-4*x^3*y^5*z^12-152*x^4*y^4*z^12-152*x^4*y^5*z^11-160*x^5*y^4*z^11-38*x^4*y^6*z^10-160*x^5*y^5*z^10-40*x^5*y^6*z^9-64*x^6*y^5*z^9-64*x^6*y^6*z^8-16*x^7*y^5*z^8-16*x^6*y^7*z^7-16*x^7*y^6*z^7-4*x^7*y^7*z^6-24*x^4*y^4*z^13-24*x^4*y^5*z^12-88*x^5*y^4*z^12-6*x^4*y^6*z^11-88*x^5*y^5*z^11-22*x^5*y^6*z^10-80*x^6*y^5*z^10-80*x^6*y^6*z^9-64*x^7*y^5*z^9-20*x^6*y^7*z^8-64*x^7*y^6*z^8-16*x^7*y^7*z^7-16*x^5*y^4*z^13-16*x^5*y^5*z^12-4*x^5*y^6*z^11-32*x^6*y^5*z^11-32*x^6*y^6*z^10-80*x^7*y^5*z^10-8*x^6*y^7*z^9-80*x^7*y^6*z^9-20*x^7*y^7*z^8-4*x^6*y^5*z^12-4*x^6*y^6*z^11-32*x^7*y^5*z^11-x^6*y^7*z^10-32*x^7*y^6*z^10-8*x^7*y^7*z^9-4*x^7*y^5*z^12-4*x^7*y^6*z^11-x^7*y^7*z^10", Z, GREVLEX, vars),
                g5 = parse("2*z^3+x^2*y*z^2-8*z^6+x^3*y*z^2-6*z^7-4*z^8-8*x^2*y*z^5-11*x^2*y*z^6-4*x^3*y*z^5-4*x^2*y*z^7-3*x^3*y*z^6-2*x^4*y^2*z^4+64*y*z^10+32*y^2*z^9-2*x^3*y*z^7-4*x^4*y^2*z^5-2*x^5*y^2*z^4+96*y*z^11+48*y^2*z^10-x^4*y^2*z^6-4*x^5*y^2*z^5+100*y*z^12+50*y^2*z^11+96*x^2*y^2*z^9+48*x^2*y^3*z^8-x^5*y^2*z^6-80*y*z^13-40*y^2*z^12+224*x^2*y^2*z^10+112*x^2*y^3*z^9+32*x^3*y^2*z^9+16*x^3*y^3*z^8-176*y*z^14-88*y^2*z^13+210*x^2*y^2*z^11+105*x^2*y^3*z^10+48*x^3*y^2*z^10+24*x^3*y^3*z^9+48*x^4*y^3*z^8+24*x^4*y^4*z^7-200*y*z^15-100*y^2*z^14-144*x^2*y^2*z^12-72*x^2*y^3*z^11+50*x^3*y^2*z^11+25*x^3*y^3*z^10+152*x^4*y^3*z^9+76*x^4*y^4*z^8+32*x^5*y^3*z^8+16*x^5*y^4*z^7-96*y*z^16-48*y^2*z^15-520*x^2*y^2*z^13-260*x^2*y^3*z^12-104*x^3*y^2*z^12-52*x^3*y^3*z^11+160*x^4*y^3*z^10+80*x^4*y^4*z^9+88*x^5*y^3*z^9+44*x^5*y^4*z^8+8*x^6*y^4*z^7+4*x^6*y^5*z^6-32*y*z^17-16*y^2*z^16-520*x^2*y^2*z^14-260*x^2*y^3*z^13-184*x^3*y^2*z^13-92*x^3*y^3*z^12-116*x^4*y^3*z^11-58*x^4*y^4*z^10+80*x^5*y^3*z^10+40*x^5*y^4*z^9+32*x^6*y^4*z^8+16*x^6*y^5*z^7+8*x^7*y^4*z^7+4*x^7*y^5*z^6-272*x^2*y^2*z^15-136*x^2*y^3*z^14-200*x^3*y^2*z^14-100*x^3*y^3*z^13-516*x^4*y^3*z^12-258*x^4*y^4*z^11-148*x^5*y^3*z^11-74*x^5*y^4*z^10+40*x^6*y^4*z^9+20*x^6*y^5*z^8+32*x^7*y^4*z^8+16*x^7*y^5*z^7-64*x^2*y^2*z^16-32*x^2*y^3*z^15-96*x^3*y^2*z^15-48*x^3*y^3*z^14-530*x^4*y^3*z^13-265*x^4*y^4*z^12-440*x^5*y^3*z^12-220*x^5*y^4*z^11-32*x^6*y^3*z^11-64*x^6*y^4*z^10-24*x^6*y^5*z^9+40*x^7*y^4*z^9+20*x^7*y^5*z^8-32*x^3*y^2*z^16-16*x^3*y^3*z^15-264*x^4*y^3*z^14-132*x^4*y^4*z^13-420*x^5*y^3*z^13-210*x^5*y^4*z^12-48*x^6*y^3*z^12-238*x^6*y^4*z^11-107*x^6*y^5*z^10-80*x^7*y^4*z^10-40*x^7*y^5*z^9-48*x^4*y^3*z^15-24*x^4*y^4*z^14-224*x^5*y^3*z^14-112*x^5*y^4*z^13-50*x^6*y^3*z^13-265*x^6*y^4*z^12-120*x^6*y^5*z^11-302*x^7*y^4*z^11-151*x^7*y^5*z^10-32*x^8*y^4*z^10-24*x^8*y^5*z^9-4*x^8*y^6*z^8-48*x^5*y^3*z^15-24*x^5*y^4*z^14-24*x^6*y^3*z^14-120*x^6*y^4*z^13-54*x^6*y^5*z^12-320*x^7*y^4*z^12-160*x^7*y^5*z^11-88*x^8*y^4*z^11-76*x^8*y^5*z^10-16*x^8*y^6*z^9-16*x^9*y^5*z^9-8*x^9*y^6*z^8-8*x^6*y^3*z^15-20*x^6*y^4*z^14-8*x^6*y^5*z^13-152*x^7*y^4*z^13-76*x^7*y^5*z^12-80*x^8*y^4*z^12-80*x^8*y^5*z^11-20*x^8*y^6*z^10-64*x^9*y^5*z^10-32*x^9*y^6*z^9-8*x^10*y^5*z^9-4*x^10*y^6*z^8-24*x^7*y^4*z^14-12*x^7*y^5*z^13-44*x^8*y^4*z^13-38*x^8*y^5*z^12-8*x^8*y^6*z^11-80*x^9*y^5*z^11-40*x^9*y^6*z^10-32*x^10*y^5*z^10-16*x^10*y^6*z^9-8*x^8*y^4*z^14-6*x^8*y^5*z^13-x^8*y^6*z^12-32*x^9*y^5*z^12-16*x^9*y^6*z^11-40*x^10*y^5*z^11-20*x^10*y^6*z^10-4*x^9*y^5*z^13-2*x^9*y^6*z^12-16*x^10*y^5*z^12-8*x^10*y^6*z^11-2*x^10*y^5*z^13-x^10*y^6*z^12", Z, GREVLEX, vars),
                g6 = parse("8*z^3+6*z^4+4*z^5+4*x^2*y*z^2+8*x^2*y*z^3+2*x^2*y*z^4-32*y*z^10-16*y^2*z^9-24*y*z^11-12*y^2*z^10-16*y*z^12-8*y^2*z^11-48*x^2*y^2*z^9-24*x^2*y^3*z^8-56*x^2*y^2*z^10-28*x^2*y^3*z^9-32*x^3*y^2*z^9-16*x^3*y^3*z^8-24*x^2*y^2*z^11-12*x^2*y^3*z^10-24*x^3*y^2*z^10-12*x^3*y^3*z^9-24*x^4*y^3*z^8-12*x^4*y^4*z^7-16*x^3*y^2*z^11-8*x^3*y^3*z^10-38*x^4*y^3*z^9-19*x^4*y^4*z^8-32*x^5*y^3*z^8-16*x^5*y^4*z^7-12*x^4*y^3*z^10-6*x^4*y^4*z^9-44*x^5*y^3*z^9-22*x^5*y^4*z^8-8*x^6*y^3*z^8-8*x^6*y^4*z^7-2*x^6*y^5*z^6-16*x^5*y^3*z^10-8*x^5*y^4*z^9-6*x^6*y^3*z^9-11*x^6*y^4*z^8-4*x^6*y^5*z^7-8*x^7*y^4*z^7-4*x^7*y^5*z^6-4*x^6*y^3*z^10-4*x^6*y^4*z^9-x^6*y^5*z^8-16*x^7*y^4*z^8-8*x^7*y^5*z^7-4*x^8*y^4*z^7-2*x^8*y^5*z^6-4*x^7*y^4*z^9-2*x^7*y^5*z^8-8*x^8*y^4*z^8-4*x^8*y^5*z^7-2*x^8*y^4*z^9-x^8*y^5*z^8", Z, GREVLEX, vars),
                g7 = parse("2*y*z+y^2+8*y^2*z^5+8*y^3*z^4+2*y^4*z^3-16*y*z^7-8*y^2*z^6-12*y*z^8-6*y^2*z^7+32*y^3*z^6+48*y^4*z^5+24*y^5*z^4+4*x^2*y^3*z^4+4*y^6*z^3+4*x^2*y^4*z^3+x^2*y^5*z^2-8*y*z^9-4*y^2*z^8-16*x^2*y^2*z^6-8*x^2*y^3*z^5+4*x^3*y^3*z^4+4*x^3*y^4*z^3+x^3*y^5*z^2-22*x^2*y^2*z^7-11*x^2*y^3*z^6-8*x^3*y^2*z^6+16*x^2*y^4*z^5-4*x^3*y^3*z^5+24*x^2*y^5*z^4+12*x^2*y^6*z^3+2*x^2*y^7*z^2-8*x^2*y^2*z^8-4*x^2*y^3*z^7-6*x^3*y^2*z^7-3*x^3*y^3*z^6+16*x^3*y^4*z^5-4*x^4*y^3*z^5+24*x^3*y^5*z^4-2*x^4*y^4*z^4+12*x^3*y^6*z^3+2*x^3*y^7*z^2-4*x^3*y^2*z^8-2*x^3*y^3*z^7-8*x^4*y^3*z^6-4*x^4*y^4*z^5-4*x^5*y^3*z^5-2*x^5*y^4*z^4-2*x^4*y^3*z^7-x^4*y^4*z^6-8*x^5*y^3*z^6-4*x^5*y^4*z^5-2*x^5*y^3*z^7-x^5*y^4*z^6", Z, GREVLEX, vars),
                g8 = parse("2*z^3+x^2*y*z^2+x^3*y*z^2-1024*y^2*z^17-1024*y^3*z^16-256*y^4*z^15-2304*y^2*z^18-2304*y^3*z^17-576*y^4*z^16-3264*y^2*z^19-3264*y^3*z^18-816*y^4*z^17-2560*x^2*y^3*z^16-2560*x^2*y^4*z^15-640*x^2*y^5*z^14-2736*y^2*z^20-2736*y^3*z^19-684*y^4*z^18-7680*x^2*y^3*z^17-7680*x^2*y^4*z^16-1024*x^3*y^3*z^16-1920*x^2*y^5*z^15-1024*x^3*y^4*z^15-256*x^3*y^5*z^14-1632*y^2*z^21-1632*y^3*z^20-408*y^4*z^19-11040*x^2*y^3*z^18-11040*x^2*y^4*z^17-2304*x^3*y^3*z^17-2760*x^2*y^5*z^16-2304*x^3*y^4*z^16-576*x^3*y^5*z^15-2560*x^4*y^4*z^15-2560*x^4*y^5*z^14-640*x^4*y^6*z^13-576*y^2*z^22-576*y^3*z^21-144*y^4*z^20-9840*x^2*y^3*z^19-9840*x^2*y^4*z^18-3264*x^3*y^3*z^18-2460*x^2*y^5*z^17-3264*x^3*y^4*z^17-816*x^3*y^5*z^16-9600*x^4*y^4*z^16-9600*x^4*y^5*z^15-2048*x^5*y^4*z^15-2400*x^4*y^6*z^14-2048*x^5*y^5*z^14-512*x^5*y^6*z^13-128*y^2*z^23-128*y^3*z^22-32*y^4*z^21-5520*x^2*y^3*z^20-5520*x^2*y^4*z^19-2736*x^3*y^3*z^19-1380*x^2*y^5*z^18-2736*x^3*y^4*z^18-684*x^3*y^5*z^17-15120*x^4*y^4*z^17-15120*x^4*y^5*z^16-6528*x^5*y^4*z^16-3780*x^4*y^6*z^15-6528*x^5*y^5*z^15-256*x^6*y^4*z^15-1632*x^5*y^6*z^14-1536*x^6*y^5*z^14-1344*x^6*y^6*z^13-320*x^6*y^7*z^12-1920*x^2*y^3*z^21-1920*x^2*y^4*z^20-1632*x^3*y^3*z^20-480*x^2*y^5*z^19-1632*x^3*y^4*z^19-408*x^3*y^5*z^18-13740*x^4*y^4*z^18-13740*x^4*y^5*z^17-9408*x^5*y^4*z^17-3435*x^4*y^6*z^16-9408*x^5*y^5*z^16-576*x^6*y^4*z^16-2352*x^5*y^6*z^15-6336*x^6*y^5*z^15-5904*x^6*y^6*z^14-1536*x^7*y^5*z^14-1440*x^6*y^7*z^13-1536*x^7*y^6*z^13-384*x^7*y^7*z^12-320*x^2*y^3*z^22-320*x^2*y^4*z^21-576*x^3*y^3*z^21-80*x^2*y^5*z^20-576*x^3*y^4*z^20-144*x^3*y^5*z^19-7560*x^4*y^4*z^19-7560*x^4*y^5*z^18-8472*x^5*y^4*z^18-1890*x^4*y^6*z^17-8472*x^5*y^5*z^17-816*x^6*y^4*z^17-2118*x^5*y^6*z^16-11016*x^6*y^5*z^16-10404*x^6*y^6*z^15-6336*x^7*y^5*z^15-2550*x^6*y^7*z^14-6336*x^7*y^6*z^14-384*x^8*y^5*z^14-1584*x^7*y^7*z^13-704*x^8*y^6*z^13-416*x^8*y^7*z^12-80*x^8*y^8*z^11-128*x^3*y^3*z^22-128*x^3*y^4*z^21-32*x^3*y^5*z^20-2400*x^4*y^4*z^20-2400*x^4*y^5*z^19-4704*x^5*y^4*z^19-600*x^4*y^6*z^18-4704*x^5*y^5*z^18-684*x^6*y^4*z^18-1176*x^5*y^6*z^17-10204*x^6*y^5*z^17-9691*x^6*y^6*z^16-10416*x^7*y^5*z^16-2380*x^6*y^7*z^15-10416*x^7*y^6*z^15-1344*x^8*y^5*z^15-2604*x^7*y^7*z^14-3024*x^8*y^6*z^14-2016*x^8*y^7*z^13-512*x^9*y^6*z^13-420*x^8*y^8*z^12-512*x^9*y^7*z^12-128*x^9*y^8*z^11-320*x^4*y^4*z^21-320*x^4*y^5*z^20-1632*x^5*y^4*z^20-80*x^4*y^6*z^19-1632*x^5*y^5*z^19-408*x^6*y^4*z^19-408*x^5*y^6*z^18-5508*x^6*y^5*z^18-5202*x^6*y^6*z^17-9504*x^7*y^5*z^17-1275*x^6*y^7*z^16-9504*x^7*y^6*z^16-1944*x^8*y^5*z^16-2376*x^7*y^7*z^15-5304*x^8*y^6*z^15-3846*x^8*y^7*z^14-2592*x^9*y^6*z^14-840*x^8*y^8*z^13-2592*x^9*y^7*z^13-192*x^10*y^6*z^13-648*x^9*y^8*z^12-224*x^10*y^7*z^12-80*x^10*y^8*z^11-8*x^10*y^9*z^10-256*x^5*y^4*z^21-256*x^5*y^5*z^20-144*x^6*y^4*z^20-64*x^5*y^6*z^19-1584*x^6*y^5*z^19-1476*x^6*y^6*z^18-5208*x^7*y^5*z^18-360*x^6*y^7*z^17-5208*x^7*y^6*z^17-1776*x^8*y^5*z^17-1302*x^7*y^7*z^16-5056*x^8*y^6*z^16-3724*x^8*y^7*z^15-4992*x^9*y^6*z^15-820*x^8*y^8*z^14-4992*x^9*y^7*z^14-912*x^10*y^6*z^14-1248*x^9*y^8*z^13-1104*x^10*y^7*z^13-420*x^10*y^8*z^12-64*x^11*y^7*z^12-48*x^10*y^9*z^11-64*x^11*y^8*z^11-16*x^11*y^9*z^10-32*x^6*y^4*z^21-192*x^6*y^5*z^20-168*x^6*y^6*z^19-1584*x^7*y^5*z^19-40*x^6*y^7*z^18-1584*x^7*y^6*z^18-972*x^8*y^5*z^18-396*x^7*y^7*z^17-2652*x^8*y^6*z^17-1923*x^8*y^7*z^16-4768*x^9*y^6*z^16-420*x^8*y^8*z^15-4768*x^9*y^7*z^15-1632*x^10*y^6*z^15-1192*x^9*y^8*z^14-2064*x^10*y^7*z^14-840*x^10*y^8*z^13-384*x^11*y^7*z^13-108*x^10*y^9*z^12-384*x^11*y^8*z^12-32*x^12*y^7*z^12-96*x^11*y^9*z^11-32*x^12*y^8*z^11-8*x^12*y^9*z^10-192*x^7*y^5*z^20-192*x^7*y^6*z^19-336*x^8*y^5*z^19-48*x^7*y^7*z^18-756*x^8*y^6*z^18-504*x^8*y^7*z^17-2496*x^9*y^6*z^17-105*x^8*y^8*z^16-2496*x^9*y^7*z^16-1488*x^10*y^6*z^16-624*x^9*y^8*z^15-1936*x^10*y^7*z^15-820*x^10*y^8*z^14-864*x^11*y^7*z^14-112*x^10*y^9*z^13-864*x^11*y^8*z^13-192*x^12*y^7*z^13-216*x^11*y^9*z^12-192*x^12*y^8*z^12-48*x^12*y^9*z^11-48*x^8*y^5*z^20-88*x^8*y^6*z^19-52*x^8*y^7*z^18-648*x^9*y^6*z^18-10*x^8*y^8*z^17-648*x^9*y^7*z^17-816*x^10*y^6*z^17-162*x^9*y^8*z^16-1032*x^10*y^7*z^16-420*x^10*y^8*z^15-896*x^11*y^7*z^15-54*x^10*y^9*z^14-896*x^11*y^8*z^14-432*x^12*y^7*z^14-224*x^11*y^9*z^13-432*x^12*y^8*z^13-108*x^12*y^9*z^12-64*x^9*y^6*z^19-64*x^9*y^7*z^18-228*x^10*y^6*z^18-16*x^9*y^8*z^17-276*x^10*y^7*z^17-105*x^10*y^8*z^16-432*x^11*y^7*z^16-12*x^10*y^9*z^15-432*x^11*y^8*z^15-448*x^12*y^7*z^15-108*x^11*y^9*z^14-448*x^12*y^8*z^14-112*x^12*y^9*z^13-24*x^10*y^6*z^19-28*x^10*y^7*z^18-10*x^10*y^8*z^17-96*x^11*y^7*z^17-x^10*y^9*z^16-96*x^11*y^8*z^16-216*x^12*y^7*z^16-24*x^11*y^9*z^15-216*x^12*y^8*z^15-54*x^12*y^9*z^14-8*x^11*y^7*z^18-8*x^11*y^8*z^17-48*x^12*y^7*z^17-2*x^11*y^9*z^16-48*x^12*y^8*z^16-12*x^12*y^9*z^15-4*x^12*y^7*z^18-4*x^12*y^8*z^17-x^12*y^9*z^16", Z, GREVLEX, vars),
                g9 = parse("4*z^3+3*z^4+2*z^5+2*x^2*y*z^2+8*z^6+4*x^2*y*z^3+6*z^7-8*y^2*z^5-8*y^3*z^4+x^2*y*z^4-2*y^4*z^3+4*z^8+8*y*z^7+4*y^2*z^6+8*x^2*y*z^5+11*x^2*y*z^6+4*x^3*y*z^5-4*x^2*y^3*z^4-4*x^2*y^4*z^3-x^2*y^5*z^2+32*y^2*z^8+32*y^3*z^7+4*x^2*y*z^7+8*y^4*z^6+8*x^2*y^2*z^6+3*x^3*y*z^6+4*x^2*y^3*z^5-4*x^3*y^3*z^4+2*x^4*y^2*z^4-4*x^3*y^4*z^3-x^3*y^5*z^2+24*y^2*z^9+24*y^3*z^8+6*y^4*z^7+2*x^3*y*z^7+8*x^3*y^2*z^6+4*x^3*y^3*z^5+4*x^4*y^2*z^5+2*x^5*y^2*z^4+16*y^2*z^10+16*y^3*z^9+4*y^4*z^8+32*x^2*y^3*z^7+32*x^2*y^4*z^6+x^4*y^2*z^6+8*x^2*y^5*z^5+2*x^4*y^3*z^5+4*x^5*y^2*z^5+x^4*y^4*z^4+44*x^2*y^3*z^8+44*x^2*y^4*z^7+16*x^3*y^3*z^7+11*x^2*y^5*z^6+16*x^3*y^4*z^6+x^5*y^2*z^6+4*x^3*y^5*z^5+4*x^5*y^3*z^5+2*x^5*y^4*z^4+16*x^2*y^3*z^9+16*x^2*y^4*z^8+12*x^3*y^3*z^8+4*x^2*y^5*z^7+12*x^3*y^4*z^7+3*x^3*y^5*z^6+8*x^4*y^4*z^6+8*x^4*y^5*z^5+2*x^6*y^3*z^5+2*x^4*y^6*z^4+x^6*y^4*z^4+8*x^3*y^3*z^9+8*x^3*y^4*z^8+2*x^3*y^5*z^7+16*x^4*y^4*z^7+16*x^4*y^5*z^6+8*x^5*y^4*z^6+4*x^4*y^6*z^5+8*x^5*y^5*z^5+2*x^5*y^6*z^4+4*x^4*y^4*z^8+4*x^4*y^5*z^7+16*x^5*y^4*z^7+x^4*y^6*z^6+16*x^5*y^5*z^6+4*x^5*y^6*z^5+4*x^5*y^4*z^8+4*x^5*y^5*z^7+x^5*y^6*z^6", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(g1, g2, g3, g4, g5, g6, g7, g8, g9);

        List<MultivariatePolynomial<BigInteger>> f4 = F4GB(gens, GREVLEX);
        // NOTE: not working with trySparse = true
        List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX, null, false);
        assertEquals(f4, mod);
    }

    @Test
    @RequiresSingular
    public void testModularGB8() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                g1 = parse("3*z^2+x*y^3*z^2", Z, GREVLEX, vars),
                g2 = parse("2*z^3+x^2*y*z+4*x^3*y^3*z", Z, GREVLEX, vars),
                g3 = parse("2*x^2*y^2*z^2+2*x^3*y^2*z+x^3*y^2*z^2", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(g1, g2, g3);

        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> f4 = F4GB(gens, GREVLEX);
        List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX, null, false);
        assertEquals(f4, mod);

        RandomGenerator rnd = getRandom();
        rnd.setSeed(-6534340102157296164L);

        List<MultivariatePolynomial<BigInteger>> shuffled = shuffleGB(gens, rnd, 2, 3);

        // NOTE: not working with trySparse = true
        mod = ModularGB(shuffled, GREVLEX, null, false);
        assertEquals(SingularGB(shuffled, GREVLEX), mod);

        // fixme both F4 and Buchberger takes too long due to intermediate expression swell
//        f4 = F4GB(shuffled, GREVLEX);
//        assertEquals(f4, mod);
    }

    @Test
    public void testModularGB9() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2*x4^3 + 11*x1*x3^3*x4 + 15*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
                b = parse("11*x1^3 + 13*x3^2*x4^2 + x1^3*x2^3*x3 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
                c = parse("7*x1^2*x2*x4 + 4*x1^3*x3 + 12*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> f4 = F4GB(gens, GREVLEX);
        List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX, null, false);
        assertEquals(f4, mod);
    }

    @Test
    @TimeConsuming
    public void testModularGB10() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("2*x3*x4+4*x1^3*x3^3+2*x1^3*x2^3+3*x2^2*x3^3*x4^2", Z, GREVLEX, vars),
                b = parse("x3^2*x4^2+x1*x2*x3*x4+4*x2^3*x4^3+2*x1^2*x2^2*x4^3", Z, GREVLEX, vars),
                c = parse("4*x1*x3*x4^2+4*x1^3*x3*x4+x1^2*x2^3+2*x1*x2*x3^2*x4^2+2*x1^2*x2^2*x3*x4+3*x1*x2^4*x4^3", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        GBResult<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> f4 = F4GB(gens, GREVLEX);
        List<MultivariatePolynomial<BigInteger>> mod = ModularGB(gens, GREVLEX, null, false);
        assertEquals(f4, mod);
    }

    @Test
    public void testSparseGB1() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z", Z, GREVLEX, vars),
                b = parse("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4", Z, GREVLEX, vars),
                c = parse("8*x^3 + 12*y^3 + x*z^2 + 3", Z, GREVLEX, vars),
                d = parse("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c, d);


        List<MultivariatePolynomial<BigInteger>> expected = new ArrayList<>(Arrays.asList(
                parse("x", Z, GREVLEX, vars),
                parse("z^2", Z, GREVLEX, vars),
                parse("1 + 4*y^3", Z, GREVLEX, vars)
        ));

        List<MultivariatePolynomial<BigInteger>> gb = solveGB(gens,
                expected.stream().map(p -> ((SortedSet<DegreeVector>) p.terms.keySet())).collect(Collectors.toList()), GREVLEX);
        assertEquals(expected, gb);
    }

    @Test
    @Ignore("too long")
    public void testSparseGB2() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2*x4^3 + 11*x1*x3^3*x4 + 15*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
                b = parse("11*x1^3 + 13*x3^2*x4^2 + x1^3*x2^3*x3 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
                c = parse("7*x1^2*x2*x4 + 4*x1^3*x3 + 12*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
//                a = parse("112312423412343253451*x1*x3^3*x4 + 11232143245*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
//                b = parse("1111*x1^3 + 13*x3^2*x4^2 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
//                c = parse("12345435135345345413457*x1^2*x2*x4 + 121234431*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        for (int i = 0; i < 1000; ++i) {
            List<MultivariatePolynomialZp64> modGens = gens.stream().map(p -> MultivariatePolynomial.asOverZp64(p.setRing(Zp(123419)))).collect(Collectors.toList());
            List<MultivariatePolynomialZp64> modGB = BuchbergerGB(modGens, GREVLEX);
            List<MultivariatePolynomial<BigInteger>> gb = solveGB(gens,
                    modGB.stream().map(p -> p.terms.keySet()).collect(Collectors.toList()), GREVLEX);
            assertTrue(isGroebnerBasis(gens, gb, GREVLEX));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testSparseGB3() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<BigInteger>
                b = parse("x1*x2*x3^2 + 123*x2*x4*x1^2*x5 + 123*x3*x2^3", Z, GREVLEX, vars),
                c = parse("17*x1^3*x2^3*x3^3*x4 - x2*x4^2*x1^4*x5 + 17*x3^3*x2 - 17*x4 - 17", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(c, b);


        List<MultivariatePolynomial<BigInteger>> expected = BuchbergerGB(gens, GREVLEX);
        List<MultivariatePolynomialZp64> modGens = gens.stream().map(p -> MultivariatePolynomial.asOverZp64(p.setRing(Zp(123419)))).collect(Collectors.toList());
        List<MultivariatePolynomialZp64> modGB = BuchbergerGB(modGens, GREVLEX);

        List<MultivariatePolynomial<BigInteger>> gb = solveGB(gens,
                modGB.stream().map(p -> p.terms.keySet()).collect(Collectors.toList()), GREVLEX);

        assertEquals(expected, gb);
    }

    @Test
    public void testSparseGB4() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                g1 = parse("2*y*z+y^2+16*y^2*z^5+16*y^3*z^4+4*y^4*z^3+12*y^2*z^6+12*y^3*z^5+3*y^4*z^4+8*y^2*z^7-24*y^3*z^6-46*y^4*z^5-24*y^5*z^4+8*x^2*y^3*z^4-4*y^6*z^3+8*x^2*y^4*z^3+2*x^2*y^5*z^2-24*y^3*z^7-36*y^4*z^6-18*y^5*z^5+16*x^2*y^3*z^5-3*y^6*z^4+16*x^2*y^4*z^4+4*x^2*y^5*z^3-16*y^3*z^8-24*y^4*z^7-12*y^5*z^6+4*x^2*y^3*z^6-2*y^6*z^5-12*x^2*y^4*z^5-23*x^2*y^5*z^4-12*x^2*y^6*z^3-2*x^2*y^7*z^2+64*y^3*z^9+96*y^4*z^8+48*y^5*z^7+8*y^6*z^6-32*x^2*y^4*z^6-48*x^2*y^5*z^5-24*x^2*y^6*z^4-4*x^2*y^7*z^3+48*y^3*z^10+72*y^4*z^9+36*y^5*z^8+6*y^6*z^7-8*x^2*y^4*z^7-12*x^2*y^5*z^6-6*x^2*y^6*z^5-x^2*y^7*z^4+32*y^3*z^11+48*y^4*z^10+24*y^5*z^9+4*y^6*z^8+64*x^2*y^4*z^8+96*x^2*y^5*z^7+48*x^2*y^6*z^6+8*x^2*y^7*z^5+88*x^2*y^4*z^9+132*x^2*y^5*z^8+32*x^3*y^4*z^8+66*x^2*y^6*z^7+48*x^3*y^5*z^7+11*x^2*y^7*z^6+24*x^3*y^6*z^6+4*x^3*y^7*z^5+32*x^2*y^4*z^10+48*x^2*y^5*z^9+24*x^3*y^4*z^9+24*x^2*y^6*z^8+36*x^3*y^5*z^8+4*x^2*y^7*z^7+18*x^3*y^6*z^7+16*x^4*y^5*z^7+3*x^3*y^7*z^6+24*x^4*y^6*z^6+12*x^4*y^7*z^5+2*x^4*y^8*z^4+16*x^3*y^4*z^10+24*x^3*y^5*z^9+12*x^3*y^6*z^8+32*x^4*y^5*z^8+2*x^3*y^7*z^7+48*x^4*y^6*z^7+16*x^5*y^5*z^7+24*x^4*y^7*z^6+24*x^5*y^6*z^6+4*x^4*y^8*z^5+12*x^5*y^7*z^5+2*x^5*y^8*z^4+8*x^4*y^5*z^9+12*x^4*y^6*z^8+32*x^5*y^5*z^8+6*x^4*y^7*z^7+48*x^5*y^6*z^7+x^4*y^8*z^6+24*x^5*y^7*z^6+4*x^5*y^8*z^5+8*x^5*y^5*z^9+12*x^5*y^6*z^8+6*x^5*y^7*z^7+x^5*y^8*z^6", Z, GREVLEX, vars),
                g2 = parse("2*y*z+y^2+6*z^3+3*z^4+2*z^5+3*x^2*y*z^2+4*x^2*y*z^3+x^3*y*z^2+x^2*y*z^4-64*y^2*z^11-64*y^3*z^10-16*y^4*z^9-48*y^2*z^12-48*y^3*z^11-12*y^4*z^10-32*y^2*z^13-32*y^3*z^12-8*y^4*z^11-96*x^2*y^3*z^10-96*x^2*y^4*z^9-24*x^2*y^5*z^8-112*x^2*y^3*z^11-112*x^2*y^4*z^10-64*x^3*y^3*z^10-28*x^2*y^5*z^9-64*x^3*y^4*z^9-16*x^3*y^5*z^8-48*x^2*y^3*z^12-48*x^2*y^4*z^11-48*x^3*y^3*z^11-12*x^2*y^5*z^10-48*x^3*y^4*z^10-12*x^3*y^5*z^9-48*x^4*y^4*z^9-48*x^4*y^5*z^8-12*x^4*y^6*z^7-32*x^3*y^3*z^12-32*x^3*y^4*z^11-8*x^3*y^5*z^10-76*x^4*y^4*z^10-76*x^4*y^5*z^9-64*x^5*y^4*z^9-19*x^4*y^6*z^8-64*x^5*y^5*z^8-16*x^5*y^6*z^7-24*x^4*y^4*z^11-24*x^4*y^5*z^10-88*x^5*y^4*z^10-6*x^4*y^6*z^9-88*x^5*y^5*z^9-16*x^6*y^4*z^9-22*x^5*y^6*z^8-24*x^6*y^5*z^8-12*x^6*y^6*z^7-2*x^6*y^7*z^6-32*x^5*y^4*z^11-32*x^5*y^5*z^10-12*x^6*y^4*z^10-8*x^5*y^6*z^9-28*x^6*y^5*z^9-19*x^6*y^6*z^8-16*x^7*y^5*z^8-4*x^6*y^7*z^7-16*x^7*y^6*z^7-4*x^7*y^7*z^6-8*x^6*y^4*z^11-12*x^6*y^5*z^10-6*x^6*y^6*z^9-32*x^7*y^5*z^9-x^6*y^7*z^8-32*x^7*y^6*z^8-8*x^8*y^5*z^8-8*x^7*y^7*z^7-8*x^8*y^6*z^7-2*x^8*y^7*z^6-8*x^7*y^5*z^10-8*x^7*y^6*z^9-16*x^8*y^5*z^9-2*x^7*y^7*z^8-16*x^8*y^6*z^8-4*x^8*y^7*z^7-4*x^8*y^5*z^10-4*x^8*y^6*z^9-x^8*y^7*z^8", Z, GREVLEX, vars),
                g3 = parse("4*z^3+3*z^4+2*z^5+2*x^2*y*z^2+4*x^2*y*z^3+16*y^2*z^5+16*y^3*z^4+x^2*y*z^4+4*y^4*z^3+12*y^2*z^6+12*y^3*z^5+3*y^4*z^4+8*y^2*z^7+8*y^3*z^6+2*y^4*z^5+8*x^2*y^3*z^4+8*x^2*y^4*z^3+2*x^2*y^5*z^2+64*y^2*z^8+64*y^3*z^7+16*y^4*z^6+16*x^2*y^3*z^5+16*x^2*y^4*z^4+4*x^2*y^5*z^3+96*y^2*z^9+96*y^3*z^8+24*y^4*z^7+4*x^2*y^3*z^6+4*x^2*y^4*z^5+x^2*y^5*z^4+100*y^2*z^10+100*y^3*z^9+25*y^4*z^8+64*x^2*y^3*z^7+64*x^2*y^4*z^6+16*x^2*y^5*z^5+48*y^2*z^11+48*y^3*z^10+12*y^4*z^9+176*x^2*y^3*z^8+176*x^2*y^4*z^7+44*x^2*y^5*z^6+16*y^2*z^12+16*y^3*z^11+4*y^4*z^10+160*x^2*y^3*z^9+160*x^2*y^4*z^8+40*x^2*y^5*z^7+16*x^4*y^4*z^6+16*x^4*y^5*z^5+4*x^4*y^6*z^4+88*x^2*y^3*z^10+88*x^2*y^4*z^9+22*x^2*y^5*z^8+64*x^4*y^4*z^7+64*x^4*y^5*z^6+16*x^4*y^6*z^5+16*x^2*y^3*z^11+16*x^2*y^4*z^10+4*x^2*y^5*z^9+80*x^4*y^4*z^8+80*x^4*y^5*z^7+20*x^4*y^6*z^6+32*x^4*y^4*z^9+32*x^4*y^5*z^8+8*x^4*y^6*z^7+4*x^4*y^4*z^10+4*x^4*y^5*z^9+x^4*y^6*z^8", Z, GREVLEX, vars),
                g4 = parse("2*y*z+y^2-16*z^6-24*z^7-25*z^8-16*x^2*y*z^5-12*z^9-44*x^2*y*z^6-4*z^10-40*x^2*y*z^7-4*x^4*y^2*z^4-22*x^2*y*z^8-16*x^4*y^2*z^5+64*y^3*z^9-4*x^2*y*z^9+96*y^4*z^8+48*y^5*z^7+8*y^6*z^6-20*x^4*y^2*z^6-128*y^2*z^11-80*y^3*z^10+40*y^4*z^9+36*y^5*z^8+6*y^6*z^7-8*x^4*y^2*z^7-192*y^2*z^12-160*y^3*z^11+24*y^5*z^9+4*y^6*z^8+64*x^2*y^4*z^8-x^4*y^2*z^8+96*x^2*y^5*z^7+48*x^2*y^6*z^6+8*x^2*y^7*z^5-200*y^2*z^13-200*y^3*z^12-50*y^4*z^11-192*x^2*y^3*z^10-104*x^2*y^4*z^9+84*x^2*y^5*z^8+32*x^3*y^4*z^8+66*x^2*y^6*z^7+48*x^3*y^5*z^7+11*x^2*y^7*z^6+24*x^3*y^6*z^6+4*x^3*y^7*z^5-96*y^2*z^14-96*y^3*z^13-24*y^4*z^12-448*x^2*y^3*z^11-416*x^2*y^4*z^10-64*x^3*y^3*z^10-64*x^2*y^5*z^9-40*x^3*y^4*z^9+24*x^2*y^6*z^8+20*x^3*y^5*z^8+4*x^2*y^7*z^7+18*x^3*y^6*z^7+16*x^4*y^5*z^7+3*x^3*y^7*z^6+24*x^4*y^6*z^6+12*x^4*y^7*z^5+2*x^4*y^8*z^4-32*y^2*z^15-32*y^3*z^14-8*y^4*z^13-420*x^2*y^3*z^12-420*x^2*y^4*z^11-96*x^3*y^3*z^11-105*x^2*y^5*z^10-80*x^3*y^4*z^10-96*x^4*y^4*z^9+12*x^3*y^6*z^8-64*x^4*y^5*z^8+2*x^3*y^7*z^7+24*x^4*y^6*z^7+16*x^5*y^5*z^7+24*x^4*y^7*z^6+24*x^5*y^6*z^6+4*x^4*y^8*z^5+12*x^5*y^7*z^5+2*x^5*y^8*z^4-224*x^2*y^3*z^13-224*x^2*y^4*z^12-100*x^3*y^3*z^12-56*x^2*y^5*z^11-100*x^3*y^4*z^11-25*x^3*y^5*z^10-304*x^4*y^4*z^10-296*x^4*y^5*z^9-64*x^5*y^4*z^9-64*x^4*y^6*z^8-32*x^5*y^5*z^8+6*x^4*y^7*z^7+32*x^5*y^6*z^7+x^4*y^8*z^6+24*x^5*y^7*z^6+4*x^5*y^8*z^5-48*x^2*y^3*z^14-48*x^2*y^4*z^13-48*x^3*y^3*z^13-12*x^2*y^5*z^12-48*x^3*y^4*z^12-12*x^3*y^5*z^11-320*x^4*y^4*z^11-320*x^4*y^5*z^10-176*x^5*y^4*z^10-80*x^4*y^6*z^9-168*x^5*y^5*z^9-32*x^5*y^6*z^8-16*x^6*y^5*z^8+6*x^5*y^7*z^7-16*x^6*y^6*z^7+x^5*y^8*z^6-4*x^6*y^7*z^6-16*x^3*y^3*z^14-16*x^3*y^4*z^13-4*x^3*y^5*z^12-152*x^4*y^4*z^12-152*x^4*y^5*z^11-160*x^5*y^4*z^11-38*x^4*y^6*z^10-160*x^5*y^5*z^10-40*x^5*y^6*z^9-64*x^6*y^5*z^9-64*x^6*y^6*z^8-16*x^7*y^5*z^8-16*x^6*y^7*z^7-16*x^7*y^6*z^7-4*x^7*y^7*z^6-24*x^4*y^4*z^13-24*x^4*y^5*z^12-88*x^5*y^4*z^12-6*x^4*y^6*z^11-88*x^5*y^5*z^11-22*x^5*y^6*z^10-80*x^6*y^5*z^10-80*x^6*y^6*z^9-64*x^7*y^5*z^9-20*x^6*y^7*z^8-64*x^7*y^6*z^8-16*x^7*y^7*z^7-16*x^5*y^4*z^13-16*x^5*y^5*z^12-4*x^5*y^6*z^11-32*x^6*y^5*z^11-32*x^6*y^6*z^10-80*x^7*y^5*z^10-8*x^6*y^7*z^9-80*x^7*y^6*z^9-20*x^7*y^7*z^8-4*x^6*y^5*z^12-4*x^6*y^6*z^11-32*x^7*y^5*z^11-x^6*y^7*z^10-32*x^7*y^6*z^10-8*x^7*y^7*z^9-4*x^7*y^5*z^12-4*x^7*y^6*z^11-x^7*y^7*z^10", Z, GREVLEX, vars),
                g5 = parse("2*z^3+x^2*y*z^2-8*z^6+x^3*y*z^2-6*z^7-4*z^8-8*x^2*y*z^5-11*x^2*y*z^6-4*x^3*y*z^5-4*x^2*y*z^7-3*x^3*y*z^6-2*x^4*y^2*z^4+64*y*z^10+32*y^2*z^9-2*x^3*y*z^7-4*x^4*y^2*z^5-2*x^5*y^2*z^4+96*y*z^11+48*y^2*z^10-x^4*y^2*z^6-4*x^5*y^2*z^5+100*y*z^12+50*y^2*z^11+96*x^2*y^2*z^9+48*x^2*y^3*z^8-x^5*y^2*z^6-80*y*z^13-40*y^2*z^12+224*x^2*y^2*z^10+112*x^2*y^3*z^9+32*x^3*y^2*z^9+16*x^3*y^3*z^8-176*y*z^14-88*y^2*z^13+210*x^2*y^2*z^11+105*x^2*y^3*z^10+48*x^3*y^2*z^10+24*x^3*y^3*z^9+48*x^4*y^3*z^8+24*x^4*y^4*z^7-200*y*z^15-100*y^2*z^14-144*x^2*y^2*z^12-72*x^2*y^3*z^11+50*x^3*y^2*z^11+25*x^3*y^3*z^10+152*x^4*y^3*z^9+76*x^4*y^4*z^8+32*x^5*y^3*z^8+16*x^5*y^4*z^7-96*y*z^16-48*y^2*z^15-520*x^2*y^2*z^13-260*x^2*y^3*z^12-104*x^3*y^2*z^12-52*x^3*y^3*z^11+160*x^4*y^3*z^10+80*x^4*y^4*z^9+88*x^5*y^3*z^9+44*x^5*y^4*z^8+8*x^6*y^4*z^7+4*x^6*y^5*z^6-32*y*z^17-16*y^2*z^16-520*x^2*y^2*z^14-260*x^2*y^3*z^13-184*x^3*y^2*z^13-92*x^3*y^3*z^12-116*x^4*y^3*z^11-58*x^4*y^4*z^10+80*x^5*y^3*z^10+40*x^5*y^4*z^9+32*x^6*y^4*z^8+16*x^6*y^5*z^7+8*x^7*y^4*z^7+4*x^7*y^5*z^6-272*x^2*y^2*z^15-136*x^2*y^3*z^14-200*x^3*y^2*z^14-100*x^3*y^3*z^13-516*x^4*y^3*z^12-258*x^4*y^4*z^11-148*x^5*y^3*z^11-74*x^5*y^4*z^10+40*x^6*y^4*z^9+20*x^6*y^5*z^8+32*x^7*y^4*z^8+16*x^7*y^5*z^7-64*x^2*y^2*z^16-32*x^2*y^3*z^15-96*x^3*y^2*z^15-48*x^3*y^3*z^14-530*x^4*y^3*z^13-265*x^4*y^4*z^12-440*x^5*y^3*z^12-220*x^5*y^4*z^11-32*x^6*y^3*z^11-64*x^6*y^4*z^10-24*x^6*y^5*z^9+40*x^7*y^4*z^9+20*x^7*y^5*z^8-32*x^3*y^2*z^16-16*x^3*y^3*z^15-264*x^4*y^3*z^14-132*x^4*y^4*z^13-420*x^5*y^3*z^13-210*x^5*y^4*z^12-48*x^6*y^3*z^12-238*x^6*y^4*z^11-107*x^6*y^5*z^10-80*x^7*y^4*z^10-40*x^7*y^5*z^9-48*x^4*y^3*z^15-24*x^4*y^4*z^14-224*x^5*y^3*z^14-112*x^5*y^4*z^13-50*x^6*y^3*z^13-265*x^6*y^4*z^12-120*x^6*y^5*z^11-302*x^7*y^4*z^11-151*x^7*y^5*z^10-32*x^8*y^4*z^10-24*x^8*y^5*z^9-4*x^8*y^6*z^8-48*x^5*y^3*z^15-24*x^5*y^4*z^14-24*x^6*y^3*z^14-120*x^6*y^4*z^13-54*x^6*y^5*z^12-320*x^7*y^4*z^12-160*x^7*y^5*z^11-88*x^8*y^4*z^11-76*x^8*y^5*z^10-16*x^8*y^6*z^9-16*x^9*y^5*z^9-8*x^9*y^6*z^8-8*x^6*y^3*z^15-20*x^6*y^4*z^14-8*x^6*y^5*z^13-152*x^7*y^4*z^13-76*x^7*y^5*z^12-80*x^8*y^4*z^12-80*x^8*y^5*z^11-20*x^8*y^6*z^10-64*x^9*y^5*z^10-32*x^9*y^6*z^9-8*x^10*y^5*z^9-4*x^10*y^6*z^8-24*x^7*y^4*z^14-12*x^7*y^5*z^13-44*x^8*y^4*z^13-38*x^8*y^5*z^12-8*x^8*y^6*z^11-80*x^9*y^5*z^11-40*x^9*y^6*z^10-32*x^10*y^5*z^10-16*x^10*y^6*z^9-8*x^8*y^4*z^14-6*x^8*y^5*z^13-x^8*y^6*z^12-32*x^9*y^5*z^12-16*x^9*y^6*z^11-40*x^10*y^5*z^11-20*x^10*y^6*z^10-4*x^9*y^5*z^13-2*x^9*y^6*z^12-16*x^10*y^5*z^12-8*x^10*y^6*z^11-2*x^10*y^5*z^13-x^10*y^6*z^12", Z, GREVLEX, vars),
                g6 = parse("8*z^3+6*z^4+4*z^5+4*x^2*y*z^2+8*x^2*y*z^3+2*x^2*y*z^4-32*y*z^10-16*y^2*z^9-24*y*z^11-12*y^2*z^10-16*y*z^12-8*y^2*z^11-48*x^2*y^2*z^9-24*x^2*y^3*z^8-56*x^2*y^2*z^10-28*x^2*y^3*z^9-32*x^3*y^2*z^9-16*x^3*y^3*z^8-24*x^2*y^2*z^11-12*x^2*y^3*z^10-24*x^3*y^2*z^10-12*x^3*y^3*z^9-24*x^4*y^3*z^8-12*x^4*y^4*z^7-16*x^3*y^2*z^11-8*x^3*y^3*z^10-38*x^4*y^3*z^9-19*x^4*y^4*z^8-32*x^5*y^3*z^8-16*x^5*y^4*z^7-12*x^4*y^3*z^10-6*x^4*y^4*z^9-44*x^5*y^3*z^9-22*x^5*y^4*z^8-8*x^6*y^3*z^8-8*x^6*y^4*z^7-2*x^6*y^5*z^6-16*x^5*y^3*z^10-8*x^5*y^4*z^9-6*x^6*y^3*z^9-11*x^6*y^4*z^8-4*x^6*y^5*z^7-8*x^7*y^4*z^7-4*x^7*y^5*z^6-4*x^6*y^3*z^10-4*x^6*y^4*z^9-x^6*y^5*z^8-16*x^7*y^4*z^8-8*x^7*y^5*z^7-4*x^8*y^4*z^7-2*x^8*y^5*z^6-4*x^7*y^4*z^9-2*x^7*y^5*z^8-8*x^8*y^4*z^8-4*x^8*y^5*z^7-2*x^8*y^4*z^9-x^8*y^5*z^8", Z, GREVLEX, vars),
                g7 = parse("2*y*z+y^2+8*y^2*z^5+8*y^3*z^4+2*y^4*z^3-16*y*z^7-8*y^2*z^6-12*y*z^8-6*y^2*z^7+32*y^3*z^6+48*y^4*z^5+24*y^5*z^4+4*x^2*y^3*z^4+4*y^6*z^3+4*x^2*y^4*z^3+x^2*y^5*z^2-8*y*z^9-4*y^2*z^8-16*x^2*y^2*z^6-8*x^2*y^3*z^5+4*x^3*y^3*z^4+4*x^3*y^4*z^3+x^3*y^5*z^2-22*x^2*y^2*z^7-11*x^2*y^3*z^6-8*x^3*y^2*z^6+16*x^2*y^4*z^5-4*x^3*y^3*z^5+24*x^2*y^5*z^4+12*x^2*y^6*z^3+2*x^2*y^7*z^2-8*x^2*y^2*z^8-4*x^2*y^3*z^7-6*x^3*y^2*z^7-3*x^3*y^3*z^6+16*x^3*y^4*z^5-4*x^4*y^3*z^5+24*x^3*y^5*z^4-2*x^4*y^4*z^4+12*x^3*y^6*z^3+2*x^3*y^7*z^2-4*x^3*y^2*z^8-2*x^3*y^3*z^7-8*x^4*y^3*z^6-4*x^4*y^4*z^5-4*x^5*y^3*z^5-2*x^5*y^4*z^4-2*x^4*y^3*z^7-x^4*y^4*z^6-8*x^5*y^3*z^6-4*x^5*y^4*z^5-2*x^5*y^3*z^7-x^5*y^4*z^6", Z, GREVLEX, vars),
                g8 = parse("2*z^3+x^2*y*z^2+x^3*y*z^2-1024*y^2*z^17-1024*y^3*z^16-256*y^4*z^15-2304*y^2*z^18-2304*y^3*z^17-576*y^4*z^16-3264*y^2*z^19-3264*y^3*z^18-816*y^4*z^17-2560*x^2*y^3*z^16-2560*x^2*y^4*z^15-640*x^2*y^5*z^14-2736*y^2*z^20-2736*y^3*z^19-684*y^4*z^18-7680*x^2*y^3*z^17-7680*x^2*y^4*z^16-1024*x^3*y^3*z^16-1920*x^2*y^5*z^15-1024*x^3*y^4*z^15-256*x^3*y^5*z^14-1632*y^2*z^21-1632*y^3*z^20-408*y^4*z^19-11040*x^2*y^3*z^18-11040*x^2*y^4*z^17-2304*x^3*y^3*z^17-2760*x^2*y^5*z^16-2304*x^3*y^4*z^16-576*x^3*y^5*z^15-2560*x^4*y^4*z^15-2560*x^4*y^5*z^14-640*x^4*y^6*z^13-576*y^2*z^22-576*y^3*z^21-144*y^4*z^20-9840*x^2*y^3*z^19-9840*x^2*y^4*z^18-3264*x^3*y^3*z^18-2460*x^2*y^5*z^17-3264*x^3*y^4*z^17-816*x^3*y^5*z^16-9600*x^4*y^4*z^16-9600*x^4*y^5*z^15-2048*x^5*y^4*z^15-2400*x^4*y^6*z^14-2048*x^5*y^5*z^14-512*x^5*y^6*z^13-128*y^2*z^23-128*y^3*z^22-32*y^4*z^21-5520*x^2*y^3*z^20-5520*x^2*y^4*z^19-2736*x^3*y^3*z^19-1380*x^2*y^5*z^18-2736*x^3*y^4*z^18-684*x^3*y^5*z^17-15120*x^4*y^4*z^17-15120*x^4*y^5*z^16-6528*x^5*y^4*z^16-3780*x^4*y^6*z^15-6528*x^5*y^5*z^15-256*x^6*y^4*z^15-1632*x^5*y^6*z^14-1536*x^6*y^5*z^14-1344*x^6*y^6*z^13-320*x^6*y^7*z^12-1920*x^2*y^3*z^21-1920*x^2*y^4*z^20-1632*x^3*y^3*z^20-480*x^2*y^5*z^19-1632*x^3*y^4*z^19-408*x^3*y^5*z^18-13740*x^4*y^4*z^18-13740*x^4*y^5*z^17-9408*x^5*y^4*z^17-3435*x^4*y^6*z^16-9408*x^5*y^5*z^16-576*x^6*y^4*z^16-2352*x^5*y^6*z^15-6336*x^6*y^5*z^15-5904*x^6*y^6*z^14-1536*x^7*y^5*z^14-1440*x^6*y^7*z^13-1536*x^7*y^6*z^13-384*x^7*y^7*z^12-320*x^2*y^3*z^22-320*x^2*y^4*z^21-576*x^3*y^3*z^21-80*x^2*y^5*z^20-576*x^3*y^4*z^20-144*x^3*y^5*z^19-7560*x^4*y^4*z^19-7560*x^4*y^5*z^18-8472*x^5*y^4*z^18-1890*x^4*y^6*z^17-8472*x^5*y^5*z^17-816*x^6*y^4*z^17-2118*x^5*y^6*z^16-11016*x^6*y^5*z^16-10404*x^6*y^6*z^15-6336*x^7*y^5*z^15-2550*x^6*y^7*z^14-6336*x^7*y^6*z^14-384*x^8*y^5*z^14-1584*x^7*y^7*z^13-704*x^8*y^6*z^13-416*x^8*y^7*z^12-80*x^8*y^8*z^11-128*x^3*y^3*z^22-128*x^3*y^4*z^21-32*x^3*y^5*z^20-2400*x^4*y^4*z^20-2400*x^4*y^5*z^19-4704*x^5*y^4*z^19-600*x^4*y^6*z^18-4704*x^5*y^5*z^18-684*x^6*y^4*z^18-1176*x^5*y^6*z^17-10204*x^6*y^5*z^17-9691*x^6*y^6*z^16-10416*x^7*y^5*z^16-2380*x^6*y^7*z^15-10416*x^7*y^6*z^15-1344*x^8*y^5*z^15-2604*x^7*y^7*z^14-3024*x^8*y^6*z^14-2016*x^8*y^7*z^13-512*x^9*y^6*z^13-420*x^8*y^8*z^12-512*x^9*y^7*z^12-128*x^9*y^8*z^11-320*x^4*y^4*z^21-320*x^4*y^5*z^20-1632*x^5*y^4*z^20-80*x^4*y^6*z^19-1632*x^5*y^5*z^19-408*x^6*y^4*z^19-408*x^5*y^6*z^18-5508*x^6*y^5*z^18-5202*x^6*y^6*z^17-9504*x^7*y^5*z^17-1275*x^6*y^7*z^16-9504*x^7*y^6*z^16-1944*x^8*y^5*z^16-2376*x^7*y^7*z^15-5304*x^8*y^6*z^15-3846*x^8*y^7*z^14-2592*x^9*y^6*z^14-840*x^8*y^8*z^13-2592*x^9*y^7*z^13-192*x^10*y^6*z^13-648*x^9*y^8*z^12-224*x^10*y^7*z^12-80*x^10*y^8*z^11-8*x^10*y^9*z^10-256*x^5*y^4*z^21-256*x^5*y^5*z^20-144*x^6*y^4*z^20-64*x^5*y^6*z^19-1584*x^6*y^5*z^19-1476*x^6*y^6*z^18-5208*x^7*y^5*z^18-360*x^6*y^7*z^17-5208*x^7*y^6*z^17-1776*x^8*y^5*z^17-1302*x^7*y^7*z^16-5056*x^8*y^6*z^16-3724*x^8*y^7*z^15-4992*x^9*y^6*z^15-820*x^8*y^8*z^14-4992*x^9*y^7*z^14-912*x^10*y^6*z^14-1248*x^9*y^8*z^13-1104*x^10*y^7*z^13-420*x^10*y^8*z^12-64*x^11*y^7*z^12-48*x^10*y^9*z^11-64*x^11*y^8*z^11-16*x^11*y^9*z^10-32*x^6*y^4*z^21-192*x^6*y^5*z^20-168*x^6*y^6*z^19-1584*x^7*y^5*z^19-40*x^6*y^7*z^18-1584*x^7*y^6*z^18-972*x^8*y^5*z^18-396*x^7*y^7*z^17-2652*x^8*y^6*z^17-1923*x^8*y^7*z^16-4768*x^9*y^6*z^16-420*x^8*y^8*z^15-4768*x^9*y^7*z^15-1632*x^10*y^6*z^15-1192*x^9*y^8*z^14-2064*x^10*y^7*z^14-840*x^10*y^8*z^13-384*x^11*y^7*z^13-108*x^10*y^9*z^12-384*x^11*y^8*z^12-32*x^12*y^7*z^12-96*x^11*y^9*z^11-32*x^12*y^8*z^11-8*x^12*y^9*z^10-192*x^7*y^5*z^20-192*x^7*y^6*z^19-336*x^8*y^5*z^19-48*x^7*y^7*z^18-756*x^8*y^6*z^18-504*x^8*y^7*z^17-2496*x^9*y^6*z^17-105*x^8*y^8*z^16-2496*x^9*y^7*z^16-1488*x^10*y^6*z^16-624*x^9*y^8*z^15-1936*x^10*y^7*z^15-820*x^10*y^8*z^14-864*x^11*y^7*z^14-112*x^10*y^9*z^13-864*x^11*y^8*z^13-192*x^12*y^7*z^13-216*x^11*y^9*z^12-192*x^12*y^8*z^12-48*x^12*y^9*z^11-48*x^8*y^5*z^20-88*x^8*y^6*z^19-52*x^8*y^7*z^18-648*x^9*y^6*z^18-10*x^8*y^8*z^17-648*x^9*y^7*z^17-816*x^10*y^6*z^17-162*x^9*y^8*z^16-1032*x^10*y^7*z^16-420*x^10*y^8*z^15-896*x^11*y^7*z^15-54*x^10*y^9*z^14-896*x^11*y^8*z^14-432*x^12*y^7*z^14-224*x^11*y^9*z^13-432*x^12*y^8*z^13-108*x^12*y^9*z^12-64*x^9*y^6*z^19-64*x^9*y^7*z^18-228*x^10*y^6*z^18-16*x^9*y^8*z^17-276*x^10*y^7*z^17-105*x^10*y^8*z^16-432*x^11*y^7*z^16-12*x^10*y^9*z^15-432*x^11*y^8*z^15-448*x^12*y^7*z^15-108*x^11*y^9*z^14-448*x^12*y^8*z^14-112*x^12*y^9*z^13-24*x^10*y^6*z^19-28*x^10*y^7*z^18-10*x^10*y^8*z^17-96*x^11*y^7*z^17-x^10*y^9*z^16-96*x^11*y^8*z^16-216*x^12*y^7*z^16-24*x^11*y^9*z^15-216*x^12*y^8*z^15-54*x^12*y^9*z^14-8*x^11*y^7*z^18-8*x^11*y^8*z^17-48*x^12*y^7*z^17-2*x^11*y^9*z^16-48*x^12*y^8*z^16-12*x^12*y^9*z^15-4*x^12*y^7*z^18-4*x^12*y^8*z^17-x^12*y^9*z^16", Z, GREVLEX, vars),
                g9 = parse("4*z^3+3*z^4+2*z^5+2*x^2*y*z^2+8*z^6+4*x^2*y*z^3+6*z^7-8*y^2*z^5-8*y^3*z^4+x^2*y*z^4-2*y^4*z^3+4*z^8+8*y*z^7+4*y^2*z^6+8*x^2*y*z^5+11*x^2*y*z^6+4*x^3*y*z^5-4*x^2*y^3*z^4-4*x^2*y^4*z^3-x^2*y^5*z^2+32*y^2*z^8+32*y^3*z^7+4*x^2*y*z^7+8*y^4*z^6+8*x^2*y^2*z^6+3*x^3*y*z^6+4*x^2*y^3*z^5-4*x^3*y^3*z^4+2*x^4*y^2*z^4-4*x^3*y^4*z^3-x^3*y^5*z^2+24*y^2*z^9+24*y^3*z^8+6*y^4*z^7+2*x^3*y*z^7+8*x^3*y^2*z^6+4*x^3*y^3*z^5+4*x^4*y^2*z^5+2*x^5*y^2*z^4+16*y^2*z^10+16*y^3*z^9+4*y^4*z^8+32*x^2*y^3*z^7+32*x^2*y^4*z^6+x^4*y^2*z^6+8*x^2*y^5*z^5+2*x^4*y^3*z^5+4*x^5*y^2*z^5+x^4*y^4*z^4+44*x^2*y^3*z^8+44*x^2*y^4*z^7+16*x^3*y^3*z^7+11*x^2*y^5*z^6+16*x^3*y^4*z^6+x^5*y^2*z^6+4*x^3*y^5*z^5+4*x^5*y^3*z^5+2*x^5*y^4*z^4+16*x^2*y^3*z^9+16*x^2*y^4*z^8+12*x^3*y^3*z^8+4*x^2*y^5*z^7+12*x^3*y^4*z^7+3*x^3*y^5*z^6+8*x^4*y^4*z^6+8*x^4*y^5*z^5+2*x^6*y^3*z^5+2*x^4*y^6*z^4+x^6*y^4*z^4+8*x^3*y^3*z^9+8*x^3*y^4*z^8+2*x^3*y^5*z^7+16*x^4*y^4*z^7+16*x^4*y^5*z^6+8*x^5*y^4*z^6+4*x^4*y^6*z^5+8*x^5*y^5*z^5+2*x^5*y^6*z^4+4*x^4*y^4*z^8+4*x^4*y^5*z^7+16*x^5*y^4*z^7+x^4*y^6*z^6+16*x^5*y^5*z^6+4*x^5*y^6*z^5+4*x^5*y^4*z^8+4*x^5*y^5*z^7+x^5*y^6*z^6", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(g1, g2, g3, g4, g5, g6, g7, g8, g9);

        GBResult<MonomialZp64, MultivariatePolynomialZp64> modGB = F4GB(mod(gens, nextPrime(1 << 25)), GREVLEX);
        List<MultivariatePolynomial<BigInteger>> sparseGB = solveGB(gens, modGB.stream().map(AMultivariatePolynomial::getSkeleton).collect(Collectors.toList()), GREVLEX);
        System.out.println(sparseGB);
    }


    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SingularResult<Term, Poly> SingularGB(List<Poly> ideal,
                                          Comparator<DegreeVector> monomialOrder) throws IOException, InterruptedException {
        return SingularGB(ideal, monomialOrder, null);
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SingularResult<Term, Poly> SingularGB(List<Poly> ideal,
                                          Comparator<DegreeVector> monomialOrder,
                                          int... eliminateVars) throws IOException, InterruptedException {
        if (!isSingularAvailable())
            throw new RuntimeException("no Singular");
        String order;
        if (monomialOrder == LEX)
            order = "lp";
        else if (monomialOrder == GREVLEX)
            order = "dp";
        else
            throw new IllegalArgumentException("bad order");

        setMonomialOrder(ideal, monomialOrder);

        Poly factory = ideal.get(0);
        IPolynomialRing<Poly> ring = Rings.PolynomialRing(factory);
        int nVars = factory.nVariables;
        String[] vars = new String[nVars];
        for (int i = 1; i <= nVars; i++)
            vars[i - 1] = "x" + i;

        // prepare tmp file
        File tmpSource = File.createTempFile("singular_" + Integer.toHexString(ideal.hashCode()), null);
        File tmpOutput = File.createTempFile("singular_" + Integer.toHexString(ideal.hashCode()), "output");
        tmpSource.deleteOnExit();
        tmpOutput.deleteOnExit();

        try {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(tmpSource))) {

                writer.write("system(\"--ticks-per-sec\",1000);");
                writer.newLine();
                writer.write("option(redSB);");
                writer.newLine();


                writer.write(String.format("ring r = %s, (%s), %s;",
                        factory.coefficientRingCardinality() == null ? "0" : factory.coefficientRingCardinality(),
                        Arrays.stream(vars).collect(Collectors.joining(",")),
                        order));
                writer.newLine();

                for (int i = 0; i < ideal.size(); i++) {
                    writer.write(String.format("poly p%s = ", i));
                    writer.write(ideal.get(i).toString(vars));
                    writer.write(";");
                    writer.newLine();
                }

                // ideal
                writer.write(String.format("ideal I = %s;",
                        IntStream.range(0, ideal.size()).mapToObj(i -> "p" + i)
                                .collect(Collectors.joining(","))));
                writer.newLine();
                writer.write("int t = timer;");
                writer.newLine();
                if (eliminateVars == null)
                    writer.write("ideal G = std(I);");
                else
                    writer.write("ideal G = std(eliminate(I," + Arrays.stream(eliminateVars).mapToObj(i -> vars[i]).collect(Collectors.joining("*")) + "));");
                writer.newLine();
                writer.write("int elapsed = timer-t;");
                writer.newLine();
                writer.write("print(\"OUTPUTSTARTSHERE\");");
                writer.newLine();
                writer.write("print(G);");
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
            if (!singularOut.contains("TIMESEPARATOR"))
                //<- there was a error in Singular
                return null;

            String[] split = singularOut.split("OUTPUTSTARTSHERE")[1].split("TIMESEPARATOR");
            // parse polynomials
            List<Poly> std = Arrays.stream(split[0].split(",")).map(str -> ring.parse(str, vars)).collect(Collectors.toList());
//            // minimize Groebner basis
//            minimizeGroebnerBases(std);
//            // reduce Groebner basis
//            removeRedundant(std);
            // canonicalize
            canonicalize(std);
            long elapsedNanos = 1000L * 1000L * Long.parseLong(split[1].trim());
            return new SingularResult<>(std, elapsedNanos);
        } finally {
            tmpSource.delete();
            tmpOutput.delete();
        }
    }

    static final class SingularResult<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
            extends ListWrapper<Poly> {
        final long nanoseconds;

        public SingularResult(List<Poly> list, long nanoseconds) {
            super(list);
            this.nanoseconds = nanoseconds;
        }

        @Override
        public String toString() {
            return list.toString() + "\n" + "Time: " + nanosecondsToString(nanoseconds);
        }
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> shuffleGB(List<Poly> basis, RandomGenerator rnd, int depth, int redundancy) {
        assert redundancy >= 2;
        List<Poly> newList = new ArrayList<>();
        for (int red = 0; red < redundancy; ++red)
            for (int i = 0; i < basis.size(); ++i) {
                Poly el = basis.get(i).clone();
                Poly tmp = el.clone();
                for (int j = 0; j < redundancy; ++j) {
                    for (int k = 0; k < depth; ++k) {
                        Poly r = basis.get(rnd.nextInt(basis.size()));
                        if (rnd.nextBoolean())
                            r = r.clone().negate();
                        if (rnd.nextBoolean())
                            tmp = tmp.add(r);
                        else
                            tmp = tmp.multiply(r);
                    }
                }
                newList.add(el.add(tmp));
            }
        return newList;
    }

    @Target(ElementType.METHOD)
    @Retention(RetentionPolicy.RUNTIME)
    @interface RequiresSingular {}
}