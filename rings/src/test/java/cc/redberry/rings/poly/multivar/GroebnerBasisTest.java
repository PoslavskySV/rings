package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.GroebnerBasis.*;
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
import static cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class GroebnerBasisTest extends AMultivariateTest {
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

        List<MultivariatePolynomial<Rational<BigInteger>>> bHom = BuchbergerHomogeneousGB(gens, GREVLEX);
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
            List<MultivariatePolynomial<Rational<BigInteger>>> hm = BuchbergerHomogeneousGB(gens, GREVLEX);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(hm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> nhm = BuchbergerGB(gens, GREVLEX);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(nhm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            List<MultivariatePolynomial<Rational<BigInteger>>> ref = SingularGB(gens, GREVLEX).std;
            Assert.assertEquals(hm, nhm);
            Assert.assertEquals(ref, nhm);
        }
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
        IntegersZp64 ring = new IntegersZp64(17);
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
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator)))
                            .map(p -> p.setOrdering(GREVLEX))
                            .collect(Collectors.toList());

            long start;

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = singular.std;
            System.out.println("   Singular  : " + TimeUnits.nanosecondsToString(singular.nanoseconds));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + TimeUnits.nanosecondsToString(f4));
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
        for (int i = 5; i <= its(7, 8); i++) {
            System.out.println(String.format("=> Cyclic%s:", i));
            List<MultivariatePolynomialZp64> ideal =
                    GroebnerBasisData.cyclic(i)
                            .stream()
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator)))
                            .map(p -> p.setOrdering(GREVLEX))
                            .collect(Collectors.toList());

            long start;

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = singular.std;
            System.out.println("   Singular  : " + TimeUnits.nanosecondsToString(singular.nanoseconds));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + TimeUnits.nanosecondsToString(f4));
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

        SingularResult<MonomialZp64, MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX);

        for (int i = 0; i < its(1, 2); ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> f4 = GroebnerBasis.F4GB(ideal, GREVLEX);
            assertEquals(expected.std, f4);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> buchberger = SingularGB(ideal, GREVLEX).std;
            assertEquals(expected.std, buchberger);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();
        }

        //1713ms
        //675ms

        //1520ms
        //731ms
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
        assertEquals(expected.std, buchberger);
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
        assertEquals(expected.std, buchberger);
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
            List<MultivariatePolynomialZp64> expected = SingularGB(ideal, GREVLEX).std;
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
            List<MultivariatePolynomialZp64> buch = GroebnerBasis.BuchbergerGB(ideal, GREVLEX);
            System.out.println("Buchberger: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> f4 = GroebnerBasis.F4GB(ideal, GREVLEX);
            System.out.println("F4: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

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
            Comparator<SyzygyPair> strategy = GroebnerBasis.normalSelectionStrategy(GREVLEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> norm = BuchbergerGB(ideal, GREVLEX, strategy, NO_MINIMIZATION);
            System.out.println("Normal strategy: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBasis.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> sugar = BuchbergerGB(ideal, GREVLEX, strategy, NO_MINIMIZATION);
            System.out.println("Sugar strategy:  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

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
            Comparator<SyzygyPair> strategy = GroebnerBasis.normalSelectionStrategy(LEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> norm = BuchbergerGB(ideal, LEX, strategy, NO_MINIMIZATION);
            System.out.println("Normal strategy: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBasis.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> sugar = BuchbergerGB(ideal, LEX, strategy, NO_MINIMIZATION);
            System.out.println("Sugar strategy:  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

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
            System.out.println("F4       : " + TimeUnits.nanosecondsToString(time));

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            tSingular.addValue(singular.nanoseconds);
            System.out.println("Singular : " + TimeUnits.nanosecondsToString(singular.nanoseconds));

            List<MultivariatePolynomialZp64> expected = singular.std;
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

            List<MultivariatePolynomialZp64> expected = singular.std;
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
            System.out.println("F4       : " + TimeUnits.nanosecondsToString(time));

            //start = System.nanoTime();
            //List<MultivariatePolynomial<BigInteger>> buch = BuchbergerGB(ideal, GREVLEX);
            //time = System.nanoTime() - start;
            //tBuch.addValue(time);
            //System.out.println("Buch     : " + TimeUnits.nanosecondsToString(time));

            SingularResult<?, MultivariatePolynomial<BigInteger>> singular = SingularGB(ideal, GREVLEX);
            tSingular.addValue(singular.nanoseconds);
            System.out.println("Singular : " + TimeUnits.nanosecondsToString(singular.nanoseconds));
            System.out.println();

            List<MultivariatePolynomial<BigInteger>> expected = singular.std;
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

        List<MultivariatePolynomial<Rational<BigInteger>>> expected = BuchbergerGB(gens, GREVLEX, NO_MINIMIZATION,
                () -> new SyzygyTreeSet<>(new TreeSet<>(normalSelectionStrategy(GREVLEX))));

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
        List<MultivariatePolynomial<Rational<BigInteger>>> expected = BuchbergerGB(gens, GREVLEX, NO_MINIMIZATION,
                () -> new SyzygyTreeSet<>(new TreeSet<>(normalSelectionStrategy(GREVLEX))));

        assertEquals(expected, actual);
        List<MultivariatePolynomial<Rational<BigInteger>>> f4 = F4GB(gens, GREVLEX);
        assertEquals(expected, f4);
    }

    @Test
    public void test17() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4"};
        MultivariatePolynomial<BigInteger>
                a = parse("6*x2*x4^3 + 11*x1*x3^3*x4 + 15*x2^3*x3^2*x4 + 13*x1^3*x2^3*x4", Z, GREVLEX, vars),
                b = parse("11*x1^3 + 13*x3^2*x4^2 + x1^3*x2^3*x3 + 10*x1^3*x2^2*x3^2*x4", Z, GREVLEX, vars),
                c = parse("7*x1^2*x2*x4 + 4*x1^3*x3 + 12*x1^2*x2^2*x3^2*x4^2", Z, GREVLEX, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c);

        List<MultivariatePolynomial<Rational<MultivariatePolynomial<BigInteger>>>> id = gens.stream()
                .map(p -> MultivariateConversions.split(p, 2, 3))
                .map(p -> p.mapCoefficients(Frac(p.ring), cf -> new Rational<>(p.ring, cf)))
                .collect(Collectors.toList());
        id = BuchbergerGB(id, GREVLEX);
        assertTrue(id.get(0).isConstant());
    }

    @Test
    @RequiresSingular
    public void test18() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2*y^2 + x*y - z", Q, GREVLEX, vars),
                f2 = parse("x^17*y - y^2 + x*z*y - 1", Q, GREVLEX, vars);

        List<MultivariatePolynomial<Rational<BigInteger>>> ideal = Arrays.asList(f1, f2);
        List<MultivariatePolynomial<Rational<BigInteger>>> sing = SingularGB(ideal, GREVLEX).std;
        List<MultivariatePolynomial<Rational<BigInteger>>> buch = canonicalize(GroebnerBasis.BuchbergerGB(ideal, GREVLEX));
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
        List<MultivariatePolynomial<Rational<BigInteger>>> sing = SingularGB(ideal, GREVLEX).std;
        List<MultivariatePolynomial<Rational<BigInteger>>> buch = canonicalize(GroebnerBasis.BuchbergerGB(ideal, GREVLEX));
        Assert.assertEquals(sing, buch);
    }

    @Test
    @Ignore // fixme: remove @Ignore when modular algorithm will be implemented!
    public void test20() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("(-807152359/1175978805)*z+(1708357903/571090061)*x^2*y^2*z^2+(39119166838761599/323038390588954371)*x^3*y^2*z^3", Q, GREVLEX, vars),
                b = parse("(-960519798/1555504243)*x*z^3-(1278846706/1239147733)*x*y^3-(62586766/904196831)*x*y*z^3-(792306301/1609075855)*x^3*y^3*z^3", Q, GREVLEX, vars),
                c = parse("(9306287/4567935)*x^2-(1422841761/1340607578)*x*y*z^3-(115093936/778347949)*x^3*z^3-(44182447/32319755)*x^2*y^3*z^3", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Arrays.asList(a, b, c);

        List<MultivariatePolynomial<Rational<BigInteger>>> actual = BuchbergerGB(gens, GREVLEX);
        assertTrue(actual.stream().allMatch(p -> p.stream().allMatch(Rational::isIntegral)));
        canonicalize(actual);
        System.out.println("s");
        List<MultivariatePolynomial<Rational<BigInteger>>> expected = BuchbergerGB(gens, GREVLEX, NO_MINIMIZATION,
                () -> new SyzygyTreeSet<>(new TreeSet<>(normalSelectionStrategy(GREVLEX))));

        assertEquals(expected, actual);
        List<MultivariatePolynomial<Rational<BigInteger>>> f4 = F4GB(gens, GREVLEX);
        assertEquals(expected, f4);
    }

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    SingularResult<Term, Poly> SingularGB(List<Poly> ideal, Comparator<DegreeVector> monomialOrder) throws IOException, InterruptedException {
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
        int nVars = factory.nVariables;
        String[] vars = new String[nVars];
        for (int i = 1; i <= nVars; i++)
            vars[i - 1] = "x" + i;

        // prepare tmp file
        File tmpSource = File.createTempFile("singular_" + Integer.toHexString(ideal.hashCode()), null);
        File tmpOutput = File.createTempFile("singular_" + Integer.toHexString(ideal.hashCode()), "output");
        tmpSource.deleteOnExit();
        tmpOutput.deleteOnExit();

        Runtime.getRuntime().addShutdownHook(new Thread(tmpSource::delete));
        Runtime.getRuntime().addShutdownHook(new Thread(tmpOutput::delete));
        try {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter(tmpSource))) {

                writer.write("system(\"--ticks-per-sec\",1000);");
                writer.newLine();
                writer.write("option(redSB);");
                writer.newLine();


                writer.write(String.format("ring r = %s, (%s), %s;",
                        factory.coefficientRingCardinality() == null ? "QQ" : factory.coefficientRingCardinality(),
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
                writer.write("ideal G = std(I);");
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
            List<Poly> std = Arrays.stream(split[0].split(",")).map(str -> factory.parsePoly(str, vars)).collect(Collectors.toList());
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

    static final class SingularResult<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        final List<Poly> std;
        final long nanoseconds;

        SingularResult(List<Poly> std, long nanoseconds) {
            this.std = std;
            this.nanoseconds = nanoseconds;
        }

        @Override
        public String toString() {
            return std.toString() + "\n" + "Time: " + TimeUnits.nanosecondsToString(nanoseconds);
        }
    }


    @Target(ElementType.METHOD)
    @Retention(RetentionPolicy.RUNTIME)
    @interface RequiresSingular {}
}