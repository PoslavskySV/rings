package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.Rings.Zp64;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import static cc.redberry.rings.poly.multivar.GroebnerBasisData.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;
import static org.junit.Assert.assertEquals;

/**
 * @since 1.0
 */
@Ignore
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
    public void test2a() throws Exception {
        String[] vars = {"x1", "x2", "x3", "x4", "x5"};
        MultivariatePolynomial<Rational<BigInteger>>
                a = parse("x1^2*x2^2*x3 + x5^4 - 1 + x2^3*x4 + x3^5 + x4", Q, GREVLEX, vars),
                b = parse("x1*x2*x3^2 + x2*x4*x1^2*x5 + x3*x2^3", Q, GREVLEX, vars),
                c = parse("x1^3*x2^3*x3^3*x4 - x2*x4^2*x1^4*x5 + x3^3*x2 - x4 - 1", Q, GREVLEX, vars);
        List<MultivariatePolynomial<Rational<BigInteger>>> gens = Stream.of(a, b, c)
                .map(p -> p.homogenize(vars.length)).collect(Collectors.toList());

        System.out.println(Arrays.toString(IntStream.range(0, vars.length + 1).map(i -> gens.stream().mapToInt(p -> p.degree(i)).sum()).toArray()));

        for (int i = 0; i < 111; i++) {
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
            System.out.println(ref.equals(hm));
            System.out.println(ref.equals(nhm));
//            Assert.assertEquals(hm, nhm);
            System.out.println();

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
    public void name() throws Exception {
        IntegersZp64 ring = new IntegersZp64(17);
        List<MultivariatePolynomialZp64> ideal = GroebnerBasisData.katsura(11)
                .stream()
                .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator)))
                .map(p -> p.setOrdering(GREVLEX))
                .collect(Collectors.toList());

        for (int i = 0; i < 100; ++i) {
            long start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb = F4GB(ideal, GREVLEX);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }

//        for (MultivariatePolynomial<Rational<BigInteger>> p : GroebnerBasisData.katsura(11)) {
//            System.out.println("polynomialArray.emplace_back(\"" + p + "\");");
//        }
    }

    @Test
    @RequiresSingular
    public void test6_katsura() throws Exception {
        IntegersZp64 ring = new IntegersZp64(65521);
        for (int i = 9; i < 11; i++) {
            System.out.println(String.format("=> Katsura%s:", i));
            int nVars = i;
            List<MultivariatePolynomialZp64> ideal =
                    GroebnerBasisData.katsura(i)
                            .stream()
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator)))
                            .map(p -> p.setOrdering(GREVLEX))
//                            .map(p -> p.homogenize(nVars))
                            .collect(Collectors.toList());

            long start;

//            start = System.nanoTime();
//            List<MultivariatePolynomialZp64> actualBuchberger = BuchbergerGB(ideal, GREVLEX);
////            long buchberger = System.nanoTime() - start;
//            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
//            List<MultivariatePolynomialZp64> expected = singular.std;
//            System.out.println("   Singular  : " + TimeUnits.nanosecondsToString(singular.nanoseconds));


            STEP0 = 0;
            STEP1 = 0;
            STEP2 = 0;
            STEP3 = 0;
            STEP4 = 0;
            STEP5 = 0;
            STEP6 = 0;
            NPLUS = 0;

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + TimeUnits.nanosecondsToString(f4));
            System.out.println("   STEP0     : " + TimeUnits.nanosecondsToString(STEP0));
            System.out.println("   STEP1     : " + TimeUnits.nanosecondsToString(STEP1));
            System.out.println("   STEP2     : " + TimeUnits.nanosecondsToString(STEP2));
            System.out.println("   STEP3     : " + TimeUnits.nanosecondsToString(STEP3));
            System.out.println("   STEP4     : " + TimeUnits.nanosecondsToString(STEP4));
            System.out.println("   STEP5     : " + TimeUnits.nanosecondsToString(STEP5));
            System.out.println("   STEP6     : " + TimeUnits.nanosecondsToString(STEP6));
            System.out.println("   NPLUS     : " + TimeUnits.nanosecondsToString(NPLUS));



//            assertEquals(expected, actualBuchberger);
//            System.out.println(expected);
//            System.out.println(actualF4);
//            System.out.println(ideal);
//            assertEquals(expected, actualF4);
//
//            System.out.println("   Buchberger: " + TimeUnits.nanosecondsToString(buchberger));
            //System.out.println("   F4 rrr    : " + TimeUnits.nanosecondsToString(REDUCTION));
            System.out.println();
        }
    }

    @Test
    @RequiresSingular
    public void test6_cyclic() throws Exception {
        IntegersZp64 ring = new IntegersZp64(17);
        for (int i = 5; i < 16; i++) {
            System.out.println(String.format("=> Cyclic%s:", i));
          List<MultivariatePolynomialZp64> ideal =
                    GroebnerBasisData.cyclic(i)
                            .stream()
                            .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator)))
                            .map(p -> p.setOrdering(GREVLEX))
                            .collect(Collectors.toList());

            long start;

//            start = System.nanoTime();
//            List<MultivariatePolynomialZp64> actualBuchberger = BuchbergerGB(ideal, GREVLEX);
//            long buchberger = System.nanoTime() - start;

            SingularResult<MonomialZp64, MultivariatePolynomialZp64> singular = SingularGB(ideal, GREVLEX);
            List<MultivariatePolynomialZp64> expected = singular.std;
            System.out.println("   Singular  : " + TimeUnits.nanosecondsToString(singular.nanoseconds));

            start = System.nanoTime();
            //REDUCTION = 0;
            List<MultivariatePolynomialZp64> actualF4 = F4GB(ideal, GREVLEX);
            long f4 = System.nanoTime() - start;
            System.out.println("   F4        : " + TimeUnits.nanosecondsToString(f4));


//            assertEquals(expected, actualBuchberger);
//            System.out.println(expected);
//            System.out.println(actualF4);
            assertEquals(expected, actualF4);

//              System.out.println("   Buchberger: " + TimeUnits.nanosecondsToString(buchberger));
            //System.out.println("   F4 rrr    : " + TimeUnits.nanosecondsToString(REDUCTION));
            System.out.println();
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

        for (int i = 0; i < 1000; ++i) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> f4 = GroebnerBasis.F4GB(ideal, GREVLEX);
            assertEquals(expected.std, f4);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> buchberger = GroebnerBasis.BuchbergerGB(ideal, GREVLEX);
            assertEquals(expected.std, buchberger);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();
        }

    }

    @Test
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
    public void test2aaaa() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2*y^2 + x*y + 5*z^3*y^2", Q, LEX, vars),
                f2 = parse("x*y^4 - y^2 - 5*z^3*x^2", Q, LEX, vars);


        List<MultivariatePolynomial<Rational<BigInteger>>> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2), GREVLEX);
        for (MultivariatePolynomial<Rational<BigInteger>> t : r)
            System.out.println(t.toString(vars));
    }

    @Test
    public void test3a() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("x^2*y^2 + x*y + 5*z^3*y^2", domain, LEX, vars)),
                f2 = asOverZp64(parse("x*y^4 - y^2 - 5*z^3*x^2", domain, LEX, vars)),
                f3 = asOverZp64(parse("x - 1", domain, LEX, vars));


        List<MultivariatePolynomialZp64> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2, f3), LEX);
        for (MultivariatePolynomialZp64 t : r)
            System.out.println(t.toString(vars));
    }

    @Test
    public void test3aa() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("11*x^3+3*x^4*y*z^2+12*x^3*y^4+7*x^2*y^6*z+16*x^3*y*z^6+13*x^4*y*z^5+7*x^5*y*z^5+9*x^4*y^3*z^5+7*x^6*y^4*z^2+10*x^4*y^4*z^5", domain, GREVLEX, vars)),
                f2 = asOverZp64(parse("5*y+14*x*y*z^6+5*x^2*y^2*z^4+16*x*y^4*z^3+9*x^5*y^3*z+7*x^2*y^5*z^3+15*x^5*y^6+16*x^5*y*z^6+9*x^2*y^5*z^6+10*x^5*y^5*z^3", domain, GREVLEX, vars)),
                f3 = asOverZp64(parse("8*y^5+14*y^3*z^3+2*x*y^2*z^4+13*x*y^3*z^3+2*x^2*y^4*z^2+16*x^3*y^3*z^4+8*x^3*y^3*z^5+11*x^6*y^4*z^2+8*x^5*y^5*z^3", domain, GREVLEX, vars));


        List<MultivariatePolynomialZp64> r;
//        r = GroebnerBasis.F4GB(Arrays.asList(f1, f2, f3), GREVLEX);
//        System.out.println(r);

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2, f3);
        r = GroebnerBasis.BuchbergerGB(ideal, GREVLEX);
        System.out.println(r);
        System.out.println(SingularGB(ideal, GREVLEX));
//        for (MultivariatePolynomialZp64 t : r)
//            System.out.println(t.toString(vars));
    }


    @Test
    public void test3Random() throws Exception {
        String[] vars = {"x", "y", "z"};
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(1001, 100);
        int nEls = 3;
        for (int i = 0; i < nIterations; i++) {
            System.out.println(i);
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
    public void test4a() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^17*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x^5*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb1 = GroebnerBasis.BuchbergerGB(ideal, LEX);
            System.out.println("GB1: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));


            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3 = GroebnerBasis.BuchbergerGB(ideal, GREVLEX);
            System.out.println("GB3: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            for (List<MultivariatePolynomialZp64> gb : Arrays.asList(gb1, gb3)) {
                gb.forEach(MultivariatePolynomialZp64::monic);
                gb.sort(MultivariatePolynomialZp64::compareTo);
            }
            System.out.println(gb1.size());

            if (!gb1.equals(gb3)) {
                System.out.println(gb1);
                System.out.println(gb3);
            }
            assert gb1.equals(gb3);
            System.out.println();
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB(gb, GREVLEX)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }

    @Test
    public void test5a() throws Exception {
        String[] vars = {"t", "x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2*t + x*y*z - y*z^4 + x", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^17*y^4 - z^4*y^2 + x*z*y^4 - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x^15*y^4*z - x*y^2 + x*z + t^3", Rings.Zp(17), GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;
            Comparator<SyzygyPair> strategy = GroebnerBasis.normalSelectionStrategy(GREVLEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3 = BuchbergerGB(
                    ideal,
                    GREVLEX,
                    strategy,//System.out.printf("%s\n", cur) != null &&
                    (prev, cur) -> (cur > prev + 11128));
            System.out.println(gb3.size());
            System.out.println("Normal strategy: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBasis.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3a = BuchbergerGB(
                    ideal,
                    GREVLEX,
                    strategy,
                    (prev, cur) -> (cur > prev + 11150));
            System.out.println("Sugar strategy:  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));


            System.out.println();

            for (List<MultivariatePolynomialZp64> gb : Arrays.asList(gb3)) {
                gb.forEach(MultivariatePolynomialZp64::monic);
                gb.sort(MultivariatePolynomialZp64::compareTo);
            }

//            assert gb2.equals(gb3);
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB(gb, GREVLEX)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }

    @Test
    public void test5b() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y*z - y*z^2 + x", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^3*y^4 - z*y^2 + x*z*y - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;
            Comparator<SyzygyPair> strategy = GroebnerBasis.normalSelectionStrategy(LEX);

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3 = BuchbergerGB(
                    ideal,
                    LEX,
                    strategy,
                    (prev, cur) -> (cur > prev + 11150));
            System.out.println("Normal strategy: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            strategy = GroebnerBasis.withSugar(strategy);
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3a = BuchbergerGB(
                    ideal,
                    LEX,
                    strategy,
                    (prev, cur) -> (cur > prev + 11150));
            System.out.println("Sugar strategy:  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));


            System.out.println();

            for (List<MultivariatePolynomialZp64> gb : Arrays.asList(gb3)) {
                gb.forEach(MultivariatePolynomialZp64::monic);
                gb.sort(MultivariatePolynomialZp64::compareTo);
            }

//            assert gb2.equals(gb3);
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB(gb, GREVLEX)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }

    @Test
    public void test5aa() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), GREVLEX, vars),
                f2 = parse("x^37*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), GREVLEX, vars),
                f3 = parse("x^15*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomial<Rational<BigInteger>>> ideal = GroebnerBasisData.katsura2();
            long start;

            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> gb3 = BuchbergerGB(
                    ideal,
                    GREVLEX,
                    GroebnerBasis.normalSelectionStrategy(GREVLEX),
                    (prev, cur) -> (cur - prev) > 10);

            System.out.println("GB3: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            System.out.println();

            for (List<MultivariatePolynomial<Rational<BigInteger>>> gb : Arrays.asList(gb3)) {
                gb.forEach(MultivariatePolynomial<Rational<BigInteger>>::monic);
                gb.sort(MultivariatePolynomial<Rational<BigInteger>>::compareTo);
            }

            if (i == 0)
                System.out.println(gb3);

//            assert gb2.equals(gb3);
        }

    }

    //    @Benchmark
    @Test
    public void testKatsura_GREVLEX() throws Exception {
        for (int i = MIN_KATSURA; i < MAX_KATSURA; i++) {
            IntegersZp64 ring = Zp64(524287);
            List<MultivariatePolynomialZp64> katsura = katsura(i)
                    .stream()
                    .map(p -> p.mapCoefficients(ring, r -> ring.modulus(r.numerator.longValueExact())))
                    .collect(Collectors.toList());

            Comparator<DegreeVector> order = GREVLEX;
            Comparator<SyzygyPair> normalStrategy = GroebnerBasis.normalSelectionStrategy(order);
            Comparator<SyzygyPair> sugarStrategy = GroebnerBasis.withSugar(normalStrategy);
            MinimizationStrategy minimizationStrategy = (prev, cur) -> false;

            long start, elapsed;

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> sugar = BuchbergerGB(katsura, order, sugarStrategy, minimizationStrategy);
            elapsed = System.nanoTime() - start;
            System.out.println(String.format("Katsura%s, sugar: %s", i, TimeUnits.nanosecondsToString(elapsed)));

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> norm = BuchbergerGB(katsura, order, normalStrategy, minimizationStrategy);
            elapsed = System.nanoTime() - start;
            System.out.println(String.format("Katsura%s, normal: %s", i, TimeUnits.nanosecondsToString(elapsed)));

            canonicalize(norm);
            canonicalize(sugar);

            assertEquals(norm, sugar);
            System.out.println();
        }
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
                        factory.coefficientRingCardinality() == null ? 0 : factory.coefficientRingCardinality(),
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