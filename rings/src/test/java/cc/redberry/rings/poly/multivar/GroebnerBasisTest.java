package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.Rings.Zp64;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.*;
import static cc.redberry.rings.poly.multivar.GroebnerBasis.GroebnerBasis;
import static cc.redberry.rings.poly.multivar.GroebnerBasisData.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariateDivision.remainder;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;

/**
 * @since 1.0
 */
@Ignore
public class GroebnerBasisTest extends AMultivariateTest {
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

        GroebnerBasis.canonicalize(expected);
        List<MultivariatePolynomial<Rational<BigInteger>>> actual = GroebnerBasis(Arrays.asList(a, b), GREVLEX);
        Assert.assertEquals(expected, actual);
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

        System.out.println(Arrays.toString(IntStream.range(0, vars.length + 1).map(i -> gens.stream().mapToInt(p -> p.degree(i)).sum()).toArray()));

        for (int i = 0; i < 1; i++) {
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> hm = BuchbergerHomogeneousGB(gens, GREVLEX);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(hm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            start = System.nanoTime();
            List<MultivariatePolynomial<Rational<BigInteger>>> nhm = BuchbergerGB(gens, GREVLEX);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertTrue(nhm.stream().allMatch(MultivariatePolynomial::isHomogeneous));

            Assert.assertEquals(hm, nhm);
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

        Assert.assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
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

        Assert.assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
    }


    @Test
    public void test5() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp64 ring = new IntegersZp64(17);
        MultivariatePolynomialZp64
                f1 = MultivariatePolynomialZp64.parse("x^2*y^2 + x*y + 5*z^13*y^2", ring, GREVLEX, vars),
                f2 = MultivariatePolynomialZp64.parse("x*y^24 - y^2*z^6 - 5*z^3*x^2 - 1", ring, GREVLEX, vars);

        List<MultivariatePolynomialZp64> ideal = Arrays.asList(f1, f2);

        List<MultivariatePolynomialZp64> fgb = F4GB(ideal, GREVLEX);
//        MultivariatePolynomialZp64 p = fgb.remove(3);
//        System.out.println(Arrays.toString(MultivariateDivision.divideAndRemainder(p, fgb.toArray(new MultivariatePolynomialZp64[0]))));

        Assert.assertEquals(BuchbergerGB(ideal, GREVLEX), F4GB(ideal, GREVLEX));
    }

    @Test
    public void test2a() throws Exception {
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
//        List<MultivariatePolynomialZp64> r = GroebnerBasis.F4GB(Arrays.asList(f1, f2, f3), GREVLEX);
//        System.out.println(r);

        r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2, f3), GREVLEX);
        System.out.println(r);
//        for (MultivariatePolynomialZp64 t : r)
//            System.out.println(t.toString(vars));
    }


    @Test
    public void test3Random() throws Exception {
        String[] vars = {"x", "y", "z"};
        RandomGenerator rnd = getRandom();
        IntegersZp64 domain = new IntegersZp64(17);
        int nIterations = its(100, 100);
        int nEls = 3;
        for (int i = 0; i < nIterations; i++) {
            List<MultivariatePolynomialZp64> ideal = new ArrayList<>();
            for (int j = 0; j < nEls; j++)
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(3, 4, 10, domain, GREVLEX, rnd));
            List<MultivariatePolynomialZp64> gb = GroebnerBasis.BuchbergerGB(ideal, GREVLEX);

            System.out.println(gb);
            System.out.println(gb.size());

            // some checks for gb
            for (MultivariatePolynomialZp64 iel : ideal)
                Assert.assertTrue(remainder(iel, gb).isZero());


            for (int j = 0; j < 10; j++) {
                MultivariatePolynomialZp64 rndPoly = RandomMultivariatePolynomials.randomPolynomial(3, 10, 15, domain, GREVLEX, rnd);
                MultivariatePolynomialZp64 rem = remainder(rndPoly, gb).monic();
                for (int k = 0; k < 5; k++) {
                    Collections.shuffle(gb, new Random(rnd.nextLong()));
                    Assert.assertEquals(rem, remainder(rndPoly, gb).monic());
                }
            }
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

            Assert.assertEquals(norm, sugar);
            System.out.println();
        }
    }
}