package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;

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
        String[] vars = {"x", "y"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2*y^2+x*y", Rings.Q, MonomialOrder.LEX, vars),
                f2 = parse("x*y^4-y^2", Rings.Q, MonomialOrder.LEX, vars);

        System.out.println(GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2)));
    }

    @Test
    public void test2() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational<BigInteger>>
                f1 = parse("x^2*y^2 + x*y + 5*z^3*y^2", Rings.Q, MonomialOrder.LEX, vars),
                f2 = parse("x*y^4 - y^2 - 5*z^3*x^2", Rings.Q, MonomialOrder.LEX, vars);


        List<MultivariatePolynomial<Rational<BigInteger>>> r = GroebnerBasis.BuchbergerGB3(Arrays.asList(f1, f2));
        for (MultivariatePolynomial<Rational<BigInteger>> t : r)
            System.out.println(t.toString(vars));
    }

    @Test
    public void test3() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("x^2*y^2 + x*y + 5*z^3*y^2", domain, MonomialOrder.LEX, vars)),
                f2 = asOverZp64(parse("x*y^4 - y^2 - 5*z^3*x^2", domain, MonomialOrder.LEX, vars)),
                f3 = asOverZp64(parse("x - 1", domain, MonomialOrder.LEX, vars));


        List<MultivariatePolynomialZp64> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2, f3));
        for (MultivariatePolynomialZp64 t : r)
            System.out.println(t.toString(vars));
    }

    @Test
    public void test3a() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomialZp64
                f1 = asOverZp64(parse("11*x^3+3*x^4*y*z^2+12*x^3*y^4+7*x^2*y^6*z+16*x^3*y*z^6+13*x^4*y*z^5+7*x^5*y*z^5+9*x^4*y^3*z^5+7*x^6*y^4*z^2+10*x^4*y^4*z^5", domain, MonomialOrder.GREVLEX, vars)),
                f2 = asOverZp64(parse("5*y+14*x*y*z^6+5*x^2*y^2*z^4+16*x*y^4*z^3+9*x^5*y^3*z+7*x^2*y^5*z^3+15*x^5*y^6+16*x^5*y*z^6+9*x^2*y^5*z^6+10*x^5*y^5*z^3", domain, MonomialOrder.GREVLEX, vars)),
                f3 = asOverZp64(parse("8*y^5+14*y^3*z^3+2*x*y^2*z^4+13*x*y^3*z^3+2*x^2*y^4*z^2+16*x^3*y^3*z^4+8*x^3*y^3*z^5+11*x^6*y^4*z^2+8*x^5*y^5*z^3", domain, MonomialOrder.GREVLEX, vars));


        List<MultivariatePolynomialZp64> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2, f3));
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
                ideal.add(RandomMultivariatePolynomials.randomPolynomial(3, 4, 10, domain, MonomialOrder.GREVLEX, rnd));
            List<MultivariatePolynomialZp64> gb = GroebnerBasis.BuchbergerGB3(ideal);

            System.out.println(gb);
            System.out.println(gb.size());

            // some checks for gb
            for (MultivariatePolynomialZp64 iel : ideal)
                Assert.assertTrue(remainder(iel, gb).isZero());


            for (int j = 0; j < 10; j++) {
                MultivariatePolynomialZp64 rndPoly = RandomMultivariatePolynomials.randomPolynomial(3, 10, 15, domain, MonomialOrder.GREVLEX, rnd);
                MultivariatePolynomialZp64 rem = remainder(rndPoly, gb).monic();
                for (int k = 0; k < 5; k++) {
                    Collections.shuffle(gb, new Random(rnd.nextLong()));
                    Assert.assertEquals(rem, remainder(rndPoly, gb).monic());
                }
            }
        }
    }

    @Test
    public void test4() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f2 = parse("x^17*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f3 = parse("x^5*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;
            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb1 = GroebnerBasis.BuchbergerGB(ideal);
            System.out.println("GB1: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));


            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3 = GroebnerBasis.BuchbergerGB3(ideal);
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
        System.out.println(GroebnerBasis.BuchbergerGB3(gb)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }

    @Test
    public void test5() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f2 = parse("x^37*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f3 = parse("x^15*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb3 = GroebnerBasis.BuchbergerGB3(ideal);
            System.out.println("GB3: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            System.out.println();

            for (List<MultivariatePolynomialZp64> gb : Arrays.asList(gb3)) {
                gb.forEach(MultivariatePolynomialZp64::monic);
                gb.sort(MultivariatePolynomialZp64::compareTo);
            }

//            assert gb2.equals(gb3);
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB3(gb)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }
}