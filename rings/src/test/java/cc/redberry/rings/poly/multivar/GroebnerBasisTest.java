package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.util.TimeUnits;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.asOverZp64;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;

/**
 * @since 1.0
 */
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


        List<MultivariatePolynomial<Rational<BigInteger>>> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2));
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

    @Ignore
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
            List<MultivariatePolynomialZp64> gb2 = GroebnerBasis.BuchbergerGB2(ideal);
            System.out.println("GB2: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            for (List<MultivariatePolynomialZp64> gb : Arrays.asList(gb1, gb2)) {
                gb.forEach(MultivariatePolynomialZp64::monic);
                gb.sort(MultivariatePolynomialZp64::compareTo);
            }
            System.out.println(gb1.size());

            if (!gb1.equals(gb2)) {
                System.out.println(gb1);
                System.out.println(gb2);
            }
            assert gb1.equals(gb2);
            System.out.println();
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB2(gb)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }

    @Ignore
    @Test
    public void test5() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                f1 = parse("x^2*y^2 + x*y - z", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f2 = parse("x^17*y^4 - y^2 + x*z*y - 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars),
                f3 = parse("x^15*y^4*z - x*y^2 + x*z + 1", Rings.Zp(17), MonomialOrder.GREVLEX, vars);

        for (int i = 0; i < 1000; i++) {
            List<MultivariatePolynomialZp64> ideal = Arrays.asList(asOverZp64(f1), asOverZp64(f2), asOverZp64(f3));
            long start;

            start = System.nanoTime();
            List<MultivariatePolynomialZp64> gb2 = GroebnerBasis.BuchbergerGB2(ideal);
            System.out.println("GB2: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));


            System.out.println(gb2.size());
        }

        List<MultivariatePolynomial<BigInteger>> gb = Arrays.asList(f1, f2, f3);
        System.out.println(GroebnerBasis.BuchbergerGB2(gb)
                .stream()
                .map(p -> p.toString(vars))
                .collect(Collectors.joining(", ")));
    }
}