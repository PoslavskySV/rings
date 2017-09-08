package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.Rational;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.Rationals;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.asLongPolyZp;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class GroebnerBasisTest {

    @Test
    public void sadd() throws Exception {
String[] variables = {"x", "y"};
MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.parse("x^2 + 2*x*y + y^2", variables);

// will use the specified mapping and print x^2 + 2*x*y + y^2
System.out.println(poly.toString(variables));

// will use the default mapping and print a^2 + 2*a*b + b^2
System.out.println(poly.toString());
    }

    @Test
    public void test1() throws Exception {
        String[] vars = {"x", "y"};
        MultivariatePolynomial<Rational>
                f1 = parse("x^2*y^2+x*y", Rationals.Rationals, MonomialOrder.LEX, vars),
                f2 = parse("x*y^4-y^2", Rationals.Rationals, MonomialOrder.LEX, vars);

        System.out.println(GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2)));
    }

    @Test
    public void test2() throws Exception {
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<Rational>
                f1 = parse("x^2*y^2 + x*y + 5*z^3*y^2", Rationals.Rationals, MonomialOrder.LEX, vars),
                f2 = parse("x*y^4 - y^2 - 5*z^3*x^2", Rationals.Rationals, MonomialOrder.LEX, vars);


        List<MultivariatePolynomial<Rational>> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2));
        for (MultivariatePolynomial<Rational> t : r)
            System.out.println(t.toString(vars));
    }

    @Test
    public void test3() throws Exception {
        String[] vars = {"x", "y", "z"};
        IntegersModulo domain = new IntegersModulo(17);
        lMultivariatePolynomialZp
                f1 = asLongPolyZp(parse("x^2*y^2 + x*y + 5*z^3*y^2", domain, MonomialOrder.LEX, vars)),
                f2 = asLongPolyZp(parse("x*y^4 - y^2 - 5*z^3*x^2", domain, MonomialOrder.LEX, vars)),
                f3 = asLongPolyZp(parse("x - 1", domain, MonomialOrder.LEX, vars));


        List<lMultivariatePolynomialZp> r = GroebnerBasis.BuchbergerGB(Arrays.asList(f1, f2, f3));
        for (lMultivariatePolynomialZp t : r)
            System.out.println(t.toString(vars));
    }
}