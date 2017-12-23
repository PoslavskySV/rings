package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.Q;

/**
 * @since 2.3
 */
public class IdealTest {
    @Test
    public void test1() throws Exception {
        String[] vars = {"x", "y", "z"};
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                I1 = Ideal.parse(new String[]{"x", "y"}, Q, MonomialOrder.LEX, vars),
                I2 = Ideal.parse(new String[]{"y^2", "z"}, Q, MonomialOrder.LEX, vars),
                I3 = Ideal.parse(new String[]{"y^2", "x*z", "y*z"}, Q, MonomialOrder.LEX, vars);
        Assert.assertEquals(I3, I1.intersection(I2));
    }


    @Test
    public void test2() throws Exception {
        String[] vars = {"x", "y", "z"};
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                I1 = Ideal.parse(new String[]{"x^2*y - 1", "x + y - z"}, Q, MonomialOrder.GREVLEX, vars),
                I2 = Ideal.parse(new String[]{"x^2*z^2 - y^2", "x^3 + y - z^2 - 1"}, Q, MonomialOrder.GREVLEX, vars),
                I3 = Ideal.parse(new String[]{"y^2", "x*z", "y*z"}, Q, MonomialOrder.GREVLEX, vars);
        System.out.println(I1);
        System.out.println(I2);
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> inter = I1.intersection(I2);
        System.out.println(inter);
    }


}