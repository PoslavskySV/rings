package cc.redberry.rings.poly;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class UtilTest {
    @Test
    public void test1() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> poly = UnivariatePolynomial.parse("(1/2) + (1/3)*x + (1/4)*x^2 + (1/5)*x^3 + (1/6)*x^7", Rings.Q);
        System.out.println(Util.toCommonDenominator(poly));

        MultivariatePolynomial<Rational<BigInteger>> mpoly = MultivariatePolynomial.parse("(1/2) + (1/3)*x + (1/4)*x^2 + (1/5)*x^3111  + 1/5*x^66", Rings.Q);
        System.out.println(Util.toCommonDenominator(mpoly));
    }
}