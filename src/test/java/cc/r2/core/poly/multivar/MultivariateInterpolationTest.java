package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import org.junit.Test;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.LEX;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateInterpolationTest {
    @Test
    @SuppressWarnings("unchecked")
    public void test1() throws Exception {
        String[] variables = {"a", "b"};
        Domain<BigInteger> domain = new ModularDomain(17);
        MultivariatePolynomial<BigInteger> val1 = MultivariatePolynomial.parse("a^2 + a^3 + 1", domain, LEX, variables);
        MultivariatePolynomial<BigInteger> val2 = MultivariatePolynomial.parse("12*a^2 + 13*a^3 + 11", domain, LEX, variables);
        MultivariatePolynomial<BigInteger> val3 = MultivariatePolynomial.parse("2*a^2 + 3*a^3 + 1", domain, LEX, variables);

        BigInteger[] points = {BigInteger.valueOf(1), BigInteger.valueOf(2), BigInteger.valueOf(3)};
        MultivariatePolynomial<BigInteger>[] values = new MultivariatePolynomial[]{val1, val2, val3};
        MultivariatePolynomial<BigInteger> res = MultivariateInterpolation.interpolateNewton(1, points, values);
        assertInterpolation(1, res, points, values);
    }

    private static <E> void assertInterpolation(int variable, MultivariatePolynomial<E> poly, E[] points, MultivariatePolynomial<E>[] values) {
        for (int i = 0; i < points.length; i++)
            assertEquals(values[i], poly.evaluate(variable, points[i]));
    }
}