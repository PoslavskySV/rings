package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.TreeMap;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariatePolynomialTest {
    @Test
    public void name() throws Exception {
        MultivariatePolynomial.DegreeVector[] arr = {
                new MultivariatePolynomial.DegreeVector(2, 0, 0),
                new MultivariatePolynomial.DegreeVector(1, 1, 0),
                new MultivariatePolynomial.DegreeVector(1, 0, 1),
                new MultivariatePolynomial.DegreeVector(0, 2, 0),
                new MultivariatePolynomial.DegreeVector(0, 1, 1),
                new MultivariatePolynomial.DegreeVector(0, 0, 2)
        };

        Arrays.sort(arr, MultivariatePolynomial.GREVLEX);
        ArraysUtil.reverse(arr);
        System.out.println(Arrays.toString(arr));

        Arrays.sort(arr, MultivariatePolynomial.GRLEX);
        ArraysUtil.reverse(arr);
        System.out.println(Arrays.toString(arr));
    }

    @Test
    public void testArith1() throws Exception {
        MultivariatePolynomial a = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", MultivariatePolynomial.LEX);
        Assert.assertEquals(BigInteger.ZERO, a.cc());
        Assert.assertEquals(BigInteger.ONE, a.lc());
        Assert.assertEquals(BigInteger.ONE, a.clone().increment().cc());
        Assert.assertEquals(BigInteger.NEGATIVE_ONE, a.clone().decrement().cc());

        MultivariatePolynomial b = MultivariatePolynomial.parse("a*b - a^2 + c^3*b^2", MultivariatePolynomial.LEX);
        Assert.assertEquals(MultivariatePolynomial.parse("2*a^2", MultivariatePolynomial.LEX, "a", "b", "c"), a.clone().subtract(b));
        Assert.assertEquals(MultivariatePolynomial.parse("2*a*b + 2*c^3*b^2", MultivariatePolynomial.LEX, "a", "b", "c"), a.clone().add(b));
        Assert.assertEquals(MultivariatePolynomial.parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", MultivariatePolynomial.LEX), a.multiply(b));
    }

    @Test
    public void sadasd() throws Exception {
        TreeMap<Integer, Integer> set = new TreeMap<>();
        set.put(1, 1);
        set.put(2, 1);
        set.put(3, 1);
        System.out.println(set);
        set.compute(3, (k, v) -> null);
        System.out.println(set);

    }
}