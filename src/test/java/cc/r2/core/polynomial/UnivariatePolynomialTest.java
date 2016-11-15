package cc.r2.core.polynomial;

import cc.r2.core.number.BigInteger;
import org.junit.Test;

import static cc.r2.core.number.BigIntegerRing.IntegerRing;

/**
 * Created by poslavsky on 02/11/2016.
 */
public class UnivariatePolynomialTest {
    @Test
    @SuppressWarnings("unchecked")
    public void name1() throws Exception {
        UnivariatePolynomial a = new UnivariatePolynomial(IntegerRing, new BigInteger[]{BigInteger.valueOf(1), BigInteger.valueOf(2), BigInteger.valueOf(3)});
        UnivariatePolynomial b = new UnivariatePolynomial(IntegerRing, new BigInteger[]{BigInteger.valueOf(-1), BigInteger.valueOf(-2), BigInteger.valueOf(-3)});

        System.out.println(a);
        System.out.println(a.multiplyByMonomial(BigInteger.valueOf(3), 2));

        System.out.println(a.add(b));
        System.out.println(a.subtract(b));
        System.out.println(a.multiply(b));

    }
}