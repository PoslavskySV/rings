package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.test.AbstractTest;
import org.junit.Test;

/**
 * @since 1.0
 */
public class RandomMultivariatePolynomialsTest extends AbstractTest {
    @Test
    public void name() throws Exception {
        System.out.println(RandomMultivariatePolynomials.randomPolynomial(5, 5, 5, BigInteger.valueOf(123), MonomialOrder.LEX, getRandom()));
    }
}