package cc.r2.core.poly.multivar;

import cc.r2.core.bigint.BigInteger;
import cc.r2.core.test.AbstractTest;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class RandomMultivariatePolynomialsTest extends AbstractTest{
    @Test
    public void name() throws Exception {
        System.out.println(RandomMultivariatePolynomials.randomPolynomial(5,5,5, BigInteger.valueOf(123), MonomialOrder.LEX, getRandom()));
    }
}