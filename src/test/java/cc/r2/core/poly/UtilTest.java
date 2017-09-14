package cc.r2.core.poly;

import cc.r2.core.poly.multivar.MultivariatePolynomial;
import cc.r2.core.poly.multivar.MultivariatePolynomialZp64;
import cc.r2.core.poly.univar.UnivariatePolynomialZ64;
import org.junit.Assert;
import org.junit.Test;

import java.io.Serializable;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class UtilTest {
    @Test
    public void testsSerialization1() throws Exception {
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3));
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2));
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2).toBigPoly());
        assertSerialization(MultivariatePolynomial.parse("a^2 - 2*b"));
        assertSerialization(MultivariatePolynomialZp64.parse("a^2 - 2*b", new IntegersZp64(2)));
    }

    private static <T extends Serializable> void assertSerialization(T object) {
        String compressed = Util.compress(object);
        T uncomressed = Util.uncompress(compressed);
        Assert.assertEquals(object, uncomressed);
    }
}