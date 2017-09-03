package cc.r2.core.poly;

import cc.r2.core.poly.multivar.MultivariatePolynomial;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import org.junit.Assert;
import org.junit.Test;

import java.io.Serializable;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class CommonUtilsTest {
    @Test
    public void testsSerialization1() throws Exception {
        assertSerialization(lUnivariatePolynomialZ.create(1, 2, 3));
        assertSerialization(lUnivariatePolynomialZ.create(1, 2, 3).modulus(2));
        assertSerialization(lUnivariatePolynomialZ.create(1, 2, 3).modulus(2).toBigPoly());
        assertSerialization(MultivariatePolynomial.parse("a^2 - 2*b"));
        assertSerialization(lMultivariatePolynomialZp.parse("a^2 - 2*b", new lIntegersModulo(2)));
    }

    private static <T extends Serializable> void assertSerialization(T object) {
        String compressed = CommonUtils.compress(object);
        T uncomressed = CommonUtils.uncompress(compressed);
        Assert.assertEquals(object, uncomressed);
    }
}