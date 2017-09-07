package cc.r2.core.poly;

import cc.r2.core.poly.multivar.MultivariatePolynomial;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class IPolynomialTest {

    @Test
    public void testComparable() throws Exception {
        assertSignum(lUnivariatePolynomialZ.create(1, 2, 3));
        lUnivariatePolynomialZ negate = lUnivariatePolynomialZ.create(1, 2, 3).negate();
        assertSignum(negate);
        assertSignum(lUnivariatePolynomialZ.create(1, 2, 3).modulus(2));
        assertSignum(lUnivariatePolynomialZ.create(1, 2, 3).modulus(2).toBigPoly());
        assertSignum(MultivariatePolynomial.parse("a^2 - 2*b"));
        assertSignum(lMultivariatePolynomialZp.parse("a^2 - 2*b", new lIntegersModulo(2)));
    }

    private static <T extends IPolynomial<T>> void assertSignum(T p) {
        T zero = p.createZero();
        if (!p.isZero())
            assertEquals(1, p.compareTo(zero));
    }
}