package cc.r2.core.poly;

import cc.r2.core.poly.multivar.MultivariatePolynomial;
import cc.r2.core.poly.multivar.MultivariatePolynomialZp64;
import cc.r2.core.poly.univar.UnivariatePolynomialZ64;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class IPolynomialTest {

    @Test
    public void testComparable() throws Exception {
        assertSignum(UnivariatePolynomialZ64.create(1, 2, 3));
        UnivariatePolynomialZ64 negate = UnivariatePolynomialZ64.create(1, 2, 3).negate();
        assertSignum(negate);
        assertSignum(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2));
        assertSignum(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2).toBigPoly());
        assertSignum(MultivariatePolynomial.parse("a^2 - 2*b"));
        assertSignum(MultivariatePolynomialZp64.parse("a^2 - 2*b", new IntegersZp64(2)));
    }

    private static <T extends IPolynomial<T>> void assertSignum(T p) {
        T zero = p.createZero();
        if (!p.isZero())
            assertEquals(1, p.compareTo(zero));
    }
}