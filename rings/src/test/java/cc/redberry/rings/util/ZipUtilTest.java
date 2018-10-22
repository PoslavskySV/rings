package cc.redberry.rings.util;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import org.junit.Assert;
import org.junit.Test;

import java.io.Serializable;

import static cc.redberry.rings.Rings.Frac;
import static cc.redberry.rings.Rings.MultivariateRing;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class ZipUtilTest {
    @Test
    public void testsSerialization1() throws Exception {
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3));
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2));
        assertSerialization(UnivariatePolynomialZ64.create(1, 2, 3).modulus(2).toBigPoly());
        assertSerialization(MultivariatePolynomial.parse("a^2 - 2*b"));
        assertSerialization(MultivariatePolynomialZp64.parse("a^2 - 2*b", new IntegersZp64(2)));
    }

    @Test
    public void testsSerialization2() throws Exception {
        MultivariatePolynomialZp64 p = MultivariatePolynomialZp64.parse("a^2 - 2*b", new IntegersZp64(2));
        Rationals<MultivariatePolynomialZp64> fracs = Frac(MultivariateRing(p));
        assertSerialization(fracs);
        assertSerialization(fracs.mkNumerator(p));
    }

    private static <T extends Serializable> void assertSerialization(T object) {
        String compressed = ZipUtil.compress(object);
        T uncomressed = ZipUtil.uncompress(compressed);
        Assert.assertEquals(object, uncomressed);
    }
}