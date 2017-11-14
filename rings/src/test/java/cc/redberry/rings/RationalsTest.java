package cc.redberry.rings;

import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.Frac;
import static cc.redberry.rings.Rings.UnivariateRingZp64;

/**
 * @since 1.0
 */
public class RationalsTest {
    @Test
    public void test1() throws Exception {
        for (int i = 19; i < 1000; i++)
            Assert.assertNotNull(Frac(UnivariateRingZp64(17)).randomElement());
    }
}