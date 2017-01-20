package cc.r2.core.polynomial;

import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

/**
 * Created by poslavsky on 20/01/2017.
 */
public class RandomPolynomialsTest {
    @Test
    public void test1() throws Exception {
        RandomPolynomials.randomPoly(10, 2, new Well1024a());
        RandomPolynomials.randomMonicPoly(10, 2, new Well1024a());
    }
}