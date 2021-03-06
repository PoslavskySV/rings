package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp64;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class RandomUnivariatePolynomialsTest extends AUnivariateTest {
    @Test
    public void test1() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(94);
        UnivariatePolynomialZp64 p = RandomUnivariatePolynomials.randomMonicPoly(30, 3816990131L, rnd);
        assertTrue(p.isMonic());
    }

    @Test
    public void test2() throws Exception {
        long modulus = 3816990131L;
        long a = 3567169893L;
        long b = 3117218385L;
        IntegersZp64 domain = new IntegersZp64(modulus);
        System.out.println(domain.multiply(a, b));
        System.out.println(domain.modulusFits32);
    }
}