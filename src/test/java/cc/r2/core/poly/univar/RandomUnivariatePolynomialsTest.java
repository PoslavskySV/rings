package cc.r2.core.poly.univar;

import cc.r2.core.poly.test.APolynomialTest;
import cc.r2.core.poly.lIntegersModulo;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class RandomUnivariatePolynomialsTest extends APolynomialTest {
    @Test
    public void test1() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(94);
        lUnivariatePolynomialZp p = RandomUnivariatePolynomials.randomMonicPoly(30, 3816990131L, rnd);
        assertTrue(p.isMonic());
    }

    @Test
    public void test2() throws Exception {
        long modulus = 3816990131L;
        long a = 3567169893L;
        long b = 3117218385L;
        lIntegersModulo domain = new lIntegersModulo(modulus);
        System.out.println(domain.multiply(a, b));
        System.out.println(domain.modulusFits32);
    }
}