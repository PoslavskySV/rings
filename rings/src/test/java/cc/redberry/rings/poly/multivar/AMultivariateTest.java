package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.poly.test.APolynomialTest;

/**
 * @since 1.0
 */
public class AMultivariateTest extends APolynomialTest {
    public AMultivariateTest() {
        // prevent to switch from generic ring to 64-bit arithmetic
        Conversions64bit.SWITCH_TO_64bit = System.getProperty("test.switch64", "false").equals("true");
    }
}
