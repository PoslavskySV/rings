package cc.redberry.rings.poly.univar;

import cc.redberry.rings.poly.test.APolynomialTest;
import org.junit.Assume;

/**
 * @since 1.0
 */
public class AUnivariateTest extends APolynomialTest {
    private final boolean run;

    public AUnivariateTest() {
        // prevent to switch from generic ring to 64-bit arithmetic
        Conversions64bit.SWITCH_TO_64bit = System.getProperty("test.switch64", "false").equals("true");
        this.run = !"true".equals(System.getProperty("skipUnivariate"));
    }

    @Override
    public void beforeMethod() throws Exception {
        super.beforeMethod();
        Assume.assumeTrue(run);
    }
}
