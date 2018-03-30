package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.util.TimeUnits;
import org.junit.After;
import org.junit.Assume;

/**
 * @since 1.0
 */
public class AMultivariateTest extends APolynomialTest {
    private final boolean run;

    long start;

    public AMultivariateTest() {
        // prevent to switch from generic ring to 64-bit arithmetic
        Conversions64bit.SWITCH_TO_64bit = System.getProperty("test.switch64", "false").equals("true");
        this.run = !"true".equals(System.getProperty("skipMultivariate"));
    }

    @Override
    public void beforeMethod() throws Exception {
        super.beforeMethod();
        start = System.nanoTime();
        Assume.assumeTrue(run);
    }

    @After
    public void afterMethod() {
        System.out.println("=====> " + name.getMethodName() + "   elapsed " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
    }
}
