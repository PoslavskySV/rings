package cc.r2.core.util;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class TimeUnits {
    private TimeUnits() {}

    public static String nanosecondsToString(long nano) {
        String pf = "ns";
        if (nano / 1000 != 0) {
            pf = "us";
            nano /= 1000;
        }

        if (nano / 1000 != 0) {
            pf = "ms";
            nano /= 1000;
        }

        if (nano / 1000 != 0) {
            pf = "s";
            nano /= 1000;
        }

        return nano + pf;
    }
}
