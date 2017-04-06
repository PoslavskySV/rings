package cc.r2.core.util;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

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

    public static String statisticsNanotime(DescriptiveStatistics stats) {
        return statisticsNanotime(stats, false);
    }

    public static String statisticsNanotime(DescriptiveStatistics stats, boolean median) {
        return nanosecondsToString((long) (median ? stats.getPercentile(0.5) : stats.getMean())) + " ± " + nanosecondsToString((long) stats.getStandardDeviation());
    }
}
