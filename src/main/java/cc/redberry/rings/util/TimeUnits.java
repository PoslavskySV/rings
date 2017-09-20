package cc.redberry.rings.util;

import org.apache.commons.math3.exception.MathIllegalStateException;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * @since 1.0
 */
public final class TimeUnits {
    private TimeUnits() {}

    public static String nanosecondsToString(long nano) {
        String pf = "ns";
        if (nano / 1000 > 1) {
            pf = "us";
            nano /= 1000;
        }

        if (nano / 1000 > 1) {
            pf = "ms";
            nano /= 1000;
        }

        if (nano / 1000 > 1) {
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

    private static String ns(double nano) { return nanosecondsToString((long) nano);}

    public static String statisticsNanotimeFull(DescriptiveStatistics stats) {
        StringBuilder outBuffer = new StringBuilder();
        String endl = "\n";
        outBuffer.append("DescriptiveStatistics:").append(endl);
        outBuffer.append("n: ").append(stats.getN()).append(endl);
        outBuffer.append("min: ").append(ns(stats.getMin())).append(endl);
        outBuffer.append("max: ").append(ns(stats.getMax())).append(endl);
        outBuffer.append("mean: ").append(ns(stats.getMean())).append(endl);
        outBuffer.append("std dev: ").append(ns(stats.getStandardDeviation()))
                .append(endl);
        try {
            // No catch for MIAE because actual parameter is valid below
            outBuffer.append("median: ").append(ns(stats.getPercentile(50))).append(endl);
        } catch (MathIllegalStateException ex) {
            outBuffer.append("median: unavailable").append(endl);
        }
        outBuffer.append("skewness: ").append(ns(stats.getSkewness())).append(endl);
        outBuffer.append("kurtosis: ").append(ns(stats.getKurtosis())).append(endl);
        return outBuffer.toString();
    }
}
