package cc.redberry.rings.util;

import cc.redberry.rings.bigint.BigInteger;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @since 1.0
 */
public final class RandomUtil {
    private RandomUtil() {}

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param length array length
     * @param min    min value
     * @param max    max value
     * @param rnd    random source
     * @return array of length {@code length} with elements bounded by {@code bound} (by absolute value)
     */
    public static int[] randomIntArray(int length, int min, int max, RandomGenerator rnd) {
        int[] data = new int[length];
        for (int i = 0; i < length; ++i)
            data[i] = min + rnd.nextInt(max - min);
        return data;
    }

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param length array length
     * @param total  total
     * @param rnd    random source
     * @return array of length {@code length} with elements bounded by {@code bound} (by absolute value)
     */
    public static int[] randomSharpIntArray(int length, int total, RandomGenerator rnd) {
        int[] data = new int[length];
        for (int i = 0; i < length; ++i) {
            if (total <= 0)
                break;
            data[i] = rnd.nextInt(total);
            total -= data[i];
        }
        return data;
    }

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param length array length
     * @param min    min value
     * @param max    max value
     * @param rnd    random source
     * @return array of length {@code length} with elements bounded by {@code bound} (by absolute value)
     */
    public static long[] randomLongArray(int length, long min, long max, RandomGenerator rnd) {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        long[] data = new long[length];
        for (int i = 0; i < length; ++i)
            data[i] = rndd.nextLong(min, max - 1);
        return data;
    }

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param length array length
     * @param min    min value
     * @param max    max value
     * @param rnd    random source
     * @return array of length {@code length} with elements bounded by {@code bound} (by absolute value)
     */
    public static BigInteger[] randomBigIntegerArray(int length, BigInteger min, BigInteger max, RandomGenerator rnd) {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        BigInteger[] data = new BigInteger[length];
        BigInteger delta = max.subtract(min);
        for (int i = 0; i < length; ++i)
            data[i] = min.add(randomInt(delta, rnd));
        return data;
    }

    /**
     * Returns random integer in range {@code [0, bound)}.
     *
     * @param bound maximal allowed value
     * @param rnd   random
     * @return a BigInteger {@code b} so that {@code 0 <= b < bound}
     */
    public static BigInteger randomInt(BigInteger bound, RandomGenerator rnd) {
        BigInteger r;
        do {
            r = new BigInteger(bound.bitLength(), rnd);
        } while (r.compareTo(bound) >= 0);
        return r;
    }
}
