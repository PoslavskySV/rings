package cc.r2.core.polynomial;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Util for generating random polynomials
 */
public final class RandomPolynomials {
    private RandomPolynomials() {
    }

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value)
     */
    public static long[] randomLongArray(int degree, int bound, RandomGenerator rnd) {
        long[] data = new long[degree + 1];
        for (int i = 0; i <= degree; ++i) {
            data[i] = rnd.nextInt(bound);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = -data[i];
        }
        while (data[degree] == 0)
            data[degree] = rnd.nextInt(bound);
        return data;
    }

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value)
     */
    public static BigInteger[] randomBigArray(int degree, long bound, RandomGenerator rnd) {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        BigInteger[] data = new BigInteger[degree + 1];
        for (int i = 0; i <= degree; ++i) {
            data[i] = BigInteger.valueOf(rndd.nextLong(0, bound));
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = data[i].negate();
        }
        while (data[degree].isZero())
            data[degree] = BigInteger.valueOf(rndd.nextLong(0, bound));
        return data;
    }

    private static final int DEFAULT_BOUND = 100;

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static MutableLongPoly randomPoly(int degree, RandomGenerator rnd) {
        return randomPoly(degree, DEFAULT_BOUND, rnd);
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static MutableLongPoly randomMonicPoly(int degree, long modulus, RandomGenerator rnd) {
        MutableLongPoly r = randomPoly(degree, (int) modulus - 1, rnd);
        while (r.data[degree] % modulus == 0) {r.data[r.degree] = rnd.nextLong();}
        r.modulus(modulus);
        r.monic(modulus);
        return r;
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static MutableLongPoly randomPoly(int degree, int bound, RandomGenerator rnd) {
        return MutableLongPoly.create(randomLongArray(degree, bound, rnd));
    }
}
