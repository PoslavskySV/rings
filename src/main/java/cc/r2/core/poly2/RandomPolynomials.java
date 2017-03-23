package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Util to generate random polynomials.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class RandomPolynomials {
    private RandomPolynomials() {}

    /**
     * Creates random array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return array of length {@code degree + 1} with elements bounded by {@code bound} (by absolute value)
     */
    public static long[] randomLongArray(int degree, long bound, RandomGenerator rnd) {
        long[] data = new long[degree + 1];
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        for (int i = 0; i <= degree; ++i) {
            data[i] = rndd.nextLong(0, bound - 1);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = -data[i];
        }
        while (data[degree] == 0)
            data[degree] = rndd.nextLong(0, bound - 1);
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
    public static BigInteger[] randomBigArray(int degree, BigInteger bound, RandomGenerator rnd) {
        long lBound = bound.isLong() ? bound.longValue() : Long.MAX_VALUE;
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        BigInteger[] data = new BigInteger[degree + 1];
        for (int i = 0; i <= degree; ++i) {
            data[i] = randomInt(bound, rnd);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = data[i].negate();
        }
        while (data[degree].equals(BigInteger.ZERO))
            data[degree] = BigInteger.valueOf(rndd.nextLong(0, lBound));
        return data;
    }

    public static BigInteger randomInt(BigInteger bound, RandomGenerator rnd) {
        BigInteger r;
        do {
            r = new BigInteger(bound.bitLength(), rnd);
        } while (r.compareTo(bound) >= 0);
        return r;
    }

    private static final int DEFAULT_BOUND = 100;

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static MutablePolynomialZ randomPoly(int degree, RandomGenerator rnd) {
        return randomPoly(degree, DEFAULT_BOUND, rnd);
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static MutablePolynomialMod randomMonicPoly(int degree, long modulus, RandomGenerator rnd) {
        MutablePolynomialZ r = randomPoly(degree, modulus, rnd);
        while (r.data[degree] % modulus == 0) {r.data[r.degree] = rnd.nextLong();}
        return r.modulus(modulus, false).monic();
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static bMutablePolynomialMod randomMonicPoly(int degree, BigInteger modulus, RandomGenerator rnd) {
        bMutablePolynomialZ r = randomPoly(degree, modulus, rnd);
        while ((r.data[degree].mod(modulus)).isZero()) {r.data[r.degree] = rnd.nextLong();}
        return r.modulus(modulus, false).monic();
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static MutablePolynomialZ randomPoly(int degree, long bound, RandomGenerator rnd) {
        return MutablePolynomialZ.create(randomLongArray(degree, bound, rnd));
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static bMutablePolynomialZ randomPoly(int degree, BigInteger bound, RandomGenerator rnd) {
        return bMutablePolynomialZ.create(randomBigArray(degree, bound, rnd));
    }
}
