package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.util.RandomUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * Methods to generate random polynomials.
 *

 * @since 1.0
 */
public final class RandomUnivariatePolynomials {
    private RandomUnivariatePolynomials() {}

    private static final int DEFAULT_BOUND = 100;

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param factory the factory (used to infer the type and the ring)
     * @param degree  polynomial degree
     * @param rnd     random source
     * @return random polynomial of specified {@code degree}
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly randomPoly(Poly factory, int degree, RandomGenerator rnd) {
        if (factory instanceof UnivariatePolynomialZ64)
            return (Poly) randomPoly(degree, rnd);
        else if (factory instanceof UnivariatePolynomialZp64)
            return (Poly) randomMonicPoly(degree, ((UnivariatePolynomialZp64) factory).modulus(), rnd);
        else if (factory instanceof UnivariatePolynomial) {
            UnivariatePolynomial p = randomPoly(degree, ((UnivariatePolynomial) factory).ring, rnd);
            if (factory.isOverField())
                p = p.monic();
            return (Poly) p;
        }
        throw new RuntimeException(factory.getClass().toString());
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static UnivariatePolynomialZ64 randomPoly(int degree, RandomGenerator rnd) {
        return randomPoly(degree, DEFAULT_BOUND, rnd);
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static UnivariatePolynomialZp64 randomMonicPoly(int degree, long modulus, RandomGenerator rnd) {
        UnivariatePolynomialZ64 r = randomPoly(degree, modulus, rnd);
        while (r.data[degree] % modulus == 0) {
            r.data[r.degree] = rnd.nextLong();
        }
        return r.modulus(modulus, false).monic();
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static UnivariatePolynomial<BigInteger> randomMonicPoly(int degree, BigInteger modulus, RandomGenerator rnd) {
        UnivariatePolynomial<BigInteger> r = randomPoly(degree, modulus, rnd);
        while ((r.data[degree].mod(modulus)).isZero()) {r.data[r.degree] = RandomUtil.randomInt(modulus, rnd);}
        return r.setRing(new IntegersZp(modulus)).monic();
    }

    /**
     * Creates random polynomial of specified {@code degree}.
     *
     * @param degree polynomial degree
     * @param ring   the ring
     * @param rnd    random source
     * @return random polynomial of specified {@code degree}
     */
    public static <E> UnivariatePolynomial<E> randomMonicPoly(int degree, Ring<E> ring, RandomGenerator rnd) {
        return randomPoly(degree, ring, rnd).monic();
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute
     * value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static UnivariatePolynomialZ64 randomPoly(int degree, long bound, RandomGenerator rnd) {
        return UnivariatePolynomialZ64.create(randomLongArray(degree, bound, rnd));
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute
     * value).
     *
     * @param degree polynomial degree
     * @param bound  absolute bound for coefficients
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static UnivariatePolynomial<BigInteger> randomPoly(int degree, BigInteger bound, RandomGenerator rnd) {
        return UnivariatePolynomial.create(Rings.Z, randomBigArray(degree, bound, rnd));
    }

    /**
     * Creates random polynomial of specified {@code degree} with elements from specified ring
     *
     * @param degree polynomial degree
     * @param ring   the ring
     * @param rnd    random source
     * @return random polynomial of specified {@code degree} with elements bounded by {@code bound} (by absolute value)
     */
    public static <E> UnivariatePolynomial<E> randomPoly(int degree, Ring<E> ring, RandomGenerator rnd) {
        return UnivariatePolynomial.create(ring, randomArray(degree, ring, rnd));
    }

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
            data[i] = RandomUtil.randomInt(bound, rnd);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = data[i].negate();
        }
        while (data[degree].equals(BigInteger.ZERO))
            data[degree] = BigInteger.valueOf(rndd.nextLong(0, lBound));
        return data;
    }

    /**
     * Creates random array of length {@code degree + 1} with elements from the specified ring
     *
     * @param degree polynomial degree
     * @param ring   the ring
     * @param rnd    random source
     * @return array of length {@code degree + 1} with elements from specified ring
     */
    public static <E> E[] randomArray(int degree, Ring<E> ring, RandomGenerator rnd) {
        E[] data = ring.createArray(degree + 1);
        for (int i = 0; i <= degree; ++i)
            data[i] = ring.randomElement(rnd);
        while (ring.isZero(data[degree]))
            data[degree] = ring.randomElement(rnd);
        return data;
    }
}
