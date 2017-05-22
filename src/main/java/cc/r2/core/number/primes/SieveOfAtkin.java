package cc.r2.core.number.primes;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.BitSet;

/**
 * Created by poslavsky on 14/11/2016.
 */
public final class SieveOfAtkin {
    private final int limit;
    private final BigInteger blimit;
    private final BitSet sieve;

    private SieveOfAtkin(int limit, BigInteger blimit, BitSet sieve) {
        this.limit = limit;
        this.blimit = blimit;
        this.sieve = sieve;
    }

    SieveOfAtkin toLimit(int newLimit) {
        return limit == newLimit ? this : new SieveOfAtkin(newLimit, BigInteger.valueOf(newLimit), sieve);
    }

    /**
     * Constructs Atkin sieve.
     *
     * @param limit limit
     */
    private SieveOfAtkin(int limit) {
        this(limit, BigInteger.valueOf(limit));
    }

    /**
     * Constructs Atkin sieve.
     *
     * @param limit  limit
     * @param blimit the same as limit (just for saving cost of int-BigInteger conversion.
     */
    private SieveOfAtkin(int limit, BigInteger blimit) {
        this.limit = limit;
        this.blimit = blimit;

        this.sieve = new BitSet(limit + 1);
        int limitSqrt = (int) Math.sqrt((double) limit);

        // the sieve works only for integers > 3, so
        // set these trivially to their proper values
        sieve.set(2);
        sieve.set(3);

        // loop through all possible integer values for x and y
        // up to the square root of the max prime for the sieve
        // we don't need any larger values for x or y since the
        // max value for x or y will be the square root of n
        // in the quadratics
        // the theorem showed that the quadratics will produce all
        // primes that also satisfy their wheel factorizations, so
        // we can produce the value of n from the quadratic first
        // and then filter n through the wheel quadratic
        // there may be more efficient ways to do this, but this
        // is the design in the Wikipedia article
        // loop through all integers for x and y for calculating
        // the quadratics
        for (int x = 1; x <= limitSqrt; x++) {
            for (int y = 1; y <= limitSqrt; y++) {
                // first quadratic using m = 12 and r in R1 = {r : 1, 5}
                int n = (4 * x * x) + (y * y);
                if (n <= limit && (n % 12 == 1 || n % 12 == 5))
                    sieve.flip(n);

                // second quadratic using m = 12 and r in R2 = {r : 7}
                n = (3 * x * x) + (y * y);
                if (n <= limit && (n % 12 == 7))
                    sieve.flip(n);

                // third quadratic using m = 12 and r in R3 = {r : 11}
                n = (3 * x * x) - (y * y);
                if (x > y && n <= limit && (n % 12 == 11))
                    sieve.flip(n);

                // note that R1 union R2 union R3 is the set R
                // R = {r : 1, 5, 7, 11}
                // which is all values 0 < r < 12 where r is
                // a relative prime of 12
                // Thus all primes become candidates
            }
        }

        // remove all perfect squares since the quadratic
        // wheel factorization filter removes only some of them
        for (int n = 5; n <= limitSqrt; n++) {
            if (sieve.get(n)) {
                int x = n * n;
                for (int i = x; i <= limit; i += x)
                    sieve.clear(i);
            }
        }
    }

    public boolean isPrime(int n) {
        if (n > limit)
            throw new IndexOutOfBoundsException("Out of sieve bounds.");
        return sieve.get(n);
    }

    /** Returns the last prime in this sieve */
    public int lastPrime() {
        for (int i = limit; i >= 0; --i)
            if (isPrime(i))
                return i;
        throw new IllegalStateException("No ant primes in the sieve");
    }

    public int randomPrime(RandomGenerator rnd) {
        int i;
        do {
            i = rnd.nextInt(limit);
        } while (!isPrime(i));
        return i;
    }

    public int getLimit() {
        return limit;
    }

    public BigInteger getLimitAsBigInteger() {
        return blimit;
    }

    //cached sieve
    static final SieveOfAtkin SmallPrimesSieve = new SieveOfAtkin(64 * 1024);

    public static SieveOfAtkin createSieve(int limit) {
        if (limit <= SmallPrimesSieve.limit)
            return SmallPrimesSieve.toLimit(limit);
        return new SieveOfAtkin(limit);
    }

    public static SieveOfAtkin createSieve(BigInteger limit) {
        if (limit.compareTo(SmallPrimesSieve.blimit) < 9)
            return SmallPrimesSieve.toLimit(limit.intValue());
        return new SieveOfAtkin(limit.intValueExact(), limit);
    }
}