package cc.r2.core.number;

import static java.lang.Math.*;

/**
 * Created by poslavsky on 06/12/2016.
 */
public final class ChineseRemainders {
    public static BigInteger CRT(final BigInteger[] coprimes, final BigInteger[] remainders) {
        if (coprimes.length != remainders.length)
            throw new IllegalArgumentException();
        BigInteger m = coprimes[0];
        for (int i = 1; i < coprimes.length; i++) {
            if (coprimes[i].signum() <= 0)
                throw new RuntimeException("Negative CRT input: " + coprimes[i]);
            m = coprimes[i].multiply(m);
        }

        BigInteger result = BigInteger.ZERO;
        for (int i = 0; i < coprimes.length; i++) {
            BigInteger mi = m.divide(coprimes[i]);
            BigInteger eea = bezout0(mi, coprimes[i]);
            result = result.add(mi.multiply(eea.multiply(remainders[i]).mod(coprimes[i])));
        }
        return result;
    }

    private static BigInteger bezout0(BigInteger a, BigInteger b) {
        BigInteger s = BigInteger.ZERO, old_s = BigInteger.ONE;
        BigInteger r = b, old_r = a;

        BigInteger q;
        BigInteger tmp;
        while (!r.isZero()) {
            q = old_r.divide(r);

            tmp = old_r;
            old_r = r;
            r = tmp.subtract(q.multiply(r));

            tmp = old_s;
            old_s = s;
            s = tmp.subtract(q.multiply(s));
        }
        assert old_r.isOne();
        return old_s;
    }


    /**
     * Runs Chinese Remainders algorithm
     *
     * @param prime1     #1 prime
     * @param prime2     #2 prime
     * @param remainder1 #1 remainder
     * @param remainder2 #2 remainder
     * @return the result
     */
    public static long ChineseRemainders(long prime1, long prime2,
                                         long remainder1, long remainder2) {
        if (prime1 <= 0 || prime2 <= 0)
            throw new RuntimeException("Negative CRT input: " + prime1 + " " + prime2);

        long modulus = multiplyExact(prime1, prime2);
        long result = 0;

        result = floorMod(addExact(result,
                floorMod(multiplyExact(prime2,
                        floorMod(multiplyExact(bezout0(prime2, prime1), remainder1), prime1)), modulus)), modulus);
        result = floorMod(addExact(result,
                floorMod(multiplyExact(prime1,
                        floorMod(multiplyExact(bezout0(prime1, prime2), remainder2), prime2)), modulus)), modulus);

        return result;
    }

    /**
     * Runs Chinese Remainders algorithm
     *
     * @param primes     list of coprime numbers
     * @param remainders remainder
     * @return the result
     */
    public static long ChineseRemainders(final long[] primes,
                                         final long[] remainders) {
        if (primes.length != remainders.length)
            throw new IllegalArgumentException();

        long modulus = primes[0];
        for (int i = 1; i < primes.length; ++i) {
            if (primes[i] <= 0)
                throw new RuntimeException("Negative CRT input: " + primes[i]);
            modulus = multiplyExact(primes[i], modulus);
        }

        long result = 0;
        for (int i = 0; i < primes.length; ++i) {
            long iModulus = modulus / primes[i];
            long bezout = bezout0(iModulus, primes[i]);
            result = floorMod(addExact(result,
                    floorMod(multiplyExact(iModulus,
                            floorMod(multiplyExact(bezout, remainders[i]), primes[i])), modulus)), modulus);
        }
        return result;
    }

    static long bezout0(long a, long b) {
        long s = 0, old_s = 1;
        long r = b, old_r = a;

        long q;
        long tmp;
        while (r != 0) {
            q = old_r / r;

            tmp = old_r;
            old_r = r;
            r = subtractExact(tmp , multiplyExact(q , r));

            tmp = old_s;
            old_s = s;
            s = subtractExact(tmp , multiplyExact(q , s));
        }
        assert old_r == 1 : "a = " + a + " b = " + b;
        return old_s;
    }
}
