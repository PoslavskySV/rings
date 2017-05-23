package cc.r2.core.number;

import cc.redberry.libdivide4j.FastDivision;

import static java.lang.Math.*;

/**
 * Created by poslavsky on 06/12/2016.
 */
public final class ChineseRemainders {

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
        long result;

        //(result + (prime2 * ((bezout0(prime2, prime1) * remainder1) % prime1)) % modulus) % modulus
        result = floorMod(multiplyExact(prime2,
                floorMod(multiplyExact(bezout0(prime2, prime1), remainder1), prime1)), modulus);
        result = floorMod(addExact(result,
                floorMod(multiplyExact(prime1,
                        floorMod(multiplyExact(bezout0(prime1, prime2), remainder2), prime2)), modulus)), modulus);

        return result;
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
    public static BigInteger ChineseRemainders(BigInteger prime1, BigInteger prime2,
                                               BigInteger remainder1, BigInteger remainder2) {
        if (prime1.signum() <= 0 || prime2.signum() <= 0)
            throw new RuntimeException("Negative CRT input: " + prime1 + " " + prime2);

        BigInteger modulus = prime1.multiply(prime2);
        BigInteger result = BigInteger.ZERO;

        //(result + (prime2 * ((bezout0(prime2, prime1) * remainder1) % prime1)) ) % modulus
        result = result.add(prime2.multiply(bezout0(prime2, prime1).multiply(remainder1).mod(prime1))).mod(modulus);
        result = result.add(prime1.multiply(bezout0(prime1, prime2).multiply(remainder2).mod(prime2))).mod(modulus);

        return result;
    }

    /**
     * Runs Chinese Remainders algorithm
     *
     * @param magic      magic structure for fast modular arithmetics ({@link #createMagic(long, long)}
     * @param remainder1 #1 remainder
     * @param remainder2 #2 remainder
     * @return the result
     */
    public static long ChineseRemainders(ChineseRemaindersMagic magic,
                                         long remainder1, long remainder2) {
        long result;

        //(result + (prime2 * ((bezout0(prime2, prime1) * remainder1) % prime1)) % modulus) % modulus
        result = magic.mulModModulus(magic.prime2, magic.mulModPrime1(magic.bezoutPrime2Prime1, remainder1));
        result = magic.addMod(result,
                magic.mulModModulus(magic.prime1, magic.mulModPrime2(magic.bezoutPrime1Prime2, remainder2)));

        return result;
    }

    /** Magic for fast repeated Chinese Remainders */
    public static ChineseRemaindersMagic createMagic(long prime1, long prime2) {
        long modulus = multiplyExact(prime1, prime2);
        return new ChineseRemaindersMagic(
                prime1, prime2, modulus,
                FastDivision.magicUnsigned(prime1), FastDivision.magicUnsigned(prime2), FastDivision.magicUnsigned(modulus),
                prime1 <= Integer.MAX_VALUE, prime2 <= Integer.MAX_VALUE, modulus <= Integer.MAX_VALUE,
                bezout0(prime1, prime2), bezout0(prime2, prime1)
        );
    }

    public static class ChineseRemaindersMagic {
        final long prime1, prime2, modulus;
        final FastDivision.Magic magic1, magic2, magicModulus;
        final boolean prime1IsInt, prime2IsInt, modulusIsInt;
        final FastDivision.Magic magic32Prime1, magic32Prime2, magic32Modulus;
        final long bezoutPrime1Prime2, bezoutPrime2Prime1;

        ChineseRemaindersMagic(long prime1, long prime2, long modulus, FastDivision.Magic magic1, FastDivision.Magic magic2, FastDivision.Magic magicModulus, boolean prime1IsInt, boolean prime2IsInt, boolean modulusIsInt, long bezoutPrime1Prime2, long bezoutPrime2Prime1) {
            this.prime1 = prime1;
            this.prime2 = prime2;
            this.modulus = modulus;
            this.magic1 = magic1;
            this.magic2 = magic2;
            this.magicModulus = magicModulus;
            this.prime1IsInt = prime1IsInt;
            this.prime2IsInt = prime2IsInt;
            this.modulusIsInt = modulusIsInt;
            this.magic32Prime1 = prime1IsInt ? null : FastDivision.magic32ForMultiplyMod(prime1);
            this.magic32Prime2 = prime2IsInt ? null : FastDivision.magic32ForMultiplyMod(prime2);
            this.magic32Modulus = modulusIsInt ? null : FastDivision.magic32ForMultiplyMod(modulus);
            this.bezoutPrime1Prime2 = Math.floorMod(bezoutPrime1Prime2, prime2);
            this.bezoutPrime2Prime1 = Math.floorMod(bezoutPrime2Prime1, prime1);
        }

        long mulModPrime1(long a, long b) {
            return prime1IsInt ? FastDivision.modUnsignedFast(a * b, magic1) : FastDivision.multiplyMod128Unsigned(a, b, prime1, magic32Prime1);
        }

        long mulModPrime2(long a, long b) {
            return prime2IsInt ? FastDivision.modUnsignedFast(a * b, magic2) : FastDivision.multiplyMod128Unsigned(a, b, prime2, magic32Prime2);
        }

        long mulModModulus(long a, long b) {
            return modulusIsInt ? FastDivision.modUnsignedFast(a * b, magicModulus) : FastDivision.multiplyMod128Unsigned(a, b, modulus, magic32Modulus);
        }

        long addMod(long a, long b) {
            long r = a + b;
            return r - modulus >= 0 ? r - modulus : r;
        }
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

    /**
     * Runs Chinese Remainders algorithm
     *
     * @param primes     list of coprime numbers
     * @param remainders remainder
     * @return the result
     */
    public static BigInteger ChineseRemainders(final BigInteger[] primes, final BigInteger[] remainders) {
        if (primes.length != remainders.length)
            throw new IllegalArgumentException();
        BigInteger m = primes[0];
        for (int i = 1; i < primes.length; i++) {
            if (primes[i].signum() <= 0)
                throw new RuntimeException("Negative CRT input: " + primes[i]);
            m = primes[i].multiply(m);
        }

        BigInteger result = BigInteger.ZERO;
        for (int i = 0; i < primes.length; i++) {
            BigInteger mi = m.divide(primes[i]);
            BigInteger eea = bezout0(mi, primes[i]);
            result = result.add(mi.multiply(eea.multiply(remainders[i]).mod(primes[i])));
        }
        return result;
    }

    private static long bezout0(long a, long b) {
        long s = 0, old_s = 1;
        long r = b, old_r = a;

        long q;
        long tmp;
        while (r != 0) {
            q = old_r / r;

            tmp = old_r;
            old_r = r;
            r = subtractExact(tmp, multiplyExact(q, r));

            tmp = old_s;
            old_s = s;
            s = subtractExact(tmp, multiplyExact(q, s));
        }
        assert old_r == 1 : "a = " + a + " b = " + b;
        return old_s;
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
}
