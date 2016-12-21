package cc.r2.core.number.primes;

import cc.r2.core.number.BigInteger;

import java.util.Arrays;

public final class PrimesIterator {
    private final int[] smallPrimes = SmallPrimes.SmallPrimes12;
    private long pointer;
    private final long from;
    private static final int largeSieveLimit = 16777216;
    private final SieveOfAtkin smallSieve = SieveOfAtkin.SmallPrimesSieve;
    private SieveOfAtkin largeSieve = null;

    public PrimesIterator() {from = 0;}

    public PrimesIterator(long from) {
        this.from = from;
        if (from < smallPrimes[smallPrimes.length - 1]) {
            pointer = Arrays.binarySearch(smallPrimes, (int) from);
            if (pointer < 0) pointer = ~pointer;
        } else pointer = from;
    }

    public long take() {
        if (pointer < smallPrimes.length)
            return smallPrimes[(int) (pointer++)];

        for (; pointer < smallSieve.getLimit(); )
            if (smallSieve.isPrime((int) (pointer)))
                return pointer - 1;

        if (pointer < largeSieveLimit) {
            if (largeSieve == null)
                largeSieve = SieveOfAtkin.createSieve(largeSieveLimit);

            for (; pointer < largeSieve.getLimit(); )
                if (largeSieve.isPrime((int) (pointer++)))
                    return pointer - 1;
        }

        if (pointer < Integer.MAX_VALUE - 1)
            return pointer = SmallPrimes.nextPrime((int) (pointer + 1));

        if (pointer < Long.MAX_VALUE - 1)
            try {
                return pointer = BigPrimes.nextPrime(BigInteger.valueOf(pointer + 1)).longValueExact();
            } catch (ArithmeticException e) {
                return -1;
            }

        return -1;
    }
}
