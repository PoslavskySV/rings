package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.redberry.libdivide4j.FastDivision;

/**
 * Helper methods for arithmetics with {@code longs}.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class LongArithmetics {
    /** Max supported modulus bits */
    public static final int MAX_SUPPORTED_MODULUS_BITS = 62;
    /** Max supported modulus */
    public static final long MAX_SUPPORTED_MODULUS = (1L << MAX_SUPPORTED_MODULUS_BITS) - 1L;
    public static final BigInteger b_MAX_SUPPORTED_MODULUS = BigInteger.valueOf(MAX_SUPPORTED_MODULUS);

    private LongArithmetics() {}

    /**
     * Returns true if {@code val} fits into 32-bit machine word (unsigned) and false otherwise
     *
     * @param val the value
     * @return true if {@code val} fits into 32-bit machine word (unsigned) and false otherwise
     */
    public static boolean fits32bitWord(long val) {
        return Long.compareUnsigned(val, (1L << 32) - 1) <= 0;
    }

    /**
     * Returns true if {@code val} fits into 32-bit machine word (unsigned) and false otherwise
     *
     * @param val the value
     * @return true if {@code val} fits into 32-bit machine word (unsigned) and false otherwise
     */
    public static boolean fits31bitWord(long val) {
        return Long.compareUnsigned(val, (1L << 31) - 1) <= 0;
    }

    /**
     * Delegates to {@link Math#multiplyExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long safeMultiply(long x, long y) {
        return Math.multiplyExact(x, y);
    }

    /**
     * Delegates to {@link Math#multiplyExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long safeMultiply(long x, long y, long z) {
        return Math.multiplyExact(Math.multiplyExact(x, y), z);
    }

    /**
     * Delegates to {@link Math#addExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long safeAdd(long x, long y) {
        return Math.addExact(x, y);
    }

    /**
     * Delegates to {@link Math#subtractExact(long, long)}
     *
     * @throws ArithmeticException if the result overflows a long
     **/
    public static long safeSubtract(long a, long b) {
        return Math.subtractExact(a, b);
    }

    public static boolean isOverflowMultiply(long x, long y) {
        long r = x * y;
        long ax = Math.abs(x);
        long ay = Math.abs(y);
        if (((ax|ay) >>> 31 != 0)) {
            // Some bits greater than 2^31 that might cause overflow
            // Check the result using the divide operator
            // and check for the special case of Long.MIN_VALUE * -1
            if (((y != 0) && (r / y != x)) ||
                    (x == Long.MIN_VALUE && y == -1)) {
                return true;
            }
        }
        return false;
    }

    public static boolean isOverflowAdd(long x, long y) {
        long r = x + y;
        // HD 2-12 Overflow iff both arguments have the opposite sign of the result
        if (((x^r)&(y^r)) < 0) {
            return true;
        }
        return false;
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e}
     * @throws ArithmeticException if the result overflows a long
     */
    public static long safePow(final long base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        long result = 1L;
        long k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = safeMultiply(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = safeMultiply(k2p, k2p);
        }
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e}
     */
    public static long unsafePow(final long base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        long result = 1L;
        long k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result *= k2p;
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p *= k2p;
        }
    }

    /**
     * Returns the greatest common divisor of two longs.
     *
     * @param p a long
     * @param q a long
     * @return greatest common divisor of {@code a} and {@code b}
     */
    public static long gcd(final long p, final long q) {
        long u = p;
        long v = q;
        if ((u == 0) || (v == 0)) {
            if ((u == Long.MIN_VALUE) || (v == Long.MIN_VALUE))
                throw new IllegalArgumentException("long overflow");
            return Math.abs(u) + Math.abs(v);
        }
        // keep u and v negative, as negative integers range down to
        // -2^63, while positive numbers can only be as large as 2^63-1
        // (i.e. we can't necessarily negate a negative number without
        // overflow)
        /* assert u!=0 && v!=0; */
        if (u > 0) {
            u = -u;
        } // make u negative
        if (v > 0) {
            v = -v;
        } // make v negative
        // B1. [Find power of 2]
        int k = 0;
        while ((u&1) == 0 && (v&1) == 0 && k < 63) { // while u and v are
            // both even...
            u /= 2;
            v /= 2;
            k++; // cast out twos.
        }
        if (k == 63) {
            throw new IllegalArgumentException("Overflow");
        }
        // B2. Initialize: u and v have been divided by 2^k and at least
        // one is odd.
        long t = ((u&1) == 1) ? v : -(u / 2)/* B3 */;
        // t negative: u was odd, v may be even (t replaces v)
        // t positive: u was even, v is odd (t replaces u)
        do {
            /* assert u<0 && v<0; */
            // B4/B3: cast out twos from t.
            while ((t&1) == 0) { // while t is even..
                t /= 2; // cast out twos
            }
            // B5 [reset max(u,v)]
            if (t > 0) {
                u = -t;
            } else {
                v = t;
            }
            // B6/B3. at this point both u and v should be odd.
            t = (v - u) / 2;
            // |u| larger: t positive (replace u)
            // |v| larger: t negative (replace v)
        } while (t != 0);
        return -u * (1L << k); // gcd is u*2^k
    }

    /**
     * Runs extended Euclidean algorithm to compute {@code [gcd(a,b), x, y]} such that {@code x * a + y * b = gcd(a, b)}
     *
     * @param a a long
     * @param b a long
     * @return array of {@code [gcd(a,b), x, y]} such that {@code x * a + y * b = gcd(a, b)}
     */
    public static long[] gcdExtended(long a, long b) {
        long s = 0, old_s = 1;
        long t = 1, old_t = 0;
        long r = b, old_r = a;

        long q;
        long tmp;
        while (r != 0) {
            q = old_r / r;

            tmp = old_r;
            old_r = r;
            r = tmp - q * r;

            tmp = old_s;
            old_s = s;
            s = tmp - q * s;

            tmp = old_t;
            old_t = t;
            t = tmp - q * t;
        }
        assert old_r == a * old_s + b * old_t;
        return new long[]{old_r, old_s, old_t};
    }

    /**
     * Returns the least common multiple of two longs
     *
     * @param a a long
     * @param b a long
     * @return least common multiple of {@code a} and {@code b}
     * @throws ArithmeticException if the result overflows a long
     */
    public static long lcm(long a, long b) {
        if (a == 0 || b == 0)
            return 0;
        return safeMultiply(a / gcd(a, b), b);
    }

    /**
     * Returns the greatest common an array of longs
     *
     * @param integers array of longs
     * @param from     from position (inclusive)
     * @param to       to position (exclusive)
     * @return greatest common divisor of array
     */
    public static long gcd(final long[] integers, int from, int to) {
        if (integers.length < 2)
            throw new IllegalArgumentException();
        long gcd = gcd(integers[from], integers[from + 1]);
        if (gcd == 1)
            return gcd;
        for (int i = from + 2; i < to; i++) {
            gcd = gcd(integers[i], gcd);
            if (gcd == 1)
                return gcd;
        }
        return gcd;
    }

    /**
     * Returns the greatest common an array of longs
     *
     * @param integers array of longs
     * @return greatest common divisor of array
     */
    public static long gcd(final long... integers) {
        return gcd(integers, 0, integers.length);
    }

    /** Delegates to {@link Math#floorMod(long, long)} */
    public static long mod(long num, long modulus) {
        if (num < 0)
            num += modulus; //<- may help
        return (num >= modulus || num < 0) ? Math.floorMod(num, modulus) : num;
    }

    /**
     * Returns {@code value mod modulus} in the symmetric representation ({@code -modulus/2 <= result <= modulus/2})
     *
     * @param value   a long
     * @param modulus modulus
     * @return {@code value mod modulus} in the symmetric representation ({@code -modulus/2 <= result <= modulus/2})
     */
    public static long symMod(long value, long modulus) {
        value = mod(value, modulus);
        return value <= modulus / 2 ? value : value - modulus;
    }

    /**
     * Returns a solution of congruence {@code num * x = 1 mod modulus}
     *
     * @param num     base
     * @param modulus modulus
     * @return {@code a^(-1) mod p}
     * @throws IllegalArgumentException {@code a} and {@code modulus} are not coprime
     */
    public static long modInverse(long num, long modulus) {
        if (num == 1)
            return num;
        if (num < 0)
            num = mod(num, modulus);

        long s = 0, old_s = 1;
        long r = modulus, old_r = num;

        long q;
        long tmp;
        while (r != 0) {
            q = old_r / r;

            tmp = old_r;
            old_r = r;
            r = tmp - q * r;

            tmp = old_s;
            old_s = s;
            s = tmp - q * s;
        }

        if (old_r != 1)
            throw new IllegalArgumentException(
                    String.format("modInverse(%s, %s) : not invertible (old_r = %s)", num, modulus, old_r));
        return mod(old_s, modulus);
    }

    /**
     * Casts {@code long} to signed {@code int} throwing exception in case of overflow.
     *
     * @param value the long
     * @return int value
     * @throws ArithmeticException if the result overflows a long
     */
    public static int safeToInt(long value) {
        if ((int) value != value) {
            throw new ArithmeticException("integer overflow: " + value);
        }
        return (int) value;
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code modulus}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param modulus  the modulus
     * @return {@code base} in a power of {@code e}
     */
    public static long powMod(final long base, long exponent, long modulus) {
        return powModSigned(base, exponent, FastDivision.magicSigned(modulus));
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code magic.modulus}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param magic    magic modulus
     * @return {@code base} in a power of {@code e}
     */
    public static long powModSigned(final long base, long exponent, FastDivision.Magic magic) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return 1;

        long result = 1;
        long k2p = FastDivision.modSignedFast(base, magic);
        for (; ; ) {
            if ((exponent&1) != 0)
                result = FastDivision.modSignedFast(result * k2p, magic);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = FastDivision.modSignedFast(k2p * k2p, magic);
        }
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code magic.modulus}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @param magic    magic modulus
     * @return {@code base} in a power of {@code e}
     */
    public static long powModUnsigned(final long base, long exponent, FastDivision.Magic magic) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return 1;

        long result = 1;
        long k2p = FastDivision.modUnsignedFast(base, magic);
        for (; ; ) {
            if ((exponent&1) != 0)
                result = FastDivision.modUnsignedFast(result * k2p, magic);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = FastDivision.modUnsignedFast(k2p * k2p, magic);
        }
    }

    /**
     * Tests whether {@code n} is a perfect power {@code n == a^b} and returns {@code {a, b}} if so and {@code null}
     * otherwise
     *
     * @param n the number
     * @return array {@code {a, b}} so that {@code n = a^b} or {@code null} is {@code n} is not a perfect power
     */
    public static long[] perfectPowerDecomposition(long n) {
        if (n < 0) {
            long[] ipp = perfectPowerDecomposition(-n);
            if (ipp == null)
                return null;
            if (ipp[1] % 2 == 0)
                return null;
            ipp[0] = -ipp[0];
            return ipp;
        }

        if ((-n&n) == n)
            return new long[]{2, 63 - Long.numberOfLeadingZeros(n)};

        int lgn = 1 + (64 - Long.numberOfLeadingZeros(n));
        for (int b = 2; b < lgn; b++) {
            //b lg a = lg n
            long lowa = 1L;
            long higha = 1L << (lgn / b + 1);
            while (lowa < higha - 1) {
                long mida = (lowa + higha) >> 1;
                long ab = LongArithmetics.unsafePow(mida, b);
                if (ab > n)
                    higha = mida;
                else if (ab < n)
                    lowa = mida;
                else {
                    long[] ipp = perfectPowerDecomposition(mida);
                    if (ipp != null)
                        return new long[]{ipp[0], ipp[1] * b};
                    else
                        return new long[]{mida, b};
                }
            }
        }
        return null;
    }
}
