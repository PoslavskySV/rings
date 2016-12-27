package cc.r2.core.polynomial;


public final class LongArithmetics {
    private LongArithmetics() {
    }

    /** dividend / divisor */
    public static long divide(long dividend, long divisor) {
        return dividend / divisor;
    }

    /** dividend % divisor */
    public static long rem(long dividend, long divisor) {
        return dividend % divisor;
    }

    /** Delegates to {@link Math#multiplyExact(long, long)} */
    public static long multiplyExact(long x, long y) {
        return Math.multiplyExact(x, y);
    }

    /** Delegates to {@link Math#addExact(long, long)} */
    public static long addExact(long x, long y) {
        return Math.addExact(x, y);
    }

    /** Delegates to {@link Math#subtractExact(long, long)} */
    public static long subtractExact(long a, long b) {
        return Math.subtractExact(a, b);
    }

    /** Delegates to {@link Math#floorMod(long, long)} */
    public static long floorMod(long num, long modulus) {
        return Math.floorMod(num, modulus);
    }

    /**
     * Returns the product of its arguments modulo {@code modulus},
     * throwing an exception if long overflow occurs at intermediate step.
     *
     * @param x       the first value
     * @param y       the second value
     * @param modulus modulus
     * @return the result
     * @throws ArithmeticException if long overflow occurs at intermediate step
     */
    public static long multiplyMod(long x, long y, long modulus) {
        return floorMod(multiplyExact(floorMod(x, modulus), floorMod(y, modulus)), modulus);
    }

    /**
     * Returns the product of its arguments modulo {@code modulus},
     * throwing an exception if long overflow occurs at intermediate step.
     *
     * @param x       the first value
     * @param y       the second value
     * @param z       the third value
     * @param modulus modulus
     * @return the result
     * @throws ArithmeticException if long overflow occurs at intermediate step
     */
    public static long multiplyMod(long x, long y, long z, long modulus) {
        return multiplyMod(multiplyMod(x, y, modulus), z, modulus);
    }

    /**
     * Returns the sum of its arguments modulo {@code modulus},
     * throwing an exception if long overflow occurs at intermediate step.
     *
     * @param x       the first value
     * @param y       the second value
     * @param modulus modulus
     * @return the result
     * @throws ArithmeticException if long overflow occurs at intermediate step
     */
    public static long addMod(long x, long y, long modulus) {
        return floorMod(addExact(floorMod(x, modulus), floorMod(y, modulus)), modulus);
    }

    /**
     * Returns the sum of its arguments modulo {@code modulus},
     * throwing an exception if long overflow occurs at intermediate step.
     *
     * @param x       the first value
     * @param y       the second value
     * @param z       the third value
     * @param modulus modulus
     * @return the result
     * @throws ArithmeticException if long overflow occurs at intermediate step
     */
    public static long addMod(long x, long y, long z, long modulus) {
        return addMod(addMod(x, y, modulus), z, modulus);
    }

    /**
     * Returns the difference of its arguments modulo {@code modulus},
     * throwing an exception if long overflow occurs at intermediate step.
     *
     * @param x       the first value
     * @param y       the second value
     * @param modulus modulus
     * @return the result
     * @throws ArithmeticException if long overflow occurs at intermediate step
     */
    public static long subtractMod(long x, long y, long modulus) {
        return floorMod(subtractExact(floorMod(x, modulus), floorMod(y, modulus)), modulus);
    }

    /**
     * Returns {@code a/b mod p}
     *
     * @param dividend a long
     * @param divider  a long
     * @param modulus  mosulus
     * @return {@code a/b mod p}
     * @throws IllegalArgumentException if {@code b} and {@code p} are not coprime
     */
    public static long divide(long dividend, long divider, long modulus) {
        return floorMod(multiplyExact(dividend, modInverse(divider, modulus)), modulus);
    }

    /**
     * Returns {@code value mod modulus} in the symmetric representation ({@code -modulus/2 <= result <= modulus/2})
     *
     * @param value   a long
     * @param modulus modulus
     * @return {@code value mod modulus} in the symmetric representation ({@code -modulus/2 <= result <= modulus/2})
     */
    public static long symMod(long value, long modulus) {
        if (modulus < 0)
            throw new IllegalArgumentException("Negative modulus");
        value = floorMod(value, modulus);
        return value <= modulus / 2 ? value : value - modulus;
    }

    /**
     * Returns {@code base} in a power of {@code e} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code e} (non negative)
     * @throws ArithmeticException if long overflow occurs at some intermediate step
     */
    public static long powExact(final long base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        long result = 1L;
        long k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = multiplyExact(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiplyExact(k2p, k2p);
        }
    }

    /**
     * Returns a solution of congruence {@code a * x = 1 mod p}
     *
     * @param num     base
     * @param modulus modulus
     * @return {@code a^(-1) mod p}
     * @throws IllegalArgumentException {@code a} and {@code modulus} are not coprime
     */
    public static long modInverse(long num, long modulus) {
        if (num < 0)
            num = floorMod(num, modulus);

        long s = 0, old_s = 1;
        long r = modulus, old_r = num;

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
        if (old_r != 1)
            throw new IllegalArgumentException("Not invertible: val = " + num + ", modulus = " + modulus + ", old_r = " + old_r);
        return floorMod(old_s, modulus);
    }

    /**
     * Returns {@code base} in a power of {@code e} (may be negative) modulus {@code p}
     *
     * @param base     base
     * @param exponent exponent (may be negative, if base and modulus are coprime)
     * @param modulus  modulus
     * @return {@code base} in a power of {@code e} (may be negative) modulus {@code p}
     * @throws IllegalArgumentException if e < 0 and {@code base} and {@code modulus} are not coprime
     */
    public static long modPow(long base, long exponent, long modulus) {
        if (exponent < 0) {
            base = modInverse(base, modulus);
            exponent = -exponent;
        }

        long result = 1L;
        long k2p = base;
        while (true) {
            if ((exponent&1) != 0)
                result = floorMod(multiplyExact(result, k2p), modulus);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = floorMod(multiplyExact(k2p, k2p), modulus);
        }
    }

    /**
     * Returns the greatest common divisor of two longs
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
     */
    public static long lcm(long a, long b) {
        if (a == 0 || b == 0)
            return 0;
        return a / gcd(a, b) * b;
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

    /**
     * Returns true if multiplication x*y does not cause long overflow, Copy-paste of {@link Math#multiplyExact(long, long)}
     *
     * @param x a long
     * @param y a long
     * @return true if multiplication x*y does not cause long overflow, false otherwise
     */
    public static boolean safeMultiply(long x, long y) {
        long r = x * y;
        long ax = Math.abs(x);
        long ay = Math.abs(y);
        // Some bits greater than 2^31 that might cause overflow
        // Check the result using the divide operator
        // and check for the special case of Long.MIN_VALUE * -1
        return !((((ax|ay) >>> 31 != 0)) && (((y != 0) && (r / y != x)) || (x == Long.MIN_VALUE && y == -1)));
    }

    /**
     * Returns true if summation x+y does not cause long overflow. Copy-paste of {@link Math#addExact(long, long)}
     *
     * @param x a long
     * @param y a long
     * @return true if summation x+y does not cause long overflow, false otherwise
     */
    public static boolean safeAdd(long x, long y) {
        long r = x + y;
        // HD 2-12 Overflow iff both arguments have the opposite sign of the result
        return !(((x^r)&(y^r)) < 0);
    }

    /**
     * Returns true if subtraction x-y does not cause long overflow. Copy-paste of {@link Math#subtractExact(long, long)}
     *
     * @param x a long
     * @param y a long
     * @return true if subtraction x-y does not cause long overflow, false otherwise
     */
    public static boolean safeSubtract(long x, long y) {
        long r = x - y;
        // HD 2-12 Overflow iff the arguments have different signs and
        // the sign of the result is different than the sign of x
        return !(((x^y)&(x^r)) < 0);
    }

    /**
     * Returns true if exponentiation base^e does not cause long overflow. Copy-paste of {@link #pow(long, long)}
     *
     * @param base base
     * @param e    exponent (non negative)
     * @return true if exponentiation base^e does not cause long overflow, false otherwise
     */
    public static boolean safePow(final long base, long e) {
        if (e < 0) return false;

        long result = 1L;
        long k2p = base;
        while (e != 0) {
            if ((e&1) != 0) {
                if (!safeMultiply(result, k2p))
                    return false;
                result *= k2p;
            }
            if ((e&1) == 0 && !safeMultiply(k2p, k2p))
                return false;
            k2p *= k2p;
            e = e >> 1;
        }
        return true;
    }
}
