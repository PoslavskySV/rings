package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;

/**
 * Created by poslavsky on 26/01/2017.
 */
public final class LongModularArithmetics {
    private LongModularArithmetics() {}

    /** Max supported modulus */
    public static final long MAX_SUPPORTED_MODULUS = (1L << 62) - 1;

    /**
     * Returns highest 64 bits of (signed) long multiplication.
     *
     * @param x the number
     * @param y the number
     * @return highest 64 bits of (signed) long multiplication.
     */
    public static long multiplyHighSigned(long x, long y) {
        long x_high = x >> 32;
        long x_low = x&0xFFFFFFFFL;
        long y_high = y >> 32;
        long y_low = y&0xFFFFFFFFL;

        long z2 = x_low * y_low;
        long t = x_high * y_low + (z2 >>> 32);
        long z1 = t&0xFFFFFFFFL;
        long z0 = t >> 32;
        z1 += x_low * y_high;
        return x_high * y_high + z0 + (z1 >> 32);
    }

    /**
     * Returns highest 64 bits of (unsigned) long multiplication.
     *
     * @param x the number
     * @param y the number
     * @return highest 64 bits of (unsigned) long multiplication.
     */
    public static long multiplyHighUnsigned(long x, long y) {
        long x_high = x >>> 32;
        long y_high = y >>> 32;
        long x_low = x&0xFFFFFFFFL;
        long y_low = y&0xFFFFFFFFL;

        long z2 = x_low * y_low;
        long t = x_high * y_low + (z2 >>> 32);
        long z1 = t&0xFFFFFFFFL;
        long z0 = t >>> 32;
        z1 += x_low * y_high;
        return x_high * y_high + z0 + (z1 >>> 32);
    }

    /**
     * Returns lowest 64 bits of either signed or unsigned long multiplication.
     *
     * @param x the number
     * @param y the number
     * @return lowest 64 bits of long multiplication.
     */
    public static long multiplyLow(long x, long y) {
        long n = x * y;
        return n;
    }

    /**
     * Return's quotient and remainder of 128 bit integer division by 64 bit integer.
     * <p>
     * Code taken from Hacker's Delight:
     * http://www.hackersdelight.org/HDcode/divlu.c.
     *
     * @param u1 highest 64 dividend bits
     * @param u0 lowest 64 dividend bits
     * @param v  the divider
     * @return {quotient, remainder}
     */
    public static long[] divideAndRemainder128(long u1, long u0, long v) {
        long b = (1L << 32); // Number base (16 bits).
        long
                un1, un0,           // Norm. dividend LSD's.
                vn1, vn0,           // Norm. divisor digits.
                q1, q0,             // Quotient digits.
                un64, un21, un10,   // Dividend digit pairs.
                rhat;               // A remainder.
        int s;              // Shift amount for norm.

        if (u1 >= v)                          // If overflow, set rem.
            return new long[]{-1L, -1L};      // possible quotient.


        // count leading zeros
        s = Long.numberOfLeadingZeros(v); // 0 <= s <= 63.
        if (s > 0) {
            v = v << s;         // Normalize divisor.
            un64 = (u1 << s)|((u0 >>> (64 - s))&(-s >> 31));
            un10 = u0 << s;     // Shift dividend left.
        } else {
            // Avoid undefined behavior.
            un64 = u1|u0;
            un10 = u0;
        }

        vn1 = v >>> 32;            // Break divisor up into
        vn0 = v&0xFFFFFFFFL;     // two 32-bit digits.

        un1 = un10 >>> 32;         // Break right half of
        un0 = un10&0xFFFFFFFFL;  // dividend into two digits.

        q1 = Long.divideUnsigned(un64, vn1);            // Compute the first
        rhat = un64 - q1 * vn1;     // quotient digit, q1.
        while (true) {
            if (Long.compareUnsigned(q1, b) >= 0 || Long.compareUnsigned(q1 * vn0, b * rhat + un1) > 0) { //if (q1 >= b || q1 * vn0 > b * rhat + un1) {
                q1 = q1 - 1;
                rhat = rhat + vn1;
                if (Long.compareUnsigned(rhat, b) < 0)
                    continue;
            }
            break;
        }

        un21 = un64 * b + un1 - q1 * v;  // Multiply and subtract.

        q0 = Long.divideUnsigned(un21, vn1);            // Compute the second
        rhat = un21 - q0 * vn1;     // quotient digit, q0.
        while (true) {
            if (Long.compareUnsigned(q0, b) >= 0 || Long.compareUnsigned(q0 * vn0, b * rhat + un0) > 0) {
                q0 = q0 - 1;
                rhat = rhat + vn1;
                if (Long.compareUnsigned(rhat, b) < 0)
                    continue;
            }
            break;
        }
        long r = (un21 * b + un0 - q0 * v) >>> s;    // return it.
        return new long[]{q1 * b + q0, r};
    }

    /**
     * Computes magic for fast unsigned integer division.
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param d the divider (must be positive)
     * @return the magic
     */
    public static MagicDivider magicUnsigned(long d) {
        return magicUnsigned(d, false);
    }

    /**
     * Computes magic for fast unsigned integer division.
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param d the divider (must be positive)
     * @return the magic
     */
    public static MagicDivider magicUnsigned(long d, boolean branchfree) {
        // 1 is not supported with branchfree algorithm
        assert (!branchfree || d != 1);

        long resultMagic;
        int resultMore;
        int floor_log_2_d = 63 - Long.numberOfLeadingZeros(d);
        if ((d&(d - 1)) == 0) {
            // Power of 2
            if (!branchfree) {
                resultMagic = 0;
                resultMore = floor_log_2_d|0x80;
            } else {
                // We want a magic number of 2**64 and a shift of floor_log_2_d
                // but one of the shifts is taken up by LIBDIVIDE_ADD_MARKER, so we
                // subtract 1 from the shift
                resultMagic = 0;
                resultMore = (floor_log_2_d - 1)|0x40;
            }
        } else {
            long proposed_m, rem;
            int more;

            long[] tmp = divideAndRemainder128(1L << floor_log_2_d, 0, d); // == (1 << (64 + floor_log_2_d)) / d
            proposed_m = tmp[0];
            rem = tmp[1];

//            assert (rem > 0 && rem < d);
            long e = d - rem;

            // This power works if e < 2**floor_log_2_d.
            if (!branchfree && e < (1L << floor_log_2_d)) {
                // This power works
                more = floor_log_2_d;
            } else {
                // We have to use the general 65-bit algorithm.  We need to compute
                // (2**power) / d. However, we already have (2**(power-1))/d and
                // its remainder. By doubling both, and then correcting the
                // remainder, we can compute the larger division.
                // don't care about overflow here - in fact, we expect it
                proposed_m += proposed_m;
                long twice_rem = rem + rem;
                if (twice_rem >= d || twice_rem < rem) proposed_m += 1;
                more = floor_log_2_d|0x40;
            }
            resultMagic = 1 + proposed_m;
            resultMore = more;
            // result.more's shift should in general be ceil_log_2_d. But if we
            // used the smaller power, we subtract one from the shift because we're
            // using the smaller power. If we're using the larger power, we
            // subtract one from the shift because it's taken care of by the add
            // indicator. So floor_log_2_d happens to be correct in both cases,
            // which is why we do it outside of the if statement.
        }
        return new MagicDivider(resultMagic, resultMore, d);
    }

    /**
     * Returns unsigned {@code dividend / divider} using fast integer division
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider }
     */
    public static long divideUnsignedFast(long dividend, MagicDivider divider) {
        int more = divider.more;
        if ((more&0x80) != 0) {
            return dividend >>> (more&0x3F);
        } else {
            long q = multiplyHighUnsigned(divider.magic, dividend);
            if ((more&0x40) != 0) {
                long t = ((dividend - q) >>> 1) + q;
                return t >>> (more&0x3F);
            } else {
                return q >>> more; // all upper bits are 0 - don't need to mask them off
            }
        }
    }


    /**
     * Computes magic for fast signed integer division.
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param d the divider (must be positive)
     * @return the magic
     */
    public static MagicDivider magicSigned(long d) {
        return magicSigned(d, true);
    }

    /**
     * Computes magic for fast integer division.
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param d          the divider (must be positive)
     * @param branchfree use branch free computation
     * @return the magic
     */
    public static MagicDivider magicSigned(long d, boolean branchfree) {
        if (d < 0)
            throw new IllegalArgumentException("Negative divider.");

        assert (!branchfree || (d != 1 && d != -1));

        long resultMagic;
        int resultMore;
        // If d is a power of 2, or negative a power of 2, we have to use a shift.
        // This is especially important because the magic algorithm fails for -1.
        // To check if d is a power of 2 or its inverse, it suffices to check
        // whether its absolute value has exactly one bit set.  This works even for
        // INT_MIN, because abs(INT_MIN) == INT_MIN, and INT_MIN has one bit set
        // and is a power of 2.
        long ud = d;
        long absD = (d < 0 ? -ud : ud); // gcc optimizes this to the fast abs trick
        int floor_log_2_d = 63 - Long.numberOfLeadingZeros(absD);
        // check if exactly one bit is set,
        // don't care if absD is 0 since that's divide by zero
        if ((absD&(absD - 1)) == 0) {
            // Branchfree and non-branchfree cases are the same
            resultMagic = 0;
            resultMore = floor_log_2_d|(d < 0 ? 0x80 : 0);
        } else {
            // the dividend here is 2**(floor_log_2_d + 63), so the low 64 bit word
            // is 0 and the high word is floor_log_2_d - 1
            int more;
            long rem, proposed_m;
            long[] tmp = divideAndRemainder128(1L << (floor_log_2_d - 1), 0, absD);
            proposed_m = tmp[0];
            rem = tmp[1];
            long e = absD - rem;

            // We are going to start with a power of floor_log_2_d - 1.
            // This works if works if e < 2**floor_log_2_d.
            if (!branchfree && e < (1L << floor_log_2_d)) {
                // This power works
                more = floor_log_2_d - 1;
            } else {
                // We need to go one higher. This should not make proposed_m
                // overflow, but it will make it negative when interpreted as an
                // int32_t.
                proposed_m += proposed_m;
                long twice_rem = rem + rem;
                if (Long.compareUnsigned(twice_rem, absD) >= 0 || Long.compareUnsigned(twice_rem, rem) < 0)
                    proposed_m += 1;
                // note that we only set the LIBDIVIDE_NEGATIVE_DIVISOR bit if we
                // also set ADD_MARKER this is an annoying optimization that
                // enables algorithm #4 to avoid the mask. However we always set it
                // in the branchfree case
                more = floor_log_2_d|0x40;
            }
            proposed_m += 1;
            long magic = proposed_m;

            // Mark if we are negative
            if (d < 0) {
                more |= 0x80;
                if (!branchfree) {
                    magic = -magic;
                }
            }

            resultMore = more;
            resultMagic = magic;
        }
        return new MagicDivider(resultMagic, resultMore, d);
    }

    /**
     * Returns signed {@code dividend / divider} using fast integer division
     * <p>
     * Code taken from C libdivide library http://libdivide.com
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider }
     */
    public static long divideSignedFast(long dividend, MagicDivider divider) {
        int more = divider.more;
        long magic = divider.magic;
        if (magic == 0) { //shift path
            int shifter = more&0x3F;
            long uq = dividend + ((dividend >> 63)&((1L << shifter) - 1));
            long q = uq;
            q = q >> shifter;
            // must be arithmetic shift and then sign-extend
            long shiftMask = more >> 7;
            q = (q^shiftMask) - shiftMask;
            return q;
        } else {
            long uq = multiplyHighSigned(magic, dividend);
            if ((more&0x40) != 0) {
                // must be arithmetic shift and then sign extend
                long sign = more >> 7;
                uq += (((long) dividend^sign) - sign);
            }
            long q = (long) uq;
            q >>= more&0x3F;
            if (q < 0)
                q += 1;
            return q;
        }
    }

    /**
     * Calculates the remainder using fast integer division
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend % divider }
     */
    public static long remainderSignedFast(long dividend, MagicDivider divider) {
        long quot = divideSignedFast(dividend, divider);
        return dividend - quot * divider.divider;
    }

    /**
     * Calculates the remainder using fast integer division
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend % divider }
     */
    public static long remainderUnsignedFast(long dividend, MagicDivider divider) {
        long quot = divideUnsignedFast(dividend, divider);
        return dividend - quot * divider.divider;
    }

    /**
     * Computes floor division of the dividend by the divider using fast integer division returning (meaningful for signed operations)
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider }
     * @see Math#floorDiv(long, long)
     */
    public static long floorDivideFast(long dividend, MagicDivider divider) {
        long r = divideSignedFast(dividend, divider);
        // if the signs are different and modulo not zero, round down
        if ((dividend^divider.divider) < 0 && (r * divider.divider != dividend)) {
            r--;
        }
        return r;
    }

    /**
     * Calculates the modulus using fast integer division
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend % divider }
     */
    public static long modSignedFast(long dividend, MagicDivider divider) {
        return dividend - floorDivideFast(dividend, divider) * divider.divider;
    }

    /**
     * Calculates the modulus using fast integer division
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend % divider }
     */
    public static long modUnsignedFast(long dividend, MagicDivider divider) {
        return dividend - divideUnsignedFast(dividend, divider) * divider.divider;
    }

    /**
     * Computes magic for fast mulmod operation.
     *
     * @param divider the divider (must be positive)
     * @return the magic
     */
    public static MagicDivider magic32ForMultiplyMod(long divider) {
        long v = divider;
        int s = Long.numberOfLeadingZeros(v); // 0 <= s <= 63.
        if (s > 0)
            v = v << s;
        return magicUnsigned(v >>> 32);
    }

    /**
     * Returns unsigned {@code (a*b)%divider}
     *
     * @param a       the first multiplier
     * @param b       the second multiplier
     * @param divider the divider
     * @param magic32 magic for fast division {@link #magic32ForMultiplyMod(long)}
     * @return {@code (a*b)%divider }
     */
    public static long multiplyMod128Unsigned(long a, long b, long divider, MagicDivider magic32) {
        return multiplyMod128Unsigned0(multiplyHighUnsigned(a, b), multiplyLow(a, b), divider, magic32);
    }

    /**
     * Returns unsigned {@code (low|(high<<64))%divider}
     *
     * @param high    the highest bits
     * @param low     the lowest bits
     * @param divider the divider
     * @param magic32 magic for fast division {@link #magic32ForMultiplyMod(long)}
     * @return {@code (low|(high<<64))%divider}
     */
    public static long multiplyMod128Unsigned0(long high, long low, long divider, MagicDivider magic32) {
        long b = (1L << 32); // Number base (16 bits).
        long
                un1, un0,           // Norm. dividend LSD's.
                vn1, vn0,           // Norm. divisor digits.
                q1, q0,             // Quotient digits.
                un64, un21, un10,   // Dividend digit pairs.
                rhat;               // A remainder.
        int s;              // Shift amount for norm.

        if (high >= divider)                          // If overflow, set rem.
            throw new IllegalArgumentException();


        // count leading zeros
        s = Long.numberOfLeadingZeros(divider); // 0 <= s <= 63.
        if (s > 0) {
            divider = divider << s;         // Normalize divisor.
            un64 = (high << s)|((low >>> (64 - s))&(-s >> 31));
            un10 = low << s;     // Shift dividend left.
        } else {
            // Avoid undefined behavior.
            un64 = high|low;
            un10 = low;
        }

        vn1 = divider >>> 32;            // Break divisor up into
        vn0 = divider&0xFFFFFFFFL;     // two 32-bit digits.

        un1 = un10 >>> 32;         // Break right half of
        un0 = un10&0xFFFFFFFFL;  // dividend into two digits.

        q1 = divideUnsignedFast(un64, magic32);            // Compute the first
        rhat = un64 - q1 * vn1;     // quotient digit, q1.
        while (true) {
            if (Long.compareUnsigned(q1, b) >= 0 || Long.compareUnsigned(q1 * vn0, b * rhat + un1) > 0) { //if (q1 >= b || q1 * vn0 > b * rhat + un1) {
                q1 = q1 - 1;
                rhat = rhat + vn1;
                if (Long.compareUnsigned(rhat, b) < 0)
                    continue;
            }
            break;
        }

        un21 = un64 * b + un1 - q1 * divider;  // Multiply and subtract.

        q0 = divideUnsignedFast(un21, magic32);            // Compute the second
        rhat = un21 - q0 * vn1;     // quotient digit, q0.
        while (true) {
            if (Long.compareUnsigned(q0, b) >= 0 || Long.compareUnsigned(q0 * vn0, b * rhat + un0) > 0) {
                q0 = q0 - 1;
                rhat = rhat + vn1;
                if (Long.compareUnsigned(rhat, b) < 0)
                    continue;
            }
            break;
        }
        long r = (un21 * b + un0 - q0 * divider) >>> s;    // return it.
        return r;
    }


//    static long modul64(long x, long y, long z) {
//        /*
//        Divides (x || y) by z, for 64-bit integers x, y,
//        and z, giving the remainder (modulus) as the result.
//        Must have x < z (to get a 64-bit result). This is
//        checked for.
//        */
//
//        long i, t;
//
//        if (x >= z)
//            throw new IllegalArgumentException();
//
//        for (i = 1; i <= 64; i++) {  // Do 64 times.
//            t = x >>> 63;       // All 1's if x(63) = 1.
//            x = (x << 1)|(y >>> 63); // Shift x || y left
//            y = y << 1;               // one bit.
//            if (Long.compareUnsigned((x|t), z) >= 0) {
//                x = x - z;
//                y = y + 1;
//            }
//        }
//        return x;                    // Quotient is y.
//    }

    public static final class MagicDivider {
        /** The magic number */
        final long magic;
        /** Shifting bits **/
        final int more;
        /** The divider **/
        final long divider;

        public MagicDivider(long magic, int more, long divider) {
            this.magic = magic;
            this.more = more;
            this.divider = divider;
        }
    }

    /**
     * Converts unsigned long to BigInteger
     *
     * @param bits unsigned bits
     * @return BigInteger value of unsigned long
     */
    static BigInteger valueOfUnsigned(long bits) {
        if (bits >= 0)
            return valueOfSigned(bits);
        BigInteger r = BigInteger.valueOf(~bits);
        for (int i = 0; i < 64; ++i)
            r = r.flipBit(i);
        return r;
    }

    /**
     * Converts signed long to BigInteger
     *
     * @param bits signed bits
     * @return BigInteger value of signed long
     */
    static BigInteger valueOfSigned(long bits) {
        return BigInteger.valueOf(bits);
    }
}
