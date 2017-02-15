package cc.r2.core.polynomial;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.Arrays;

import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 25/01/2017.
 */
public class HackersDelight {

    static final int LIBDIVIDE_U64_SHIFT_PATH = 0x80;
    static final int LIBDIVIDE_ADD_MARKER = 0x40;
    static final int LIBDIVIDE_64_SHIFT_MASK = 0x3F;
    static final int LIBDIVIDE_NEGATIVE_DIVISOR = 0x80;
    static final int INT_BITS = 30;
    static final int LONG_BITS = 60;
    static final int INT_MASK = 0x3FFFFFFF;//0xFFFFFFFF;

    static final class libdivide_s64_t {
        final long magic;
        final int more;
        final long divider;

        public libdivide_s64_t(long magic, int more, long divider) {
            this.magic = magic;
            this.more = more;
            this.divider = divider;
        }
    }

    static int libdivide__count_leading_zeros64(long v) {
        return Long.numberOfLeadingZeros(v);
    }

    static long libdivide__mullhi_s64(long x, long y) {
        final long x_hi = x >> 32;
        final long y_hi = y >> 32;
        final long x_lo = x&0xFFFFFFFFL;
        final long y_lo = y&0xFFFFFFFFL;
        long result = x_lo * y_lo;
        result >>= 32;

        result += x_hi * y_lo + x_lo * y_hi;
        result >>= 32;
        result += x_hi * y_hi;
        return result;
    }


    public static long[] libdivide_128_div_64_to_64(long u1, long u0, long v) {
        long b = (1L << 32); // Number base (16 bits).
        long un1, un0,  // Norm. dividend LSD's.
                vn1, vn0,           // Norm. divisor digits.
                q1, q0,             // Quotient digits.
                un64, un21, un10,   // Dividend digit pairs.
                rhat;               // A remainder.
        int s;              // Shift amount for norm.

        if (u1 >= v) {                  // If overflow, set rem.
            return new long[]{-1L, -1L};      // possible quotient.
        }

        // count leading zeros
        s = libdivide__count_leading_zeros64(v); // 0 <= s <= 63.
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


    static libdivide_s64_t libdivide_internal_s64_gen(long d, boolean branchfree) {

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
        int floor_log_2_d = 63 - libdivide__count_leading_zeros64(absD);
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
            long[] tmp = libdivide_128_div_64_to_64(1L << (floor_log_2_d - 1), 0, absD);
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
        return new libdivide_s64_t(resultMagic, resultMore, d);
    }


    static libdivide_s64_t libdivide_s64_gen(long d) {
        return libdivide_internal_s64_gen(d, false);
    }

    static long libdivide_s64_do(long numer, libdivide_s64_t denom) {
        int more = denom.more;
        long magic = denom.magic;
        if (magic == 0) { //shift path
            int shifter = more&0x3F;
            long uq = numer + ((numer >> 63)&((1L << shifter) - 1));
            long q = uq;
            q = q >> shifter;
            // must be arithmetic shift and then sign-extend
            long shiftMask = more >> 7;
            q = (q^shiftMask) - shiftMask;
            return q;
        } else {
            long uq = libdivide__mullhi_s64(magic, numer);
            if ((more&0x40) != 0) {
                // must be arithmetic shift and then sign extend
                long sign = more >> 7;
                uq += (((long) numer^sign) - sign);
            }
            long q = (long) uq;
            q >>= more&0x3F;
            if (q < 0)
                q += 1;
            return q;
        }
    }

    static long libdivide_s64_do_rem(long numer, libdivide_s64_t denom) {
        long quot = libdivide_s64_do(numer, denom);
        return numer - quot * denom.divider;
    }


    static long[] reduceModBenchFast(long[] arr, libdivide_s64_t magic) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = libdivide_s64_do_rem(tmp[j], magic);
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    static long[] reduceModBenchPlain(long[] arr, long modulus) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = tmp[j] % modulus;
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    @Test
    public void name() throws Exception {
        long divider = 71232;
        long dividend = 1L << 55;
        libdivide_s64_t rem = libdivide_internal_s64_gen(divider, false);

        System.out.println(libdivide_s64_do(dividend, rem));
        System.out.println(dividend / divider);

    }

    @Test
    public void div() throws Exception {

        long u1 = 1 << 4;
        long u0 = 0;//0x20FaFFa0L;
        long div = 111232L;

        BigInteger bigInteger = BigInteger.valueOf(u1).shiftLeft(32).or(BigInteger.valueOf(u0));
        System.out.println(bigInteger.toString(16));
        System.out.println(bigInteger.divide(BigInteger.valueOf(div)));
        System.out.println(bigInteger.remainder(BigInteger.valueOf(div)));

        System.out.println(Arrays.toString(libdivide_128_div_64_to_64(u1, u0, div)));

    }


    @Test
    public void testd7_bench() throws Exception {
//        long[] modulus = {3};
        RandomGenerator rnd = new Well1024a();
        DescriptiveStatistics plain = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        for (int i = 0; i < 100000; i++) {
            if (i == 10000) {
                fast.clear();
                plain.clear();
            }
            long[] arr = new long[200];
            for (int j = 0; j < arr.length; j++) {
                arr[j] = rnd.nextInt();
                if (arr[j] < 0) arr[j] = -arr[j];
            }

            long modulus;
            do {
                modulus = rnd.nextLong();
                modulus = modulus % 100;
                if (modulus < 0)
                    modulus = -modulus;
            } while (modulus == 0);


            long[] f = reduceModBenchFast(arr, libdivide_internal_s64_gen(modulus, false));
            long[] p = reduceModBenchPlain(arr, modulus);

            assertEquals(f[1], p[1]);

            fast.addValue(f[0]);
            plain.addValue(p[0]);
        }

        System.out.println("==== fast ====");
        System.out.println(fast);
        System.out.println("==== plain ====");
        System.out.println(plain);
    }
}
