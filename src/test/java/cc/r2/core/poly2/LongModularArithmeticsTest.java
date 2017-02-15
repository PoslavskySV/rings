package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly2.LongModularArithmetics.*;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly2.LongModularArithmetics.*;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 26/01/2017.
 */
public class LongModularArithmeticsTest {

    static long lowBits(BigInteger num) {
        return num.and(BigInteger.valueOf(0xFFFFFFFFFFFFFFFFL)).longValue();
    }

    static long highBits(BigInteger num) {
        return num.shiftRight(64).longValue();
    }

    @Test
    public void testMulHighRandom1_unsigned() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10000; i++) {

            long x = rnd.nextLong();
            long y = rnd.nextLong();

            long high = multiplyHighUnsigned(x, y);
            long low = multiplyLow(x, y);

            BigInteger num = valueOfUnsigned(y).multiply(valueOfUnsigned(x));
            long highExpected = highBits(num);
            long lowExpected = lowBits(num);

            String errMsg = Long.toHexString(x) + "  " + Long.toHexString(y);
            assertEquals(errMsg, highExpected, high);
            assertEquals(errMsg, lowExpected, low);
        }
    }

    @Test
    public void testMulHighRandom1_signed() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10000; i++) {

            long x = rnd.nextLong();
            long y = rnd.nextLong();

            long high = multiplyHighSigned(x, y);
            long low = multiplyLow(x, y);

            BigInteger num = valueOfSigned(y).multiply(valueOfSigned(x));
            long highExpected = highBits(num);
            long lowExpected = lowBits(num);

            String errMsg = Long.toHexString(x) + "  " + Long.toHexString(y);
            assertEquals(errMsg, highExpected, high);
            assertEquals(errMsg, lowExpected, low);
        }
    }

    @Test
    public void testMulHighUnsigned1() throws Exception {
        for (long x : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long y : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE}) {
                long high = multiplyHighUnsigned(x, y);
                long low = multiplyLow(x, y);

                BigInteger num = valueOfUnsigned(x).multiply(valueOfUnsigned(y));
                long highExpected = highBits(num);
                long lowExpected = lowBits(num);
                assertEquals(highExpected, high);
                assertEquals(lowExpected, low);
            }
    }

    @Test
    public void testMulHighSigned1() throws Exception {
        for (long x : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long y : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE}) {
                long high = multiplyHighSigned(x, y);
                long low = multiplyLow(x, y);

                BigInteger num = valueOfSigned(x).multiply(valueOfSigned(y));
                long highExpected = highBits(num);
                long lowExpected = lowBits(num);
                assertEquals(highExpected, high);
                assertEquals(lowExpected, low);
            }
    }

    @Test
    public void testDivideFastRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 1000_000; i++) {
            long dividend = rnd.nextLong();
            long divider = rnd.nextLong();

            MagicDivider magic;
            magic = magicUnsigned(divider, true);
            Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));
            magic = magicUnsigned(divider, false);
            Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));

            magic = magicSigned(divider, true);
            Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));
            magic = magicSigned(divider, false);
            Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));
        }
    }

    @Test
    public void testDivideFast1() throws Exception {
        long dividend = -1;
        long divider = 324234234L;
        MagicDivider magic = magicUnsigned(divider);

        BigInteger bDividend = valueOfUnsigned(dividend);
        BigInteger bDivider = valueOfUnsigned(divider);

        Assert.assertEquals(bDividend.divide(bDivider).longValue(), divideUnsignedFast(dividend, magic));
    }

    @Test
    public void testDivideFast2() throws Exception {
        for (long dividend : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long divider : new long[]{Long.MAX_VALUE}) {
                MagicDivider magic = magicSigned(divider, true);
                Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));
                magic = magicSigned(divider, false);
                Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));

            }
    }

    @Test
    public void testDivideFast3() throws Exception {
        for (long dividend : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long divider : new long[]{-1L, Long.MAX_VALUE, Long.MIN_VALUE}) {
                MagicDivider magic = magicUnsigned(divider, true);
                Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));
                magic = magicUnsigned(divider, false);
                Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));
            }
    }

    @Test
    public void testDivideFast4() throws Exception {
        assertEquals(1432344, divideSignedFast(1432344, magicSigned(1)));
        assertEquals(-1432344, divideSignedFast(1432344, magicSigned(-1)));
    }

    @Test
    public void testDivideFast5() throws Exception {
        long dividend = -941192;
        long divider = -8;
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, true)));
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, false)));
    }

    @Test
    public void testDivideFast6() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10_000; i++) {
            long dividend = rnd.nextLong();
            for (long sign : new long[]{-1, 1})
                for (int p = 1; p < 55; p++) {
                    long divider = sign * (1L << p);

                    MagicDivider magic;
                    magic = magicUnsigned(divider, true);
                    Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));
                    magic = magicUnsigned(divider, false);
                    Assert.assertEquals(dividend + "/" + divider, valueOfUnsigned(dividend).divide(valueOfUnsigned(divider)).longValue(), divideUnsignedFast(dividend, magic));

                    magic = magicSigned(divider, true);
                    Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));
                    magic = magicSigned(divider, false);
                    Assert.assertEquals(dividend + "/" + divider, dividend / divider, divideSignedFast(dividend, magic));
                }
        }
    }

    @Test
    public void testDivideFast7() throws Exception {
        long dividend = -3188646;
        long divider = -6;
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, false)));
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, true)));
    }

    @Test
    public void testDivideFast8() throws Exception {
        long dividend = -108969;
        long divider = -3;
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, false)));
        assertEquals(dividend / divider, divideSignedFast(dividend, magicSigned(divider, true)));
    }

    @Test
    public void testMulMod128() {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 1000_000; i++) {
            long a = rnd.nextLong(), b = rnd.nextLong(), modulus;
            do {
                modulus = rnd.nextLong();
            } while (Long.compareUnsigned(modulus, MAX_SUPPORTED_MODULUS) >= 0);

            BigInteger expected = valueOfUnsigned(a).multiply(valueOfUnsigned(b)).mod(valueOfUnsigned(modulus));
            MagicDivider magic = magicUnsigned(modulus);
            a = modUnsignedFast(a, magic);
            b = modUnsignedFast(b, magic);
            MagicDivider magic32 = magic32ForMultiplyMod(modulus);
            Assert.assertEquals(a + "*" + b + " mod " + modulus, expected.longValue(), multiplyMod128Unsigned(a, b, modulus, magic32));
        }
    }

    static long[] modulusBenchmarkFast(long[] arr, MagicDivider magic) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = modSignedFast(tmp[j], magic);
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    static long[] modulusBenchmarkPlain(long[] arr, long modulus) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = Math.floorMod(tmp[j], modulus);
                r += tmp[j];
            }
            timing += System.nanoTime() - start;
        }
        return new long[]{timing, r};
    }

    @Test
    public void fastDivisionBenchmark1() throws Exception {
        boolean small = true;
        RandomGenerator rnd = new Well1024a();
        DescriptiveStatistics plain = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        for (int i = 0; i < 100000; i++) {
            if (i == 10000) {
                fast.clear();
                plain.clear();
            }
            long[] arr = new long[1000];
            for (int j = 0; j < arr.length; j++)
                arr[j] = rnd.nextLong();

            long modulus;
            if (small)
                do {
                    modulus = rnd.nextLong();
                    modulus = modulus % 100;
                    if (modulus < 0)
                        modulus = -modulus;
                } while (modulus == 0);
            else {
                modulus = rnd.nextLong();
                if (modulus < 0)
                    modulus = -modulus;
            }
            if (modulus == 0 || modulus == 1)
                modulus = 123;


            long[] f = modulusBenchmarkFast(arr, magicSigned(modulus));
            long[] p = modulusBenchmarkPlain(arr, modulus);

            assertEquals(f[1], p[1]);

            fast.addValue(f[0]);
            plain.addValue(p[0]);
        }

        System.out.println("==== Fast long modulus ====");
        System.out.println(fast);
        System.out.println("==== Plain long modulus ====");
        System.out.println(plain);
    }
}