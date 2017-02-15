package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 26/01/2017.
 */
public class LongModularArithmeticsTest {
    @Test
    public void asdasdasd() throws Exception {
        long a = -123124L;
        long b = 123L;
        System.out.println(1.0 * a / b);
        System.out.println(Math.floorDiv(a, b));
        System.out.println(a / b);
        System.out.println(LongModularArithmetics.divideSignedFast(a, LongModularArithmetics.magicSigned(b)));

        System.out.println(a);
        System.out.println(-1001L * b);
        System.out.println(-1002L * b);


    }

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

            long high = LongModularArithmetics.multiplyHighUnsigned(x, y);
            long low = LongModularArithmetics.multiplyLow(x, y);

            BigInteger num = LongModularArithmetics.valueOfUnsigned(y).multiply(LongModularArithmetics.valueOfUnsigned(x));
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

            long high = LongModularArithmetics.multiplyHighSigned(x, y);
            long low = LongModularArithmetics.multiplyLow(x, y);

            BigInteger num = LongModularArithmetics.valueOfSigned(y).multiply(LongModularArithmetics.valueOfSigned(x));
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
                long high = LongModularArithmetics.multiplyHighUnsigned(x, y);
                long low = LongModularArithmetics.multiplyLow(x, y);

                BigInteger num = LongModularArithmetics.valueOfUnsigned(x).multiply(LongModularArithmetics.valueOfUnsigned(y));
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
                long high = LongModularArithmetics.multiplyHighSigned(x, y);
                long low = LongModularArithmetics.multiplyLow(x, y);

                BigInteger num = LongModularArithmetics.valueOfSigned(x).multiply(LongModularArithmetics.valueOfSigned(y));
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

            LongModularArithmetics.MagicDivider magic;
            magic = LongModularArithmetics.magicUnsigned(divider);
            Assert.assertEquals(dividend + "/" + divider, LongModularArithmetics.valueOfUnsigned(dividend).divide(LongModularArithmetics.valueOfUnsigned(divider)).longValue(), LongModularArithmetics.divideUnsignedFast(dividend, magic));

            if (divider < 0) divider = -divider;
            magic = LongModularArithmetics.magicSigned(divider);
            Assert.assertEquals(dividend + "/" + divider, dividend / divider, LongModularArithmetics.divideSignedFast(dividend, magic));
        }
    }

    @Test
    public void testDivideFast1() throws Exception {
        long dividend = -1;
        long divider = 324234234L;
        LongModularArithmetics.MagicDivider magic = LongModularArithmetics.magicUnsigned(divider);

        BigInteger bDividend = LongModularArithmetics.valueOfUnsigned(dividend);
        BigInteger bDivider = LongModularArithmetics.valueOfUnsigned(divider);

        Assert.assertEquals(bDividend.divide(bDivider).longValue(), LongModularArithmetics.divideUnsignedFast(dividend, magic));
    }

    @Test
    public void testDivideFast2() throws Exception {
        for (long dividend : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long divider : new long[]{Long.MAX_VALUE}) {
                LongModularArithmetics.MagicDivider magic = LongModularArithmetics.magicSigned(divider);
                Assert.assertEquals(dividend + "/" + divider, dividend / divider, LongModularArithmetics.divideSignedFast(dividend, magic));

            }
    }

    @Test
    public void testDivideFast3() throws Exception {
        for (long dividend : new long[]{-1L, 1L, Long.MAX_VALUE, Long.MIN_VALUE})
            for (long divider : new long[]{-1L, Long.MAX_VALUE, Long.MIN_VALUE}) {
                LongModularArithmetics.MagicDivider magic = LongModularArithmetics.magicUnsigned(divider);
                Assert.assertEquals(dividend + "/" + divider, LongModularArithmetics.valueOfUnsigned(dividend).divide(LongModularArithmetics.valueOfUnsigned(divider)).longValue(), LongModularArithmetics.divideUnsignedFast(dividend, magic));
            }
    }

    @Test
    public void testMulMod128() {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 1000_000; i++) {
            long a = rnd.nextLong(), b = rnd.nextLong(), modulus;
            do {
                modulus = rnd.nextLong();
            } while (Long.compareUnsigned(modulus, LongModularArithmetics.MAX_SUPPORTED_MODULUS) >= 0);

            BigInteger expected = LongModularArithmetics.valueOfUnsigned(a).multiply(LongModularArithmetics.valueOfUnsigned(b)).mod(LongModularArithmetics.valueOfUnsigned(modulus));
            LongModularArithmetics.MagicDivider magic = LongModularArithmetics.magicUnsigned(modulus);
            a = LongModularArithmetics.modUnsignedFast(a, magic);
            b = LongModularArithmetics.modUnsignedFast(b, magic);
            LongModularArithmetics.MagicDivider magic32 = LongModularArithmetics.magic32ForMultiplyMod(modulus);
            Assert.assertEquals(a + "*" + b + " mod " + modulus, expected.longValue(), LongModularArithmetics.multiplyMod128Unsigned(a, b, modulus, magic32));
        }
    }

    static long[] modulusBenchmarkFast(long[] arr, LongModularArithmetics.MagicDivider magic) {
        long r = 0;
        long timing = 0;
        for (int i = 0; i < 10; i++) {
            long[] tmp = arr.clone();
            long start = System.nanoTime();
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = LongModularArithmetics.modSignedFast(tmp[j], magic);
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


            long[] f = modulusBenchmarkFast(arr, LongModularArithmetics.magicSigned(modulus));
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