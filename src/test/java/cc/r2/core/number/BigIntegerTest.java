package cc.r2.core.number;

import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

/**
 * Created by poslavsky on 01/11/2016.
 */
public class BigIntegerTest {
    @Test
    public void test1() throws Exception {

        benSmall(20000, new BigInteger("1"), new BigInteger("1000"));
    }

    @Test
    public void test2() throws Exception {
        int nTries = 500_000;
        int start = Integer.MAX_VALUE - nTries;
        int stop = Integer.MAX_VALUE;

        BigInteger zero = BigInteger.ZERO;

        DescriptiveStatistics stats1 = new DescriptiveStatistics();
        DescriptiveStatistics stats2 = new DescriptiveStatistics();
        for (int i = 0; i < 25; i++) {
            R r1 = gcdTime1(zero, start, stop);
            zero = r1.res;
            R r2 = gcdTime2(zero, start, stop);
            zero = r2.res;

            stats1.addValue(r1.time);
            stats2.addValue(r2.time);
        }

        System.out.println(stats1);
        System.out.println(stats2);
        System.out.println(zero.mod(new BigInteger("2")));
    }


    static R gcdTime1(BigInteger zero, int start, int stop) {
        long t1 = System.nanoTime();
        for (int i = start; i < stop; i++) {
            BigInteger a = BigInteger.valueOf(i);
            BigInteger b = BigInteger.valueOf(stop - i);
            BigInteger gcd1 = a.gcd(b);

            zero = zero.add(gcd1);
        }
        t1 = System.nanoTime() - t1;
        return new R(t1, zero);
    }

    static R gcdTime2(BigInteger zero, int start, int stop) {
        long t2 = System.nanoTime();
        for (int i = start; i < stop; i++) {
            BigInteger a = BigInteger.valueOf(i);
            BigInteger b = BigInteger.valueOf(stop - i);

            BigInteger gcd2 = b.gcd(a);

            zero = zero.add(gcd2);
        }
        t2 = System.nanoTime() - t2;
        return new R(t2, zero);
    }

    static void benSmall(final int nTries, BigInteger a, BigInteger b) {
        TLongList bres = new TLongArrayList();
        DescriptiveStatistics bStats = new DescriptiveStatistics();
        for (int i = 0; i < nTries; i++) {
            R r = benMulBigInt(a.add(i), b.add(i));
            bStats.addValue(r.time);
            bres.add(r.res.hashCode());
        }

        TLongList lres = new TLongArrayList();
        DescriptiveStatistics lStats = new DescriptiveStatistics();
        for (int i = 0; i < nTries; i++) {
            R r = benMulLong(a.add(i).longValue(), b.add(i).longValue());
            lStats.addValue(r.time);
            lres.add(r.res.hashCode());
        }

        System.out.println(bres.size() > 10 ? bres.subList(0, 10) : bres);
        System.out.println(lres.size() > 10 ? lres.subList(0, 10) : lres);

        System.out.println("BigInteger");
        System.out.println(bStats);


        System.out.println("Long");
        System.out.println(lStats);


        System.out.println("\n\nMean");
        System.out.println("BigInteger: " + bStats.getMean());
        System.out.println("Long:       " + lStats.getMean());
    }

    static R benMulBigInt(BigInteger start, BigInteger stop) {
        long tStart = System.nanoTime();
        int d = stop.subtract(start).abs().intValue();
        BigInteger r = BigInteger.ONE;
        for (int i = 0; i < d; i++) {
            BigInteger tmp = start.add(i);
            r = tmp.add(r.multiply(tmp));
        }
        return new R(((System.nanoTime() - tStart)), r);
    }

    static R benMulLong(long start, long stop) {
        long tStart = System.nanoTime();
        long d = Math.abs(stop - start);
        long r = 1L;
        for (int i = 0; i < d; i++) {
            long tmp = start + i;
            r = tmp + r * tmp;
        }
        long t = ((System.nanoTime() - tStart));
        return new R(t, new BigInteger(Long.toString(r)));
    }

    static final class R {
        final long time;
        final BigInteger res;

        public R(long time, BigInteger res) {
            this.time = time;
            this.res = res;
        }
    }

    @Test
    public void asdasdas() throws Exception {
        System.out.println(pow(BigInteger.valueOf(3), 4));
    }

    static <R extends RingElement<R>> R pow(R val, int p) {
        if (p < 0)
            throw new IllegalArgumentException();
        if (p == 0)
            return val.getOne();
        if (p == 1)
            return val;

        R result = val.getOne();
        R k2p = val;
        while (p != 0) {
            if ((p & 0x1) != 0) {
                result = result.multiply(k2p);
            }
            k2p = k2p.multiply(k2p);
            p = p >> 1;
        }

        return result;
    }
}