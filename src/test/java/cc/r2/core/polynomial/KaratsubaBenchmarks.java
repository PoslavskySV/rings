package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by poslavsky on 10/01/2017.
 */
public class KaratsubaBenchmarks {

    static long[] multiplyKaratsuba(
            final long[] f, final int fFrom, final int fTo,
            final long[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return new long[0];

        // single element in f
        if (fTo - fFrom == 1) {
            long[] result = new long[gTo - gFrom];
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = LongArithmetics.multiply(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            long[] result = new long[fTo - fFrom];
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = LongArithmetics.multiply(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = LongArithmetics.multiply(f[fFrom], g[gFrom]);
            result[1] = LongArithmetics.add(LongArithmetics.multiply(f[fFrom], g[gFrom + 1]), LongArithmetics.multiply(f[fFrom + 1], g[gFrom]));
            result[2] = LongArithmetics.multiply(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if ((fTo - fFrom) * (gTo - gFrom) < 1024 /* MutableLongPoly.KARATSUBA_THRESHOLD */)
            return MutableLongPoly.multiplyClassical(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyKaratsuba(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            long[] f0g = multiplyKaratsuba(f, fFrom, fFrom + split, g, gFrom, gTo);
            long[] f1g = multiplyKaratsuba(f, fFrom + split, fTo, g, gFrom, gTo);

            long[] result = Arrays.copyOf(f0g, fTo - fFrom + gTo - gFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = LongArithmetics.add(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        long[] f0g0 = multiplyKaratsuba(f, fFrom, fMid, g, gFrom, gMid);
        long[] f1g1 = multiplyKaratsuba(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(fMid - fFrom, fTo - fMid)];
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = LongArithmetics.add(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        long[] g0_plus_g1 = new long[Math.max(gMid - gFrom, gTo - gMid)];
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = LongArithmetics.add(g0_plus_g1[i - gMid], g[i]);

        long[] mid = multiplyKaratsuba(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = LongArithmetics.subtract(mid[i], f1g1[i]);


        long[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = LongArithmetics.add(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = LongArithmetics.add(result[i + 2 * split], f1g1[i]);

        return result;
    }


    static long[] a_cache = null, b_cache = null, r_cache = null;

    static long[] get_a_cache(int len) {
        if (a_cache == null || a_cache.length < len)
            return (a_cache = new long[len]);
        return a_cache;
    }

    static long[] get_b_cache(int len) {
        if (b_cache == null || b_cache.length < len)
            return (b_cache = new long[len]);
        return b_cache;
    }

    static long[] get_r_cache(int len) {
        if (r_cache == null || r_cache.length < len)
            return (r_cache = new long[len]);
        return r_cache;
    }

    static long[] multiplyKaratsuba4e(final long[] a, final long[] b) {
        if (a.length < b.length)
            return multiplyKaratsuba4e(b, a);
        return multiplyKaratsuba4e(a, 0, a.length, b, 0, b.length, get_a_cache(a.length), 0, get_b_cache(a.length), 0);
    }

    static long[] multiplyKaratsuba4e(
            final long[] a, final int aFrom, final int aTo,
            final long[] b, final int bFrom, final int bTo,
            final long[] aCache, final int aCacheFrom,
            final long[] bCache, final int bCacheFrom) {
        if (aFrom >= aTo)
            return new long[0];
        if (bFrom >= bTo)
            return new long[0];

        if (aTo - aFrom == 1) {
            long[] result = new long[bTo - bFrom];
            //single element in a
            for (int i = bFrom; i < bTo; ++i)
                result[i - bFrom] = a[aFrom] * b[i];
            return result;
        }
        if (bTo - bFrom == 1) {
            long[] result = new long[aTo - aFrom];
            //single element in b
            for (int i = aFrom; i < aTo; ++i)
                result[i - aFrom] = b[bFrom] * a[i];
            return result;
        }
        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
            long[] result = new long[3];
            //both a and b are linear
            result[0] = a[aFrom] * b[bFrom];
            result[1] = a[aFrom] * b[aFrom + 1] + a[aFrom + 1] * b[bFrom]; ;
            result[2] = a[aFrom + 1] * b[bFrom + 1];
            return result;
        }
        //classical
        if ((aTo - aFrom) * (bTo - bFrom) < 10)
            return MutableLongPoly.multiplyClassical(a, aFrom, aTo, b, bFrom, bTo);
        if (aTo - aFrom < bTo - bFrom)
            return multiplyKaratsuba4e(b, bFrom, bTo, a, aFrom, aTo, bCache, bCacheFrom, aCache, aCacheFrom);


        //we now split a and b into 2 parts:
        int split = (aTo - aFrom + 1) / 2;
        //if we can't split b
        if (bFrom + split >= bTo) {
            long[] f0g = multiplyKaratsuba4e(a, aFrom, aFrom + split, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom);
            long[] f1g = multiplyKaratsuba4e(a, aFrom + split, aTo, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom);

            long[] result = Arrays.copyOf(f0g, aTo - aFrom + bTo - bFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] += f1g[i];
            return result;
        }

        int aMid = aFrom + split, bMid = bFrom + split;


        // f0 + f1

        int f0_plus_f1 = Math.max(aMid - aFrom, aTo - aMid);
        int g0_plus_g1 = Math.max(bMid - bFrom, bTo - bMid);

        for (int i = aFrom; i < aMid; i++)
            aCache[aCacheFrom + i - aFrom] = a[i];
        for (int i = aMid; i < aTo; ++i)
            aCache[aCacheFrom + i - aMid] += a[i];

        for (int i = bFrom; i < bMid; i++)
            bCache[bCacheFrom + i - bFrom] = b[i];
        for (int i = bMid; i < bTo; ++i)
            bCache[bCacheFrom + i - bMid] += b[i];

        long[] mid = multiplyKaratsuba4e(aCache, aCacheFrom, aCacheFrom + f0_plus_f1, bCache, bCacheFrom, bCacheFrom + g0_plus_g1, aCache, aCacheFrom + split, bCache, bCacheFrom + split);

        long[] f0g0 = multiplyKaratsuba4e(a, aFrom, aMid, b, bFrom, bMid, aCache, aCacheFrom + split, bCache, bCacheFrom + split);
        long[] f1g1 = multiplyKaratsuba4e(a, aMid, aTo, b, bMid, bTo, aCache, aCacheFrom + split, bCache, bCacheFrom + split);

        if (mid.length < f0g0.length)
            mid = Arrays.copyOf(mid, f0g0.length);
        if (mid.length < f1g1.length)
            mid = Arrays.copyOf(mid, f1g1.length);

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; i++)
            mid[i] -= f0g0[i];

        for (int i = 0; i < f1g1.length; i++)
            mid[i] -= f1g1[i];


        long[] result = Arrays.copyOf(f0g0, (aTo - aFrom) + (bTo - bFrom) - 1);
        for (int i = 0; i < mid.length; i++)
            result[i + split] += mid[i];

        for (int i = 0; i < f1g1.length; i++)
            result[i + 2 * split] += f1g1[i];

        return result;
    }


    static void multiplyClassical(final long[] result, final int resultFrom, final long[] a, final int aFrom, final int aTo, final long[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassical(result, resultFrom, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
//        Arrays.fill(result, resultFrom, resultFrom + aTo - aFrom + bTo - bFrom - 1, 0);
        for (int i = 0; i < aTo - aFrom; ++i) {
            long c = a[aFrom + i];
            if (c != 0)
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[resultFrom + i + j] = LongArithmetics.add(result[resultFrom + i + j], LongArithmetics.multiply(c, b[bFrom + j]));
        }
    }

    static long[] multiplyKaratsuba4f(final long[] a, final long[] b) {
        if (a.length < b.length)
            return multiplyKaratsuba4f(b, a);
        long[] result = get_r_cache(5 * a.length);
        Arrays.fill(result, 0);
        int pSize = multiplyKaratsuba4f(a, 0, a.length, b, 0, b.length, get_a_cache(2 * a.length), 0, get_b_cache(2 * a.length), 0, result, 0);
        return Arrays.copyOf(result, pSize);
    }


    static long ops = 0;

    static int multiplyKaratsuba4f(
            long[] a, int aFrom, int aTo,
            long[] b, int bFrom, int bTo,
            long[] aCache, int aCacheFrom,
            long[] bCache, int bCacheFrom,
            long[] result, int resultFrom) {
        if (aFrom >= aTo || bFrom >= bTo) {
            result[resultFrom] = 0;
            return 0;
        }

        if (aTo - aFrom == 1) {
            //single element in a
            for (int i = bFrom; i < bTo; ++i)
                result[resultFrom + (i - bFrom)] = a[aFrom] * b[i];
            return bTo - bFrom;
        }
        if (bTo - bFrom == 1) {
            //single element in b
            for (int i = aFrom; i < aTo; ++i)
                result[resultFrom + (i - aFrom)] = b[bFrom] * a[i];
            return aTo - aFrom;
        }
        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
            ops = LongArithmetics.add(ops, 3);
            //both a and b are linear
            result[resultFrom] = a[aFrom] * b[bFrom];
            result[resultFrom + 1] = a[aFrom] * b[bFrom + 1] + a[aFrom + 1] * b[bFrom]; ;
            result[resultFrom + 2] = a[aFrom + 1] * b[bFrom + 1];
            return 3;
        }
        //classical
        if ((0L + aTo - aFrom) * (0L + bTo - bFrom) < 2048 * 2) {
            ops = LongArithmetics.add(ops, (aTo - aFrom) * (bTo - bFrom));
            Arrays.fill(result, resultFrom, resultFrom + aTo - aFrom + bTo - bFrom - 1, 0);
            multiplyClassical(result, resultFrom, a, aFrom, aTo, b, bFrom, bTo);
            return aTo - aFrom + bTo - bFrom - 1;
        }

        if (aTo - aFrom < bTo - bFrom)
            return multiplyKaratsuba4f(b, bFrom, bTo, a, aFrom, aTo, bCache, bCacheFrom, aCache, aCacheFrom, result, resultFrom);

        //we now split a and b into 2 parts:
        int split = (aTo - aFrom + 1) / 2;
        //if we can't split b
        if (bFrom + split >= bTo) {

            //f0g
            Arrays.fill(result, resultFrom + 4 * split, resultFrom + 4 * split + 2 * split, 0);
            int f0g = multiplyKaratsuba4f(a, aFrom, aFrom + split, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom, result, resultFrom + 4 * split);
            assert f0g < 2 * split;
            System.arraycopy(result, resultFrom + 4 * split, result, resultFrom, f0g);


            Arrays.fill(result, resultFrom + 4 * split, resultFrom + 4 * split + 2 * split, 0);
            int f1g = multiplyKaratsuba4f(a, aFrom + split, aTo, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom, result, resultFrom + 4 * split);
            assert f1g < 2 * split;

            for (int i = 0; i < f1g; ++i)
                result[resultFrom + i + split] += result[resultFrom + i + 4 * split];
            return aTo - aFrom + bTo - bFrom - 1;
        }

        int aMid = aFrom + split, bMid = bFrom + split;


        int f0g0From = resultFrom;
        Arrays.fill(result, resultFrom + 4 * split, resultFrom + 4 * split + 2 * split, 0);
        int f0g0 = multiplyKaratsuba4f(
                a, aFrom, aMid,
                b, bFrom, bMid,
                aCache, aCacheFrom + split,
                bCache, bCacheFrom + split,
                result, resultFrom + 4 * split);
        assert f0g0 < 2 * split;
        System.arraycopy(result, resultFrom + 4 * split, result, f0g0From, f0g0);


        // f0 + f1

        int f0_plus_f1 = Math.max(aMid - aFrom, aTo - aMid);
        int g0_plus_g1 = Math.max(bMid - bFrom, bTo - bMid);

        System.arraycopy(a, aFrom, aCache, aCacheFrom, split);
        for (int i = aMid; i < aTo; ++i)
            aCache[aCacheFrom + i - aMid] += a[i];

        System.arraycopy(b, bFrom, bCache, bCacheFrom, split);
        for (int i = bMid; i < bTo; ++i)
            bCache[bCacheFrom + i - bMid] += b[i];

        int midFrom = resultFrom + 4 * split;
        Arrays.fill(result, midFrom, midFrom + 2 * split, 0);
        int mid = multiplyKaratsuba4f(
                aCache, aCacheFrom, aCacheFrom + f0_plus_f1,
                bCache, bCacheFrom, bCacheFrom + g0_plus_g1,
                aCache, aCacheFrom + split,
                bCache, bCacheFrom + split,
                result, midFrom);
        assert mid < 2 * split;

        // (f0+g0)*(f1+g1) - f0g0
        for (int i = 0; i < f0g0; i++)
            result[midFrom + i] -= result[f0g0From + i];

        mid = Math.max(mid, f0g0);
        for (int i = 0; i < mid; ++i)
            result[resultFrom + i + split] += result[midFrom + i];

        int f1g1From = resultFrom + 4 * split;
        Arrays.fill(result, f1g1From, f1g1From + 2 * split, 0);
        int f1g1 = multiplyKaratsuba4f(
                a, aMid, aTo,
                b, bMid, bTo,
                aCache, aCacheFrom + split,
                bCache, bCacheFrom + split,
                result, f1g1From);
        assert f1g1 < 2 * split;

        // (f0+g0)*(f1+g1) - f1g1
        for (int i = 0; i < f1g1; ++i)
            result[resultFrom + i + split] -= result[f1g1From + i];

        //write f1g1
        for (int i = 0; i < f1g1; ++i)
            result[resultFrom + i + 2 * split] += result[f1g1From + i];

        return (aTo - aFrom) + (bTo - bFrom) - 1;
    }

    static void squareClassical(final long[] result, int resultFrom, long[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            long c = data[from + i];
            if (c != 0)
                for (int j = 0; j < len; ++j)
                    result[resultFrom + i + j] = LongArithmetics.add(result[resultFrom + i + j], LongArithmetics.multiply(c, data[from + j]));
        }
    }

    static long[] squareKaratsuba4e(final long[] a) {
        long[] result = get_r_cache(5 * a.length);
        Arrays.fill(result, 0);
        int pSize = squareKaratsuba4e(a, 0, a.length, get_a_cache(2 * a.length), 0, result, 0);
        return Arrays.copyOf(result, pSize);
    }


    /**
     * Karatsuba squaring
     *
     * @param f     the data
     * @param fFrom begin in f
     * @param fTo   end in f
     * @return the result
     */
    static int squareKaratsuba4e(final long[] f, final int fFrom, final int fTo,
                                 long[] aCache, int aCacheFrom,
                                 long[] result, int resultFrom) {
        if (fFrom >= fTo) {
            result[resultFrom] = 0;
            return 0;
        }
        if (fTo - fFrom == 1) {
            result[resultFrom] = LongArithmetics.multiply(f[fFrom], f[fFrom]);
            return 1;
        }
        if (fTo - fFrom == 2) {
            result[resultFrom] = LongArithmetics.multiply(f[fFrom], f[fFrom]);
            result[resultFrom + 1] = LongArithmetics.multiply(2L, f[fFrom], f[fFrom + 1]);
            result[resultFrom + 2] = LongArithmetics.multiply(f[fFrom + 1], f[fFrom + 1]);
            return 3;
        }
        //switch to classical
        if (LongArithmetics.multiply(fTo - fFrom, fTo - fFrom) < 10) {
            Arrays.fill(result, resultFrom, resultFrom + 2 * (fTo - fFrom) - 1, 0);
            squareClassical(result, resultFrom, f, fFrom, fTo);
            return 2 * (fTo - fFrom) - 1;
        }


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;

        Arrays.fill(result, resultFrom, resultFrom + 2 * split, 0);
        int f0g0 = squareKaratsuba4e(f, fFrom, fMid, aCache, aCacheFrom + split, result, resultFrom);
        assert f0g0 < 2 * split;

        Arrays.fill(result, resultFrom + 2 * split, resultFrom + 4 * split, 0);
        int f1g1 = squareKaratsuba4e(f, fMid, fTo, aCache, aCacheFrom + split, result, resultFrom + 2 * split);
        assert f1g1 < 2 * split;

        // f0 + f1
        int f0_plus_f1 = Math.max(fMid - fFrom, fTo - fMid);
        System.arraycopy(f, fFrom, aCache, aCacheFrom, split);
        for (int i = fMid; i < fTo; ++i)
            aCache[aCacheFrom + i - fMid] += f[i];

        Arrays.fill(result, resultFrom + 4 * split, resultFrom + 6 * split, 0);
        int mid = squareKaratsuba4e(aCache, 0, f0_plus_f1, aCache, aCacheFrom + split, result, resultFrom + 4 * split);
        assert mid < 2 * split;


        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0; ++i)
            result[resultFrom + 4 * split + i] -= result[resultFrom + i];
        for (int i = 0; i < f1g1; ++i)
            result[resultFrom + 4 * split + i] -= result[resultFrom + 2 * split + i];


        mid = Math.max(Math.max(mid, f0g0), f1g1);
        for (int i = 0; i < mid; ++i)
            result[resultFrom + i + split] = result[resultFrom + i + 4 * split];

        return 2 * (fTo - fFrom) - 1;
    }

    @Test
    public void name() throws Exception {
        long[] a = {1, 1, 1, 1, 1, 1, 1};

        long[] r1 = MutableLongPoly.multiplyClassical(a, 0, a.length, a, 0, a.length);
        System.out.println(Arrays.toString(r1));
        long[] r2 = squareKaratsuba4e(a);
        System.out.println(Arrays.toString(r2));
        Assert.assertArrayEquals(r1, r2);
    }

    @Test
    public void name1() throws Exception {
        DescriptiveStatistics karaf = new DescriptiveStatistics(), clas = new DescriptiveStatistics();
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 15; i++) {
            System.out.println(i);
            if (i == 10) {
                karaf.clear();
                clas.clear();
            }
            long[] a = RandomPolynomials.randomLongArray(55000, 2, rnd);
            for (int j = 0; j < a.length; j++) {
                if (a[j] == 0)
                    a[j] = 1;
            }
            long[] b = RandomPolynomials.randomLongArray(55000, 2, rnd);
            for (int j = 0; j < b.length; j++) {
                if (b[j] == 0)
                    b[j] = 1;
            }

            long start = System.nanoTime();
            long[] r1 = MutableLongPoly.multiplyClassical(a, 0, a.length, b, 0, b.length);
            long cTime = System.nanoTime() - start;
            System.out.println(cTime);
            clas.addValue(cTime);

            start = System.nanoTime();
            ops = 0;
            long[] r2 = multiplyKaratsuba4f(a, b);
            System.out.println(ops);
            long kTime = System.nanoTime() - start;
            System.out.println(kTime);
            karaf.addValue(kTime);

            Assert.assertArrayEquals(r1, r2);

        }

        System.out.println("KARAT");
        System.out.println(karaf);

        System.out.println("Class");
        System.out.println(clas);

    }
}
