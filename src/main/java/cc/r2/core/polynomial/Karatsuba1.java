package cc.r2.core.polynomial;

import java.util.Arrays;

/**
 * Created by poslavsky on 06/01/2017.
 */
public class Karatsuba1 {

    static long[] multiplyNaive(final long[] a, final long[] b, final int aFrom, final int aTo, final int bFrom, final int bTo) {
        long[] result = new long[aTo - aFrom + bTo - bFrom - 1];
        for (int i = 0; i < aTo - aFrom; i++) {
            for (int j = 0; j < bTo - bFrom; j++) {
                result[i + j] += a[aFrom + i] * b[bFrom + j];
            }
        }
        return result;
    }

    static long[] multiplyKaratsuba4a(
            final long[] a, final int aFrom, final int aTo,
            final long[] b, final int bFrom, final int bTo) {
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
        if ((aTo - aFrom) * (bTo - bFrom) < 1000)
            return multiplyNaive(a, b, aFrom, aTo, bFrom, bTo);

        int imid = (Math.max(aTo - aFrom, bTo - bFrom)) / 2;

        //f0*g0
        int aMid = Math.min(aFrom + imid, aTo);
        int bMid = Math.min(bFrom + imid, bTo);

        long[] f0g0 = multiplyKaratsuba4a(a, aFrom, aMid, b, bFrom, bMid);
        long[] f1g1 = multiplyKaratsuba4a(a, aMid, aTo, b, bMid, bTo);


        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(aMid - aFrom, aTo - aMid)];
        long[] g0_plus_g1 = new long[Math.max(bMid - bFrom, bTo - bMid)];

        for (int i = aFrom; i < aMid; i++)
            f0_plus_f1[i - aFrom] = a[i];
        for (int i = aMid; i < aTo; ++i)
            f0_plus_f1[i - aMid] += a[i];

        for (int i = bFrom; i < bMid; i++)
            g0_plus_g1[i - bFrom] = b[i];
        for (int i = bMid; i < bTo; ++i)
            g0_plus_g1[i - bMid] += b[i];

        long[] mid = multiplyKaratsuba4a(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

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
            result[i + imid] += mid[i];

        for (int i = 0; i < f1g1.length; i++)
            result[i + 2 * imid] += f1g1[i];

        return result;
    }

    static long[] multiplyKaratsuba4a(final long[] a, final long[] b) {
        return multiplyKaratsuba4a(a, 0, a.length, b, 0, b.length);
    }

    static long[] multiplyKaratsuba4c(final long[] a, final long[] b) {
        if (a.length < b.length)
            return multiplyKaratsuba4c(b, a);
        return multiplyKaratsuba4c(a, 0, a.length, b, 0, b.length);
    }

    static long[] multiplyKaratsuba4c(
            final long[] a, final int aFrom, final int aTo,
            final long[] b, final int bFrom, final int bTo) {
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
        if ((aTo - aFrom) * (bTo - bFrom) < 1000)
            return multiplyNaive(a, b, aFrom, aTo, bFrom, bTo);
        if (aTo - aFrom < bTo - bFrom)
            return multiplyKaratsuba4c(b, bFrom, bTo, a, aFrom, aTo);


        //we now split a into 2 parts:
        int split = (aTo - aFrom + 1) / 2;
        //if we can't split b
        if (bFrom + split >= bTo) {
            long[] f0g = multiplyKaratsuba4c(a, aFrom, aFrom + split, b, bFrom, bTo);
            long[] f1g = multiplyKaratsuba4c(a, aFrom + split, aTo, b, bFrom, bTo);

            long[] result = Arrays.copyOf(f0g, aTo - aFrom + bTo - bFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] += f1g[i];
            return result;
        }


        int imid = (Math.max(aTo - aFrom, bTo - bFrom) + 1) / 2;

        //f0*g0
        int aMid = Math.min(aFrom + imid, aTo);
        int bMid = Math.min(bFrom + imid, bTo);

        long[] f0g0 = multiplyKaratsuba4c(a, aFrom, aMid, b, bFrom, bMid);
        long[] f1g1 = multiplyKaratsuba4c(a, aMid, aTo, b, bMid, bTo);


        // f0 + f1
        long[] f0_plus_f1 = new long[Math.max(aMid - aFrom, aTo - aMid)];
        long[] g0_plus_g1 = new long[Math.max(bMid - bFrom, bTo - bMid)];

        for (int i = aFrom; i < aMid; i++)
            f0_plus_f1[i - aFrom] = a[i];
        for (int i = aMid; i < aTo; ++i)
            f0_plus_f1[i - aMid] += a[i];

        for (int i = bFrom; i < bMid; i++)
            g0_plus_g1[i - bFrom] = b[i];
        for (int i = bMid; i < bTo; ++i)
            g0_plus_g1[i - bMid] += b[i];

        long[] mid = multiplyKaratsuba4c(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

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
            result[i + imid] += mid[i];

        for (int i = 0; i < f1g1.length; i++)
            result[i + 2 * imid] += f1g1[i];

        return result;
    }


    static long[] multiplyKaratsuba4d(final long[] a, final long[] b) {
        if (a.length < b.length)
            return multiplyKaratsuba4d(b, a);
        return multiplyKaratsuba4d(a, 0, a.length, b, 0, b.length);
    }

    static long[] multiplyKaratsuba4d(
            final long[] a, final int aFrom, final int aTo,
            final long[] b, final int bFrom, final int bTo) {
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
        if ((aTo - aFrom) * (bTo - bFrom) < 1000)
            return multiplyNaive(a, b, aFrom, aTo, bFrom, bTo);
        if (aTo - aFrom < bTo - bFrom)
            return multiplyKaratsuba4d(b, bFrom, bTo, a, aFrom, aTo);


        //we now split a and b into 2 parts:
        int split = (aTo - aFrom + 1) / 2;
        //if we can't split b
        if (bFrom + split >= bTo) {
            long[] f0g = multiplyKaratsuba4d(a, aFrom, aFrom + split, b, bFrom, bTo);
            long[] f1g = multiplyKaratsuba4d(a, aFrom + split, aTo, b, bFrom, bTo);

            long[] result = Arrays.copyOf(f0g, aTo - aFrom + bTo - bFrom - 1);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] += f1g[i];
            return result;
        }

        int aMid = aFrom + split, bMid = bFrom + split;

        long[] f0g0 = multiplyKaratsuba4d(a, aFrom, aMid, b, bFrom, bMid);
        long[] f1g1 = multiplyKaratsuba4d(a, aMid, aTo, b, bMid, bTo);


        // f0 + f1

        long[] f0_plus_f1 = new long[Math.max(aMid - aFrom, aTo - aMid)];
        long[] g0_plus_g1 = new long[Math.max(bMid - bFrom, bTo - bMid)];

        for (int i = aFrom; i < aMid; i++)
            f0_plus_f1[i - aFrom] = a[i];
        for (int i = aMid; i < aTo; ++i)
            f0_plus_f1[i - aMid] += a[i];

        for (int i = bFrom; i < bMid; i++)
            g0_plus_g1[i - bFrom] = b[i];
        for (int i = bMid; i < bTo; ++i)
            g0_plus_g1[i - bMid] += b[i];

        long[] mid = multiplyKaratsuba4d(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

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


    static long[] a_cache = null;
    static long[] get_a_cache(int len){
        if(a_cache == null || a_cache.length < len)
            return (a_cache = new long[len]);
        return a_cache;
    }

    static long[] b_cache = null;
    static long[] get_b_cache(int len){
        if(b_cache == null || b_cache.length < len)
            return (b_cache = new long[len]);
        return b_cache;
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
        if ((aTo - aFrom) * (bTo - bFrom) < 1000)
            return multiplyNaive(a, b, aFrom, aTo, bFrom, bTo);
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


//    static long[] multiplyKaratsuba4f(final long[] a, final long[] b) {
//        if (a.length < b.length)
//            return multiplyKaratsuba4f(b, a);
//        long[] result = new long[2 * a.length];
//        multiplyKaratsuba4f(a, 0, a.length, b, 0, b.length, new long[a.length], 0, new long[a.length], 0, result, 0);
//        return result;
//    }
//
//    static void multiplyNaive(final long[] result, int resultFrom, final long[] a, final long[] b, final int aFrom, final int aTo, final int bFrom, final int bTo) {
//        for (int i = 0; i < aTo - aFrom; i++) {
//            for (int j = 0; j < bTo - bFrom; j++) {
//                result[resultFrom + i + j] += a[aFrom + i] * b[bFrom + j];
//            }
//        }
//    }
//
//    static void multiplyKaratsuba4f(
//            final long[] a, final int aFrom, final int aTo,
//            final long[] b, final int bFrom, final int bTo,
//            final long[] aCache, final int aCacheFrom,
//            final long[] bCache, final int bCacheFrom,
//            final long[] result, final int resultFrom) {
//        if (aFrom >= aTo) {
//            result[resultFrom] = 0;
//            return;
//        }
//        if (bFrom >= bTo) {
//            result[resultFrom] = 0;
//            return;
//        }
//
//        if (aTo - aFrom == 1) {
//            //single element in a
//            for (int i = bFrom; i < bTo; ++i)
//                result[resultFrom + i - bFrom] = a[aFrom] * b[i];
//            return;
//        }
//        if (bTo - bFrom == 1) {
//            //single element in b
//            for (int i = aFrom; i < aTo; ++i)
//                result[resultFrom + i - aFrom] = b[bFrom] * a[i];
//            return;
//        }
//        if (aTo - aFrom == 2 && bTo - bFrom == 2) {
//            //both a and b are linear
//            result[resultFrom + 0] = a[aFrom] * b[bFrom];
//            result[resultFrom + 1] = a[aFrom] * b[aFrom + 1] + a[aFrom + 1] * b[bFrom]; ;
//            result[resultFrom + 2] = a[aFrom + 1] * b[bFrom + 1];
//            return;
//        }
//        //classical
//        if ((aTo - aFrom) * (bTo - bFrom) < 1000) {
//            multiplyNaive(result, resultFrom, a, b, aFrom, aTo, bFrom, bTo);
//            return;
//        }
//        if (aTo - aFrom < bTo - bFrom) {
//            multiplyKaratsuba4f(b, bFrom, bTo, a, aFrom, aTo, bCache, bCacheFrom, aCache, aCacheFrom, result, resultFrom);
//            return;
//        }
//
//        //we now split a and b into 2 parts:
//        int split = (aTo - aFrom + 1) / 2;
//        //if we can't split b
//        if (bFrom + split >= bTo) {
//            int f0gLen = split + bTo - bFrom - 1;
//            int f1gLen = (aTo - split) + bTo - bFrom - 1;
//            multiplyKaratsuba4f(a, aFrom, aFrom + split, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom, result, resultFrom);
//            multiplyKaratsuba4f(a, aFrom + split, aTo, b, bFrom, bTo, aCache, aCacheFrom, bCache, bCacheFrom, result, resultFrom + f0gLen);
//
//            for (int i = 0; i < f1gLen; i++)
//                result[resultFrom + i + split] += result[resultFrom + f0gLen + i];
//            return;
//        }
//
//        int aMid = aFrom + split, bMid = bFrom + split;
//
//
//        // f0 + f1
//
//        int f0_plus_f1 = Math.max(aMid - aFrom, aTo - aMid);
//        int g0_plus_g1 = Math.max(bMid - bFrom, bTo - bMid);
//
//        for (int i = aFrom; i < aMid; i++)
//            aCache[aCacheFrom + i - aFrom] = a[i];
//        for (int i = aMid; i < aTo; ++i)
//            aCache[aCacheFrom + i - aMid] += a[i];
//
//        for (int i = bFrom; i < bMid; i++)
//            bCache[bCacheFrom + i - bFrom] = b[i];
//        for (int i = bMid; i < bTo; ++i)
//            bCache[bCacheFrom + i - bMid] += b[i];
//
//        int written = 0;
//        multiplyKaratsuba4f(aCache, aCacheFrom, aCacheFrom + f0_plus_f1, bCache, bCacheFrom, bCacheFrom + g0_plus_g1, aCache, aCacheFrom + split, bCache, bCacheFrom + split, result, resultFrom);
//
//        long[] f0g0 = multiplyKaratsuba4e(a, aFrom, aMid, b, bFrom, bMid, aCache, aCacheFrom + split, bCache, bCacheFrom + split);
//        long[] f1g1 = multiplyKaratsuba4e(a, aMid, aTo, b, bMid, bTo, aCache, aCacheFrom + split, bCache, bCacheFrom + split);
//
//        if (mid.length < f0g0.length)
//            mid = Arrays.copyOf(mid, f0g0.length);
//        if (mid.length < f1g1.length)
//            mid = Arrays.copyOf(mid, f1g1.length);
//
//        //subtract f0g0, f1g1
//        for (int i = 0; i < f0g0.length; i++)
//            mid[i] -= f0g0[i];
//
//        for (int i = 0; i < f1g1.length; i++)
//            mid[i] -= f1g1[i];
//
//
//        long[] result = Arrays.copyOf(f0g0, (aTo - aFrom) + (bTo - bFrom) - 1);
//        for (int i = 0; i < mid.length; i++)
//            result[i + split] += mid[i];
//
//        for (int i = 0; i < f1g1.length; i++)
//            result[i + 2 * split] += f1g1[i];
//
//        return result;
//    }
}
