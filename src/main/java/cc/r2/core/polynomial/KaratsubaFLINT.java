package cc.r2.core.polynomial;

import java.util.Arrays;

class KaratsubaFLINT {

    static int n_revbin(int in, int bits) {
        return Integer.reverse(in) >>> (32 - bits);
    }

    static long[] poly_mul_karatsuba(long[] poly1, long[] poly2) {
        long[] res = new long[poly1.length + poly2.length - 1];
        poly_mul_karatsuba(res, poly1, poly2);
        return res;
    }

    static void poly_mul_karatsuba(long[] res, long[] poly1, long[] poly2) {
        if (poly1.length == 0 || poly2.length == 0)
            return;

        assert res.length == poly1.length + poly2.length - 1;

        if (poly1.length >= poly2.length)
            poly_mul_karatsuba(res, poly1, poly1.length, poly2, poly2.length);
        else
            poly_mul_karatsuba(res, poly2, poly2.length, poly1, poly1.length);

    }

    private static long[] cache;

    static long[] dataCache(int length) {
        if (cache == null || cache.length < length)
            return (cache = new long[length]);
        Arrays.fill(cache, 0);
        return cache;
    }

    private static long[] tmpcache;

    static long[] tempCache(int length) {
        if (tmpcache == null || tmpcache.length < length)
            return (tmpcache = new long[length]);
        Arrays.fill(tmpcache, 0);
        return tmpcache;
    }

    static void poly_mul_karatsuba(long[] res, long[] poly1, int len1, long[] poly2, int len2) {
        assert len1 >= len2;
        if (len1 == 1) {
            res[0] = poly1[0] * poly2[0];
            return;
        }

        int loglen = 0;
        while ((1 << loglen) < len1)
            loglen++;

        int length = 1 << loglen;
        long[] data = dataCache(4 * length);//new long[4 * length];

        revbin1(data, 0, poly1, 0, len1, loglen);
        revbin1(data, length, poly2, 0, len2, loglen);

        long[] temp = tempCache(2 * length);//new long[2 * length];
        poly_mul_kara_recursive(data, 2 * length, data, 0, data, length, temp, 0, loglen);

        revbin2(res, 0, data, 2 * length, len1 + len2 - 1, loglen + 1);
    }

    static void poly_mul_kara_recursive(
            long[] out, int outShift,
            long[] rev1, int rev1Shift,
            long[] rev2, int rev2Shift,
            long[] temp, int tempShift,
            final int bits) {
        int length = 1 << bits;
        int m = length / 2;

        if (length == 2) {
            out[outShift] = rev1[rev1Shift] * rev2[rev2Shift];
            out[outShift + 2] = rev1[rev1Shift] * rev2[rev2Shift + 1] + rev1[rev1Shift + 1] * rev2[rev2Shift];
            out[outShift + 1] = rev1[rev1Shift + 1] * rev2[rev2Shift + 1];
            out[outShift + 3] = 0;
            return;
        }

        if (length == 1) {
            out[outShift] = rev1[rev1Shift] * rev2[rev2Shift];
            out[outShift + 1] = 0;
            return;
        }

        vec_add(temp, tempShift, rev1, rev1Shift, rev1, rev1Shift + m, m);
        vec_add(temp, tempShift + m, rev2, rev2Shift, rev2, rev2Shift + m, m);


        poly_mul_kara_recursive(out, outShift, rev1, rev1Shift, rev2, rev2Shift, temp, tempShift + 2 * m, bits - 1);
        poly_mul_kara_recursive(out, outShift + 2 * m, temp, tempShift, temp, tempShift + m, temp, tempShift + 2 * m, bits - 1);
        poly_mul_kara_recursive(temp, tempShift, rev1, rev1Shift + m, rev2, rev2Shift + m, temp, tempShift + 2 * m, bits - 1);

        vec_sub(out, outShift + length, out, outShift + length, out, outShift, length);
        vec_sub(out, outShift + length, out, outShift + length, temp, tempShift, length);

        vec_add_rev(out, outShift, temp, tempShift, bits);
    }

    static void vec_add(
            long[] res, int resShift,
            long[] a, int aShift,
            long[] b, int bShift,
            final int length) {
        for (int i = 0; i < length; i++)
            res[resShift + i] = a[aShift + i] + b[bShift + i];
    }

    static void vec_sub(
            long[] res, int resShift,
            long[] a, int aShift,
            long[] b, int bShift,
            final int length) {
        for (int i = 0; i < length; i++)
            res[resShift + i] = a[aShift + i] - b[bShift + i];
    }

    static void vec_add_rev(
            long[] in1, int in1Shift,
            long[] in2, int in2Shift,
            int bits) {
        for (int i = 0, to = (1 << bits) - 1; i < to; i++) {
            int j = n_revbin(n_revbin(i, bits) + 1, bits);
            in1[in1Shift + j] = in1[in1Shift + j] + in2[in2Shift + i];
        }
    }

    static void revbin1(
            long[] out, int outShift,
            long[] in, int inShift,
            int len, int bits) {
        for (int i = 0; i < len; i++)
            out[outShift + n_revbin(i, bits)] = in[inShift + i];

    }

    static void revbin2(
            long[] out, int outShift,
            long[] in, int inShift,
            int len, int bits) {
        for (int i = 0; i < len; i++)
            out[outShift + i] = in[inShift + n_revbin(i, bits)];
    }
}
