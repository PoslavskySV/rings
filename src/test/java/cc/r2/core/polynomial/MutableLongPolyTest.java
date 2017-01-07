package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * Created by poslavsky on 05/01/2017.
 */
public class MutableLongPolyTest {

    static void subtr(long[] a, int lf, int lt, int rf, int rt, int point) {
        for (int i = point, l = lf, r = rf; l < lt && r < rt; ++i, ++l, ++r)
            a[i] -= a[l] + a[r];
    }

    static long[] subtr2(long[] a, int lf, int lt, int rf, int rt, int point) {
        long[] la = a.clone();
        a = a.clone();
        for (int i = lf; i < lt; i++)
            a[point + i - lf] -= la[i];
        for (int i = rf; i < rt; i++)
            a[point + i - rf] -= la[i];
        return a;
    }


    @Test
    public void sadsd() throws Exception {
        long[] r = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int lf = 0, lt = 3, rf = 4, rt = 7, point = 1;
        System.out.println(Arrays.toString(subtr2(r, lf, lt, rf, rt, point)));
        subtr(r, lf, lt, rf, rt, point);
        System.out.println(Arrays.toString(r));

    }

    @Test
    public void name() throws Exception {
//        MutableLongPoly a = MutableLongPoly.create(1, 1, 1, 1, 1);
//        MutableLongPoly b = MutableLongPoly.create(1, 1, 1);

        MutableLongPoly a = MutableLongPoly.create(1, 1);
        MutableLongPoly b = MutableLongPoly.create(1, 1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1);

////
//        MutableLongPoly a = MutableLongPoly.create(1, 1, 1, 1, 1, 1);
//        MutableLongPoly b = MutableLongPoly.create(1, 1, 1);
//        MutableLongPoly a = MutableLongPoly.create(1, 1, 1, 1, 1, 1, 1, 1);
//        MutableLongPoly b = MutableLongPoly.create(1, 1, 1, 1, 1, 1, 1, 1);
        MutableLongPoly naive = a.clone().multiply(b.clone());
        System.out.println(naive);
        MutableLongPoly kara = a.clone().multiplyKaratsuba4a(b.clone());
        System.out.println(kara);
        System.out.println(kara.subtract(naive));
    }

    @Test
    public void sadasd() throws Exception {
        RandomGenerator rnd = new Well1024a();


        long modulus = 5659;
        int nnn = 10000;
        for (int i = 0; i < nnn; i++) {
            MutableLongPoly a = SmallPolynomialsTest.randomPoly(250 + rnd.nextInt(300), 2, rnd);
            MutableLongPoly b = SmallPolynomialsTest.randomPoly(250 + rnd.nextInt(300), 2, rnd);
            if (i > nnn - 1000)
                System.out.println("---");

            long start = System.nanoTime();
            MutableLongPoly kara = a.clone().multiplyKaratsuba4a(b.clone());
            if (i > nnn - 1000)
                System.out.println(System.nanoTime() - start);

            start = System.nanoTime();
            MutableLongPoly ord = a.clone().multiply(b.clone());
            if (i > nnn - 1000)
                System.out.println(System.nanoTime() - start);
            Assert.assertEquals(ord, kara);
        }
    }

    @Test
    public void asaasname() throws Exception {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        RandomGenerator rnd = new Well1024a();
        int r = 0;
        for (int i = 0; i < 20000; i++) {
            if (i == 10000) {
                System.out.println(stats);
                stats.clear();
            }
            MutableLongPoly a = SmallPolynomialsTest.randomPoly(300 + rnd.nextInt(300), 2, rnd);
            MutableLongPoly b = SmallPolynomialsTest.randomPoly(300 + rnd.nextInt(300), 2, rnd);
            long start = System.nanoTime();
            MutableLongPoly c = a.multiply(b);
            stats.addValue(System.nanoTime() - start);
            r += c.degree;
        }
        System.out.println(r);
        System.out.println(stats);
    }

    @Test
    public void flint1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int nnn = 1101;
        int mmm = 100;

        DescriptiveStatistics karaStats = new DescriptiveStatistics();
        DescriptiveStatistics karaStats2 = new DescriptiveStatistics();
        DescriptiveStatistics naiveStats = new DescriptiveStatistics();
        for (int i = 0; i < nnn; i++) {
            if (i == 10000) {
                karaStats.clear();
                karaStats2.clear();
                naiveStats.clear();
            }
            if (i > nnn - mmm)
                System.out.println("---");

            MutableLongPoly aP = SmallPolynomialsTest.randomPoly(5000 + rnd.nextInt(132), 3, rnd);
            MutableLongPoly bP = SmallPolynomialsTest.randomPoly(5000 + rnd.nextInt(132), 3, rnd);
            long start = System.nanoTime();
            MutableLongPoly kara = MutableLongPoly.create(Karatsuba1.multiplyKaratsuba4d(aP.data, bP.data));
            long kTime = System.nanoTime() - start;
            karaStats.addValue(kTime);


            start = System.nanoTime();
            MutableLongPoly kara2 = MutableLongPoly.create(Karatsuba1.multiplyKaratsuba4e(aP.data, bP.data));
            long kTime2 = System.nanoTime() - start;
            karaStats2.addValue(kTime2);


            start = System.nanoTime();
            MutableLongPoly j = aP.clone().multiply(bP);
            long nTime = System.nanoTime() - start;
            naiveStats.addValue(nTime);
            if (i > nnn - mmm) {
                System.out.println(kTime);
                System.out.println(kTime2);
                System.out.println(nTime);
            }

            if (!kara.subtract(j).isZero()) {
                System.out.println("KARA1");
                System.out.println(aP);
                System.out.println(bP);
                System.out.println();
                return;
            }
            if (!kara2.subtract(j).isZero()) {
                System.out.println("KARA2");
                System.out.println(aP);
                System.out.println(bP);
                System.out.println();
                return;
            }
//            System.out.println(k.subtract(aP.multiply(bP)));
        }
        System.out.println(karaStats);
        System.out.println(karaStats2);
        System.out.println(naiveStats);
    }

    @Test
    public void asd() throws Exception {
        RandomGenerator rnd = new Well1024a();
//        MutableLongPoly aP = SmallPolynomialsTest.randomPoly(300 + rnd.nextInt(432), 3, rnd);
//        MutableLongPoly bP = SmallPolynomialsTest.randomPoly(300 + rnd.nextInt(432), 3, rnd);

        MutableLongPoly aP = MutableLongPoly.create(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1);
        MutableLongPoly bP = MutableLongPoly.create(1, 1, 1, 1, 1, 1, 1);
        MutableLongPoly k = MutableLongPoly.create(KaratsubaFLINT.poly_mul_karatsuba(aP.data, bP.data));
        System.out.println(k);
        MutableLongPoly.create(Karatsuba1.multiplyKaratsuba4a(aP.data, bP.data));
//        System.out.println(karakara);
        System.out.println(aP.multiply(bP));

    }
}