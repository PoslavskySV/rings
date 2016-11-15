package cc.r2.core.polynomial;

import cc.r2.core.number.*;
import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.Well512a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.number.BigIntegerRing.IntegerRing;
import static cc.r2.core.number.BigRationalField.BigRationalField;

public class PolynomialsTest {

    static UnivariatePolynomial<BigRational> makeQx(int... cfx) {
        return new UnivariatePolynomial<>(BigRationalField, Arrays.stream(cfx).mapToObj(BigRational::new).toArray(BigRational[]::new));
    }

    static UnivariatePolynomial<BigInteger> makeZx(int... cfx) {
        return new UnivariatePolynomial<>(IntegerRing, Arrays.stream(cfx).mapToObj(BigInteger::valueOf).toArray(BigInteger[]::new));
    }

    static UnivariatePolynomial<ModPrimeBigInteger> makeZpx(ModPrimeBigIntegerField field, int... cfx) {
        return new UnivariatePolynomial<>(field, Arrays.stream(cfx).mapToObj(x -> new ModPrimeBigInteger(field, BigInteger.valueOf(x))).toArray(ModPrimeBigInteger[]::new));
    }


    @Test
    public void pDivideAndRemainder() throws Exception {
        UnivariatePolynomial<BigInteger> a = makeZx(7, -7, 0, 0, 1);
        UnivariatePolynomial<BigInteger> b = makeZx(-7, 0, 3);

        UnivariatePolynomial<BigInteger>[] qr = Polynomials.pseudoDivideAndRemainder(IntegerRing, a, b);
        Assert.assertEquals(makeZx(21, 0, 9), qr[0]);
        Assert.assertEquals(makeZx(336, -189), qr[1]);
    }

    @Test
    public void test1() throws Exception {
        UnivariatePolynomial<BigRational> a = makeQx(5, 1, 0, 2, 3);
        UnivariatePolynomial<BigRational> b = makeQx(3, 2, 1);

        UnivariatePolynomial<BigRational>[] qr = Polynomials.divideAndRemainder(BigRationalField, a, b);
        UnivariatePolynomial<BigRational> q = qr[0], r = qr[1];


        Assert.assertEquals(makeQx(-1, -4, 3), q);
        Assert.assertEquals(makeQx(8, 15), r);
    }

    @Test
    public void test1a() throws Exception {
        UnivariatePolynomial<BigInteger> a = makeZx(5, 1, 0, 2, 3);
        UnivariatePolynomial<BigInteger> b = makeZx(3, 2, 1);

        UnivariatePolynomial<BigInteger>[] qr = Polynomials.divideAndRemainder(IntegerRing, a, b);
        UnivariatePolynomial<BigInteger> q = qr[0], r = qr[1];


        Assert.assertEquals(makeZx(-1, -4, 3), q);
        Assert.assertEquals(makeZx(8, 15), r);
    }

    @Test
    public void test1c() throws Exception {
        ModPrimeBigIntegerField fi = new ModPrimeBigIntegerField(BigInteger.valueOf(113));
        UnivariatePolynomial<ModPrimeBigInteger> a = makeZpx(fi, 5, 1, 0, 2, 3);
        UnivariatePolynomial<ModPrimeBigInteger> b = makeZpx(fi, 3, 2, 1);

        UnivariatePolynomial<ModPrimeBigInteger>[] qr = Polynomials.divideAndRemainder(fi, a, b);
        UnivariatePolynomial<ModPrimeBigInteger> q = qr[0], r = qr[1];

        Assert.assertEquals(makeZpx(fi, 112, 109, 3), q);
        Assert.assertEquals(makeZpx(fi, 8, 15), r);
    }

    @Test
    public void test1b() throws Exception {
        UnivariatePolynomial<BigInteger> a = makeZx(5, 1, 0, 2, 1);
        UnivariatePolynomial<BigInteger> b = makeZx(3, 2, 1);

        UnivariatePolynomial<BigInteger>[] qr = Polynomials.divideAndRemainder(IntegerRing, a, b);
        UnivariatePolynomial<BigInteger> q = qr[0], r = qr[1];
        System.out.println(q);
        System.out.println(r);
    }

    @Test
    public void performanceTest1() throws Exception {
        TDoubleArrayList l = new TDoubleArrayList();
        for (int m = 10; m < 150; m++) {
            if (m % 5 == 0) System.out.println(m);
            l.add(performanceDivideAndRemainder(BigRationalField, m + 10, m, 10).getMean() / 1000);
        }

        int Up = 2000;
        l.clear();
        for (int m = 10; m < Up; m++) {
            if (m % 5 == 0) System.out.println(m);
            l.add(performanceDivideAndRemainder(BigRationalField, m + 10, m, 2).getMean() / 1000);
        }


        int i = 0;
        System.out.print("{");
        for (int m = 10; m < Up; m++) {
            System.out.print("{" + m + "," + l.get(i++) + "},");
        }
    }

    @Test
    public void performanceTest2() throws Exception {
        ModPrimeBigIntegerField Zp = new ModPrimeBigIntegerField(BigInteger.valueOf(224199919));
        for (int i = 0; i < 100; i++) {
            DescriptiveStatistics zp = performanceDivideAndRemainder(Zp, 50, 30, 100);
            DescriptiveStatistics q = performanceDivideAndRemainder(BigRationalField, 50, 30, 100);
            System.out.println();
            System.out.println("Zp: " + (zp.getMean()));
            System.out.println("Q:  " + (q.getMean()));
        }
    }

    static <R extends EuclideanRingElement<R>> DescriptiveStatistics performanceDivideAndRemainder(
            Ring<R> ring, int n, int m, int N) {
        DescriptiveStatistics stats = new DescriptiveStatistics();
        BigInteger minCoeff = new BigInteger("12414324");
        BigInteger maxCoeff = new BigInteger("32414324");
        RandomDataGenerator rnd = new RandomDataGenerator(new Well512a(System.currentTimeMillis()));


        for (int i = 0; i < N; i++) {

            UnivariatePolynomial<R> a = randomPoly(ring, n, minCoeff, maxCoeff, rnd);
            UnivariatePolynomial<R> b = randomPoly(ring, m, minCoeff, maxCoeff, rnd);

            long start = System.currentTimeMillis();
            UnivariatePolynomial<R>[] qr = Polynomials.divideAndRemainder(ring, a, b);
            stats.addValue((System.currentTimeMillis() - start));

            UnivariatePolynomial<R> q = qr[0], r = qr[1];
            Assert.assertEquals(a, b.multiply(q).add(r));
        }
        return stats;
    }

    @Test
    public void testRandom1() throws Exception {
        RandomDataGenerator rnd = new RandomDataGenerator(new Well512a());

        BigInteger minCoeff = new BigInteger("12414324");
        BigInteger maxCoeff = new BigInteger("32414324");

        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (int i = 0; i < 10000; i++) {
            if (i % 10 == 0) System.out.println(i);
            UnivariatePolynomial<BigRational> a = randomPoly(BigRationalField, 30, minCoeff, maxCoeff, rnd);
            UnivariatePolynomial<BigRational> b = randomPoly(BigRationalField, 20, minCoeff, maxCoeff, rnd);

            long start = System.nanoTime();
            UnivariatePolynomial<BigRational>[] qr = Polynomials.divideAndRemainder(BigRationalField, a, b);
            UnivariatePolynomial<BigRational> q = qr[0], r = qr[1];
            stats.addValue((System.nanoTime() - start));

//            System.out.println(a);
//            System.out.println(b);
//            System.out.println(q);
//            System.out.println(r);
//            System.out.println(q.degree());
//            System.out.println(r.degree());
//            System.out.println();
            Assert.assertEquals(a, b.multiply(q).add(r));
        }
//
        System.out.println(stats);
    }

    @SuppressWarnings("unchecked")
    static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R> randomPoly(
            Ring<R> ring, int degree, BigInteger minCoeff, BigInteger maxCoeff, RandomDataGenerator rnd) {
        R[] cfx = (R[]) new EuclideanRingElement[degree];
        long delta = maxCoeff.subtract(minCoeff).longValueExact();
        for (int i = 0; i < degree; i++) {
            cfx[i] = ring.parse(minCoeff.add(BigInteger.valueOf(rnd.nextLong(0, delta))));
            if (rnd.getRandomGenerator().nextBoolean() && rnd.getRandomGenerator().nextBoolean())
                cfx[i] = cfx[i].negate();
        }
        return new UnivariatePolynomial<>(ring, cfx);
    }


    @Test
    public void prs1() throws Exception {
        UnivariatePolynomial<BigRational> p1 = makeQx(7, -7, 0, 1);
        UnivariatePolynomial<BigRational> p2 = makeQx(-7, 0, 3);

        List<UnivariatePolynomial<BigRational>> prs = Polynomials.EuclideanPRS(BigRationalField, p1, p2);
        prs.forEach(System.out::println);
    }

    @Test
    public void prs123() throws Exception {
        UnivariatePolynomial<BigRational> p1 = makeQx(7, -7, 0, 1);
        UnivariatePolynomial<BigRational> p2 = makeQx(-7, 0, 3);

        List<UnivariatePolynomial<BigRational>> prs1 = Polynomials.EuclideanPRS(BigRationalField, p1, p2);
        prs1.forEach(System.out::println);
        System.out.println("\n--------\n");
        List<UnivariatePolynomial<BigRational>> prs2 = Polynomials.naivePRS(BigRationalField, p1, p2);
        prs2.forEach(System.out::println);
    }

    @Test
    public void test12312123() throws Exception {
        UnivariatePolynomial<BigInteger> p1 = makeZx(-5, 2, 8, -3, -3, 0, 1, 0, 1);
        UnivariatePolynomial<BigInteger> p2 = makeZx(21, -9, -4, 0, 5, 0, 3);
        List<UnivariatePolynomial<BigInteger>> prs = Polynomials.subResultantPRS(IntegerRing, p1, p2);
        prs.forEach(System.out::println);
    }

    @Test
    public void gcdTimings() throws Exception {
        RandomDataGenerator rnd = new RandomDataGenerator(new Well512a());

        List<long[]> data = new ArrayList<>();
        for (int i = 5; i < 100; i++) {
            long l = gcdTime(i, 3, rnd);
            System.out.println(i + " : " + l);
            data.add(new long[]{i, l});
        }

        System.out.print("{");
        for (long[] tt : data) {
            System.out.print("{" + tt[0] + "," + tt[1] + "},");
        }
    }

    @Test
    public void sasadasdasdas() throws Exception {
        for (int i = 0; i < 1000; i++) {
            System.out.println( new String(new char[i]));
        }

//
//        for (int i = 0; i < 1000000; i++) {
//            if(i%1000==0)System.out.println(i);
//            Assert.assertEquals(BigInteger.valueOf(i).isPrime(), isPrime(i));
//        }
    }

    public static boolean isPrime(int n) {
        if (n % 2 == 0)
            return false;
        return !new String(new char[n]).matches(".?|(..+?)\\1+");
    }

    static long gcdTime(int n, int tr, RandomDataGenerator rnd) {
        BigInteger minCoeff = BigInteger.valueOf(10);
        BigInteger maxCoeff = BigInteger.valueOf(1000000);


        long tot = 0;
        for (int i = 0; i < tr; i++) {
            UnivariatePolynomial<BigInteger> a = randomPoly(IntegerRing, n, minCoeff, maxCoeff, rnd);
            UnivariatePolynomial<BigInteger> b = randomPoly(IntegerRing, n - 2, minCoeff, maxCoeff, rnd);
            UnivariatePolynomial<BigInteger> common = randomPoly(IntegerRing, n / 3, minCoeff, maxCoeff, rnd);
            a = a.multiply(common);
            b = b.multiply(common);
            long start = System.nanoTime();
            List<UnivariatePolynomial<BigInteger>> t = Polynomials.subResultantPRS(IntegerRing, a, b);
            tot += System.nanoTime() - start;
        }
        tot /= tr;
        tot /= 1000_000;
        return tot;
    }

    @Test
    public void test123123() throws Exception {
        UnivariatePolynomial<BigInteger> p1 = makeZx(1, 2, 3, 4, 23, 12, 423, 2, 42, 523, 12, 321, 312, 3423, 425, 45, 23, 12, 131, 2423, 432, 412, 3123, 12);
        UnivariatePolynomial<BigInteger> p2 = makeZx(312, 312, 3, 123, 124, 546, 88, 2, 3, 6, 7, 8, 7, 878, 65, 7, 64435, 3, 4534, 53, 53, 65445, 1, 321);
        UnivariatePolynomial<BigInteger> gcd = makeZx(1, 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1);
        p1 = p1.multiply(gcd);
        p2 = p2.multiply(gcd);

        UnivariatePolynomial<BigInteger>[] qr = Polynomials.pseudoDivideAndRemainder(IntegerRing, p1, p2);
//        System.out.println(qr[0]);
//        System.out.println(qr[1]);

//        System.out.println("XXX");


        for (int i = 0; i < 100; i++) {
            System.out.println();
            List<UnivariatePolynomial<BigInteger>> prs;

            long start = System.currentTimeMillis();
            prs = Polynomials.primitivePRS(IntegerRing, p1, p2);
            Assert.assertEquals(gcd, prs.get(prs.size() - 1));
            System.out.println(System.currentTimeMillis() - start);

            start = System.currentTimeMillis();
            prs = Polynomials.subResultantPRS(IntegerRing, p1, p2);
            System.out.println(System.currentTimeMillis() - start);

            Assert.assertEquals(gcd, Polynomials.primitivePart(prs.get(prs.size() - 1)));
        }


        List<UnivariatePolynomial<BigInteger>> prs = Polynomials.subResultantPRS(IntegerRing, p1, p2);
        System.out.println(gcd);
        System.out.println(Polynomials.primitivePart(prs.get(prs.size() - 1)));
//        System.out.println();
//        prs.forEach(System.out::println);
//        System.out.println(prs.get(prs.size() - 1));
//        System.out.println(prs.get(prs.size() - 2));
//        System.out.println(prs.get(prs.size() - 3));
//        System.out.println(prs.get(prs.size() - 4));
//        System.out.println(prs);
    }
}