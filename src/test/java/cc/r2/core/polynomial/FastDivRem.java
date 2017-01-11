package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by poslavsky on 11/01/2017.
 */
public class FastDivRem {

    private static int log2(int l) {
        return (Integer.bitCount(l) == 1 ? 32 : 33) - Integer.numberOfLeadingZeros(l);
    }

    static MutableLongPoly RemainderMonomial(MutableLongPoly poly, int xDegree) {
        return poly.cut(xDegree - 1);
    }

    static final class InverseModMonomial {
        final long modulus;
        final MutableLongPoly poly;

        InverseModMonomial(MutableLongPoly poly, long modulus) {
            if (poly.cc() != 1)
                throw new IllegalArgumentException();
            this.modulus = modulus;
            this.poly = poly;
        }

        private final ArrayList<MutableLongPoly> inverses = new ArrayList<>();

        MutableLongPoly getInverse(int xDegree) {
            if (xDegree < 1)
                return null;
            int r = log2(xDegree);
            if (inverses.size() > r)
                return inverses.get(r - 1);
            int currentSize = inverses.size();
            MutableLongPoly gPrev = currentSize == 0 ? MutableLongPoly.one() : inverses.get(inverses.size() - 1);
            for (int i = currentSize; i < r; ++i) {
                MutableLongPoly tmp = gPrev.clone().multiply(2, modulus).subtract(gPrev.clone().square(modulus).multiply(poly, modulus), modulus);
                inverses.add(gPrev = RemainderMonomial(tmp, 1 << i));
            }
            return gPrev;
        }
    }

    static MutableLongPoly InverseModMonomial0(MutableLongPoly poly, int xDegree, long modulus) {
        if (xDegree < 1)
            return null;
        if (poly.cc() != 1)
            throw new IllegalArgumentException();
        int r = log2(xDegree);
        MutableLongPoly gPrev = MutableLongPoly.one();
        for (int i = 0; i < r; ++i) {
            MutableLongPoly tmp = gPrev.clone().multiply(2, modulus).subtract(gPrev.square(modulus).multiply(poly, modulus), modulus);
            gPrev = RemainderMonomial(tmp, 1 << i);
        }
        return gPrev;
    }

    @Test
    public void testInverseModStructure() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long modulus = 101;
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly p = RandomPolynomials.randomMonicPoly(2 + rnd.nextInt(100), modulus, rnd);
            p.data[0] = 1;

            InverseModMonomial invMod = new InverseModMonomial(p, modulus);
            for (int j = 0; j < 30; j++) {
                int xDegree = 1 + rnd.nextInt(1025);
                Assert.assertEquals(invMod.getInverse(xDegree), InverseModMonomial0(p, xDegree, modulus));
            }
        }
    }

    static MutableLongPoly[] FastDivideAndRemainder(MutableLongPoly a, MutableLongPoly b, MutableLongPoly invRevMod, long modulus) {
        if (a.degree < b.degree)
            return new MutableLongPoly[]{MutableLongPoly.zero(), a};
        int m = a.degree - b.degree;
        MutableLongPoly q = RemainderMonomial(a.clone().reverse().multiply(invRevMod, modulus), m + 1).reverse();
        if (q.degree < m)
            q.shiftRight(m - q.degree);
        return new MutableLongPoly[]{q, a.clone().subtract(b.clone().multiply(q, modulus), modulus)};
    }

    static MutableLongPoly[] FastDivideAndRemainder(MutableLongPoly a, MutableLongPoly b, InverseModMonomial invRevMod, long modulus) {
        if (a.degree < b.degree)
            return new MutableLongPoly[]{MutableLongPoly.zero(), a};
        int m = a.degree - b.degree;
        MutableLongPoly q = RemainderMonomial(a.clone().reverse().multiply(invRevMod.getInverse(m + 1), modulus), m + 1).reverse();
        if (q.degree < m)
            q.shiftRight(m - q.degree);
        return new MutableLongPoly[]{q, a.clone().subtract(b.clone().multiply(q, modulus), modulus)};
    }

    @Test
    public void testFastRemainderRandom() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long modulus = 17;
        for (int i = 0; i < 300; i++) {
            MutableLongPoly b = RandomPolynomials.randomMonicPoly(30, modulus, rnd);
            MutableLongPoly a = RandomPolynomials.randomMonicPoly(rnd.nextInt(30), modulus, rnd);

            MutableLongPoly invMod = InverseModMonomial0(b.clone().reverse(), a.degree - b.degree + 1, modulus);
            MutableLongPoly[] fast = FastDivideAndRemainder(a, b, invMod, modulus);
            MutableLongPoly[] plain = SmallPolynomials.divideAndRemainder(a, b, modulus, true);
            if (!Arrays.equals(fast, plain)) {
                System.out.println(a.toStringForCopy());
                System.out.println(b.toStringForCopy());
            }
            Assert.assertArrayEquals(fast, plain);
        }
    }

    @Test
    public void testFastRemainder1() throws Exception {
        long modulus = 7;
        MutableLongPoly a = MutableLongPoly.create(5, 1, 4, 6, 4, 3, 5, 5, 3, 4, 2, 2, 5, 2, 5, 6, 1, 1, 2, 5, 1, 0, 0, 6, 6, 5, 5, 1, 0, 1, 4, 1, 1);
        MutableLongPoly b = MutableLongPoly.create(2, 5, 3, 1, 1, 5, 6, 3, 4, 0, 0, 5, 4, 0, 2, 1);
        MutableLongPoly invMod = InverseModMonomial0(b.clone().reverse(), a.degree - b.degree + 1, modulus);
        MutableLongPoly[] fast = FastDivideAndRemainder(a, b, invMod, modulus);
        MutableLongPoly[] plain = SmallPolynomials.divideAndRemainder(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void testFastRemainder2() throws Exception {
        long modulus = 7;
        MutableLongPoly a = MutableLongPoly.create(5, 3, 3, 3, 5, 3, 1, 4, -3, 1, 4, 5, 0, 2, 2, -5, 1).modulus(modulus);
        MutableLongPoly b = MutableLongPoly.create(0, 4, 6, 1, 2, 4, 0, 0, 6, 5, 2, 3, 1, 4, 0, 1);
        MutableLongPoly invMod = InverseModMonomial0(b.clone().reverse(), a.degree - b.degree + 1, modulus);
        MutableLongPoly[] fast = FastDivideAndRemainder(a, b, invMod, modulus);
        MutableLongPoly[] plain = SmallPolynomials.divideAndRemainder(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void testFastRemainder3() throws Exception {
        long modulus = 17;
        MutableLongPoly a = MutableLongPoly.create(0, 6, 2, 1, 10, 15, 16, 15, 2, 11, 13, 0, 1, 15, 5, 13, 8, 14, 13, 14, 15, 1, 1);
        MutableLongPoly b = MutableLongPoly.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1);
        MutableLongPoly invMod = InverseModMonomial0(b.clone().reverse(), a.degree - b.degree + 1, modulus);
        MutableLongPoly[] fast = FastDivideAndRemainder(a, b, invMod, modulus);
        MutableLongPoly[] plain = SmallPolynomials.divideAndRemainder(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void testFastRemainder4() throws Exception {
        long modulus = 17;
        MutableLongPoly a = MutableLongPoly.create(5, 9, 4, 9, 8, 12, 11, 9, 1, 6, 15, 7, 11, 2, 11, 13, 11, 10, 5, 1);
        MutableLongPoly b = MutableLongPoly.create(11, 15, 9, 5, 11, 5, 14, 9, 1, 0, 16, 12, 11, 5, 15, 10, 15, 2, 14, 3, 1, 16, 16, 12, 13, 1, 12, 11, 1, 15, 1);
        MutableLongPoly invMod = InverseModMonomial0(b.clone().reverse(), a.degree - b.degree + 1, modulus);
        MutableLongPoly[] fast = FastDivideAndRemainder(a, b, invMod, modulus);
        MutableLongPoly[] plain = SmallPolynomials.divideAndRemainder(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void testInverseModRandom() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int modulus = 17;
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly f = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(100), modulus, rnd);
            f.data[0] = 1;
            int modDegree = 1 + rnd.nextInt(2 * f.degree);
            MutableLongPoly invmod = InverseModMonomial0(f, modDegree, modulus);
            assertInverseModMonomial(f, invmod, modDegree, modulus);
        }
    }

    static void assertInverseModMonomial(MutableLongPoly poly, MutableLongPoly invMod, int monomialDegree, long modulus) {
        Assert.assertTrue(SmallPolynomialArithmetics.polyMultiplyMod(poly, invMod, MutableLongPoly.createMonomial(1, monomialDegree), modulus, true).isOne());
    }

    @Test
    public void testInverseMod1() throws Exception {
        long modulus = 7;
        MutableLongPoly f = MutableLongPoly.create(0, 2, 3, 4, -5, 1).reverse();
        int modDegree = f.degree;
        MutableLongPoly invmod = InverseModMonomial0(f, modDegree, modulus);
        MutableLongPoly r = SmallPolynomialArithmetics.polyMultiplyMod(f, invmod, MutableLongPoly.createMonomial(1, modDegree), modulus, true);
        Assert.assertTrue(r.isOne());
    }

    @Test
    public void testInverseMod2() throws Exception {
        long modulus = 17;
        MutableLongPoly f = MutableLongPoly.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1).reverse();
        int modDegree = 9;
        MutableLongPoly invmod = InverseModMonomial0(f, modDegree, modulus);
        MutableLongPoly r = SmallPolynomialArithmetics.polyMultiplyMod(f, invmod, MutableLongPoly.createMonomial(1, modDegree), modulus, true);
        Assert.assertTrue(r.isOne());
    }


    @Test
    public void testPerformance() throws Exception {
        long modulus = LargeDDFTest.bigModulus;
        RandomGenerator rnd = new Well1024a();
        MutableLongPoly polyModulus = RandomPolynomials.randomPoly(300, (int) modulus, rnd);
        polyModulus.monic(modulus);

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        InverseModMonomial invRev = new InverseModMonomial(polyModulus.clone().reverse(), modulus);
        int nIterations = 1500;
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000) {
                classic.clear();
                fast.clear();
            }
            MutableLongPoly a = RandomPolynomials.randomPoly(2 * polyModulus.degree + rnd.nextInt(100), (int) modulus, rnd);
            a = a.modulus(modulus);

            long start = System.nanoTime();
            MutableLongPoly[] qdPlain = SmallPolynomials.divideAndRemainder(a, polyModulus, modulus, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            MutableLongPoly[] qdNewton = FastDivideAndRemainder(a, polyModulus, invRev, modulus);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            if (i > nIterations - 10) {
                System.out.println("====");
                System.out.println(plain);
                System.out.println(newton);
            }
            Assert.assertArrayEquals(qdPlain, qdNewton);
        }

        System.out.println("==== Plain ====");
        System.out.println(classic);

        System.out.println("==== Fast ====");
        System.out.println(fast);

    }
}
