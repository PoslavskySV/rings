package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static cc.r2.core.polynomial.SmallPolynomials.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class SmallPolynomialsTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950);
        MutableLongPoly b = MutableLongPoly.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216);

        long modulus = 11;
        SmallPolynomials.PolynomialRemainders prs = SmallPolynomials.Euclid(a, b, modulus);
        MutableLongPoly gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, modulus)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, modulus)[1].isZero());
    }

    @Test(expected = AssertionError.class)
    public void test2() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(new long[]{28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575, 0}, 7);
        MutableLongPoly divider = new MutableLongPoly(new long[]{24487310, 38204421, 12930314, 41553770, -1216266, 7382581, 15631547, 0, 0}, 6);
        pseudoDivideAndRemainder(dividend, divider);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(0, 14, 50, 93, 108, 130, 70);
        MutableLongPoly divider = MutableLongPoly.create(63, 92, 143, 245, 146, 120, 90);
        MutableLongPoly gcd = PolynomialEuclid(dividend, divider, true).gcd();
        MutableLongPoly expected = MutableLongPoly.create(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test4() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(7, -7, 0, 1);
        MutableLongPoly divider = MutableLongPoly.create(-7, 0, 3);
        SmallPolynomials.PolynomialRemainders naive = PolynomialEuclid(dividend.clone(), divider.clone(), false);
        List<MutableLongPoly> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(MutableLongPoly.create(63, -42));
        expectedNaive.add(MutableLongPoly.create(-1L));
        assertEquals(expectedNaive, naive.remainders);

        SmallPolynomials.PolynomialRemainders primitive = PolynomialEuclid(dividend.clone(), divider.clone(), true);
        List<MutableLongPoly> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(MutableLongPoly.create(-3, 2));
        expectedPrimitive.add(MutableLongPoly.create(-1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        SmallPolynomials.PolynomialRemainders subresultant = SubresultantEuclid(dividend.clone(), divider.clone());
        List<MutableLongPoly> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(MutableLongPoly.create(63, -42));
        expectedSubresultant.add(MutableLongPoly.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }


    @Test
    public void test5() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            MutableLongPoly dividend = MutableLongPoly.create(7, 4, 10, 10, 0, 2).multiply(sign);
            MutableLongPoly divider = MutableLongPoly.create(6, 7, 9, 8, 3).multiply(sign);
            SmallPolynomials.PolynomialRemainders subresultant = SubresultantEuclid(dividend, divider);
            List<MutableLongPoly> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(MutableLongPoly.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(MutableLongPoly.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(MutableLongPoly.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(MutableLongPoly.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test5a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                MutableLongPoly dividend = MutableLongPoly.create(1, 1, 6, -3, 5).multiply(sign1);
                MutableLongPoly divider = MutableLongPoly.create(-9, -3, -10, 2, 8).multiply(sign2);
                SmallPolynomials.PolynomialRemainders subresultant = SubresultantEuclid(dividend, divider);
                List<MutableLongPoly> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(MutableLongPoly.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(MutableLongPoly.create(4344, 3818, 9774));
                expectedSubresultant.add(MutableLongPoly.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(MutableLongPoly.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(7, 4, 10, 10, 0, 2);
        SmallPolynomials.PolynomialRemainders prs = PolynomialEuclid(dividend.clone(), dividend.clone(), false);

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test7() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutableLongPoly divider = MutableLongPoly.create(7, 4, 10, 10, 0, 2);
        SmallPolynomials.PolynomialRemainders prs;

        for (boolean prim : new boolean[]{true, false}) {
            prs = PolynomialEuclid(dividend.clone(), divider.clone(), prim);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());

            prs = PolynomialEuclid(divider.clone(), dividend.clone(), prim);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        prs = SubresultantEuclid(dividend.clone(), divider.clone());
        assertEquals(2, prs.remainders.size());
        assertEquals(divider, prs.gcd());
    }

    @Test
    public void test8() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutableLongPoly divider = MutableLongPoly.create(7, 4, 10, 10, 0, 2);

        for (SmallPolynomials.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        for (SmallPolynomials.PolynomialRemainders prs : runAlgorithms(divider, dividend)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }
    }

    @Test
    public void test9() throws Exception {
        //test long overflow
        MutableLongPoly dividend = MutableLongPoly.create(1, 1, 1, 1).multiply(2);
        MutableLongPoly divider = MutableLongPoly.create(1, 0, 2);

        for (SmallPolynomials.PolynomialRemainders prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test10() throws Exception {
        MutableLongPoly dividend = MutableLongPoly.create(28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575);
        MutableLongPoly divider = MutableLongPoly.one();
        MutableLongPoly[] qr = divideAndRemainder(dividend, divider);
        assertEquals(dividend, qr[0]);
        assertTrue(qr[1].isZero());
    }

    @Test
    public void test11() throws Exception {
        RandomGenerator rnd = new Well1024a();
        MutableLongPoly[] qd;
        long[] primes = {3, 5, 7, 11, 13, 67, 97, 113, 127};
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly dividend = randomPoly(rnd.nextInt(20), rnd);
            MutableLongPoly divider = randomPoly(rnd.nextInt(20), rnd);

            for (long prime : primes) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                qd = divideAndRemainder(dividend, divider, prime);
                assertQuotientRemainder(dividend, divider, qd, prime);
            }
        }
    }

    @Test
    public void test12() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        MutableLongPoly b = MutableLongPoly.create(1, 2, 1);
        MutableLongPoly expected = MutableLongPoly.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            MutableLongPoly dividend = randomPoly(5, rnd);
            MutableLongPoly divider = randomPoly(0, rnd);
            for (SmallPolynomials.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10000; i++) {
            MutableLongPoly dividend = randomPoly(5, 5, rnd);
            MutableLongPoly divider = randomPoly(5, 5, rnd);
            for (SmallPolynomials.PolynomialRemainders prs :
                    runAlgorithms(dividend, divider,
                            GCDAlgorithm.PolynomialPrimitiveEuclid,
                            GCDAlgorithm.SubresultantEuclid)) {
                assertPolynomialRemainders(dividend, divider, prs);
            }
        }
    }

    private static void assertQuotientRemainder(MutableLongPoly dividend, MutableLongPoly divider, MutableLongPoly[] qr) {
        if (qr == null) return;
        assertEquals(dividend, divider.clone().multiply(qr[0]).add(qr[1]));
    }

    private static void assertQuotientRemainder(MutableLongPoly dividend, MutableLongPoly divider, MutableLongPoly[] qr, long modulus) {
        if (qr == null) return;
        assertEquals(dividend.clone().modulus(modulus), divider.clone().multiply(qr[0], modulus).add(qr[1], modulus));
    }

    private static void assertPolynomialRemainders(MutableLongPoly a, MutableLongPoly b, PolynomialRemainders prs) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        a = a.clone().divide(content(a));
        b = b.clone().divide(content(b));
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        MutableLongPoly gcd = primitivePart(prs.gcd().clone());
        assertTrue(pseudoDivideAndRemainder(a, gcd)[1].isZero());
        assertTrue(pseudoDivideAndRemainder(b, gcd)[1].isZero());
    }

    private static List<PolynomialRemainders> runAlgorithms(MutableLongPoly dividend, MutableLongPoly divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        PolynomialRemainders gcd(MutableLongPoly dividend, MutableLongPoly divider) {
            switch (this) {
                case PolynomialEuclid:
                    return PolynomialEuclid(dividend.clone(), divider.clone(), false);
                case PolynomialPrimitiveEuclid:
                    return PolynomialEuclid(dividend.clone(), divider.clone(), true);
                case SubresultantEuclid:
                    return SubresultantEuclid(dividend.clone(), divider.clone());
            }
            throw new IllegalArgumentException();
        }
    }

    private static final int CBOUND = 100;

    public static MutableLongPoly randomPoly(int degree, RandomGenerator rnd) {
        return randomPoly(degree, CBOUND, rnd);
    }

    public static MutableLongPoly randomPoly(int degree, int bound, RandomGenerator rnd) {
        long[] data = new long[degree + 1];
        for (int i = 0; i <= degree; ++i) {
            data[i] = rnd.nextInt(bound);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                data[i] = -data[i];
        }
        while (data[degree] == 0)
            data[degree] = rnd.nextInt(bound);
        return new MutableLongPoly(data, degree);
    }
}