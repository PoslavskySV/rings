package cc.r2.core.poly2;

import cc.r2.core.poly2.PolynomialGCD.*;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly2.DivideAndRemainder.divideAndRemainder;
import static cc.r2.core.poly2.PolynomialGCD.*;
import static cc.r2.core.poly2.RandomPolynomials.randomPoly;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 15/02/2017.
 */
public class PolynomialGCDTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        long modulus = 11;
        MutablePolynomialMod a = MutablePolynomialMod.create(modulus, 3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950);
        MutablePolynomialMod b = MutablePolynomialMod.create(modulus, 2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216);

        PolynomialRemainders<MutablePolynomialMod> prs = PolynomialGCD.Euclid(a, b);
        MutablePolynomialMod gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test
    public void test2() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(0, 14, 50, 93, 108, 130, 70);
        MutablePolynomialZ divider = MutablePolynomialZ.create(63, 92, 143, 245, 146, 120, 90);
        MutablePolynomialZ gcd = PolynomialEuclid(dividend, divider, true).gcd();
        MutablePolynomialZ expected = MutablePolynomialZ.create(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(7, -7, 0, 1);
        MutablePolynomialZ divider = MutablePolynomialZ.create(-7, 0, 3);
        PolynomialRemainders<MutablePolynomialZ> naive = PolynomialEuclid(dividend.clone(), divider.clone(), false);
        List<MutablePolynomialZ> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(MutablePolynomialZ.create(63, -42));
        expectedNaive.add(MutablePolynomialZ.create(-1L));
        assertEquals(expectedNaive, naive.remainders);

        PolynomialRemainders<MutablePolynomialZ> primitive = PolynomialEuclid(dividend.clone(), divider.clone(), true);
        List<MutablePolynomialZ> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(MutablePolynomialZ.create(-3, 2));
        expectedPrimitive.add(MutablePolynomialZ.create(-1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        PolynomialRemainders<MutablePolynomialZ> subresultant = SubresultantEuclid(dividend.clone(), divider.clone());
        List<MutablePolynomialZ> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(MutablePolynomialZ.create(63, -42));
        expectedSubresultant.add(MutablePolynomialZ.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }

    @Test
    public void test4() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            MutablePolynomialZ dividend = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(sign);
            MutablePolynomialZ divider = MutablePolynomialZ.create(6, 7, 9, 8, 3).multiply(sign);
            PolynomialRemainders<MutablePolynomialZ> subresultant = SubresultantEuclid(dividend, divider);
            List<MutablePolynomialZ> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(MutablePolynomialZ.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(MutablePolynomialZ.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(MutablePolynomialZ.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(MutablePolynomialZ.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test4a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                MutablePolynomialZ dividend = MutablePolynomialZ.create(1, 1, 6, -3, 5).multiply(sign1);
                MutablePolynomialZ divider = MutablePolynomialZ.create(-9, -3, -10, 2, 8).multiply(sign2);
                PolynomialRemainders<MutablePolynomialZ> subresultant = SubresultantEuclid(dividend, divider);
                List<MutablePolynomialZ> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(MutablePolynomialZ.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(MutablePolynomialZ.create(4344, 3818, 9774));
                expectedSubresultant.add(MutablePolynomialZ.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(MutablePolynomialZ.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test5() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2);
        PolynomialGCD.PolynomialRemainders prs = PolynomialEuclid(dividend.clone(), dividend.clone(), false);

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutablePolynomialZ divider = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2);
        PolynomialRemainders<MutablePolynomialZ> prs;

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
    public void test7() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutablePolynomialZ divider = MutablePolynomialZ.create(7, 4, 10, 10, 0, 2);

        for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(divider, dividend)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }
    }

    @Test
    public void test8() throws Exception {
        //test long overflow
        MutablePolynomialZ dividend = MutablePolynomialZ.create(1, 1, 1, 1).multiply(2);
        MutablePolynomialZ divider = MutablePolynomialZ.create(1, 0, 2);

        for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test9() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 2, 1);
        MutablePolynomialZ expected = MutablePolynomialZ.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            MutablePolynomialZ dividend = randomPoly(5, rnd);
            MutablePolynomialZ divider = randomPoly(0, rnd);
            for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10000; i++) {
            MutablePolynomialZ dividend = randomPoly(5, 5, rnd);
            MutablePolynomialZ divider = randomPoly(5, 5, rnd);
            for (PolynomialGCD.PolynomialRemainders prs :
                    runAlgorithms(dividend, divider,
                            GCDAlgorithm.PolynomialPrimitiveEuclid,
                            GCDAlgorithm.SubresultantEuclid)) {
                assertPolynomialRemainders(dividend, divider, prs);
            }
        }
    }

    @Test
    public void testRandom3() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long[] primes = {3, 5, 7, 11, 13, 67, 97, 113, 127};
        for (int i = 0; i < 1000; i++) {
            MutablePolynomialZ dividend = randomPoly(10, 500, rnd);
            MutablePolynomialZ divider = randomPoly(10, 500, rnd);
            for (long prime : primes) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                MutablePolynomialMod a = dividend.modulus(prime);
                MutablePolynomialMod b = divider.modulus(prime);
                PolynomialRemainders euclid = Euclid(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }

    @Test
    public void test10() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        MutablePolynomialZ b = MutablePolynomialZ.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(MutablePolynomialZ.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test11() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(0, 1, 1, -6, 17, -18, 14);
        MutablePolynomialZ b = MutablePolynomialZ.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(MutablePolynomialZ.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test12() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(1, 2, 3, 3, 6, 9);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 3, 6, 5, 3);
        assertEquals(MutablePolynomialZ.create(1, 2, 3), ModularGCD(a, b));

        a = MutablePolynomialZ.create(0, 0, 1, 2);
        b = MutablePolynomialZ.create(-1, 0, 4);
        assertEquals(MutablePolynomialZ.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test13() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(-1, 0, 4);
        MutablePolynomialZ b = MutablePolynomialZ.create(-1, 0, 0, 0, 16);
        MutablePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(MutablePolynomialZ.create(-1, 0, 4), gcd);
    }

    @Test
    public void test14() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(-1, 0, 4);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, -5, 6);
        MutablePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(MutablePolynomialZ.create(-1, 2), gcd);
    }

    @Test
    public void testRandom5() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int overflow = 0;
        int larger = 0;
        DescriptiveStatistics timings = null;
        for (int k = 0; k < 2; k++) {
            timings = new DescriptiveStatistics();
            for (int i = 0; i < 10000; i++) {
                try {
                    MutablePolynomialZ a = randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    MutablePolynomialZ b = randomPoly(1 + rnd.nextInt(6), 100, rnd);
                    MutablePolynomialZ gcd = randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    MutablePolynomialZ mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    MutablePolynomialZ[] qr = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    MutablePolynomialZ[] qr1 = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
                    assertArrayEquals(qr, qr1);
                    if (!qr[0].isConstant()) ++larger;

                    assertGCD(a, b, mgcd);
                } catch (ArithmeticException e) {++overflow;}
            }
        }
        System.out.println("Overflows: " + overflow);
        System.out.println("Larger gcd: " + larger);
        System.out.println("\nTiming statistics:\n" + timings);
    }


    @Test
    public void test15() throws Exception {
        MutablePolynomialZ poly = MutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test16() throws Exception {
        MutablePolynomialZ poly = MutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.multiply(LongArithmetics.safePow(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divideOrNull(LongArithmetics.safePow(5, poly.degree));
        poly.multiply(512);
        assertEquals(-654063683625L, poly.evaluateAtRational(17, 2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void test17() throws Exception {
        MutablePolynomialZ poly = MutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }

    @Test
    public void test18() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        MutablePolynomialZ b = MutablePolynomialZ.zero();
        assertEquals(a, SubresultantEuclid(a, b).gcd());
        assertEquals(a, SubresultantEuclid(b, a).gcd());
        assertEquals(a, PolynomialEuclid(a, b, true).gcd());
        assertEquals(a, PolynomialEuclid(b, a, true).gcd());
        assertEquals(a, Euclid(a, b).gcd());
        assertEquals(a, Euclid(b, a).gcd());
        assertEquals(a, PolynomialGCD(a, b));
        assertEquals(a, PolynomialGCD(b, a));
    }

    @Test
    public void test19() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        MutablePolynomialZ b = MutablePolynomialZ.create(2);
        assertEquals(b, SubresultantEuclid(a, b).gcd());
        assertEquals(b, SubresultantEuclid(b, a).gcd());
        assertEquals(b, PolynomialEuclid(a, b, true).gcd());
        assertEquals(b, PolynomialEuclid(b, a, true).gcd());
        assertEquals(b, PolynomialEuclid(a, b, false).gcd());
        assertEquals(b, PolynomialEuclid(b, a, false).gcd());
        assertEquals(b, Euclid(a, b).gcd());
        assertEquals(b, Euclid(b, a).gcd());
        assertEquals(b, PolynomialGCD(a, b));
        assertEquals(b, PolynomialGCD(b, a));
    }

    @Test
    public void test20() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(32);
        MutablePolynomialZ b = MutablePolynomialZ.create(24);
        MutablePolynomialZ gcd = MutablePolynomialZ.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, SubresultantEuclid(a, b).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, true).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, false).gcd());
    }

    private static void assertGCD(MutablePolynomialZ a, MutablePolynomialZ b, MutablePolynomialZ gcd) {
        MutablePolynomialZ[] qr;
        for (MutablePolynomialZ dividend : Arrays.asList(a, b)) {
            qr = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(dividend, gcd, true);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutablePolynomialZ a, MutablePolynomialZ b, PolynomialRemainders<MutablePolynomialZ> prs) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        a = a.clone().divideOrNull(a.content());
        b = b.clone().divideOrNull(b.content());
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        MutablePolynomialZ gcd = prs.gcd().clone().primitivePart();
        assertTrue(DivideAndRemainder.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivideAndRemainder.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutablePolynomialMod a, MutablePolynomialMod b, PolynomialRemainders<MutablePolynomialMod> prs) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        MutablePolynomialMod gcd = prs.gcd().clone();
        assertTrue(DivideAndRemainder.divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivideAndRemainder.divideAndRemainder(b, gcd, true)[1].isZero());
    }


    private static List<PolynomialRemainders> runAlgorithms(MutablePolynomialZ dividend, MutablePolynomialZ divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        PolynomialRemainders<MutablePolynomialZ> gcd(MutablePolynomialZ dividend, MutablePolynomialZ divider) {
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

    @Test
    public void test21() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(8, 2, -1, -2, -7);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }
}