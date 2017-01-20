package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.polynomial.DivideAndRemainder.divideAndRemainder;
import static cc.r2.core.polynomial.PolynomialGCD.*;
import static cc.r2.core.polynomial.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.polynomial.RandomPolynomials.randomPoly;
import static org.junit.Assert.*;

public class PolynomialGCDTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950);
        MutablePolynomial b = MutablePolynomial.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216);

        long modulus = 11;
        PolynomialGCD.PolynomialRemainders prs = PolynomialGCD.Euclid(a, b, modulus);
        MutablePolynomial gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, modulus, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, modulus, true)[1].isZero());
    }

    @Test
    public void test2() throws Exception {
        //test long overflow
        MutablePolynomial dividend = MutablePolynomial.create(0, 14, 50, 93, 108, 130, 70);
        MutablePolynomial divider = MutablePolynomial.create(63, 92, 143, 245, 146, 120, 90);
        MutablePolynomial gcd = PolynomialEuclid(dividend, divider, true).gcd();
        MutablePolynomial expected = MutablePolynomial.create(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        MutablePolynomial dividend = MutablePolynomial.create(7, -7, 0, 1);
        MutablePolynomial divider = MutablePolynomial.create(-7, 0, 3);
        PolynomialGCD.PolynomialRemainders naive = PolynomialEuclid(dividend.clone(), divider.clone(), false);
        List<MutablePolynomial> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(MutablePolynomial.create(63, -42));
        expectedNaive.add(MutablePolynomial.create(-1L));
        assertEquals(expectedNaive, naive.remainders);

        PolynomialGCD.PolynomialRemainders primitive = PolynomialEuclid(dividend.clone(), divider.clone(), true);
        List<MutablePolynomial> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(MutablePolynomial.create(-3, 2));
        expectedPrimitive.add(MutablePolynomial.create(-1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        PolynomialGCD.PolynomialRemainders subresultant = SubresultantEuclid(dividend.clone(), divider.clone());
        List<MutablePolynomial> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(MutablePolynomial.create(63, -42));
        expectedSubresultant.add(MutablePolynomial.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }

    @Test
    public void test4() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            MutablePolynomial dividend = MutablePolynomial.create(7, 4, 10, 10, 0, 2).multiply(sign);
            MutablePolynomial divider = MutablePolynomial.create(6, 7, 9, 8, 3).multiply(sign);
            PolynomialGCD.PolynomialRemainders subresultant = SubresultantEuclid(dividend, divider);
            List<MutablePolynomial> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(MutablePolynomial.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(MutablePolynomial.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(MutablePolynomial.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(MutablePolynomial.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test4a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                MutablePolynomial dividend = MutablePolynomial.create(1, 1, 6, -3, 5).multiply(sign1);
                MutablePolynomial divider = MutablePolynomial.create(-9, -3, -10, 2, 8).multiply(sign2);
                PolynomialGCD.PolynomialRemainders subresultant = SubresultantEuclid(dividend, divider);
                List<MutablePolynomial> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(MutablePolynomial.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(MutablePolynomial.create(4344, 3818, 9774));
                expectedSubresultant.add(MutablePolynomial.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(MutablePolynomial.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test5() throws Exception {
        //test long overflow
        MutablePolynomial dividend = MutablePolynomial.create(7, 4, 10, 10, 0, 2);
        PolynomialGCD.PolynomialRemainders prs = PolynomialEuclid(dividend.clone(), dividend.clone(), false);

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        MutablePolynomial dividend = MutablePolynomial.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutablePolynomial divider = MutablePolynomial.create(7, 4, 10, 10, 0, 2);
        PolynomialGCD.PolynomialRemainders prs;

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
        MutablePolynomial dividend = MutablePolynomial.create(7, 4, 10, 10, 0, 2).multiply(2);
        MutablePolynomial divider = MutablePolynomial.create(7, 4, 10, 10, 0, 2);

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
        MutablePolynomial dividend = MutablePolynomial.create(1, 1, 1, 1).multiply(2);
        MutablePolynomial divider = MutablePolynomial.create(1, 0, 2);

        for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test9() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        MutablePolynomial b = MutablePolynomial.create(1, 2, 1);
        MutablePolynomial expected = MutablePolynomial.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 100; i++) {
            MutablePolynomial dividend = randomPoly(5, rnd);
            MutablePolynomial divider = randomPoly(0, rnd);
            for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 10000; i++) {
            MutablePolynomial dividend = randomPoly(5, 5, rnd);
            MutablePolynomial divider = randomPoly(5, 5, rnd);
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
            MutablePolynomial dividend = randomPoly(10, 500, rnd);
            MutablePolynomial divider = randomPoly(10, 500, rnd);
            for (long prime : primes) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                PolynomialRemainders euclid = Euclid(dividend, divider, prime);
                assertPolynomialRemainders(dividend, divider, euclid, prime);
            }
        }
    }

    @Test
    public void test10() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        MutablePolynomial b = MutablePolynomial.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(MutablePolynomial.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test11() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(0, 1, 1, -6, 17, -18, 14);
        MutablePolynomial b = MutablePolynomial.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(MutablePolynomial.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test12() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 2, 3, 3, 6, 9);
        MutablePolynomial b = MutablePolynomial.create(1, 3, 6, 5, 3);
        assertEquals(MutablePolynomial.create(1, 2, 3), ModularGCD(a, b));

        a = MutablePolynomial.create(0, 0, 1, 2);
        b = MutablePolynomial.create(-1, 0, 4);
        assertEquals(MutablePolynomial.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test13() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(-1, 0, 4);
        MutablePolynomial b = MutablePolynomial.create(-1, 0, 0, 0, 16);
        MutablePolynomial gcd = ModularGCD(a, b);
        assertEquals(MutablePolynomial.create(-1, 0, 4), gcd);
    }

    @Test
    public void test14() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(-1, 0, 4);
        MutablePolynomial b = MutablePolynomial.create(1, -5, 6);
        MutablePolynomial gcd = ModularGCD(a, b);
        assertEquals(MutablePolynomial.create(-1, 2), gcd);
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
                    MutablePolynomial a = randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    MutablePolynomial b = randomPoly(1 + rnd.nextInt(6), 100, rnd);
                    MutablePolynomial gcd = randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    MutablePolynomial mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    MutablePolynomial[] qr = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    MutablePolynomial[] qr1 = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test16() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.multiply(LongArithmetics.pow(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divide(LongArithmetics.pow(5, poly.degree));
        poly.multiply(512);
        assertEquals(-654063683625L, poly.evaluateAtRational(17, 2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void test17() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }


    @Test
    public void test18() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(8, -2 * 8, 8, 8 * 2);
        MutablePolynomial b = MutablePolynomial.zero();
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
        MutablePolynomial a = MutablePolynomial.create(8, -2 * 8, 8, 8 * 2);
        MutablePolynomial b = MutablePolynomial.create(2);
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
        MutablePolynomial a = MutablePolynomial.create(32);
        MutablePolynomial b = MutablePolynomial.create(24);
        MutablePolynomial gcd = MutablePolynomial.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, SubresultantEuclid(a, b).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, true).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, false).gcd());
    }

    private static void assertGCD(MutablePolynomial a, MutablePolynomial b, MutablePolynomial gcd) {
        MutablePolynomial[] qr;
        for (MutablePolynomial dividend : Arrays.asList(a, b)) {
            qr = DivideAndRemainder.pseudoDivideAndRemainderAdaptive(dividend, gcd, true);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutablePolynomial a, MutablePolynomial b, PolynomialRemainders prs) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        a = a.clone().divide(a.content());
        b = b.clone().divide(b.content());
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        MutablePolynomial gcd = prs.gcd().clone().primitivePart();
        assertTrue(DivideAndRemainder.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivideAndRemainder.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutablePolynomial a, MutablePolynomial b, PolynomialRemainders prs, long modulus) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        assertEquals(a.clone().modulus(modulus), prs.remainders.get(0));
        assertEquals(b.clone().modulus(modulus), prs.remainders.get(1));

        MutablePolynomial gcd = prs.gcd().clone();
        assertTrue(DivideAndRemainder.divideAndRemainder(a, gcd, modulus, true)[1].isZero());
        assertTrue(DivideAndRemainder.divideAndRemainder(b, gcd, modulus, true)[1].isZero());
    }


    private static List<PolynomialRemainders> runAlgorithms(MutablePolynomial dividend, MutablePolynomial divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        PolynomialRemainders gcd(MutablePolynomial dividend, MutablePolynomial divider) {
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
        MutablePolynomial a = MutablePolynomial.create(8, 2, -1, -2, -7);
        MutablePolynomial b = MutablePolynomial.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }
}