package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly2.PolynomialGCD.*;
import cc.r2.core.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly2.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly2.PolynomialGCD.*;
import static cc.r2.core.poly2.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.poly2.RandomPolynomials.randomPoly;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 15/02/2017.
 */
public class PolynomialGCDTest extends AbstractPolynomialTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        long modulus = 11;
        MutablePolynomialMod a = MutablePolynomialZ.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950).modulus(modulus);
        MutablePolynomialMod b = MutablePolynomialZ.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216).modulus(modulus);

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

        for (PolynomialGCD.PolynomialRemainders<MutablePolynomialZ> prs : runAlgorithms(dividend, divider))
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
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            MutablePolynomialZ dividend = randomPoly(5, rnd);
            MutablePolynomialZ divider = randomPoly(0, rnd);
            for (PolynomialGCD.PolynomialRemainders<MutablePolynomialZ> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(10000, 10000); i++) {
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
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(500, 3000); i++) {
            MutablePolynomialZ dividend = randomPoly(10, 500, rnd);
            MutablePolynomialZ divider = randomPoly(10, 500, rnd);
            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                MutablePolynomialMod a = dividend.modulus(prime);
                MutablePolynomialMod b = divider.modulus(prime);
                PolynomialRemainders<MutablePolynomialMod> euclid = Euclid(a, b);
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
    @Benchmark(runAnyway = true)
    public void testRandom5() throws Exception {
        RandomGenerator rnd = getRandom();
        int overflow = 0;
        int larger = 0;
        DescriptiveStatistics timings = null;
        for (int k = 0; k < 2; k++) {
            timings = new DescriptiveStatistics();
            for (int i = 0; i < its(5000, 10000); i++) {
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

                    MutablePolynomialZ[] qr = DivisionWithRemainder.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    MutablePolynomialZ[] qr1 = DivisionWithRemainder.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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

    @SuppressWarnings("unchecked")
    private static void assertGCD(IMutablePolynomialZ a, IMutablePolynomialZ b, IMutablePolynomialZ gcd) {
        IMutablePolynomial[] qr;
        for (IMutablePolynomialZ dividend : Arrays.asList(a, b)) {
            qr = DivisionWithRemainder.pseudoDivideAndRemainder(dividend, gcd, true);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    @SuppressWarnings("ConstantConditions")
    private static <T extends IMutablePolynomialZ<T>> void assertPolynomialRemainders(T a, T b, PolynomialRemainders<T> prs) {
        if (a.degree() < b.degree()) {
            assertPolynomialRemainders(b, a, prs);
            return;
        }

        a = a.clone().primitivePartSameSign();
        b = b.clone().primitivePartSameSign();
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        T gcd = prs.gcd().clone().primitivePart();
        assertTrue(DivisionWithRemainder.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivisionWithRemainder.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }

    @SuppressWarnings("ConstantConditions")
    private static <T extends IMutablePolynomialZp<T>> void assertPolynomialRemainders(T a, T b, PolynomialRemainders<T> prs) {
        if (a.degree() < b.degree()) {
            assertPolynomialRemainders(b, a, prs);
            return;
        }

        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        T gcd = prs.gcd().clone();
        assertTrue(DivisionWithRemainder.divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivisionWithRemainder.divideAndRemainder(b, gcd, true)[1].isZero());
    }

    private static <T extends IMutablePolynomialZ<T>> List<PolynomialRemainders<T>> runAlgorithms(T dividend, T divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders<T>> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        <T extends IMutablePolynomialZ<T>> PolynomialRemainders<T> gcd(T dividend, T divider) {
            switch (this) {
                case PolynomialEuclid:
                    return PolynomialEuclid(dividend, divider, false);
                case PolynomialPrimitiveEuclid:
                    return PolynomialEuclid(dividend, divider, true);
                case SubresultantEuclid:
                    return SubresultantEuclid(dividend, divider);
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

    @Test
    public void test22() throws Exception {
        MutablePolynomialMod a = MutablePolynomialZ.create(1, 2, 3, 4, 3, 2, 1).modulus(25);
        MutablePolynomialMod b = MutablePolynomialZ.create(1, 2, 3, 1, 3, 2, 1).modulus(25);
        assertExtendedGCD(a, b);
    }

    @Test
    public void test23() throws Exception {
        //test long overflow
        bMutablePolynomialZ a = bMutablePolynomialZ.create(0, 14, 50, 11233219232222L, 108, 130, 70);
        bMutablePolynomialZ b = bMutablePolynomialZ.create(63, 92, 143, 1245222, 146, 120, 90);
        bMutablePolynomialZ gcd1 = bMutablePolynomialZ.create(1, 2, 3, 4, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -4, -3, -2, -1, 999);
        bMutablePolynomialZ gcd2 = bMutablePolynomialZ.create(999999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, -123124342345L);
        bMutablePolynomialZ gcd3 = bMutablePolynomialZ.create(991999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, 123124342345L);
        bMutablePolynomialZ gcd4 = bMutablePolynomialZ.create(Long.MAX_VALUE, Long.MAX_VALUE - 1, Long.MAX_VALUE - 2, Long.MAX_VALUE - 3, Long.MAX_VALUE - 4, Long.MAX_VALUE - 5);
        bMutablePolynomialZ gcd = gcd1.multiply(gcd2).multiply(gcd3).multiply(gcd4);

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        bMutablePolynomialZ gcdActual = PolynomialEuclid(a, b, false).gcd();
        assertGCD(a, b, gcdActual);

        PolynomialRemainders<bMutablePolynomialZ> prs = SubresultantEuclid(a, b);
        assertPolynomialRemainders(a, b, prs);

        bMutablePolynomialZ gcdSubresultant = prs.gcd();
        assertEquals(gcdActual.degree, gcdSubresultant.degree);

        System.out.println(gcdActual.normMax().bitLength());
    }

    @Test
    public void testRandom1_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            bMutablePolynomialZ dividend = randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
            bMutablePolynomialZ divider = randomPoly(0, BigInteger.LONG_MAX_VALUE, rnd);
            for (PolynomialGCD.PolynomialRemainders<bMutablePolynomialZ> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(20, 100); i++) {
            bMutablePolynomialZ dividend = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            bMutablePolynomialZ divider = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            for (PolynomialGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider, GCDAlgorithm.SubresultantEuclid, GCDAlgorithm.PolynomialPrimitiveEuclid)) {
                assertPolynomialRemainders(dividend, divider, prs);
            }
        }
    }


    @Test
    public void testRandom3_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(50, 500); i++) {
            System.out.println(i);
            bMutablePolynomialZ dividend = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            bMutablePolynomialZ divider = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            for (BigInteger prime : getProbablePrimesArray(BigInteger.LONG_MAX_VALUE.shiftLeft(10), 10)) {
                if (dividend.lc().mod(prime).isZero() || divider.lc().mod(prime).isZero())
                    continue;
                bMutablePolynomialMod a = dividend.modulus(prime);
                bMutablePolynomialMod b = divider.modulus(prime);
                PolynomialRemainders<bMutablePolynomialMod> euclid = Euclid(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }


    static <T extends MutablePolynomialAbstract<T>> void assertExtendedGCD(T a, T b) {
        assertExtendedGCD(ExtendedEuclid(a, b), a, b);
    }

    static <T extends MutablePolynomialAbstract<T>> void assertExtendedGCD(T[] eea, T a, T b) {
        Assert.assertEquals(eea[0], a.clone().multiply(eea[1]).add(b.clone().multiply(eea[2])));
        assertEquals(eea[0].degree, PolynomialGCD(a, b).degree);
    }

}