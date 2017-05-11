package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.LongArithmetics;
import cc.r2.core.poly.ModularDomain;
import cc.r2.core.poly.univar2.PolynomialGCD.*;
import cc.r2.core.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly.IntegersDomain.IntegersDomain;
import static cc.r2.core.poly.univar2.DivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly.univar2.PolynomialGCD.*;
import static cc.r2.core.poly.univar2.PolynomialGCD.PolynomialGCD;
import static cc.r2.core.poly.univar2.RandomPolynomials.randomPoly;
import static cc.r2.core.poly.univar2.gMutablePolynomial.asLongPolyZ;
import static cc.r2.core.poly.univar2.gMutablePolynomial.asLongPolyZp;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 15/02/2017.
 */
public class PolynomialGCDTest extends AbstractPolynomialTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        long modulus = 11;
        lMutablePolynomialZp a = lMutablePolynomialZ.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950).modulus(modulus);
        lMutablePolynomialZp b = lMutablePolynomialZ.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216).modulus(modulus);

        PolynomialRemainders<lMutablePolynomialZp> prs = PolynomialGCD.Euclid(a, b);
        lMutablePolynomialZp gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test
    public void test2() throws Exception {
        //test long overflow
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(0, 14, 50, 93, 108, 130, 70);
        lMutablePolynomialZ divider = lMutablePolynomialZ.create(63, 92, 143, 245, 146, 120, 90);
        lMutablePolynomialZ gcd = PolynomialEuclid(dividend, divider, true).gcd();
        lMutablePolynomialZ expected = lMutablePolynomialZ.create(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(7, -7, 0, 1);
        lMutablePolynomialZ divider = lMutablePolynomialZ.create(-7, 0, 3);
        PolynomialRemainders<lMutablePolynomialZ> naive = PolynomialEuclid(dividend.clone(), divider.clone(), false);
        List<lMutablePolynomialZ> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(lMutablePolynomialZ.create(63, -42));
        expectedNaive.add(lMutablePolynomialZ.create(-1L));
        assertEquals(expectedNaive, naive.remainders);

        PolynomialRemainders<lMutablePolynomialZ> primitive = PolynomialEuclid(dividend.clone(), divider.clone(), true);
        List<lMutablePolynomialZ> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(lMutablePolynomialZ.create(-3, 2));
        expectedPrimitive.add(lMutablePolynomialZ.create(-1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        PolynomialRemainders<lMutablePolynomialZ> subresultant = SubresultantEuclid(dividend.clone(), divider.clone());
        List<lMutablePolynomialZ> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(lMutablePolynomialZ.create(63, -42));
        expectedSubresultant.add(lMutablePolynomialZ.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }

    @Test
    public void test4() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            lMutablePolynomialZ dividend = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(sign);
            lMutablePolynomialZ divider = lMutablePolynomialZ.create(6, 7, 9, 8, 3).multiply(sign);
            PolynomialRemainders<lMutablePolynomialZ> subresultant = SubresultantEuclid(dividend, divider);
            List<lMutablePolynomialZ> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(lMutablePolynomialZ.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(lMutablePolynomialZ.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(lMutablePolynomialZ.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(lMutablePolynomialZ.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test4a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                lMutablePolynomialZ dividend = lMutablePolynomialZ.create(1, 1, 6, -3, 5).multiply(sign1);
                lMutablePolynomialZ divider = lMutablePolynomialZ.create(-9, -3, -10, 2, 8).multiply(sign2);
                PolynomialRemainders<lMutablePolynomialZ> subresultant = SubresultantEuclid(dividend, divider);
                List<lMutablePolynomialZ> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(lMutablePolynomialZ.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(lMutablePolynomialZ.create(4344, 3818, 9774));
                expectedSubresultant.add(lMutablePolynomialZ.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(lMutablePolynomialZ.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test5() throws Exception {
        //test long overflow
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2);
        PolynomialGCD.PolynomialRemainders prs = PolynomialEuclid(dividend.clone(), dividend.clone(), false);

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        lMutablePolynomialZ divider = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2);
        PolynomialRemainders<lMutablePolynomialZ> prs;

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
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        lMutablePolynomialZ divider = lMutablePolynomialZ.create(7, 4, 10, 10, 0, 2);

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
        lMutablePolynomialZ dividend = lMutablePolynomialZ.create(1, 1, 1, 1).multiply(2);
        lMutablePolynomialZ divider = lMutablePolynomialZ.create(1, 0, 2);

        for (PolynomialGCD.PolynomialRemainders<lMutablePolynomialZ> prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test9() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, 2, 1);
        lMutablePolynomialZ expected = lMutablePolynomialZ.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            lMutablePolynomialZ dividend = randomPoly(5, rnd);
            lMutablePolynomialZ divider = randomPoly(0, rnd);
            for (PolynomialGCD.PolynomialRemainders<lMutablePolynomialZ> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(10000, 10000); i++) {
            lMutablePolynomialZ dividend = randomPoly(5, 5, rnd);
            lMutablePolynomialZ divider = randomPoly(5, 5, rnd);
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
            lMutablePolynomialZ dividend = randomPoly(10, 500, rnd);
            lMutablePolynomialZ divider = randomPoly(10, 500, rnd);
            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                lMutablePolynomialZp a = dividend.modulus(prime);
                lMutablePolynomialZp b = divider.modulus(prime);
                PolynomialRemainders<lMutablePolynomialZp> euclid = Euclid(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }

    @Test
    public void test10() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(lMutablePolynomialZ.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test11() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(0, 1, 1, -6, 17, -18, 14);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(lMutablePolynomialZ.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test12() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(1, 2, 3, 3, 6, 9);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, 3, 6, 5, 3);
        assertEquals(lMutablePolynomialZ.create(1, 2, 3), ModularGCD(a, b));

        a = lMutablePolynomialZ.create(0, 0, 1, 2);
        b = lMutablePolynomialZ.create(-1, 0, 4);
        assertEquals(lMutablePolynomialZ.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test13() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(-1, 0, 4);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(-1, 0, 0, 0, 16);
        lMutablePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(lMutablePolynomialZ.create(-1, 0, 4), gcd);
    }

    @Test
    public void test14() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(-1, 0, 4);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, -5, 6);
        lMutablePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(lMutablePolynomialZ.create(-1, 2), gcd);
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
                    lMutablePolynomialZ a = randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    lMutablePolynomialZ b = randomPoly(1 + rnd.nextInt(6), 100, rnd);
                    lMutablePolynomialZ gcd = randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    lMutablePolynomialZ mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    lMutablePolynomialZ[] qr = DivisionWithRemainder.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    lMutablePolynomialZ[] qr1 = DivisionWithRemainder.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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
        lMutablePolynomialZ poly = lMutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test16() throws Exception {
        lMutablePolynomialZ poly = lMutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
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
        lMutablePolynomialZ poly = lMutablePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }

    @Test
    public void test18() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        lMutablePolynomialZ b = lMutablePolynomialZ.zero();
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
        lMutablePolynomialZ a = lMutablePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(2);
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
        lMutablePolynomialZ a = lMutablePolynomialZ.create(32);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(24);
        lMutablePolynomialZ gcd = lMutablePolynomialZ.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, SubresultantEuclid(a, b).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, true).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, false).gcd());
    }

    @SuppressWarnings("unchecked")
    private static void assertGCD(IMutablePolynomial a, IMutablePolynomial b, IMutablePolynomial gcd) {
        IMutablePolynomial[] qr;
        for (IMutablePolynomial dividend : Arrays.asList(a, b)) {
            qr = DivisionWithRemainder.pseudoDivideAndRemainder(dividend, gcd, true);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    @SuppressWarnings("ConstantConditions")
    private static <T extends IMutablePolynomial<T>> void assertPolynomialRemainders(T a, T b, PolynomialRemainders<T> prs) {
        if (a.degree() < b.degree()) {
            assertPolynomialRemainders(b, a, prs);
            return;
        }

        if (!a.isOverField()) {
            a = a.clone().primitivePartSameSign();
            b = b.clone().primitivePartSameSign();
        }
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        T gcd = prs.gcd().clone();
        if (!a.isOverField())
            gcd = gcd.primitivePart();
        assertTrue(DivisionWithRemainder.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(DivisionWithRemainder.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }
//
//    @SuppressWarnings("ConstantConditions")
//    private static <T extends IMutablePolynomialZp<T>> void assertPolynomialRemainders(T a, T b, PolynomialRemainders<T> prs) {
//        if (a.degree() < b.degree()) {
//            assertPolynomialRemainders(b, a, prs);
//            return;
//        }
//
//        assertEquals(a, prs.remainders.get(0));
//        assertEquals(b, prs.remainders.get(1));
//
//        T gcd = prs.gcd().clone();
//        assertTrue(DivisionWithRemainder.divideAndRemainder(a, gcd, true)[1].isZero());
//        assertTrue(DivisionWithRemainder.divideAndRemainder(b, gcd, true)[1].isZero());
//    }

    private static <T extends IMutablePolynomial<T>> List<PolynomialRemainders<T>> runAlgorithms(T dividend, T divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders<T>> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        <T extends IMutablePolynomial<T>> PolynomialRemainders<T> gcd(T dividend, T divider) {
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
        lMutablePolynomialZ a = lMutablePolynomialZ.create(8, 2, -1, -2, -7);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }

    @Test
    public void test22() throws Exception {
        lMutablePolynomialZp a = lMutablePolynomialZ.create(1, 2, 3, 4, 3, 2, 1).modulus(25);
        lMutablePolynomialZp b = lMutablePolynomialZ.create(1, 2, 3, 1, 3, 2, 1).modulus(25);
        assertExtendedGCD(a, b);
    }

    @Test
    public void test23() throws Exception {
        //test long overflow
        gMutablePolynomial<BigInteger> a = gMutablePolynomial.create(0, 14, 50, 11233219232222L, 108, 130, 70);
        gMutablePolynomial<BigInteger> b = gMutablePolynomial.create(63, 92, 143, 1245222, 146, 120, 90);
        gMutablePolynomial<BigInteger> gcd1 = gMutablePolynomial.create(1, 2, 3, 4, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -4, -3, -2, -1, 999);
        gMutablePolynomial<BigInteger> gcd2 = gMutablePolynomial.create(999999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, -123124342345L);
        gMutablePolynomial<BigInteger> gcd3 = gMutablePolynomial.create(991999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, 123124342345L);
        gMutablePolynomial<BigInteger> gcd4 = gMutablePolynomial.create(Long.MAX_VALUE, Long.MAX_VALUE - 1, Long.MAX_VALUE - 2, Long.MAX_VALUE - 3, Long.MAX_VALUE - 4, Long.MAX_VALUE - 5);
        gMutablePolynomial<BigInteger> gcd = gcd1.multiply(gcd2).multiply(gcd3).multiply(gcd4);

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        gMutablePolynomial<BigInteger> gcdActual = PolynomialEuclid(a, b, false).gcd();
        assertGCD(a, b, gcdActual);

        PolynomialRemainders<gMutablePolynomial<BigInteger>> prs = SubresultantEuclid(a, b);
        assertPolynomialRemainders(a, b, prs);

        gMutablePolynomial<BigInteger> gcdSubresultant = prs.gcd();
        assertEquals(gcdActual.degree, gcdSubresultant.degree);

        gMutablePolynomial<BigInteger> gcdModular = ModularGCD(a, b);
        assertEquals(gcdActual.degree, gcdModular.degree);

        System.out.println(gMutablePolynomial.normMax(gcdActual).bitLength());
    }

    @Test
    public void testRandom1_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            gMutablePolynomial<BigInteger> dividend = randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
            gMutablePolynomial<BigInteger> divider = randomPoly(0, BigInteger.LONG_MAX_VALUE, rnd);
            for (PolynomialGCD.PolynomialRemainders<gMutablePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(20, 100); i++) {
            gMutablePolynomial<BigInteger> dividend = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            gMutablePolynomial<BigInteger> divider = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            for (PolynomialGCD.PolynomialRemainders<gMutablePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider, GCDAlgorithm.SubresultantEuclid, GCDAlgorithm.PolynomialPrimitiveEuclid)) {
                assertPolynomialRemainders(dividend, divider, prs);
            }
        }
    }

    @Test
    public void testRandom3_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(50, 500); i++) {
            gMutablePolynomial<BigInteger> dividend = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            gMutablePolynomial<BigInteger> divider = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            for (BigInteger prime : getProbablePrimesArray(BigInteger.LONG_MAX_VALUE.shiftLeft(10), 10)) {
                if (dividend.lc().mod(prime).isZero() || divider.lc().mod(prime).isZero())
                    continue;
                ModularDomain domain = new ModularDomain(prime);
                gMutablePolynomial<BigInteger> a = dividend.setDomain(domain);
                gMutablePolynomial<BigInteger> b = divider.setDomain(domain);
                PolynomialRemainders<gMutablePolynomial<BigInteger>> euclid = Euclid(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }

    @Test
    @Benchmark(runAnyway = true)
    public void testRandom5_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        int overflow = 0;
        int larger = 0;
        DescriptiveStatistics timings = null;
        for (int kk = 0; kk < 2; kk++) {
            timings = new DescriptiveStatistics();
            for (int i = 0; i < its(300, 1000); i++) {
                try {
                    gMutablePolynomial<BigInteger> a = randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    gMutablePolynomial<BigInteger> b = randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    gMutablePolynomial<BigInteger> gcd = randomPoly(2 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    gMutablePolynomial<BigInteger> mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    gMutablePolynomial<BigInteger>[] qr = DivisionWithRemainder.pseudoDivideAndRemainder(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    gMutablePolynomial<BigInteger>[] qr1 = DivisionWithRemainder.pseudoDivideAndRemainder(mgcd.clone(), gcd, false);
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
    public void test24() throws Exception {
        long modulus = 419566991;

        gMutablePolynomial<BigInteger> poly = gMutablePolynomial.create(IntegersDomain, new BigInteger("-4914"), new BigInteger("6213"), new BigInteger("3791"), new BigInteger("996"), new BigInteger("-13304"), new BigInteger("-1567"), new BigInteger("2627"), new BigInteger("15845"), new BigInteger("-12626"), new BigInteger("-6383"), new BigInteger("294"), new BigInteger("26501"), new BigInteger("-17063"), new BigInteger("-14635"), new BigInteger("9387"), new BigInteger("-7141"), new BigInteger("-8185"), new BigInteger("17856"), new BigInteger("4431"), new BigInteger("-13075"), new BigInteger("-7050"), new BigInteger("14672"), new BigInteger("3690"), new BigInteger("-3990"));
        gMutablePolynomial<BigInteger> a = gMutablePolynomial.create(IntegersDomain, new BigInteger("419563715"), new BigInteger("419566193"), new BigInteger("3612"), new BigInteger("3444"), new BigInteger("419563127"), new BigInteger("419564681"), new BigInteger("419565017"), new BigInteger("419564387"), new BigInteger("419563463"), new BigInteger("3192"), new BigInteger("419563841"), new BigInteger("419563001"));
        gMutablePolynomial<BigInteger> b = gMutablePolynomial.create(IntegersDomain, new BigInteger("209783497"), new BigInteger("9989688"), new BigInteger("379608231"), new BigInteger("399587609"), new BigInteger("59938143"), new BigInteger("29969072"), new BigInteger("99896901"), new BigInteger("359628849"), new BigInteger("329659781"), new BigInteger("239752567"), new BigInteger("19979379"), new BigInteger("179814423"), new BigInteger("1"));

        ModularDomain domain = new ModularDomain(modulus);
        gMutablePolynomial<BigInteger> aMod = a.setDomain(domain).monic(poly.lc());
        gMutablePolynomial<BigInteger> bMod = b.setDomain(domain).monic();
        gMutablePolynomial<BigInteger>[] xgcd = PolynomialGCD.ExtendedEuclid(aMod, bMod);
        lMutablePolynomialZp[] lxgcd = PolynomialGCD.ExtendedEuclid(asLongPolyZp(aMod), asLongPolyZp(bMod));
        assertEquals(asLongPolyZp(xgcd[0]), lxgcd[0]);
        assertEquals(asLongPolyZp(xgcd[1]), lxgcd[1]);
        assertEquals(asLongPolyZp(xgcd[2]), lxgcd[2]);
    }

    static <T extends lMutablePolynomialAbstract<T>> void assertExtendedGCD(T a, T b) {
        assertExtendedGCD(ExtendedEuclid(a, b), a, b);
    }

    static <T extends lMutablePolynomialAbstract<T>> void assertExtendedGCD(T[] eea, T a, T b) {
        assertEquals(eea[0], a.clone().multiply(eea[1]).add(b.clone().multiply(eea[2])));
        assertEquals(eea[0].degree, PolynomialGCD(a, b).degree);
    }
}