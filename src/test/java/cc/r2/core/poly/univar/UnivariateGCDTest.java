package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.Rational;
import cc.r2.core.poly.test.APolynomialTest;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.MachineArithmetic;
import cc.r2.core.poly.univar.UnivariateGCD.*;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.poly.Integers.Integers;
import static cc.r2.core.poly.Rationals.Rationals;
import static cc.r2.core.poly.univar.UnivariateDivision.divideAndRemainder;
import static cc.r2.core.poly.univar.RandomUnivariatePolynomials.randomPoly;
import static cc.r2.core.poly.univar.UnivariateGCD.*;
import static cc.r2.core.poly.univar.UnivariatePolynomial.asLongPolyZp;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 15/02/2017.
 */
public class UnivariateGCDTest extends APolynomialTest {
    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        long modulus = 11;
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950).modulus(modulus);
        lUnivariatePolynomialZp b = lUnivariatePolynomialZ.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216).modulus(modulus);

        PolynomialRemainders<lUnivariatePolynomialZp> prs = UnivariateGCD.EuclidRemainders(a, b);
        lUnivariatePolynomialZp gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test
    public void test2() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(0, 14, 50, 93, 108, 130, 70);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(63, 92, 143, 245, 146, 120, 90);
        lUnivariatePolynomialZ gcd = EuclidPseudoRemainders(dividend, divider, true).gcd();
        lUnivariatePolynomialZ expected = lUnivariatePolynomialZ.create(7, 4, 10, 10);
        assertEquals(expected, gcd);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(7, -7, 0, 1);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(-7, 0, 3);
        PolynomialRemainders<lUnivariatePolynomialZ> naive = EuclidPseudoRemainders(dividend.clone(), divider.clone(), false);
        List<lUnivariatePolynomialZ> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(lUnivariatePolynomialZ.create(63, -42));
        expectedNaive.add(lUnivariatePolynomialZ.create(1L));
        assertEquals(expectedNaive, naive.remainders);

        PolynomialRemainders<lUnivariatePolynomialZ> primitive = EuclidPseudoRemainders(dividend.clone(), divider.clone(), true);
        List<lUnivariatePolynomialZ> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(lUnivariatePolynomialZ.create(-3, 2));
        expectedPrimitive.add(lUnivariatePolynomialZ.create(1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        PolynomialRemainders<lUnivariatePolynomialZ> subresultant = EuclidSubresultantRemainders(dividend.clone(), divider.clone());
        List<lUnivariatePolynomialZ> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(lUnivariatePolynomialZ.create(63, -42));
        expectedSubresultant.add(lUnivariatePolynomialZ.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }

    @Test
    public void test4() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(sign);
            lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(6, 7, 9, 8, 3).multiply(sign);
            PolynomialRemainders<lUnivariatePolynomialZ> subresultant = EuclidSubresultantRemainders(dividend, divider);
            List<lUnivariatePolynomialZ> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(lUnivariatePolynomialZ.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(lUnivariatePolynomialZ.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(lUnivariatePolynomialZ.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(lUnivariatePolynomialZ.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test4a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(1, 1, 6, -3, 5).multiply(sign1);
                lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(-9, -3, -10, 2, 8).multiply(sign2);
                PolynomialRemainders<lUnivariatePolynomialZ> subresultant = EuclidSubresultantRemainders(dividend, divider);
                List<lUnivariatePolynomialZ> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(lUnivariatePolynomialZ.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(lUnivariatePolynomialZ.create(4344, 3818, 9774));
                expectedSubresultant.add(lUnivariatePolynomialZ.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(lUnivariatePolynomialZ.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test5() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2);
        UnivariateGCD.PolynomialRemainders prs = EuclidPseudoRemainders(dividend.clone(), dividend.clone(), false);

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2);
        PolynomialRemainders<lUnivariatePolynomialZ> prs;

        for (boolean prim : new boolean[]{true, false}) {
            prs = EuclidPseudoRemainders(dividend.clone(), divider.clone(), prim);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());

            prs = EuclidPseudoRemainders(divider.clone(), dividend.clone(), prim);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        prs = EuclidSubresultantRemainders(dividend.clone(), divider.clone());
        assertEquals(2, prs.remainders.size());
        assertEquals(divider, prs.gcd());
    }

    @Test
    public void test7() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2).multiply(2);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(7, 4, 10, 10, 0, 2);

        for (UnivariateGCD.PolynomialRemainders prs : runAlgorithms(dividend, divider)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        for (UnivariateGCD.PolynomialRemainders prs : runAlgorithms(divider, dividend)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }
    }

    @Test
    public void test8() throws Exception {
        //test long overflow
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(1, 1, 1, 1).multiply(2);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(1, 0, 2);

        for (UnivariateGCD.PolynomialRemainders<lUnivariatePolynomialZ> prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test9() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, 2, 1);
        lUnivariatePolynomialZ expected = lUnivariatePolynomialZ.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            lUnivariatePolynomialZ dividend = randomPoly(5, rnd);
            lUnivariatePolynomialZ divider = randomPoly(0, rnd);
            for (UnivariateGCD.PolynomialRemainders<lUnivariatePolynomialZ> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2() throws Exception {
        RandomGenerator rnd = getRandom();
        int overflows = 0;
        int nIterations = its(10000, 10000);
        for (int i = 0; i < nIterations; i++) {
            lUnivariatePolynomialZ dividend = randomPoly(5, 5, rnd);
            lUnivariatePolynomialZ divider = randomPoly(5, 5, rnd);
            try {
                for (UnivariateGCD.PolynomialRemainders prs :
                        runAlgorithms(dividend, divider,
                                GCDAlgorithm.PolynomialPrimitiveEuclid,
                                GCDAlgorithm.SubresultantEuclid)) {
                    assertPolynomialRemainders(dividend, divider, prs);
                }
            } catch (ArithmeticException e) {
                if (!e.getMessage().equals("long overflow"))
                    throw e;
                ++overflows;
            }
        }
        assertTrue(overflows <= Math.ceil(nIterations / 1000.));
        if (overflows > 0)
            System.out.println("Overflows: " + overflows);
    }

    @Test
    public void testRandom2a() throws Exception {
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(4, -1, -4, -3, 2, 4);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(-4, 1, 4, 3, -2, -4);
        assertPolynomialRemainders(dividend, divider, UnivariateGCD.EuclidPseudoRemainders(dividend, divider, true));
    }

    @Test
    public void testRandom2b() throws Exception {
        lUnivariatePolynomialZ dividend = lUnivariatePolynomialZ.create(0, 0, 4, 0, 4, 4);
        lUnivariatePolynomialZ divider = lUnivariatePolynomialZ.create(0, 0, 4, 0, 4, 4);
        assertPolynomialRemainders(dividend, divider, UnivariateGCD.EuclidPseudoRemainders(dividend, divider, true));
    }

    @Test
    public void testRandom3() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(500, 3000); i++) {
            lUnivariatePolynomialZ dividend = randomPoly(10, 500, rnd);
            lUnivariatePolynomialZ divider = randomPoly(10, 500, rnd);
            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                lUnivariatePolynomialZp a = dividend.modulus(prime);
                lUnivariatePolynomialZp b = divider.modulus(prime);
                PolynomialRemainders<lUnivariatePolynomialZp> euclid = EuclidRemainders(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }

    @Test
    public void test10() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(lUnivariatePolynomialZ.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test11() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(0, 1, 1, -6, 17, -18, 14);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(lUnivariatePolynomialZ.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test12() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(1, 2, 3, 3, 6, 9);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, 3, 6, 5, 3);
        assertEquals(lUnivariatePolynomialZ.create(1, 2, 3), ModularGCD(a, b));

        a = lUnivariatePolynomialZ.create(0, 0, 1, 2);
        b = lUnivariatePolynomialZ.create(-1, 0, 4);
        assertEquals(lUnivariatePolynomialZ.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test13() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(-1, 0, 4);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(-1, 0, 0, 0, 16);
        lUnivariatePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(lUnivariatePolynomialZ.create(-1, 0, 4), gcd);
    }

    @Test
    public void test14() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(-1, 0, 4);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, -5, 6);
        lUnivariatePolynomialZ gcd = ModularGCD(a, b);
        assertEquals(lUnivariatePolynomialZ.create(-1, 2), gcd);
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
                    lUnivariatePolynomialZ a = randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    lUnivariatePolynomialZ b = randomPoly(1 + rnd.nextInt(6), 100, rnd);
                    lUnivariatePolynomialZ gcd = randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    lUnivariatePolynomialZ mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    lUnivariatePolynomialZ[] qr = UnivariateDivision.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    lUnivariatePolynomialZ[] qr1 = UnivariateDivision.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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
        lUnivariatePolynomialZ poly = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test16() throws Exception {
        lUnivariatePolynomialZ poly = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.multiply(MachineArithmetic.safePow(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divideOrNull(MachineArithmetic.safePow(5, poly.degree));
        poly.multiply(512);
        assertEquals(-654063683625L, poly.evaluateAtRational(17, 2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void test17() throws Exception {
        lUnivariatePolynomialZ poly = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }

    @Test
    public void test18() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.zero();
        assertEquals(a, EuclidSubresultantRemainders(a, b).gcd());
        assertEquals(a, EuclidSubresultantRemainders(b, a).gcd());
        assertEquals(a, EuclidPseudoRemainders(a, b, true).gcd());
        assertEquals(a, EuclidPseudoRemainders(b, a, true).gcd());
        assertEquals(a, EuclidGCD(a, b));
        assertEquals(a, EuclidGCD(b, a));
        assertEquals(a, PolynomialGCD(a, b));
        assertEquals(a, PolynomialGCD(b, a));
    }

    @Test
    public void test19() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(8, -2 * 8, 8, 8 * 2);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(2);
        assertEquals(b, EuclidSubresultantRemainders(a, b).gcd());
        assertEquals(b, EuclidSubresultantRemainders(b, a).gcd());
        assertEquals(b, EuclidPseudoRemainders(a, b, true).gcd());
        assertEquals(b, EuclidPseudoRemainders(b, a, true).gcd());
        assertEquals(b, EuclidPseudoRemainders(a, b, false).gcd());
        assertEquals(b, EuclidPseudoRemainders(b, a, false).gcd());
        assertEquals(b, EuclidGCD(a, b));
        assertEquals(b, EuclidGCD(b, a));
        assertEquals(b, PolynomialGCD(a, b));
        assertEquals(b, PolynomialGCD(b, a));
    }

    @Test
    public void test20() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(32);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(24);
        lUnivariatePolynomialZ gcd = lUnivariatePolynomialZ.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, EuclidSubresultantRemainders(a, b).gcd());
        assertEquals(gcd, EuclidPseudoRemainders(a, b, true).gcd());
        assertEquals(gcd, EuclidPseudoRemainders(a, b, false).gcd());
    }

    @SuppressWarnings("unchecked")
    private static void assertGCD(IUnivariatePolynomial a, IUnivariatePolynomial b, IUnivariatePolynomial gcd) {
        IUnivariatePolynomial[] qr;
        for (IUnivariatePolynomial dividend : Arrays.asList(a, b)) {
            qr = UnivariateDivision.pseudoDivideAndRemainder(dividend, gcd, true);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    @SuppressWarnings("ConstantConditions")
    private static <T extends IUnivariatePolynomial<T>> void assertPolynomialRemainders(T a, T b, PolynomialRemainders<T> prs) {
        if (a.degree() < b.degree()) {
            assertPolynomialRemainders(b, a, prs);
            return;
        }

        if (!a.isOverField()) {
            a = a.clone().primitivePartSameSign();
            b = b.clone().primitivePartSameSign();
        }
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.size() == 2 ? prs.gcd().primitivePartSameSign() : prs.remainders.get(1));

        T gcd = prs.gcd().clone();
        if (!a.isOverField())
            gcd = gcd.primitivePart();
        assertTrue(UnivariateDivision.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(UnivariateDivision.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
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
//        assertTrue(UnivariateDivision.divideAndRemainder(a, gcd, true)[1].isZero());
//        assertTrue(UnivariateDivision.divideAndRemainder(b, gcd, true)[1].isZero());
//    }

    private static <T extends IUnivariatePolynomial<T>> List<PolynomialRemainders<T>> runAlgorithms(T dividend, T divider, GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainders<T>> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        <T extends IUnivariatePolynomial<T>> PolynomialRemainders<T> gcd(T dividend, T divider) {
            switch (this) {
                case PolynomialEuclid:
                    return EuclidPseudoRemainders(dividend, divider, false);
                case PolynomialPrimitiveEuclid:
                    return EuclidPseudoRemainders(dividend, divider, true);
                case SubresultantEuclid:
                    return EuclidSubresultantRemainders(dividend, divider);
            }
            throw new IllegalArgumentException();
        }
    }

    @Test
    public void test21() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(8, 2, -1, -2, -7);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }

    @Test
    public void test22() throws Exception {
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(1, 2, 3, 4, 3, 2, 1).modulus(25);
        lUnivariatePolynomialZp b = lUnivariatePolynomialZ.create(1, 2, 3, 1, 3, 2, 1).modulus(25);
        assertExtendedGCD(a, b);
    }

    @Test
    public void test23() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(0, 14, 50, 11233219232222L, 108, 130, 70);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(63, 92, 143, 1245222, 146, 120, 90);
        UnivariatePolynomial<BigInteger> gcd1 = UnivariatePolynomial.create(1, 2, 3, 4, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -4, -3, -2, -1, 999);
        UnivariatePolynomial<BigInteger> gcd2 = UnivariatePolynomial.create(999999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, -123124342345L);
        UnivariatePolynomial<BigInteger> gcd3 = UnivariatePolynomial.create(991999L, 123L, 123L, 342425L, 312L, -12312432423L, 13212123123123L, 123124342345L);
        UnivariatePolynomial<BigInteger> gcd4 = UnivariatePolynomial.create(Long.MAX_VALUE, Long.MAX_VALUE - 1, Long.MAX_VALUE - 2, Long.MAX_VALUE - 3, Long.MAX_VALUE - 4, Long.MAX_VALUE - 5);
        UnivariatePolynomial<BigInteger> gcd = gcd1.multiply(gcd2).multiply(gcd3).multiply(gcd4);

        a = a.multiply(gcd);
        b = b.multiply(gcd);

        UnivariatePolynomial<BigInteger> gcdActual = EuclidPseudoRemainders(a, b, false).gcd();
        assertGCD(a, b, gcdActual);

        PolynomialRemainders<UnivariatePolynomial<BigInteger>> prs = EuclidSubresultantRemainders(a, b);
        assertPolynomialRemainders(a, b, prs);

        UnivariatePolynomial<BigInteger> gcdSubresultant = prs.gcd();
        assertEquals(gcdActual.degree, gcdSubresultant.degree);

        UnivariatePolynomial<BigInteger> gcdModular = ModularGCD(a, b);
        assertEquals(gcdActual.degree, gcdModular.degree);

        System.out.println(UnivariatePolynomial.normMax(gcdActual).bitLength());
    }

    @Test
    public void testRandom1_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomial<BigInteger> dividend = randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
            UnivariatePolynomial<BigInteger> divider = randomPoly(0, BigInteger.LONG_MAX_VALUE, rnd);
            for (UnivariateGCD.PolynomialRemainders<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider)) {
                assertEquals(0, prs.gcd().degree);
            }
        }
    }

    @Test
    public void testRandom2_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(20, 100); i++) {
            UnivariatePolynomial<BigInteger> dividend = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            UnivariatePolynomial<BigInteger> divider = randomPoly(getRandomData().nextInt(10, 50), BigInteger.LONG_MAX_VALUE.shiftRight(1), rnd);
            for (UnivariateGCD.PolynomialRemainders<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider, GCDAlgorithm.SubresultantEuclid, GCDAlgorithm.PolynomialPrimitiveEuclid)) {
                assertPolynomialRemainders(dividend, divider, prs);
            }
        }
    }

    @Test
    public void testRandom3_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(50, 500); i++) {
            UnivariatePolynomial<BigInteger> dividend = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            UnivariatePolynomial<BigInteger> divider = randomPoly(getRandomData().nextInt(10, 30), BigInteger.LONG_MAX_VALUE, rnd);
            for (BigInteger prime : getProbablePrimesArray(BigInteger.LONG_MAX_VALUE.shiftLeft(10), 10)) {
                if (dividend.lc().mod(prime).isZero() || divider.lc().mod(prime).isZero())
                    continue;
                IntegersModulo domain = new IntegersModulo(prime);
                UnivariatePolynomial<BigInteger> a = dividend.setDomain(domain);
                UnivariatePolynomial<BigInteger> b = divider.setDomain(domain);
                PolynomialRemainders<UnivariatePolynomial<BigInteger>> euclid = EuclidRemainders(a, b);
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
                    UnivariatePolynomial<BigInteger> a = randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    UnivariatePolynomial<BigInteger> b = randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    UnivariatePolynomial<BigInteger> gcd = randomPoly(2 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    UnivariatePolynomial<BigInteger> mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    UnivariatePolynomial<BigInteger>[] qr = UnivariateDivision.pseudoDivideAndRemainder(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    UnivariatePolynomial<BigInteger>[] qr1 = UnivariateDivision.pseudoDivideAndRemainder(mgcd.clone(), gcd, false);
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

        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(Integers, new BigInteger("-4914"), new BigInteger("6213"), new BigInteger("3791"), new BigInteger("996"), new BigInteger("-13304"), new BigInteger("-1567"), new BigInteger("2627"), new BigInteger("15845"), new BigInteger("-12626"), new BigInteger("-6383"), new BigInteger("294"), new BigInteger("26501"), new BigInteger("-17063"), new BigInteger("-14635"), new BigInteger("9387"), new BigInteger("-7141"), new BigInteger("-8185"), new BigInteger("17856"), new BigInteger("4431"), new BigInteger("-13075"), new BigInteger("-7050"), new BigInteger("14672"), new BigInteger("3690"), new BigInteger("-3990"));
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(Integers, new BigInteger("419563715"), new BigInteger("419566193"), new BigInteger("3612"), new BigInteger("3444"), new BigInteger("419563127"), new BigInteger("419564681"), new BigInteger("419565017"), new BigInteger("419564387"), new BigInteger("419563463"), new BigInteger("3192"), new BigInteger("419563841"), new BigInteger("419563001"));
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(Integers, new BigInteger("209783497"), new BigInteger("9989688"), new BigInteger("379608231"), new BigInteger("399587609"), new BigInteger("59938143"), new BigInteger("29969072"), new BigInteger("99896901"), new BigInteger("359628849"), new BigInteger("329659781"), new BigInteger("239752567"), new BigInteger("19979379"), new BigInteger("179814423"), new BigInteger("1"));

        IntegersModulo domain = new IntegersModulo(modulus);
        UnivariatePolynomial<BigInteger> aMod = a.setDomain(domain).monic(poly.lc());
        UnivariatePolynomial<BigInteger> bMod = b.setDomain(domain).monic();
        UnivariatePolynomial<BigInteger>[] xgcd = UnivariateGCD.ExtendedEuclidGCD(aMod, bMod);
        lUnivariatePolynomialZp[] lxgcd = UnivariateGCD.ExtendedEuclidGCD(asLongPolyZp(aMod), asLongPolyZp(bMod));
        assertEquals(asLongPolyZp(xgcd[0]), lxgcd[0]);
        assertEquals(asLongPolyZp(xgcd[1]), lxgcd[1]);
        assertEquals(asLongPolyZp(xgcd[2]), lxgcd[2]);
    }

    @Test
    public void test25() throws Exception {
        UnivariatePolynomial<Rational> a = UnivariatePolynomial.create(Rationals,
                Rationals.parse("2/3"),
                Rationals.parse("4/5"),
                Rationals.parse("1/2"),
                Rationals.parse("-31/2"));
        UnivariatePolynomial<Rational> b = UnivariatePolynomial.create(Rationals,
                Rationals.parse("7/3"),
                Rationals.parse("4/7"),
                Rationals.parse("3/2"),
                Rationals.parse("-31/12"));
        UnivariatePolynomial<Rational> gcd = UnivariatePolynomial.create(Rationals,
                Rationals.parse("4/3"),
                Rationals.parse("-4/7"),
                Rationals.parse("-1/2"),
                Rationals.parse("-1/12"));
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertGCD(a, b, PolynomialGCD(a, b));
    }

    @Test
    public void test26() throws Exception {
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(1, 2, 3, 1, 2, 3, 4, 3, 2, 1).modulus(29);//.square
        // ().square().square().square().square().square().square().square().square().square();
        lUnivariatePolynomialZp b = lUnivariatePolynomialZ.create(1, 2, 3, 1, 1, 2, 3, 3).modulus(29);
        assertArrayEquals(ExtendedEuclidGCD(a, b), ExtendedHalfGCD(a, b));
    }

    @Test
    public void test27() throws Exception {
        // assert java heap
        lUnivariatePolynomialZp a = RandomUnivariatePolynomials.randomMonicPoly(15_000, 19, getRandom());
        lUnivariatePolynomialZp b = RandomUnivariatePolynomials.randomMonicPoly(15_000, 19, getRandom());
        assertExtendedGCD(ExtendedEuclidGCD(a, b), a, b);
    }

    @Test
    public void test29_randomHalfGCD() throws Exception {
        int
                minimalDegree = UnivariateGCD.SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE,
                maximalDegree = minimalDegree * 4,
                nIterations = its(1000, 10000);
        testHalfGCDRandom(minimalDegree, maximalDegree, nIterations);
    }

    @Benchmark(runAnyway = true)
    @Test
    public void test30_randomHalfGCD() throws Exception {
        int
                minimalDegree = 5 * UnivariateGCD.SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE,
                maximalDegree = minimalDegree * 4,
                nIterations = its(100, 100);
        testHalfGCDRandom(minimalDegree, maximalDegree, nIterations);
    }

    private static void testHalfGCDRandom(int minimalDegree, int maximalDegree, int nIterations) throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                euclid = new DescriptiveStatistics(),
                half = new DescriptiveStatistics();
        int nonTrivGCD = 0;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                euclid.clear();
                half.clear();
            }

            long modulus = getModulusRandom(rndd.nextInt(2, 20));
            lUnivariatePolynomialZp a = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree), modulus, rnd).multiply(1 + rnd.nextLong());
            lUnivariatePolynomialZp b = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree), modulus, rnd).multiply(1 + rnd.nextLong());
            try {
                long start;

                start = System.nanoTime();
                lUnivariatePolynomialZp expected = EuclidGCD(a, b).monic();
                euclid.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                lUnivariatePolynomialZp actual = HalfGCD(a, b).monic();
                half.addValue(System.nanoTime() - start);

                assertEquals(expected, actual);

                if (!expected.isConstant())
                    ++nonTrivGCD;
            } catch (Throwable tr) {
                System.out.println("lUnivariatePolynomialZ." + a.toStringForCopy() + ".modulus(" + modulus + ");");
                System.out.println("lUnivariatePolynomialZ." + b.toStringForCopy() + ".modulus(" + modulus + ");");
                throw tr;
            }
        }
        System.out.println("Non-trivial gcds: " + nonTrivGCD);
        System.out.println("Euclid:  " + TimeUnits.statisticsNanotime(euclid));
        System.out.println("HalfGCD: " + TimeUnits.statisticsNanotime(half));
    }

    @Test
    public void test31_randomExtendedHalfGCD() throws Exception {
        int
                minimalDegree = UnivariateGCD.SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE,
                maximalDegree = minimalDegree * 4,
                nIterations = its(1000, 2000);
        testExtendedHalfGCDRandom(minimalDegree, maximalDegree, nIterations);
    }

    @Benchmark(runAnyway = true)
    @Test
    public void test32_randomExtendedHalfGCD() throws Exception {
        int
                minimalDegree = 5 * UnivariateGCD.SWITCH_TO_HALF_GCD_ALGORITHM_DEGREE,
                maximalDegree = minimalDegree * 4,
                nIterations = its(100, 100);
        testExtendedHalfGCDRandom(minimalDegree, maximalDegree, nIterations);
    }

    private static void testExtendedHalfGCDRandom(int minimalDegree, int maximalDegree, int nIterations) throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                euclid = new DescriptiveStatistics(),
                half = new DescriptiveStatistics();
        int nonTrivGCD = 0;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                euclid.clear();
                half.clear();
            }

            long modulus = getModulusRandom(rndd.nextInt(2, 20));
            lUnivariatePolynomialZp a = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree),
                    modulus, rnd).multiply(1 + rnd.nextLong());
            lUnivariatePolynomialZp b = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree),
                    modulus, rnd).multiply(1 + rnd.nextLong());
            try {
                long start;

                start = System.nanoTime();
                lUnivariatePolynomialZp[] expected = ExtendedEuclidGCD(a, b);
                euclid.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                lUnivariatePolynomialZp[] actual = ExtendedHalfGCD(a, b);
                half.addValue(System.nanoTime() - start);

                assertArrayEquals(expected, actual);

                if (!expected[0].isConstant())
                    ++nonTrivGCD;
            } catch (Throwable tr) {
                System.out.println("lUnivariatePolynomialZ." + a.toStringForCopy() + ".modulus(" + modulus + ");");
                System.out.println("lUnivariatePolynomialZ." + b.toStringForCopy() + ".modulus(" + modulus + ");");
                throw tr;
            }
        }
        System.out.println("Non-trivial gcds: " + nonTrivGCD);
        System.out.println("Euclid:  " + TimeUnits.statisticsNanotime(euclid));
        System.out.println("HalfGCD: " + TimeUnits.statisticsNanotime(half));
    }

    static <T extends IUnivariatePolynomial<T>> void assertExtendedGCD(T a, T b) {
        assertExtendedGCD(ExtendedEuclidGCD(a, b), a, b);
    }

    static <T extends IUnivariatePolynomial<T>> void assertExtendedGCD(T[] eea, T a, T b) {
        assertEquals(eea[0], a.clone().multiply(eea[1]).add(b.clone().multiply(eea[2])));
        assertEquals(eea[0].degree(), PolynomialGCD(a, b).degree());
    }
}