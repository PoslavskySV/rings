package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static cc.r2.core.polynomial.LongArithmetics.gcd;
import static cc.r2.core.polynomial.LongArithmetics.powExact;
import static cc.r2.core.polynomial.SmallPolynomials.*;
import static org.junit.Assert.*;

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

    @Test(expected = ArithmeticException.class)
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

    @Test
    public void testRandom3() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long[] primes = {3, 5, 7, 11, 13, 67, 97, 113, 127};
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly dividend = randomPoly(10, 500, rnd);
            MutableLongPoly divider = randomPoly(10, 500, rnd);
            for (long prime : primes) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                PolynomialRemainders euclid = Euclid(dividend, divider, prime);
                assertPolynomialRemainders(dividend, divider, euclid, prime);
            }
        }
    }

    @Test
    public void test13() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        MutableLongPoly b = MutableLongPoly.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(MutableLongPoly.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test14() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(0, 1, 1, -6, 17, -18, 14);
        MutableLongPoly b = MutableLongPoly.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(MutableLongPoly.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test15() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1, 2, 3, 3, 6, 9);
        MutableLongPoly b = MutableLongPoly.create(1, 3, 6, 5, 3);
        assertEquals(MutableLongPoly.create(1, 2, 3), ModularGCD(a, b));

        a = MutableLongPoly.create(0, 0, 1, 2);
        b = MutableLongPoly.create(-1, 0, 4);
        assertEquals(MutableLongPoly.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test16() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(-1, 0, 4);
        MutableLongPoly b = MutableLongPoly.create(-1, 0, 0, 0, 16);
        MutableLongPoly gcd = ModularGCD(a, b);
        assertEquals(MutableLongPoly.create(-1, 0, 4), gcd);
    }

    @Test
    public void test17() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(-1, 0, 4);
        MutableLongPoly b = MutableLongPoly.create(1, -5, 6);
        MutableLongPoly gcd = ModularGCD(a, b);
        assertEquals(MutableLongPoly.create(-1, 2), gcd);
    }

    @Test
    public void testRandom4() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int passed = 0;
        int wins = 0;
        for (int i = 0; i < 10000; i++) {
            MutableLongPoly dividend = randomPoly(15, 1000, rnd);
            MutableLongPoly divider = randomPoly(10, 1000, rnd);
            double norm = -1;
            try {
                MutableLongPoly[] qr = pseudoDivideAndRemainder(dividend, divider);
                assertPseudoQuotientRemainder(dividend, divider, qr);
                norm = qr[0].norm();
                ++passed;
            } catch (ArithmeticException e) {}

            double normAdaptive = -1;
            try {
                MutableLongPoly[] qr = pseudoDivideAndRemainderAdaptive(dividend, divider);
                assertPseudoQuotientRemainder(dividend, divider, qr);
                normAdaptive = qr[0].norm();
                ++passed;
            } catch (ArithmeticException e) {}

            if (norm != -1) {
                assertTrue(normAdaptive != -1);
                assertTrue(normAdaptive <= norm);
                if (normAdaptive < norm)
                    ++wins;
            }
        }
        System.out.println(passed);
        System.out.println(wins);
    }

    @Test
    public void testRandom5() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int overflow = 0;
        int larger = 0;
        DescriptiveStatistics timings = null;
        for (int k = 0; k < 2; k++) {
            timings = new DescriptiveStatistics();
            for (int i = 0; i < 1000; i++) {
                try {
                    MutableLongPoly a = randomPoly(7, 100, rnd);
                    MutableLongPoly b = randomPoly(6, 100, rnd);
                    MutableLongPoly gcd = randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    MutableLongPoly mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    MutableLongPoly[] qr = pseudoDivideAndRemainderAdaptive(mgcd, gcd);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
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
    public void test18() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test19() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.multiply(powExact(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divide(powExact(5, poly.degree));
        poly.multiply(512);
        assertEquals(-654063683625L, poly.evaluateAtRational(17, 2));
    }

    @Test(expected = IllegalArgumentException.class)
    public void test20() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }

    @Test
    public void test21() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics
                fast = new DescriptiveStatistics(), fastPseudo = new DescriptiveStatistics(),
                gen = new DescriptiveStatistics(), genPseudo = new DescriptiveStatistics();

        for (int i = 0; i < 10_000; i++) {
            MutableLongPoly dividend = randomPoly(rndd.nextInt(1, 10), 10, rnd);
            MutableLongPoly divider;
            do {
                divider = MutableLongPoly.create(rndd.nextInt(-10, 10), 1);
            } while (divider.degree == 0);

            if (i == 100) {
                fast.clear();
                gen.clear();
            }

            long start = System.nanoTime();
            MutableLongPoly[] actual = divideAndRemainderLinearDivider(dividend, divider);
            fast.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            MutableLongPoly[] expected = divideAndRemainderGeneral0(dividend, divider, 1);
            gen.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);


            for (long modulus : new long[]{2, 3, 5, 7, 11, 13, 19, 101}) {
                do {
                    divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
                } while (divider.degree == 0 || gcd(divider.lc(), modulus) != 1);
                start = System.nanoTime();
                actual = divideAndRemainderLinearDividerModulus(dividend, divider, modulus);
                fast.addValue(System.nanoTime() - start);
                start = System.nanoTime();
                expected = divideAndRemainderModulus(dividend, divider, modulus);
                gen.addValue(System.nanoTime() - start);
                assertArrayEquals(expected, actual);
            }

            do {
                divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = pseudoDivideAndRemainderLinearDivider(dividend, divider);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            long factor = powExact(divider.lc(), dividend.degree - divider.degree + 1);
            expected = divideAndRemainderGeneral0(dividend, divider, factor);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);

            do {
                divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            expected = pseudoDivideAndRemainderAdaptive0(dividend, divider);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);
        }
        System.out.println("Fast:");
        System.out.println(fast);
        System.out.println("General:");
        System.out.println(gen);

        System.out.println("Fast:");
        System.out.println(fastPseudo);
        System.out.println("General:");
        System.out.println(genPseudo);
    }

    @Test
    public void test22() throws Exception {
        MutableLongPoly dividend = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly divider = MutableLongPoly.create(1, 2);
        System.out.println(Arrays.toString(divideAndRemainderLinearDivider(dividend, divider)));
        System.out.println(Arrays.toString(divideAndRemainder(dividend, divider)));
    }

    private static void assertGCD(MutableLongPoly a, MutableLongPoly b, MutableLongPoly gcd) {
        MutableLongPoly[] qr;
        for (MutableLongPoly dividend : Arrays.asList(a, b)) {
            qr = pseudoDivideAndRemainderAdaptive(dividend, gcd);
            assertNotNull(qr);
            assertTrue(qr[1].isZero());
        }
    }

    private static void assertQuotientRemainder(MutableLongPoly dividend, MutableLongPoly divider, MutableLongPoly[] qr) {
        if (qr == null) return;
        assertEquals(dividend, divider.clone().multiply(qr[0]).add(qr[1]));
    }

    private static void assertPseudoQuotientRemainder(MutableLongPoly dividend, MutableLongPoly divider, MutableLongPoly[] qr) {
        if (qr == null) return;
        MutableLongPoly d = divider.clone().multiply(qr[0]).add(qr[1]);
        MutableLongPoly[] factor = divideAndRemainder(d, dividend);
        assertNotNull(factor);
        assertTrue(factor[1].isZero());
        assertTrue(factor[0].isConstant());
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

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutableLongPoly a, MutableLongPoly b, PolynomialRemainders prs, long modulus) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        assertEquals(a.clone().modulus(modulus), prs.remainders.get(0));
        assertEquals(b.clone().modulus(modulus), prs.remainders.get(1));

        MutableLongPoly gcd = prs.gcd().clone();
        assertTrue(divideAndRemainder(a, gcd, modulus)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, modulus)[1].isZero());
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