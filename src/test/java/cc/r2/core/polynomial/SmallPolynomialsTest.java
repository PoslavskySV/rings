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
import static cc.r2.core.polynomial.LongArithmetics.pow;
import static cc.r2.core.polynomial.SmallPolynomialArithmetics.pow;
import static cc.r2.core.polynomial.SmallPolynomialArithmetics.powMod;
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
            for (int i = 0; i < 10000; i++) {
                try {
                    MutableLongPoly a = randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    MutableLongPoly b = randomPoly(1 + rnd.nextInt(6), 100, rnd);
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
        poly.multiply(pow(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divide(pow(5, poly.degree));
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
            long factor = pow(divider.lc(), dividend.degree - divider.degree + 1);
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
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.zero();
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
    public void test23() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(2);
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

    @Test(expected = ArithmeticException.class)
    public void test24() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(0);
        divideAndRemainder(a, b);
    }

    @Test
    public void test25() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(0);
        MutableLongPoly[] zeros = {MutableLongPoly.zero(), MutableLongPoly.zero()};
        assertArrayEquals(zeros, divideAndRemainder(b, a));
        assertArrayEquals(zeros, pseudoDivideAndRemainder(b, a));
        assertArrayEquals(zeros, pseudoDivideAndRemainderAdaptive(b, a));
        assertArrayEquals(zeros, divideAndRemainder(b, a, 13));
    }

    @Test
    public void test26() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(32);
        MutableLongPoly b = MutableLongPoly.create(24);
        MutableLongPoly gcd = MutableLongPoly.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, SubresultantEuclid(a, b).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, true).gcd());
        assertEquals(gcd, PolynomialEuclid(a, b, false).gcd());
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


    @Test
    public void test27() throws Exception {
        MutableLongPoly poly = pow(MutableLongPoly.create(1, 3).multiply(2), 3).multiply(pow(MutableLongPoly.create(-3, -5, 7), 2));
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutableLongPoly.create(1, 3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutableLongPoly.create(3);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutableLongPoly.create(33);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
        poly = MutableLongPoly.create(22, 22).multiply(MutableLongPoly.create(12, 12, 12)).multiply(12);
        assertFactorization(poly, SquareFreeFactorizationYun(poly));
    }

    @Test
    public void test28() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics yun = new DescriptiveStatistics(), musser = new DescriptiveStatistics();

        long overflows = 0;
        for (int i = 0; i < 1000; i++) {
            if (i == 100) {
                yun.clear();
                musser.clear();
            }
            int nbase = rndd.nextInt(1, 5);
            MutableLongPoly poly;
            try {
                poly = MutableLongPoly.create(rndd.nextLong(1, 10));
                for (int j = 0; j < nbase; j++) {
                    MutableLongPoly factor = randomPoly(rndd.nextInt(1, 3), 10, rnd);
                    int exponent = rndd.nextInt(1, 5);
                    poly = poly.multiply(pow(factor, exponent));
                }
            } catch (ArithmeticException e) {
                --i;
                continue;
            }
            try {
                long start = System.nanoTime();
                Factorization yunFactorization = SquareFreeFactorizationYun(poly);
                yun.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                Factorization musserFactorization = SquareFreeFactorizationMusser(poly);
                musser.addValue(System.nanoTime() - start);

                assertEquals(yunFactorization.factors.length, musserFactorization.factors.length);
                assertFactorization(poly, musserFactorization);
                assertFactorization(poly, yunFactorization);
            } catch (ArithmeticException exc) {
                if (!exc.getMessage().contains("overflow"))
                    throw exc;
                ++overflows;
            }
        }

        System.out.println("Overflows: " + overflows);
        System.out.println("Timings Yun:\n" + yun);
        System.out.println("\nTimings Musser:\n" + musser);
    }

    @Test
    public void test29() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, -1458, 6561, -6561);
        Factorization factorization = SquareFreeFactorizationYun(poly);
        assertFactorization(poly, factorization);
    }

    @Test
    public void test30() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, 2, -1, -2, -7);
        MutableLongPoly b = MutableLongPoly.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }

    static void assertFactorization(MutableLongPoly poly, Factorization factorization) {
        MutableLongPoly r = MutableLongPoly.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(pow(factorization.factors[i], factorization.exponents[i]));
        assertEquals(poly, r);
    }

    static void assertFactorization(MutableLongPoly poly, Factorization factorization, long modulus) {
        MutableLongPoly r = MutableLongPoly.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(powMod(factorization.factors[i], factorization.exponents[i], modulus), modulus);
        assertEquals(poly, r);
    }

    @Test
    public void test31() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 0, 1, 0, 2);
        poly = poly.multiply(MutableLongPoly.create(1, 0, 1, 0, 2), 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
    }

    @Test
    public void test32() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 0, 0, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
    }

    @Test
    public void test33() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics timings = new DescriptiveStatistics();
        long overflows = 0;
        long[] primes = {2, 3, 5, 7, 11, 13, 101};
        for (int i = 0; i < 1000; i++) {
            if (i == 100)
                timings.clear();
            int nbase = rndd.nextInt(1, 3);
            for (long modulus : primes) {
                MutableLongPoly poly;
                poly = MutableLongPoly.create(rndd.nextLong(1, 10)).modulus(modulus);
                for (int j = 0; j < nbase; j++) {
                    MutableLongPoly f = randomPoly(rndd.nextInt(1, 3), 10, rnd);
                    poly = poly.multiply(powMod(f, rndd.nextInt(1, 3), modulus), modulus);
                }
                try {
                    long start = System.nanoTime();
                    System.out.println(modulus);
                    System.out.println(poly);
                    Factorization factorization = SquareFreeFactorization(poly, modulus);
                    timings.addValue(System.nanoTime() - start);
                    assertFactorization(poly, factorization, modulus);
                } catch (ArithmeticException exc) {
                    if (!exc.getMessage().contains("overflow"))
                        throw exc;
                    ++overflows;
                }

            }
        }

        System.out.println("Overflows: " + overflows);
        System.out.println("Timings:\n" + timings);
    }

    @Test
    public void test33a() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, 1, 3, 4, 3, 3, 2, 3);
        assertFactorization(poly, SquareFreeFactorization(poly, 5), 5);
    }

    @Test
    public void test33b() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, 0, 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test33c() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, 0, 1, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 2), 2);
    }

    @Test
    public void test33d() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test33e() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(2, 3, 2, 1, 3, 3, 3);
        System.out.println(poly);
        System.out.println(SquareFreeFactorization(poly, 5));
        assertFactorization(poly, SquareFreeFactorization(poly, 5), 5);
    }

    @Test
    public void test33f() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1);
        assertFactorization(poly, SquareFreeFactorization(poly, 3), 3);
    }

    @Test
    public void test33g() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, 0, 0, 8, 20, 67, 55);
        System.out.println(poly.clone().monic(101));
        System.out.println(poly);
        System.out.println(SquareFreeFactorization(poly, 101));
        assertFactorization(poly, SquareFreeFactorization(poly, 101), 101);
    }
}