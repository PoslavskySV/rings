package cc.r2.core.polynomial;

import cc.r2.core.number.primes.SieveOfAtkin;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

import static cc.r2.core.polynomial.LongArithmetics.gcd;
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
        assertTrue(divideAndRemainder(a, gcd, modulus, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, modulus, true)[1].isZero());
    }

    @Test(expected = ArithmeticException.class)
    public void test2() throws Exception {
        //test long overflow
        MutableLongPoly dividend = new MutableLongPoly(new long[]{28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575, 0}, 7);
        MutableLongPoly divider = new MutableLongPoly(new long[]{24487310, 38204421, 12930314, 41553770, -1216266, 7382581, 15631547, 0, 0}, 6);
        pseudoDivideAndRemainder(dividend, divider, true);
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
        MutableLongPoly[] qr = divideAndRemainder(dividend, divider, true);
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
                qd = divideAndRemainder(dividend, divider, prime, true);
                assertQuotientRemainder(dividend, divider, qd, prime);
                qd = divideAndRemainder(dividend.clone(), divider, prime, false);
                try {
                    assertQuotientRemainder(dividend, divider, qd, prime);
                } catch (AssertionError err) {
                    System.out.println(dividend.toStringForCopy());
                    System.out.println(divider.toStringForCopy());
                    System.out.println(prime);
                    throw err;
                }
            }
        }
    }

    @Test
    public void test11a() throws Exception {
        MutableLongPoly dividend = MutableLongPoly.create(95, 45, 67, 5, -2, 65, 24, 24, 60);
        MutableLongPoly divider = MutableLongPoly.create(94, 86);
        long prime = 7;

        MutableLongPoly[] qd = divideAndRemainder(dividend.clone(), divider, prime, false);
        assertQuotientRemainder(dividend, divider, qd, prime);
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
                MutableLongPoly[] qr = pseudoDivideAndRemainder(dividend, divider, true);
                assertPseudoQuotientRemainder(dividend, divider, qr);
                qr = pseudoDivideAndRemainder(dividend.clone(), divider, false);
                assertPseudoQuotientRemainder(dividend, divider, qr);
                norm = qr[0].norm();
                ++passed;
            } catch (ArithmeticException e) {}

            double normAdaptive = -1;
            try {
                MutableLongPoly[] qr = pseudoDivideAndRemainderAdaptive(dividend, divider, true);
                assertPseudoQuotientRemainder(dividend, divider, qr);
                qr = pseudoDivideAndRemainderAdaptive(dividend.clone(), divider, false);
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

                    MutableLongPoly[] qr = pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    MutableLongPoly[] qr1 = pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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
    public void testRemainderRandom1() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long[] primes = {2, 3, 5, 7, 11, 17, 41, 43, 101};
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly dividend = randomPoly(5 + rnd.nextInt(20), 100, rnd);
            MutableLongPoly divider = randomPoly(rnd.nextInt(20), 100, rnd);
            if (dividend.degree < divider.degree) {
                MutableLongPoly tmp = dividend;
                dividend = divider;
                divider = tmp;
            }
            for (long prime : primes) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                MutableLongPoly expected = divideAndRemainder(dividend, divider, prime, true)[1];
                assertEquals(expected, remainder(dividend, divider, prime, true));
                assertEquals(expected, remainder(dividend.clone(), divider, prime, false));
            }
        }
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
        poly.multiply(LongArithmetics.pow(5, poly.degree));

        assertEquals(3045900, poly.evaluateAtRational(1, 5));
        assertEquals(-74846713560L, poly.evaluateAtRational(13, 5));
        assertEquals(40779736470L, poly.evaluateAtRational(13, -5));

        poly.divide(LongArithmetics.pow(5, poly.degree));
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
            MutableLongPoly[] actual = divideAndRemainderLinearDivider(dividend, divider, true);
            fast.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            MutableLongPoly[] expected = divideAndRemainderGeneral0(dividend, divider, 1, true);
            gen.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);
            MutableLongPoly[] expectedNoCopy = divideAndRemainderGeneral0(dividend.clone(), divider, 1, false);
            MutableLongPoly[] actualNoCopy = divideAndRemainderLinearDivider(dividend.clone(), divider, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);


            for (long modulus : new long[]{2, 3, 5, 7, 11, 13, 19, 101}) {
                do {
                    divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
                } while (divider.degree == 0 || gcd(divider.lc(), modulus) != 1);
                start = System.nanoTime();
                actual = divideAndRemainderLinearDividerModulus(dividend, divider, modulus, true);
                fast.addValue(System.nanoTime() - start);
                start = System.nanoTime();
                expected = divideAndRemainderModulus(dividend, divider, modulus, true);
                gen.addValue(System.nanoTime() - start);
                assertArrayEquals(expected, actual);

                actualNoCopy = divideAndRemainderLinearDividerModulus(dividend.clone(), divider, modulus, false);
                expectedNoCopy = divideAndRemainderModulus(dividend.clone(), divider, modulus, false);
                assertArrayEquals(expected, actualNoCopy);
                assertArrayEquals(expected, expectedNoCopy);
            }

            do {
                divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = pseudoDivideAndRemainderLinearDivider(dividend, divider, true);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            long factor = LongArithmetics.pow(divider.lc(), dividend.degree - divider.degree + 1);
            expected = divideAndRemainderGeneral0(dividend, divider, factor, true);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);

            actualNoCopy = pseudoDivideAndRemainderLinearDivider(dividend.clone(), divider, false);
            expectedNoCopy = divideAndRemainderGeneral0(dividend.clone(), divider, factor, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);

            do {
                divider = MutableLongPoly.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, true);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            expected = pseudoDivideAndRemainderAdaptive0(dividend, divider, true);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);

            actualNoCopy = pseudoDivideAndRemainderLinearDividerAdaptive(dividend.clone(), divider, false);
            expectedNoCopy = pseudoDivideAndRemainderAdaptive0(dividend.clone(), divider, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);
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
        divideAndRemainder(a, b, false);
    }

    @Test
    public void test25() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(0);
        MutableLongPoly[] zeros = {MutableLongPoly.zero(), MutableLongPoly.zero()};
        assertArrayEquals(zeros, divideAndRemainder(b, a, true));
        assertArrayEquals(zeros, pseudoDivideAndRemainder(b, a, true));
        assertArrayEquals(zeros, pseudoDivideAndRemainderAdaptive(b, a, true));
        assertArrayEquals(zeros, divideAndRemainder(b, a, 13, true));
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
            qr = pseudoDivideAndRemainderAdaptive(dividend, gcd, true);
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
        MutableLongPoly[] factor = divideAndRemainder(d, dividend, true);
        assertNotNull(factor);
        assertTrue(factor[1].isZero());
        assertTrue(factor[0].isConstant());
    }

    private static void assertQuotientRemainder(MutableLongPoly dividend, MutableLongPoly divider, MutableLongPoly[] qr, long modulus) {
        if (qr == null) return;
        assertEquals(dividend.clone().modulus(modulus), divider.clone().multiply(qr[0], modulus).add(qr[1], modulus));
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutableLongPoly a, MutableLongPoly b, PolynomialRemainders prs) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        a = a.clone().divide(content(a));
        b = b.clone().divide(content(b));
        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        MutableLongPoly gcd = primitivePart(prs.gcd().clone());
        assertTrue(pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }

    @SuppressWarnings("ConstantConditions")
    private static void assertPolynomialRemainders(MutableLongPoly a, MutableLongPoly b, PolynomialRemainders prs, long modulus) {
        if (a.degree < b.degree)
            assertPolynomialRemainders(b, a, prs);

        assertEquals(a.clone().modulus(modulus), prs.remainders.get(0));
        assertEquals(b.clone().modulus(modulus), prs.remainders.get(1));

        MutableLongPoly gcd = prs.gcd().clone();
        assertTrue(divideAndRemainder(a, gcd, modulus, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, modulus, true)[1].isZero());
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
        MutableLongPoly poly = pow(MutableLongPoly.create(1, 3).multiply(2), 3, false).multiply(pow(MutableLongPoly.create(-3, -5, 7), 2, false));
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
                    poly = poly.multiply(pow(factor, exponent, true));
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
            r = r.multiply(pow(factorization.factors[i], factorization.exponents[i], true));
        assertEquals(poly, r);
    }

    static void assertFactorization(MutableLongPoly poly, Factorization factorization, long modulus) {
        MutableLongPoly r = MutableLongPoly.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(powMod(factorization.factors[i], factorization.exponents[i], modulus, true), modulus);
        assertEquals(poly.clone().modulus(modulus), r);
    }

    static void assertDistinctDegreeFactorization(MutableLongPoly poly, Factorization factorization, long modulus) {
        MutableLongPoly r = MutableLongPoly.create(factorization.factor);
        for (int i = 0; i < factorization.factors.length; i++)
            r = r.multiply(factorization.factors[i], modulus);
        assertEquals(poly.clone().modulus(modulus), r);
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
        final class test {
            long overflows = 0;
            final DescriptiveStatistics arithmetics = new DescriptiveStatistics();
            final DescriptiveStatistics timings = new DescriptiveStatistics();

            void run(int bound, int maxDegree, int maxNBase, int maxExponent, int nIterations) {
                RandomGenerator rnd = new Well1024a();
                RandomDataGenerator rndd = new RandomDataGenerator(rnd);
                long[] primes = {2, 3, 5, 7, 11, 13, 101};
                long start;
                for (int i = 0; i < nIterations; i++) {
                    if (i == 100)
                        timings.clear();
                    int nbase = rndd.nextInt(1, maxNBase);
                    for (long modulus : primes) {
                        MutableLongPoly poly;
                        start = System.nanoTime();
                        poly = MutableLongPoly.create(rndd.nextLong(1, 1000)).modulus(modulus);
                        for (int j = 0; j < nbase; j++) {
                            MutableLongPoly f = randomPoly(rndd.nextInt(1, maxDegree), bound, rnd);
                            f = f.modulus(modulus);
                            poly = poly.multiply(powMod(f, rndd.nextInt(1, maxExponent), modulus, true), modulus);
                        }
                        arithmetics.addValue(System.nanoTime() - start);
                        try {
                            start = System.nanoTime();
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
            }
        }
        test smallPolys = new test();
        smallPolys.run(10, 3, 5, 5, 1000);
        System.out.println("Overflows: " + smallPolys.overflows);
        System.out.println("Timings:\n" + smallPolys.timings);
        System.out.println("Arithmetics:\n" + smallPolys.arithmetics);

        test largePolys = new test();
        largePolys.run(1000, 15, 10, 20, 10);
        System.out.println("Overflows: " + largePolys.overflows);
        System.out.println("Timings:\n" + largePolys.timings);
        System.out.println("Arithmetics:\n" + largePolys.arithmetics);
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
        assertFactorization(poly, SquareFreeFactorization(poly, 101), 101);
    }

    @Test
    public void test34() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(0, -1, -1, -1, 0, 1, -1, 1, 1);
        Factorization f = DistinctDegreeFactorization(poly, 3);
        assertDistinctDegreeFactorization(poly, f, 3);
    }

    @Test
    public void test35() throws Exception {
        final class test {
            long overflows = 0;
            long nonSquareFree = 0;
            long iterations = 0;
            final DescriptiveStatistics arithmetics = new DescriptiveStatistics();
            final DescriptiveStatistics timings = new DescriptiveStatistics();
            final DescriptiveStatistics nFactors = new DescriptiveStatistics();
            long[] primes = {2, 3, 5, 7, 11, 13, 101};


            void run(int bound, int minNBase, int maxNBase, int nIterations, boolean ensureSquareFree) {

                RandomGenerator rnd = new Well1024a();
                RandomDataGenerator rndd = new RandomDataGenerator(rnd);

                long start;
                out:
                for (int i = 0; i < nIterations; i++) {
                    if (i == 100)
                        timings.clear();
                    int nFactors = rndd.nextInt(minNBase, maxNBase);
                    for (long modulus : primes) {
                        MutableLongPoly poly;
                        start = System.nanoTime();
                        poly = MutableLongPoly.create(rndd.nextLong(1, 1000)).modulus(modulus);
                        List<MutableLongPoly> factors = new ArrayList<>();
                        int nBases = rndd.nextInt(minNBase, maxNBase);
                        for (int j = 1; j <= nBases; ++j) {
                            MutableLongPoly f;
                            do { f = randomPoly(j, bound, rnd).modulus(modulus);} while (f.degree != j);
                            factors.add(f);
                            poly = poly.multiply(f, modulus);
                        }
                        arithmetics.addValue(System.nanoTime() - start);
                        if (poly.isConstant()) {
                            --i;
                            continue out;
                        }
                        boolean squareFree = false;
                        List<MutableLongPoly> polynomials = new ArrayList<>();
                        if (!(squareFree = isSquareFree(poly, modulus))) {
                            ++nonSquareFree;
                            if (ensureSquareFree) {
                                Factorization factorization = SquareFreeFactorization(poly, modulus);
                                polynomials.addAll(Arrays.asList(factorization.factors));
                            }
                        }
                        if (polynomials.isEmpty()) polynomials.add(poly);

                        for (MutableLongPoly toTest : polynomials) {
                            try {
                                toTest = toTest.multiply(rndd.nextInt(1, 1000), modulus);
                                start = System.nanoTime();
                                Factorization factorization = DistinctDegreeFactorization(toTest, modulus);
                                timings.addValue(System.nanoTime() - start);
                                this.nFactors.addValue(factorization.factors.length);
                                assertDistinctDegreeFactorization(toTest, factorization, modulus);
                            } catch (ArithmeticException exc) {
                                if (!exc.getMessage().contains("overflow"))
                                    throw exc;
                                ++overflows;
                            }
                        }
                        ++iterations;
                    }
                }
            }

            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();
                sb.append("Total iterations: " + iterations).append("\n");
                sb.append("Overflows: " + overflows).append("\n");
                sb.append("Not square-free: " + nonSquareFree).append("\n");
                sb.append("\nNumber of factors:\n" + nFactors).append("\n");
                sb.append("\nTimings:\n" + timings).append("\n");
                sb.append("\nArithmetics:\n" + arithmetics).append("\n");
                return sb.toString();
            }
        }

        test smallPolys = new test();
        smallPolys.run(10, 1, 5, 10000, false);
        System.out.println(smallPolys);

        test largePolys = new test();
        largePolys.primes = new long[]{29, 31, 89, 101, 107, 139, 223};
        largePolys.run(1000, 8, 10, 20, true);
        System.out.println(largePolys);
    }

    @Test
    public void test36() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        class test {
            final DescriptiveStatistics timings = new DescriptiveStatistics();
            DescriptiveStatistics nFactors = new DescriptiveStatistics();

            void run(int maxDegree, int nIterations, int modBound) {
                SieveOfAtkin sieve = SieveOfAtkin.createSieve(modBound);
                for (int i = 0; i < nIterations; i++) {
                    if (i == 101)
                        timings.clear();
                    int modulus = sieve.randomPrime(rnd);
                    MutableLongPoly poly = randomPoly(rndd.nextInt(0, maxDegree), 10000, rnd).modulus(modulus);
                    long start = System.nanoTime();
                    Factorization factorization = DistinctDegreeFactorization(poly, modulus);
                    timings.addValue(System.nanoTime() - start);
                    nFactors.addValue(factorization.factors.length);
                    assertDistinctDegreeFactorization(poly, factorization, modulus);
                }
            }

            @Override
            public String toString() {
                StringBuilder sb = new StringBuilder();
                sb.append("\nNumber of factors:\n" + nFactors).append("\n");
                sb.append("\nTimings:\n" + timings).append("\n");
                return sb.toString();
            }
        }
        test small = new test();
        small.run(10, 1000, 1000);
        System.out.println(small);
        test big = new test();
        big.run(100, 100, 20);
        System.out.println(big);
    }

    @Test
    public void test35a() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(19, 20, 13, 10, 26, 19, 6, 29, 13, 20, 10, 12, 20, 3, 21, 16, 25, 10, 26, 22, 25, 2, 23, 29, 21, 14, 8, 26, 16, 7, 7, 1);
        Factorization factorization = DistinctDegreeFactorization(poly, 31);
        assertEquals(5, factorization.factors.length);
        assertDistinctDegreeFactorization(poly, factorization, 31);
    }

    @Test
    public void test37a() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(9, 7, 2, 3, 10, 1, 1);
        Factorization factorization = DistinctDegreeFactorization(poly, 11);
        assertDistinctDegreeFactorization(poly, factorization, 11);
    }

    @Test
    public void test38() throws Exception {
        MutableLongPoly poly = MutableLongPoly.create(172, 85, 84, 151, 122, 53, 107, 117, 82, 152, 133, 151, 178, 1);
        Factorization fct = DistinctDegreeFactorizationComplete(poly, 181);
        assertFactorization(poly, fct, 181);
    }

    @Test
    public void test39() throws Exception {
        assertEquals(1, SquareFreeFactorization(MutableLongPoly.create(0, 0, 0, 0, 1)).factors.length);
        assertEquals(1, SquareFreeFactorization(MutableLongPoly.create(0, 0, 0, 0, 1), 31).factors.length);
    }

    static void assertMMATest(int nLines, String file) throws Exception {
        InputStream resource = new GZIPInputStream(FactorizationTestDataTest.class.getClassLoader()
                .getResourceAsStream(file));
        int nEntries = 0;
        DescriptiveStatistics stats = new DescriptiveStatistics();
        for (FactorizationTestData fct : FactorizationTestData.allMod(resource)) {
            if (nEntries == 100)
                stats.clear();
            SmallPolynomials.Factorization expected = fct.factorization;
            long start = System.nanoTime();
            System.out.println(Arrays.toString(fct.poly.data));
            SmallPolynomials.Factorization actual = SmallPolynomials.DistinctDegreeFactorizationComplete(fct.poly, fct.modulus);
            stats.addValue(System.nanoTime() - start);
            Assert.assertEquals("Modulus: " + fct.modulus, expected.canonical(), actual.canonical());
            ++nEntries;
        }
        assertTrue(nEntries > 0);
        if (nLines != -1)
            assertEquals(nLines, nEntries);
        System.out.println(stats);
    }
//
//    @Ignore
//    @Test
//    public void test40() throws Exception {
//        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
//    }
//
//    @Ignore
//    @Test
//    public void test41() throws Exception {
//        assertMMATest(100_000, "cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.gz");
//        assertMMATest(1000, "cc/r2/core/polynomial/DistinctDegreeFactorizationLarge.gz");
//        assertMMATest(100, "cc/r2/core/polynomial/DistinctDegreeFactorizationHuge.gz");
//    }
//
//    @Test
//    public void testasdasd() throws Exception {
//        MutableLongPoly poly = MutableLongPoly.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1);
//        int modulus = 5659;
//        long start = System.nanoTime();
//        System.out.println(DistinctDegreeFactorization(poly, modulus));
//        System.out.println(System.nanoTime() - start);
//
//        start = System.nanoTime();
//        System.out.println(DistinctDegreeFactorization(poly, modulus));
//        System.out.println(System.nanoTime() - start);
//    }
//
//
//    static long[][] qMatrix(MutableLongPoly poly, long modulus) {
//        int pDegree = poly.degree;
//        long[][] matrix = new long[pDegree][pDegree];
//        long[] prevRow = new long[pDegree], nextRow = new long[pDegree];
//        prevRow[0] = 1;
//        matrix[0] = prevRow.clone();
//        for (int i = 1; i <= (pDegree - 1) * modulus; i++) {
//            nextRow[0] = symMod(-prevRow[pDegree - 1] * poly.data[0], modulus);
//            for (int j = 1; j < poly.degree; j++) {
//                nextRow[j] = symMod(prevRow[j - 1] - prevRow[pDegree - 1] * poly.data[j], modulus);
//            }
//            if (i % modulus == 0)
//                matrix[i / (int) modulus] = nextRow.clone();
//            long[] tmp = prevRow;
//            prevRow = nextRow;
//            nextRow = tmp;
//        }
//        return matrix;
//    }
//
//    static String toStringMatrix(long[][] matrix) {
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < matrix.length; i++) {
//            sb.append(Arrays.toString(matrix[i])).append("\n");
//        }
//        return sb.toString();
//    }
//
//    @Test
//    public void name() throws Exception {
//        long[][] expected = {
//                {1, 0, 0, 0, 0, 0},
//                {3, 5, -3, -3, -5, 5},
//                {3, -5, -5, 1, -1, 0},
//                {-2, 4, -1, 3, -4, -2},
//                {-4, -3, -1, 0, 0, -3},
//                {-3, -1, -4, -3, -1, -3}
//        };
//        MutableLongPoly poly = MutableLongPoly.create(1, -3, -1, -3, 1, -3, 1);
//        int modulus = 11;
//        long[][] qMatrix = qMatrix(poly, modulus);
//        System.out.println(toStringMatrix(qMatrix));
//        assertArrayEquals(expected, qMatrix);
//    }
//
//    @Test
//    public void wrer() throws Exception {
//        MutableLongPoly poly = MutableLongPoly.create(6650, 68859, 22275, 45078, 86304, 9759, 77160, 70073, 41899, 52881, 62889, 58468, 35826, 60356, 67213, 66957, 48370, 17669, 9933, 85458, 3134, 76771, 30441, 33067, 35939, 15710, 2403, 8585, 55218, 72652, 23952, 85278, 92366, 81522, 47437, 32453, 19760, 5051, 84527, 55625, 38211, 18165, 38887, 94661, 4046, 88205, 91932, 42789, 41182, 33497, 57403, 82501, 35133, 2346, 35376, 92459, 69637, 50572, 31966, 5279, 33814, 11215, 30244, 39497, 82716, 36040, 25972, 16361, 88885, 89514, 66641, 78008, 88470, 51393, 5626, 54147, 24953, 48299, 77990, 74869, 22067, 94204, 11658, 30396, 61221, 28882, 24978, 11737, 79083, 52379, 45547, 7482, 89156, 84783, 13140, 38412, 10110, 72974, 74516, 75284, 25327, 66808, 54726, 3462, 53452, 56885, 5921, 68793, 33047, 39883, 49840, 67584, 13360, 43291, 19317, 39530, 5922, 39463, 86786, 15846, 21785, 40463, 83277, 74177, 41218, 14196, 51191, 43599, 23830, 87613, 1414, 27672, 32990, 81745, 52957, 27855, 71616, 93334, 65928, 8242, 92984, 8345, 17228, 59512, 35349, 28330, 19021, 39366, 85001, 22699, 10186, 27312, 42484, 62155, 65370, 14172, 68282, 61633, 10726, 84239, 66430, 15752, 90164, 81410, 79784, 5751, 45762, 78313, 27020, 37809, 2897, 15129, 14970, 24014, 81092, 53643, 88663, 42889, 84295, 18189, 59806, 91795, 88777, 50017, 38189, 41721, 50622, 89687, 54431, 54986, 20530, 68806, 44449, 62479, 34149, 55409, 59757, 54592, 3636, 22578, 36217, 22896, 38901, 38573, 68767, 38291, 13457, 64421, 28767, 16589, 51589, 12948, 45939, 26680, 48003, 43471, 7013, 37294, 25587, 51820, 65336, 25703, 93341, 59022, 76069, 48653, 41795, 41392, 48562, 26240, 76813, 76274, 3876, 56383, 57752, 24556, 76413, 87271, 84231, 67364, 49656, 59996, 20438, 66506, 43313, 57650, 80206, 36887, 17852, 77602, 81316, 61562, 33687, 78752, 43969, 73349, 65202, 10234, 10062, 51956, 87369, 66206, 82705, 70217, 74172, 34581, 94543, 7664, 24364, 18110, 66448, 1);
//        int modulus = 5659;
//        System.out.println(poly.modulus(modulus));
//        MutableLongPoly xq = mod(MutableLongPoly.createMonomial(1, modulus), poly, modulus);
//        MutableLongPoly xq2 = composition(xq, xq, poly, modulus);
//        System.out.println(xq2);
//        System.out.println(mod(MutableLongPoly.createMonomial(1, modulus * modulus), poly, modulus));
//        MutableLongPoly xq3 = composition(xq2, xq, poly, modulus);
//        System.out.println(xq3);
//        System.out.println(mod(MutableLongPoly.createMonomial(1, modulus * modulus * modulus), poly, modulus));
//        MutableLongPoly xq4 = composition(xq2, xq2, poly, modulus);
//        System.out.println(xq4);
//        System.out.println(mod(MutableLongPoly.createMonomial(1, modulus * modulus * modulus * modulus), poly, modulus));
//
//    }
//
//    static MutableLongPoly composition(MutableLongPoly a, MutableLongPoly b, MutableLongPoly polyModulus, long modulus) {
//        MutableLongPoly res = MutableLongPoly.zero();
////        for (int i = a.degree; i >= 0; --i) {
////            res = LongArithmetics.add(LongArithmetics.multiply(res, point), data[i]);
////        }
////        return res;
//
//        for (int i = 0; i <= a.degree; i++) {
//            res.add(powMod(b, i, polyModulus, modulus).multiply(a.data[i], modulus), modulus);
//        }
//        return mod(res, polyModulus, modulus);
//    }

    //
//    @Test
//    public void test39() throws Exception {
//        InputStream resource = FactorizationTestDataTest.class.getClassLoader().getResourceAsStream("cc/r2/core/polynomial/DistinctDegreeFactorizationSmall.txt");
//        assertNotNull(resource);
//        for (FactorizationTestData fct : FactorizationTestData.allMod(resource)) {
//            Factorization factorization = DistinctDegreeFactorizationComplete(fct.poly, fct.modulus);
//
//        }
//    }
}