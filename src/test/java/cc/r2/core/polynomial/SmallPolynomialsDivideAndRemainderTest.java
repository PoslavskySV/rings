package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.polynomial.RandomPolynomials.randomPoly;
import static cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.*;
import static org.junit.Assert.*;

public class SmallPolynomialsDivideAndRemainderTest {
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
        MutableLongPoly dividend = MutableLongPoly.create(28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575, 0);
        MutableLongPoly divider = MutableLongPoly.create(24487310, 38204421, 12930314, 41553770, -1216266, 7382581, 15631547, 0, 0);
        pseudoDivideAndRemainder(dividend, divider, true);
    }

    @Test
    public void test3() throws Exception {
        MutableLongPoly dividend = MutableLongPoly.create(28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575);
        MutableLongPoly divider = MutableLongPoly.one();
        MutableLongPoly[] qr = divideAndRemainder(dividend, divider, true);
        assertEquals(dividend, qr[0]);
        assertTrue(qr[1].isZero());
    }

    @Test
    public void test4_ModularSmallPolynomialsRandom() throws Exception {
        // polynomials
        RandomGenerator rnd = new Well1024a();
        MutableLongPoly[] qd;
        long[] primes = {3, 5, 7, 11, 13, 67, 97, 113, 127};
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly dividend = randomPoly(rnd.nextInt(DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD_MONIC), rnd);
            MutableLongPoly divider = randomPoly(rnd.nextInt(DIVIDEND_DEGREE_FAST_DIVISION_THRESHOLD_MONIC), rnd);

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
    public void test4a() throws Exception {
        MutableLongPoly dividend = MutableLongPoly.create(95, 45, 67, 5, -2, 65, 24, 24, 60);
        MutableLongPoly divider = MutableLongPoly.create(94, 86);
        long prime = 7;

        MutableLongPoly[] qd = divideAndRemainder(dividend.clone(), divider, prime, false);
        assertQuotientRemainder(dividend, divider, qd, prime);
    }

    @Test
    public void test5_SmallPolynomialsRandom() throws Exception {
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
    public void test6_ModularSmallPolynomialsRemainderRandom() throws Exception {
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
    public void test7_LinearDividerRandom() throws Exception {
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
                } while (divider.degree == 0 || LongArithmetics.gcd(divider.lc(), modulus) != 1);
                start = System.nanoTime();
                actual = divideAndRemainderLinearDividerModulus(dividend, divider, modulus, true);
                fast.addValue(System.nanoTime() - start);
                start = System.nanoTime();
                expected = divideAndRemainderClassic0(dividend, divider, modulus, true);
                gen.addValue(System.nanoTime() - start);
                assertArrayEquals(expected, actual);

                actualNoCopy = divideAndRemainderLinearDividerModulus(dividend.clone(), divider, modulus, false);
                expectedNoCopy = divideAndRemainderClassic0(dividend.clone(), divider, modulus, false);
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


    @Test(expected = ArithmeticException.class)
    public void test8() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(0);
        divideAndRemainder(a, b, false);
    }


    @Test
    public void test9() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(8, -2 * 8, 8, 8 * 2);
        MutableLongPoly b = MutableLongPoly.create(0);
        MutableLongPoly[] zeros = {MutableLongPoly.zero(), MutableLongPoly.zero()};
        assertArrayEquals(zeros, divideAndRemainder(b, a, true));
        assertArrayEquals(zeros, pseudoDivideAndRemainder(b, a, true));
        assertArrayEquals(zeros, pseudoDivideAndRemainderAdaptive(b, a, true));
        assertArrayEquals(zeros, divideAndRemainder(b, a, 13, true));
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
        assertEquals(dividend.clone().modulus(modulus), divider.clone().multiply(qr[0], modulus).add(qr[1], modulus).modulus(modulus));
    }

    private static MutableLongPoly inverseModMonomial0(MutableLongPoly poly, int xDegree, long modulus) {
        if (xDegree < 1)
            return null;
        if (poly.cc() != 1)
            throw new IllegalArgumentException();
        int r = log2(xDegree);
        MutableLongPoly gPrev = MutableLongPoly.one();
        for (int i = 0; i < r; ++i) {
            MutableLongPoly tmp = gPrev.clone().multiply(2, modulus).subtract(gPrev.square(modulus).multiply(poly, modulus), modulus);
            gPrev = remainderMonomial(tmp, 1 << i, false);
        }
        return gPrev;
    }

    @Test
    public void test10_InverseModRandom() throws Exception {
        RandomGenerator rnd = new Well1024a();
        int modulus = 17;
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly f = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(100), modulus, rnd);
            f.data[0] = 1;
            int modDegree = 1 + rnd.nextInt(2 * f.degree);
            MutableLongPoly invmod = inverseModMonomial0(f, modDegree, modulus);
            assertInverseModMonomial(f, invmod, modDegree, modulus);
        }
    }

    static void assertInverseModMonomial(MutableLongPoly poly, MutableLongPoly invMod, int monomialDegree, long modulus) {
        Assert.assertTrue(SmallPolynomialArithmetics.polyMultiplyMod(poly, invMod, MutableLongPoly.createMonomial(1, monomialDegree), modulus, true).isOne());
    }

    @Test
    public void test11_InverseModStructureRandom() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long modulus = 101;
        for (int i = 0; i < 1000; i++) {
            MutableLongPoly p = RandomPolynomials.randomMonicPoly(2 + rnd.nextInt(100), modulus, rnd);
            p.data[0] = 1;

            InverseModMonomial invMod = fastDivisionPreConditioning(p, modulus);
            for (int j = 0; j < 30; j++) {
                int xDegree = 1 + rnd.nextInt(1025);
                Assert.assertEquals(invMod.getInverse(xDegree), inverseModMonomial0(p.clone().reverse(), xDegree, modulus));
            }
        }
    }

    @Test
    public void test12_FastDivisionRandom() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long modulus = 17;
        for (int i = 0; i < 300; i++) {
            MutableLongPoly b = RandomPolynomials.randomMonicPoly(30, modulus, rnd);
            MutableLongPoly a = RandomPolynomials.randomMonicPoly(rnd.nextInt(30), modulus, rnd);

            InverseModMonomial invMod = fastDivisionPreConditioning(b, modulus);
            MutableLongPoly[] fast = divideAndRemainderFast(a, b, invMod, modulus, true);
            MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
            Assert.assertArrayEquals(fast, plain);
        }
    }

    @Test
    public void test13() throws Exception {
        long modulus = 7;
        MutableLongPoly a = MutableLongPoly.create(5, 1, 4, 6, 4, 3, 5, 5, 3, 4, 2, 2, 5, 2, 5, 6, 1, 1, 2, 5, 1, 0, 0, 6, 6, 5, 5, 1, 0, 1, 4, 1, 1);
        MutableLongPoly b = MutableLongPoly.create(2, 5, 3, 1, 1, 5, 6, 3, 4, 0, 0, 5, 4, 0, 2, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(b, modulus);
        MutableLongPoly[] fast = divideAndRemainderFast(a, b, invMod, modulus, true);
        MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void test14() throws Exception {
        long modulus = 7;
        MutableLongPoly a = MutableLongPoly.create(5, 3, 3, 3, 5, 3, 1, 4, -3, 1, 4, 5, 0, 2, 2, -5, 1).modulus(modulus);
        MutableLongPoly b = MutableLongPoly.create(0, 4, 6, 1, 2, 4, 0, 0, 6, 5, 2, 3, 1, 4, 0, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(b, modulus);
        MutableLongPoly[] fast = divideAndRemainderFast(a, b, invMod, modulus, true);
        MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void test15() throws Exception {
        long modulus = 17;
        MutableLongPoly a = MutableLongPoly.create(0, 6, 2, 1, 10, 15, 16, 15, 2, 11, 13, 0, 1, 15, 5, 13, 8, 14, 13, 14, 15, 1, 1);
        MutableLongPoly b = MutableLongPoly.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(b, modulus);
        MutableLongPoly[] fast = divideAndRemainderFast(a, b, invMod, modulus, true);
        MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void test16() throws Exception {
        long modulus = 17;
        MutableLongPoly a = MutableLongPoly.create(5, 9, 4, 9, 8, 12, 11, 9, 1, 6, 15, 7, 11, 2, 11, 13, 11, 10, 5, 1);
        MutableLongPoly b = MutableLongPoly.create(11, 15, 9, 5, 11, 5, 14, 9, 1, 0, 16, 12, 11, 5, 15, 10, 15, 2, 14, 3, 1, 16, 16, 12, 13, 1, 12, 11, 1, 15, 1);
        InverseModMonomial invMod = fastDivisionPreConditioning(b, modulus);
        MutableLongPoly[] fast = divideAndRemainderFast(a, b, invMod, modulus, true);
        MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }

    @Test
    public void test17() throws Exception {
        long modulus = 7;
        MutableLongPoly f = MutableLongPoly.create(0, 2, 3, 4, -5, 1).reverse();
        int modDegree = f.degree;
        MutableLongPoly invmod = inverseModMonomial0(f, modDegree, modulus);
        MutableLongPoly r = SmallPolynomialArithmetics.polyMultiplyMod(f, invmod, MutableLongPoly.createMonomial(1, modDegree), modulus, true);
        Assert.assertTrue(r.isOne());
    }

    @Test
    public void test18() throws Exception {
        long modulus = 17;
        MutableLongPoly f = MutableLongPoly.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1).reverse();
        int modDegree = 9;
        MutableLongPoly invmod = inverseModMonomial0(f, modDegree, modulus);
        MutableLongPoly r = SmallPolynomialArithmetics.polyMultiplyMod(f, invmod, MutableLongPoly.createMonomial(1, modDegree), modulus, true);
        Assert.assertTrue(r.isOne());
    }

    @Test
    public void test19_FastDivisionPerformance() throws Exception {
        long modulus = LargeDDFTest.bigModulus;
        RandomGenerator rnd = new Well1024a();
        MutableLongPoly divider = RandomPolynomials.randomMonicPoly(300, modulus, rnd);

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        InverseModMonomial invRev = fastDivisionPreConditioning(divider, modulus);
        int nIterations = 1500;
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000) {
                classic.clear();
                fast.clear();
            }
            MutableLongPoly dividend = RandomPolynomials.randomPoly(3 * divider.degree / 2, (int) modulus, rnd);
            dividend = dividend.modulus(modulus);

            long start = System.nanoTime();
            MutableLongPoly[] qdPlain = divideAndRemainderClassic(dividend, divider, modulus, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            MutableLongPoly[] qdNewton = divideAndRemainderFast(dividend, divider, invRev, modulus, true);
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

    @Test
    public void test20_FastDivisionPerformance() throws Exception {
        long modulus = LargeDDFTest.bigModulus;
        RandomGenerator rnd = new Well1024a();

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        int nIterations = 15000;
        int dividerDegree = 36;
        int dividendDegree = 56;
        for (int i = 0; i < nIterations; i++) {
            if (i == 1000) {
                classic.clear();
                fast.clear();
            }

            MutableLongPoly divider = RandomPolynomials.randomMonicPoly(dividerDegree, modulus, rnd);
            MutableLongPoly dividend = RandomPolynomials.randomPoly(dividendDegree, (int) modulus, rnd);
            dividend = dividend.modulus(modulus);
            divider.multiply(3, modulus);

            long start = System.nanoTime();
            MutableLongPoly[] qdPlain = divideAndRemainderClassic(dividend, divider, modulus, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            MutableLongPoly[] qdNewton = divideAndRemainderFast(dividend, divider, modulus, true);
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

    @Test
    public void test21() throws Exception {
        long modulus = 17;
        MutableLongPoly a = MutableLongPoly.create(5, 9, 4, 9, 8, 12, 11, 9, 1, 6, 15, 7, 11, 2, 11, 13, 11, 10, 5, 1);
        MutableLongPoly b = MutableLongPoly.create(11, 15, 9, 5, 11, 5, 14, 9, 1, 0, 16, 12, 11, 5, 15, 10, 15, 2, 14, 3, 1, 16, 16, 12, 13, 1, 12, 11, 1, 15, 13);
        MutableLongPoly[] fast = divideAndRemainderFast(a, b, modulus, true);
        MutableLongPoly[] plain = divideAndRemainderClassic(a, b, modulus, true);
        Assert.assertArrayEquals(fast, plain);
    }
}