package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FiniteField;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.test.Benchmark;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.Arrays;

import static cc.redberry.rings.poly.univar.RandomUnivariatePolynomials.randomPoly;
import static cc.redberry.rings.poly.univar.UnivariateDivision.*;
import static cc.redberry.rings.poly.univar.UnivariatePolynomial.asOverZp64;
import static org.junit.Assert.*;


/**
 * Created by poslavsky on 15/02/2017.
 */
public class UnivariateDivisionTest extends AUnivariateTest {
    @Test
    @SuppressWarnings("ConstantConditions")
    public void test1() throws Exception {
        long modulus = 11;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216).modulus(modulus);

        UnivariateGCD.PolynomialRemainders<UnivariatePolynomialZp64> prs = UnivariateGCD.EuclidRemainders(a, b);
        UnivariatePolynomialZp64 gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(UnivariateDivision.divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(UnivariateDivision.divideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test(expected = ArithmeticException.class)
    public void test2() throws Exception {
        //test long overflow
        UnivariatePolynomialZ64 dividend = UnivariatePolynomialZ64.create(28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575, 0);
        UnivariatePolynomialZ64 divider = UnivariatePolynomialZ64.create(24487310, 38204421, 12930314, 41553770, -1216266, 7382581, 15631547, 0, 0);
        pseudoDivideAndRemainder(dividend, divider, true);
    }

    @Test
    public void test3() throws Exception {
        UnivariatePolynomialZ64 dividend = UnivariatePolynomialZ64.create(28130, 95683, 93697, 176985, 135507, 101513, 75181, 17575);
        UnivariatePolynomialZ64 divider = UnivariatePolynomialZ64.one();
        UnivariatePolynomialZ64[] qr = divideAndRemainder(dividend, divider, true);
        assertEquals(dividend, qr[0]);
        assertTrue(qr[1].isZero());
    }

    @Test
    public void test4_ModularSmallPolynomialsRandom() throws Exception {
        int thr = 81;
        // polynomials
        RandomGenerator rnd = getRandom();
        UnivariatePolynomialZp64[] qd;
        for (int i = 0; i < its(1000, 10_000); i++) {
            UnivariatePolynomialZ64 dividend = randomPoly(rnd.nextInt(thr), rnd);
            UnivariatePolynomialZ64 divider = randomPoly(rnd.nextInt(thr), rnd);

            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;

                UnivariatePolynomialZp64 a = dividend.modulus(prime, true);
                UnivariatePolynomialZp64 b = divider.modulus(prime, true);
                try {
                    qd = UnivariateDivision.divideAndRemainder(a, b, true);
                    assertQuotientRemainder(a, b, qd);
                    qd = UnivariateDivision.divideAndRemainder(a.clone(), b, false);

                    assertQuotientRemainder(a, b, qd);
                } catch (Exception err) {
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
        long prime = 7;
        UnivariatePolynomialZp64 dividend = UnivariatePolynomialZ64.create(95, 45, 67, 5, -2, 65, 24, 24, 60).modulus(prime);
        UnivariatePolynomialZp64 divider = UnivariatePolynomialZ64.create(94, 86).modulus(prime);

        UnivariatePolynomialZp64[] qd = UnivariateDivision.divideAndRemainder(dividend.clone(), divider, false);
        assertQuotientRemainder(dividend, divider, qd);
    }

    @Test
    public void test5_SmallPolynomialsRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        int passed = 0;
        int wins = 0;
        for (int i = 0; i < its(1000, 10_000); i++) {
            UnivariatePolynomialZ64 dividend = randomPoly(15, 1000, rnd);
            UnivariatePolynomialZ64 divider = randomPoly(10, 1000, rnd);
            double norm = -1;
            try {
                try {
                    UnivariatePolynomialZ64[] qr = pseudoDivideAndRemainder(dividend, divider, true);
                    assertPseudoQuotientRemainder(dividend, divider, qr);
                    qr = UnivariateDivision.pseudoDivideAndRemainder(dividend.clone(), divider, false);
                    assertPseudoQuotientRemainder(dividend, divider, qr);
                    norm = qr[0].norm2();
                    ++passed;
                } catch (ArithmeticException e) {}

                double normAdaptive = -1;
                try {
                    UnivariatePolynomialZ64[] qr = UnivariateDivision.pseudoDivideAndRemainderAdaptive(dividend, divider, true);
                    assertPseudoQuotientRemainder(dividend, divider, qr);
                    qr = UnivariateDivision.pseudoDivideAndRemainderAdaptive(dividend.clone(), divider, false);
                    assertPseudoQuotientRemainder(dividend, divider, qr);
                    normAdaptive = qr[0].norm2();
                    ++passed;
                } catch (ArithmeticException e) {}


                if (norm != -1) {
                    assertTrue(normAdaptive != -1);
                    assertTrue(normAdaptive <= norm);
                    if (normAdaptive < norm)
                        ++wins;
                }
            } catch (Exception | AssertionError e) {
                System.out.println(dividend.toStringForCopy());
                System.out.println(divider.toStringForCopy());
                throw e;
            }
        }
        System.out.println(passed);
        System.out.println(wins);
    }

    @Test
    public void test5_SmallPolynomialsRandom_a() throws Exception {
        UnivariatePolynomialZ64 dividend = UnivariatePolynomialZ64.create(1, 4, -5, -2, 9, 4, -5, 7, 5, -5, 6, 3, 9, 8, 9, -8);
        UnivariatePolynomialZ64 divider = UnivariatePolynomialZ64.create(7, 6, -1, 5, 0, 1, 0, 0, 8, 3, 7);
        UnivariatePolynomialZ64[] qr = pseudoDivideAndRemainder(dividend, divider, true);
        assertPseudoQuotientRemainder(dividend, divider, qr);
    }

    @Test
    public void test5_SmallPolynomialsRandom_b() throws Exception {
        UnivariatePolynomialZ64 dividend = UnivariatePolynomialZ64.create(6, 9, 9, 3, 2, 6, -6, 7, -2, 8, 4, -8, 7, 3, 3, -6);
        UnivariatePolynomialZ64 divider = UnivariatePolynomialZ64.create(0, 0, 0, 3, 1, 8, 8, -5, 6, -7, 9);
        UnivariatePolynomialZ64[] qr = pseudoDivideAndRemainder(dividend, divider, true);
        assertPseudoQuotientRemainder(dividend, divider, qr);
    }

    @Test
    public void test5_SmallPolynomialsRandom_с() throws Exception {
        UnivariatePolynomialZ64 dividend = UnivariatePolynomialZ64.create(5, -6, 0, -9, 3, -1, -5, 8, 0, 1, -1, -8, 8, 2, 0, 2);
        UnivariatePolynomialZ64 divider = UnivariatePolynomialZ64.create(1, 7, 6, -7, 9, 3, 3, 6, 3, -7, -3);
        UnivariatePolynomialZ64[] qr = pseudoDivideAndRemainder(dividend, divider, true);
        UnivariatePolynomialZ64 d = divider.clone().multiply(qr[0]).add(qr[1]);
        assertPseudoQuotientRemainder(dividend, divider, qr);
    }

    @Test
    public void test6_ModularSmallPolynomialsRemainderRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(1000, 10_000); i++) {
            UnivariatePolynomialZ64 dividend = randomPoly(5 + rnd.nextInt(20), 100, rnd);
            UnivariatePolynomialZ64 divider = randomPoly(rnd.nextInt(20), 100, rnd);
            if (dividend.degree < divider.degree) {
                UnivariatePolynomialZ64 tmp = dividend;
                dividend = divider;
                divider = tmp;
            }
            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                UnivariatePolynomialZp64 a = dividend.clone().modulus(prime);
                UnivariatePolynomialZp64 b = divider.clone().modulus(prime);
                UnivariatePolynomialZp64 expected = UnivariateDivision.divideAndRemainder(a, b, true)[1];
                assertEquals(expected, UnivariateDivision.remainder(a, b, true));
                assertEquals(expected, UnivariateDivision.remainder(a.clone(), b, false));
            }
        }
    }

    @Test
    public void test7_LinearDividerRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        DescriptiveStatistics
                fast = new DescriptiveStatistics(), fastPseudo = new DescriptiveStatistics(),
                gen = new DescriptiveStatistics(), genPseudo = new DescriptiveStatistics();

        long nIterations = its(1000, 15_000);
        out:
        for (int i = 0; i < nIterations; i++) {
            UnivariatePolynomialZ64 dividend = randomPoly(rndd.nextInt(1, 10), 10, rnd);
            UnivariatePolynomialZ64 divider;
            do {
                divider = UnivariatePolynomialZ64.create(rndd.nextInt(-10, 10), 1);
            } while (divider.degree == 0);

            if (i == nIterations / 10)
                Arrays.asList(fast, fastPseudo, gen, genPseudo).forEach(DescriptiveStatistics::clear);

            long start = System.nanoTime();
            UnivariatePolynomialZ64[] actual = UnivariateDivision.divideAndRemainderLinearDivider(dividend, divider, true);
            fast.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            UnivariatePolynomialZ64[] expected = UnivariateDivision.divideAndRemainderClassic0(dividend, divider, 1, true);
            gen.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);
            UnivariatePolynomialZ64[] expectedNoCopy = divideAndRemainderClassic0(dividend.clone(), divider, 1, false);
            UnivariatePolynomialZ64[] actualNoCopy = divideAndRemainderLinearDivider(dividend.clone(), divider, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);


            for (long modulus : getSmallModulusArray(10)) {
                do {
                    divider = UnivariatePolynomialZ64.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
                } while (divider.degree == 0 || MachineArithmetic.gcd(divider.lc(), modulus) != 1);

                UnivariatePolynomialZp64 dividendMod = dividend.modulus(modulus), dividerMod = divider.modulus(modulus);
                if (dividendMod.degree == 0) {
                    --i;
                    continue out;
                }
                start = System.nanoTime();
                UnivariatePolynomialZp64[] actualMod = UnivariateDivision.divideAndRemainderLinearDividerModulus(dividendMod, dividerMod, true);
                fast.addValue(System.nanoTime() - start);
                start = System.nanoTime();
                UnivariatePolynomialZp64[] expectedMod = UnivariateDivision.divideAndRemainderClassic0(dividendMod, dividerMod, true);
                gen.addValue(System.nanoTime() - start);
                assertArrayEquals(expectedMod, actualMod);

                UnivariatePolynomialZp64[] actualNoCopyMod = UnivariateDivision.divideAndRemainderLinearDividerModulus(dividendMod.clone(), dividerMod, false);
                UnivariatePolynomialZp64[] expectedNoCopyMod = UnivariateDivision.divideAndRemainderClassic0(dividendMod.clone(), dividerMod, false);
                assertArrayEquals(expectedMod, actualNoCopyMod);
                assertArrayEquals(expectedMod, expectedNoCopyMod);
            }

            do {
                divider = UnivariatePolynomialZ64.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = UnivariateDivision.pseudoDivideAndRemainderLinearDivider(dividend, divider, true);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            long factor = MachineArithmetic.safePow(divider.lc(), dividend.degree - divider.degree + 1);
            expected = UnivariateDivision.divideAndRemainderClassic0(dividend, divider, factor, true);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);

            actualNoCopy = pseudoDivideAndRemainderLinearDivider(dividend.clone(), divider, false);
            expectedNoCopy = divideAndRemainderClassic0(dividend.clone(), divider, factor, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);

            do {
                divider = UnivariatePolynomialZ64.create(rndd.nextLong(-10, 10), rndd.nextLong(-10, 10));
            } while (divider.degree == 0);
            start = System.nanoTime();
            actual = UnivariateDivision.pseudoDivideAndRemainderLinearDividerAdaptive(dividend, divider, true);
            fastPseudo.addValue(System.nanoTime() - start);
            start = System.nanoTime();
            expected = UnivariateDivision.pseudoDivideAndRemainderAdaptive0(dividend, divider, true);
            genPseudo.addValue(System.nanoTime() - start);
            assertArrayEquals(expected, actual);

            actualNoCopy = UnivariateDivision.pseudoDivideAndRemainderLinearDividerAdaptive(dividend.clone(), divider, false);
            expectedNoCopy = UnivariateDivision.pseudoDivideAndRemainderAdaptive0(dividend.clone(), divider, false);
            assertArrayEquals(expected, actualNoCopy);
            assertArrayEquals(expected, expectedNoCopy);
        }
        System.out.println("Fast:    " + fast.getMean());
        System.out.println("General: " + gen.getMean());

        System.out.println("       pseudo ");
        System.out.println("Fast:    " + fastPseudo.getMean());
        System.out.println("General: " + genPseudo.getMean());
    }


    @Test(expected = ArithmeticException.class)
    public void test8() throws Exception {
        UnivariatePolynomialZ64 a = UnivariatePolynomialZ64.create(8, -2 * 8, 8, 8 * 2);
        UnivariatePolynomialZ64 b = UnivariatePolynomialZ64.create(0);
        divideAndRemainder(a, b, false);
    }


    @Test
    public void test9() throws Exception {
        UnivariatePolynomialZ64 a = UnivariatePolynomialZ64.create(8, -2 * 8, 8, 8 * 2);
        UnivariatePolynomialZ64 b = UnivariatePolynomialZ64.create(0);
        UnivariatePolynomialZ64[] zeros = {UnivariatePolynomialZ64.zero(), UnivariatePolynomialZ64.zero()};
        assertArrayEquals(zeros, divideAndRemainder(b, a, true));
        assertArrayEquals(zeros, pseudoDivideAndRemainder(b, a, true));
        assertArrayEquals(zeros, UnivariateDivision.pseudoDivideAndRemainderAdaptive(b, a, true));
        assertArrayEquals(Arrays.stream(zeros).map(x -> x.modulus(13)).toArray(size -> new UnivariatePolynomialZp64[size]),
                UnivariateDivision.divideAndRemainder(b.modulus(13), a.modulus(13), true));
    }

    private static <T extends IUnivariatePolynomial<T>> void assertQuotientRemainder(T dividend, T divider, T[] qr) {
        if (qr == null) return;
        assertEquals(dividend, divider.clone().multiply(qr[0]).add(qr[1]));
    }

    private static <T extends IUnivariatePolynomial<T>> void assertPseudoQuotientRemainder(T dividend, T divider, T[] qr) {
        if (qr == null) return;
        T d = divider.clone().multiply(qr[0]).add(qr[1]);
        T[] factor = divideAndRemainder(d, dividend, true);
        assertNotNull(factor);
        assertTrue(factor[1].isZero());
        assertTrue(factor[0].isConstant());
    }

    private static UnivariatePolynomialZp64 inverseModMonomial0(UnivariatePolynomialZp64 poly, int xDegree) {
        if (xDegree < 1)
            return null;
        if (poly.cc() != 1)
            throw new IllegalArgumentException();
        int r = UnivariateDivision.log2(xDegree);
        UnivariatePolynomialZp64 gPrev = poly.createOne();
        for (int i = 0; i < r; ++i) {
            UnivariatePolynomialZp64 tmp = gPrev.clone().multiply(2).subtract(gPrev.square().multiply(poly));
            gPrev = UnivariateDivision.remainderMonomial(tmp, 1 << i, false);
        }
        return gPrev;
    }

    @Test
    public void test10_InverseModRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        long modulus = getModulusRandom(10);
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomialZp64 f = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(100), modulus, rnd);
            f.data[0] = 1;
            int modDegree = 1 + rnd.nextInt(2 * f.degree);
            UnivariatePolynomialZp64 invMod = inverseModMonomial0(f, modDegree);
            assertInverseModMonomial(f, invMod, modDegree);
        }
    }

    static void assertInverseModMonomial(UnivariatePolynomialZp64 poly, UnivariatePolynomialZp64 invMod, int monomialDegree) {
        assertTrue(polyMultiplyMod(poly, invMod, UnivariatePolynomialZp64.monomial(poly.ring.modulus, 1, monomialDegree), true).isOne());
    }

    /**
     * Returns the remainder of dividing the product {@code (m1 * m2)} by {@code polyModulus}.
     *
     * @param m1          the first multiplier
     * @param m2          the second multiplier
     * @param polyModulus the modulus
     * @param copy        whether to clone {@code m1}; if not, the result will be placed directly to the data structure
     *                    of the first multiplier {@code m1} and the original data of {@code m1} will be lost
     * @return {@code (m1 * m2) % polyModulus}
     */
    public static <T extends IUnivariatePolynomial<T>> T polyMultiplyMod(T m1, T m2, T polyModulus, boolean copy) {
        return remainder((copy ? m1.clone() : m1).multiply(m2), polyModulus, false);
    }

    @Test
    public void test11_InverseModStructureRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            long modulus = getModulusRandom(20);
            UnivariatePolynomialZp64 p = RandomUnivariatePolynomials.randomMonicPoly(2 + rnd.nextInt(100), modulus, rnd);
            p.data[0] = 1;

            UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(p);
            for (int j = 0; j < 30; j++) {
                int xDegree = 1 + rnd.nextInt(1025);
                assertEquals(invMod.getInverse(xDegree), inverseModMonomial0(p.clone().reverse(), xDegree));
            }
        }
    }


    @Test
    public void test12_FastDivisionRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 500); i++) {
            long modulus = getModulusRandom(getRandomData().nextInt(29, 33));
            UnivariatePolynomialZp64 b = RandomUnivariatePolynomials.randomMonicPoly(30, modulus, rnd);
            UnivariatePolynomialZp64 a = RandomUnivariatePolynomials.randomMonicPoly(rnd.nextInt(30), modulus, rnd);

            UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(b);
            UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, invMod, true);
            UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
            assertArrayEquals(fast, plain);
        }
    }

    @Test
    public void test13() throws Exception {
        long modulus = 7;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(5, 1, 4, 6, 4, 3, 5, 5, 3, 4, 2, 2, 5, 2, 5, 6, 1, 1, 2, 5, 1, 0, 0, 6, 6, 5, 5, 1, 0, 1, 4, 1, 1).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(2, 5, 3, 1, 1, 5, 6, 3, 4, 0, 0, 5, 4, 0, 2, 1).modulus(modulus);
        UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(b);
        UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, invMod, true);
        UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
        assertArrayEquals(fast, plain);
    }

    @Test
    public void test14() throws Exception {
        long modulus = 7;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(5, 3, 3, 3, 5, 3, 1, 4, -3, 1, 4, 5, 0, 2, 2, -5, 1).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(0, 4, 6, 1, 2, 4, 0, 0, 6, 5, 2, 3, 1, 4, 0, 1).modulus(modulus);
        UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(b);
        UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, invMod, true);
        UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
        assertArrayEquals(fast, plain);
    }

    @Test
    public void test15() throws Exception {
        long modulus = 17;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(0, 6, 2, 1, 10, 15, 16, 15, 2, 11, 13, 0, 1, 15, 5, 13, 8, 14, 13, 14, 15, 1, 1).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1).modulus(modulus);
        UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(b);
        UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, invMod, true);
        UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
        assertArrayEquals(fast, plain);
    }

    @Test
    public void test16() throws Exception {
        long modulus = 17;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(5, 9, 4, 9, 8, 12, 11, 9, 1, 6, 15, 7, 11, 2, 11, 13, 11, 10, 5, 1).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(11, 15, 9, 5, 11, 5, 14, 9, 1, 0, 16, 12, 11, 5, 15, 10, 15, 2, 14, 3, 1, 16, 16, 12, 13, 1, 12, 11, 1, 15, 1).modulus(modulus);
        UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(b);
        UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, invMod, true);
        UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
        assertArrayEquals(fast, plain);
    }

    @Test
    public void test17() throws Exception {
        long modulus = 7;
        UnivariatePolynomialZp64 f = UnivariatePolynomialZ64.create(0, 2, 3, 4, -5, 1).modulus(modulus).reverse();
        int modDegree = f.degree;
        UnivariatePolynomialZp64 invmod = inverseModMonomial0(f, modDegree);
        UnivariatePolynomialZp64 r = polyMultiplyMod(f, invmod, UnivariatePolynomialZp64.monomial(modulus, 1, modDegree), true);
        assertTrue(r.isOne());
    }

    @Test
    public void test18() throws Exception {
        long modulus = 17;
        UnivariatePolynomialZp64 f = UnivariatePolynomialZ64.create(7, 12, 12, 12, 13, 2, 7, 10, 7, 15, 13, 1, 10, 16, 6, 1).modulus(modulus).reverse();
        int modDegree = 9;
        UnivariatePolynomialZp64 invmod = inverseModMonomial0(f, modDegree);
        UnivariatePolynomialZp64 r = polyMultiplyMod(f, invmod, UnivariatePolynomialZp64.monomial(modulus, 1, modDegree), true);
        assertTrue(r.isOne());
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test19_FastDivisionPerformance() throws Exception {
        long modulus = 5659;
        RandomGenerator rnd = getRandom();
        UnivariatePolynomialZp64 divider = RandomUnivariatePolynomials.randomMonicPoly(118, modulus, rnd);

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        UnivariateDivision.InverseModMonomial invRev = UnivariateDivision.fastDivisionPreConditioning(divider);
        long nIterations = its(1000, 15000);
        for (int i = 0; i < nIterations; i++) {
            if (i * 10 == nIterations) {
                classic.clear();
                fast.clear();
            }
            UnivariatePolynomialZ64 dividendZ = RandomUnivariatePolynomials.randomPoly(3 * divider.degree / 2, (int) modulus, rnd);
            UnivariatePolynomialZp64 dividend = dividendZ.modulus(modulus);

            long start = System.nanoTime();
            UnivariatePolynomialZp64[] qdPlain = UnivariateDivision.divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            UnivariatePolynomialZp64[] qdNewton = UnivariateDivision.divideAndRemainderFast(dividend, divider, invRev, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals(qdPlain, qdNewton);
        }

        System.out.println("==== Plain ====");
        System.out.println(classic.getPercentile(50));

        System.out.println("==== Fast ====");
        System.out.println(fast.getPercentile(50));
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test20_FastDivisionPerformance() throws Exception {
        long modulus = BigPrimes.nextPrime(124987324L);
        RandomGenerator rnd = getRandom();

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        long nIterations = its(1000, 15000);
        int dividerDegree = 156;
        int dividendDegree = 256;
        for (int i = 0; i < nIterations; i++) {
            if (i * 10 == nIterations) {
                classic.clear();
                fast.clear();
            }

            UnivariatePolynomialZp64 divider = RandomUnivariatePolynomials.randomMonicPoly(dividerDegree, modulus, rnd);
            UnivariatePolynomialZ64 dividendZ = RandomUnivariatePolynomials.randomPoly(dividendDegree, (int) modulus, rnd);
            UnivariatePolynomialZp64 dividend = dividendZ.modulus(modulus);
            divider.multiply(3);

            long start = System.nanoTime();
            UnivariatePolynomialZp64[] qdPlain = UnivariateDivision.divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            UnivariatePolynomialZp64[] qdNewton = UnivariateDivision.divideAndRemainderFast(dividend, divider, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals(qdPlain, qdNewton);
        }


        System.out.println("==== Plain ====");
        System.out.println(classic.getMean());

        System.out.println("==== Fast ====");
        System.out.println(fast.getMean());
    }

    @Test
    public void test21() throws Exception {
        long modulus = 17;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(5, 9, 4, 9, 8, 12, 11, 9, 1, 6, 15, 7, 11, 2, 11, 13, 11, 10, 5, 1).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(11, 15, 9, 5, 11, 5, 14, 9, 1, 0, 16, 12, 11, 5, 15, 10, 15, 2, 14, 3, 1, 16, 16, 12, 13, 1, 12, 11, 1, 15, 13).modulus(modulus);
        UnivariatePolynomialZp64[] fast = UnivariateDivision.divideAndRemainderFast(a, b, true);
        UnivariatePolynomialZp64[] plain = UnivariateDivision.divideAndRemainderClassic(a, b, true);
        assertArrayEquals(fast, plain);
    }

    @Test
    public void test22() throws Exception {
        IntegersZp domain = new IntegersZp(7);
        UnivariatePolynomial<BigInteger> bDividend = UnivariatePolynomial.create(domain, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1);
        UnivariatePolynomialZp64 lDividend = asOverZp64(bDividend);


        UnivariatePolynomial<BigInteger> bDivider = UnivariatePolynomial.create(domain, 1, 2, 3, 3, 2, 1);
        UnivariatePolynomialZp64 lDivider = asOverZp64(bDivider);

        UnivariatePolynomial<BigInteger>[] bqd = UnivariateDivision.divideAndRemainderFast(bDividend, bDivider, true);
        UnivariatePolynomialZp64[] lqd = UnivariateDivision.divideAndRemainderFast(lDividend, lDivider, true);
        assertArrayEquals(new UnivariatePolynomialZp64[]{asOverZp64(bqd[0]), asOverZp64(bqd[1])}, lqd);
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test19_BigInteger_FastDivisionPerformance() throws Exception {
        BigInteger modulus = new BigInteger("1247842098624308285367648396697");//BigPrimes.nextPrime(new BigInteger(100, rnd));
        RandomGenerator rnd = getRandom();
        UnivariatePolynomial<BigInteger> divider = RandomUnivariatePolynomials.randomMonicPoly(128, modulus, rnd);

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        UnivariateDivision.InverseModMonomial<UnivariatePolynomial<BigInteger>> invRev = UnivariateDivision.fastDivisionPreConditioning(divider);
        long nIterations = its(1000, 5000);
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                classic.clear();
                fast.clear();
            }
            UnivariatePolynomial<BigInteger> dividendZ = RandomUnivariatePolynomials.randomPoly(3 * divider.degree / 2, modulus, rnd);
            UnivariatePolynomial<BigInteger> dividend = dividendZ.setRing(new IntegersZp(modulus));

            long start = System.nanoTime();
            UnivariatePolynomial<BigInteger>[] qdPlain = UnivariateDivision.divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            UnivariatePolynomial<BigInteger>[] qdNewton = UnivariateDivision.divideAndRemainderFast(dividend, divider, invRev, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals(qdPlain, qdNewton);
        }

        System.out.println("==== Plain ====");
        System.out.println(classic.getPercentile(50));

        System.out.println("==== Fast ====");
        System.out.println(fast.getPercentile(50));
    }

    @Test
    @Benchmark(runAnyway = true)
    public void test20_BigInteger_FastDivisionPerformance() throws Exception {
        RandomGenerator rnd = getRandom();
        BigInteger modulus = new BigInteger("1247842098624308285367648396697");//BigPrimes.nextPrime(new BigInteger(100, rnd));

        DescriptiveStatistics classic = new DescriptiveStatistics(), fast = new DescriptiveStatistics();
        long nIterations = its(1000, 5000);
        int dividerDegree = 126;
        int dividendDegree = 256;
        for (int i = 0; i < nIterations; i++) {
            if (nIterations / 10 == i) {
                classic.clear();
                fast.clear();
            }

            UnivariatePolynomial<BigInteger> divider = RandomUnivariatePolynomials.randomMonicPoly(dividerDegree, modulus, rnd);
            UnivariatePolynomial<BigInteger> dividendZ = RandomUnivariatePolynomials.randomPoly(dividendDegree, modulus, rnd);
            UnivariatePolynomial<BigInteger> dividend = dividendZ.setRing(new IntegersZp(modulus));
            divider.multiply(BigInteger.THREE);

            long start = System.nanoTime();
            UnivariatePolynomial<BigInteger>[] qdPlain = UnivariateDivision.divideAndRemainderClassic(dividend, divider, true);
            long plain = System.nanoTime() - start;
            classic.addValue(plain);

            start = System.nanoTime();
            UnivariatePolynomial<BigInteger>[] qdNewton = UnivariateDivision.divideAndRemainderFast(dividend, divider, true);
            long newton = System.nanoTime() - start;
            fast.addValue(newton);

            assertArrayEquals("dividend = " + dividend + ";\ndivider = " + divider + ";\n", qdPlain, qdNewton);
        }


        System.out.println("==== Plain ====");
        System.out.println(classic.getMean());

        System.out.println("==== Fast ====");
        System.out.println(fast.getMean());
    }

//    @Test
//    public void test23() throws Exception {
//        BigInteger modulus = new BigInteger("998427238390739620139");
//        bMutablePolynomialZp dividend = bMutablePolynomialZ.parse("989441076315244786644+174683251098354358x^1+2939699558711223765x^2+993164729241539182424x^3+8652504087827847685x^4+2978039521215483585x^5+5687372540827878771x^6+3684693598277313443x^7+3034113231916032517x^8+1842720927561159970x^9+1401489172494884190x^10").modulus(modulus);
//        bMutablePolynomialZp divider = bMutablePolynomialZ.parse("718119058879299323824+59748620370951943044x^1+27715597040703811206x^2+3x^3").modulus(modulus);
//
//        bMutablePolynomialZp[] classic = UnivariateDivision.divideAndRemainderClassic(dividend, divider, true);
//        bMutablePolynomialZp[] fast = UnivariateDivision.divideAndRemainderFast(dividend, divider, true);
//
//        System.out.println(Arrays.toString(classic));
//        System.out.println(Arrays.toString(fast));
//        assertQuotientRemainder(dividend, divider, classic);
//    }

    @Test(expected = ArithmeticException.class)
    public void testDivideByZero() throws Exception {
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4).modulus(3);
        System.out.println(Arrays.toString(divideAndRemainder(poly, poly.createZero(), false)));
    }

    @Test
    public void test24() throws Exception {
        UnivariatePolynomialZ64
                a = UnivariatePolynomialZ64.create(1, 2, 3),
                b = UnivariatePolynomialZ64.create(3, 4, 5, 3),
                c = UnivariatePolynomialZ64.create(5, 6, 7);
        a = a.multiply(b).add(c);
        UnivariatePolynomialZ64 r = remainder(a, b, true);
        assertNotNull(r);
        assertEquals(divideAndRemainder(a, b, true)[1], r);
    }

    @Test
    public void test25() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = 10000;
        for (int i = 0; i < nIterations; i++) {
            UnivariatePolynomialZ64
                    a = randomPoly(rndd.nextInt(1, 5), rnd),
                    b = randomPoly(rndd.nextInt(a.degree + 1, a.degree + 1 + 5), rnd),
                    rem = randomPoly(rndd.nextInt(0, a.degree), rnd);
            a = a.multiply(b).add(rem);
            UnivariatePolynomialZ64 r = remainder(a, b, true);
            assertNotNull(r);
            assertEquals(rem, r);

            long modulus = getModulusRandom(10);
            UnivariatePolynomialZp64 rp = remainder(a.modulus(modulus), b.modulus(modulus), true);
            assertNotNull(rp);
            assertEquals(rem.modulus(modulus), rp);

            UnivariatePolynomial<BigInteger> br = remainder(a.toBigPoly(), b.toBigPoly(), true);
            assertNotNull(br);
            assertEquals(rem.toBigPoly(), br);

            IntegersZp mDomain = new IntegersZp(modulus);
            UnivariatePolynomial<BigInteger> brp = remainder(a.toBigPoly().setRing(mDomain), b.toBigPoly().setRing(mDomain), true);
            assertNotNull(brp);
            assertEquals(rem.toBigPoly().setRing(mDomain), brp);
        }
    }

    @Test
    public void test26() throws Exception {
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(Rings.Z,
                new BigInteger("50924915532280396921059912058744139974274010464717486445566429090906004335144553740099697325823981482217595251233173036221148816419907025031109037285812602964126981398549972404972445712495606246251612602646515107556469916968928792288517114629524761056619784125538889985573728321588316873829401388110009323502747325559537492401603348781184777517075447383315534987333005141430297794452769344710334113504161517191471620923923009168593251562951163703046204657823731073238889133009776679947959284882493347742931374480069087134327139957626655815839165341149468676605096098488269737054120888284410186940922588811865217977201434671980392818129023662348422616219595408007164847429098690018036788941569325028261778694448305169116366300968797299392607269512834399626743351451294227407026614474070168043362619155452474803802317594566620282821580866690460901450722523940253763809526916168510380458823716227207456429213838428363313858223438249583518220909107257207453781168662470638622229018698708892857265124137406565985565106519169754356798964383232"),
                new BigInteger("37940726704358973981646894376520845717379643722668546596604774191462649371799436327770509765417751410853020097236259181340723780640646980528823132305918233917348431603978910702974168356786878991952039076697435685153720389954147196105732666546398882808075467883408078172773943074952881821883968276970586716085556521926276724965832984764113660527315402786928901376759667964550095121161582001082154197929523015687629358088978002050194666997140380633120191985503069336469087766567723745827536389920145897766163207961925615721265621132693244136269578827145639329792340827479965479742934874844926129035539562075005079119861286286379018264676225498300503424412683235789124363971763488749051910815230910592202224639481529849460050198630242102474952014158371176492350931799358619802630663775654806908014463551827044421713475118514891517003720144816583416199973842512734641894953396023836643821945946850868732474395966931084754653075573569398089513199914298230954227062809869737803168417328622780936780114654437560463636354964403162931061261855232"));
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(Rings.Z,
                new BigInteger("-425784855629707210539690417657285596009748771201595246777393523226748147711074288140558389550043516613530500261798739170010747708483239470487926534197999795185337107301428405197650111401424401453281285342023506707519100774028143103154883688933695208217191334712115038212326510433566424650523293261999555565022871454564811453312766524934990529599572227820236569578196841210873189542906771110582428737577172320473113774712560826543606283857910494481715951013601677939588946090257941915294685671582476036525511342773894137847945807902831920328130096765066331366254244760254056538237885678798173308462659097285980035599674879837185559523111215792005945149771714839595770409009983341486144352570200746450022078747332821616059006090021772872988207037534463289442687034943036344354336745002436373485142219064592871177691253579923036133374333408048223917924387929466882231114740060638812996586535143719347485716399431436091507892979143097224033509341961137283167473128691774034159819704359664495382150875212794275998540147497752601895426653284077738548432459176682797681127140157440"));
        assertPseudoQuotientRemainder(dividend, divider, pseudoDivideAndRemainder(dividend, divider, true));
    }

    @Test
    public void test27_finiteFields() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int
                minDegree = 5,
                maxDegree = 15;
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomial<UnivariatePolynomialZp64>
                    dividend = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minDegree, maxDegree), FiniteField.GF17p5, rnd),
                    divider = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(1, dividend.degree + 1), FiniteField.GF17p5, rnd);

            assertQuotientRemainder(dividend, divider, divideAndRemainder(dividend, divider, true));
        }
    }
}