package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.AlgebraicNumberField;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariateResultants.*;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.RandomDataGenerator;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.univar.RandomUnivariatePolynomials.randomPoly;
import static cc.redberry.rings.poly.univar.UnivariateDivision.*;
import static cc.redberry.rings.poly.univar.UnivariateGCD.*;
import static cc.redberry.rings.poly.univar.UnivariateResultants.*;
import static cc.redberry.rings.util.TimeUnits.nanosecondsToString;
import static org.junit.Assert.*;

/**
 * Created by poslavsky on 15/02/2017.
 */
public class UnivariateGCDTest extends AUnivariateTest {

    @SuppressWarnings("ConstantConditions")
    @Test
    public void test1() throws Exception {
        long modulus = 11;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(3480, 8088, 8742, 13810, 12402, 10418, 8966, 4450, 950).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(2204, 2698, 3694, 3518, 5034, 5214, 5462, 4290, 1216).modulus(modulus);

        APolynomialRemainderSequence<UnivariatePolynomialZp64> prs = ClassicalPRS(a, b);
        UnivariatePolynomialZp64 gcd = prs.gcd();
        assertEquals(3, gcd.degree);
        assertTrue(divideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(divideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test
    public void test2() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(0, 14, 50, 93, 108, 130, 70);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(63, 92, 143, 245, 146, 120, 90);
        UnivariatePolynomial<BigInteger> gcd = PrimitivePRS(dividend, divider).gcd();
        UnivariatePolynomial<BigInteger> expected = UnivariatePolynomial.create(7, 4, 10, 10).negate();
        assertEquals(expected, gcd);
    }

    @Test
    public void test3() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(7, -7, 0, 1);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(-7, 0, 3);
        APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> naive = PseudoPRS(dividend.clone(), divider.clone());
        List<UnivariatePolynomial<BigInteger>> expectedNaive = new ArrayList<>();
        expectedNaive.add(dividend);
        expectedNaive.add(divider);
        expectedNaive.add(UnivariatePolynomial.create(63, -42));
        expectedNaive.add(UnivariatePolynomial.create(-441L));
        assertEquals(expectedNaive, naive.remainders);

        APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> primitive = PrimitivePRS(dividend.clone(), divider.clone());
        List<UnivariatePolynomial<BigInteger>> expectedPrimitive = new ArrayList<>();
        expectedPrimitive.add(dividend);
        expectedPrimitive.add(divider);
        expectedPrimitive.add(UnivariatePolynomial.create(3, -2));
        expectedPrimitive.add(UnivariatePolynomial.create(1L));
        assertEquals(expectedPrimitive, primitive.remainders);

        APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> subresultant = SubresultantPRS(dividend.clone(), divider.clone());
        List<UnivariatePolynomial<BigInteger>> expectedSubresultant = new ArrayList<>();
        expectedSubresultant.add(dividend);
        expectedSubresultant.add(divider);
        expectedSubresultant.add(UnivariatePolynomial.create(63, -42));
        expectedSubresultant.add(UnivariatePolynomial.create(-49L));
        assertEquals(expectedSubresultant, subresultant.remainders);
    }

    @Test
    public void test4() throws Exception {
        for (long sign = -1; sign <= 2; sign += 2) {
            UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2).multiply(sign);
            UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(6, 7, 9, 8, 3).multiply(sign);
            APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> subresultant = SubresultantPRS(dividend, divider);
            List<UnivariatePolynomial<BigInteger>> expectedSubresultant = new ArrayList<>();
            expectedSubresultant.add(dividend);
            expectedSubresultant.add(divider);
            expectedSubresultant.add(UnivariatePolynomial.create(159, 112, 192, 164).multiply(sign));
            expectedSubresultant.add(UnivariatePolynomial.create(4928, 3068, 5072).multiply(sign));
            expectedSubresultant.add(UnivariatePolynomial.create(65840, -98972).multiply(sign));
            expectedSubresultant.add(UnivariatePolynomial.create(3508263).multiply(sign));
            assertEquals(expectedSubresultant, subresultant.remainders);
        }
    }

    @Test
    public void test4a() throws Exception {
        for (long sign1 = -1; sign1 <= 2; sign1 += 2)
            for (long sign2 = -1; sign2 <= 2; sign2 += 2) {
                UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(1, 1, 6, -3, 5).multiply(sign1);
                UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(-9, -3, -10, 2, 8).multiply(sign2);
                APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> subresultant = SubresultantPRS(dividend, divider);
                List<UnivariatePolynomial<BigInteger>> expectedSubresultant = new ArrayList<>();
                expectedSubresultant.add(dividend);
                expectedSubresultant.add(divider);
                expectedSubresultant.add(UnivariatePolynomial.create(-53, -23, -98, 34).multiply(sign1 * sign2));
                expectedSubresultant.add(UnivariatePolynomial.create(4344, 3818, 9774));
                expectedSubresultant.add(UnivariatePolynomial.create(-292677, 442825).multiply(sign1 * sign2));
                expectedSubresultant.add(UnivariatePolynomial.create(22860646));
                assertEquals(expectedSubresultant, subresultant.remainders);
            }
    }

    @Test
    public void test5() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2);
        APolynomialRemainderSequence prs = PseudoPRS(dividend.clone(), dividend.clone());

        assertEquals(2, prs.remainders.size());
        assertEquals(dividend, prs.gcd());
    }

    @Test
    public void test6() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2).multiply(2);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2);
        APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs;

        for (boolean prim : new boolean[]{true, false}) {
            prs = prim ? PrimitivePRS(dividend, divider) : PseudoPRS(dividend, divider);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());

            prs = prim ? PrimitivePRS(divider, dividend) : PseudoPRS(divider, dividend);
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        prs = SubresultantPRS(dividend.clone(), divider.clone());
        assertEquals(2, prs.remainders.size());
        assertEquals(divider, prs.gcd());
    }

    @Test
    public void test7() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2).multiply(2);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(7, 4, 10, 10, 0, 2);

        for (APolynomialRemainderSequence prs : runAlgorithms(dividend, divider)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }

        for (APolynomialRemainderSequence prs : runAlgorithms(divider, dividend)) {
            assertEquals(2, prs.remainders.size());
            assertEquals(divider, prs.gcd());
        }
    }

    @Test
    public void test8() throws Exception {
        //test long overflow
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(1, 1, 1, 1).multiply(2);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(1, 0, 2);

        for (APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider))
            assertEquals(0, prs.gcd().degree);
    }

    @Test
    public void test9() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(-58, 81, -29, -77, 81, 42, -38, 48, 94, 6, 55);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(1, 2, 1);
        UnivariatePolynomial<BigInteger> expected = UnivariatePolynomial.create(-58, -35, 75, -54, -102, 127, 127, 14, 152, 242, 161, 116, 55);
        assertEquals(expected, a.multiply(b));
    }

    @Test
    public void testRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomial<BigInteger> dividend = RandomUnivariatePolynomials.randomPoly(5, rnd).toBigPoly();
            UnivariatePolynomial<BigInteger> divider = RandomUnivariatePolynomials.randomPoly(0, rnd).toBigPoly();
            for (APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider)) {
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
            UnivariatePolynomial<BigInteger> dividend = RandomUnivariatePolynomials.randomPoly(5, 5, rnd).toBigPoly();
            UnivariatePolynomial<BigInteger> divider = RandomUnivariatePolynomials.randomPoly(5, 5, rnd).toBigPoly();
            try {
                for (APolynomialRemainderSequence prs :
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
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(4, -1, -4, -3, 2, 4);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(-4, 1, 4, 3, -2, -4);
        assertPolynomialRemainders(dividend, divider, PrimitivePRS(dividend, divider));
    }

    @Test
    public void testRandom2b() throws Exception {
        UnivariatePolynomial<BigInteger> dividend = UnivariatePolynomial.create(0, 0, 4, 0, 4, 4);
        UnivariatePolynomial<BigInteger> divider = UnivariatePolynomial.create(0, 0, 4, 0, 4, 4);
        assertPolynomialRemainders(dividend, divider, PseudoPRS(dividend, divider));
    }

    @Test
    public void testRandom3() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(500, 3000); i++) {
            UnivariatePolynomialZ64 dividend = RandomUnivariatePolynomials.randomPoly(10, 500, rnd);
            UnivariatePolynomialZ64 divider = RandomUnivariatePolynomials.randomPoly(10, 500, rnd);
            for (long prime : getModulusArray(9, 1, 40)) {
                if (dividend.lc() % prime == 0 || divider.lc() % prime == 0)
                    continue;
                UnivariatePolynomialZp64 a = dividend.modulus(prime);
                UnivariatePolynomialZp64 b = divider.modulus(prime);
                APolynomialRemainderSequence<UnivariatePolynomialZp64> euclid = ClassicalPRS(a, b);
                assertPolynomialRemainders(a, b, euclid);
            }
        }
    }

    @Test
    public void test10() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(1740, 4044, 4371, 6905, 6201, 5209, 4483, 2225, 475);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(1102, 1349, 1847, 1759, 2517, 2607, 2731, 2145, 608);
        assertEquals(UnivariatePolynomial.create(29, 21, 32, 19), ModularGCD(a, b));
    }

    @Test
    public void test11() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(0, 1, 1, -6, 17, -18, 14);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(0, 4, -3, 0, 8, 0, 4);
        assertEquals(UnivariatePolynomial.create(0, 1, -2, 2), ModularGCD(a, b));
    }

    @Test
    public void test12() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(1, 2, 3, 3, 6, 9);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(1, 3, 6, 5, 3);
        assertEquals(UnivariatePolynomial.create(1, 2, 3), ModularGCD(a, b));

        a = UnivariatePolynomial.create(0, 0, 1, 2);
        b = UnivariatePolynomial.create(-1, 0, 4);
        assertEquals(UnivariatePolynomial.create(1, 2), ModularGCD(a, b));
    }

    @Test
    public void test13() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(-1, 0, 4);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(-1, 0, 0, 0, 16);
        UnivariatePolynomial<BigInteger> gcd = ModularGCD(a, b);
        assertEquals(UnivariatePolynomial.create(-1, 0, 4), gcd);
    }

    @Test
    public void test14() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(-1, 0, 4);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(1, -5, 6);
        UnivariatePolynomial<BigInteger> gcd = ModularGCD(a, b);
        assertEquals(UnivariatePolynomial.create(-1, 2), gcd);
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
                    UnivariatePolynomialZ64 a = RandomUnivariatePolynomials.randomPoly(1 + rnd.nextInt(7), 100, rnd);
                    UnivariatePolynomialZ64 b = RandomUnivariatePolynomials.randomPoly(1 + rnd.nextInt(6), 100, rnd);
                    UnivariatePolynomialZ64 gcd = RandomUnivariatePolynomials.randomPoly(2 + rnd.nextInt(5), 30, rnd);
                    a = a.multiply(gcd);
                    b = b.multiply(gcd);

                    long start = System.nanoTime();
                    UnivariatePolynomialZ64 mgcd = ModularGCD(a, b);
                    timings.addValue(System.nanoTime() - start);

                    assertFalse(mgcd.isConstant());

                    UnivariatePolynomialZ64[] qr = UnivariateDivision.pseudoDivideAndRemainderAdaptive(mgcd, gcd, true);
                    assertNotNull(qr);
                    assertTrue(qr[1].isZero());
                    UnivariatePolynomialZ64[] qr1 = UnivariateDivision.pseudoDivideAndRemainderAdaptive(mgcd.clone(), gcd, false);
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
        UnivariatePolynomialZ64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        assertEquals(0, poly.evaluate(1));
        assertEquals(-3999, poly.evaluate(2));
        assertEquals(1881, poly.evaluate(-2));
        assertEquals(566220912246L, poly.evaluate(-17));
    }

    @Test
    public void test16() throws Exception {
        UnivariatePolynomialZ64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
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
        UnivariatePolynomialZ64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, -1, -2, -3, -4, -5);
        poly.evaluateAtRational(1, 5);
    }

    @Test
    public void test18() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(8, -2 * 8, 8, 8 * 2);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.zero(Rings.Z);
        assertEquals(a, SubresultantPRS(a, b).gcd());
        assertEquals(a, SubresultantPRS(b, a).gcd());
        assertEquals(a, PrimitivePRS(a, b).gcd());
        assertEquals(a, PrimitivePRS(b, a).gcd());
        assertEquals(a, EuclidGCD(a, b));
        assertEquals(a, EuclidGCD(b, a));
        assertEquals(a, PolynomialGCD(a, b));
        assertEquals(a, PolynomialGCD(b, a));
    }

    @Test
    public void test19() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(8, -2 * 8, 8, 8 * 2);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(2);
        assertEquals(b, SubresultantPRS(a, b).gcd());
        assertEquals(b, SubresultantPRS(b, a).gcd());
        assertEquals(b, PrimitivePRS(a, b).gcd());
        assertEquals(b, PrimitivePRS(b, a).gcd());
        assertEquals(b, PseudoPRS(a, b).gcd());
        assertEquals(b, PseudoPRS(b, a).gcd());
        assertEquals(b, EuclidGCD(a, b));
        assertEquals(b, EuclidGCD(b, a));
        assertEquals(b, PolynomialGCD(a, b));
        assertEquals(b, PolynomialGCD(b, a));
    }

    @Test
    public void test20() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(32);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(24);
        UnivariatePolynomial<BigInteger> gcd = UnivariatePolynomial.create(8);
        assertEquals(gcd, PolynomialGCD(a, b));
        assertEquals(gcd, SubresultantPRS(a, b).gcd());
        assertEquals(gcd, PrimitivePRS(a, b).gcd());
        assertEquals(gcd, PseudoPRS(a, b).gcd());
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
    private static <T extends IUnivariatePolynomial<T>> void assertPolynomialRemainders(T a, T b, APolynomialRemainderSequence<T> prs) {
        if (a.degree() < b.degree()) {
            assertPolynomialRemainders(b, a, prs);
            return;
        }

        assertEquals(a, prs.remainders.get(0));
        assertEquals(b, prs.remainders.get(1));

        T gcd = prs.gcd().clone();
        if (!a.isOverField())
            gcd = gcd.primitivePart();
        assertTrue(UnivariateDivision.pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(UnivariateDivision.pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }
//
//    @SuppressWarnings("ConstantConditions")
//    private static <T extends IMutablePolynomialZp<T>> void assertPolynomialRemainders(T a, T b, APolynomialRemainderSequence<T> prs) {
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

    private static List<PolynomialRemainderSequence<BigInteger>> runAlgorithms(UnivariatePolynomial<BigInteger> dividend,
                                                                               UnivariatePolynomial<BigInteger> divider,
                                                                               GCDAlgorithm... algorithms) {
        ArrayList<PolynomialRemainderSequence<BigInteger>> r = new ArrayList<>();
        for (GCDAlgorithm algorithm : algorithms)
            r.add(algorithm.gcd(dividend, divider));
        return r;
    }

    private enum GCDAlgorithm {
        PolynomialEuclid,
        PolynomialPrimitiveEuclid,
        SubresultantEuclid;

        PolynomialRemainderSequence<BigInteger> gcd(UnivariatePolynomial<BigInteger> dividend, UnivariatePolynomial<BigInteger> divider) {
            switch (this) {
                case PolynomialEuclid:
                    return PseudoPRS(dividend, divider);
                case PolynomialPrimitiveEuclid:
                    return PrimitivePRS(dividend, divider);
                case SubresultantEuclid:
                    return SubresultantPRS(dividend, divider);
            }
            throw new IllegalArgumentException();
        }
    }

    @Test
    public void test21() throws Exception {
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(8, 2, -1, -2, -7);
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(1, -9, -5, -21);
        assertTrue(ModularGCD(a, b).isOne());
    }

    @Test
    public void test22() throws Exception {
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(1, 2, 3, 4, 3, 2, 1).modulus(25);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(1, 2, 3, 1, 3, 2, 1).modulus(25);
        assertExtendedEuclidGCD(a, b);
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

        UnivariatePolynomial<BigInteger> gcdActual = PseudoPRS(a, b).gcd();
        assertGCD(a, b, gcdActual);

        APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs = SubresultantPRS(a, b);
        assertPolynomialRemainders(a, b, prs);

        UnivariatePolynomial<BigInteger> gcdSubresultant = prs.gcd();
        assertEquals(gcdActual.degree, gcdSubresultant.degree);

        UnivariatePolynomial<BigInteger> gcdModular = ModularGCD(a, b);
        assertEquals(gcdActual.degree, gcdModular.degree);

        System.out.println(gcdActual.normMax().bitLength());
    }

    @Test
    public void testRandom1_bigPoly() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomial<BigInteger> dividend = RandomUnivariatePolynomials.randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
            UnivariatePolynomial<BigInteger> divider = RandomUnivariatePolynomials.randomPoly(0, BigInteger.LONG_MAX_VALUE, rnd);
            for (APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider)) {
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
            for (APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> prs : runAlgorithms(dividend, divider, GCDAlgorithm.SubresultantEuclid, GCDAlgorithm.PolynomialPrimitiveEuclid)) {
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
                IntegersZp domain = new IntegersZp(prime);
                UnivariatePolynomial<BigInteger> a = dividend.setRing(domain);
                UnivariatePolynomial<BigInteger> b = divider.setRing(domain);
                APolynomialRemainderSequence<UnivariatePolynomial<BigInteger>> euclid = ClassicalPRS(a, b);
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
                    UnivariatePolynomial<BigInteger> a = RandomUnivariatePolynomials.randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    UnivariatePolynomial<BigInteger> b = RandomUnivariatePolynomials.randomPoly(1 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
                    UnivariatePolynomial<BigInteger> gcd = RandomUnivariatePolynomials.randomPoly(2 + rnd.nextInt(30), BigInteger.LONG_MAX_VALUE, rnd);
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

        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(Rings.Z, new BigInteger("-4914"), new BigInteger("6213"), new BigInteger("3791"), new BigInteger("996"), new BigInteger("-13304"), new BigInteger("-1567"), new BigInteger("2627"), new BigInteger("15845"), new BigInteger("-12626"), new BigInteger("-6383"), new BigInteger("294"), new BigInteger("26501"), new BigInteger("-17063"), new BigInteger("-14635"), new BigInteger("9387"), new BigInteger("-7141"), new BigInteger("-8185"), new BigInteger("17856"), new BigInteger("4431"), new BigInteger("-13075"), new BigInteger("-7050"), new BigInteger("14672"), new BigInteger("3690"), new BigInteger("-3990"));
        UnivariatePolynomial<BigInteger> a = UnivariatePolynomial.create(Rings.Z, new BigInteger("419563715"), new BigInteger("419566193"), new BigInteger("3612"), new BigInteger("3444"), new BigInteger("419563127"), new BigInteger("419564681"), new BigInteger("419565017"), new BigInteger("419564387"), new BigInteger("419563463"), new BigInteger("3192"), new BigInteger("419563841"), new BigInteger("419563001"));
        UnivariatePolynomial<BigInteger> b = UnivariatePolynomial.create(Rings.Z, new BigInteger("209783497"), new BigInteger("9989688"), new BigInteger("379608231"), new BigInteger("399587609"), new BigInteger("59938143"), new BigInteger("29969072"), new BigInteger("99896901"), new BigInteger("359628849"), new BigInteger("329659781"), new BigInteger("239752567"), new BigInteger("19979379"), new BigInteger("179814423"), new BigInteger("1"));

        IntegersZp domain = new IntegersZp(modulus);
        UnivariatePolynomial<BigInteger> aMod = a.setRing(domain).monic(poly.lc());
        UnivariatePolynomial<BigInteger> bMod = b.setRing(domain).monic();
        UnivariatePolynomial<BigInteger>[] xgcd = UnivariateGCD.ExtendedEuclidGCD(aMod, bMod);
        UnivariatePolynomialZp64[] lxgcd = UnivariateGCD.ExtendedEuclidGCD(UnivariatePolynomial.asOverZp64(aMod), UnivariatePolynomial.asOverZp64(bMod));
        assertEquals(UnivariatePolynomial.asOverZp64(xgcd[0]), lxgcd[0]);
        assertEquals(UnivariatePolynomial.asOverZp64(xgcd[1]), lxgcd[1]);
        assertEquals(UnivariatePolynomial.asOverZp64(xgcd[2]), lxgcd[2]);
    }

    @Test
    public void test25() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.create(Rings.Q,
                Rings.Q.parse("2/3"),
                Rings.Q.parse("4/5"),
                Rings.Q.parse("1/2"),
                Rings.Q.parse("-31/2"));
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.create(Rings.Q,
                Rings.Q.parse("7/3"),
                Rings.Q.parse("4/7"),
                Rings.Q.parse("3/2"),
                Rings.Q.parse("-31/12"));
        UnivariatePolynomial<Rational<BigInteger>> gcd = UnivariatePolynomial.create(Rings.Q,
                Rings.Q.parse("4/3"),
                Rings.Q.parse("-4/7"),
                Rings.Q.parse("-1/2"),
                Rings.Q.parse("-1/12"));
        a = a.clone().multiply(gcd);
        b = b.clone().multiply(gcd);
        assertGCD(a, b, PolynomialGCD(a, b));
    }

    @Test
    public void test26() throws Exception {
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(1, 2, 3, 1, 2, 3, 4, 3, 2, 1).modulus(29);//.square
        // ().square().square().square().square().square().square().square().square().square();
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(1, 2, 3, 1, 1, 2, 3, 3).modulus(29);
        assertArrayEquals(ExtendedEuclidGCD(a, b), ExtendedHalfGCD(a, b));
    }

    @Test
    public void test27() throws Exception {
        // assert java heap
        UnivariatePolynomialZp64 a = RandomUnivariatePolynomials.randomMonicPoly(15_000, 19, getRandom());
        UnivariatePolynomialZp64 b = RandomUnivariatePolynomials.randomMonicPoly(15_000, 19, getRandom());
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
            UnivariatePolynomialZp64 a = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree), modulus, rnd).multiply(1 + rnd.nextLong());
            UnivariatePolynomialZp64 b = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree), modulus, rnd).multiply(1 + rnd.nextLong());
            try {
                long start;

                start = System.nanoTime();
                UnivariatePolynomialZp64 expected = EuclidGCD(a, b).monic();
                euclid.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                UnivariatePolynomialZp64 actual = HalfGCD(a, b).monic();
                half.addValue(System.nanoTime() - start);

                assertEquals(expected, actual);

                if (!expected.isConstant())
                    ++nonTrivGCD;
            } catch (Throwable tr) {
                System.out.println("UnivariatePolynomial<BigInteger>." + a.toStringForCopy() + ".modulus(" + modulus + ");");
                System.out.println("UnivariatePolynomial<BigInteger>." + b.toStringForCopy() + ".modulus(" + modulus + ");");
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

    @Test
    public void test33() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.parse("1 + 23123*x^7 + 2344*x^15", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.parse("1 + 23*x - 23454*x^4", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>>[] xgcd = ModularExtendedRationalGCD(a, b);
        assertExtendedGCD(xgcd, a, b);
    }

    @Test
    public void test33a() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.parse("1 + 23123*x^7 + 2344*x^15", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.parse("1 + 23*x - 23454*x^4", Rings.Q);

        a.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(32)));
        b.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(12332)));

        assertExtendedGCD(ModularExtendedRationalGCD(a, b), a, b);
        assertExtendedGCD(ModularExtendedRationalGCD(b, a), b, a);
    }

    @Test
    public void test34() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.parse("1 + 23123*x^7 + 2344*x^15", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.parse("1 + 23*x - 23454*x^4", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> g = UnivariatePolynomial.parse("1 + (23/2)*x - 23454*x^3", Rings.Q);

        a.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(32)));
        b.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(12332)));

        a.multiply(g);
        b.multiply(g);

        assertExtendedGCD(ModularExtendedRationalGCD(a, b), a, b);
        assertExtendedGCD(ModularExtendedRationalGCD(b, a), b, a);
    }

    @Test
    @Benchmark
    public void test35_performance() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.parse("1 + 23123*x^7 + 2344*x^15", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.parse("1 + 23*x + 23454*x^4", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> g = UnivariatePolynomial.parse("1 + (23/2)*x + 23454*x^3", Rings.Q);

        a.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(32)));
        b.multiply(new Rational<>(Rings.Z, BigInteger.valueOf(123), BigInteger.valueOf(12332)));

        a.multiply(g);
        b.multiply(g);

        System.out.println(a);
        System.out.println(b);
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            assertExtendedGCD(ModularExtendedRationalGCD(a, b), a, b);
            assertExtendedGCD(ModularExtendedRationalGCD(b, a), b, a);
            System.out.println(nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @Benchmark
    public void test36_performance() throws Exception {
        UnivariatePolynomial<Rational<BigInteger>> a = UnivariatePolynomial.parse("(296/15) + (874/9)*x + (2083/20)*x^2 + ((-11819)/90)*x^3 + ((-147)/8)*  x^4 + (152461/360)*x^5 + (223567/1440)*x^6 + (22223/432)*  x^7 + ((-583021)/2880)*x^8 + (45407/240)*x^9 + (235373/1260)*  x^10 + ((-58349)/378)*x^11 + (269417/2520)*x^12 + (2402/45)*  x^13 + (206113/420)*x^14 + ((-218167)/1890)*x^15 + ((-62221)/5040)*  x^16 + ((-59279)/2520)*x^17 + (164803/630)*x^18 + (1027/54)*  x^19 + ((-539)/30)*x^20 + ((-97)/3)*x^21 + (64/3)*x^22", Rings.Q);
        UnivariatePolynomial<Rational<BigInteger>> b = UnivariatePolynomial.parse("(388/15) + 221*x + (76253/120)*x^2 + (73661/120)*x^3 + ((-21007)/240)* x^4 + (58939/720)*x^5 + (3215/8)*x^6 + (2599/6)*x^7 + (29683/105)* x^8 + ((-7141)/105)*x^9 + (16021/84)*x^10 + (8807/240)* x^11 + (20747/168)*x^12 + ((-1597627)/10080)*x^13 + (1846219/3360)* x^14 + (334471/6720)*x^15 + ((-644489)/6720)*x^16 + ((-551)/20)* x^17 + (17611/120)*x^18 + (3127/30)*x^19 + ((-4591)/120)* x^20 + (229/30)*x^21 + ((-34)/3)*x^22 + (26/3)*x^23", Rings.Q);

        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            UnivariatePolynomial<Rational<BigInteger>>[] r = ModularExtendedRationalGCD(a, b);
            System.out.println(nanosecondsToString(System.nanoTime() - start));
            assertExtendedGCD(r, a, b);
        }
    }

    @Test
    public void test36_modularExtendedGCDRandom() throws Exception {
        DescriptiveStatistics
                ratRec = new DescriptiveStatistics(),
                resultant = new DescriptiveStatistics();
        int nIterations = its(100, 300);
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        BigInteger bound = BigInteger.valueOf(100);
        for (int i = 0; i < nIterations; i++) {
            if (i % 50 == 0)
                System.out.println("=> " + i);
            UnivariatePolynomial<BigInteger>
                    a = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(0, 15), bound, rnd),
                    b = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(0, 15), bound, rnd),
                    g = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(0, 15), bound, rnd);

            UnivariatePolynomial<Rational<BigInteger>> ar = a.mapCoefficients(Rings.Q, c -> new Rational<>(Rings.Z, c, BigInteger.valueOf(rndd.nextInt(2, 10))));
            UnivariatePolynomial<Rational<BigInteger>> br = b.mapCoefficients(Rings.Q, c -> new Rational<>(Rings.Z, c, BigInteger.valueOf(rndd.nextInt(2, 10))));
            UnivariatePolynomial<Rational<BigInteger>> gr = g.mapCoefficients(Rings.Q, c -> new Rational<>(Rings.Z, c, BigInteger.valueOf(rndd.nextInt(2, 10))));

            ar.multiply(gr);
            br.multiply(gr);

            long start = System.nanoTime();
            UnivariatePolynomial<Rational<BigInteger>>[] xgcd1 = ModularExtendedRationalGCD(ar, br);
            ratRec.addValue(System.nanoTime() - start);
            assertExtendedGCD(xgcd1, ar, br);

            start = System.nanoTime();
            UnivariatePolynomial<Rational<BigInteger>>[] xgcd2 = ModularExtendedResultantGCDInQ(ar, br);
            resultant.addValue(System.nanoTime() - start);
            assertExtendedGCD(xgcd2, ar, br);

            normalizeExtendedGCD(xgcd1);
            normalizeExtendedGCD(xgcd2);
            assertArrayEquals(xgcd1, xgcd2);

            assertTrue(g.degree <= xgcd1[0].degree);
            assertExtendedGCD(ar.increment(), br);
        }

        System.out.println("Rational reconstruction: " + ratRec);
        System.out.println("Resultant: " + resultant);
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
            UnivariatePolynomialZp64 a = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree),
                    modulus, rnd).multiply(1 + rnd.nextLong());
            UnivariatePolynomialZp64 b = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(minimalDegree, maximalDegree),
                    modulus, rnd).multiply(1 + rnd.nextLong());
            try {
                long start;

                start = System.nanoTime();
                UnivariatePolynomialZp64[] expected = ExtendedEuclidGCD(a, b);
                euclid.addValue(System.nanoTime() - start);

                start = System.nanoTime();
                UnivariatePolynomialZp64[] actual = ExtendedHalfGCD(a, b);
                half.addValue(System.nanoTime() - start);

                assertArrayEquals(expected, actual);

                if (!expected[0].isConstant())
                    ++nonTrivGCD;
            } catch (Throwable tr) {
                System.out.println("UnivariatePolynomial<BigInteger>." + a.toStringForCopy() + ".modulus(" + modulus + ");");
                System.out.println("UnivariatePolynomial<BigInteger>." + b.toStringForCopy() + ".modulus(" + modulus + ");");
                throw tr;
            }
        }
        System.out.println("Non-trivial gcds: " + nonTrivGCD);
        System.out.println("Euclid:  " + TimeUnits.statisticsNanotime(euclid));
        System.out.println("HalfGCD: " + TimeUnits.statisticsNanotime(half));
    }

    @Test
    public void test37_algext() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(UnivariatePolynomial.create(Q, Q.valueOf(-2), Q.valueOf(0), Q.valueOf(1)));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + s*x^5) * ( 1 + s*x^2 + 12*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + s*x^5) * ( 14 - s*x + 2*x^17)");
        assertEquals(coder.parse("-1+1/2*s+x^5"), UnivariateGCD.HalfGCD(a, b));
    }

    @Test
    public void test38_algext() {
        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicExtension(UnivariatePolynomial.create(Z, Z.valueOf(-2), Z.valueOf(0), Z.valueOf(1)));
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> a = coder.parse("(1 - s + x^5) * ( 1 + s*x^2 + 12*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> b = coder.parse("(1 - s + x^5) * ( 14 - s*x + 2*x^17)");

        assertEquals(coder.parse("1 - s + x^5"), UnivariateGCD.gcdAssociateInNumberField0(a, b));
    }

    @Test
    public void test39_algext() {
        UnivariatePolynomial<BigInteger> minimalPoly = UnivariatePolynomial.create(-2, 0, 0, 0, 0, 0, 1);
        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> a = coder.parse("(1 - s + (12 - 3*s^5) * x^5) * ( 1 + s*x^2 + 12*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> b = coder.parse("(1 - s + (12 - 3*s^5) * x^5) * ( 14 - s*x + 2*s*x^17)");

        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcd = UnivariateGCD.gcdAssociateInNumberField(a, b);

        assertTrue(pseudoDivideAndRemainder(a, gcd, true)[1].isZero());
        assertTrue(pseudoDivideAndRemainder(b, gcd, true)[1].isZero());
    }

    @Test
    @Benchmark
    public void test40_algext() {
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly = UnivariatePolynomial.create(-2, 0, 0, 0, 0, 0, 1).mapCoefficients(Q, Q::valueOfBigInteger);
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + (12133  - 3*s^5) * x^5) * ( 143 + s*x^2 + 12*s^6*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + (12133- 3*s^5) * x^5) * ( 14 - s*x + 2213*s*x^17)^2");

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> mod = UnivariateGCD.PolynomialGCDInNumberField(a, b);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> half = UnivariateGCD.HalfGCD(a, b);
            System.out.println("Half GCD: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(mod, half);
            System.out.println();
        }
    }

    @Test
    public void test41algext() {
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly =
                UnivariatePolynomial.create(2, 0, 0, 0, 0, 0, 5).mapCoefficients(Q, Q::valueOfBigInteger);

        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + (1 - 3*s^5) * x^5) * ( 3 + s*x^2 + 12*s^6*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + (1 - 3*s^5) * x^5) * ( 14 - s*x + 2*s*x^17)^2");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> gcd = UnivariateGCD.PolynomialGCDInNumberField(a, b);
        assertGCD(a, b, gcd);
    }

    @Test
    @Benchmark
    public void test42_algext() {
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly = UnivariatePolynomial.create(-2, 0, 0, 0, 0, 0, 5).mapCoefficients(Q, Q::valueOfBigInteger);
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("(1 - s + (12133  - 3*s^5) * x^5) * ( 143 + s*x^2 + 12*s^6*x^5)");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("(1 - s + (12133- 3*s^5) * x^5) * ( 14 - s*x + 2213*s*x^17)^2");

        for (int i = 0; i < 1; ++i) {
            long start;
            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> mod = UnivariateGCD.PolynomialGCDInNumberField(a, b);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> half = UnivariateGCD.HalfGCD(a, b);
            System.out.println("Half GCD: " + nanosecondsToString(System.nanoTime() - start));

            assertEquals(mod, half);
            System.out.println();
        }
    }

    @Test
    public void test43_algext_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(30, 100); i++) {
            UnivariatePolynomial<Rational<BigInteger>> minimalPoly =
                    IrreduciblePolynomials.randomIrreduciblePolynomialOverZ(rndd.nextInt(1, 5), rnd)
                            .mapCoefficients(Q, Q::valueOfBigInteger);
            AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = new AlgebraicNumberField<>(minimalPoly);

//            Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");
//            UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
//            Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

            for (int j = 0; j < 3; ++j) {
                UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                        a = randomPolyOverAlgExt(field, 10, 30, rndd),
                        b = randomPolyOverAlgExt(field, 10, 30, rndd),
                        gcd = randomPolyOverAlgExt(field, 10, 30, rndd);

                a = a.multiply(gcd);
                b = b.multiply(gcd);

                long start = System.nanoTime();
                UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> mgcd = PolynomialGCDInNumberField(a, b);
                System.out.print(i + ": " + nanosecondsToString(System.nanoTime() - start) + "      ");

                UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>[] qr = UnivariateDivision.divideAndRemainder(mgcd, gcd, true);
                assertNotNull(qr);
                assertTrue(qr[1].isZero());
                assertGCD(a, b, mgcd);
            }
            System.out.println();
        }
    }

    @Test
    @Benchmark
    public void test44_algext_performance() {
        UnivariatePolynomial<Rational<BigInteger>> minimalPoly =
                UnivariatePolynomial.create(31229703,
                        31584466, 9500649, 9480702,
                        23265262, 5568454, 3392530,
                        30401154, 15298203, 25411939,
                        30401154, 15298203, 25411939,
                        31584466, 9500649, 9480702,
                        30401154, 15298203, 25411939,
                        30401154, 15298203, 25411939, 1).mapCoefficients(Q, Q::valueOfBigInteger);
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(minimalPoly);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> a = coder.parse("-2629984078 - 2747848492*s - 826556509*s^2 - 824821066*s^3 - 2024077758*s^4 - 484455432*s^5 - 295150100*s^6 - 2644900377*s^7 - 1330943630*s^8 - 2210838750*s^9 + (-2539295050 - 2653095078*s - 798054533*s^2 - 796378949*s^3 - 1954282008*s^4 - 467750096*s^5 - 284972456*s^6 - 2553696896*s^7 - 1285048960*s^8 - 2134602853*s^9)*x + (2750903023 + 2874186377*s + 864559081*s^2 + 862743871*s^3 + 2117138752*s^4 + 506729319*s^5 + 308720296*s^6 + 2766505047*s^7 + 1392136388*s^8 + 2312486542*s^9)*x^2 + (-2902051419 - 3032108714*s - 912062266*s^2 - 910147339*s^3 - 2233465130*s^4 - 534571658*s^5 - 325682847*s^6 - 2918510751*s^7 - 1468627461*s^8 - 2439546109*s^9)*x^3 + (-2025390071 - 2116159152*s - 636543408*s^2 - 635207023*s^3 - 1558772527*s^4 - 373086362*s^5 - 227299470*s^6 - 2036877303*s^7 - 1024979692*s^8 - 1702599819*s^9)*x^4 + (-2962510945 - 3095277661*s - 931063590*s^2 - 929108773*s^3 - 2279995709*s^4 - 545708558*s^5 - 332467948*s^6 - 2979313000*s^7 - 1499223872*s^8 - 2490370012*s^9)*x^5 + (-1934700903 - 2021405804*s - 608041537*s^2 - 606764926*s^3 - 1488976751*s^4 - 356381007*s^5 - 217121886*s^6 - 1945673813*s^7 - 979084946*s^8 - 1626364082*s^9)*x^6 + (-997580106 - 1042287350*s - 313521458*s^2 - 312863069*s^3 - 767753632*s^4 - 183758950*s^5 - 111953411*s^6 - 1003238154*s^7 - 504840626*s^8 - 838594076*s^9)*x^7 + (-1178958351 - 1231794140*s - 370525402*s^2 - 369747315*s^3 - 907345216*s^4 - 217169661*s^5 - 132308571*s^6 - 1185645024*s^7 - 596629988*s^8 - 991065545*s^9)*x^8 + (-2025390077 - 2116159255*s - 636543412*s^2 - 635207025*s^3 - 1558772532*s^4 - 373086332*s^5 - 227299454*s^6 - 2036877243*s^7 - 1024979650*s^8 - 1702599890*s^9)*x^9 + (1118499087 + 1168625309*s + 351524046*s^2 + 350785991*s^3 + 860814787*s^4 + 206032829*s^5 + 125523644*s^6 + 1124842769*s^7 + 566033562*s^8 + 940241738*s^9)*x^10 + (-2267227821 - 2368834898*s - 712548607*s^2 - 711052633*s^3 - 1744894553*s^4 - 417633996*s^5 - 254439748*s^6 - 2280086614*s^7 - 1147365191*s^8 - 1905895490*s^9)*x^11 + (1239417854 + 1294963126*s + 389526531*s^2 + 388708754*s^3 + 953875698*s^4 + 228306652*s^5 + 139093768*s^6 + 1246447343*s^7 + 627226377*s^8 + 1041889549*s^9)*x^12 + (-906891056 - 947533928*s - 285019558*s^2 - 284421059*s^3 - 697957802*s^4 - 167053574*s^5 - 101775989*s^6 - 912034527*s^7 - 458946168*s^8 - 762358095*s^9)*x^13 + (-1844011957 - 1926652508*s - 579539643*s^2 - 578322774*s^3 - 1419181047*s^4 - 339675694*s^5 - 206944325*s^6 - 1854470337*s^7 - 933190317*s^8 - 1550128270*s^9)*x^14 + (1602174351 + 1673976642*s + 503534338*s^2 + 502477269*s^3 + 1233058976*s^4 + 295128130*s^5 + 179804110*s^6 + 1611261193*s^7 + 810804831*s^8 + 1346832799*s^9)*x^15 + (-453445531 - 473766944*s - 142509715*s^2 - 142210541*s^3 - 348978962*s^4 - 83526711*s^5 - 50887854*s^6 - 456017293*s^7 - 229473009*s^8 - 381179039*s^9)*x^16 + (-2509065293 - 2621510610*s - 788553921*s^2 - 786898292*s^3 - 1931016741*s^4 - 462181736*s^5 - 281579982*s^6 - 2523295772*s^7 - 1269750839*s^8 - 2109190927*s^9)*x^17 + (-1964930601 - 2052990286*s - 617542158*s^2 - 616245716*s^3 - 1512241977*s^4 - 361949565*s^5 - 220514529*s^6 - 1976074947*s^7 - 994383246*s^8 - 1651775975*s^9)*x^18 + (2327687048 + 2432003909*s + 731549963*s^2 + 730014071*s^3 + 1791425258*s^4 + 428771016*s^5 + 261224799*s^6 + 2340888912*s^7 + 1177961640*s^8 + 1956719266*s^9)*x^19");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> b = coder.parse("-2297457369 - 2400419408*s - 722049354*s^2 - 720533283*s^3 - 1768159825*s^4 - 423202465*s^5 - 257832326*s^6 - 2310487767*s^7 - 1162663504*s^8 - 1931307309*s^9 + (-1118498918 - 1168625170*s - 351523998*s^2 - 350785953*s^3 - 860814620*s^4 - 206032844*s^5 - 125523543*s^6 - 1124842734*s^7 - 566033484*s^8 - 940241840*s^9)*x + (-634823841 - 663273730*s - 199513585*s^2 - 199094740*s^3 - 488570500*s^4 - 116937491*s^5 - 71243056*s^6 - 638424217*s^7 - 321262269*s^8 - 533650647*s^9)*x^2 + (-120918804 - 126337886*s - 38002569*s^2 - 37922791*s^3 - 93061011*s^4 - 22273824*s^5 - 13570040*s^6 - 121604571*s^7 - 61192741*s^8 - 101647722*s^9)*x^3 + (-2629984238 - 2747848535*s - 826556405*s^2 - 824821002*s^3 - 2024077765*s^4 - 484455463*s^5 - 295150103*s^6 - 2644900316*s^7 - 1330943621*s^8 - 2210838616*s^9)*x^4 + (60459502 + 63168844*s + 19001371*s^2 + 18961478*s^3 + 46530559*s^4 + 11136967*s^5 + 6785128*s^6 + 60802294*s^7 + 30596476*s^8 + 50823949*s^9)*x^5 + (-2599754450 - 2716264170*s - 817055796*s^2 - 815340331*s^3 - 2000812621*s^4 - 478886953*s^5 - 291757547*s^6 - 2614499288*s^7 - 1315645534*s^8 - 2185426692*s^9)*x^6 + (-1964930626 - 2052990342*s - 617542144*s^2 - 616245697*s^3 - 1512242053*s^4 - 361949415*s^5 - 220514397*s^6 - 1976074963*s^7 - 994383244*s^8 - 1651775954*s^9)*x^7 + (-1571944495 - 1642392288*s - 494033670*s^2 - 492996488*s^3 - 1209793576*s^4 - 289559689*s^5 - 176411479*s^6 - 1580860073*s^7 - 795506501*s^8 - 1321420729*s^9)*x^8 + (2811362414 + 2937355383*s + 883560380*s^2 + 881705200*s^3 + 2163669280*s^4 + 517866158*s^5 + 315505366*s^6 + 2827307413*s^7 + 1422732958*s^8 + 2363310334*s^9)*x^9 + (2116079261 + 2210912639*s + 665045462*s^2 + 663649187*s^3 + 1628568308*s^4 + 389791803*s^5 + 237477108*s^6 + 2128080876*s^7 + 1070874245*s^8 + 1778835730*s^9)*x^10 + (-2236998012 - 2337250400*s - 703047951*s^2 - 701571937*s^3 - 1721629439*s^4 - 412065534*s^5 - 251047253*s^6 - 2249685310*s^7 - 1132067031*s^8 - 1880483482*s^9)*x^11 + (-2629984108 - 2747848541*s - 826556482*s^2 - 824821041*s^3 - 2024077781*s^4 - 484455457*s^5 - 295150030*s^6 - 2644900382*s^7 - 1330943563*s^8 - 2210838627*s^9)*x^12");
        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> g = coder.parse("-2811362312 - 2937355319*s - 883560337*s^2 - 881705356*s^3 - 2163669316*s^4 - 517866198*s^5 - 315505348*s^6 - 2827307397*s^7 - 1422732833*s^8 - 2363310305*s^9 + (-1027809863 - 1073871832*s - 323021979*s^2 - 322343852*s^3 - 791018874*s^4 - 189327442*s^5 - 115345990*s^6 - 1033639203*s^7 - 520138994*s^8 - 864005935*s^9)*x + (755742670 + 789611730*s + 237516324*s^2 + 237017568*s^3 + 581631609*s^4 + 139211399*s^5 + 84813301*s^6 + 760028909*s^7 + 382455125*s^8 + 635298486*s^9)*x^2 + (-2569524725 - 2684679604*s - 807555125*s^2 - 805859617*s^3 - 1977547175*s^4 - 473318570*s^5 - 288364973*s^6 - 2584098178*s^7 - 1300347246*s^8 - 2160014731*s^9)*x^3 + (-90689044 - 94753334*s - 28502018*s^2 - 28442203*s^3 - 69795734*s^4 - 16705297*s^5 - 10177621*s^6 - 91203421*s^7 - 45894548*s^8 - 76235851*s^9)*x^4 + (-60459337 - 63168954*s - 19001279*s^2 - 18961469*s^3 - 46530495*s^4 - 11136820*s^5 - 6785055*s^6 - 60802219*s^7 - 30596432*s^8 - 50823814*s^9)*x^5 + (-2206768278 - 2305665984*s - 693547285*s^2 - 692091191*s^3 - 1698364119*s^4 - 406497048*s^5 - 247654596*s^6 - 2219284280*s^7 - 1116768829*s^8 - 1855071634*s^9)*x^6 + (-2992740680 - 3126862134*s - 940564345*s^2 - 938589464*s^3 - 2303260921*s^4 - 551276902*s^5 - 335860390*s^6 - 3009714187*s^7 - 1514522192*s^8 - 2515781892*s^9)*x^7 + (-937120739 - 979118398*s - 294520061*s^2 - 293901791*s^3 - 721223167*s^4 - 172622016*s^5 - 105168421*s^6 - 942435722*s^7 - 474244292*s^8 - 787770068*s^9)*x^8 + (-2569524839 - 2684679579*s - 807555086*s^2 - 805859625*s^3 - 1977547366*s^4 - 473318589*s^5 - 288365123*s^6 - 2584098088*s^7 - 1300347180*s^8 - 2160014805*s^9)*x^9 + (-302297052 - 315844685*s - 95006420*s^2 - 94806952*s^3 - 232652706*s^4 - 55684503*s^5 - 33925202*s^6 - 304011605*s^7 - 152982033*s^8 - 254119332*s^9)*x^10 + (362756428 + 379013624*s + 114007873*s^2 + 113768522*s^3 + 279183228*s^4 + 66821400*s^5 + 40710454*s^6 + 364813756*s^7 + 183578437*s^8 + 304943317*s^9)*x^11 + (-2327687116 - 2432003792*s - 731549932*s^2 - 730014098*s^3 - 1791425234*s^4 - 428770976*s^5 - 261224810*s^6 - 2340888824*s^7 - 1177961606*s^8 - 1956719264*s^9)*x^12 + (-2962510886 - 3095277709*s - 931063515*s^2 - 929108855*s^3 - 2279995592*s^4 - 545708444*s^5 - 332467892*s^6 - 2979313007*s^7 - 1499223814*s^8 - 2490369949*s^9)*x^13 + (1904471341 + 1989821396*s + 598540924*s^2 + 597284317*s^3 + 1465711520*s^4 + 350812607*s^5 + 213729418*s^6 + 1915272618*s^7 + 963786874*s^8 + 1600952190*s^9)*x^14 + (1420796032 + 1484469924*s + 446530582*s^2 + 445593084*s^3 + 1093467219*s^4 + 261717435*s^5 + 159448870*s^6 + 1428854284*s^7 + 719015466*s^8 + 1194361180*s^9)*x^15 + (-2327687165 - 2432003807*s - 731549938*s^2 - 730013976*s^3 - 1791425083*s^4 - 428770908*s^5 - 261224823*s^6 - 2340888813*s^7 - 1177961603*s^8 - 1956719254*s^9)*x^16 + (1390566435 + 1452885413*s + 437029859*s^2 + 436112351*s^3 + 1070202123*s^4 + 256148955*s^5 + 156056436*s^6 + 1398453101*s^7 + 703717376*s^8 + 1168949284*s^9)*x^17 + (-272067241 - 284260115*s - 85505836*s^2 - 85326368*s^3 - 209387338*s^4 - 50116022*s^5 - 30532678*s^6 - 273610327*s^7 - 137683742*s^8 - 228707356*s^9)*x^18 + (-1027809875 - 1073871903*s - 323021967*s^2 - 322343855*s^3 - 791018854*s^4 - 189327350*s^5 - 115345977*s^6 - 1033639175*s^7 - 520138874*s^8 - 864005908*s^9)*x^19 + (-2055619749 - 2147743611*s - 646044049*s^2 - 644687698*s^3 - 1582037829*s^4 - 378654787*s^5 - 230691998*s^6 - 2067278570*s^7 - 1040277772*s^8 - 1728011857*s^9)*x^20 + (-1783552468 - 1863483492*s - 560538209*s^2 - 559361408*s^3 - 1372650543*s^4 - 328538714*s^5 - 200159198*s^6 - 1793668083*s^7 - 902593935*s^8 - 1499304498*s^9)*x^21 + (-876661420 - 915949489*s - 275518819*s^2 - 274940409*s^3 - 674692521*s^4 - 161485158*s^5 - 98383273*s^6 - 881633387*s^7 - 443647831*s^8 - 736946230*s^9)*x^22 + (-2116079189 - 2210912556*s - 665045429*s^2 - 663649110*s^3 - 1628568244*s^4 - 389791730*s^5 - 237477099*s^6 - 2128080829*s^7 - 1070874203*s^8 - 1778835671*s^9)*x^23 + (-1239417745 - 1294963053*s - 389526527*s^2 - 388708825*s^3 - 953875650*s^4 - 228306522*s^5 - 139093680*s^6 - 1246447329*s^7 - 627226263*s^8 - 1041889490*s^9)*x^24 + (-2509065417 - 2621510653*s - 788553920*s^2 - 786898228*s^3 - 1931016745*s^4 - 462181615*s^5 - 281579891*s^6 - 2523295782*s^7 - 1269750822*s^8 - 2109190863*s^9)*x^25");

        UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
                ag = a.clone().multiply(g).multiply(a.clone().decrement()),
                bg = b.clone().multiply(g).multiply(b.clone().decrement());
        for (int i = 0; i < 2; ++i) {
            long start;
            start = System.nanoTime();
            UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>> mod = UnivariateGCD.PolynomialGCDInNumberField(ag, bg);
            System.out.println("Modular: " + nanosecondsToString(System.nanoTime() - start));
            assertTrue(remainder(mod, g, true).isZero());
        }

        // Modular: 7s
        // Modular: 5s
    }

    @SuppressWarnings("unchecked")
    private UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>
    randomPolyOverAlgExt(AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field,
                         int minDeg, int maxDeg,
                         RandomDataGenerator rndd) {
        return UnivariatePolynomial.create(field, IntStream.rangeClosed(0, rndd.nextInt(minDeg, maxDeg))
                .mapToObj(__ -> RandomUnivariatePolynomials.randomPoly(field.getMinimalPoly().degree, 100, rndd.getRandomGenerator()).toBigPoly()
                        .mapCoefficients(Q, Q::valueOfBigInteger)).toArray(UnivariatePolynomial[]::new));
    }

    static <T extends IUnivariatePolynomial<T>> void assertExtendedEuclidGCD(T a, T b) {
        assertExtendedGCD(ExtendedEuclidGCD(a, b), a, b);
    }

    static <T extends IUnivariatePolynomial<T>> T[] assertExtendedGCD(T a, T b) {
        T[] eea = PolynomialExtendedGCD(a, b);
        assertExtendedGCD(eea, a, b);
        return eea;
    }

    static <T extends IUnivariatePolynomial<T>> void assertExtendedGCD(T[] eea, T a, T b) {
        assertEquals(eea[0], a.clone().multiply(eea[1]).add(b.clone().multiply(eea[2])));
        assertEquals(eea[0].degree(), PolynomialGCD(a, b).degree());
    }
}