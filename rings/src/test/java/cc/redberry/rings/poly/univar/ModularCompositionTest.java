package cc.redberry.rings.poly.univar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;

import static cc.redberry.rings.poly.univar.ModularComposition.*;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 26/02/2017.
 */
public class ModularCompositionTest extends AUnivariateTest {
    @Test
    public void testXPowers1() throws Exception {
        long modulus = 43;
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertXPowers(polyModulus);
    }

    @Test
    public void testXPowers2() throws Exception {
        long modulus = 43;
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1).modulus(modulus);
        assertXPowers(poly);
    }

    @Test
    public void testXPowers3Random() throws Exception {
        RandomGenerator rnd = AbstractTest.getRandom();
        for (int i = 0; i < AbstractTest.its(100, 300); i++)
            assertXPowers(RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), getModulusRandom(AbstractTest.getRandomData().nextInt(2, 6)), rnd));
    }

    static void assertXPowers(UnivariatePolynomialZp64 polyModulus) {
        ArrayList<UnivariatePolynomialZp64> xPowers = ModularComposition.xPowers(polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
        for (int i = 0; i < xPowers.size(); i++) {
            UnivariatePolynomialZp64 expected = UnivariatePolynomialArithmetic.polyPowMod(
                    UnivariatePolynomialZp64.monomial(polyModulus.ring.modulus, 1, (int) polyModulus.ring.modulus),
                    i, polyModulus, false);
            assertEquals(expected, xPowers.get(i));
        }
    }

    @Test
    public void testPolyPowers1() throws Exception {
        long modulus = 43;
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4).modulus(modulus);
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertPolyPowers(poly, polyModulus, 10);
    }

    @Test
    public void testPolyPowers2() throws Exception {
        long modulus = 43;
        assertPolyPowers(UnivariatePolynomialZ64.create(1, 2, 3, 4).modulus(modulus), UnivariatePolynomialZ64.create(1).modulus(modulus), 10);
        assertPolyPowers(UnivariatePolynomialZ64.create(1).modulus(modulus), UnivariatePolynomialZ64.create(1, 2, 3, 1).modulus(modulus), 10);
    }


    @Test
    public void testPolyPowers3Random() throws Exception {
        RandomGenerator rnd = AbstractTest.getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < AbstractTest.its(100, 300); i++) {
                UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                UnivariatePolynomialZp64 polyModulus = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPolyPowers(poly, polyModulus, 10);
            }
        }
    }

    static void assertPolyPowers(UnivariatePolynomialZp64 poly, UnivariatePolynomialZp64 polyModulus, int nIterations) {
        ArrayList<UnivariatePolynomialZp64> hPowers = polyPowers(poly, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus), nIterations);
        for (int i = 0; i < hPowers.size(); i++) {
            UnivariatePolynomialZp64 expected = UnivariatePolynomialArithmetic.polyPowMod(poly, i, polyModulus, true);
            assertEquals(expected, hPowers.get(i));
        }
    }

    @Test
    public void testPowModulusMod1() throws Exception {
        long modulus = 43;
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4, 4).modulus(modulus);
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 4, 1, 1).modulus(modulus);
        assertPowModulusMod(poly, polyModulus);
    }

    @Test
    public void testPowModulusMod2() throws Exception {
        long modulus = 43;
        assertPowModulusMod(UnivariatePolynomialZ64.create(1, 2, 3, 4, 1).modulus(modulus), UnivariatePolynomialZ64.create(1).modulus(modulus));
    }

    @Test
    public void testPowModulusMod3Random() throws Exception {
        RandomGenerator rnd = AbstractTest.getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < AbstractTest.its(50, 300); i++) {
                UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                UnivariatePolynomialZp64 polyModulus = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPowModulusMod(poly, polyModulus);
            }
        }
    }

    static void assertPowModulusMod(UnivariatePolynomialZp64 poly, UnivariatePolynomialZp64 polyModulus) {
        UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(polyModulus);
        assertEquals(UnivariatePolynomialArithmetic.polyPowMod(poly, poly.ring.modulus, polyModulus, true), powModulusMod(poly, polyModulus, invMod, xPowers(polyModulus, invMod)));
    }

    @Test
    public void testComposition1() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1, 2, 3, 4, 45, 2, 1).modulus(modulus);
            UnivariatePolynomialZp64 point = UnivariatePolynomialZ64.create(1, 2, 3).modulus(modulus);
            UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 1).modulus(modulus);
            assertComposition(poly, point, polyModulus);
        }
    }

    @Test
    public void testComposition2() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1).modulus(modulus);
            UnivariatePolynomialZp64 point = UnivariatePolynomialZ64.create(1, 2, 1).modulus(modulus);
            UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 1).modulus(modulus);
            assertComposition(poly, point, polyModulus);
            assertComposition(poly, polyModulus, point);
            assertComposition(point, poly, polyModulus);
            assertComposition(point, polyModulus, poly);
            assertComposition(polyModulus, point, poly);
            assertComposition(polyModulus, poly, point);
        }
    }

    @Test
    public void testComposition3Random() throws Exception {
        RandomGenerator rnd = AbstractTest.getRandom();
        for (long modulus : getSmallModulusArray(5)) {
            for (int i = 0; i < 100; i++) {
                UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                UnivariatePolynomialZp64 point = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                UnivariatePolynomialZp64 polyModulus = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertComposition(poly, point, polyModulus);
            }
        }
    }

    @Test
    public void testComposition3a() throws Exception {
        long modulus = 43;
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(0, 0, 0, 1, 1).modulus(modulus);
        UnivariatePolynomialZp64 point = UnivariatePolynomialZ64.create(1, 1, 0, 1, 0, 1, 0, 1, 1).modulus(modulus);
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(0, 1).modulus(modulus);
        assertComposition(poly, point, polyModulus);
    }

    @Test
    @Benchmark
    public void testComposition4Random_performance() throws Exception {
        DescriptiveStatistics hornerStats = new DescriptiveStatistics(), brenKungStats = new DescriptiveStatistics();
        RandomGenerator rnd = AbstractTest.getRandom();
        long modulus = SmallPrimes.nextPrime(10000);
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                hornerStats.clear();
                brenKungStats.clear();
            }
            UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            UnivariatePolynomialZp64 point = RandomUnivariatePolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            UnivariatePolynomialZp64 polyModulus = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);

//            System.out.println("====");
            long start = System.nanoTime();
            UnivariatePolynomialZp64 horner = compositionHorner(poly, point, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
            hornerStats.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            UnivariatePolynomialZp64 brentKung = compositionBrentKung(poly, point, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
            brenKungStats.addValue(System.nanoTime() - start);

            assertEquals(horner, brentKung);
        }

        System.out.println("Horner");
        System.out.println(hornerStats.getMean());

        System.out.println("BrentKung");
        System.out.println(brenKungStats.getMean());
    }

    @Test
    public void testComposition3Random_big_poly() throws Exception {
        RandomGenerator rnd = AbstractTest.getRandom();
        for (long modulus : getSmallModulusArray(50)) {
            BigInteger bModulus = BigInteger.valueOf(modulus);
            for (int i = 0; i < 5; i++) {
                UnivariatePolynomial<BigInteger> poly = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                UnivariatePolynomial<BigInteger> point = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                UnivariatePolynomial<BigInteger> polyModulus = RandomUnivariatePolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                UnivariatePolynomial<BigInteger> bComp = compositionBrentKung(poly, point, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
                UnivariatePolynomialZp64 lPolyModulus = UnivariatePolynomial.asOverZp64(polyModulus);
                UnivariatePolynomialZp64 lComp = compositionBrentKung(UnivariatePolynomial.asOverZp64(poly), UnivariatePolynomial.asOverZp64(point), lPolyModulus, UnivariateDivision.fastDivisionPreConditioning(lPolyModulus));
                Assert.assertEquals(lComp, UnivariatePolynomial.asOverZp64(bComp));
            }
        }
    }

    static void assertComposition(UnivariatePolynomialZp64 poly, UnivariatePolynomialZp64 point, UnivariatePolynomialZp64 polyModulus) {
        UnivariatePolynomialZp64 horner = compositionHorner(poly, point, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
        UnivariatePolynomialZp64 brentKung = compositionBrentKung(poly, point, polyModulus, UnivariateDivision.fastDivisionPreConditioning(polyModulus));
        UnivariatePolynomialZp64 plain = plainComposition(poly, point, polyModulus);
        assertEquals(plain, horner);
        assertEquals(plain, brentKung);
    }

    private static UnivariatePolynomialZp64 plainComposition(UnivariatePolynomialZp64 poly, UnivariatePolynomialZp64 point, UnivariatePolynomialZp64 polyModulus) {
        UnivariatePolynomialZp64 res = poly.createZero();
        for (int i = 0; i <= poly.degree; i++)
            res = UnivariatePolynomialArithmetic.polyMod(res.add(UnivariatePolynomialArithmetic.polyPowMod(point, i, polyModulus, true).multiply(poly.data[i])), polyModulus, false);
        return res;
    }
}