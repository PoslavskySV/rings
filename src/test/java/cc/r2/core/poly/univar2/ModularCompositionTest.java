package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.univar2.DivisionWithRemainder.InverseModMonomial;
import cc.r2.core.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;

import static cc.r2.core.poly.univar2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly.univar2.ModularComposition.*;
import static cc.r2.core.poly.univar2.gMutablePolynomial.asLongPolyZp;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 26/02/2017.
 */
public class ModularCompositionTest extends AbstractPolynomialTest {
    @Test
    public void testXPowers1() throws Exception {
        long modulus = 43;
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertXPowers(polyModulus);
    }

    @Test
    public void testXPowers2() throws Exception {
        long modulus = 43;
        lMutablePolynomialZp poly = lMutablePolynomialZ.create(1).modulus(modulus);
        assertXPowers(poly);
    }

    @Test
    public void testXPowers3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 300); i++)
            assertXPowers(RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), getModulusRandom(getRandomData().nextInt(2, 6)), rnd));
    }

    static void assertXPowers(lMutablePolynomialZp polyModulus) {
        ArrayList<lMutablePolynomialZp> xPowers = ModularComposition.xPowers(polyModulus, fastDivisionPreConditioning(polyModulus));
        for (int i = 0; i < xPowers.size(); i++) {
            lMutablePolynomialZp expected = PolynomialArithmetics.polyPowMod(
                    lMutablePolynomialZp.createMonomial(polyModulus.domain.modulus, 1, (int) polyModulus.domain.modulus),
                    i, polyModulus, false);
            assertEquals(expected, xPowers.get(i));
        }
    }

    @Test
    public void testPolyPowers1() throws Exception {
        long modulus = 43;
        lMutablePolynomialZp poly = lMutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus);
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertPolyPowers(poly, polyModulus, 10);
    }

    @Test
    public void testPolyPowers2() throws Exception {
        long modulus = 43;
        assertPolyPowers(lMutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus), lMutablePolynomialZ.create(1).modulus(modulus), 10);
        assertPolyPowers(lMutablePolynomialZ.create(1).modulus(modulus), lMutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus), 10);
    }


    @Test
    public void testPolyPowers3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < its(100, 300); i++) {
                lMutablePolynomialZp poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                lMutablePolynomialZp polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPolyPowers(poly, polyModulus, 10);
            }
        }
    }

    static void assertPolyPowers(lMutablePolynomialZp poly, lMutablePolynomialZp polyModulus, int nIterations) {
        ArrayList<lMutablePolynomialZp> hPowers = polyPowers(poly, polyModulus, fastDivisionPreConditioning(polyModulus), nIterations);
        for (int i = 0; i < hPowers.size(); i++) {
            lMutablePolynomialZp expected = PolynomialArithmetics.polyPowMod(poly, i, polyModulus, true);
            assertEquals(expected, hPowers.get(i));
        }
    }

    @Test
    public void testPowModulusMod1() throws Exception {
        long modulus = 43;
        lMutablePolynomialZp poly = lMutablePolynomialZ.create(1, 2, 3, 4, 4).modulus(modulus);
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 4, 1, 1).modulus(modulus);
        assertPowModulusMod(poly, polyModulus);
    }

    @Test
    public void testPowModulusMod2() throws Exception {
        long modulus = 43;
        assertPowModulusMod(lMutablePolynomialZ.create(1, 2, 3, 4, 1).modulus(modulus), lMutablePolynomialZ.create(1).modulus(modulus));
    }

    @Test
    public void testPowModulusMod3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < its(50, 300); i++) {
                lMutablePolynomialZp poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                lMutablePolynomialZp polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPowModulusMod(poly, polyModulus);
            }
        }
    }

    static void assertPowModulusMod(lMutablePolynomialZp poly, lMutablePolynomialZp polyModulus) {
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
        assertEquals(PolynomialArithmetics.polyPowMod(poly, poly.domain.modulus, polyModulus, true), powModulusMod(poly, polyModulus, invMod, xPowers(polyModulus, invMod)));
    }

    @Test
    public void testComposition1() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            lMutablePolynomialZp poly = lMutablePolynomialZ.create(1, 2, 3, 4, 45, 2, 1).modulus(modulus);
            lMutablePolynomialZp point = lMutablePolynomialZ.create(1, 2, 3).modulus(modulus);
            lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus);
            assertComposition(poly, point, polyModulus);
        }
    }

    @Test
    public void testComposition2() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            lMutablePolynomialZp poly = lMutablePolynomialZ.create(1).modulus(modulus);
            lMutablePolynomialZp point = lMutablePolynomialZ.create(1, 2, 1).modulus(modulus);
            lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus);
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
        RandomGenerator rnd = getRandom();
        for (long modulus : getSmallModulusArray(5)) {
            for (int i = 0; i < 100; i++) {
                lMutablePolynomialZp poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                lMutablePolynomialZp point = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                lMutablePolynomialZp polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertComposition(poly, point, polyModulus);
            }
        }
    }

    @Test
    public void testComposition3a() throws Exception {
        long modulus = 43;
        lMutablePolynomialZp poly = lMutablePolynomialZ.create(0, 0, 0, 1, 1).modulus(modulus);
        lMutablePolynomialZp point = lMutablePolynomialZ.create(1, 1, 0, 1, 0, 1, 0, 1, 1).modulus(modulus);
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(0, 1).modulus(modulus);
        assertComposition(poly, point, polyModulus);
    }

    @Test
    @Benchmark
    public void testComposition4Random_performance() throws Exception {
        DescriptiveStatistics hornerStats = new DescriptiveStatistics(), brenKungStats = new DescriptiveStatistics();
        RandomGenerator rnd = getRandom();
        long modulus = SmallPrimes.nextPrime(10000);
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                hornerStats.clear();
                brenKungStats.clear();
            }
            lMutablePolynomialZp poly = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            lMutablePolynomialZp point = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            lMutablePolynomialZp polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);

//            System.out.println("====");
            long start = System.nanoTime();
            lMutablePolynomialZp horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
            hornerStats.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            lMutablePolynomialZp brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
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
        RandomGenerator rnd = getRandom();
        for (long modulus : getSmallModulusArray(50)) {
            BigInteger bModulus = BigInteger.valueOf(modulus);
            for (int i = 0; i < 5; i++) {
                gMutablePolynomial<BigInteger> poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                gMutablePolynomial<BigInteger> point = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                gMutablePolynomial<BigInteger> polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), bModulus, rnd);
                gMutablePolynomial<BigInteger> bComp = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
                lMutablePolynomialZp lPolyModulus = asLongPolyZp(polyModulus);
                lMutablePolynomialZp lComp = compositionBrentKung(asLongPolyZp(poly), asLongPolyZp(point), lPolyModulus, fastDivisionPreConditioning(lPolyModulus));
                assertEquals(lComp, asLongPolyZp(bComp));
            }
        }
    }

    static void assertComposition(lMutablePolynomialZp poly, lMutablePolynomialZp point, lMutablePolynomialZp polyModulus) {
        lMutablePolynomialZp horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
        lMutablePolynomialZp brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
        lMutablePolynomialZp plain = plainComposition(poly, point, polyModulus);
        assertEquals(plain, horner);
        assertEquals(plain, brentKung);
    }

    private static lMutablePolynomialZp plainComposition(lMutablePolynomialZp poly, lMutablePolynomialZp point, lMutablePolynomialZp polyModulus) {
        lMutablePolynomialZp res = poly.createZero();
        for (int i = 0; i <= poly.degree; i++)
            res = PolynomialArithmetics.polyMod(res.add(PolynomialArithmetics.polyPowMod(point, i, polyModulus, true).multiply(poly.data[i])), polyModulus, false);
        return res;
    }
}