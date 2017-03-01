package cc.r2.core.poly2;

import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly2.DivisionWithRemainder.InverseModMonomial;
import cc.r2.core.test.Benchmark;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Test;

import java.util.ArrayList;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.ModularComposition.*;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 26/02/2017.
 */
public class ModularCompositionTest extends AbstractPolynomialTest {
    @Test
    public void testXPowers1() throws Exception {
        long modulus = 43;
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertXPowers(polyModulus);
    }

    @Test
    public void testXPowers2() throws Exception {
        long modulus = 43;
        MutablePolynomialMod poly = MutablePolynomialZ.create(1).modulus(modulus);
        assertXPowers(poly);
    }

    @Test
    public void testXPowers3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 300); i++)
            assertXPowers(RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), getModulusRandom(getRandomData().nextInt(2, 6)), rnd));
    }

    static void assertXPowers(MutablePolynomialMod polyModulus) {
        ArrayList<MutablePolynomialMod> xPowers = ModularComposition.xPowers(polyModulus, fastDivisionPreConditioning(polyModulus));
        for (int i = 0; i < xPowers.size(); i++) {
            MutablePolynomialMod expected = PolynomialArithmetics.polyPowMod(
                    MutablePolynomialMod.createMonomial(polyModulus.modulus, 1, (int) polyModulus.modulus),
                    i, polyModulus, false);
            assertEquals(expected, xPowers.get(i));
        }
    }

    @Test
    public void testPolyPowers1() throws Exception {
        long modulus = 43;
        MutablePolynomialMod poly = MutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus);
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
        assertPolyPowers(poly, polyModulus, 10);
    }

    @Test
    public void testPolyPowers2() throws Exception {
        long modulus = 43;
        assertPolyPowers(MutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus), MutablePolynomialZ.create(1).modulus(modulus), 10);
        assertPolyPowers(MutablePolynomialZ.create(1).modulus(modulus), MutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus), 10);
    }


    @Test
    public void testPolyPowers3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < its(100, 300); i++) {
                MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomialMod polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPolyPowers(poly, polyModulus, 10);
            }
        }
    }

    static void assertPolyPowers(MutablePolynomialMod poly, MutablePolynomialMod polyModulus, int nIterations) {
        ArrayList<MutablePolynomialMod> hPowers = polyPowers(poly, polyModulus, fastDivisionPreConditioning(polyModulus), nIterations);
        for (int i = 0; i < hPowers.size(); i++) {
            MutablePolynomialMod expected = PolynomialArithmetics.polyPowMod(poly, i, polyModulus, true);
            assertEquals(expected, hPowers.get(i));
        }
    }

    @Test
    public void testPowModulusMod1() throws Exception {
        long modulus = 43;
        MutablePolynomialMod poly = MutablePolynomialZ.create(1, 2, 3, 4, 4).modulus(modulus);
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 4, 1, 1).modulus(modulus);
        assertPowModulusMod(poly, polyModulus);
    }

    @Test
    public void testPowModulusMod2() throws Exception {
        long modulus = 43;
        assertPowModulusMod(MutablePolynomialZ.create(1, 2, 3, 4, 1).modulus(modulus), MutablePolynomialZ.create(1).modulus(modulus));
    }

    @Test
    public void testPowModulusMod3Random() throws Exception {
        RandomGenerator rnd = getRandom();
        for (long modulus : getModulusArray(3, 2, 45)) {
            for (int i = 0; i < its(50, 300); i++) {
                MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomialMod polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPowModulusMod(poly, polyModulus);
            }
        }
    }

    static void assertPowModulusMod(MutablePolynomialMod poly, MutablePolynomialMod polyModulus) {
        InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
        assertEquals(PolynomialArithmetics.polyPowMod(poly, poly.modulus, polyModulus, true), powModulusMod(poly, polyModulus, invMod, xPowers(polyModulus, invMod)));
    }

    @Test
    public void testComposition1() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            MutablePolynomialMod poly = MutablePolynomialZ.create(1, 2, 3, 4, 45, 2, 1).modulus(modulus);
            MutablePolynomialMod point = MutablePolynomialZ.create(1, 2, 3).modulus(modulus);
            MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus);
            assertComposition(poly, point, polyModulus);
        }
    }

    @Test
    public void testComposition2() throws Exception {
        for (long modulus : getModulusArray(1, 1, 50)) {
            MutablePolynomialMod poly = MutablePolynomialZ.create(1).modulus(modulus);
            MutablePolynomialMod point = MutablePolynomialZ.create(1, 2, 1).modulus(modulus);
            MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 1).modulus(modulus);
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
                MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomialMod point = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomialMod polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertComposition(poly, point, polyModulus);
            }
        }
    }

    @Test
    public void testComposition3a() throws Exception {
        long modulus = 43;
        MutablePolynomialMod poly = MutablePolynomialZ.create(0, 0, 0, 1, 1).modulus(modulus);
        MutablePolynomialMod point = MutablePolynomialZ.create(1, 1, 0, 1, 0, 1, 0, 1, 1).modulus(modulus);
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(0, 1).modulus(modulus);
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
            MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            MutablePolynomialMod point = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            MutablePolynomialMod polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);

//            System.out.println("====");
            long start = System.nanoTime();
            MutablePolynomialMod horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
            hornerStats.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            MutablePolynomialMod brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
            brenKungStats.addValue(System.nanoTime() - start);

            assertEquals(horner, brentKung);
        }

        System.out.println("Horner");
        System.out.println(hornerStats.getMean());

        System.out.println("BrentKung");
        System.out.println(brenKungStats.getMean());
    }

    static void assertComposition(MutablePolynomialMod poly, MutablePolynomialMod point, MutablePolynomialMod polyModulus) {
        MutablePolynomialMod horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
        MutablePolynomialMod brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus));
        MutablePolynomialMod plain = plainComposition(poly, point, polyModulus);
        assertEquals(plain, horner);
        assertEquals(plain, brentKung);
    }

    private static MutablePolynomialMod plainComposition(MutablePolynomialMod poly, MutablePolynomialMod point, MutablePolynomialMod polyModulus) {
        MutablePolynomialMod res = poly.createZero();
        for (int i = 0; i <= poly.degree; i++)
            res = PolynomialArithmetics.polyMod(res.add(PolynomialArithmetics.polyPowMod(point, i, polyModulus, true).multiply(poly.data[i])), polyModulus, false);
        return res;
    }
}