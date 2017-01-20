package cc.r2.core.polynomial;

import cc.r2.core.number.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Ignore;
import org.junit.Test;

import java.util.ArrayList;

import static cc.r2.core.polynomial.DivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.ModularComposition.*;
import static cc.r2.core.polynomial.MutablePolynomial.createMonomial;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyMod;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyPowMod;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 20/01/2017.
 */
public class ModularCompositionTest {
    @Test
    public void testXPowers1() throws Exception {
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 1);
        long modulus = 43;
        assertXPowers(polyModulus, modulus);
    }

    @Test
    public void testXPowers2() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1);
        long modulus = 43;
        assertXPowers(poly, modulus);
    }

    @Test
    public void testXPowers3Random() throws Exception {
        RandomGenerator rnd = new Well1024a();
        long modulus = 7;
        for (int i = 0; i < 100; i++)
            assertXPowers(RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd), modulus);
    }

    static void assertXPowers(MutablePolynomial polyModulus, long modulus) {
        ArrayList<MutablePolynomial> xPowers = xPowers(polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus);
        for (int i = 0; i < xPowers.size(); i++) {
            MutablePolynomial expected = polyPowMod(createMonomial(1, (int) modulus), i, polyModulus, modulus, false);
            assertEquals(expected, xPowers.get(i));
        }
    }

    @Test
    public void testPolyPowers1() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4);
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 4, 5, 6, 1);
        long modulus = 43;
        assertPolyPowers(poly, polyModulus, modulus, 10);
    }

    @Test
    public void testPolyPowers2() throws Exception {
        long modulus = 43;
        assertPolyPowers(MutablePolynomial.create(1, 2, 3, 4), MutablePolynomial.create(1), modulus, 10);
        assertPolyPowers(MutablePolynomial.create(1), MutablePolynomial.create(1, 2, 3, 1), modulus, 10);
    }

    @Test
    public void testPolyPowers3Random() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (long modulus : new long[]{2, 3, 7, 43}) {
            for (int i = 0; i < 100; i++) {
                MutablePolynomial poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomial polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPolyPowers(poly, polyModulus, modulus, 10);
            }
        }
    }

    static void assertPolyPowers(MutablePolynomial poly, MutablePolynomial polyModulus, long modulus, int nIterations) {
        ArrayList<MutablePolynomial> hPowers = polyPowers(poly, polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus, nIterations);
        for (int i = 0; i < hPowers.size(); i++) {
            MutablePolynomial expected = polyPowMod(poly, i, polyModulus, modulus, true);
            assertEquals(expected, hPowers.get(i));
        }
    }

    @Test
    public void testPowModulusMod1() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 4);
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 4, 1, 1);
        long modulus = 43;
        assertPowModulusMod(poly, polyModulus, modulus);
    }

    @Test
    public void testPowModulusMod2() throws Exception {
        long modulus = 43;
//        assertPowModulusMod(MutablePolynomial.create(1), MutablePolynomial.create(1, 2, 3, 4, 1), modulus);
        assertPowModulusMod(MutablePolynomial.create(1, 2, 3, 4, 1), MutablePolynomial.create(1), modulus);
    }

    @Test
    public void testPowModulusMod3Random() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (long modulus : new long[]{2, 3, 7, 43}) {
            for (int i = 0; i < 100; i++) {
                MutablePolynomial poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomial polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertPowModulusMod(poly, polyModulus, modulus);
            }
        }
    }

    static void assertPowModulusMod(MutablePolynomial poly, MutablePolynomial polyModulus, long modulus) {
        DivideAndRemainder.InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus, modulus);
        assertEquals(polyPowMod(poly, modulus, polyModulus, modulus, true), powModulusMod(poly, polyModulus, invMod, modulus, xPowers(polyModulus, invMod, modulus)));
    }

    @Test
    public void testComposition1() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1, 2, 3, 4, 45, 2, 1);
        MutablePolynomial point = MutablePolynomial.create(1, 2, 3);
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 1);
        long modulus = 43;
        assertComposition(poly, point, polyModulus, modulus);
    }

    @Test
    public void testComposition2() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(1);
        MutablePolynomial point = MutablePolynomial.create(1, 2, 1);
        MutablePolynomial polyModulus = MutablePolynomial.create(1, 2, 3, 1);
        long modulus = 43;
        assertComposition(poly, point, polyModulus, modulus);
        assertComposition(poly, polyModulus, point, modulus);
        assertComposition(point, poly, polyModulus, modulus);
        assertComposition(point, polyModulus, poly, modulus);
        assertComposition(polyModulus, point, poly, modulus);
        assertComposition(polyModulus, poly, point, modulus);
    }

    @Test
    public void testComposition3Random() throws Exception {
        RandomGenerator rnd = new Well1024a();
        for (long modulus : new long[]{2, 3, 7, 43}) {
            for (int i = 0; i < 100; i++) {
                MutablePolynomial poly = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomial point = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                MutablePolynomial polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);
                assertComposition(poly, point, polyModulus, modulus);
            }
        }
    }

    @Test
    public void testComposition3a() throws Exception {
        MutablePolynomial poly = MutablePolynomial.create(0, 0, 0, 1, 1);
        MutablePolynomial point = MutablePolynomial.create(1, 1, 0, 1, 0, 1, 0, 1, 1);
        MutablePolynomial polyModulus = MutablePolynomial.create(0, 1);
        long modulus = 43;
        assertComposition(poly, point, polyModulus, modulus);
    }


    @Ignore
    @Test
    public void testComposition4Random_performace() throws Exception {
        DescriptiveStatistics hornerStats = new DescriptiveStatistics(), brenKungStats = new DescriptiveStatistics();
        RandomGenerator rnd = new Well1024a();
        long modulus = SmallPrimes.nextPrime(10000);
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                hornerStats.clear();
                brenKungStats.clear();
            }
            MutablePolynomial poly = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            MutablePolynomial point = RandomPolynomials.randomMonicPoly(55 + rnd.nextInt(10), modulus, rnd);
            MutablePolynomial polyModulus = RandomPolynomials.randomMonicPoly(1 + rnd.nextInt(10), modulus, rnd);

//            System.out.println("====");
            long start = System.nanoTime();
            MutablePolynomial horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus);
            hornerStats.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            MutablePolynomial brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus);
            brenKungStats.addValue(System.nanoTime() - start);

            assertEquals(horner, brentKung);
        }

        System.out.println("Horner");
        System.out.println(hornerStats);

        System.out.println("BrentKung");
        System.out.println(brenKungStats);
    }

    static void assertComposition(MutablePolynomial poly, MutablePolynomial point, MutablePolynomial polyModulus, long modulus) {
        MutablePolynomial horner = compositionHorner(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus);
        MutablePolynomial brentKung = compositionBrentKung(poly, point, polyModulus, fastDivisionPreConditioning(polyModulus, modulus), modulus);
        MutablePolynomial plain = plainComposition(poly, point, polyModulus, modulus);
        assertEquals(plain, horner);
        assertEquals(plain, brentKung);
    }

    private static MutablePolynomial plainComposition(MutablePolynomial poly, MutablePolynomial point, MutablePolynomial polyModulus, long modulus) {
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = 0; i <= poly.degree; i++)
            res = polyMod(res.add(polyPowMod(point, i, polyModulus, modulus, true).multiply(poly.data[i], modulus), modulus), polyModulus, modulus, false);
        return res;
    }
}