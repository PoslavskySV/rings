package cc.r2.core.poly2;

import cc.r2.core.poly2.Factorization.HenselData;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Collections;
import java.util.List;

import static cc.r2.core.poly2.Factorization.factor;
import static cc.r2.core.poly2.Factorization.liftFactorization;
import static cc.r2.core.poly2.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreePart;
import static cc.r2.core.poly2.SquareFreeFactorization.isSquareFree;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class FactorizationTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        Assert.assertTrue(factor(MutablePolynomialZ.create(3, 7).modulus(19)).get(0).isMonic());
    }

    static void assertHenselLift(HenselData lift) {
        Assert.assertEquals(lift.base.modulus(lift.modulus), lift.aFactor.clone().multiply(lift.bFactor));
        Assert.assertEquals(lift.base.modulus(lift.modulus).createOne(),
                lift.aFactor.clone().multiply(lift.aCoFactor)
                        .add(lift.bFactor.clone().multiply(lift.bCoFactor)));
    }

    void testHenselStepRandomModulus(
            MutablePolynomialZ poly,
            MutablePolynomialZ aFactor,
            MutablePolynomialZ bFactor,
            int nPrimes) {

        for (long modulus : getModulusArray(nPrimes, 0, 5, 0)) {
            MutablePolynomialMod aMod = aFactor.modulus(modulus);
            MutablePolynomialMod bMod = bFactor.modulus(modulus);
            if (!PolynomialGCD.PolynomialGCD(aMod, bMod).isConstant())
                continue;
            HenselData hensel = Factorization.createHenselInput(modulus,
                    poly, aMod, bMod);
            assertHenselLift(hensel);

            for (int i = 0; i < 10; i++) {
                try {
                    hensel = hensel.liftQuadratic();
                } catch (ArithmeticException ex) {
                    if (!"long overflow".equals(ex.getMessage()))
                        throw ex;
                    break;
                }
                assertHenselLift(hensel);
            }
        }
    }

    @Test
    public void testHensel1() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(-2, -1, 2, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(-2, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, 20);
    }

    @Test
    public void testHensel2() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(3, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(4, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, 20);
    }

    static void testMultiFactorHenselLifting(MutablePolynomialZ base, long modulus, int nIterations) {
        testMultiFactorHenselLifting(base, modulus, nIterations, true);
    }


    static void testMultiFactorHenselLifting(MutablePolynomialZ base, long modulus, int nIterations, boolean rethrow) {
        base = SquareFreePart(base);

        MutablePolynomialMod baseMod = base.modulus(modulus);
        if (!isSquareFree(baseMod))
            return;

        List<MutablePolynomialMod> modularFactors = factor(baseMod).factors;
        assertFactorization(baseMod, baseMod.lc(), modularFactors);

        long liftedModulus = modulus;
        for (int i = 1; i < nIterations; i++) {
            try {
                List<MutablePolynomialMod> tmp = liftFactorization(modulus, i, base, modularFactors);
                MutablePolynomialMod bm = base.modulus(liftedModulus = liftedModulus * liftedModulus);
                assertFactorization(bm, bm.lc(), tmp);
            } catch (ArithmeticException ex) {
                if (rethrow || !"long overflow".equals(ex.getMessage()))
                    throw ex;
                break;
            }
        }
    }

    @Test
    public void testHensel3() throws Exception {
        MutablePolynomialZ base = MutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 5, 5);
    }

    @Test
    public void testHensel4() throws Exception {
        MutablePolynomialZ base = MutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 11, 5);
    }

    @Test
    public void testHensel5_random() throws Exception {
        RandomDataGenerator rnd = getRandomData();
        for (int i = 0; i < its(10000, 1000); ++i) {
            MutablePolynomialZ poly = RandomPolynomials.randomPoly(rnd.nextInt(1, 20), 1000, rnd.getRandomGenerator());
            for (long modulus : getModulusArray((int) its(5, 15), 0, 5, 0)) {
                MutablePolynomialMod polyMod = poly.modulus(modulus);
                if (polyMod.isConstant() || polyMod.isMonomial())
                    continue;

                while (poly.lc() % modulus == 0)
                    poly.data[poly.degree] = rnd.nextInt(0, 100);

                try {
                    testMultiFactorHenselLifting(poly, modulus, 10, false);
                } catch (AssertionError err) {
                    System.out.println(poly.toStringForCopy());
                    System.out.println(modulus);
                    throw err;
                }
            }
        }
    }


    @Test
    public void testHensel5() throws Exception {
        MutablePolynomialZ base = MutablePolynomialZ.create(69, 30);
        long modulus = 29;

        MutablePolynomialMod baseMod = base.modulus(modulus);
        System.out.println(factor(baseMod));
        List<MutablePolynomialMod> fact = liftFactorization(modulus, 1, base, Collections.singletonList(baseMod.monic()));
        MutablePolynomialMod sin = fact.get(0);
        System.out.println(sin.multiply(base.lc()));


        System.out.println(factor(base.modulus(modulus)));
        testMultiFactorHenselLifting(base, modulus, 1);

    }

//
//    @Test
//    public void test1() throws Exception {
//        long modulus = 127;
//        MutablePolynomialMod poly =
//                MutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus)
//                        .multiply(MutablePolynomialZ.create(1, 2, 3, 4, 5).modulus(modulus))
//                        .square()
//                        .multiply(MutablePolynomialZ.create(2, 3, 4).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 34).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 31).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, -32, 42, 12, 23, 31).modulus(modulus));
//        poly = poly.multiply(poly.shiftLeft(5)).multiply(poly.shiftRight(5));
//
//        FactorDecomposition<MutablePolynomialMod> factorization = Factorization.factor(poly);
//        System.out.println(factorization);
//        System.out.println(poly);
//        System.out.println(factorization.toPolynomial(poly));
//
//    }
//
//
//    @Test
//    public void name() throws Exception {
////        System.out.println(factorBigPrime(MutablePolynomialZ.create(-1, 0, 1)));
////        System.out.println(factorBigPrime(MutablePolynomialZ.create(-4, 0, 1)));
//
//        MutablePolynomialZ a = MutablePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
//        MutablePolynomialZ b = MutablePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
//        MutablePolynomialZ poly = a.multiply(b).primitivePart();
//        Assert.assertTrue(SquareFreeFactorization.isSquareFree(poly));
//
////        MutablePolynomialMod moduloImage = poly.modulus(6101).monic();
////        System.out.println(moduloImage);
//////        System.out.println(poly);
////        System.out.println(factor(moduloImage));
////        FactorDecomposition<MutablePolynomialMod> modFact = factor(moduloImage.monic());
////        System.out.println(modFact);
//
////        System.out.println(modFact.get(0).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
////        System.out.println(modFact.get(0).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
////        System.out.println(modFact.get(2).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
////        System.out.println(modFact.get(2).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
////
////        System.out.println(modFact);
//
//        poly = MutablePolynomialZ.create(46225, 0, -5596840, 0, 13950764, 0, -7453176, 0, 1513334, 0, -141912, 0, 6476, 0, -136, 0, 1);
//        for (int i = 0; i < 1000; i++) {
//            long start = System.nanoTime();
//            System.out.println(factorBigPrime(poly));
//            System.out.println(System.nanoTime() - start);
//        }
//
//    }
}