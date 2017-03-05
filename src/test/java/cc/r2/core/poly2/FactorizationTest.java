package cc.r2.core.poly2;

import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly2.Factorization.factorBigPrime;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class FactorizationTest {

    @Test
    public void test1() throws Exception {
        long modulus = 127;
        MutablePolynomialMod poly =
                MutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus)
                        .multiply(MutablePolynomialZ.create(1, 2, 3, 4, 5).modulus(modulus))
                        .square()
                        .multiply(MutablePolynomialZ.create(2, 3, 4).modulus(modulus))
                        .multiply(MutablePolynomialZ.create(22, 32, 42).modulus(modulus))
                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 34).modulus(modulus))
                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 31).modulus(modulus))
                        .multiply(MutablePolynomialZ.create(22, -32, 42, 12, 23, 31).modulus(modulus));
        poly = poly.multiply(poly.shiftLeft(5)).multiply(poly.shiftRight(5));

        FactorDecomposition<MutablePolynomialMod> factorization = Factorization.factor(poly);
        System.out.println(factorization);
        System.out.println(poly);
        System.out.println(factorization.toPolynomial(poly));

    }


    @Test
    public void name() throws Exception {
//        System.out.println(factorBigPrime(MutablePolynomialZ.create(-1, 0, 1)));
//        System.out.println(factorBigPrime(MutablePolynomialZ.create(-4, 0, 1)));

        MutablePolynomialZ a = MutablePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
        MutablePolynomialZ poly = a.multiply(b).primitivePart();
        Assert.assertTrue(SquareFreeFactorization.isSquareFree(poly));

//        MutablePolynomialMod moduloImage = poly.modulus(6101).monic();
//        System.out.println(moduloImage);
////        System.out.println(poly);
//        System.out.println(factor(moduloImage));
//        FactorDecomposition<MutablePolynomialMod> modFact = factor(moduloImage.monic());
//        System.out.println(modFact);

//        System.out.println(modFact.get(0).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(0).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(2).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(2).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
//
//        System.out.println(modFact);

        poly = MutablePolynomialZ.create(46225, 0, -5596840, 0, 13950764, 0, -7453176, 0, 1513334, 0, -141912, 0, 6476, 0, -136, 0, 1);
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            System.out.println(factorBigPrime(poly));
            System.out.println(System.nanoTime() - start);
        }

    }
}