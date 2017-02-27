package cc.r2.core.poly2;

import org.junit.Test;

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

        for (int i = 0; i < 19999; i++) {
            long start = System.nanoTime();
            System.out.println(i + Factorization.factor(DistinctDegreeFactorizationTest.bigPoly).factors.size());
            System.out.println(System.nanoTime() - start);
            System.out.println("-====");
        }


    }
}