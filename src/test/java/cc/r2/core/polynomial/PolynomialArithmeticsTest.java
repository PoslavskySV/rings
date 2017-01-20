package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import static cc.r2.core.polynomial.PolynomialArithmetics.*;
import static org.junit.Assert.assertEquals;

public class PolynomialArithmeticsTest {
    @Test
    public void test1() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 4);
        MutablePolynomial b = MutablePolynomial.create(0, 2, 3);
        MutablePolynomial polyModulus = MutablePolynomial.create(0, 4, 0, 1);
        long modulus = 5;
        assertEquals(MutablePolynomial.create(0, 4, 1), polyMultiplyMod(a, b, polyModulus, modulus, true));
        assertEquals(MutablePolynomial.create(0, 4, 1), polyMultiplyMod(a.clone(), b, polyModulus, modulus, false));
    }

    @Test
    public void test3() throws Exception {
        assertEquals(MutablePolynomial.create(1, 2, 1), polyPow(MutablePolynomial.create(1, 1), 2, true));
        assertEquals(MutablePolynomial.create(1, 2, 1), polyPow(MutablePolynomial.create(1, 1), 2, false));
    }

    @Test
    public void test4() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(1, 0, 1, 0, 1);
        MutablePolynomial b = MutablePolynomial.create(1, 1, 1);
        PolynomialArithmetics.polyPowMod(b, 2, 2, true);
    }

    @Test
    public void test5() throws Exception {
        MutablePolynomial a = MutablePolynomial.create(0, 0, 0, 1);
        MutablePolynomial polyModulus = MutablePolynomial.create(0, -1, -1, -1, 0, 1, -1, 1, 1);
        long modulus = 3;
        assertEquals(MutablePolynomial.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), PolynomialArithmetics.polyPowMod(a, modulus, polyModulus, modulus, true));
        assertEquals(MutablePolynomial.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), PolynomialArithmetics.polyPowMod(a, modulus, polyModulus, modulus, false));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        long[] primes = {2, 3, 5, 7, 11, 17, 67, 29, 31, 89, 101, 107, 139, 223};
        for (int i = 0; i < 100; i++) {
            MutablePolynomial poly = RandomPolynomials.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            MutablePolynomial polyModulus = RandomPolynomials.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : primes) {
                MutablePolynomial base = poly.clone().monic(prime);
                MutablePolynomial modulus = polyModulus.clone().monic(prime);
                assertEquals(polyMod(polyPowMod(base, exponent, prime, true), modulus, prime, false), PolynomialArithmetics.polyPowMod(base, exponent, polyModulus, prime, true));
            }
        }
    }
}