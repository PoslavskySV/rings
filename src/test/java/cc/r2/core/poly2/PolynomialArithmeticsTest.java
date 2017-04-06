package cc.r2.core.poly2;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.PolynomialArithmetics.*;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class PolynomialArithmeticsTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        long modulus = 5;
        lMutablePolynomialZp a = lMutablePolynomialZ.create(1, 4).modulus(modulus);
        lMutablePolynomialZp b = lMutablePolynomialZ.create(0, 2, 3).modulus(modulus);
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(0, 4, 0, 1).modulus(modulus);

        assertEquals(lMutablePolynomialZ.create(0, 4, 1).modulus(modulus),
                polyMultiplyMod(a, b, polyModulus, true));
        assertEquals(lMutablePolynomialZ.create(0, 4, 1).modulus(modulus),
                polyMultiplyMod(a.clone(), b, polyModulus, false));
    }

    @Test
    public void test3() throws Exception {
        assertEquals(lMutablePolynomialZ.create(1, 2, 1), polyPow(lMutablePolynomialZ.create(1, 1), 2, true));
        assertEquals(lMutablePolynomialZ.create(1, 2, 1), polyPow(lMutablePolynomialZ.create(1, 1), 2, false));
    }

    @Test
    public void test4() throws Exception {
        lMutablePolynomialZ a = lMutablePolynomialZ.create(1, 0, 1, 0, 1);
        lMutablePolynomialZ b = lMutablePolynomialZ.create(1, 1, 1);
        polyPow(b.modulus(2), 2, true);
    }

    @Test
    public void test5() throws Exception {
        long modulus = 3;
        lMutablePolynomialZp a = lMutablePolynomialZ.create(0, 0, 0, 1).modulus(modulus);
        lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(0, -1, -1, -1, 0, 1, -1, 1, 1).modulus(modulus);
        assertEquals(lMutablePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), polyPowMod(a, modulus, polyModulus, true));
        assertEquals(lMutablePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), polyPowMod(a, modulus, polyModulus, false));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(100, 1000); i++) {
            lMutablePolynomialZ poly = RandomPolynomials.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            lMutablePolynomialZ polyModulus = RandomPolynomials.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : getModulusArray(9, 1, 40)) {
                lMutablePolynomialZp base = poly.modulus(prime).monic();
                lMutablePolynomialZp modulus = polyModulus.modulus(prime).monic();
                assertEquals(polyMod(polyPow(base, exponent, true), modulus, false), polyPowMod(base, exponent, modulus, true));
            }
        }
    }

    @Test
    public void test7() throws Exception {
        for (long modulus : new long[]{13, SmallPrimes.nextPrime(1 << 10), SmallPrimes.nextPrime(1 << 13), BigPrimes.nextPrime(1L << 43)}) {
            lMutablePolynomialZp polyModulus = lMutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
            DivisionWithRemainder.InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
            for (int exp = 0; exp <= 2500; exp++) {
                assertEquals(smallMonomial(exp, polyModulus, invMod), createMonomialMod(exp, polyModulus, invMod));
            }
        }
    }
}