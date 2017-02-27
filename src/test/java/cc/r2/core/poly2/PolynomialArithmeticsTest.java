package cc.r2.core.poly2;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import static cc.r2.core.poly2.DivisionWithRemainder.fastDivisionPreConditioning;
import static cc.r2.core.poly2.PolynomialArithmetics.*;
import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class PolynomialArithmeticsTest {
    @Test
    public void test1() throws Exception {
        long modulus = 5;
        MutablePolynomialMod a = MutablePolynomialZ.create(1, 4).modulus(modulus);
        MutablePolynomialMod b = MutablePolynomialZ.create(0, 2, 3).modulus(modulus);
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(0, 4, 0, 1).modulus(modulus);

        assertEquals(MutablePolynomialZ.create(0, 4, 1).modulus(modulus),
                polyMultiplyMod(a, b, polyModulus, true));
        assertEquals(MutablePolynomialZ.create(0, 4, 1).modulus(modulus),
                polyMultiplyMod(a.clone(), b, polyModulus, false));
    }

    @Test
    public void test3() throws Exception {
        assertEquals(MutablePolynomialZ.create(1, 2, 1), polyPow(MutablePolynomialZ.create(1, 1), 2, true));
        assertEquals(MutablePolynomialZ.create(1, 2, 1), polyPow(MutablePolynomialZ.create(1, 1), 2, false));
    }

    @Test
    public void test4() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(1, 0, 1, 0, 1);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 1, 1);
        polyPow(b.modulus(2), 2, true);
    }

    @Test
    public void test5() throws Exception {
        long modulus = 3;
        MutablePolynomialMod a = MutablePolynomialZ.create(0, 0, 0, 1).modulus(modulus);
        MutablePolynomialMod polyModulus = MutablePolynomialZ.create(0, -1, -1, -1, 0, 1, -1, 1, 1).modulus(modulus);
        assertEquals(MutablePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), polyPowMod(a, modulus, polyModulus, true));
        assertEquals(MutablePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), polyPowMod(a, modulus, polyModulus, false));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        long[] primes = {2, 3, 5, 7, 11, 17, 67, 29, 31, 89, 101, 107, 139, 223};
        for (int i = 0; i < 100; i++) {
            MutablePolynomialZ poly = RandomPolynomials.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            MutablePolynomialZ polyModulus = RandomPolynomials.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : primes) {
                MutablePolynomialMod base = poly.modulus(prime).monic();
                MutablePolynomialMod modulus = polyModulus.modulus(prime).monic();
                assertEquals(polyMod(polyPow(base, exponent, true), modulus, false), polyPowMod(base, exponent, modulus, true));
            }
        }
    }

    @Test
    public void test7() throws Exception {
        for (long modulus : new long[]{13, SmallPrimes.nextPrime(1 << 10), SmallPrimes.nextPrime(1 << 13), BigPrimes.nextPrime(1L << 43)}) {
            MutablePolynomialMod polyModulus = MutablePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
            DivisionWithRemainder.InverseModMonomial invMod = fastDivisionPreConditioning(polyModulus);
            for (int exp = 0; exp <= 2500; exp++) {
                assertEquals(smallMonomial(exp, polyModulus, invMod), createMonomialMod(exp, polyModulus, invMod));
            }
        }
    }
}