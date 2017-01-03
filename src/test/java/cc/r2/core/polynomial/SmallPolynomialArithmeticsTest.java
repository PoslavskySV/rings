package cc.r2.core.polynomial;

import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.junit.Test;

import static cc.r2.core.polynomial.SmallPolynomialArithmetics.*;
import static org.junit.Assert.assertEquals;

public class SmallPolynomialArithmeticsTest {
    @Test
    public void test1() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1, 4);
        MutableLongPoly b = MutableLongPoly.create(0, 2, 3);
        MutableLongPoly polyModulus = MutableLongPoly.create(0, 4, 0, 1);
        long modulus = 5;
        assertEquals(MutableLongPoly.create(0, 4, 1), multiplyMod(a, b, polyModulus, modulus));
    }

    @Test
    public void test2() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(1, 2, 3, 4, 0);
        assertEquals(MutableLongPoly.create(2, 6, 12), derivative(a));
    }

    @Test
    public void test3() throws Exception {
        assertEquals(MutableLongPoly.create(1, 2, 1), pow(MutableLongPoly.create(1, 1), 2));
    }

    @Test
    public void test4() throws Exception {
//        System.out.println(LongArithmetics.modInverse(7, 7));
//        assertEquals(MutableLongPoly.create(1, 2, 1), pow(MutableLongPoly.create(1, 1), 2));
        MutableLongPoly a = MutableLongPoly.create(1, 0, 1, 0, 1);
        MutableLongPoly b = MutableLongPoly.create(1, 1, 1);
        System.out.println(SmallPolynomialArithmetics.powMod(b, 2, 2));
    }

    @Test
    public void test5() throws Exception {
        MutableLongPoly a = MutableLongPoly.create(0, 0, 0, 1);
        MutableLongPoly polyModulus = MutableLongPoly.create(0, -1, -1, -1, 0, 1, -1, 1, 1);
        long modulus = 3;
        assertEquals(MutableLongPoly.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), powMod(a, modulus, polyModulus, modulus));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = new Well1024a();
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        long[] primes = {2, 3, 5, 7, 11, 17, 67, 29, 31, 89, 101, 107, 139, 223};
        for (int i = 0; i < 100; i++) {
            MutableLongPoly poly = SmallPolynomialsTest.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            MutableLongPoly polyModulus = SmallPolynomialsTest.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : primes) {
                MutableLongPoly base = poly.clone().monic(prime);
                MutableLongPoly modulus = polyModulus.clone().monic(prime);
                assertEquals(mod(powMod(base, exponent, prime), modulus, prime), powMod(base, exponent, polyModulus, prime));
            }
        }
    }

    @Test
    public void name() throws Exception {
        MutableLongPoly r = MutableLongPoly.create(0,0,1);
        System.out.println(powMod(r, 23, 23));
        System.out.println(powMod0(r, 23));

    }
}