package cc.r2.core.poly.univar2;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class PolynomialArithmeticsTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        long modulus = 5;
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(1, 4).modulus(modulus);
        lUnivariatePolynomialZp b = lUnivariatePolynomialZ.create(0, 2, 3).modulus(modulus);
        lUnivariatePolynomialZp polyModulus = lUnivariatePolynomialZ.create(0, 4, 0, 1).modulus(modulus);

        Assert.assertEquals(lUnivariatePolynomialZ.create(0, 4, 1).modulus(modulus),
                PolynomialArithmetics.polyMultiplyMod(a, b, polyModulus, true));
        Assert.assertEquals(lUnivariatePolynomialZ.create(0, 4, 1).modulus(modulus),
                PolynomialArithmetics.polyMultiplyMod(a.clone(), b, polyModulus, false));
    }

    @Test
    public void test3() throws Exception {
        Assert.assertEquals(lUnivariatePolynomialZ.create(1, 2, 1), PolynomialArithmetics.polyPow(lUnivariatePolynomialZ.create(1, 1), 2, true));
        Assert.assertEquals(lUnivariatePolynomialZ.create(1, 2, 1), PolynomialArithmetics.polyPow(lUnivariatePolynomialZ.create(1, 1), 2, false));
    }

    @Test
    public void test4() throws Exception {
        lUnivariatePolynomialZ a = lUnivariatePolynomialZ.create(1, 0, 1, 0, 1);
        lUnivariatePolynomialZ b = lUnivariatePolynomialZ.create(1, 1, 1);
        PolynomialArithmetics.polyPow(b.modulus(2), 2, true);
    }

    @Test
    public void test5() throws Exception {
        long modulus = 3;
        lUnivariatePolynomialZp a = lUnivariatePolynomialZ.create(0, 0, 0, 1).modulus(modulus);
        lUnivariatePolynomialZp polyModulus = lUnivariatePolynomialZ.create(0, -1, -1, -1, 0, 1, -1, 1, 1).modulus(modulus);
        Assert.assertEquals(lUnivariatePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), PolynomialArithmetics.polyPowMod(a, modulus, polyModulus, true));
        Assert.assertEquals(lUnivariatePolynomialZ.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), PolynomialArithmetics.polyPowMod(a, modulus, polyModulus, false));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(100, 1000); i++) {
            lUnivariatePolynomialZ poly = RandomPolynomials.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            lUnivariatePolynomialZ polyModulus = RandomPolynomials.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : getModulusArray(9, 1, 40)) {
                lUnivariatePolynomialZp base = poly.modulus(prime).monic();
                lUnivariatePolynomialZp modulus = polyModulus.modulus(prime).monic();
                Assert.assertEquals(PolynomialArithmetics.polyMod(PolynomialArithmetics.polyPow(base, exponent, true), modulus, false), PolynomialArithmetics.polyPowMod(base, exponent, modulus, true));
            }
        }
    }

    @Test
    public void test7() throws Exception {
        for (long modulus : new long[]{13, SmallPrimes.nextPrime(1 << 10), SmallPrimes.nextPrime(1 << 13), BigPrimes.nextPrime(1L << 43)}) {
            lUnivariatePolynomialZp polyModulus = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
            DivisionWithRemainder.InverseModMonomial invMod = DivisionWithRemainder.fastDivisionPreConditioning(polyModulus);
            for (int exp = 0; exp <= 2500; exp++) {
                Assert.assertEquals(PolynomialArithmetics.smallMonomial(exp, polyModulus, invMod), PolynomialArithmetics.createMonomialMod(exp, polyModulus, invMod));
            }
        }
    }
}