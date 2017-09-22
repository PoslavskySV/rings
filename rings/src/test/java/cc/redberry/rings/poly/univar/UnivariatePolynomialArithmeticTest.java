package cc.redberry.rings.poly.univar;

import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.primes.SmallPrimes;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class UnivariatePolynomialArithmeticTest extends AUnivariateTest {
    @Test
    public void test1() throws Exception {
        long modulus = 5;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(1, 4).modulus(modulus);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(0, 2, 3).modulus(modulus);
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(0, 4, 0, 1).modulus(modulus);

        Assert.assertEquals(UnivariatePolynomialZ64.create(0, 4, 1).modulus(modulus),
                UnivariatePolynomialArithmetic.polyMultiplyMod(a, b, polyModulus, true));
        Assert.assertEquals(UnivariatePolynomialZ64.create(0, 4, 1).modulus(modulus),
                UnivariatePolynomialArithmetic.polyMultiplyMod(a.clone(), b, polyModulus, false));
    }

    @Test
    public void test3() throws Exception {
        Assert.assertEquals(UnivariatePolynomialZ64.create(1, 2, 1), UnivariatePolynomialArithmetic.polyPow(UnivariatePolynomialZ64.create(1, 1), 2, true));
        Assert.assertEquals(UnivariatePolynomialZ64.create(1, 2, 1), UnivariatePolynomialArithmetic.polyPow(UnivariatePolynomialZ64.create(1, 1), 2, false));
    }

    @Test
    public void test4() throws Exception {
        UnivariatePolynomialZ64 a = UnivariatePolynomialZ64.create(1, 0, 1, 0, 1);
        UnivariatePolynomialZ64 b = UnivariatePolynomialZ64.create(1, 1, 1);
        UnivariatePolynomialArithmetic.polyPow(b.modulus(2), 2, true);
    }

    @Test
    public void test5() throws Exception {
        long modulus = 3;
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(0, 0, 0, 1).modulus(modulus);
        UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(0, -1, -1, -1, 0, 1, -1, 1, 1).modulus(modulus);
        Assert.assertEquals(UnivariatePolynomialZ64.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), UnivariatePolynomialArithmetic.polyPowMod(a, modulus, polyModulus, true));
        Assert.assertEquals(UnivariatePolynomialZ64.create(0, -1, 0, 0, 1, 1, 1, -1).modulus(modulus), UnivariatePolynomialArithmetic.polyPowMod(a, modulus, polyModulus, false));
    }

    @Test
    public void test6() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomialZ64 poly = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(1, 5), 100, rnd);
            UnivariatePolynomialZ64 polyModulus = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(poly.degree == 1 ? 0 : 1, poly.degree), 100, rnd);
            poly.data[poly.degree] = 1;
            polyModulus.data[polyModulus.degree] = 1;
            int exponent = 2 + rnd.nextInt(20);
            for (long prime : getModulusArray(9, 1, 40)) {
                UnivariatePolynomialZp64 base = poly.modulus(prime).monic();
                UnivariatePolynomialZp64 modulus = polyModulus.modulus(prime).monic();
                Assert.assertEquals(UnivariatePolynomialArithmetic.polyMod(UnivariatePolynomialArithmetic.polyPow(base, exponent, true), modulus, false), UnivariatePolynomialArithmetic.polyPowMod(base, exponent, modulus, true));
            }
        }
    }

    @Test
    public void test7() throws Exception {
        for (long modulus : new long[]{13, SmallPrimes.nextPrime(1 << 10), SmallPrimes.nextPrime(1 << 13), BigPrimes.nextPrime(1L << 43)}) {
            UnivariatePolynomialZp64 polyModulus = UnivariatePolynomialZ64.create(1, 2, 3, 4, 5, 6, 1).modulus(modulus);
            UnivariateDivision.InverseModMonomial invMod = UnivariateDivision.fastDivisionPreConditioning(polyModulus);
            for (int exp = 0; exp <= 2500; exp++) {
                Assert.assertEquals(UnivariatePolynomialArithmetic.smallMonomial(exp, polyModulus, invMod), UnivariatePolynomialArithmetic.createMonomialMod(exp, polyModulus, invMod));
            }
        }
    }
}