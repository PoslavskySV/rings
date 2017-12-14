package cc.redberry.rings;

import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

/**
 * @since 2.3
 */
public class RationalReconstructionTest extends AbstractTest {
    @Ignore
    @Test
    public void test2() throws Exception {
        UnivariatePolynomialZp64 num = UnivariatePolynomialZp64.create(17, new long[]{1, 2, 3, 4, 5, 6});
        UnivariatePolynomialZp64 den = UnivariatePolynomialZp64.create(17, new long[]{1, 2, 3, 4, 5, 6, 8, 9, 8, 7, 6, 5, 4, 3, 2, 1});
        UnivariatePolynomialZp64 mod = IrreduciblePolynomials.randomIrreduciblePolynomial(17, 1 + num.degree() * den.degree() * 2, getRandom());
        Ring<UnivariatePolynomialZp64> r = Rings.GF(mod);

        UnivariatePolynomialZp64 n = r.multiply(num, r.reciprocal(den));

        UnivariatePolynomialZp64[] rec = RationalReconstruction.reconstruct(n, mod, 20, 20);
        System.out.println(Arrays.toString(rec));
        // to within a sign
        //Assert.assertArrayEquals(new UnivariatePolynomialZp64[]{num, den}, rec);
    }

    @Test
    public void testReconstructionRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < 1000; i++) {
            long num = rndd.nextLong(10, 1000);
            long den = rndd.nextLong(num + 1, 1000);
            long gcd = MachineArithmetic.gcd(num, den);
            num = num / gcd;
            den = den / gcd;
            long modulus = BigPrimes.nextPrime(2 * num * den);
            //if (rnd.nextBoolean())
            //  num = -num;
            long n = MachineArithmetic.mod(num * MachineArithmetic.modInverse(den, modulus), modulus);
            assert num != 0 && den != 0 && n != 0;
            long[] re = RationalReconstruction.reconstruct(n, modulus, Math.abs(num) + 10, Math.abs(den) + 10);
            Assert.assertEquals(n, MachineArithmetic.mod(re[0] * MachineArithmetic.modInverse(re[1], modulus), modulus));
        }
    }
}