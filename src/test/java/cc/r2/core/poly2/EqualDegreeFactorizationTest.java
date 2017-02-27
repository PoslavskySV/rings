package cc.r2.core.poly2;

import cc.r2.core.number.primes.BigPrimes;
import org.junit.Test;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class EqualDegreeFactorizationTest {
    @Test
    public void test1() throws Exception {
        long modulus = BigPrimes.nextPrime(1L << 40);
        MutablePolynomialMod p =
                MutablePolynomialZ.create(171, 234, 31, 1241, 123435).modulus(modulus)
                        .multiply(MutablePolynomialZ.create(13435, 42143241, 13, 22342341, 19123).modulus(modulus));


        for (int i = 0; i < 10000; i++) {
            long start = System.nanoTime();
            System.out.println(EqualDegreeFactorization.CantorZassenhaus(p, 4).factors.size());
            System.out.println(System.nanoTime() - start);
            System.out.println("=====");

        }


    }
}