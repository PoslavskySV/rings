package cc.r2.core.number;

import org.junit.Test;

import java.util.ArrayList;
import java.util.Random;

/**
 * Created by poslavsky on 11/11/2016.
 */
public class PrimesTest {
    @Test
    public void name1() throws Exception {
        BigInteger a = new BigInteger("37545145100657");
        ArrayList<BigInteger> factors = new ArrayList<>();
        Primes.factor(a, factors, new Random());
        System.out.println(factors);

        System.out.println(a.isPrime());
    }
}