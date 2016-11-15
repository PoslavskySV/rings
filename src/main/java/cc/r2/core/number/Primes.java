package cc.r2.core.number;

import java.util.ArrayList;
import java.util.Random;

import static cc.r2.core.number.BigInteger.ONE;
import static cc.r2.core.number.BigInteger.TWO;

/**
 * Created by poslavsky on 11/11/2016.
 */
public class Primes {
    public static BigInteger rho(BigInteger N, Random rn) {
        BigInteger divisor;
        BigInteger c = new BigInteger(N.bitLength(), rn);
        BigInteger x = new BigInteger(N.bitLength(), rn);
        BigInteger xx = x;

        // check divisibility by 2
        if (N.mod(TWO).isZero()) return TWO;

        do {
            x = x.multiply(x).mod(N).add(c).mod(N);
            xx = xx.multiply(xx).mod(N).add(c).mod(N);
            xx = xx.multiply(xx).mod(N).add(c).mod(N);
            divisor = x.subtract(xx).gcd(N);
        } while ((divisor.compareTo(ONE)) == 0);

        return divisor;
    }

    public static void factor(BigInteger N, ArrayList<BigInteger> factors, Random random) {
        if (N.compareTo(ONE) == 0) return;
        if (N.isProbablePrime(20)) {
            factors.add(N);
            return;
        }
        BigInteger divisor = rho(N, random);
        factor(divisor, factors, random);
        factor(N.divide(divisor), factors, random);
    }

}
