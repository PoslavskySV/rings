package cc.r2.core.number.qsi.factpor;


import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.ArrayList;
import java.util.List;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.number.qsi.factpor.QuadraticSieve.sqrtBigInt;

public final class Primes {
    private Primes() {
    }

    public static boolean isPrime(int j) {
        long k = (long) Math.sqrt(j);
        if (j % 2 == 0)
            return false;
        for (long p = 3; p <= k; p += 2)
            if (j % p == 0)
                return false;
        return true;
    }

    /**
     * Plain trial division algorithm
     *
     * @param integer integer to factor
     * @return list of factors
     */
    public static List<BigInteger> TrialDivision(BigInteger integer, BigInteger upperBound) {
        if (integer.isZero())
            return new ArrayList<>();
        if (integer.isOne()) {
            ArrayList<BigInteger> f = new ArrayList<>();
            f.add(integer);
            return f;
        }
        ArrayList<BigInteger> factors = new ArrayList<>();
        BigInteger temp = integer;

        //rule out all factors = 2 (reduce integer to odd number)
        BigInteger mod2 = temp.mod(TWO);
        while (mod2.isZero()) {
            factors.add(TWO);
            temp = temp.divide(TWO);
            mod2 = temp.mod(TWO);
        }
        integer = temp;

        //main loop over all odd multipliers
        BigInteger i = THREE;
        for (; i.compareTo(upperBound) <= 0
                && temp.compareTo(i) >= 0
                && integer.compareTo(i.multiply(i)) >= 0; i = i.add(TWO)) {
            BigInteger[] qr = temp.divideAndRemainder(i);
            if (qr[1].isZero()) {
                factors.add(i);
                temp = qr[0];
                i = i.subtract(TWO);
            }
        }
        if (temp.compareTo(ONE) > 0)
            factors.add(temp);
        if (factors.isEmpty())
            factors.add(integer);
        return factors;
    }

    /**
     * Fermat's factoring algorithm works like trial division, but walks in the opposite
     * direction. Thus, it can be used to factor a number that we know has a factor in the
     * interval [Sqrt(n) - upperBound, Sqrt(n) + upperBound].
     *
     * @param n          number to factor
     * @param upperBound upper bound
     * @return a single factor
     */
    public static BigInteger fermat(BigInteger n, long upperBound) {
        long cnt = 0;
        BigInteger x = sqrtBigInt(n).add(ONE);
        BigInteger u = x.multiply(TWO).add(ONE);
        BigInteger v = ONE;
        BigInteger r = x.multiply(x).subtract(n);
        while (!r.isZero()) {
            cnt++;
            if (cnt > upperBound)
                return ZERO;
            while (r.compareTo(ZERO) > 0) {
                r = r.subtract(v);
                v = v.add(TWO);
            }
            if (r.compareTo(ZERO) < 0) {
                r = r.add(u);
                u = u.add(TWO);
            }
        }
        return u.subtract(v).divide(TWO);
    }

    /**
     * Pollards's rho algorithm (random search version).
     *
     * @param n        integer to factor
     * @param attempts number of random attempts
     * @return a single factor of {@code n} or null if no factors found
     */
    public static BigInteger PollardRho(BigInteger n, int attempts, RandomGenerator rn) {
        BigInteger divisor;
        BigInteger c = new BigInteger(n.bitLength(), rn);
        BigInteger x = new BigInteger(n.bitLength(), rn);
        BigInteger xx = x;

        // check divisibility by 2
        if (n.mod(TWO).isZero()) return TWO;

        do {
            x = x.multiply(x).mod(n).add(c).mod(n);
            xx = xx.multiply(xx).mod(n).add(c).mod(n);
            xx = xx.multiply(xx).mod(n).add(c).mod(n);
            divisor = x.subtract(xx).gcd(n);
        } while (attempts-- > 0 && divisor.isOne());

        return divisor.isOne() ? null : divisor;
    }

    /**
     * Pollards's rho algorithm.
     *
     * @param n          integer to factor
     * @param upperBound expected B-smoothness
     * @return a single factor of {@code n} or null if no factors found
     */
    public static BigInteger PollardRho(BigInteger n, long upperBound) {
        long range = 1;
        long terms = 0;
        BigInteger x1 = TWO;
        BigInteger x2 = FIVE;
        BigInteger product = ONE;
        while (terms <= upperBound) {
            for (long j = 1; j <= range; j++) {
                x2 = x2.multiply(x2).add(ONE).mod(n);
                product = product.multiply(x1.subtract(x2)).mod(n);
                if (terms++ > upperBound)
                    break;
                if (terms % 5 == 0) {
                    BigInteger g = n.gcd(product);
                    if (g.compareTo(ONE) > 0) {
                        return g;
                    }
                    product = ONE;
                }
            }
            x1 = x2;
            range *= 2;
            for (long j = 1; j <= range; j++)
                x2 = x2.multiply(x2).add(ONE).mod(n);
        }
        return null;
    }


    /**
     * Pollards's p-1 algorithm.
     *
     * @param n          integer to factor
     * @param upperBound expected B-smoothness
     * @return a single factor of {@code n} or null if no factors found
     */
    public static BigInteger PollardP1(BigInteger n, long upperBound) {
        BigInteger g, i, m;
        for (int outerCnt = 0; outerCnt < 5; outerCnt++) {
            switch (outerCnt) {
                case 0:
                    m = TWO;
                    break;
                case 1:
                    m = THREE;
                    break;
                case 2:
                    m = FOUR;
                    break;
                case 3:
                    m = FIVE;
                    break;
                case 4:
                    m = SEVEN;
                    break;
                default:
                    m = TWO;
            }
            i = ONE;
            for (long cnt = 2; cnt <= upperBound; cnt++) {
                i = i.add(ONE);
                m = m.modPow(i, n);
                if (cnt % 5 == 0) {
                    g = n.gcd(m.subtract(ONE));
                    if ((g.compareTo(ONE) > 0) && (g.compareTo(n) < 0)) {
                        return g;
                    }
                }
            }
        }
        return null;
    }

    public static BigInteger QuadraticSieve(BigInteger n, int bound) {
        return new QuadraticSieve(n).quadraticSieve(bound);
    }
}
