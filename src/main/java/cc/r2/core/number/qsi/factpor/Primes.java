package cc.r2.core.number.qsi.factpor;


import cc.r2.core.number.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;

import java.util.ArrayList;
import java.util.List;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.number.qsi.factpor.QuadraticSieve.sqrtBigInt;
import static cc.r2.core.number.qsi.factpor.SieveOfAtkin.SmallPrimesSieve;

public final class Primes {
    private Primes() {
    }

    public static boolean isPrime(int j) {
        if (j == 2)
            return true;
        long k = (long) Math.sqrt(j);
        if (j % 2 == 0)
            return false;
        for (long p = 3; p <= k; p += 2)
            if (j % p == 0)
                return false;
        return true;
    }

    public static boolean isPrime(BigInteger b) {
        return findFactor(b).equals(b);
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

    private static final RandomGenerator PollardRnd = new Well1024a();


    static BigInteger findFactor(BigInteger n) {
        int numBits = n.bitCount();
        BigInteger r;

        //switching between algorithms

        if (numBits < 20) {
            // t = 1e4 - 3e6
            r = PollardRho(n, 131_072);
            if (r != null)
                return r;
        }

        if (numBits < 30) {
            // t = 5e6
            r = PollardRho(n, 1024, PollardRnd);
            if (r != null)
                return r;
            // t = 2e6 - 5e7
            r = PollardRho(n, 131_072);
            if (r != null)
                return r;
        }

        if (numBits < 60) {
            // t = 2e5
            r = PollardRho(n, 128);
            if (r != null)
                return r;

            // t = 5e5
            r = PollardRho(n, 128, PollardRnd);
            if (r != null)
                return r;

            // t = 1e6
            r = PollardP1(n, 128);
            if (r != null)
                return r;

            // t = 2e7
            r = PollardRho(n, 131_072);
            if (r != null)
                return r;
        }

        //<-really large number with large primes

        // t = 5e5
        r = PollardRho(n, 128);
        if (r != null)
            return r;

        // t = 5e5
        r = PollardP1(n, 128);
        if (r != null)
            return r;

        // t = 5e6
        r = PollardRho(n, 1032, PollardRnd);
        if (r != null)
            return r;

        // t = 1e8
        r = PollardRho(n, 131_072);
        if (r != null)
            return r;


        // t = 1e9
        r = PollardP1(n, 131_072);
        if (r != null)
            return r;

        // t = 1e9 -> oo
        // be sure that trial division is done
        assert n.compareTo(BigInteger.valueOf(32768)) > 0;
        r = QuadraticSieve(n, 32768);
        assert r != null;
        return r;
    }

    private static boolean checkKnownSmallPrime(BigInteger b) {
        return b.compareTo(SmallPrimesSieve.getLimitAsBigInteger()) < 0
                && SmallPrimesSieve.isPrime(b.intValue());
    }

    public static ArrayList<BigInteger> factor(BigInteger num) {
        ArrayList<BigInteger> factors = new ArrayList<>();

        //fast check for small prime
        if (checkKnownSmallPrime(num)) {
            factors.add(num);
            return factors;
        }

        //start with trial divisions
        num = TrialDivision(num, factors);

        //switch to hard algorithms
        HardFactor(num, factors);

        return factors;
    }

    static BigInteger TrialDivision(BigInteger num, ArrayList<BigInteger> factors) {
        for (int p : PrimesList.SmallPrimes12) {
            BigInteger prime = BigInteger.valueOf(p);
            BigInteger[] qr = num.divideAndRemainder(prime);
            while (qr[1].isZero()) {
                num = qr[0];
                factors.add(prime);
                qr = num.divideAndRemainder(prime);
            }
        }
        return num;
    }

    static void HardFactor(BigInteger num, ArrayList<BigInteger> factors) {
        BigInteger factor;
        while (true) {
            factor = findFactor(num);
            if (factor.isOne() || factor.equals(num)) {
                factors.add(num);
                return;
            } else
                //TODO: ? check whether factor is prime ?
                HardFactor(factor, factors);
            num = num.divide(factor);
        }
    }
}
