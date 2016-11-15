package cc.r2.core.number;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static cc.r2.core.number.BigInteger.*;

/**
 * Created by poslavsky on 11/11/2016.
 */
public final class FactorInteger {
    private FactorInteger() {
    }

    /**
     * Plain trial division algorithm
     *
     * @param integer integer to factor
     * @return list of factors
     */
    public static List<BigInteger> TrialDivision(BigInteger integer) {
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
        for (; temp.compareTo(i) >= 0
                && integer.compareTo(i.multiply(i)) >= 0; i = i.add(2)) {
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

    private static final BigInteger[] SmallPrimes = {
            BigInteger.valueOf(2), BigInteger.valueOf(3),
            BigInteger.valueOf(5), BigInteger.valueOf(7),
            BigInteger.valueOf(11), BigInteger.valueOf(13),
            BigInteger.valueOf(17), BigInteger.valueOf(19)
    };

    private static final BigInteger PollardB = BigInteger.valueOf(100);

    public static BigInteger PollardP1(BigInteger integer) {
        return PollardP1(integer, TWO, BigInteger.valueOf(1000));
    }

    public static BigInteger PollardP1(BigInteger integer,
                                       BigInteger primeBase,
                                       BigInteger upperBound) {
        BigInteger mb = ONE;
        BigInteger prime = TWO;
        while (prime.compareTo(upperBound) < 0) {
            mb = mb.multiply(prime.pow(log(prime, upperBound)));
            prime = prime.nextProbablePrime();
        }

        BigInteger b = primeBase.modPow(mb, integer);
        BigInteger gcd = b.decrement().gcd(integer);
        if (!gcd.isOne() && !gcd.equals(integer))
            return gcd;

        System.out.println("NOT FOUND");
        upperBound = upperBound.multiply(upperBound);

        BigInteger nextPrime = prime = prime.nextProbablePrime();
        BigInteger cmod, nextCmod = b.modPow(prime, integer);
        do {
            prime = nextPrime;
            cmod = nextCmod;

            gcd = integer.gcd(cmod.decrement());
            if (!gcd.isOne() && !gcd.equals(integer))
                return gcd;

            nextPrime = prime.nextProbablePrime();
//            nextCmod = cmod.multiply(b.modPow(nextPrime.subtract(prime), integer)).mod(integer);
            nextCmod = cmod.modPow(nextPrime.subtract(prime), integer);
        } while (prime.compareTo(upperBound) < 0);

        return null;
    }

    public static BigInteger PollardP1x(BigInteger integer) {
        return PollardP1x(integer, ONE, TWO, ONE, BigInteger.valueOf(100));
    }

    public static BigInteger PollardP1x(BigInteger integer,
                                        final BigInteger mb,
                                        BigInteger primeBase,
                                        BigInteger lowerBound,
                                        BigInteger upperBound) {
        if (lowerBound.compareTo(integer) > 0)
            return null;
        BigInteger prime = lowerBound.nextProbablePrime();
        BigInteger newMb = mb;
        while (prime.compareTo(upperBound) < 0) {
            if(!prime.isPrime())
                continue;
            int log = log(prime, upperBound);
            if(log == 1) log = 2;
            newMb = newMb.multiply(prime.pow(log));
//            System.out.println(prime);
            prime = prime.nextProbablePrime();
        }

        BigInteger gcd = primeBase.modPow(newMb, integer).decrement().gcd(integer);
        if (gcd.isOne()) {
            System.out.println("GCD ONE");
            return PollardP1x(integer, newMb, primeBase, upperBound, upperBound.add(BigInteger.valueOf(1000)));
        } else if (gcd.equals(integer)) {
            System.out.println("GCD SAME");
            return PollardP1x(integer, mb, primeBase.nextProbablePrime(), lowerBound, upperBound);
        } else
            return gcd;
    }

    public static int log(BigInteger base, BigInteger n) {
        double l = (Math.log(n.doubleValue()) / Math.log(base.doubleValue()));
        if (l > Integer.MAX_VALUE / 2)
            throw new IllegalArgumentException();
        return (int) l;
    }

    public static BigInteger PollardRho(BigInteger val, int constant) {

        BigInteger toFactor = val;

        BigInteger divisor;
        BigInteger c  = BigInteger.ONE; //new BigInteger(toFactor.bitLength(), random);
        BigInteger x  = BigInteger.valueOf(2 + constant);//new BigInteger(toFactor.bitLength(), random);
        BigInteger xx = x;

        // check divisibility by 2
        if (toFactor.mod(TWO).compareTo(ZERO) == 0) return TWO;

        do {
            x  =  x.multiply(x).mod(toFactor).add(c).mod(toFactor);
            xx = xx.multiply(xx).mod(toFactor).add(c).mod(toFactor);
            xx = xx.multiply(xx).mod(toFactor).add(c).mod(toFactor);
            divisor = x.subtract(xx).gcd(toFactor);
        } while((divisor.compareTo(ONE)) == 0);

        return divisor;
    }
}
