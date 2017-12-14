package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.primes.BigPrimes;

import java.util.Iterator;

/**
 * The ring of integers (Z).
 *
 * @since 1.0
 */
public final class Integers extends AIntegers {
    private static final long serialVersionUID = 1L;

    /** The ring of integers (Z) */
    public static final Integers Integers = new Integers();

    private Integers() {}

    @Override
    public boolean isField() {return false;}

    @Override
    public boolean isEuclideanRing() {return true;}

    @Override
    public BigInteger cardinality() {return null;}

    @Override
    public BigInteger characteristic() {return BigInteger.ZERO;}

    @Override
    public boolean isUnit(BigInteger element) {return isOne(element);}

    @Override
    public BigInteger add(BigInteger a, BigInteger b) {return a.add(b);}

    @Override
    public BigInteger subtract(BigInteger a, BigInteger b) {return a.subtract(b);}

    @Override
    public BigInteger negate(BigInteger element) {return element.negate();}

    @Override
    public BigInteger multiply(BigInteger a, BigInteger b) {return a.multiply(b);}

    @Override
    public BigInteger[] divideAndRemainder(BigInteger a, BigInteger b) {return a.divideAndRemainder(b);}

    @Override
    public BigInteger remainder(BigInteger a, BigInteger b) {return a.mod(b);}

    @Override
    public BigInteger reciprocal(BigInteger element) {
        if (isOne(element) || isMinusOne(element))
            return element;
        throw new UnsupportedOperationException();
    }

    @Override
    public BigInteger pow(BigInteger base, int exponent) {
        return base.pow(exponent);
    }

    @Override
    public BigInteger pow(BigInteger base, long exponent) {
        if (exponent < Integer.MAX_VALUE)
            return pow(base, (int) exponent);
        return super.pow(base, exponent);
    }

    @Override
    public BigInteger pow(BigInteger base, BigInteger exponent) {
        if (exponent.isLong())
            return pow(base, exponent.longValueExact());
        return super.pow(base, exponent);
    }

    @Override
    public final BigInteger gcd(BigInteger a, BigInteger b) {
        return a.gcd(b);
    }

    @Override
    public FactorDecomposition<BigInteger> factorSquareFree(BigInteger element) {
        return factor(element);
    }

    @Override
    public FactorDecomposition<BigInteger> factor(BigInteger element) {
        return FactorDecomposition.of(this, BigPrimes.primeFactors(element));
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return val;}

    @Override
    public BigInteger valueOf(long val) {return BigInteger.valueOf(val);}

    @Override
    public BigInteger getNegativeOne() {return BigInteger.NEGATIVE_ONE;}

    @Override
    public boolean isMinusOne(BigInteger bigInteger) {return bigInteger.isMinusOne();}

    @Override
    public final int signum(BigInteger element) {
        return element.signum();
    }

    @Override
    public BigInteger abs(BigInteger el) {
        return el.abs();
    }

    @Override
    public String toString() {return "Z";}

    @Override
    public Iterator<BigInteger> iterator() {
        throw new UnsupportedOperationException("Ring of infinite cardinality.");
    }
}
