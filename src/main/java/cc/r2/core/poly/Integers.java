package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;

import java.util.Iterator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Integers extends AbstractIntegers {
    public static final Integers Integers = new Integers();

    private Integers() {}

    @Override
    public boolean isField() {return false;}

    @Override
    public BigInteger cardinality() {return null;}

    @Override
    public BigInteger characteristics() {return BigInteger.ZERO;}

    @Override
    public boolean isUnit(BigInteger a) {
        return isOne(a);
    }

    @Override
    public BigInteger add(BigInteger a, BigInteger b) {return a.add(b);}

    @Override
    public BigInteger subtract(BigInteger a, BigInteger b) {return a.subtract(b);}

    @Override
    public BigInteger negate(BigInteger val) {return val.negate();}

    @Override
    public BigInteger multiply(BigInteger a, BigInteger b) {return a.multiply(b);}

    @Override
    public BigInteger[] divideAndRemainder(BigInteger a, BigInteger b) {return a.divideAndRemainder(b);}

    @Override
    public BigInteger remainder(BigInteger a, BigInteger b) {
        return a.mod(b);
    }

    @Override
    public BigInteger reciprocal(BigInteger a) {
        if (isOne(a))
            return a;
        throw new UnsupportedOperationException();
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return val;}

    @Override
    public BigInteger valueOf(long val) {return BigInteger.valueOf(val);}

    @Override
    public BigInteger getNegativeOne() {return BigInteger.NEGATIVE_ONE;}

    @Override
    public boolean isMinusOne(BigInteger bigInteger) {
        return bigInteger.isMinusOne();
    }

    @Override
    public String toString() {return "Z";}

    @Override
    public Iterator<BigInteger> iterator() {
        throw new UnsupportedOperationException("Domain of infinite cardinality.");
    }
}
