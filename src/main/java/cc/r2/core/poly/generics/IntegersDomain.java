package cc.r2.core.poly.generics;

import cc.r2.core.number.BigInteger;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntegersDomain extends AbstractIntegersDomain {
    public static final IntegersDomain IntegersDomain = new IntegersDomain();

    private IntegersDomain() {}

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
    public BigInteger reciprocal(BigInteger a) {
        throw new UnsupportedOperationException();
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return val;}

    @Override
    public BigInteger valueOf(long val) {return BigInteger.valueOf(val);}
}
