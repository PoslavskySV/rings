package cc.r2.core.number;


public final class BigIntegerRing implements Ring<BigInteger> {
    public static final BigIntegerRing IntegerRing = new BigIntegerRing();

    @Override
    public BigInteger getOne() {
        return BigInteger.ONE;
    }

    @Override
    public BigInteger getZero() {
        return BigInteger.ZERO;
    }

    @Override
    public BigInteger parse(Object o) {
        if (o instanceof BigInteger)
            return (BigInteger) o;
        return new BigInteger((String) o);
    }

    @Override
    public Class<BigInteger> getElementType() {
        return BigInteger.class;
    }
}
