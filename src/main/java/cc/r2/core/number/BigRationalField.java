package cc.r2.core.number;


public final class BigRationalField implements Ring<BigRational> {
    public static final BigRationalField BigRationalField = new BigRationalField();

    @Override
    public BigRational getOne() {
        return BigRational.ONE;
    }

    @Override
    public BigRational getZero() {
        return BigRational.ZERO;
    }

    @Override
    public BigRational parse(Object o) {
        if (o instanceof BigInteger)
            return new BigRational((BigInteger) o, BigInteger.ONE);
        else throw new RuntimeException();
    }

    @Override
    public Class<BigRational> getElementType() {
        return BigRational.class;
    }
}
