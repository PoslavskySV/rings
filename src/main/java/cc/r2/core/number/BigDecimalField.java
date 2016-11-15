package cc.r2.core.number;

public final class BigDecimalField implements Ring<BigDecimal> {
    public static final BigDecimalField BigDecimalField = new BigDecimalField();

    private BigDecimalField() {
    }

    @Override
    public BigDecimal getOne() {
        return BigDecimal.ONE;
    }

    @Override
    public BigDecimal getZero() {
        return BigDecimal.ZERO;
    }

    @Override
    public BigDecimal parse(Object o) {
        if (o instanceof BigDecimal)
            return (BigDecimal) o;
        if (o instanceof String)
            return new BigDecimal((String) o);
        throw new IllegalArgumentException(o.toString());
    }

    @Override
    public Class<BigDecimal> getElementType() {
        return BigDecimal.class;
    }
}
