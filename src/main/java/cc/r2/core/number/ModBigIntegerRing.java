package cc.r2.core.number;

public final class ModBigIntegerRing extends ModBigIntegerRingAbstract<ModBigInteger> {
    public ModBigIntegerRing(BigInteger mod) {
        super(mod);
    }

    @Override
    ModBigInteger once_createOne() {
        return new ModBigInteger(this, BigInteger.ONE);
    }

    @Override
    ModBigInteger once_createZero() {
        return new ModBigInteger(this, BigInteger.ZERO);
    }

    @Override
    public ModBigInteger parse(Object o) {
        return new ModBigInteger(this, BigIntegerRing.IntegerRing.parse(o));
    }

    @Override
    public Class<ModBigInteger> getElementType() {
        return ModBigInteger.class;
    }
}
