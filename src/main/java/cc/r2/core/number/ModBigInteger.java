package cc.r2.core.number;

public final class ModBigInteger extends ModBigIntegerAbstract
        implements RingElement<ModBigInteger> {
    final ModBigIntegerRing ring;

    ModBigInteger(boolean asis, ModBigIntegerRing ring, BigInteger val) {
        super(asis, ring.mod, val);
        this.ring = ring;
    }

    public ModBigInteger(ModBigIntegerRing ring, BigInteger val) {
        super(ring.mod, val);
        this.ring = ring;
    }

    @Override
    public ModBigInteger add(ModBigInteger a) {
        return new ModBigInteger(true, ring, add0(a.val));
    }

    @Override
    public ModBigInteger subtract(ModBigInteger a) {
        return new ModBigInteger(true, ring, subtract0(a.val));
    }

    @Override
    public ModBigInteger multiply(ModBigInteger a) {
        return new ModBigInteger(true, ring, multiply0(a.val));
    }

    @Override
    public ModBigInteger negate() {
        if (val.isZero()) return this;
        return new ModBigInteger(ring, val.negate());
    }

    @Override
    public ModBigInteger getZero() {
        return ring.getZero();
    }

    @Override
    public ModBigInteger getOne() {
        return ring.getOne();
    }

    @Override
    public Ring<ModBigInteger> getRing() {
        return ring;
    }
}
