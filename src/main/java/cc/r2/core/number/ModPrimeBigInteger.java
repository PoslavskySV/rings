package cc.r2.core.number;

public final class ModPrimeBigInteger extends ModBigIntegerAbstract
        implements FieldElement<ModPrimeBigInteger> {
    private ModPrimeBigIntegerField ring;

    ModPrimeBigInteger(boolean asis, ModPrimeBigIntegerField ring, BigInteger val) {
        super(asis, ring.mod, val);
        this.ring = ring;
    }

    public ModPrimeBigInteger(ModPrimeBigIntegerField ring, BigInteger val) {
        super(ring.mod, val);
        this.ring = ring;
    }

    @Override
    public ModPrimeBigInteger add(ModPrimeBigInteger a) {
        return new ModPrimeBigInteger(true, ring, add(a.val));
    }

    @Override
    public ModPrimeBigInteger subtract(ModPrimeBigInteger a) {
        return new ModPrimeBigInteger(true, ring, subtract(a.val));
    }

    @Override
    public ModPrimeBigInteger multiply(ModPrimeBigInteger a) {
        return new ModPrimeBigInteger(true, ring, multiply(a.val));
    }

    @Override
    public ModPrimeBigInteger divide(ModPrimeBigInteger a) {
        if (a.isZero())
            throw new IllegalArgumentException("Division by zero");
        //solve: this = q * a (mod p)
        return new ModPrimeBigInteger(ring, val.multiply(a.val.modPow(BigInteger.NEGATIVE_ONE, mod)));
    }

    @Override
    public ModPrimeBigInteger getZero() {
        return ring.getZero();
    }

    @Override
    public ModPrimeBigInteger getOne() {
        return ring.getOne();
    }

    @Override
    public ModPrimeBigInteger negate() {
        if (val.isZero()) return this;
        return new ModPrimeBigInteger(ring, val.negate());
    }

    @Override
    public boolean isZero() {
        return val.isZero();
    }

    @Override
    public boolean isOne() {
        return val.isOne();
    }

    @Override
    public Ring<ModPrimeBigInteger> getRing() {
        return ring;
    }

    @Override
    public ModPrimeBigInteger gcd(ModPrimeBigInteger oth) {
        return getOne();
    }

    @Override
    public ModPrimeBigInteger[] divideAndRemainder(ModPrimeBigInteger oth) {
        return new ModPrimeBigInteger[]{divide(oth), getZero()};
    }

    @Override
    public int compareTo(ModPrimeBigInteger o) {
        return val.compareTo(o.val);
    }
}
