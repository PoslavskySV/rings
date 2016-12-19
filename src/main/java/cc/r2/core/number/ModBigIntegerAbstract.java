package cc.r2.core.number;

abstract class ModBigIntegerAbstract {
    protected final BigInteger mod;
    protected final BigInteger val;

    ModBigIntegerAbstract(BigInteger mod, BigInteger val) {
        this.mod = mod;
        this.val = val.mod(mod);
    }

    ModBigIntegerAbstract(boolean asis, BigInteger mod, BigInteger val) {
        this.mod = mod;
        this.val = val;
    }

    public BigInteger value() {
        return val;
    }

    BigInteger mod0(BigInteger a) {
        return a.mod(mod);
    }

    BigInteger add0(BigInteger a) {
        return mod0(val.add(a));
    }

    BigInteger subtract0(BigInteger a) {
        return mod0(val.subtract(a));
    }

    BigInteger multiply0(BigInteger a) {
        return mod0(val.multiply(a));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o)
            return true;
        if (!(o instanceof ModBigIntegerAbstract))
            return false;
        ModBigIntegerAbstract oth = (ModBigIntegerAbstract) o;
        return mod.equals(oth.mod) && val.equals(oth.val);

    }

    @Override
    public int hashCode() {
        return mod.hashCode() + 11 * val.hashCode();
    }

    @Override
    public String toString() {
        return val.toString();
    }
}
