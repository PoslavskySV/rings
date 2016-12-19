package cc.r2.core.number;

import cc.r2.core.number.primes.BigPrimes;

public final class ModPrimeBigIntegerField extends ModBigIntegerRingAbstract<ModPrimeBigInteger> {
    public ModPrimeBigIntegerField(BigInteger mod) {
        super(mod);
        if (!BigPrimes.isPrime(mod))
            throw new IllegalArgumentException("Not a prime.");
    }

    @Override
    public ModPrimeBigInteger parse(Object o) {
        return new ModPrimeBigInteger(this, BigIntegerRing.IntegerRing.parse(o));
    }

    @Override
    ModPrimeBigInteger once_createOne() {
        return new ModPrimeBigInteger(this, BigInteger.ONE);
    }

    @Override
    ModPrimeBigInteger once_createZero() {
        return new ModPrimeBigInteger(this, BigInteger.ZERO);
    }

    @Override
    public Class<ModPrimeBigInteger> getElementType() {
        return ModPrimeBigInteger.class;
    }
}
