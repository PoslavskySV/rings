package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
abstract class AbstractIntegers implements Domain<BigInteger> {
    @Override
    public final BigInteger getZero() {
        return BigInteger.ZERO;
    }

    @Override
    public final BigInteger getOne() {
        return BigInteger.ONE;
    }

    @Override
    public final boolean isZero(BigInteger bigInteger) {
        return bigInteger.isZero();
    }

    @Override
    public final boolean isOne(BigInteger bigInteger) {
        return bigInteger.isOne();
    }

    @Override
    public final BigInteger parse(String string) {
        return valueOf(new BigInteger(string));
    }

    @Override
    public final BigInteger gcd(BigInteger a, BigInteger b) {
        return a.gcd(b);
    }

    @Override
    public final int signum(BigInteger val) {
        return val.signum();
    }

    @Override
    public final int compare(BigInteger o1, BigInteger o2) {return o1.compareTo(o2);}
}
