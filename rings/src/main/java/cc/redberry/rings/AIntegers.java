package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;

/**
 * @since 1.0
 */
abstract class AIntegers extends ARing<BigInteger> {
    private static final long serialVersionUID = 1L;

    @Override
    public final BigInteger getZero() {
        return BigInteger.ZERO;
    }

    @Override
    public final BigInteger getOne() {
        return BigInteger.ONE;
    }

    @Override
    public final boolean isZero(BigInteger element) {
        return element.isZero();
    }

    @Override
    public final boolean isOne(BigInteger element) {
        return element.isOne();
    }

    @Override
    public final BigInteger parse(String string) {
        return valueOf(new BigInteger(string.trim()));
    }

    @Override
    public final int compare(BigInteger o1, BigInteger o2) {return o1.compareTo(o2);}

    @Override
    public final BigInteger[] createArray(int length) {
        return new BigInteger[length];
    }

    @Override
    public final BigInteger[][] createArray2d(int length) {
        return new BigInteger[length][];
    }

    @Override
    public final BigInteger[][] createArray2d(int m, int n) {
        return new BigInteger[m][n];
    }

    @Override
    public final BigInteger valueOfBigInteger(BigInteger val) {
        return valueOf(val);
    }

    @Override
    public BigInteger copy(BigInteger element) {
        return element;
    }
}
