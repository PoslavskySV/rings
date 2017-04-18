package cc.r2.core.poly.generics;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class ModularDomain extends AbstractIntegersDomain {
    public final BigInteger modulus;

    public ModularDomain(BigInteger modulus) {
        this.modulus = modulus;
    }

    public ModularDomain(long modulus) {
        this(BigInteger.valueOf(modulus));
    }

    @Override
    public boolean isField() {return true;}

    @Override
    public BigInteger size() {return modulus;}

    /**
     * Returns {@code val mod this.modulus}
     *
     * @param val some integer
     * @return {@code val mod this.modulus}
     */
    public BigInteger modulus(BigInteger val) {
        return (val.signum() >= 0 && val.compareTo(modulus) < 0) ? val : val.mod(modulus);
    }

    /**
     * Converts {@code value} to a symmetric representation of Zp
     *
     * @param value field element
     * @return {@code value} in a symmetric representation of Zp
     */
    public BigInteger symMod(BigInteger value) {
        return value.compareTo(modulus.shiftRight(1)) <= 0 ? value : value.subtract(modulus);
    }

    @Override
    public BigInteger add(BigInteger a, BigInteger b) {
        BigInteger r = a.add(b), rm = r.subtract(modulus);
        return rm.signum() >= 0 ? rm : r;
    }

    @Override
    public BigInteger subtract(BigInteger a, BigInteger b) {
        BigInteger r = a.subtract(b);
        return r.signum() < 0 ? r.add(modulus) : r;
    }

    @Override
    public BigInteger negate(BigInteger val) {return val.isZero() ? val : modulus.subtract(val);}

    @Override
    public BigInteger multiply(BigInteger a, BigInteger b) {return modulus(a.multiply(b));}

    @Override
    public BigInteger[] divideAndRemainder(BigInteger a, BigInteger b) {
        return new BigInteger[]{multiply(a, b.modInverse(modulus)), BigInteger.ZERO};
    }

    @Override
    public BigInteger reciprocal(BigInteger a) {
        return a.modInverse(modulus);
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return modulus(val);}

    @Override
    public BigInteger valueOf(long val) {return valueOf(BigInteger.valueOf(val));}

    @Override
    public BigInteger randomElement(RandomGenerator rnd) {return RandomUtil.randomInt(modulus, rnd);}
}
