package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * Ring of integers modulo some {@code modulus}.
 *
 * @since 1.0
 */
public final class IntegersZp extends AIntegers {
    private static final long serialVersionUID = 1L;
    /**
     * The modulus.
     */
    public final BigInteger modulus;

    /**
     * Creates Zp ring for specified modulus.
     *
     * @param modulus the modulus
     */
    public IntegersZp(BigInteger modulus) {
        this.modulus = modulus;
    }

    /**
     * Creates Zp ring for specified modulus.
     *
     * @param modulus the modulus
     */
    public IntegersZp(long modulus) {
        this(BigInteger.valueOf(modulus));
    }

    @Override
    public boolean isField() {return true;}

    @Override
    public boolean isEuclideanRing() {return true;}

    @Override
    public BigInteger cardinality() {return modulus;}

    @Override
    public BigInteger characteristic() {return modulus;}

    @Override
    public boolean isUnit(BigInteger element) {
        return !element.isZero() && !modulus.divideAndRemainder(element)[1].isZero();
    }

    /**
     * Returns {@code val mod this.modulus}
     *
     * @param val the integer
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
    public BigInteger symmetricForm(BigInteger value) {
        return value.compareTo(modulus.shiftRight(1)) <= 0 ? value : value.subtract(modulus);
    }

    /**
     * Converts to a {@link IntegersZp64 }
     */
    public IntegersZp64 asMachineRing() {
        return new IntegersZp64(modulus.longValueExact());
    }

    @Override
    public BigInteger add(BigInteger a, BigInteger b) {
        a = valueOf(a); b = valueOf(b);
        BigInteger r = a.add(b), rm = r.subtract(modulus);
        return rm.signum() >= 0 ? rm : r;
    }

    @Override
    public BigInteger subtract(BigInteger a, BigInteger b) {
        a = valueOf(a); b = valueOf(b);
        BigInteger r = a.subtract(b);
        return r.signum() < 0 ? r.add(modulus) : r;
    }

    @Override
    public BigInteger negate(BigInteger element) {return element.isZero() ? element : modulus.subtract(valueOf(element));}

    @Override
    public BigInteger multiply(BigInteger a, BigInteger b) {return modulus(a.multiply(b));}

    @Override
    public BigInteger[] divideAndRemainder(BigInteger a, BigInteger b) {
        return new BigInteger[]{divide(a, b), BigInteger.ZERO};
    }

    public BigInteger divide(BigInteger a, BigInteger b) {
        return multiply(a, b.modInverse(modulus));
    }

    @Override
    public BigInteger remainder(BigInteger a, BigInteger b) {
        return getZero();
    }

    @Override
    public BigInteger reciprocal(BigInteger element) {
        return element.modInverse(modulus);
    }

    @Override
    public FactorDecomposition<BigInteger> factorSquareFree(BigInteger element) {
        return factor(element);
    }

    @Override
    public FactorDecomposition<BigInteger> factor(BigInteger element) {
        return FactorDecomposition.of(this, element);
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return modulus(val);}

    @Override
    public BigInteger valueOf(long val) {return valueOf(BigInteger.valueOf(val));}

    @Override
    public BigInteger randomElement(RandomGenerator rnd) {return RandomUtil.randomInt(modulus, rnd);}

    @Override
    public Iterator<BigInteger> iterator() {
        return new It();
    }

    private final class It implements Iterator<BigInteger> {
        private BigInteger val = BigInteger.ZERO;

        @Override
        public boolean hasNext() {
            return val.compareTo(modulus) < 0;
        }

        @Override
        public BigInteger next() {
            BigInteger r = val;
            val = val.increment();
            return r;
        }
    }

    /** ring for perfectPowerBase() */
    private IntegersZp ppBaseDomain = null;

    /**
     * Returns ring for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     *
     * @return ring for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     */
    public IntegersZp perfectPowerBaseDomain() {
        if (ppBaseDomain == null) {
            synchronized (this) {
                if (ppBaseDomain == null) {
                    BigInteger base = perfectPowerBase();
                    if (base == null)
                        ppBaseDomain = this;
                    else
                        ppBaseDomain = new IntegersZp(base);
                }
            }
        }

        return ppBaseDomain;
    }

    private IntegersZp64 lDomain;

    /**
     * Returns machine integer ring or null if modulus is larger than {@code long}
     *
     * @return machine integer ring or null if modulus is larger than {@code long}
     */
    public IntegersZp64 asZp64() {
        if (!modulus.isLong())
            return null;
        if (lDomain == null)
            synchronized (this) {
                if (lDomain == null)
                    lDomain = new IntegersZp64(modulus.longValueExact());
            }

        return lDomain;
    }

    @Override
    public String toString() {return "Z/" + modulus;}

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntegersZp that = (IntegersZp) o;

        return modulus.equals(that.modulus);
    }

    @Override
    public int hashCode() {
        return modulus.hashCode();
    }
}
