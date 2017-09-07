package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * Domain of integers modulo some {@code modulus}.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntegersModulo extends AIntegers {
    private static final long serialVersionUID = 1L;
    /**
     * The modulus.
     */
    public final BigInteger modulus;

    /**
     * Creates Zp domain for specified modulus.
     *
     * @param modulus the modulus
     */
    public IntegersModulo(BigInteger modulus) {
        this.modulus = modulus;
    }

    /**
     * Creates Zp domain for specified modulus.
     *
     * @param modulus the modulus
     */
    public IntegersModulo(long modulus) {
        this(BigInteger.valueOf(modulus));
    }

    @Override
    public boolean isField() {return true;}

    @Override
    public BigInteger cardinality() {return modulus;}

    @Override
    public BigInteger characteristics() {return modulus;}

    @Override
    public boolean isUnit(BigInteger element) {
        return !modulus.divideAndRemainder(element)[1].isZero();
    }

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
    public BigInteger symmetricForm(BigInteger value) {
        return value.compareTo(modulus.shiftRight(1)) <= 0 ? value : value.subtract(modulus);
    }

    /**
     * Converts to a {@link lIntegersModulo }
     */
    public lIntegersModulo asLong() {
        return new lIntegersModulo(modulus.longValueExact());
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
    public BigInteger negate(BigInteger element) {return element.isZero() ? element : modulus.subtract(element);}

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

    /** domain for perfectPowerBase() */
    private IntegersModulo ppBaseDomain = null;

    /**
     * Returns domain for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     *
     * @return domain for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     */
    public IntegersModulo perfectPowerBaseDomain() {
        if (ppBaseDomain == null) {
            synchronized (this) {
                if (ppBaseDomain == null) {
                    BigInteger base = perfectPowerBase();
                    if (base == null)
                        ppBaseDomain = this;
                    else
                        ppBaseDomain = new IntegersModulo(base);
                }
            }
        }

        return ppBaseDomain;
    }

    lIntegersModulo lDomain;

    /**
     * Returns machine integer domain or null if modulus is larger than {@code long}
     *
     * @return machine integer domain or null if modulus is larger than {@code long}
     */
    public lIntegersModulo asMachineSizedDomain() {
        if (!modulus.isLong())
            return null;
        if (lDomain == null)
            synchronized (this) {
                if (lDomain == null)
                    lDomain = new lIntegersModulo(modulus.longValueExact());
            }

        return lDomain;
    }

    @Override
    public String toString() {return "Z/" + modulus;}

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntegersModulo that = (IntegersModulo) o;

        return modulus.equals(that.modulus);
    }

    @Override
    public int hashCode() {
        return modulus.hashCode();
    }
}
