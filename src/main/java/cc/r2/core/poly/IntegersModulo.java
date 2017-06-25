package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IntegersModulo extends AbstractIntegers {
    public final BigInteger modulus;

    public IntegersModulo(BigInteger modulus) {
        this.modulus = modulus;
    }

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
    public boolean isUnit(BigInteger a) {
        return !modulus.divideAndRemainder(a)[1].isZero();
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
    public BigInteger negate(BigInteger val) {return val.isZero() ? val : modulus.subtract(val);}

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
    public BigInteger reciprocal(BigInteger a) {
        return a.modInverse(modulus);
    }

    @Override
    public BigInteger valueOf(BigInteger val) {return modulus(val);}

    @Override
    public BigInteger valueOf(long val) {return valueOf(BigInteger.valueOf(val));}

    @Override
    public BigInteger randomElement(RandomGenerator rnd) {return RandomUtil.randomInt(modulus, rnd);}

    /**
     * if modulus = a^b, a and b are stored in this array
     * if perfectPowerDecomposition[0] ==  null   => the data is not yet initialized
     * if perfectPowerDecomposition[1] ==  null   => modulus is not a perfect power
     */
    private final BigInteger[] perfectPowerDecomposition = new BigInteger[2];

    private void checkPerfectPower() {
        // lazy initialization
        if (perfectPowerDecomposition[0] == null) {
            synchronized ( perfectPowerDecomposition ){
                if (perfectPowerDecomposition[0] != null)
                    return;

                BigInteger[] ipp = BigIntegerArithmetics.perfectPowerDecomposition(modulus);
                if (ipp == null) {
                    // not a perfect power
                    perfectPowerDecomposition[0] = BigInteger.NEGATIVE_ONE;
                    perfectPowerDecomposition[1] = null;
                    return;
                }
                perfectPowerDecomposition[0] = ipp[0];
                perfectPowerDecomposition[1] = ipp[1];
            }
        }
    }

    /**
     * Returns whether the modulus is a perfect power
     *
     * @return whether the modulus is a perfect power
     */
    public boolean isPerfectPower() {
        checkPerfectPower();
        return perfectPowerDecomposition[1] != null;
    }

    /**
     * Returns {@code base} if {@code modulus == base^exponent}, and {@code null} otherwise
     *
     * @return {@code base} if {@code modulus == base^exponent}, and {@code null} otherwise
     */
    public BigInteger perfectPowerBase() {
        if (!isPerfectPower())
            return null;
        return perfectPowerDecomposition[0];
    }

    /**
     * Returns {@code exponent} if {@code modulus == base^exponent}, and {@code null} otherwise
     *
     * @return {@code exponent} if {@code modulus == base^exponent}, and {@code null} otherwise
     */
    public BigInteger perfectPowerExponent() {
        if (!isPerfectPower())
            return null;
        return perfectPowerDecomposition[1];
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
            synchronized ( this ){
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
            synchronized ( this ){
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
