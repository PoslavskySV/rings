package cc.redberry.rings;

import cc.redberry.libdivide4j.FastDivision.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MachineArithmetic;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well44497b;

import static cc.redberry.libdivide4j.FastDivision.*;

/**
 * Zp ring over machine numbers which provides fast modular arithmetic.
 *
 * @see cc.redberry.libdivide4j.FastDivision
 * @since 1.0
 */
public final class IntegersZp64 implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    /** the modulus */
    public final long modulus;
    /** magic **/
    public final Magic magic, magic32MulMod;
    /** whether modulus less then 2^32 (if so, faster mulmod available) **/
    public final boolean modulusFits32;

    public IntegersZp64(long modulus, Magic magic, Magic magic32MulMod, boolean modulusFits32) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32MulMod;
        this.modulusFits32 = modulusFits32;
    }

    /**
     * Creates the ring.
     *
     * @param modulus the modulus
     */
    public IntegersZp64(long modulus) {
        this(modulus, magicSigned(modulus), magic32ForMultiplyMod(modulus), MachineArithmetic.fits31bitWord(modulus));
    }

    /** Returns {@code val % this.modulus} */
    public long modulus(long val) {
        return modSignedFast(val, magic);
    }

    /** Returns {@code val % this.modulus} */
    public long modulus(BigInteger val) {
        return val.isLong() ? modSignedFast(val.longValue(), magic) : val.mod(BigInteger.valueOf(modulus)).longValue();
    }

    /** Inplace sets elements of {@code data} to {@code data % this.modulus} */
    public void modulus(long[] data) {
        for (int i = 0; i < data.length; ++i)
            data[i] = modulus(data[i]);
    }

    /** Multiply mod operation */
    public long multiply(long a, long b) {
        return modulusFits32 ? modulus(a * b) : multiplyMod128Unsigned(a, b, modulus, magic32MulMod);
    }

    /** Add mod operation */
    public long add(long a, long b) {
        long r = a + b;
        return r - modulus >= 0 ? r - modulus : r;
    }

    /** Subtract mod operation */
    public long subtract(long a, long b) {
        long r = a - b;
        return r + ((r >> 63) & modulus);
    }

    /** Subtract mod operation */
    public long divide(long a, long b) {
        return multiply(a, reciprocal(b));
    }


    /** cached modular inverses */
    private volatile int[] cachedReciprocals = null;

    /** builds a table of cached reciprocals */
    public void buildCachedReciprocals() {
        if (cachedReciprocals != null)
            return;
        synchronized (this) {
            if (cachedReciprocals == null) {
                int[] cachedReciprocals = new int[MachineArithmetic.safeToInt(modulus)];
                for (int val = 1; val < cachedReciprocals.length; ++val)
                    cachedReciprocals[val] = (int) MachineArithmetic.modInverse(val, modulus);
                this.cachedReciprocals = cachedReciprocals;
            }
        }
    }

    /** Returns modular inverse of {@code val} */
    public long reciprocal(long val) {
        return cachedReciprocals == null
                ? MachineArithmetic.modInverse(val, modulus)
                : (val < cachedReciprocals.length ? cachedReciprocals[(int) val] : MachineArithmetic.modInverse(val, modulus));
    }

    /** Negate mod operation */
    public long negate(long val) {
        return val == 0 ? val : modulus - val;
    }

    /** to symmetric modulus */
    public long symmetricForm(long value) {
        return value <= modulus / 2 ? value : value - modulus;
    }

    /**
     * Converts this to a generic ring over big integers
     *
     * @return generic ring
     */
    public IntegersZp asGenericRing() {
        return new IntegersZp(modulus);
    }

    /**
     * Returns {@code base} in a power of non-negative {@code e} modulo {@code magic.modulus}
     *
     * @param base     the base
     * @param exponent the non-negative exponent
     * @return {@code base} in a power of {@code e}
     */
    public long powMod(final long base, long exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();
        if (exponent == 0)
            return 1;

        long result = 1;
        long k2p = base;
        for (; ; ) {
            if ((exponent & 1) != 0)
                result = multiply(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiply(k2p, k2p);
        }
    }

    /**
     * Returns a random element from this ring
     *
     * @param rnd the source of randomness
     * @return random element from this ring
     */
    public long randomElement(RandomGenerator rnd) {
        return modulus(rnd.nextLong());
    }

    /**
     * Returns a random element from this ring
     *
     * @return random element from this ring
     */
    public long randomElement() {
        return randomElement(new Well44497b(System.nanoTime()));
    }

    /**
     * Returns a random non zero element from this ring
     *
     * @param rnd the source of randomness
     * @return random non zero element from this ring
     */
    public long randomNonZeroElement(RandomGenerator rnd) {
        long el;
        do {
            el = randomElement(rnd);
        } while (el == 0);
        return el;
    }

    /**
     * Gives value!
     *
     * @param value the number
     * @return value!
     */
    public long factorial(int value) {
        long result = 1;
        for (int i = 2; i <= value; ++i)
            result = multiply(result, i);
        return result;
    }

    /**
     * if modulus = a^b, a and b are stored in this array
     */
    private final long[] perfectPowerDecomposition = {-1, -1};

    private void checkPerfectPower() {
        // lazy initialization
        if (perfectPowerDecomposition[0] == -1) {
            synchronized (perfectPowerDecomposition) {
                if (perfectPowerDecomposition[0] != -1)
                    return;

                long[] ipp = MachineArithmetic.perfectPowerDecomposition(modulus);
                if (ipp == null) {
                    // not a perfect power
                    perfectPowerDecomposition[0] = modulus;
                    perfectPowerDecomposition[1] = 1;
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
        return perfectPowerExponent() > 1;
    }

    /**
     * Returns {@code base} if {@code modulus == base^exponent}, and {@code -1} otherwisec
     *
     * @return {@code base} if {@code modulus == base^exponent}, and {@code -1} otherwisec
     */
    public long perfectPowerBase() {
        checkPerfectPower();
        return perfectPowerDecomposition[0];
    }

    /**
     * Returns {@code exponent} if {@code modulus == base^exponent}, and {@code -1} otherwisec
     *
     * @return {@code exponent} if {@code modulus == base^exponent}, and {@code -1} otherwisec
     */
    public long perfectPowerExponent() {
        checkPerfectPower();
        return perfectPowerDecomposition[1];
    }

    /** ring for perfectPowerBase() */
    private IntegersZp64 ppBaseDomain = null;

    /**
     * Returns ring for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     *
     * @return ring for {@link #perfectPowerBase()} or {@code this} if modulus is not a perfect power
     */
    public IntegersZp64 perfectPowerBaseDomain() {
        if (ppBaseDomain == null) {
            synchronized (this) {
                if (ppBaseDomain == null) {
                    long base = perfectPowerBase();
                    if (base == -1)
                        ppBaseDomain = this;
                    else
                        ppBaseDomain = new IntegersZp64(base);
                }
            }
        }

        return ppBaseDomain;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        IntegersZp64 that = (IntegersZp64) o;

        return modulus == that.modulus;
    }

    @Override
    public String toString() {return "Z/" + modulus;}

    @Override
    public int hashCode() {
        return (int) (modulus ^ (modulus >>> 32));
    }
}
