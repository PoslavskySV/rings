package cc.r2.core.poly;

import cc.redberry.libdivide4j.FastDivision.*;
import org.apache.commons.math3.random.RandomGenerator;

import static cc.redberry.libdivide4j.FastDivision.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lIntegersModulo {
    /** the modulus */
    public final long modulus;
    /** magic **/
    public final Magic magic, magic32MulMod;
    /** whether modulus less then 2^32 (if so, faster mulmod available) **/
    public final boolean modulusFits32;

    public lIntegersModulo(long modulus, Magic magic, Magic magic32MulMod, boolean modulusFits32) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32MulMod;
        this.modulusFits32 = modulusFits32;
    }

    public lIntegersModulo(long modulus) {
        this(modulus, magicSigned(modulus), magic32ForMultiplyMod(modulus), LongArithmetics.fits31bitWord(modulus));
    }

    /** Returns {@code val % this.modulus} */
    public long modulus(long val) {
        return modSignedFast(val, magic);
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
        return r + ((r >> 63)&modulus);
    }

    /** Subtract mod operation */
    public long divide(long a, long b) {
        return multiply(a, reciprocal(b));
    }

    /** Returns modular inverse of {@code val} */
    public long reciprocal(long val) {
        return LongArithmetics.modInverse(val, modulus);
    }

    /** Negate mod operation */
    public long negate(long val) {
        return val == 0 ? val : modulus - val;
    }

    /** to symmetric modulus */
    public long symmetricForm(long value) {
        return value <= modulus / 2 ? value : value - modulus;
    }

    public IntegersModulo asDomain() {
        return new IntegersModulo(modulus);
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
            if ((exponent&1) != 0)
                result = multiply(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiply(k2p, k2p);
        }
    }

    /**
     * Returns a random element from this domain
     *
     * @param rnd the source of randomness
     * @return random element from this domain
     */
    public long randomElement(RandomGenerator rnd) {
        return modulus(rnd.nextLong());
    }

    /**
     * Returns a random non zero element from this domain
     *
     * @param rnd the source of randomness
     * @return random non zero element from this domain
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        lIntegersModulo that = (lIntegersModulo) o;

        return modulus == that.modulus;
    }

    @Override
    public String toString() {return "Z/" + modulus;}

    @Override
    public int hashCode() {
        return (int) (modulus^(modulus >>> 32));
    }
}
