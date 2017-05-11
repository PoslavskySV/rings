package cc.r2.core.poly;

import cc.redberry.libdivide4j.FastDivision.*;

import static cc.redberry.libdivide4j.FastDivision.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lModularDomain {
    /** the modulus */
    public final long modulus;
    /** magic **/
    public final Magic magic, magic32MulMod;
    /** whether modulus less then 2^32 (if so, faster mulmod available) **/
    public final boolean modulusFits32;

    public lModularDomain(long modulus, Magic magic, Magic magic32MulMod, boolean modulusFits32) {
        this.modulus = modulus;
        this.magic = magic;
        this.magic32MulMod = magic32MulMod;
        this.modulusFits32 = modulusFits32;
    }

    public lModularDomain(long modulus) {
        this(modulus, magicSigned(modulus), magic32ForMultiplyMod(modulus), LongArithmetics.fits31bitWord(modulus));
    }

    /** modulus operation */
    public long mod(long val) {
        return modSignedFast(val, magic);
    }

    /** reduce data mod modulus **/
    public void mod(long[] data) {
        for (int i = 0; i < data.length; ++i)
            data[i] = mod(data[i]);
    }

    /** multiplyMod operation */
    public long multiplyMod(long a, long b) {
        return modulusFits32 ? mod(a * b) : multiplyMod128Unsigned(a, b, modulus, magic32MulMod);
    }

    /** addMod operation */
    public long addMod(long a, long b) {
        long r = a + b;
        return r - modulus >= 0 ? r - modulus : r;
    }

    /** subtractMod operation */
    public long subtractMod(long a, long b) {
        long r = a - b;
        return r + ((r >> 63)&modulus);
    }

    public long reciprocal(long val) {
        return LongArithmetics.modInverse(val, modulus);
    }

    /** negateMod operation */
    public long negateMod(long val) {
        return val == 0 ? val : modulus - val;
    }

    /** to symmetric modulus */
    public long symmetricForm(long value) {
        return value <= modulus / 2 ? value : value - modulus;
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
                result = multiplyMod(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiplyMod(k2p, k2p);
        }
    }
}
