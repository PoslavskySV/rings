package cc.redberry.rings.poly.test;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.test.TimeConsuming;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.lang.annotation.ElementType;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

/**
 * @since 1.0
 */
public class APolynomialTest extends AbstractTest {

    /** random prime number with specified number of bits */
    public static long getModulusRandom(int nBits) {
        if (nBits <= 1 || nBits > MachineArithmetic.MAX_SUPPORTED_MODULUS_BITS)
            throw new IllegalArgumentException();
        return BigPrimes.nextPrime(getRandomData().nextLong(1L << (nBits - 1), (1L << nBits) - 1));
    }

    public static long[] getSmallModulusArray(int n) {
        return getModulusArray(n, 0, 0);
    }

    public static long[] getLargeModulusArray(int n, int maxModulusBits) {
        return getModulusArray(0, n, maxModulusBits);
    }

    public static long[] getModulusArray(int nSmall, int nLarge, int smallModulusBits, int maxModulusBits) {
        long[] res = new long[nSmall + nLarge];
        int i = 0;
        for (; i < nSmall; i++)
            res[i] = getModulusRandom(getRandomData().nextInt(2, smallModulusBits));
        for (; i < res.length; i++)
            res[i] = getModulusRandom(getRandomData().nextInt(32, maxModulusBits));
        return res;
    }

    public static long[] getModulusArray(int nSmall, int nLarge, int maxModulusBits) {
        return getModulusArray(nSmall, nLarge, 31, maxModulusBits);
    }

    public static BigInteger[] getProbablePrimesArray(BigInteger from, int n) {
        BigInteger[] res = new BigInteger[n];
        res[0] = from.nextProbablePrime();
        for (int i = 1; i < n; i++)
            res[i] = res[i - 1].nextProbablePrime();
        return res;
    }

    public static long[] getOneSmallOneLargeModulus(int maxModulusBits) {
        return getModulusArray(1, 1, maxModulusBits);
    }

    @Test
    @TimeConsuming
    public void test0() throws Exception {
        for (int nBits = 2; nBits < 60; nBits++) {
            for (int i = 0; i < 10; i++) {
                int modulusBits = 64 - Long.numberOfLeadingZeros(getModulusRandom(nBits));
                Assert.assertTrue(nBits <= modulusBits && modulusBits <= nBits + 1);
            }
        }
    }

    protected static String getSingularPath() {
        return System.getProperty("singularPath", "/Applications/Singular.app/Contents/bin/Singular");
    }

    protected static boolean isSingularAvailable() {
        return new File(getSingularPath()).exists();
    }

    @Before
    @Override
    public void beforeMethod() throws Exception {
        if (getClass().getMethod(name.getMethodName()).isAnnotationPresent(RequiresSingular.class))
            Assume.assumeTrue(isSingularAvailable());
        super.beforeMethod();
    }

    @Target(ElementType.METHOD)
    @Retention(RetentionPolicy.RUNTIME)
    protected @interface RequiresSingular {}
}
