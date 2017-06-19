package cc.r2.core.poly.univar;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.lIntegersModulo;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lUnivariatePolynomialZpTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        lUnivariatePolynomialZp aL = lUnivariatePolynomialZ.create(1, 2, 3, 4, 5, 6).modulus(59);
        for (int i = 0; i < 5; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            Assert.assertTrue(check(aL));
        }
    }

    @Test
    public void test2() throws Exception {
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(3);
        Assert.assertEquals(0, factory.domain.negate(0));
        Assert.assertEquals(0, factory.negate().lc());
    }

    @Test
    public void test4() throws Exception {
        System.out.println(lUnivariatePolynomialZ.create(0).firstNonZeroCoefficientPosition());
    }

    private static boolean check(lUnivariatePolynomialZp poly) {
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] >= poly.domain.modulus)
                return false;
        }
        return true;
    }

    @Ignore
    @Test
    public void test_performance() throws Exception {
        int nIterations = 10000;
        DescriptiveStatistics
                mul_a = new DescriptiveStatistics(),
                mul_p = new DescriptiveStatistics(),
                mul_k = new DescriptiveStatistics();
        RandomGenerator rnd = getRandom();
        int arrSize = 100;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 100) {
                mul_a.clear(); mul_p.clear(); mul_k.clear();
            }
            long[] a = new long[arrSize + rnd.nextInt(arrSize)];
            long[] b = new long[arrSize + rnd.nextInt(arrSize)];
            long modulus = (1 << 10) + rnd.nextInt(1 << 10);
            for (int j = 0; j < a.length; j++) {
                a[j] = rnd.nextInt((int) modulus);
            }

            for (int j = 0; j < b.length; j++) {
                b[j] = rnd.nextInt((int) modulus);
            }

            lIntegersModulo domain = new lIntegersModulo(modulus);

            long start;
            start = System.nanoTime();
            long[] a_r = mul_a(a, b, modulus, domain);
            mul_a.addValue(System.nanoTime() - start);

            lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(domain);
            start = System.nanoTime();
            long[] c_r = new long[a.length + b.length - 1];
            factory.multiplyClassicalSafeNoTrick(c_r, a, 0, a.length, b, 0, b.length);
            mul_p.addValue(System.nanoTime() - start);
            assert Arrays.equals(c_r, a_r);

//            start = System.nanoTime();
//            long[] d_r = factory.multiplyKaratsubaSafe(a, 0, a.length, b, 0, b.length);
//            mul_k.addValue(System.nanoTime() - start);
//            assert Arrays.equals(d_r, a_r);
        }

        System.out.println(mul_a.getMean() / 1e6);
        System.out.println(mul_p.getMean() / 1e6);
        System.out.println(mul_k.getMean() / 1e6);
    }

    private static long[] mul_a(long[] a, long[] b, long p, lIntegersModulo domain) {
        long p2 = p * p;
        int
                d_a = a.length - 1,
                d_b = b.length - 1,
                d_c = d_a + d_b;

        long result[] = new long[d_c + 1];
        for (int k = 0; k <= d_c; ++k) {
            long t = 0;
            for (int i = Math.max(0, k - d_b), to = Math.min(k, d_a); i <= to; ++i) {
                if (t < 0) ;
                else t = t - p2;
                t = t + a[i] * b[k - i];
            }
            t = domain.modulus(t);
            if (t < 0)
                t = t + p;
            result[k] = t;
        }
        return result;
    }
}