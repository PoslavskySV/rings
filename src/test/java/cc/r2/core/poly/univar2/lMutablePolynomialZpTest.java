package cc.r2.core.poly.univar2;

import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lMutablePolynomialZpTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        lMutablePolynomialZp aL = lMutablePolynomialZ.create(1, 2, 3, 4, 5, 6).modulus(59);
        for (int i = 0; i < 5; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            Assert.assertTrue(check(aL));
        }
    }

    @Test
    public void test2() throws Exception {
        lMutablePolynomialZp factory = lMutablePolynomialZp.zero(3);
        Assert.assertEquals(0, factory.domain.negateMod(0));
        Assert.assertEquals(0, factory.negate().lc());
    }

    static int LIM = 3;

    private static long test(lMutablePolynomialZp aL) {
        for (int i = 0; i < LIM; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
        }
        return aL.evaluate(aL.degree);
    }

    private static long test0(cc.r2.core.poly.univar.lMutablePolynomialZp aL) {
        for (int i = 0; i < LIM; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree() * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
        }
        return aL.evaluate(aL.degree());
    }

    @Test
    public void test3() throws Exception {
        DescriptiveStatistics old = new DescriptiveStatistics(), nev = new DescriptiveStatistics();
        int nIterations = 3000;
        LIM = 2;
        for (int i = 0; i < nIterations; i++) {
            if (i == nIterations / 10) {
                old.clear();
                nev.clear();
            }
            long modulus = getModulusRandom(15);
            long[] arr = RandomUtil.randomLongArray(100, 0, modulus, getRandom());

            long start = System.nanoTime();
            long X = test(lMutablePolynomialZ.create(arr).modulus(modulus));
            nev.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            long Y = test0(cc.r2.core.poly.univar.lMutablePolynomialZ.create(arr).modulus(modulus));
            old.addValue(System.nanoTime() - start);

            assert X == Y;
        }

        System.out.println(old.getPercentile(0.5));
        System.out.println(nev.getPercentile(0.5));
    }

    private static boolean check(lMutablePolynomialZp poly) {
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] >= poly.domain.modulus)
                return false;
        }
        return true;
    }
}