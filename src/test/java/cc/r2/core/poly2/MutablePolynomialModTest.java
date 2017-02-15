package cc.r2.core.poly2;

import cc.r2.core.number.primes.BigPrimes;
import cc.r2.core.polynomial.MutablePolynomial;
import cc.r2.core.polynomial.RandomPolynomials;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by poslavsky on 27/01/2017.
 */
public class MutablePolynomialModTest {
    private static long subtractMod(long a, long b, long modulus) {
        long r = a - b;
        return r + ((r >> 63)&modulus);
    }

    @Test
    public void name() throws Exception {
        long a = 1345, b = 142445, m = 12314324L;

        System.out.println(subtractMod(a, b, m));
        System.out.println(Math.floorMod(a - b, m));

    }

    //
//    @Test
//    public void asdasd() throws Exception {
//        System.out.println(Long.toBinaryString(1L << 63));
//        System.out.println(Long.numberOfTrailingZeros(1L << 63));
//
//    }
//
//    @Test
//    public void test1() throws Exception {
//        long[] data = {66, 513, 304, -894, 1000, 163, 796, 219, 319, -919, 396};
//        long modulus = SmallPrimes.nextPrime(1000);
//        MutablePolynomial oldPoly = MutablePolynomial.create(data.clone()).modulus(modulus);
//        MutablePolynomialMod newPoly = MutablePolynomialMod.create(modulus, data.clone());
//
//        System.out.println(oldBench(oldPoly, modulus));
//        System.out.println(newBench(newPoly));
//
//        System.out.println(oldPoly);
//        System.out.println(newPoly);
//
//    }
//
    static long oldBench(MutablePolynomial m, long modulus) {
        MutablePolynomial old = m.clone();
        int deg = m.degree;
        for (int i = 0; i < 5; i++) {
            MutablePolynomial oldM = m.clone();
            m = m.multiply(m.add(m.shiftLeft(m.degree / 3), modulus), modulus).subtract(oldM, modulus);
            m = m.cut(deg);
            m = m.subtract(MutablePolynomial.derivative(m, modulus), modulus);
            m = m.multiply(old, modulus);
        }
        return m.evaluate(modulus / 3, modulus);
    }

    //
    static MutablePolynomialMod newBench0(MutablePolynomialMod m) {
        MutablePolynomialMod old = m.clone();
        int deg = m.degree;
        for (int i = 0; i < 5; i++) {
            MutablePolynomialMod oldM = m.clone();
            m = m.multiply(m.add(m.shiftLeft(m.degree / 3))).subtract(oldM);
            m = m.cut(deg);
            m = m.subtract(m.derivative());
            m = m.multiply(old);
        }
        return m;
    }

    static MutablePolynomialMod newBench1(MutablePolynomialMod tmp) {
        MutablePolynomialMod old = tmp.clone();
        for (int i = 0; i < 5; i++) {
            MutablePolynomialMod oldM = tmp.clone();
            tmp = tmp.add(old);
            tmp = tmp.multiply(old);
            tmp = tmp.subtract(oldM);
            tmp = tmp.subtract(oldM.derivative());
            tmp = tmp.multiply(old);
        }
        return tmp;
    }

    static long newBench(MutablePolynomialMod m) {
        return newBench0(m).evaluate(m.modulus / 3);
    }

    @Test
    public void test1132() throws Exception {
        long modulus = BigPrimes.nextPrime(1L << 59);
        RandomGenerator rnd = new Well1024a();
        int degree = 130;
        long[] data = RandomPolynomials.randomLongArray(degree, modulus, rnd);
        MutablePolynomialMod polyMod = MutablePolynomialMod.createSigned(modulus, data.clone());
        System.out.println(newBench1(polyMod));
    }

    @Test
    public void test1sad132() throws Exception {
        long modulus = BigPrimes.nextPrime(1L << 59);
        RandomGenerator rnd = new Well1024a();
        int degree = 130;
        long[] data = RandomPolynomials.randomLongArray(degree, modulus, rnd);
        MutablePolynomialMod p = MutablePolynomialMod.createSigned(modulus, data);
        System.out.println("mod=" + modulus + ";");
        System.out.println("poly= " + p + ";");
        System.out.println("red = " + newBench1(p) + ";");
    }

    @Test
    public void testPerformance() throws Exception {
        DescriptiveStatistics oldf = new DescriptiveStatistics(), newf = new DescriptiveStatistics();

        int degree = 30;
        long modulus = BigPrimes.nextPrime(1L << 31);

        RandomGenerator rnd = new Well1024a();
        for (int i = 0; i < 15000; i++) {
            if (i == 10000) {
                System.out.println("clear");
                oldf.clear();
                newf.clear();
            }

            long[] data = RandomPolynomials.randomLongArray(degree, modulus, rnd);

            MutablePolynomial poly = MutablePolynomial.create(data.clone());
            poly.modulus(modulus);
            long start;
            start = System.nanoTime();
            long oldR = oldBench(poly, modulus);
            oldf.addValue(System.nanoTime() - start);

            MutablePolynomialMod polyMod = MutablePolynomialMod.createSigned(modulus, data.clone());
            start = System.nanoTime();
            long newR = newBench(polyMod);
            newf.addValue(System.nanoTime() - start);

            Assert.assertEquals(oldR, newR);
        }

        System.out.println("==== old ==== ");
        System.out.println(oldf);

        System.out.println("==== new ==== ");
        System.out.println(newf);
    }
}