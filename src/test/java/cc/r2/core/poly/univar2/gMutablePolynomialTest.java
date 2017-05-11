package cc.r2.core.poly.univar2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.IntegersDomain;
import cc.r2.core.poly.ModularDomain;
import cc.r2.core.poly.univar.bMutablePolynomialZ;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.univar2.gMutablePolynomial.asLongPolyZp;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class gMutablePolynomialTest extends AbstractPolynomialTest {


    @Test
    public void test1() throws Exception {
        BigInteger modulus = BigInteger.valueOf(59);
        gMutablePolynomial<BigInteger> aZ = gMutablePolynomial.create(IntegersDomain.IntegersDomain, 1, 2, 3, 4, 5, 6);
        ModularDomain domain = new ModularDomain(modulus);
        gMutablePolynomial<BigInteger> aZp = aZ.setDomain(domain);
        lMutablePolynomialZp aL = asLongPolyZp(aZp);

        for (int i = 0; i < 5; i++) {
//            a = (a.clone() * a.clone().decrement() - a.clone().derivative() + (a.clone().square())) * a.clone();
            aZp = (aZp.clone().multiply(aZp.clone().decrement()).subtract(aZp.clone().derivative()).add(aZp.clone().square())).multiply(aZp.clone());
            aZp = aZp.truncate(aZp.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            aZp = aZp.subtract(aZp.derivative()).decrement();
//            a = a.multiply(a.derivative().increment().truncate(10));

            aZ = (aZ.clone().multiply(aZ.clone().decrement()).subtract(aZ.clone().derivative()).add(aZ.clone().square())).multiply(aZ.clone());
            aZ = aZ.truncate(aZ.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            aZ = aZ.subtract(aZ.derivative()).decrement();
//            aZ = aZ.multiply(aZ.derivative().increment().truncate(10));

            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
            aL = aL.subtract(aL.derivative()).decrement();
//            aL = aL.multiply(aL.derivative().increment().truncate(10));
        }

        System.out.println(aZp.degree);
        Assert.assertEquals(aL, asLongPolyZp(aZp));
        Assert.assertEquals(aZp, aZ.setDomain(domain));
    }

    @Test
    public void test2() throws Exception {
        Assert.assertEquals(
                bMutablePolynomialZ.parse("1 + x+ x^23"),
                bMutablePolynomialZ.parse("x+1+x^23"));
    }
//
//    @Test
//    public void test3() throws Exception {
//        BigInteger modulus = new BigInteger("998427238390739620139");
//        bMutablePolynomialZp factory = bMutablePolynomialZp.one(modulus);
//
//        BigInteger a = new BigInteger("213471654376351273471236215473").mod(modulus);
//        BigInteger b = new BigInteger("41982734698213476213918476921834").mod(modulus);
//        Assert.assertEquals(b.subtract(a).mod(modulus), factory.subtractMod(b, a));
//        Assert.assertEquals(a.subtract(b).mod(modulus), factory.subtractMod(a, b));
//    }
//
//    @Test
//    public void test4() throws Exception {
//        bMutablePolynomialZp factory = bMutablePolynomialZp.zero(BigInteger.valueOf(3));
//        Assert.assertEquals(BigInteger.ZERO, factory.negateMod(BigInteger.ZERO));
//        Assert.assertEquals(BigInteger.ZERO, factory.negate().lc());
//    }


    static int LIM = 3;

    private static BigInteger test(gMutablePolynomial<BigInteger> aL) {
        for (int i = 0; i < LIM; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
        }
        return aL.evaluate(aL.degree);
    }

    private static BigInteger test(bMutablePolynomialZp aL) {
        for (int i = 0; i < LIM; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree() * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
        }
        return aL.evaluate(BigInteger.valueOf(aL.degree()));
    }

    private static BigInteger test(bMutablePolynomialZ aL) {
        for (int i = 0; i < LIM; i++) {
            aL = (aL.clone().multiply(aL.clone().decrement()).subtract(aL.clone().derivative()).add(aL.clone().square())).multiply(aL.clone());
            aL = aL.truncate(aL.degree() * 3 / 2).shiftRight(2).shiftLeft(2).increment().negate();
        }
        return aL.evaluate(BigInteger.valueOf(aL.degree()));
    }


    @Test
    public void test4() throws Exception {
        DescriptiveStatistics old = new DescriptiveStatistics(), nev = new DescriptiveStatistics();
        int nIterations = 800;
        LIM = 2;
        BigInteger modulus = BigInteger.valueOf(getModulusRandom(50));
        modulus = modulus.nextProbablePrime();

        System.out.println(modulus);

        for (int i = 0; i < nIterations; i++) {
            if (i % 10 == 0)
                System.out.println(i);
            if (i == nIterations / 10) {
                old.clear();
                nev.clear();
            }
            BigInteger[] arr = RandomUtil.randomBigIntegerArray(100, BigInteger.ZERO, modulus, getRandom());

            long start = System.nanoTime();
            BigInteger X = test(bMutablePolynomialZ.create(arr).modulus(modulus));
            old.addValue(System.nanoTime() - start);

            ModularDomain domain = new ModularDomain(modulus);
            start = System.nanoTime();
            BigInteger Y = test(gMutablePolynomial.create(domain, arr));
            nev.addValue(System.nanoTime() - start);

            assert X.equals(Y);
        }

        System.out.println(old.getPercentile(0.5));
        System.out.println(nev.getPercentile(0.5));
    }

}