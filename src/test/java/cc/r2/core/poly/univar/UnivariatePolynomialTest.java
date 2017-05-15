package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.Integers;
import cc.r2.core.poly.IntegersModulo;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.univar.UnivariatePolynomial.asLongPolyZp;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class UnivariatePolynomialTest extends AbstractPolynomialTest {

    @Test
    public void test1() throws Exception {
        BigInteger modulus = BigInteger.valueOf(59);
        UnivariatePolynomial<BigInteger> aZ = UnivariatePolynomial.create(Integers.Integers, 1, 2, 3, 4, 5, 6);
        IntegersModulo domain = new IntegersModulo(modulus);
        UnivariatePolynomial<BigInteger> aZp = aZ.setDomain(domain);
        lUnivariatePolynomialZp aL = asLongPolyZp(aZp);

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
    public void test3() throws Exception {
        Assert.assertEquals(-1, UnivariatePolynomial.create(0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(-1, UnivariatePolynomial.create(0, 0, 0).firstNonZeroCoefficientPosition());
        Assert.assertEquals(0, UnivariatePolynomial.create(1).firstNonZeroCoefficientPosition());
        Assert.assertEquals(1, UnivariatePolynomial.create(0, 1).firstNonZeroCoefficientPosition());
    }
//
//    @Test
//    public void test3() throws Exception {
//        BigInteger modulus = new BigInteger("998427238390739620139");
//        bMutablePolynomialZp factory = bMutablePolynomialZp.one(modulus);
//
//        BigInteger a = new BigInteger("213471654376351273471236215473").mod(modulus);
//        BigInteger b = new BigInteger("41982734698213476213918476921834").mod(modulus);
//        Assert.assertEquals(b.subtract(a).mod(modulus), factory.subtract(b, a));
//        Assert.assertEquals(a.subtract(b).mod(modulus), factory.subtract(a, b));
//    }
//
//    @Test
//    public void test4() throws Exception {
//        bMutablePolynomialZp factory = bMutablePolynomialZp.zero(BigInteger.valueOf(3));
//        Assert.assertEquals(BigInteger.ZERO, factory.negate(BigInteger.ZERO));
//        Assert.assertEquals(BigInteger.ZERO, factory.negate().lc());
//    }
}