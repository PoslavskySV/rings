package cc.redberry.rings.poly;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.Rings.UnivariateRing;

/**
 *
 */
public class AlgebraicNumberFieldTest extends APolynomialTest {
    @Test
    public void test1() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicExtension(UnivariatePolynomial.create(Q, Q.valueOf(-2), Q.valueOf(0), Q.valueOf(1)));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        Assert.assertEquals(coder.parse("2"), coder.parse("s^2"));
        Assert.assertEquals(coder.parse("s/2"), coder.parse("1/s"));
        Assert.assertEquals(coder.parse("s - 1"), coder.parse("1/(1 + s)"));
        Assert.assertEquals(coder.parse("-s"), coder.parse("-s"));
    }

    @Test
    public void test2() {
        UnivariatePolynomial<BigInteger> mp = UnivariatePolynomial.create(-2, 0, 0, 0, 0, 0, 1);
        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicExtension(mp);
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        UnivariatePolynomial<BigInteger> element = coder.parse("- 2 + 13*s + s^2 - 13*s^4");
        UnivariatePolynomial<BigInteger>[] ann = field.cancellingMultiplier2(element);
        UnivariatePolynomial<BigInteger> expected = ann[1];
        UnivariatePolynomial<BigInteger> actual = field.multiply(element, ann[0]);
        Assert.assertEquals(expected, actual);
        Assert.assertTrue(expected.isConstant());
    }

    @Test
    public void test3_random() {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < 30; ++i) {
            UnivariatePolynomial<BigInteger> minPoly = IrreduciblePolynomials.randomIrreduciblePolynomialOverZ(rndd.nextInt(2, 10), rnd);
            minPoly.setLC(Z.getOne());
            if (!IrreduciblePolynomials.irreducibleQ(minPoly)) {
                --i;
                continue;
            }

            AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicExtension(minPoly);
            for (int j = 0; j < 10; ++j) {
                UnivariatePolynomial<BigInteger> element = field.randomElement(rnd);
                element = element.setRing(Zp(100)).setRingUnsafe(Z);

                if (element.isZero()) {
                    --j;
                    continue;
                }

                UnivariatePolynomial<BigInteger>[] ann2 = field.cancellingMultiplier2(element);

                UnivariatePolynomial<BigInteger> expected = ann2[1];
                UnivariatePolynomial<BigInteger> actual = field.multiply(element, ann2[0]);
                Assert.assertEquals(expected, actual);
                Assert.assertTrue(expected.isConstant());
            }
        }
    }

    @Test
    public void test4() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = new AlgebraicNumberField<>(UnivariatePolynomial.create(-2, 0, 1).mapCoefficients(Q, Q::mkNumerator));

        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        Assert.assertEquals(coder.parse("-1"), field.norm(coder.parse("1 + s")));
        Assert.assertEquals(coder.parse("-14"), field.norm(coder.parse("-2 + 3*s")));
    }

    @Test
    public void test5() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = new AlgebraicNumberField<>(Coder.mkUnivariateCoder(UnivariateRing(Q), "x").parse("-12 -x + 13*x^2 - 11*x^3 + x^7"));

        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        Assert.assertEquals(coder.parse("2077536306"), field.norm(coder.parse("-2 - 2*s + 3*s^5")));
    }
}