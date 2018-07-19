package cc.redberry.rings.poly;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.AlgebraicNumberField;
import static cc.redberry.rings.Rings.MultivariateRing;
import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.Rings.UnivariateRing;

/**
 *
 */
public class AlgebraicNumberFieldTest extends APolynomialTest {
    @Test
    public void test1() {
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field = AlgebraicNumberField(UnivariatePolynomial.create(Q, Q.valueOf(-2), Q.valueOf(0), Q.valueOf(1)));
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        Assert.assertEquals(coder.parse("2"), coder.parse("s^2"));
        Assert.assertEquals(coder.parse("s/2"), coder.parse("1/s"));
        Assert.assertEquals(coder.parse("s - 1"), coder.parse("1/(1 + s)"));
        Assert.assertEquals(coder.parse("-s"), coder.parse("-s"));
    }

    @Test
    public void test2() {
        UnivariatePolynomial<BigInteger> mp = UnivariatePolynomial.create(-22, 12, 13, -112343242, 3, 1, 0, 0, 0, 123, 1234134, 324, 123423442, 0, 1).square().increment();
        AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicNumberField(mp);
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkUnivariateCoder(field, "s");
        UnivariatePolynomial<BigInteger> element = coder.parse("- 2 + 13*s + s^2 - 13*s^4");
        UnivariatePolynomial<BigInteger>[] ann = field.normalizer2(element);
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

            AlgebraicNumberField<UnivariatePolynomial<BigInteger>> field = AlgebraicNumberField(minPoly);
            for (int j = 0; j < 10; ++j) {
                UnivariatePolynomial<BigInteger> element = field.randomElement(rnd);
                element = element.setRing(Zp(100)).setRingUnsafe(Z);

                if (element.isZero()) {
                    --j;
                    continue;
                }

                UnivariatePolynomial<BigInteger>[] ann2 = field.normalizer2(element);

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

    @Test
    public void test6() {
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> pCoder = Coder.mkUnivariateCoder(UnivariateRing(Q), "x");
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = new AlgebraicNumberField<>(pCoder.parse("-12 -x + 13*x^2 - 11*x^3 + x^7"));

        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");
        Coder<UnivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(UnivariateRing(field), cfCoder, "x");

        Assert.assertEquals(
                pCoder.parse("144 - 313*x^3 + 436*x^6 - 3237*x^9 + 4581*x^12 - 5994*x^15 + 17561*x^18 - 7875*x^21 + 11443*x^24 - 11268*x^27 + 10009*x^30 - 4236*x^33 + 703*x^36 + 22*x^39 - 12*x^42"),
                field.normOfPolynomial(coder.parse("s^2 - (1-s^3)*x^3 + (1 + s)*x^6")));
    }

    @Test
    public void test7() {
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> pCoder = Coder.mkUnivariateCoder(UnivariateRing(Q), "x");
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = new AlgebraicNumberField<>(pCoder.parse("-12 -x + 13*x^2 - 11*x^3 + x^7"));

        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(MultivariateRing(3, field), cfCoder, "x", "y", "z");

        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> mCoder = Coder.mkMultivariateCoder(MultivariateRing(3, Q), "x", "y", "z");

        Assert.assertEquals(
                mCoder.parse("144 - 313*y*z^2 + 301*x^2*y^2*z^2 + 135*y^2*z^4 - 1699*x^2*y^3*z^4 + 400*x^4*y^4*z^4 - 1538*y^3*z^6 + 3898*x^2*y^4*z^6 - 2563*x^4*y^5*z^6 + 143*x^6*y^6*z^6 + 283*y^4*z^8 - 4422*x^2*y^5*z^8 + 6958*x^4*y^6*z^8 - 1525*x^6*y^7*z^8 - 82*x^8*y^8*z^8 + 991*y^5*z^10 + 5089*x^2*y^6*z^10 - 13592*x^4*y^7*z^10 + 6547*x^6*y^8*z^10 - 427*x^8*y^9*z^10 + 72*x^10*y^10*z^10 + 5371*y^6*z^12 - 2728*x^2*y^7*z^12 + 11346*x^4*y^8*z^12 - 11581*x^6*y^9*z^12 + 2727*x^8*y^10*z^12 - 474*x^10*y^11*z^12 + 31*x^12*y^12*z^12 + 9970*y^7*z^14 - 6368*x^2*y^8*z^14 + 740*x^4*y^9*z^14 + 7210*x^6*y^10*z^14 - 3762*x^8*y^11*z^14 + 672*x^10*y^12*z^14 + 22*x^12*y^13*z^14 - 12*x^14*y^14*z^14"),
                field.normOfPolynomial(coder.parse("s^2 - (1-s^3)*y*z^2 + (1 + s)*x^2*y^2*z^2")));
    }

    @Test
    public void test8() {
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> pCoder = Coder.mkUnivariateCoder(UnivariateRing(Q), "x");
        AlgebraicNumberField<UnivariatePolynomial<Rational<BigInteger>>> field
                = new AlgebraicNumberField<>(pCoder.parse("-12 -x + 13*x^2 - 11*x^3 + x^7"));

        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkUnivariateCoder(field, "s");

        Assert.assertEquals(
                pCoder.parse("4153 - 6474*x - 1314*x^2 + 2957*x^3 - 818*x^4 + 62*x^5 - 7*x^6 + x^7"),
                field.minimalPolynomial(cfCoder.parse("s - s^3 + 1")));
        Assert.assertEquals(cfCoder.parse("7"), field.trace(cfCoder.parse("s - s^3 + 1")));
        Assert.assertEquals(cfCoder.parse("21154/123"), field.trace(cfCoder.parse("12*s^2 - s^13/123 + 11*s/12 + 1")));
        Assert.assertEquals(cfCoder.parse("21154/123"), field.trace(cfCoder.parse("-12*s^2 - s^13/123 + 11*s/12 + 1")));
    }
}