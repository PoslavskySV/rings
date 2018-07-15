package cc.redberry.rings.poly;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.multivar.Monomial;
import cc.redberry.rings.poly.multivar.MultivariateFactorization;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariateFactorization;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;
import org.junit.Test;

import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.Rings.UnivariateRing;
import static org.junit.Assert.assertEquals;

/**
 *
 */
public class MultipleFieldExtensionTest {

    @Test
    public void test1() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomial<Rational<BigInteger>>
                minPolySqrt2 = auxCoder.parse("x^2 - 2"),
                minPolySqrt3 = auxCoder.parse("x^2 + 3");

        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                > field = MultipleFieldExtension.mkMultipleExtension(minPolySqrt2, minPolySqrt3);

        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(field, "s2", "s3");
        assertEquals(coder.parse("-6"), field.pow(coder.parse("s2 * s3"), 2));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void test2() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomial<Rational<BigInteger>>
                minPolySqrt2 = auxCoder.parse("x^2 - 2"),
                minPolySqrt3 = auxCoder.parse("x^2 + 3"),
                minPolySqrt5 = auxCoder.parse("x^2 + 5");

        for (int i = 0; i < 2; ++i) {
            long start = System.nanoTime();
            MultipleFieldExtension<
                    Monomial<Rational<BigInteger>>,
                    MultivariatePolynomial<Rational<BigInteger>>,
                    UnivariatePolynomial<Rational<BigInteger>>
                    > field = MultipleFieldExtension.mkMultipleExtension(minPolySqrt2, minPolySqrt3, minPolySqrt5);
            System.out.println("Create extension: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(field, "s2", "s3", "s5");

            assertEquals(coder.parse("(2 * (-3) * (-5)) + 1 / (2 * (-3) * (-5))"), coder.parse("(s2 * s3 * s5)^2 + 1 / (s2 * s3 * s5)^2"));
            System.out.println("Arithmetic      : " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println();
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void test3() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomial<Rational<BigInteger>>
                minPolySqrt2 = auxCoder.parse("x^2 - 2"),
                minPolySqrt3 = auxCoder.parse("x^2 - 3"),
                minPolySqrt5 = auxCoder.parse("x^2 - 5");
        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                > field = MultipleFieldExtension.mkMultipleExtension(minPolySqrt2, minPolySqrt3, minPolySqrt5);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkPolynomialCoder(field, "s2", "s3", "s5");

        UnivariateRing<UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> uRing = UnivariateRing(field);
        Coder<UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkUnivariateCoder(uRing, cfCoder, "x");

        UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> poly = coder.parse("((x - s2) * (x - s3) * (x - s5))^2");
        PolynomialFactorDecomposition<UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> factors = UnivariateFactorization.Factor(poly);

        assertEquals(3, factors.size());
        assertEquals(6, factors.sumExponents());
        assertEquals(poly, factors.multiply());
    }

    @Test
    @SuppressWarnings("unchecked")
    public void test4() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomial<Rational<BigInteger>>
                minPolySqrt2 = auxCoder.parse("x^2 - 2"),
                minPolySqrt3 = auxCoder.parse("x^2 - 3"),
                minPolySqrt5 = auxCoder.parse("x^2 - 5");
        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                > field = MultipleFieldExtension.mkMultipleExtension(minPolySqrt2, minPolySqrt3, minPolySqrt5);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> cfCoder = Coder.mkPolynomialCoder(field, "s2", "s3", "s5");

        MultivariateRing<MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> mRing = Rings.MultivariateRing(3, field);
        Coder<MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?> coder = Coder.mkMultivariateCoder(mRing, cfCoder, "x", "y", "z");

        MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> poly = coder.parse("((x*y - s2) * (y*z - s3) * (x*z - s5))^2");
        PolynomialFactorDecomposition<MultivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> factors = MultivariateFactorization.Factor(poly);

        assertEquals(3, factors.size());
        assertEquals(6, factors.sumExponents());
        assertEquals(poly, factors.multiply());
    }
}