package cc.redberry.rings.poly;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.multivar.*;
import cc.redberry.rings.poly.univar.UnivariateFactorization;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.TimeUnits;
import org.junit.Test;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.Rings.UnivariateRing;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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

    @Test
    public void test5() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomial<Rational<BigInteger>> poly = auxCoder.parse("2*x^3 - 5");
        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                > splittingField = MultipleFieldExtension.mkSplittingField(poly);
        assertEquals(6, splittingField.getSimpleExtension().degree());

        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(splittingField, "s1", "s2", "s3");
        assertEquals(coder.parse("5/2"), coder.parse("s1 * s2 * s3"));
        assertEquals(coder.parse("0"), coder.parse("s1 * s2  +  s1 * s3 + s2 * s3"));
        assertEquals(coder.parse("0"), coder.parse("s1 + s2 + s3"));
    }

    @Test
    public void test6() {
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> auxRing = UnivariateRing(Q);
        Coder<UnivariatePolynomial<Rational<BigInteger>>, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");


        // some irreducible polynomial
        UnivariatePolynomial<Rational<BigInteger>> poly = auxCoder.parse("2*x^3 - 3*x^2 + 4*x + 5");
        // create splitting field as multiple field extension
        // s1,s2,s3 are roots of specified poly
        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                >
                splittingField = MultipleFieldExtension.mkSplittingField(poly);
        // check the degree of this extension (6 = 3!)
        assertEquals(6, splittingField.getSimpleExtension().degree());

        // assert Vieta's identities
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> coder = Coder.mkPolynomialCoder(splittingField, "s1", "s2", "s3");
        assertEquals(coder.parse("-5/2"), coder.parse("s1 * s2 * s3"));
        assertEquals(coder.parse("2"), coder.parse("s1 * s2  +  s1 * s3 + s2 * s3"));
        assertEquals(coder.parse("3/2"), coder.parse("s1 + s2 + s3"));

        UnivariateRing<UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>> pRing = UnivariateRing(splittingField);
        Coder<UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?> pCoder = Coder.mkUnivariateCoder(pRing, coder, "x");

        UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> p = pCoder.parse("2*x^3 - 3*x^2 + 4*x + 5");
        assertEquals(3, UnivariateFactorization.Factor(p).size());

        MultivariatePolynomial<Rational<BigInteger>> el = splittingField.variable(0);
        System.out.println(splittingField.subtract(splittingField.image(splittingField.inverse(el)), el));
        assertTrue(splittingField.subtract(splittingField.image(splittingField.inverse(el)), el).isZero());
    }

    @Test
    public void test6a() {
        UnivariateRing<UnivariatePolynomialZp64> auxRing = UnivariateRingZp64(SmallPrimes.nextPrime(1 << 12));
        Coder<UnivariatePolynomialZp64, ?, ?> auxCoder = Coder.mkPolynomialCoder(auxRing, "x");

        UnivariatePolynomialZp64 poly = auxCoder.parse("17*x^3 - 14*x^2 + 25*x +  15");
        MultipleFieldExtension<
                MonomialZp64,
                MultivariatePolynomialZp64,
                UnivariatePolynomialZp64
                > splittingField = MultipleFieldExtension.mkSplittingField(poly);
        assertEquals(3, splittingField.getSimpleExtension().degree());

        Coder<MultivariatePolynomialZp64, ?, ?> coder = Coder.mkPolynomialCoder(splittingField, "s1", "s2", "s3");
        assertEquals(coder.parse("-15/17"), coder.parse("s1 * s2 * s3"));
        assertEquals(coder.parse("25/17"), coder.parse("s1 * s2  +  s1 * s3 + s2 * s3"));
        assertEquals(coder.parse("14/17"), coder.parse("s1 + s2 + s3"));
    }

    @Test
    public void test7() {
        UnivariatePolynomial<Rational<BigInteger>> minPoly1 = UnivariatePolynomial.create(3, 0, 0, 1).mapCoefficients(Q, Q::mkNumerator);
        MultipleFieldExtension<
                Monomial<Rational<BigInteger>>,
                MultivariatePolynomial<Rational<BigInteger>>,
                UnivariatePolynomial<Rational<BigInteger>>
                > field = MultipleFieldExtension.mkMultipleExtension(minPoly1);

        MultivariatePolynomial<Rational<BigInteger>> alpha1 = field.variable(0);
        UnivariatePolynomial<MultivariatePolynomial<Rational<BigInteger>>> minPoly2 = UnivariatePolynomial.create(field, alpha1, field.valueOf(3), field.pow(alpha1, 2));

        field = field.joinAlgebraicElement(minPoly2);
        MultivariatePolynomial<Rational<BigInteger>> el = field.variable(0);
        assertTrue(field.inverse(field.subtract(field.image(field.inverse(el)), el)).isZero());
    }
}