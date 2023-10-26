package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecompositionTest;
import cc.redberry.rings.poly.FiniteField;
import cc.redberry.rings.poly.PolynomialFactorDecomposition;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.util.TimeUnits;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static cc.redberry.rings.poly.multivar.MultivariateGCDTest.createMonomial;
import static cc.redberry.rings.poly.multivar.MultivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics;
import static cc.redberry.rings.poly.multivar.MultivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics;

/**
 * @since 1.0
 */
public class MultivariateSquareFreeFactorizationTest extends AMultivariateTest {

    @Test
    public void test1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("x^2 + y^2*x + z"),
                b = MultivariatePolynomial.parse("x^2 - 2*y^2*x - z"),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13"),
                poly = a.square().multiply(b.square().square()).multiply(c.square().square());

        for (int i = 0; i < AbstractTest.its(1, 5); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> yun = MultivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics(poly);
            System.out.println("Yun: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(poly, yun);

            start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> mus = SquareFreeFactorizationMusserZeroCharacteristics(poly);
            System.out.println("Musser: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            FactorDecompositionTest.assertFactorization(poly, mus);
            Assert.assertEquals(yun.size(), mus.size());
        }
    }

    @Test
    public void test2() throws Exception {
        IntegersZp domain = new IntegersZp(7);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("11 + x^7*y^14 + z^7*x^14 + y^7", domain),
                b = MultivariatePolynomial.parse("11 + 3*y^7*z^14 + 4*x^7*y^14 + 5*z^7", domain),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13", domain),
                poly = a.square().multiply(b.square()).multiply(c.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asOverZp64(poly);
        for (int i = 0; i < its(1, 5); i++) {
            long start = System.nanoTime();
            FactorDecompositionTest.assertFactorization(lPoly, MultivariateSquareFreeFactorization.SquareFreeFactorizationBernardin(lPoly));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }

//        905ms
//        877ms
//        911ms
//        694ms
//        684ms
//        729ms
//        854ms
//        891ms
//        805ms
//        613ms
    }

    @Test
    public void test3() throws Exception {
        IntegersZp domain = new IntegersZp(2);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("11 + x^7*y^14 + z^7*x^14 + y^7", domain),
                b = MultivariatePolynomial.parse("11 + 3*y^7*z^14 + 4*x^7*y^14 + 5*z^7", domain),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13", domain),
                poly = a.square().multiply(b.square()).multiply(c.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asOverZp64(poly);
        FactorDecompositionTest.assertFactorization(lPoly, MultivariateSquareFreeFactorization.SquareFreeFactorizationBernardin(lPoly));
    }

    @Test
    public void test4() throws Exception {
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("11 + x^7*y^14 + z^7*x^14 + y^7", domain),
                b = MultivariatePolynomial.parse("11 + 3*y^7*z^14 + 4*x^7*y^14 + 5*z^7", domain),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13", domain),
                poly = a.square().multiply(b.square()).multiply(c.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asOverZp64(poly);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationBernardin(lPoly);
        FactorDecompositionTest.assertFactorization(lPoly, decomposition);
    }

    @Test
    public void test5() throws Exception {
        IntegersZp domain = new IntegersZp(2);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("1 + a^6*b^14 + a^2*b^4 + a^7", domain),
                b = MultivariatePolynomial.parse("1 + a^3*b^4 + a + b", domain),
                poly = a.square().multiply(b.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asOverZp64(poly);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationBernardin(lPoly);
        FactorDecompositionTest.assertFactorization(lPoly, decomposition);
    }

    @Test
    public void test5_finiteField() throws Exception {
        long modulus = 2;

        FiniteField<UnivariatePolynomialZp64> field = new FiniteField<>(IrreduciblePolynomials.randomIrreduciblePolynomial(modulus, 4, AbstractTest.getRandom()));
        MultivariatePolynomial<UnivariatePolynomialZp64>
                a = MultivariatePolynomial.zero(3, field, MonomialOrder.LEX)
                .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 2, 3, 4, 5).modulus(modulus)), 1, 1, 3))
                .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 3, 2, 13).modulus(modulus)), 3, 2, 1))
                .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 13, 12, 13).modulus(modulus)), 0, 2, 1)),
                b = MultivariatePolynomial.zero(3, field, MonomialOrder.LEX)
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5).modulus(modulus)), 1, 1, 13))
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 1, 2, 13).modulus(modulus)), 2, 2, 1))
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 113, 112, 13).modulus(modulus)), 10, 2, 1)),
                c = MultivariatePolynomial.one(3, field, MonomialOrder.LEX)
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5, 12).modulus(modulus)), 11, 1, 13))
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(11, 2, 1, 1, 2, 13).modulus(modulus)), 21, 2, 1))
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(modulus)), 10, 12, 1))
                        .add(createMonomial(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(modulus)), 0, 0, 1)),
                poly = a.square().multiply(b.square()).multiply(c.square());


        PolynomialFactorDecomposition<MultivariatePolynomial<UnivariatePolynomialZp64>> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationBernardin(poly);
        FactorDecompositionTest.assertFactorization(poly, decomposition);
    }

    @Test
    public void test6() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e", "f", "g", "h"};
        for (Ring<BigInteger> ring : Arrays.<Ring<BigInteger>>asList(new IntegersZp(2), Rings.Z)) {
            MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.parse("a^2*b^4*c*e^5", ring, vars);
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> expected = PolynomialFactorDecomposition.empty(poly);
            expected.addFactor(MultivariatePolynomial.parse("a", ring, vars), 2);
            expected.addFactor(MultivariatePolynomial.parse("b", ring, vars), 4);
            expected.addFactor(MultivariatePolynomial.parse("c", ring, vars), 1);
            expected.addFactor(MultivariatePolynomial.parse("e", ring, vars), 5);

            Assert.assertEquals(expected, MultivariateSquareFreeFactorization.SquareFreeFactorization(poly));
        }
    }

    @Test
    public void test7() {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("2*y^3-3*x*y^2+x^3");

        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>>
                r = SquareFreeFactorizationYunZeroCharacteristics(a);

        FactorDecompositionTest.assertFactorization(a, r);
    }
}
