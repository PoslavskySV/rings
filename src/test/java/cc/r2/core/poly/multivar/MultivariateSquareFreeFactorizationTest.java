package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.test.APolynomialTest;
import cc.r2.core.poly.univar.IrreduciblePolynomials;
import cc.r2.core.poly.univar.UnivariatePolynomialZ64;
import cc.r2.core.poly.univar.UnivariatePolynomialZp64;
import cc.r2.core.util.TimeUnits;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static cc.r2.core.poly.FactorDecompositionTest.assertFactorization;
import static cc.r2.core.poly.multivar.MonomialOrder.LEX;
import static cc.r2.core.poly.multivar.MultivariateSquareFreeFactorization.SquareFreeFactorization;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateSquareFreeFactorizationTest extends APolynomialTest {

    @Test
    public void test1() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("x^2 + y^2*x + z"),
                b = MultivariatePolynomial.parse("x^2 - 2*y^2*x - z"),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13"),
                poly = a.square().multiply(b.square().square()).multiply(c.square().square());

        for (int i = 0; i < its(1, 5); i++) {
            long start = System.nanoTime();
            FactorDecomposition<MultivariatePolynomial<BigInteger>> yun = MultivariateSquareFreeFactorization.SquareFreeFactorizationYunZeroCharacteristics(poly);
            System.out.println("Yun: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            assertFactorization(poly, yun);

            start = System.nanoTime();
            FactorDecomposition<MultivariatePolynomial<BigInteger>> mus = MultivariateSquareFreeFactorization.SquareFreeFactorizationMusserZeroCharacteristics(poly);
            System.out.println("Musser: " + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            assertFactorization(poly, mus);
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

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asLongPolyZp(poly);
        for (int i = 0; i < its(1, 5); i++) {
            long start = System.nanoTime();
            assertFactorization(lPoly, MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser(lPoly));
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

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asLongPolyZp(poly);
        assertFactorization(lPoly, MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser(lPoly));
    }

    @Test
    public void test4() throws Exception {
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("11 + x^7*y^14 + z^7*x^14 + y^7", domain),
                b = MultivariatePolynomial.parse("11 + 3*y^7*z^14 + 4*x^7*y^14 + 5*z^7", domain),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13", domain),
                poly = a.square().multiply(b.square()).multiply(c.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asLongPolyZp(poly);
        FactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser(lPoly);
        assertFactorization(lPoly, decomposition);
    }

    @Test
    public void test5() throws Exception {
        IntegersZp domain = new IntegersZp(2);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("1 + a^6*b^14 + a^2*b^4 + a^7", domain),
                b = MultivariatePolynomial.parse("1 + a^3*b^4 + a + b", domain),
                poly = a.square().multiply(b.square());

        MultivariatePolynomialZp64 lPoly = MultivariatePolynomial.asLongPolyZp(poly);
        FactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser(lPoly);
        assertFactorization(lPoly, decomposition);
    }

    @Test
    public void test5_finiteField() throws Exception {
        long modulus = 2;

        FiniteField<UnivariatePolynomialZp64> field = new FiniteField<>(IrreduciblePolynomials.randomIrreduciblePolynomial(modulus, 4, getRandom()));
        MultivariatePolynomial<UnivariatePolynomialZp64>
                a = MultivariatePolynomial.zero(3, field, LEX)
                .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(1, 2, 3, 4, 5).modulus(modulus)), 1, 1, 3))
                .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 3, 2, 13).modulus(modulus)), 3, 2, 1))
                .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 13, 12, 13).modulus(modulus)), 0, 2, 1)),
                b = MultivariatePolynomial.zero(3, field, LEX)
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5).modulus(modulus)), 1, 1, 13))
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 1, 1, 2, 13).modulus(modulus)), 2, 2, 1))
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 11, 113, 112, 13).modulus(modulus)), 10, 2, 1)),
                c = MultivariatePolynomial.one(3, field, LEX)
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(1, 1, 3, 4, 5, 12).modulus(modulus)), 11, 1, 13))
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(11, 2, 1, 1, 2, 13).modulus(modulus)), 21, 2, 1))
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(modulus)), 10, 12, 1))
                        .add(Monomial.create(field.valueOf(UnivariatePolynomialZ64.create(2, 111, 113, 112, 13, 12).modulus(modulus)), 0, 0, 1)),
                poly = a.square().multiply(b.square()).multiply(c.square());


        FactorDecomposition<MultivariatePolynomial<UnivariatePolynomialZp64>> decomposition = MultivariateSquareFreeFactorization.SquareFreeFactorizationMusser(poly);
        assertFactorization(poly, decomposition);
    }

    @Test
    public void test6() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e", "f", "g", "h"};
        for (Domain<BigInteger> domain : Arrays.<Domain<BigInteger>>asList(new IntegersZp(2), Domains.Z)) {
            MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.parse("a^2*b^4*c*e^5", domain, vars);
            FactorDecomposition<MultivariatePolynomial<BigInteger>> expected = FactorDecomposition.empty(poly);
            expected.addFactor(MultivariatePolynomial.parse("a", domain, vars), 2);
            expected.addFactor(MultivariatePolynomial.parse("b", domain, vars), 4);
            expected.addFactor(MultivariatePolynomial.parse("c", domain, vars), 1);
            expected.addFactor(MultivariatePolynomial.parse("e", domain, vars), 5);

            Assert.assertEquals(expected, SquareFreeFactorization(poly));
        }
    }
}