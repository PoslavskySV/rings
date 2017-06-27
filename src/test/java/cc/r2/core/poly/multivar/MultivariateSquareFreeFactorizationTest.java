package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.util.TimeUnits;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.FactorDecompositionTest.assertFactorization;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateSquareFreeFactorizationTest extends AbstractPolynomialTest {

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
        IntegersModulo domain = new IntegersModulo(7);
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("11 + x^7*y^14 + z^7*x^14 + y^7", domain),
                b = MultivariatePolynomial.parse("11 + 3*y^7*z^14 + 4*x^7*y^14 + 5*z^7", domain),
                c = MultivariatePolynomial.parse("z*y^2*x^2 - 2*y^3*x - 1234*z^7*x^12*y^13", domain),
                poly = a.square().multiply(b.square()).multiply(c.square());

        lMultivariatePolynomialZp lPoly = MultivariatePolynomial.asLongPolyZp(poly);
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
}