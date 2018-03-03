package cc.redberry.rings.poly.multivar;

import cc.redberry.combinatorics.Combinatorics;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import org.junit.Test;

import java.util.Arrays;

import static cc.redberry.rings.poly.multivar.MultivariateConversions.*;
import static org.junit.Assert.assertEquals;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateConversionsTest extends AMultivariateTest {
    @Test
    public void test1() throws Exception {
        MultivariatePolynomial<BigInteger>
                poly = MultivariatePolynomial.parse("x1*x2 + 2*x2*x3 + 3*x3*x4 + 4*x4*x5 + 5*x5*x6");
        for (int i = 1; i < poly.nVariables; i++)
            for (int[] c : Combinatorics.combinationsWithPermutations(poly.nVariables, i)) {
                MultivariatePolynomial<MultivariatePolynomial<BigInteger>> split = split(poly, c);
                assertEquals(poly, merge(split, c));
            }
    }

    @Test
    public void test2() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("x1*x2 + 12*x2*x3 + 3*x3*x4^5 + 4*x4*x5 + 5*x5*x6"),
                b = MultivariatePolynomial.parse("x1*x2 + 2*x2^2*x3 + 3*x3^6*x4^7 + 41*x4*x5 + 1*x5*x6");

        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(a);

        int var = 3;
        IPolynomialRing<UnivariatePolynomial<MultivariatePolynomial<BigInteger>>> uRing = MultivariateConversions.asUnivariate(ring, var);

        Rationals<MultivariatePolynomial<BigInteger>> frac = Rings.Frac(uRing.factory().ring);
        UnivariatePolynomial<Rational<MultivariatePolynomial<BigInteger>>>
                ra = asUnivariate(a, var).mapCoefficients(frac, p -> new Rational<>(frac.ring, p)),
                rb = asUnivariate(b, var).mapCoefficients(frac, p -> new Rational<>(frac.ring, p));

        System.out.println(Arrays.toString(PolynomialMethods.divideAndRemainder(rb, ra)));
    }
}