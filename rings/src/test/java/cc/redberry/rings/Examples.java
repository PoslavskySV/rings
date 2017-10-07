package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariateDivision;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import org.junit.Test;

import java.util.Arrays;

import static cc.redberry.rings.Rings.UnivariateRingZp64;
import static cc.redberry.rings.Rings.Z;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class Examples {

    @Test
    public void name() throws Exception {
// when parsing "x" will be considered as the "first variable"
// and "y" as "the second" => in the result the particular
// names "x" and "y" are erased
        MultivariatePolynomial<BigInteger> poly1 = MultivariatePolynomial.parse("x^2 + x*y", "x", "y");
// parse the same polynomial but using "a" and "b" instead of "x" and "y"
        MultivariatePolynomial<BigInteger> poly2 = MultivariatePolynomial.parse("a^2 + a*b", "a", "b");
// polynomials are equal (no matter which variable names were used when parsing)
        assert poly1.equals(poly2);
// degree in the first variable
        assert poly1.degree(0) == 2;
// degree in the second variable
        assert poly1.degree(1) == 1;

// this poly differs from poly2 since now "a" is "the second"
// variable and "b" is "the first"
        MultivariatePolynomial<BigInteger> poly3 = MultivariatePolynomial.parse("a^2 + a*b", "b", "a");
        assert !poly3.equals(poly2);
// swap the first and the second variables and the result is equal to poly2
        assert AMultivariatePolynomial.swapVariables(poly3, 0, 1).equals(poly2);


// the default toString() will use the default
// variables "a", "b", "c"  and so on (alphabetical)
// the result will be "a*b + a^2"
        System.out.println(poly1);
// specify which variable names use for printing
// the result will be "x*y + x^2"
        System.out.println(poly1.toString(new String[]{"x", "y"}));
// the result will be "y*x + y^2"
        System.out.println(poly1.toString(new String[]{"y", "x"}));
    }

    @Test
    public void test3() throws Exception {
// parse polynomials
        UnivariatePolynomial
                p1 = UnivariatePolynomial.parse("x", Z),
                p2 = UnivariatePolynomial.parse("x^2", Z),
                p3 = UnivariatePolynomial.parse("x^3", Z);

// this WILL modify poly1
        p1.add(p2);
// this will NOT modify poly2
        p2.copy().add(p3);

    }

    @Test
    public void test4() throws Exception {
UnivariateRing<UnivariatePolynomialZp64> ring = UnivariateRingZp64(17);
// some random divider
UnivariatePolynomialZp64 divider = ring.randomElement();
// some random dividend
UnivariatePolynomialZp64 dividend = ring.add(
        ring.valueOf(1),
        ring.multiply(divider, ring.valueOf(2)),
        ring.pow(divider, 2));

// quotient and remainder using built-in methods
UnivariatePolynomialZp64[] divRemPlain
        = UnivariateDivision.divideAndRemainder(dividend, divider, true);

// precomputed Newton inverses, need to calculate it only once
UnivariateDivision.InverseModMonomial<UnivariatePolynomialZp64> invMod
        = UnivariateDivision.fastDivisionPreConditioning(divider);
// quotient and remainder computed using fast
// algorithm with precomputed Newton inverses

UnivariatePolynomialZp64[] divRemFast
        = UnivariateDivision.divideAndRemainderFast(dividend, divider, invMod, true);

// results are the same
assert Arrays.equals(divRemPlain, divRemFast);
    }
}
