package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.poly.univar.UnivariateInterpolation.InterpolationZp64;
import org.apache.commons.math3.random.Well44497b;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.PolynomialMethods.Factor;
import static cc.redberry.rings.poly.PolynomialMethods.polyPow;
import static cc.redberry.rings.poly.univar.IrreduciblePolynomials.*;
import static cc.redberry.rings.poly.univar.UnivariateGCD.*;
import static cc.redberry.rings.poly.univar.UnivariateSquareFreeFactorization.SquareFreeFactorization;

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

    @Test
    public void test5() throws Exception {
// Polynomials over field
        UnivariatePolynomialZp64 a = UnivariatePolynomialZ64.create(1, 3, 2).modulus(17);
        UnivariatePolynomialZp64 b = UnivariatePolynomialZ64.create(1, 0, -1).modulus(17);
// Euclid and Half-GCD algorithms for polynomials over field
        assert EuclidGCD(a, b).equals(HalfGCD(a, b));
// Extended Euclidean algorithm
        UnivariatePolynomialZp64[] xgcd = ExtendedEuclidGCD(a, b);
        assert a.copy().multiply(xgcd[1]).add(b.copy().multiply(xgcd[2])).equals(xgcd[0]);
// Extended Half-GCD algorithm
        UnivariatePolynomialZp64[] xgcd1 = ExtendedHalfGCD(a, b);
        assert Arrays.equals(xgcd, xgcd1);


// Polynomials over Z
        UnivariatePolynomial<BigInteger> aZ = UnivariatePolynomial.create(1, 3, 2);
        UnivariatePolynomial<BigInteger> bZ = UnivariatePolynomial.create(1, 0, -1);

// GCD for polynomials over Z
        assert ModularGCD(aZ, bZ).equals(UnivariatePolynomial.create(1, 1));

// Bivariate polynomials represented as Z[y][x]
        UnivariateRing<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>>
                ringXY = UnivariateRing(UnivariateRing(Z));
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>>
                aXY = ringXY.parse("(1 + y) + (1 + y^2)*x + (y - y^2)*x^2"),
                bXY = ringXY.parse("(3 + y) + (3 + 2*y + y^2)*x + (3*y - y^2)*x^2");
//// Subresultant sequence
        PolynomialRemainders<UnivariatePolynomial<UnivariatePolynomial<BigInteger>>>
                subResultants = SubresultantRemainders(aXY, bXY);
// The GCD
        UnivariatePolynomial<UnivariatePolynomial<BigInteger>> gcdXY = subResultants.gcd().primitivePart();
        assert UnivariateDivision.remainder(aXY, gcdXY, true).isZero();
        assert UnivariateDivision.remainder(bXY, gcdXY, true).isZero();
    }

    @Test
    public void test6() throws Exception {


// ring GF(13^5)[x] (coefficient domain is finite field)
        UnivariateRing<UnivariatePolynomial<UnivariatePolynomialZp64>> ringF = UnivariateRing(GF(13, 5));
// some random polynomial composed from some factors
        UnivariatePolynomial<UnivariatePolynomialZp64> polyF = ringF.randomElement().multiply(ringF.randomElement().multiply(polyPow(ringF.randomElement(), 10)));

// perform square-free factorization
        System.out.println(SquareFreeFactorization(polyF));
// perform complete factorization
        System.out.println(Factor(polyF));


// ring Q[x]
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> ringQ = UnivariateRing(Q);
// some random polynomial composed from some factors
        UnivariatePolynomial<Rational<BigInteger>> polyQ = ringQ.randomElement().multiply(ringQ.randomElement().multiply(polyPow(ringQ.randomElement(), 10)));
// perform square-free factorization
        System.out.println(SquareFreeFactorization(polyQ));
// perform complete factorization
        System.out.println(Factor(polyQ));
    }

    @Test
    public void test7() throws Exception {
        Well44497b random = new Well44497b();

// random irreducible polynomial in Z/2[x] of degree 10
        UnivariatePolynomialZp64 poly1 = randomIrreduciblePolynomial(2, 10, random);
        assert poly1.degree() == 10;
        assert irreducibleQ(poly1);

// random irreducible polynomial in Z/2[x] of degree 10
        UnivariatePolynomial<BigInteger> poly2 = randomIrreduciblePolynomial(Zp(2), 10, random);
        assert poly2.degree() == 10;
        assert irreducibleQ(poly2);

// random irreducible polynomial in GF(11^15)[x] of degree 10
        UnivariatePolynomial<UnivariatePolynomialZp64> poly3 = randomIrreduciblePolynomial(GF(11, 15), 10, random);
        assert poly3.degree() == 10;
        assert irreducibleQ(poly3);

// random irreducible polynomial in Z[x] of degree 10
        UnivariatePolynomial<BigInteger> poly4 = randomIrreduciblePolynomialOverZ(10, random);
        assert poly4.degree() == 10;
        assert irreducibleQ(poly4);
    }

    @Test
    public void test8() throws Exception {
// points
long[] points = {1L, 2L, 3L, 12L};
// values
long[] values = {3L, 2L, 1L, 6L};

// interpolate using Newton method
UnivariatePolynomialZp64 result = new InterpolationZp64(Zp64(17))
        .update(points, values)
        .getInterpolatingPolynomial();

// result.evaluate(points(i)) = values(i)
assert IntStream.range(0, points.length).allMatch(i -> result.evaluate(points[i]) == values[i]);
    }
}