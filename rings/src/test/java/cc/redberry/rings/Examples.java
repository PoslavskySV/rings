package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.QuotientRing;
import cc.redberry.rings.poly.UnivariateQuotientRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.*;
import cc.redberry.rings.poly.multivar.GroebnerBases.HilbertSeries;
import cc.redberry.rings.poly.univar.UnivariateDivision;
import cc.redberry.rings.poly.univar.UnivariateInterpolation.InterpolationZp64;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.primes.SmallPrimes;
import org.apache.commons.math3.random.Well44497b;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.IntStream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.PolynomialMethods.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariateGCD.*;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;
import static cc.redberry.rings.poly.univar.IrreduciblePolynomials.irreducibleQ;
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
        System.out.println(poly1.toString("x", "y"));
        // the result will be "y*x + y^2"
        System.out.println(poly1.toString("y", "x"));
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
        //    // Subresultant sequence
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


    /**
     * @param <Poly> polynomial type
     */
    static <Poly extends IPolynomial<Poly>> Poly genericFunc(Poly poly) {
        return poly.createOne().add(
                poly.copy().multiply(2),
                polyPow(poly, 2).multiply(3));
    }

    @Test
    public void test9() throws Exception {
        System.out.println(genericFunc(UnivariatePolynomialZ64.create(1, 2, 3).modulus(17)));
        System.out.println(genericFunc(MultivariatePolynomial.parse("1 + x + y + z")));
    }

    /**
     * @param <Poly> polynomial type
     */
    static <Poly extends IPolynomial<Poly>> Poly genericFuncWithRing(Poly poly, Ring<Poly> ring) {
        return ring.add(
                ring.getOne(),
                ring.multiply(poly, ring.valueOf(2)),
                ring.multiply(ring.pow(poly, 2), ring.valueOf(3)));
    }

    @Test
    public void test10() throws Exception {
        UnivariateRing<UnivariatePolynomialZp64> uRing = UnivariateRingZp64(17);
        System.out.println(genericFuncWithRing(uRing.parse("1 + 2*x + 3*x^2"), uRing));

        MultivariateRing<MultivariatePolynomial<BigInteger>> mRing = MultivariateRing(3, Z);
        System.out.println(genericFuncWithRing(mRing.parse("1 + x + y + z"), mRing));
    }

    @Test
    public void test11() throws Exception {
        FiniteField<UnivariatePolynomialZp64> gf = GF(13, 4);
        UnivariatePolynomialZp64 poly = gf.pow(gf.parse("1 + z + z^2 + z^3 + z^4"), 10);

        UnivariatePolynomialZp64 noRing = genericFunc(poly);
        System.out.println(noRing);

        UnivariatePolynomialZp64 withRing = genericFuncWithRing(poly, gf);
        System.out.println(withRing);

        assert !noRing.equals(withRing);
    }


    /**
     * @param <Monomial> type of monomials
     * @param <Poly>     type of multivariate polynomials
     */
    static <Monomial extends AMonomial<Monomial>,
            Poly extends AMultivariatePolynomial<Monomial, Poly>>
    Poly genericFunc(Poly poly) { return null; }

    /**
     * @param <Monomial> type of monomials
     * @param <Poly>     type of multivariate polynomials
     */
    static <Monomial extends AMonomial<Monomial>,
            Poly extends AMultivariatePolynomial<Monomial, Poly>>
    Poly genericFuncWithRing(Poly poly, IPolynomialRing<Poly> ring) { return null; }

    @Test
    public void test12() throws Exception {
        genericFunc(MultivariatePolynomial.parse("a + b"));

        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        genericFuncWithRing(ring.parse("a + b"), ring);
    }

    @Test
    public void test13() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring
                = MultivariateRing(2, Z, GREVLEX);

        // poly in GREVLEX
        MultivariatePolynomial<BigInteger> poly = ring.parse("x + x^2*y^2 + x*y");
        assert poly.ordering == GREVLEX;

        // poly in LEX
        MultivariatePolynomial<BigInteger> poly2 = poly.setOrdering(MonomialOrder.LEX);
        assert poly2.ordering == MonomialOrder.LEX;

        // poly in GREVLEX (ordering of lhs is used)
        MultivariatePolynomial<BigInteger> add = ring.add(poly, poly2);
        assert add.ordering == GREVLEX;

        // poly in LEX (ordering of lhs is used)
        MultivariatePolynomial<BigInteger> add2 = ring.add(poly2, poly);
        assert add2.ordering == MonomialOrder.LEX;
    }

    @Test
    public void test14() throws Exception {

        String[] variables = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                dividend = MultivariatePolynomial.parse("x - x^2*y^2 + 2*x*y + 1 - z*y^2*x^2 + z", variables),
                divider1 = MultivariatePolynomial.parse("x + y", variables),
                divider2 = MultivariatePolynomial.parse("x + z", variables),
                divider3 = MultivariatePolynomial.parse("y + z", variables);

        dividend = polyPow(dividend, 3);

        {
            MultivariatePolynomial<BigInteger>[] divRem
                    = MultivariateDivision.divideAndRemainder(dividend, divider1, divider2);

            MultivariatePolynomial<BigInteger>
                    quot1 = divRem[0],
                    quot2 = divRem[1],
                    rem = divRem[2];

            assert dividend.equals(rem.copy().add(
                    quot1.copy().multiply(divider1),
                    quot2.copy().multiply(divider2)));
        }

        {
            MultivariatePolynomial<BigInteger>[] divRem
                    = MultivariateDivision.divideAndRemainder(dividend, divider1, divider2, divider3);

            MultivariatePolynomial<BigInteger>
                    quot1 = divRem[0],
                    quot2 = divRem[1],
                    quot3 = divRem[2],
                    rem = divRem[3];

            assert dividend.equals(rem.copy().add(
                    quot1.copy().multiply(divider1),
                    quot2.copy().multiply(divider2),
                    quot3.copy().multiply(divider3)));
        }
    }

    @Test
    public void test15() throws Exception {


        // some large finite field
        IntegersZp64 zpRing = Zp64(SmallPrimes.nextPrime(1 << 15));
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("x^2 - x*y + z^5", zpRing),
                b = MultivariatePolynomialZp64.parse("x^2 + x*y^7 + x*y*z^2", zpRing);

        MultivariatePolynomialZp64
                gcd = MultivariatePolynomialZp64.parse("x + y + z", zpRing),
                poly1 = a.copy().multiply(gcd),
                poly2 = b.copy().multiply(gcd);

        // EZGCD in finite field
        MultivariatePolynomialZp64 ez = EZGCD(poly1, poly2);
        assert ez.equals(gcd);

        // EEZGCD in finite field
        MultivariatePolynomialZp64 eez = EEZGCD(poly1, poly2);
        assert eez.equals(gcd);

        // ZippelGCD in finite field
        MultivariatePolynomialZp64 zippel = ZippelGCD(poly1, poly2);
        assert zippel.equals(gcd);

        // some very small finite field (Z/2)
        IntegersZp64 z2 = Zp64(2);
        MultivariatePolynomialZp64
                z2GCD = gcd.setRing(z2),
                z2Poly1 = a.setRing(z2).multiply(z2GCD),
                z2Poly2 = b.setRing(z2).multiply(z2GCD);

        // Kaltofen’s & Monagan’s generic modular GCD
        MultivariatePolynomialZp64 modGF = MultivariateGCD.KaltofenMonaganSparseModularGCDInGF(z2Poly1, z2Poly2);
        assert modGF.equals(z2GCD);

        // Z
        MultivariatePolynomial<BigInteger>
                zGCD = gcd.setRing(Z),
                zPoly1 = a.setRing(Z).multiply(zGCD),
                zPoly2 = b.setRing(Z).multiply(zGCD);

        // Modular GCD in Z with sparse interpolation
        MultivariatePolynomial<BigInteger> mod = ZippelGCDInZ(zPoly1, zPoly2);
        assert mod.equals(zGCD);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void test16() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        int rndDegree = 5, rndSize = 5;

        // some random gcd
        MultivariatePolynomial<BigInteger> gcd = ring.randomElement(rndDegree, rndSize);
        // array of random polynomials which have gcd
        MultivariatePolynomial<BigInteger>[] polys = IntStream.range(0, 10)
                .mapToObj(i -> ring.randomElement(rndDegree, rndSize).multiply(gcd))
                .toArray(MultivariatePolynomial[]::new);

        // fast algorithm for array of polynomials will be used
        MultivariatePolynomial<BigInteger> fastGCD = MultivariateGCD.PolynomialGCD(polys);
        // slow step-by-step gcd calculation
        MultivariatePolynomial<BigInteger> slowGCD = Arrays.stream(polys)
                .reduce(ring.getZero(), MultivariateGCD::PolynomialGCD);
        // result the same
        assert fastGCD.equals(slowGCD) || fastGCD.equals(slowGCD.negate());
    }

    @Test
    public void test17() throws Exception {
        // ring GF(13^5)[x, y, z] (coefficient domain is finite field)
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>>
                ringF = MultivariateRing(3, GF(13, 5));

        // generate random poly of degree 5 and size 5
        Supplier<MultivariatePolynomial<UnivariatePolynomialZp64>> randomPolyF
                = () -> ringF.randomElement(5, 5).increment();

        // some random polynomial composed from some factors
        MultivariatePolynomial<UnivariatePolynomialZp64> polyF =
                randomPolyF.get().multiply(
                        randomPolyF.get(), ringF.pow(randomPolyF.get(), 2));
        // perform square-free factorization
        System.out.println(FactorSquareFree(polyF));
        // perform complete factorization
        System.out.println(Factor(polyF));


        // ring Q[x, y, z]
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ringQ = MultivariateRing(3, Q);

        Supplier<MultivariatePolynomial<Rational<BigInteger>>> randomPolyQ
                = () -> ringQ.randomElement(5, 5).increment();
        // some random polynomial composed from some factors
        MultivariatePolynomial<Rational<BigInteger>> polyQ =
                randomPolyQ.get().multiply(
                        randomPolyQ.get(), ringQ.pow(randomPolyQ.get(), 2));
        // perform square-free factorization
        System.out.println(FactorSquareFree(polyQ));
        // perform complete factorization
        System.out.println(Factor(polyQ));
    }

    @Test
    public void test18() throws Exception {

        // ring GF(13^6)[x, y, z]
        FiniteField<UnivariatePolynomialZp64> cfRing = GF(13, 6);
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> ring = MultivariateRing(3, cfRing);


        UnivariatePolynomialZp64[] points = {
                cfRing.parse("1 + t"),
                cfRing.parse("2 + t"),
                cfRing.parse("3 + t"),
                cfRing.parse("12 + t")
        };

        String[] vars = {"x", "y", "z"};
        // some values for interpolation
        MultivariatePolynomial[] values = {
                ring.parse("x + y", vars),
                ring.parse(" x^2 + (t) * y", vars),
                ring.parse("y^3", vars),
                ring.parse("(t) * x^4 + y", vars)
        };

        // interpolation polynomial values for variable z
        MultivariatePolynomial<UnivariatePolynomialZp64> result =
                new MultivariateInterpolation.Interpolation(2, ring)
                        .update(points, values)
                        .getInterpolatingPolynomial();

        assert IntStream.range(0, points.length).allMatch(i -> result.evaluate(2, points[i]).equals(values[i]));
    }

    @Test
    public void test19() throws Exception {

        // Galois field GF(7^10)
        // (irreducible polynomial will be generated automatically)
        FiniteField<UnivariatePolynomialZp64> gf7_10 = GF(7, 10);
        assert gf7_10.characteristic().intValue() == 7;
        assert gf7_10.cardinality().equals(BigInteger.valueOf(7).pow(10));

        // GF(7^3) generated by irreducible polynomial "1 + 3*z + z^2 + z^3"
        FiniteField<UnivariatePolynomialZp64> gf7_3 = GF(UnivariatePolynomialZ64.create(1, 3, 1, 1).modulus(7));
        assert gf7_3.characteristic().intValue() == 7;
        assert gf7_3.cardinality().intValue() == 7 * 7 * 7;
    }

    @Test
    public void test20() throws Exception {
        // ring Z
        Ring<BigInteger> notField = Z;
        // it is not a fielf
        assert !notField.isField();
        // this is OK
        assert 1 == notField.reciprocal(Z.valueOf(1)).intValue();
        // this will throw ArithmeticException
        notField.reciprocal(Z.valueOf(10));
    }

    @Test
    public void test21() throws Exception {
        int intNumber = 1234567;
        // false
        boolean primeQ = SmallPrimes.isPrime(intNumber);
        // 1234577
        int intPrime = SmallPrimes.nextPrime(intNumber);
        // [127, 9721]
        int[] intFactors = SmallPrimes.primeFactors(intNumber);

        long longNumber = 12345671234567123L;
        // false
        primeQ = BigPrimes.isPrime(longNumber);
        // 12345671234567149
        long longPrime = BigPrimes.nextPrime(longNumber);
        // [1323599, 9327350077]
        long[] longFactors = BigPrimes.primeFactors(longNumber);

        BigInteger bigNumber = Z.parse("321536584276145124691487234561832756183746531874567");
        // false
        primeQ = BigPrimes.isPrime(bigNumber);
        // 321536584276145124691487234561832756183746531874827
        BigInteger bigPrime = BigPrimes.nextPrime(bigNumber);
        // [3, 29, 191, 797359, 1579057, 14916359, 1030298906727233717673336103]
        List<BigInteger> bigFactors = BigPrimes.primeFactors(bigNumber);
    }

    @Test
    public void test22() throws Exception {
        Rationals<BigInteger> field = Frac(Z);     // the same as Q

        Rational<BigInteger>
                a = field.parse("13/6"),
                b = field.parse("2/3"),
                c = field.parse("3/2");

        assert field.parse("13/6")
                .equals(field.add(
                        field.parse("2/3"),
                        field.parse("3/2")));

        assert field.parse("5/6")
                .equals(field.add(
                        field.parse("2/3"),
                        field.parse("1/6")));

    }

    @Test
    public void test23() throws Exception {
        // Ring Z/3[x]
        UnivariateRing<UnivariatePolynomialZp64> zp3x = UnivariateRingZp64(3);
        // parse univariate poly from string
        UnivariatePolynomialZp64
                p1 = zp3x.parse("4 + 8*x + 13*x^2"),
                p2 = zp3x.parse("4 - 8*x + 13*x^2");
        assert zp3x.add(p1, p2).equals(zp3x.parse("2 - x^2"));


        // GF(7^3)
        FiniteField<UnivariatePolynomialZp64> cfRing = GF(UnivariateRingZp64(7).parse("1 + 3*z + z^2 + z^3"));
        // GF(7^3)[x]
        UnivariateRing<UnivariatePolynomial<UnivariatePolynomialZp64>> gfx = UnivariateRing(cfRing);
        // parse univariate poly from string
        UnivariatePolynomial<UnivariatePolynomialZp64>
                r1 = gfx.parse("4 + (8 + z)*x + (13 - z^43)*x^2"),
                r2 = gfx.parse("4 - (8 + z)*x + (13 + z^43)*x^2");
        assert gfx.add(r1, r2).equals(gfx.parse("1 - 2*x^2"));
        UnivariatePolynomial<UnivariatePolynomialZp64>
                divRem[] = divideAndRemainder(r1, r2),
                div = divRem[0],
                rem = divRem[1];
        assert r1.equals(gfx.add(gfx.multiply(r2, div), rem));
    }

    @Test
    public void test24() throws Exception {
        String[] vars = {"x", "y", "z"};
        // Ring Z/3[x, y, z]
        MultivariateRing<MultivariatePolynomialZp64> zp3xyz = MultivariateRingZp64(3, 3);
        // parse univariate poly from string
        MultivariatePolynomialZp64
                p1 = zp3xyz.parse("4 + 8*x*y + 13*x^2*z^5", vars),
                p2 = zp3xyz.parse("4 - 8*x*y + 13*x^2*z^5", vars);
        assert zp3xyz.add(p1, p2).equals(zp3xyz.parse("2 - x^2*z^5", vars));


        // GF(7^3)
        FiniteField<UnivariatePolynomialZp64> cfRing = GF(UnivariateRingZp64(7).parse("1 + 3*z + z^2 + z^3"));
        // GF(7^3)[x]
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> gfxyz = MultivariateRing(3, cfRing);
        // parse univariate poly from string
        MultivariatePolynomial<UnivariatePolynomialZp64>
                r1 = gfxyz.parse("4 + (8 + z)*x*y + (13 - z^43)*x^2*z^5", vars),
                r2 = gfxyz.parse("4 - (8 + z)*x*y + (13 + z^43)*x^2*z^5", vars);
        assert gfxyz.add(r1, r2).equals(gfxyz.parse("1 - 2*x^2*z^5", vars));
        MultivariatePolynomial<UnivariatePolynomialZp64>
                divRem[] = divideAndRemainder(r1, r2),
                div = divRem[0],
                rem = divRem[1];
        assert r1.equals(gfxyz.add(gfxyz.multiply(r2, div), rem));
    }

    @Test
    public void test25() throws Exception {
        FiniteField<UnivariatePolynomialZp64> ring = GF(17, 9);

        UnivariatePolynomialZp64 a = ring.randomElement();
        UnivariatePolynomialZp64 b = ring.pow(a, 1000);
        UnivariatePolynomialZp64 c = ring.reciprocal(b);

        assert ring.multiply(b, c).isOne();

        UnivariatePolynomialZp64 some = ring.add(
                ring.divideExact(a, ring.add(b, c)),
                ring.pow(a, 6),
                ring.negate(ring.multiply(a, b, c)));
    }

    @Test
    public void test26() throws Exception {
        // Z[x, y, z]
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z, MonomialOrder.LEX);

        MultivariatePolynomial<BigInteger>
                x = ring.variable(0),
                y = ring.variable(1),
                z = ring.variable(2);

        // do some math
        MultivariatePolynomial<BigInteger> a = ring.decrement(ring.pow(ring.add(x, y, z), 2));
        MultivariatePolynomial<BigInteger> b = ring.add(
                ring.pow(ring.add(x, ring.negate(y), ring.negate(z), ring.getNegativeOne()), 2),
                x, y, z, ring.getNegativeOne());
        MultivariatePolynomial<BigInteger> c = ring.add(
                ring.pow(ring.add(a, b, ring.getOne()), 9),
                ring.negate(a), ring.negate(b), ring.getNegativeOne());

        // reduce c modulo a and b (multivariate division with remainder)
        MultivariatePolynomial<BigInteger>[] divRem = MultivariateDivision.divideAndRemainder(c, a, b);
        MultivariatePolynomial<BigInteger>
                div1 = divRem[0],
                div2 = divRem[1],
                rem = divRem[2];
    }

    @Test
    public void test27() throws Exception {
        // base ring
        UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> baseRing = UnivariateRing(Q);
        // poly in base ring
        UnivariatePolynomial<Rational<BigInteger>> basePoly = baseRing.parse("123 * x^31 + 123 * x^2 + (1/2) * x + 1");

        UnivariatePolynomial<Rational<BigInteger>> modulus = baseRing.parse("x^2 + 1");
        // quotient ring
        UnivariateQuotientRing<UnivariatePolynomial<Rational<BigInteger>>> quotRing = UnivariateQuotientRing(baseRing, modulus);
        // same poly in quotient ring
        UnivariatePolynomial<Rational<BigInteger>> quotPoly = quotRing.parse("123 * x^31 + 123 * x^2 + (1/2) * x + 1");

        assert basePoly.degree() == 31;
        assert quotPoly.degree() == 1;
        assert quotPoly.equals(remainder(basePoly, modulus));
    }

    @Test
    public void test38() throws Exception {
        // base ring Q[x,y,z]
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> baseRing = MultivariateRing(3, Q);

        MultivariatePolynomial<Rational<BigInteger>>
                generator1 = baseRing.parse("x^2 + y^12 - z"),
                generator2 = baseRing.parse("x^2*z + y^2 - 1");
        // ideal in a base ring generated by two polys <x^2 + y^12 - z, x^2*z + y^2 - 1>
        // a proper Groebner basis will be constructed automatically
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                ideal = Ideal.create(Arrays.asList(generator1, generator2));
        // quotient ring Q[x,y,z]/I
        QuotientRing<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                quotRing = QuotientRing(baseRing, ideal);

        // do some math in a quotient ring
        MultivariatePolynomial<Rational<BigInteger>>
                q1 = quotRing.parse("10 * x^12 + 11 * y^11 + 12 * z^10"),
                q2 = quotRing.parse("x * y - y * z - z * x"),
                polyQuot = quotRing.add(
                        quotRing.multiply(q1, 11),
                        quotRing.multiply(q1, q1, q2));

        // do the same math in a base ring
        MultivariatePolynomial<Rational<BigInteger>>
                b1 = baseRing.parse("10 * x^12 + 11 * y^11 + 12 * z^10"),
                b2 = baseRing.parse("x * y - y * z - z * x"),
                polyBase = baseRing.add(
                        baseRing.multiply(b1, 11),
                        baseRing.multiply(b1, b1, b2));

        assert !polyQuot.equals(polyBase);
        assert polyQuot.equals(ideal.normalForm(polyBase));
    }

    @Test
    public void test39() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(3, 17);

        // create ideal with two generators using GREVLEX monomial order for underlying Groebner basis
        Ideal<MonomialZp64, MultivariatePolynomialZp64> I = Ideal.create(Arrays.asList(
                ring.parse("x^2 + y^12 - z"),
                ring.parse("x^2 * z + y^2 - 1")), GREVLEX);
        // I is proper ideal
        assert I.isProper();

        // get computed Groebner basis
        List<MultivariatePolynomialZp64> gb = I.getGroebnerBasis();
        System.out.println(gb);

        // check some ideal properties
        assert I.dimension() == 1;
        assert I.degree() == 36;

        // create another ideal with only one generator
        Ideal<MonomialZp64, MultivariatePolynomialZp64> J = Ideal.create(Arrays.asList(
                ring.parse("x^4 * y^4 + 1")), GREVLEX);
        // J is principal ideal
        assert J.isPrincipal();
        assert J.dimension() == 2;
        assert J.degree() == 8;


        Ideal<MonomialZp64, MultivariatePolynomialZp64> union = I.union(J);
        // union is zero dimensional ideal
        assert union.dimension() == 0;
        // change order to LEX (elimination order)
        Ideal<MonomialZp64, MultivariatePolynomialZp64> eliminated = union.changeOrder(LEX);
        // system can now be solved easily
        System.out.println(eliminated);


        Ideal<MonomialZp64, MultivariatePolynomialZp64> intersection = I.intersection(J);
        // intersection is still 2-dimensional
        assert intersection.dimension() == 2;
        // multiplication in this case is equal to intersection
        Ideal<MonomialZp64, MultivariatePolynomialZp64> times = I.multiply(J);
        assert times.equals(intersection);


        // yet another ideal
        Ideal<MonomialZp64, MultivariatePolynomialZp64> K = Ideal.create(Arrays.asList(
                ring.parse("z * x^4 - z * y^14 + y * z^16"),
                ring.pow(ring.parse("x + y + z"), 4)), GREVLEX);
        // compute complicated quotient ideal
        Ideal<MonomialZp64, MultivariatePolynomialZp64> quot = (I.multiply(J).multiply(K)).quotient(times);
        assert quot.equals(K);
    }

    @Test
    public void test40() throws Exception {

        // base ring in LEX order
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> ring = MultivariateRing(4, Q, LEX);
        String[] variables = {"x", "y", "z", "t"};

        // some polynomial in a base ring order (LEX)
        MultivariatePolynomial<Rational<BigInteger>> poly =
                ring.parse("x + y^2 * z + z^3 * y * t + t^4 * z * y", variables);
        assert poly.ordering == LEX;

        // some ideal with Groebner basis computed in GREVLEX
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                idealGrevLex = Ideal.create(Arrays.asList(
                ring.parse("y * x^3 + z * t^3 - 1", variables),
                ring.parse("x * y - y * z - z * x + t^3", variables)),
                GREVLEX);
        assert idealGrevLex.ordering == GREVLEX;

        // normal form of poly will be computed with respect to GREVLEX
        // then the result will be re-sorted according to the base ring order (LEX)
        MultivariatePolynomial<Rational<BigInteger>> nfGrevLex = idealGrevLex.normalForm(poly);
        assert nfGrevLex.ordering == LEX;

        // the same ideal with Groebner basis in LEX order
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                idealLex = idealGrevLex.changeOrder(LEX);
        assert idealLex.ordering == LEX;

        // normal form of poly will be computed with respect to LEX
        MultivariatePolynomial<Rational<BigInteger>> nfLex = idealLex.normalForm(poly);
        assert nfLex.ordering == LEX;

        // Normal forms computed against LEX basis and GREVLEX basis
        // are different (although both polynomials are sorted in LEX)
        assert !nfGrevLex.equals(nfLex);
    }

    @Test
    public void test41() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> ring = MultivariateRingZp64(3, 32003);
        Ideal<MonomialZp64, MultivariatePolynomialZp64> ideal
                = Ideal.create(Arrays.asList(ring.parse("x^2 - x*y"), ring.parse("y^2 - z - 1"), ring.parse("z^2 + x^4")));
        // get Hilbert-Poincare series
        HilbertSeries hps = ideal.hilbertSeries();

        assert hps.dimension() == 0;
        assert hps.degree() == 8;

        // series numerator
        System.out.println(hps.initialNumerator);
        // reduced series numerator
        System.out.println(hps.numerator);

        // integer Hilbert polynomial
        System.out.println(hps.hilbertPolynomialZ());
        // rational Hilbert polynomial
        System.out.println(hps.hilbertPolynomial());
    }

    @Test
    public void test42() throws Exception {
        // some ideal in Z[x,y,z] with very simple Groebner basis
        String[] vars = {"x", "y", "z"};
        MultivariatePolynomial<BigInteger>
                a = parse("8*x^2*y^2 + 5*x*y^3 + 3*x^3*z + x^2*y*z", Z, vars),
                b = parse("x^5 + 2*y^3*z^2 + 13*y^2*z^3 + 5*y*z^4", Z, vars),
                c = parse("8*x^3 + 12*y^3 + x*z^2 + 3", Z, vars),
                d = parse("7*x^2*y^4 + 18*x*y^3*z^2 + y^3*z^3", Z, vars);
        List<MultivariatePolynomial<BigInteger>> gens = Arrays.asList(a, b, c, d);


        // The default method will use modular algorithm in this case
        List<MultivariatePolynomial<BigInteger>> gb = GroebnerBases.GroebnerBasis(gens, GREVLEX);
        // Groebner bases is very simple: <x, z^2, 1 + 4*y^3>
        System.out.println(gb);

        // Modular algorithm will take few milliseconds
        List<MultivariatePolynomial<BigInteger>> mod = GroebnerBases.ModularGB(gens, GREVLEX);
        assert (mod.equals(gb));

        // F4 algorithm will also take few milliseconds
        List<MultivariatePolynomial<BigInteger>> f4 = GroebnerBases.F4GB(gens, GREVLEX);
        assert (f4.equals(gb));

        // But Buchberger algorithm will take several minutes
        List<MultivariatePolynomial<BigInteger>> buch = GroebnerBases.BuchbergerGB(gens, GREVLEX);
        assert (buch.equals(gb));
    }

    @Test
    public void test43() {
        // polynomial ring Z[x,y,z]
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);

        // Coder for Z[x,y,z]
        Coder<MultivariatePolynomial<BigInteger>, ?, ?> coder = Coder.mkMultivariateCoder(ring, "x", "y", "z");

        // parse some element from string
        MultivariatePolynomial<BigInteger> p = coder.parse("x^2 + y^2 + z^2");

        // add special string binding
        coder.bind("p", p);

        // parse another poly
        MultivariatePolynomial<BigInteger> pSquared = coder.parse("p^2");
        assert pSquared.equals(p.square());
    }

    @Test
    public void test44() {
        // Parser for rationals
        Coder<Rational<BigInteger>, ?, ?> qCoder = Coder.mkCoder(Q);
        // parse some rational number
        Rational<BigInteger> el = qCoder.parse("1/2/3 + (1-3/5)^3 + 1");
        System.out.println(el);
    }

    @Test
    public void test45() {
        // univariate ring Z/2[t]
        UnivariateRing<UnivariatePolynomialZp64> uRing = UnivariateRingZp64(2);
        // coder for polynomials from Z/2[t]
        Coder<UnivariatePolynomialZp64, ?, ?> uCoder = Coder.mkUnivariateCoder(uRing, "t");

        // rational functions over Z/2[t]
        Rationals<UnivariatePolynomialZp64> cfRing = Frac(uRing);
        // coder for rational functions from Frac(Z/2[t])
        Coder<Rational<UnivariatePolynomialZp64>, ?, ?>
                cfCoder = Coder.mkRationalsCoder(cfRing, uCoder);

        // ring Frac(Z/2[t])[a,b,c]
        MultivariateRing<MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>>
                ring = MultivariateRing(3, cfRing);
        // coder for polynomials from Frac(Z/2[t])[a,b,c]
        Coder<MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>, ?, ?>
                coder = Coder.mkMultivariateCoder(ring, cfCoder, "a", "b", "c");

        // parse some element
        MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>
                el = coder.parse("(1 + t)*a^2 - c^3 + b/t^2 + (a + b)/(1 + t)^3");

        System.out.println(coder.stringify(el));

        // associate variable "E" with polynomial el in parser
        coder.bind("E", el);

        // "E" will be replaced with el by the parser
        MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>
                el2 = coder.parse("(a+b) * E^2 + 1");

    }

    @Test
    public void test46() {
        // univariate ring Z/2[t]
        UnivariateRing<UnivariatePolynomialZp64> uRing = UnivariateRingZp64(2);
        // coder for polynomials from Z/2[t]
        Coder<UnivariatePolynomialZp64, ?, ?> uCoder = Coder.mkUnivariateCoder(uRing, "t");

        // rational functions over Z/2[t]
        Rationals<UnivariatePolynomialZp64> cfRing = Frac(uRing);
        // coder for rational functions from Frac(Z/2[t])
        Coder<Rational<UnivariatePolynomialZp64>, ?, ?>
                cfCoder = Coder.mkRationalsCoder(cfRing, uCoder);

        // ring Frac(Z/2[t])[a,b,c]
        MultivariateRing<MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>>
                ring = MultivariateRing(3, cfRing);
        // coder for polynomials from Frac(Z/2[t])[a,b,c]
        Coder<MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>, ?, ?>
                coder = Coder.mkMultivariateCoder(ring, cfCoder, "a", "b", "c");

        // parse some element
        MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>
                el = coder.parse("(1 + t)*a^2 - c^3 + b/t^2 + (a + b)/(1 + t)^3");

        System.out.println(coder.stringify(el));

        // associate variable "E" with polynomial el in parser
        coder.bind("E", el);

        // "E" will be replaced with el by the parser
        MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>
                el2 = coder.parse("(a+b) * E^2 + 1");
    }

    @Test
    public void test47() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(3, Z);
        Rationals<MultivariatePolynomial<BigInteger>> field = Frac(ring);

        // Parser/stringifier of rational functions
        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?> coder = Coder.mkRationalsCoder(field,
                Coder.mkMultivariateCoder(ring, "x", "y", "z"));

        // parse some math expression from string
        // it will be automatically reduced to a common denominator
        // with the gcd being automatically cancelled
        Rational<MultivariatePolynomial<BigInteger>> expr1 = coder.parse("(x/y/(x - z) + (x + z)/(y - z))^2 - 1");

        // do some math ops programmatically
        Rational<MultivariatePolynomial<BigInteger>>
                x = new Rational<>(ring, ring.variable(0)),
                y = new Rational<>(ring, ring.variable(1)),
                z = new Rational<>(ring, ring.variable(2));

        Rational<MultivariatePolynomial<BigInteger>> expr2 = field.add(
                field.pow(expr1, 2),
                field.divideExact(x, y),
                field.negate(z));


        // bind expr1 and expr2 to variables to use them further in parser
        coder.bind("expr1", expr1);
        coder.bind("expr2", expr2);

        // parse some complicated expression from string
        // it will be automatically reduced to a common denominator
        // with the gcd being automatically cancelled
        Rational<MultivariatePolynomial<BigInteger>> expr3 = coder.parse(
                " expr1 / expr2 - (x*y - z)/(x-y)/expr1"
                        + " + x / expr2 - (x*z - y)/(x-y)/expr1/expr2"
                        + "+ x^2*y^2 - z^3 * (x - y)^2");

        // export expression to string
        System.out.println(coder.stringify(expr3));

        // take numerator and denominator
        MultivariatePolynomial<BigInteger> num = expr3.numerator();
        MultivariatePolynomial<BigInteger> den = expr3.denominator();

        // common GCD is always cancelled automatically
        assert field.ring.gcd(num, den).isOne();

        // compute unique factor decomposition of expression
        FactorDecomposition<Rational<MultivariatePolynomial<BigInteger>>> factors = field.factor(expr3);
        System.out.println(factors.toString(coder));
    }
}