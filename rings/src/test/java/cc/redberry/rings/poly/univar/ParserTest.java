package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FiniteField;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * @since 1.0
 */
public class ParserTest extends AUnivariateTest {
    @Test
    public void test1() throws Exception {
        assertEquals(
                UnivariatePolynomial.create(3),
                Parser.parse(Rings.Z, "3"));
        assertEquals(
                UnivariatePolynomial.create(3, 0, 0, 1, 0, 4),
                Parser.parse(Rings.Z, "3 + x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, 1, 0, 4),
                Parser.parse(Rings.Z, "x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, 1, 0, 4),
                Parser.parse(Rings.Z, "0*x + x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, -1, 0, 4),
                Parser.parse(Rings.Z, "0*x - x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, -1, 0, 4),
                Parser.parse(Rings.Z, "0*x + - x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, -1, 0, 4),
                Parser.parse(Rings.Z, "0*x - + x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, -1, 0, 4),
                Parser.parse(Rings.Z, "-x^3 + 4*x^5"));
        assertEquals(
                UnivariatePolynomial.create(0, 0, 0, -1, 0, -4),
                Parser.parse(Rings.Z, "-x^3 - 4*x^5"));
    }

    @Test
    public void test2() throws Exception {
        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(20, 4, 48, 86, -25, 93, 91, 93, -3, 3, 38, 20, 38, 40, 7);
        assertEquals(poly, Parser.parse(Rings.Z, poly.toString()));
    }

    @Test
    public void test3() throws Exception {
        UnivariatePolynomial<BigInteger> poly = UnivariatePolynomial.create(1, 97, -92, 43);
        UnivariatePolynomial<BigInteger> r = Parser.parse(Rings.Z, poly.toString());
        assertEquals(poly, r);
    }

    @Test
    public void testParseRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(1000, 1000); i++) {
            UnivariatePolynomialZ64 lPoly = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(1, 20), rnd);
            UnivariatePolynomial<BigInteger> bPoly = lPoly.toBigPoly();

            for (String str : new String[]{lPoly.toString(), bPoly.toString()})
                assertEquals(bPoly, Parser.parse(Rings.Z, str));
        }
    }

    @Test
    public void testParseRandom2() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(100, 1000); i++) {
            UnivariatePolynomial<UnivariatePolynomialZp64> poly = RandomUnivariatePolynomials.randomPoly(rndd.nextInt(1, 20), FiniteField.GF17p5, rnd);
            assertEquals(poly, Parser.parse(poly.ring, poly.toString()));
        }
    }

    @Test
    public void testFiniteField1() throws Exception {
        UnivariatePolynomialZp64
                c0 = UnivariatePolynomialZ64.create(1, 2, 3).modulus(17),
                c1 = UnivariatePolynomialZ64.create(1).modulus(17),
                c2 = UnivariatePolynomialZ64.create(0, 2, 3).modulus(17),
                c3 = UnivariatePolynomialZ64.create(-1, -2, -3).modulus(17);

        UnivariatePolynomial<UnivariatePolynomialZp64> poly =
                UnivariatePolynomial.create(FiniteField.GF17p5, c0, c1, c2, c3);
        assertEquals(poly, Parser.parse(poly.ring, poly.toString()));
    }

    @Test
    public void test4() throws Exception {
        System.out.println(UnivariatePolynomialZp64.parse("12312341231412423142342343125234234321423 + x", new IntegersZp64(3)));
    }
}