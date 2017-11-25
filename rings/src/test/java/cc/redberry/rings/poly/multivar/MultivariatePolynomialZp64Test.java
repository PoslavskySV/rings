package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class MultivariatePolynomialZp64Test extends AMultivariateTest {

    private IntegersZp64 Integers = new IntegersZp64(Integer.MAX_VALUE);

    private MultivariatePolynomialZp64 parse0(String str) {
        return MultivariatePolynomialZp64.parse(str, Integers);
    }

    private MultivariatePolynomialZp64 parse0(String str, String... vars) {
        return MultivariatePolynomialZp64.parse(str, Integers, vars);
    }

    @Test
    public void testArithmetics1() throws Exception {
        MultivariatePolynomialZp64 a = parse("a*b + a^2 + c^3*b^2", Integers);
        assertEquals(0, a.cc());
        assertEquals(1, a.lc());
        assertEquals(1, a.clone().increment().cc());
        assertEquals(Integers.modulus(-1), a.clone().decrement().cc());
        MultivariatePolynomialZp64 b = parse("a*b - a^2 + c^3*b^2", Integers);
        assertEquals(parse("2*a^2", Integers, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", Integers, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", Integers), a.multiply(b));
    }

    @Test
    public void testZero1() throws Exception {
        MultivariatePolynomialZp64 p = parse("a*b + a^2 + c^3*b^2", Integers);

        MultivariatePolynomialZp64 a = p.clone();
        a.subtract(a.clone());
        assertZero(a);

        a = p.clone();
        a.subtract(a.clone());
        assertZero(a);

        a = p.clone();
        a.add(a.clone().negate());
        assertZero(a);

        a = p.clone().toZero();
        a.add(a.clone().negate());
        assertZero(a);

        a = p.clone().multiply(p.createZero());
        assertZero(a);

        a = p.createZero();
        a.subtractLt();
        assertZero(a);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testZero2() throws Exception {
        MultivariatePolynomialZp64 poly = MultivariatePolynomialZp64.create(3,
                Integers, LEX,
                new MonomialZp64(new int[]{1, 2, 3}, 0),
                new MonomialZp64(new int[]{0, 1, 2}, 5),
                new MonomialZp64(new int[]{0, 1, 2}, 5),
                new MonomialZp64(new int[]{3, 43, 1}, 10));
        assertEquals(2, poly.size());
        assertEquals(parse("10*b*c^2 + 10*a^3*b^43*c", Integers), poly);
    }

    private static void assertZero(MultivariatePolynomialZp64 a) {
        assertEquals(0, a.size());
        assertTrue(a.isZero());
        assertTrue(a.isConstant());
        assertTrue(a.lt().isZeroVector());
        assertEquals(a.lc(), 0);
        assertEquals(a.cc(), 0);
        assertEquals(a.degree(), 0);
    }

    @Test
    public void testArithmetics2() throws Exception {
        IntegersZp64 algebra = new IntegersZp64(17);
        MultivariatePolynomialZp64 a = parse("a*b + a^2 + c^3*b^2", algebra, LEX);
        assertEquals(0, a.cc());
        assertEquals(1, a.lc());
        assertEquals(1, a.clone().increment().cc());
        assertEquals(algebra.modulus(-1), a.clone().decrement().cc());
        MultivariatePolynomialZp64 b = parse("a*b - a^2 + c^3*b^2", algebra, LEX);
        assertEquals(parse("2*a^2", algebra, LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", algebra, LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", algebra, LEX), a.multiply(b));
    }

    @Test
    public void testZeroVariables() throws Exception {
        MultivariatePolynomialZp64 p0 = parse("23", Integers);
        assertEquals(0, p0.nVariables);
        assertEquals(1, p0.size());
        assertEquals(0, p0.clone().subtract(p0).size());
    }

    @Test
    public void testCreateLinear() throws Exception {
        MultivariatePolynomialZp64 p0 = MultivariatePolynomialZp64.zero(3, Integers, LEX);
        String[] vars = {"a", "b", "c"};
        assertEquals(parse("-1+2*a", Integers, vars), p0.createLinear(0, -1, 2));
        assertEquals(parse("-1+2*b", Integers, vars), p0.createLinear(1, -1, 2));
        assertEquals(parse("-1+2*c", Integers, vars), p0.createLinear(2, -1, 2));
    }

    @Test
    public void testEliminate1() throws Exception {
        assertEquals(parse0("2^14*b"), parse0("a^14*b").eliminate(0, 2));
        assertEquals(parse0("2*a^14"), parse0("a^14*b").eliminate(1, 2));
        String str = "2^14*b - 7*2^9*b^4 + 19*2^9*b^4";
        assertEquals(parse0(str), parse0(str.replace("2^", "a^"), "a", "b").eliminate(0, 2));

        str = "-5*a^22*c*d^13 + 5*a^32*b^24*c*d + a^31*c*d^42 + c^66";
        MultivariatePolynomialZp64 poly = parse0(str);
        assertEquals(parse0(str.replace("d", "3")), poly.eliminate(3, 3));
        assertEquals(parse0(str.replace("c", "3"), "a", "b", "d"), poly.eliminate(2, 3));
        assertEquals(parse0(str.replace("b", "3"), "a", "c", "d"), poly.eliminate(1, 3));
        assertEquals(parse0(str.replace("a", "3"), "b", "c", "d"), poly.eliminate(0, 3));
    }

    @Test
    public void testEvaluate1() throws Exception {
        String[] vars = {"a", "b"};
        assertEquals(parse0("2^14*b", vars), parse0("a^14*b", vars).evaluate(0, 2));
        assertEquals(parse0("2*a^14", vars), parse0("a^14*b", vars).evaluate(1, 2));
        String str = "2^14*b - 7*2^9*b^4 + 19*2^9*b^4";
        assertEquals(parse0(str, vars), parse0(str.replace("2^", "a^"), vars).evaluate(0, 2));

        str = "-5*a^22*c*d^13 + 5*a^32*b^24*c*d + a^31*c*d^42 + c^66";
        vars = new String[]{"a", "b", "c", "d"};
        MultivariatePolynomialZp64 poly = parse0(str, vars);
        assertEquals(parse0(str.replace("d", "3"), vars), poly.evaluate(3, 3));
        assertEquals(parse0(str.replace("c", "3"), vars), poly.evaluate(2, 3));
        assertEquals(parse0(str.replace("b", "3"), vars), poly.evaluate(1, 3));
        assertEquals(parse0(str.replace("a", "3"), vars), poly.evaluate(0, 3));
    }

    @Test
    public void testEvaluate2() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64 poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", Integers, LEX, vars);
        assertEquals(parse0("18 + 18*a^2 + 18*a^3", vars), poly.evaluate(1, 1));

        IntegersZp64 pDomain = new IntegersZp64(17);
        assertEquals(parse("1 + a^2 + a^3", pDomain, LEX, vars), poly.setRing(pDomain).evaluate(1, 1));
    }

    @Test
    public void testEvaluate3() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp64 domain = new IntegersZp64(5642359);
        MultivariatePolynomialZp64 poly = parse(" b^2*c^2+b^3+a*c^4+a*b*c^2+a*b^2+a*b^2*c+2*a^2*c^2+a^2*c^3+a^2*b+a^2*b^2+a^3+a^3*c+a^3*c^2+a^4", domain, LEX, vars);

        int[] evalVars = {1, 2};
        int[] raiseFactors = {2, 1};
        MultivariatePolynomialZp64 r = poly.evaluate(new lPrecomputedPowersHolder(poly.nVariables, evalVars, new long[]{4229599, 9}, domain), evalVars, raiseFactors);
        assertEquals(parse("1694989 + 336131*a + 4996260*a^2 + 91*a^3 + a^4", domain, LEX, vars), r);
    }

    @Test
    public void testUnivar1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp64 domain = new IntegersZp64(17);
        MultivariatePolynomialZp64 poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        poly = poly.evaluate(0, 1);
        assertEquals(UnivariatePolynomialZ64.create(9, 0, 11).modulus(domain.modulus), poly.asUnivariate());

        assertEquals(poly.setRing(Integers), asMultivariate(UnivariatePolynomialZ64.create(9, 0, 11).modulus(Integers.modulus), 2, 1, poly.ordering));
        assertEquals(poly, asMultivariate(UnivariatePolynomialZ64.create(9, 0, 11).modulus(domain.modulus), 2, 1, poly.ordering));
    }

    @Test
    public void testConversion() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64 poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", new IntegersZp64(17), LEX, vars);
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
    }

    @Test
    public void testParse1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp64 domain = new IntegersZp64(17);
        MultivariatePolynomialZp64 poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        MultivariatePolynomialZp64 parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse2() throws Exception {
        String[] vars = {"a"};
        IntegersZp64 domain = new IntegersZp64(17);
        MultivariatePolynomialZp64 poly = parse("8+14*a+16*a^2+11*a^3+12*a^4+a^5", domain, LEX, vars);
        MultivariatePolynomialZp64 parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(1000, 5000);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        IntegersZp64 domain;
        for (int i = 0; i < nIterations; i++) {
            domain = rnd.nextBoolean() ? Integers : new IntegersZp64(getModulusRandom(8));
            MultivariatePolynomialZp64 poly =
                    MultivariatePolynomial.asOverZp64(RandomMultivariatePolynomials.randomPolynomial(
                            rndd.nextInt(1, 4),
                            rndd.nextInt(1, 5),
                            rndd.nextInt(1, 10),
                            domain.asGenericRing(), LEX, rnd));
            MultivariatePolynomialZp64 parsed = parse(poly.toString(vars), domain, LEX, Arrays.copyOf(vars, poly.nVariables));
            assertEquals(poly, parsed);
        }
    }

    @Test
    public void testCoefficient1() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64 poly = parse0("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", vars);
        assertEquals(parse0("7+15*a^2+6*a^3", vars), poly.coefficientOf(1, 2));
        assertEquals(parse0("1+11*b+6*b^2", vars), poly.coefficientOf(0, 3));

        for (int v = 0; v < poly.nVariables; v++) {
            final int var = v;
            MultivariatePolynomialZp64 r = IntStream.rangeClosed(0, poly.degree(var))
                    .mapToObj(i -> poly.coefficientOf(var, i).multiply(new MonomialZp64(poly.nVariables, var, i, 1)))
                    .reduce(poly.createZero(), AMultivariatePolynomial::add);
            assertEquals(poly, r);
        }
    }

    @Test
    public void testRenameVariables1() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 poly = parse0("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
        for (int i = 1; i < poly.nVariables; i++)
            for (int j = 0; j < i; j++)
                assertEquals(poly, swapVariables(swapVariables(poly, i, j), i, j));
    }


    @Test
    public void testRenameVariables2() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 poly = parse0("1 + a^4 + b^6 + c^5 + d^3 + e^2", vars);
        int[] variables = ArraysUtil.sequence(poly.nVariables);
        int[] degrees = poly.degrees();
        int lastVar = 2;

        ArraysUtil.quickSort(ArraysUtil.negate(degrees), variables);
        ArraysUtil.negate(degrees);

        assertEquals(0, variables[lastVar]);
        MultivariatePolynomialZp64 renamed = renameVariables(poly, variables);
        assertDescending(renamed.degrees());
        assertEquals(poly, renameVariables(renamed, inverse(variables)));
    }

    @Test
    public void testSubstitute1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp64 domain = new IntegersZp64(17);
        MultivariatePolynomialZp64 poly = parse("1 + a^2*b^2 + a^3*b^3 + a*b^3 + b^3 + a^2 + 2", domain, vars);
        assertEquals(
                parse("7 + 4*a + a^2 + 4*b^2 + 4*a*b^2 + a^2*b^2 + 11*b^3 + 13*a*b^3 + 6*a^2*b^3 + a^3*b^3", domain, vars),
                poly.shift(0, 2));
        assertEquals(
                parse("2 + 16*a + 2*a^2 + 16*a^3 + 3*b + 3*a*b + 15*a^2*b + 3*a^3*b + 14*b^2 + 14*a*b^2 + a^2*b^2 + 14*a^3*b^2 + b^3 + a*b^3 + a^3*b^3", domain, vars),
                poly.shift(1, -1));
        assertEquals(
                parse("1 + 9*a + 12*a^2 + 16*a^3 + 7*b + 3*a*b + 2*a^2*b + 3*a^3*b + 14*b^2 + 15*a*b^2 + 14*a^2*b^2 + 14*a^3*b^2 + 11*b^3 + 6*a*b^3 + 7*a^2*b^3 + a^3*b^3", domain, vars),
                poly.shift(new int[]{0, 1}, new long[]{-9, -1}));

        assertEquals(
                parse("3 + a^2 + 2*a*b + b^2 + a^2*b^2 + b^3 + 3*a*b^3 + a^3*b^3 + 2*b^4 + 3*a^2*b^4 + 3*a*b^5 + b^6", domain, vars),
                poly.substitute(0, parse("a + b", domain, vars)));
        assertEquals(
                parse("3 + a^2 + 11*a^3*b + a^2*b^2 + 9*a^4*b^2 + b^3 + 16*a*b^3 + 12*a^3*b^3 + 3*a^2*b^4 + 15*a*b^5 + 10*a^5*b^5 + 3*a^2*b^6 + 7*a^6*b^6 + a^3*b^7 + b^8 + 7*a^4*b^8 + 3*a*b^9 + 8*a^2*b^10 + 16*b^12", domain, vars),
                poly.substitute(0, parse("14*a^2*b - b^3 + a", domain, vars)));
    }


    @Test
    public void testSubstitute2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1321349);
        MultivariatePolynomialZp64
                a = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain);

        int[] vars = {1, 2, 3, 4};
        long[] shifts = {762555, 207901, 752954, 112652};
        long[] bShifts = ArraysUtil.negate(shifts.clone());

        assertEquals(a, a.shift(vars, shifts).shift(vars, bShifts));
    }

    @Test
    public void testSubstitute3() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1321349);
        MultivariatePolynomialZp64
                poly = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain),
                expected = parse("8*a^2*d^3*e^4 + 8*a^2*c*d^3*e^4 + 2*a^2*c^2*d^3*e^4 + 224*a*b^3*d^5*e^4 + 560*a*b^3*c*d^5*e^4 + 560*a*b^3*c^2*d^5*e^4 + 280*a*b^3*c^3*d^5*e^4 + 70*a*b^3*c^4*d^5*e^4 + 7*a*b^3*c^5*d^5*e^4", domain);

        assertEquals(expected, poly.shift(2, 2));
    }

    @Test
    public void testSubstitute4() throws Exception {
        IntegersZp64 domain = new IntegersZp64(67);
        MultivariatePolynomialZp64
                a = parse("b*a + b^2 + 1", domain),
                b = parse("a*b + 66*b + 1", domain),
                c = parse("b*a^2 + a^2 + b", domain),
                d = parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);


        assertEquals(parse("64 + 45*a + 18*a^2 + 34*a^3 + 17*a^4 + 5*a^5 + 59*a^6 + 37*a^7 + 55*a^8 + 61*a^9 + 59*a*b + 33*a^2*b + 27*a^3*b + 26*a^4*b + 49*a^5*b + 60*a^6*b + 53*a^7*b + 25*a^8*b + 14*a^9*b + 33*b^2 + 55*a^2*b^2 + 16*a^3*b^2 + 58*a^4*b^2 + 16*a^5*b^2 + 26*a^6*b^2 + 60*a^7*b^2 + 23*a^8*b^2 + 59*a^9*b^2 + 53*b^3 + 23*a*b^3 + 60*a^2*b^3 + 60*a^3*b^3 + 51*a^4*b^3 + 29*a^5*b^3 + 21*a^6*b^3 + 33*a^8*b^3 + 8*a^9*b^3 + 13*b^4 + 55*a*b^4 + 14*a^2*b^4 + 64*a^3*b^4 + 11*a^4*b^4 + 16*a^5*b^4 + 56*a^6*b^4 + 45*a^7*b^4 + 11*a^8*b^4 + 60*a^9*b^4 + 56*b^5 + 11*a*b^5 + 56*a^2*b^5 + 12*a^3*b^5 + a^4*b^5 + 35*a^5*b^5 + 61*a^6*b^5 + 41*a^7*b^5 + 47*a^8*b^5 + 63*a^9*b^5 + 66*b^6 + a*b^6 + 66*a^2*b^6 + a^3*b^6 + 30*a^5*b^6 + 52*a^6*b^6 + 43*a^7*b^6 + 18*a^8*b^6 + 59*a^9*b^6 + 5*a^5*b^7 + 63*a^6*b^7 + 6*a^7*b^7 + 11*a^8*b^7 + 17*a^9*b^7 + 50*a^5*b^8 + 17*a^6*b^8 + 50*a^7*b^8 + 18*a^8*b^8 + a^9*b^8 + 66*a^5*b^9 + a^6*b^9 + 66*a^7*b^9 + a^8*b^9", domain),
                base.shift(1, 2));
    }

    //    @Test
//    public void testMonomialContent1() throws Exception {
//        String[] vars = {"a", "b", "c", "d", "e"};
//        MultivariatePolynomialZp64 poly = parse("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
//        poly = poly.multiply(parse("a*b^2*c^3*d^4*e^5", vars));
//        DegreeVector mc = poly.monomialContent();
//        assertEquals(new DegreeVector(new int[]{1, 2, 3, 4, 5}), mc);
//    }

    @Test
    public void testSeries1() throws Exception {
        int nIterations = its(1000, 1000);
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int n = 0; n < nIterations; n++) {
            long modulus = getModulusRandom(rndd.nextInt(2, 4));
            if (modulus < 10) {
                if (rnd.nextBoolean())
                    modulus = modulus * modulus;
                else if (rnd.nextBoolean() && rnd.nextBoolean())
                    modulus = modulus * modulus * modulus;
                else if (rnd.nextBoolean() && rnd.nextBoolean() && rnd.nextBoolean() && rnd.nextBoolean())
                    modulus = modulus * modulus * modulus * modulus;
            }

            int order = rndd.nextInt((int) modulus + 1, (int) (3 * modulus) + 1);
            int exponent = rndd.nextInt(order, 5 * order);

            assertEquals("" + exponent + "   " + order + "   " + modulus, combinatorialFactorExpected(exponent, order, modulus),
                    MultivariatePolynomialZp64.seriesCoefficientFactor(exponent, order, new IntegersZp64(modulus)));
        }
    }

    @Test
    public void testSeries2() throws Exception {
        long modulus = 9;
        int order = 18;
        int exponent = 21;
        assertEquals("" + exponent + "   " + order + "   " + modulus, combinatorialFactorExpected(exponent, order, modulus),
                MultivariatePolynomialZp64.seriesCoefficientFactor(exponent, order, new IntegersZp64(modulus)));
    }

    private static long combinatorialFactorExpected(int exponent, int order, long modulus) {
        BigInteger result = BigInteger.ONE;
        for (int i = 1; i <= order; ++i)
            result = result.multiply(BigInteger.valueOf(exponent - i + 1));
        for (int i = 1; i <= order; ++i)
            result = result.divide(BigInteger.valueOf(i));
        return result.mod(BigInteger.valueOf(modulus)).longValueExact();
    }

    @Test
    public void testOverMultivariate1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(Integer.MAX_VALUE);
        MultivariatePolynomialZp64 poly = parse("2*a*b*c*d*e + 3*b*c + a^2*b*c", domain);
        assertEquals(poly, MultivariatePolynomialZp64.asNormalMultivariate(poly.asOverMultivariate(0, 3, 4)));
        assertEquals(poly, MultivariatePolynomialZp64.asNormalMultivariate(poly.asOverMultivariateEliminate(0, 3, 4),
                new int[]{0, 3, 4}, new int[]{1, 2}));
    }

    @Test
    public void testOverMultivariate2() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        for (int i = 0; i < its(100, 100); i++) {
            IntegersZp64 domain = new IntegersZp64(getModulusRandom(rndd.nextInt(2, 31)));
            MultivariatePolynomialZp64 poly = RandomMultivariatePolynomials.randomPolynomial(10, 10, 50, domain, MonomialOrder.LEX, rnd);
            for (int j = 0; j < its(10, 10); j++) {
                int[] select = rndd.nextPermutation(poly.nVariables, rndd.nextInt(1, poly.nVariables));
                assertEquals(poly, MultivariatePolynomialZp64.asNormalMultivariate(poly.asOverMultivariate(select)));
                Arrays.sort(select);
                assertEquals(poly, MultivariatePolynomialZp64.asNormalMultivariate(poly.asOverMultivariateEliminate(select),
                        select, ArraysUtil.intSetDifference(ArraysUtil.sequence(poly.nVariables), select)));
            }
        }
    }

    @Test
    public void testSetLC1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(Integer.MAX_VALUE);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                poly = parse("a^2*b^2*c + a*c + 1", domain, vars),
                lc = parse("c^2*b + 1", domain, vars);
        assertEquals(poly.setLC(0, lc).lc(0), lc);
    }

    private static void assertDescending(int[] arr) {
        for (int i = 1; i < arr.length; i++)
            assertTrue(arr[i - 1] >= arr[i]);
    }

    public static int[] inverse(int[] permutation) {
        final int[] inv = new int[permutation.length];
        for (int i = permutation.length - 1; i >= 0; --i)
            inv[permutation[i]] = i;
        return inv;
    }


    @Test
    public void testKroneckerMultiplication1() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> pRing = Rings.MultivariateRing(MultivariatePolynomialZp64.zero(3, Rings.Zp64(17), LEX));
        RandomGenerator rnd = getRandom();
        for (int i = 0; i < 10000; i++) {
            MultivariatePolynomialZp64
                    a = pRing.randomElement(5, 5 + rnd.nextInt(300), rnd),
                    b = pRing.randomElement(5, 5 + rnd.nextInt(300), rnd);

            MultivariatePolynomialZp64.KRONECKER_THRESHOLD = Integer.MAX_VALUE;
            MultivariatePolynomialZp64 plain = a.clone().multiply(b);

            MultivariatePolynomialZp64.KRONECKER_THRESHOLD = 0;
            MultivariatePolynomialZp64 kronecker = a.clone().multiply(b);
            assertEquals(plain, kronecker);
        }
    }

    @Benchmark
    @Test
    public void testKroneckerMultiplication2() throws Exception {
        MultivariateRing<MultivariatePolynomialZp64> pRing =
                Rings.MultivariateRing(MultivariatePolynomialZp64.zero(4, Rings.Zp64(17), LEX));

        for (int i = 0; i < 1000; i++) {
            MultivariatePolynomialZp64
                    a = pRing.randomElement(5, 550),
                    b = pRing.randomElement(5, 550);
            System.out.println(a.sparsity());
            System.out.println(b.sparsity());
            int nIts = 1;
            long size = 0;

            MultivariatePolynomialZp64.KRONECKER_THRESHOLD = Integer.MAX_VALUE;
            long start = System.nanoTime();
            MultivariatePolynomialZp64 ev = a.clone().multiply(b);
            for (int j = 0; j < nIts; j++)
                size += a.clone().multiply(b).size();
            double plain = System.nanoTime() - start;


            MultivariatePolynomialZp64.KRONECKER_THRESHOLD = 10;
            start = System.nanoTime();
            MultivariatePolynomialZp64 ev2 = a.clone().multiply(b);
            for (int j = 0; j < nIts; j++)
                size += a.clone().multiply(b).size();
            double kronecker = System.nanoTime() - start;

            System.out.println(1.0 * a.size() * b.size() / (size / 2 / (nIts + 1)));
            System.out.println(plain / kronecker);
            assertEquals(ev, ev2);
            System.out.println();
        }
    }
}