package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static cc.redberry.rings.bigint.BigInteger.*;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @since 1.0
 */
public class MultivariatePolynomialTest extends AMultivariateTest {
    @Test
    public void testArithmetic1() throws Exception {
        MultivariatePolynomial<BigInteger> a = parse("a*b + a^2 + c^3*b^2", LEX);
        assertEquals(ZERO, a.cc());
        assertEquals(ONE, a.lc());
        assertEquals(ONE, a.clone().increment().cc());
        assertEquals(BigInteger.NEGATIVE_ONE, a.clone().decrement().cc());
        MultivariatePolynomial<BigInteger> b = parse("a*b - a^2 + c^3*b^2", LEX);
        assertEquals(parse("2*a^2", LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", LEX), a.multiply(b));
    }

    @Test
    public void testZero1() throws Exception {
        MultivariatePolynomial<BigInteger> p = parse("a*b + a^2 + c^3*b^2", LEX);

        MultivariatePolynomial<BigInteger> a = p.clone();
        a.subtract(a);
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
        MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.create(3,
                Rings.Z, LEX,
                new Monomial<>(new int[]{1, 2, 3}, ZERO),
                new Monomial<>(new int[]{0, 1, 2}, FIVE),
                new Monomial<>(new int[]{0, 1, 2}, FIVE),
                new Monomial<>(new int[]{3, 43, 1}, TEN));
        assertEquals(2, poly.size());
        assertEquals(parse("10*b*c^2 + 10*a^3*b^43*c", LEX), poly);
    }

    private static <E> void assertZero(MultivariatePolynomial<E> a) {
        assertEquals(a.size(), 0);
        assertTrue(a.isZero());
        assertTrue(a.isConstant());
        assertTrue(a.lt().isZeroVector());
        assertEquals(a.lc(), a.ring.getZero());
        assertEquals(a.cc(), a.ring.getZero());
        assertEquals(a.degree(), 0);
    }

    @Test
    public void testArithmetics2() throws Exception {
        IntegersZp algebra = new IntegersZp(17);
        MultivariatePolynomial<BigInteger> a = parse("a*b + a^2 + c^3*b^2", algebra, LEX);
        assertEquals(ZERO, a.cc());
        assertEquals(ONE, a.lc());
        assertEquals(ONE, a.clone().increment().cc());
        assertEquals(algebra.getNegativeOne(), a.clone().decrement().cc());
        MultivariatePolynomial<BigInteger> b = parse("a*b - a^2 + c^3*b^2", algebra, LEX);
        assertEquals(parse("2*a^2", algebra, LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", algebra, LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", algebra, LEX), a.multiply(b));
    }

    @Test
    public void testZeroVariables() throws Exception {
        MultivariatePolynomial<BigInteger> p0 = parse("23", LEX);
        assertEquals(0, p0.nVariables);
        assertEquals(1, p0.size());
        assertEquals(0, p0.clone().subtract(p0).size());
    }

    @Test
    public void testCreateLinear() throws Exception {
        MultivariatePolynomial<BigInteger> p0 = MultivariatePolynomial.zero(3, Rings.Z, LEX);
        String[] vars = {"a", "b", "c"};
        assertEquals(parse("-1+2*a", vars), p0.createLinear(0, NEGATIVE_ONE, TWO));
        assertEquals(parse("-1+2*b", vars), p0.createLinear(1, NEGATIVE_ONE, TWO));
        assertEquals(parse("-1+2*c", vars), p0.createLinear(2, NEGATIVE_ONE, TWO));
    }

    @Test
    public void testEliminate1() throws Exception {
        assertEquals(parse("2^14*b"), parse("a^14*b").eliminate(0, 2));
        assertEquals(parse("2*a^14"), parse("a^14*b").eliminate(1, 2));
        String str = "2^14*b - 7*2^9*b^4 + 19*2^9*b^4";
        assertEquals(parse(str), parse(str.replace("2^", "a^"), "a", "b").eliminate(0, 2));

        str = "-5*a^22*c*d^13 + 5*a^32*b^24*c*d + a^31*c*d^42 + c^66";
        MultivariatePolynomial<BigInteger> poly = parse(str);
        assertEquals(parse(str.replace("d", "3")), poly.eliminate(3, 3));
        assertEquals(parse(str.replace("c", "3"), "a", "b", "d"), poly.eliminate(2, 3));
        assertEquals(parse(str.replace("b", "3"), "a", "c", "d"), poly.eliminate(1, 3));
        assertEquals(parse(str.replace("a", "3"), "b", "c", "d"), poly.eliminate(0, 3));
    }

    @Test
    public void testEvaluate1() throws Exception {
        String[] vars = {"a", "b"};
        assertEquals(parse("2^14*b", vars), parse("a^14*b", vars).evaluate(0, 2));
        assertEquals(parse("2*a^14", vars), parse("a^14*b", vars).evaluate(1, 2));
        String str = "2^14*b - 7*2^9*b^4 + 19*2^9*b^4";
        assertEquals(parse(str, vars), parse(str.replace("2^", "a^"), vars).evaluate(0, 2));

        str = "-5*a^22*c*d^13 + 5*a^32*b^24*c*d + a^31*c*d^42 + c^66";
        vars = new String[]{"a", "b", "c", "d"};
        MultivariatePolynomial<BigInteger> poly = parse(str, vars);
        assertEquals(parse(str.replace("d", "3"), vars), poly.evaluate(3, 3));
        assertEquals(parse(str.replace("c", "3"), vars), poly.evaluate(2, 3));
        assertEquals(parse(str.replace("b", "3"), vars), poly.evaluate(1, 3));
        assertEquals(parse(str.replace("a", "3"), vars), poly.evaluate(0, 3));
    }

    @Test
    public void testEvaluate2() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", Rings.Z, LEX, vars);
        assertEquals(parse("18 + 18*a^2 + 18*a^3", vars), poly.evaluate(1, 1));

        IntegersZp pDomain = new IntegersZp(17);
        assertEquals(parse("1 + a^2 + a^3", pDomain, LEX, vars), poly.setRing(pDomain).evaluate(1, 1));
    }

    @Test
    public void testEvaluate3() throws Exception {
        String[] vars = {"a", "b", "c"};
        IntegersZp domain = new IntegersZp(5642359);
        MultivariatePolynomial<BigInteger> poly = parse(" b^2*c^2+b^3+a*c^4+a*b*c^2+a*b^2+a*b^2*c+2*a^2*c^2+a^2*c^3+a^2*b+a^2*b^2+a^3+a^3*c+a^3*c^2+a^4", domain, LEX, vars);

        int[] evalVars = {1, 2};
        int[] raiseFactors = {2, 1};
        MultivariatePolynomial<BigInteger> r = poly.evaluate(new PrecomputedPowersHolder<>(poly.nVariables, evalVars, new BigInteger[]{BigInteger.valueOf(4229599), BigInteger.valueOf(9)}, domain), evalVars, raiseFactors);
        assertEquals(parse("1694989 + 336131*a + 4996260*a^2 + 91*a^3 + a^4", domain, LEX, vars), r);
    }

    @Test
    public void testUnivar1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        poly = poly.evaluate(0, 1);
        assertEquals(UnivariatePolynomial.create(9, 0, 11), poly.asUnivariate());
        assertEquals(UnivariatePolynomial.create(9, 0, 11).setRing(domain), poly.asUnivariate());

        assertEquals(poly.setRing(Rings.Z), asMultivariate(UnivariatePolynomial.create(9, 0, 11), 2, 1, poly.ordering));
        assertEquals(poly, asMultivariate(UnivariatePolynomial.create(9, 0, 11).setRing(domain), 2, 1, poly.ordering));
    }

    @Test
    public void testConversion() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", new IntegersZp(17), LEX, vars);
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
    }

    @Test
    public void testParse1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        MultivariatePolynomial<BigInteger> parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse2() throws Exception {
        String[] vars = {"a"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger> poly = parse("8+14*a+16*a^2+11*a^3+12*a^4+a^5", domain, LEX, vars);
        MultivariatePolynomial<BigInteger> parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(1000, 5000);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        Ring<BigInteger> ring;
        for (int i = 0; i < nIterations; i++) {
            ring = rnd.nextBoolean() ? Rings.Z : new IntegersZp(getModulusRandom(8));
            MultivariatePolynomial<BigInteger> poly =
                    RandomMultivariatePolynomials.randomPolynomial(
                            rndd.nextInt(1, 4),
                            rndd.nextInt(1, 5),
                            rndd.nextInt(1, 10),
                            ring, LEX, rnd);
            MultivariatePolynomial<BigInteger> parsed = parse(poly.toString(vars), ring, LEX, Arrays.copyOf(vars, poly.nVariables));
            assertEquals(poly, parsed);
        }
    }

    @Test
    public void testCoefficient1() throws Exception {
        String[] vars = {"a", "b"};
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", vars);
        assertEquals(parse("7+15*a^2+6*a^3", vars), poly.coefficientOf(1, 2));
        assertEquals(parse("1+11*b+6*b^2", vars), poly.coefficientOf(0, 3));

        for (int v = 0; v < poly.nVariables; v++) {
            final int var = v;
            MultivariatePolynomial<BigInteger> r = IntStream.rangeClosed(0, poly.degree(var))
                    .mapToObj(i -> poly.coefficientOf(var, i).multiply(new Monomial<BigInteger>(poly.nVariables, var, i, BigInteger.ONE)))
                    .reduce(poly.createZero(), (a, b) -> a.add(b));
            assertEquals(poly, r);
        }
    }

    @Test
    public void testRenameVariables1() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomial<BigInteger> poly = parse("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
        for (int i = 1; i < poly.nVariables; i++)
            for (int j = 0; j < i; j++)
                assertEquals(poly, swapVariables(swapVariables(poly, i, j), i, j));
    }


    @Test
    public void testRenameVariables2() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomial<BigInteger> poly = parse("1 + a^4 + b^6 + c^5 + d^3 + e^2", vars);
        int[] variables = ArraysUtil.sequence(poly.nVariables);
        int[] degrees = poly.degrees();
        int lastVar = 2;

        ArraysUtil.quickSort(ArraysUtil.negate(degrees), variables);
        ArraysUtil.negate(degrees);

        assertEquals(0, variables[lastVar]);
        MultivariatePolynomial<BigInteger> renamed = renameVariables(poly, variables);
        assertDescending(renamed.degrees());
        assertEquals(poly, renameVariables(renamed, inverse(variables)));
    }

    @Test
    public void testDerivative1() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomial<BigInteger> poly =
                parse("1 + a^4*b^2 + c*b^6 + a*c^5 + a*b*d^3 + a^3*b^3*e^2", vars);
        assertEquals(
                parse("4*a^3*b^2 + c^5 + b*d^3 + 3*a^2*b^3*e^2", vars),
                poly.derivative(0));
        assertEquals(
                parse("2*a^4*b + 6*c*b^5 + a*d^3 + 3*a^3*b^2*e^2", vars),
                poly.derivative(1));
        assertEquals(
                parse("b^6 + 5*a*c^4", vars),
                poly.derivative(2));
        assertEquals(
                parse("3*a*b*d^2", vars),
                poly.derivative(3));
        assertEquals(
                parse("2*a^3*b^3*e", vars),
                poly.derivative(4));
        assertEquals(
                parse("120*b^3*c + 6*a^3*e^2", vars),
                poly.derivative(1, 3));
        assertEquals(
                parse("0", vars),
                poly.derivative(4).derivative(2));

        for (int var = 0; var < poly.nVariables; var++) {
            MultivariatePolynomial<BigInteger> tmp = poly;
            int degree = poly.degree(var);
            for (int i = 0; i <= degree; ++i)
                tmp = tmp.derivative(var);
            assertTrue(tmp.isZero());
        }
    }

    @Test
    public void testSubstitute1() throws Exception {
        String[] vars = {"a", "b"};
        IntegersZp domain = new IntegersZp(17);
        MultivariatePolynomial<BigInteger> poly = parse("1 + a^2*b^2 + a^3*b^3 + a*b^3 + b^3 + a^2 + 2", domain, vars);
        assertEquals(
                parse("7 + 4*a + a^2 + 4*b^2 + 4*a*b^2 + a^2*b^2 + 11*b^3 + 13*a*b^3 + 6*a^2*b^3 + a^3*b^3", domain, vars),
                poly.shift(0, 2));
        assertEquals(
                parse("2 + 16*a + 2*a^2 + 16*a^3 + 3*b + 3*a*b + 15*a^2*b + 3*a^3*b + 14*b^2 + 14*a*b^2 + a^2*b^2 + 14*a^3*b^2 + b^3 + a*b^3 + a^3*b^3", domain, vars),
                poly.shift(1, -1));
        assertEquals(
                parse("1 + 9*a + 12*a^2 + 16*a^3 + 7*b + 3*a*b + 2*a^2*b + 3*a^3*b + 14*b^2 + 15*a*b^2 + 14*a^2*b^2 + 14*a^3*b^2 + 11*b^3 + 6*a*b^3 + 7*a^2*b^3 + a^3*b^3", domain, vars),
                poly.shift(new int[]{0, 1}, new BigInteger[]{BigInteger.valueOf(-9), BigInteger.valueOf(-1)}));

        assertEquals(
                parse("3 + a^2 + 2*a*b + b^2 + a^2*b^2 + b^3 + 3*a*b^3 + a^3*b^3 + 2*b^4 + 3*a^2*b^4 + 3*a*b^5 + b^6", domain, vars),
                poly.substitute(0, parse("a + b", domain, vars)));
        assertEquals(
                parse("3 + a^2 + 11*a^3*b + a^2*b^2 + 9*a^4*b^2 + b^3 + 16*a*b^3 + 12*a^3*b^3 + 3*a^2*b^4 + 15*a*b^5 + 10*a^5*b^5 + 3*a^2*b^6 + 7*a^6*b^6 + a^3*b^7 + b^8 + 7*a^4*b^8 + 3*a*b^9 + 8*a^2*b^10 + 16*b^12", domain, vars),
                poly.substitute(0, parse("14*a^2*b - b^3 + a", domain, vars)));
    }

    @Test
    public void testSubstitute2() throws Exception {
        IntegersZp domain = new IntegersZp(1321349);
        MultivariatePolynomial<BigInteger>
                a = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain);

        int[] vars = {1, 2, 3, 4};
        BigInteger[] shifts = {BigInteger.valueOf(762555), BigInteger.valueOf(207901), BigInteger.valueOf(752954), BigInteger.valueOf(112652)};
        BigInteger[] bShifts = ArraysUtil.negate(shifts.clone());

        assertEquals(a, a.shift(vars, shifts).shift(vars, bShifts));
    }

    @Test
    public void testSubstitute3() throws Exception {
        IntegersZp domain = new IntegersZp(1321349);
        MultivariatePolynomial<BigInteger>
                poly = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain),
                expected = parse("8*a^2*d^3*e^4 + 8*a^2*c*d^3*e^4 + 2*a^2*c^2*d^3*e^4 + 224*a*b^3*d^5*e^4 + 560*a*b^3*c*d^5*e^4 + 560*a*b^3*c^2*d^5*e^4 + 280*a*b^3*c^3*d^5*e^4 + 70*a*b^3*c^4*d^5*e^4 + 7*a*b^3*c^5*d^5*e^4", domain);

        assertEquals(expected, poly.shift(2, 2));
    }

    @Test
    public void testSubstitute4() throws Exception {
        IntegersZp domain = new IntegersZp(67);
        MultivariatePolynomial<BigInteger>
                a = parse("b*a + b^2 + 1", domain),
                b = parse("a*b + 66*b + 1", domain),
                c = parse("b*a^2 + a^2 + b", domain),
                d = parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);

        assertEquals(parse("64 + 45*a + 18*a^2 + 34*a^3 + 17*a^4 + 5*a^5 + 59*a^6 + 37*a^7 + 55*a^8 + 61*a^9 + 59*a*b + 33*a^2*b + 27*a^3*b + 26*a^4*b + 49*a^5*b + 60*a^6*b + 53*a^7*b + 25*a^8*b + 14*a^9*b + 33*b^2 + 55*a^2*b^2 + 16*a^3*b^2 + 58*a^4*b^2 + 16*a^5*b^2 + 26*a^6*b^2 + 60*a^7*b^2 + 23*a^8*b^2 + 59*a^9*b^2 + 53*b^3 + 23*a*b^3 + 60*a^2*b^3 + 60*a^3*b^3 + 51*a^4*b^3 + 29*a^5*b^3 + 21*a^6*b^3 + 33*a^8*b^3 + 8*a^9*b^3 + 13*b^4 + 55*a*b^4 + 14*a^2*b^4 + 64*a^3*b^4 + 11*a^4*b^4 + 16*a^5*b^4 + 56*a^6*b^4 + 45*a^7*b^4 + 11*a^8*b^4 + 60*a^9*b^4 + 56*b^5 + 11*a*b^5 + 56*a^2*b^5 + 12*a^3*b^5 + a^4*b^5 + 35*a^5*b^5 + 61*a^6*b^5 + 41*a^7*b^5 + 47*a^8*b^5 + 63*a^9*b^5 + 66*b^6 + a*b^6 + 66*a^2*b^6 + a^3*b^6 + 30*a^5*b^6 + 52*a^6*b^6 + 43*a^7*b^6 + 18*a^8*b^6 + 59*a^9*b^6 + 5*a^5*b^7 + 63*a^6*b^7 + 6*a^7*b^7 + 11*a^8*b^7 + 17*a^9*b^7 + 50*a^5*b^8 + 17*a^6*b^8 + 50*a^7*b^8 + 18*a^8*b^8 + a^9*b^8 + 66*a^5*b^9 + a^6*b^9 + 66*a^7*b^9 + a^8*b^9", domain),
                base.shift(1, 2));
    }

    @Test
    public void testOverMultivariate1() throws Exception {
        IntegersZp domain = new IntegersZp(Integer.MAX_VALUE);
        MultivariatePolynomial<BigInteger> poly = parse("2*a*b*c*d*e + 3*b*c + a^2*b*c", domain);
        assertEquals(poly, MultivariatePolynomial.asNormalMultivariate(poly.asOverMultivariate(0, 3, 4)));
        assertEquals(poly, MultivariatePolynomial.asNormalMultivariate(poly.asOverMultivariateEliminate(0, 3, 4),
                new int[]{0, 3, 4}, new int[]{1, 2}));
    }

    @Test
    public void testOverMultivariate2() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        for (int i = 0; i < its(100, 100); i++) {
            IntegersZp domain = new IntegersZp(getModulusRandom(rndd.nextInt(2, 31)));
            MultivariatePolynomial<BigInteger> poly = RandomMultivariatePolynomials.randomPolynomial(10, 10, 50, domain, MonomialOrder.LEX, rnd);
            for (int j = 0; j < its(3, 10); j++) {
                int[] select = rndd.nextPermutation(poly.nVariables, rndd.nextInt(1, poly.nVariables));
                assertEquals(poly, MultivariatePolynomial.asNormalMultivariate(poly.asOverMultivariate(select)));
                Arrays.sort(select);
                assertEquals(poly, MultivariatePolynomial.asNormalMultivariate(poly.asOverMultivariateEliminate(select),
                        select, ArraysUtil.intSetDifference(ArraysUtil.sequence(poly.nVariables), select)));
            }
        }
    }

//    @Test
//    public void testMonomialContent1() throws Exception {
//        String[] vars = {"a", "b", "c", "d", "e"};
//        MultivariatePolynomial<BigInteger> poly = parse("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
//        poly = poly.multiply(parse("a*b^2*c^3*d^4*e^5", vars));
//        DegreeVector mc = poly.monomialContent();
//        assertEquals(new DegreeVector(new int[]{1, 2, 3, 4, 5}), mc);
//    }

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
}