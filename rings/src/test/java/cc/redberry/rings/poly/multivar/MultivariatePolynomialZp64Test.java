package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.TimeUnits;
import com.carrotsearch.hppc.LongObjectHashMap;
import com.carrotsearch.hppc.cursors.ObjectCursor;
import com.koloboke.compile.KolobokeMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;
import uk.co.omegaprime.btreemap.LongObjectBTreeMap;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Supplier;
import java.util.stream.IntStream;

import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64.*;
import static cc.redberry.rings.util.ArraysUtil.sum;
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
    public void perf() throws Exception {
        IntegersZp64 domain = new IntegersZp64(Integer.MAX_VALUE);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                a = parse("1 + a + b + c + d + a*b + a^2*d", domain, vars),
                b = parse("1 - a*d + b^2 + c - d + a*c + c^2*d", domain, vars);

        a = PolynomialMethods.polyPow(parse("1 + a + b + c + d + a*b + a^2*d", domain, vars), 5);
        b = PolynomialMethods.polyPow(parse("a*d*c + b^2 + c - d*a^4 + a*c*d + c^2*d", domain, vars), 5);
        a = a.subtract(b);
        b = b.multiply(a);
        a = a.square();

        System.out.println(a.sparsity());
        System.out.println(a.size());
        System.out.println(b.sparsity());
        System.out.println(b.size());

        MultivariatePolynomialZp64 p = a.clone().multiply(b);

        System.out.println(p.sparsity());
        System.out.println(a.size() * b.size());
        System.out.println(p.size());


        LongObjectBTreeMap<PackedMonomial> aPacked = fromPoly(a);
        LongObjectBTreeMap<PackedMonomial> bPacked = fromPoly(b);

        TLongObjectHashMap<PackedMonomial> aPacked2 = fromPoly2(a);
        TLongObjectHashMap<PackedMonomial> bPacked2 = fromPoly2(b);

        KMap aPacked3 = fromPoly3(a);
        KMap bPacked3 = fromPoly3(b);

        LongObjectHashMap<PackedMonomial> aPacked4 = fromPoly4(a);
        LongObjectHashMap<PackedMonomial> bPacked4 = fromPoly4(b);


        System.out.println();
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            Map<DegreeVector, MonomialZp64> tree
                    = multiplyA(() -> new TreeMap<>(p.ordering), domain, a.terms, b.terms);
            System.out.println("TreeMap:" + TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            LongObjectHashMap<PackedMonomial> fast4 = multiplyC(domain, aPacked4, bPacked4);
            System.out.println("Fast4" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            MultivariatePolynomialZp64 btree4 = toPoly4(a.nVariables, domain, fast4);
            System.out.println(btree4.equals(p));

            start = System.nanoTime();
            KMap fast3 = multiplyC(domain, aPacked3, bPacked3);
            System.out.println("Fast3 " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            MultivariatePolynomialZp64 btree3 = toPoly3(a.nVariables, domain, fast3);
            System.out.println(btree3.equals(p));



            start = System.nanoTime();
            TLongObjectHashMap<PackedMonomial> fast2 = multiplyC(domain, aPacked2, bPacked2);
            System.out.println("Fast2 " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            MultivariatePolynomialZp64 btree2 = toPoly2(a.nVariables, domain, fast2);
            System.out.println(btree2.equals(p));



            start = System.nanoTime();
            LongObjectBTreeMap<PackedMonomial> fast = multiplyB(domain, aPacked, bPacked);
            System.out.println("Fast" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            MultivariatePolynomialZp64 btree = toPoly(a.nVariables, domain, fast);
            System.out.println(btree.equals(p));

            start = System.nanoTime();
            Map<DegreeVector, MonomialZp64> hashTable
                    = multiplyA(() -> new Hashtable<>(3 * p.size() / 2, 1.0f), domain, a.terms, b.terms);
            System.out.println("Hashtable:" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println(hashTable.equals(tree));

            start = System.nanoTime();
            Map<DegreeVector, MonomialZp64> hash
                    = multiplyA(() -> new HashMap<>(3 * p.size() / 2, 1.0f), domain, a.terms, b.terms);
            System.out.println("HashMap:" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println(hash.equals(tree));

            start = System.nanoTime();
            Map<DegreeVector, MonomialZp64> chash
                    = multiplyA(() -> new ConcurrentHashMap<>(3 * p.size() / 2, 1.0f), domain, a.terms, b.terms);
            System.out.println("HashMap:" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println(chash.equals(tree));

//            start = System.nanoTime();
//            Map<DegreeVector, MonomialZp64> hash2
//                    = multiplyA(() -> new LinkedHashMap<>(a.size() * b.size(), 1.0f), domain, a.terms, b.terms);
//            System.out.println("LinkedHashMap:" + TimeUnits.nanosecondsToString(System.nanoTime() - start));
//            System.out.println(hash2.equals(tree));

            System.out.println();
            System.out.println();
//            System.out.println(tree.equals(hash));
        }
    }


    static Map<DegreeVector, MonomialZp64> multiplyA(
            Supplier<Map<DegreeVector, MonomialZp64>> factory,
            IntegersZp64 ring,
            Map<DegreeVector, MonomialZp64> a,
            Map<DegreeVector, MonomialZp64> b) {
        Map<DegreeVector, MonomialZp64> newMap = factory.get();
        for (MonomialZp64 othElement : a.values())
            for (MonomialZp64 thisElement : b.values())
                add(newMap, thisElement.multiply(othElement, ring.multiply(thisElement.coefficient, othElement.coefficient)), ring);
        return newMap;
    }

    private static void add(Map<DegreeVector, MonomialZp64> polynomial, MonomialZp64 term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        polynomial.merge(term, term, (o, n) -> {
            long r = ring.add(o.coefficient, n.coefficient);
            if (r == 0)
                return null;
            else {
                return o.setCoefficient(r);
            }
        });
    }

    static LongObjectBTreeMap<PackedMonomial> fromPoly(MultivariatePolynomialZp64 poly) {
        LongObjectBTreeMap<PackedMonomial> map = LongObjectBTreeMap.create();
        for (MonomialZp64 m : poly) {
            long hash = pack(m.exponents, m.totalDegree);
            assert !map.containsKey(hash);
            map.put(hash, new PackedMonomial(hash, m.coefficient));
        }
        return map;
    }

    static MultivariatePolynomialZp64 toPoly(int nVars, IntegersZp64 ring, LongObjectBTreeMap<PackedMonomial> poly) {
        MultivariatePolynomialZp64 p = MultivariatePolynomialZp64.zero(nVars, ring, LEX);
        for (PackedMonomial m : poly.values()) {
            int[] exponents = unpack(m.hash, nVars);
            p.add(new MonomialZp64(exponents, m.coefficient));
        }
        return p;
    }

    static LongObjectBTreeMap<PackedMonomial> multiplyB(
            IntegersZp64 ring,
            LongObjectBTreeMap<PackedMonomial> a,
            LongObjectBTreeMap<PackedMonomial> b) {
        LongObjectBTreeMap<PackedMonomial> newMap = LongObjectBTreeMap.create();
        for (PackedMonomial othElement : a.values())
            for (PackedMonomial thisElement : b.values())
                add(newMap, thisElement.multiply(othElement, ring.multiply(thisElement.coefficient, othElement.coefficient)), ring);
        return newMap;
    }


    static KMap fromPoly3(MultivariatePolynomialZp64 poly) {
        KMap map  = KMap.withExpectedSize(poly.size());
        for (MonomialZp64 m : poly) {
            long hash = pack(m.exponents, m.totalDegree);
            assert !map.containsKey(hash);
            map.put(hash, new PackedMonomial(hash, m.coefficient));
        }
        return map;
    }

    static MultivariatePolynomialZp64 toPoly3(int nVars, IntegersZp64 ring, KMap poly) {
        MultivariatePolynomialZp64 p = MultivariatePolynomialZp64.zero(nVars, ring, LEX);
        for (PackedMonomial m : poly.values()) {
            int[] exponents = unpack(m.hash, nVars);
            p.add(new MonomialZp64(exponents, m.coefficient));
        }
        return p;
    }

    static KMap multiplyC(
            IntegersZp64 ring,
            KMap a,
            KMap b) {
        KMap newMap = KMap.withExpectedSize(a.size() * b.size());
        for (PackedMonomial othElement : a.values())
            for (PackedMonomial thisElement : b.values())
                add(newMap, thisElement.multiply(othElement, ring.multiply(thisElement.coefficient, othElement.coefficient)), ring);
        return newMap;
    }

    static TLongObjectHashMap<PackedMonomial> fromPoly2(MultivariatePolynomialZp64 poly) {
        TLongObjectHashMap<PackedMonomial> map = new TLongObjectHashMap(poly.size() + 1, 1.0f);
        for (MonomialZp64 m : poly) {
            long hash = pack(m.exponents, m.totalDegree);
            assert !map.containsKey(hash);
            map.put(hash, new PackedMonomial(hash, m.coefficient));
        }
        return map;
    }

    static MultivariatePolynomialZp64 toPoly2(int nVars, IntegersZp64 ring, TLongObjectHashMap<PackedMonomial> poly) {
        MultivariatePolynomialZp64 p = MultivariatePolynomialZp64.zero(nVars, ring, LEX);
        for (PackedMonomial m : poly.valueCollection()) {
            int[] exponents = unpack(m.hash, nVars);
            p.add(new MonomialZp64(exponents, m.coefficient));
        }
        return p;
    }


    static LongObjectHashMap<PackedMonomial> fromPoly4(MultivariatePolynomialZp64 poly) {
        LongObjectHashMap<PackedMonomial> map = new LongObjectHashMap<>(poly.size() + 1, 0.98f);
        for (MonomialZp64 m : poly) {
            long hash = pack(m.exponents, m.totalDegree);
            assert !map.containsKey(hash);
            map.put(hash, new PackedMonomial(hash, m.coefficient));
        }
        return map;
    }

    static MultivariatePolynomialZp64 toPoly4(int nVars, IntegersZp64 ring, LongObjectHashMap<PackedMonomial> poly) {
        MultivariatePolynomialZp64 p = MultivariatePolynomialZp64.zero(nVars, ring, LEX);

        for (ObjectCursor<PackedMonomial> mc : poly.values()) {
            PackedMonomial m = mc.value;
            int[] exponents = unpack(m.hash, nVars);
            p.add(new MonomialZp64(exponents, m.coefficient));
        }
        return p;
    }


    static LongObjectHashMap<PackedMonomial> multiplyC(
            IntegersZp64 ring,
            LongObjectHashMap<PackedMonomial> a,
            LongObjectHashMap<PackedMonomial> b) {
        LongObjectHashMap<PackedMonomial> newMap = new LongObjectHashMap<>(a.size() * b.size(), 0.98f);
        for (ObjectCursor<PackedMonomial > othElement : a.values())
            for (ObjectCursor<PackedMonomial > thisElement : b.values())
                add(newMap, thisElement.value.multiply(othElement.value, ring.multiply(thisElement.value.coefficient, othElement.value.coefficient)), ring);
        return newMap;
    }


    static TLongObjectHashMap<PackedMonomial> multiplyC(
            IntegersZp64 ring,
            TLongObjectHashMap<PackedMonomial> a,
            TLongObjectHashMap<PackedMonomial> b) {
        TLongObjectHashMap<PackedMonomial> newMap = new TLongObjectHashMap<>(a.size() * b.size(), 1.0f);
        for (PackedMonomial othElement : a.valueCollection())
            for (PackedMonomial thisElement : b.valueCollection())
                add(newMap, thisElement.multiply(othElement, ring.multiply(thisElement.coefficient, othElement.coefficient)), ring);
        return newMap;
    }

    private static void add(TLongObjectHashMap<PackedMonomial> polynomial, PackedMonomial term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        PackedMonomial old = polynomial.putIfAbsent(term.hash, term);
        if (old != null) {
            long r = ring.add(old.coefficient, term.coefficient);
            if (r == 0)
                polynomial.remove(term.hash);
            else
                polynomial.put(term.hash, new PackedMonomial(term.hash, r));

        }
    }
    private static void add(LongObjectHashMap<PackedMonomial> polynomial, PackedMonomial term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        PackedMonomial old = polynomial.put(term.hash, term);
        if (old != null) {
            long r = ring.add(old.coefficient, term.coefficient);
            if (r == 0)
                polynomial.remove(term.hash);
            else
                polynomial.put(term.hash, new PackedMonomial(term.hash, r));

        }
    }

    private static void add(LongObjectBTreeMap<PackedMonomial> polynomial, PackedMonomial term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        polynomial.merge(term.hash, term, (o, n) -> {
            long r = ring.add(o.coefficient, n.coefficient);
            if (r == 0)
                return null;
            else {
                return new PackedMonomial(o.hash, r);
            }
        });
    }

    private static void add(KMap polynomial, PackedMonomial term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        polynomial.merge(term.hash, term, (o, n) -> {
            long r = ring.add(o.coefficient, n.coefficient);
            if (r == 0)
                return null;
            else {
                return new PackedMonomial(o.hash, r);
            }
        });
    }

    static final Comparator<Long> packedComparator = Long::compareUnsigned;

    @Test
    public void ssss() throws Exception {
        System.out.println(Long.toBinaryString(-1 >>> (64 - 3)));
        int[] arr = {12, 22, 13, 14, 25};
        long hash = pack(arr, sum(arr));
        System.out.println(Long.toBinaryString(hash));
        System.out.println(Arrays.toString(unpack(hash, arr.length)));
        Assert.assertArrayEquals(arr, unpack(hash, arr.length));
    }

    static int[] bitTable, bitDegreeTable;

    static {
        bitTable = new int[12];
        for (int i = 0; i < bitTable.length; i++)
            bitTable[i] = (int) Math.floor(64. / (i + 1));

        bitDegreeTable = new int[12];
        Arrays.fill(bitDegreeTable, 0, 3, Integer.MAX_VALUE);
        for (int i = 3; i < bitDegreeTable.length; i++)
            bitDegreeTable[i] = 1 << bitTable[i];
    }


    static long pack(int[] degrees, int totalDegree) {
        if (bitDegreeTable.length < degrees.length
                || totalDegree >= bitDegreeTable[degrees.length])
            return -1L;

        int nVariables = degrees.length;
        int bitsPerDeg = bitTable[nVariables];

        int i = 0;
        long hash = 0L;
        for (; i < nVariables; ++i)
            hash = hash | ((long) degrees[i]) << ((nVariables - i - 1) * bitsPerDeg);
        hash = hash | ((long) totalDegree) << ((nVariables) * bitsPerDeg);

        return hash;
    }

    static int[] unpack(long hash, int nVariables) {
        assert hash != -1;
        int[] exponents = new int[nVariables];
        int bitsPerDeg = bitTable[nVariables];

        int i = 0;
        for (; i < nVariables; ++i)
            exponents[i] = (int) ((-1 >>> (64 - bitsPerDeg)) & (hash >> ((nVariables - i - 1) * bitsPerDeg)));

        return exponents;
    }

    static final class PackedMonomial {
        final long hash;
        final long coefficient;

        public PackedMonomial(long hash, long coefficient) {
            this.hash = hash;
            this.coefficient = coefficient;
        }

        PackedMonomial multiply(PackedMonomial oth, long newCoeff) {
            return new PackedMonomial(hash + oth.hash, newCoeff);
        }
    }
}