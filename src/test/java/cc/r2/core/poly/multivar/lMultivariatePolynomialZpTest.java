package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.Arrays;
import java.util.stream.IntStream;

import static cc.r2.core.poly.multivar.DegreeVector.LEX;
import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class lMultivariatePolynomialZpTest extends AbstractPolynomialTest {

    private lIntegersModulo Integers = new lIntegersModulo(Integer.MAX_VALUE);

    private lMultivariatePolynomialZp parse0(String str) {
        return lMultivariatePolynomialZp.parse(str, Integers);
    }

    private lMultivariatePolynomialZp parse0(String str, String... vars) {
        return lMultivariatePolynomialZp.parse(str, Integers, vars);
    }

    @Test
    public void testArithmetics1() throws Exception {
        lMultivariatePolynomialZp a = parse("a*b + a^2 + c^3*b^2", Integers);
        assertEquals(0, a.cc());
        assertEquals(1, a.lc());
        assertEquals(1, a.clone().increment().cc());
        assertEquals(Integers.modulus(-1), a.clone().decrement().cc());
        lMultivariatePolynomialZp b = parse("a*b - a^2 + c^3*b^2", Integers);
        assertEquals(parse("2*a^2", Integers, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", Integers, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", Integers), a.multiply(b));
    }

    @Test
    public void testZero1() throws Exception {
        lMultivariatePolynomialZp p = parse("a*b + a^2 + c^3*b^2", Integers);

        lMultivariatePolynomialZp a = p.clone();
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
        lMultivariatePolynomialZp poly = lMultivariatePolynomialZp.create(3,
                Integers, LEX,
                new lMonomialTerm(new int[]{1, 2, 3}, 0),
                new lMonomialTerm(new int[]{0, 1, 2}, 5),
                new lMonomialTerm(new int[]{0, 1, 2}, 5),
                new lMonomialTerm(new int[]{3, 43, 1}, 10));
        assertEquals(2, poly.size());
        assertEquals(parse("10*b*c^2 + 10*a^3*b^43*c", Integers), poly);
    }

    private static void assertZero(lMultivariatePolynomialZp a) {
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
        lIntegersModulo algebra = new lIntegersModulo(17);
        lMultivariatePolynomialZp a = parse("a*b + a^2 + c^3*b^2", algebra, LEX);
        assertEquals(0, a.cc());
        assertEquals(1, a.lc());
        assertEquals(1, a.clone().increment().cc());
        assertEquals(algebra.modulus(-1), a.clone().decrement().cc());
        lMultivariatePolynomialZp b = parse("a*b - a^2 + c^3*b^2", algebra, LEX);
        assertEquals(parse("2*a^2", algebra, LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(parse("2*a*b + 2*c^3*b^2", algebra, LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", algebra, LEX), a.multiply(b));
    }

    @Test
    public void testZeroVariables() throws Exception {
        lMultivariatePolynomialZp p0 = parse("23", Integers);
        assertEquals(0, p0.nVariables);
        assertEquals(1, p0.size());
        assertEquals(0, p0.clone().subtract(p0).size());
    }

    @Test
    public void testCreateLinear() throws Exception {
        lMultivariatePolynomialZp p0 = lMultivariatePolynomialZp.zero(3, Integers, LEX);
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
        lMultivariatePolynomialZp poly = parse0(str);
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
        lMultivariatePolynomialZp poly = parse0(str, vars);
        assertEquals(parse0(str.replace("d", "3"), vars), poly.evaluate(3, 3));
        assertEquals(parse0(str.replace("c", "3"), vars), poly.evaluate(2, 3));
        assertEquals(parse0(str.replace("b", "3"), vars), poly.evaluate(1, 3));
        assertEquals(parse0(str.replace("a", "3"), vars), poly.evaluate(0, 3));
    }

    @Test
    public void testEvaluate2() throws Exception {
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", Integers, LEX, vars);
        assertEquals(parse0("18 + 18*a^2 + 18*a^3", vars), poly.evaluate(1, 1));

        lIntegersModulo pDomain = new lIntegersModulo(17);
        assertEquals(parse("1 + a^2 + a^3", pDomain, LEX, vars), poly.setDomain(pDomain).evaluate(1, 1));
    }

    @Test
    public void testEvaluate3() throws Exception {
        String[] vars = {"a", "b", "c"};
        lIntegersModulo domain = new lIntegersModulo(5642359);
        lMultivariatePolynomialZp poly = parse(" b^2*c^2+b^3+a*c^4+a*b*c^2+a*b^2+a*b^2*c+2*a^2*c^2+a^2*c^3+a^2*b+a^2*b^2+a^3+a^3*c+a^3*c^2+a^4", domain, LEX, vars);

        int[] evalVars = {1, 2};
        int[] raiseFactors = {2, 1};
        lMultivariatePolynomialZp r = poly.evaluate(new lPrecomputedPowersHolder(new long[]{4229599, 9}, domain), evalVars, raiseFactors);
        assertEquals(parse("1694989 + 336131*a + 4996260*a^2 + 91*a^3 + a^4", domain, LEX, vars), r);
    }

    @Test
    public void testUnivar1() throws Exception {
        String[] vars = {"a", "b"};
        lIntegersModulo domain = new lIntegersModulo(17);
        lMultivariatePolynomialZp poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        poly = poly.evaluate(0, 1);
        assertEquals(lUnivariatePolynomialZ.create(9, 0, 11).modulus(domain.modulus), poly.asUnivariate());

        assertEquals(poly.setDomain(Integers), asMultivariate(lUnivariatePolynomialZ.create(9, 0, 11).modulus(Integers.modulus), 2, 1, poly.ordering));
        assertEquals(poly, asMultivariate(lUnivariatePolynomialZ.create(9, 0, 11).modulus(domain.modulus), 2, 1, poly.ordering));
    }

    @Test
    public void testConversion() throws Exception {
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", new lIntegersModulo(17), LEX, vars);
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
        assertEquals(poly, asNormalMultivariate(poly.asOverUnivariateEliminate(1), 1));
    }

    @Test
    public void testParse1() throws Exception {
        String[] vars = {"a", "b"};
        lIntegersModulo domain = new lIntegersModulo(17);
        lMultivariatePolynomialZp poly = parse("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", domain, LEX, vars);
        lMultivariatePolynomialZp parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse2() throws Exception {
        String[] vars = {"a"};
        lIntegersModulo domain = new lIntegersModulo(17);
        lMultivariatePolynomialZp poly = parse("8+14*a+16*a^2+11*a^3+12*a^4+a^5", domain, LEX, vars);
        lMultivariatePolynomialZp parsed = parse(poly.toString(vars), domain, LEX, vars);
        assertEquals(poly, parsed);
    }

    @Test
    public void testParse_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(1000, 5000);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        lIntegersModulo domain;
        for (int i = 0; i < nIterations; i++) {
            domain = rnd.nextBoolean() ? Integers : new lIntegersModulo(getModulusRandom(8));
            lMultivariatePolynomialZp poly =
                    MultivariatePolynomial.asLongPolyZp(RandomMultivariatePolynomial.randomPolynomial(
                            rndd.nextInt(1, 4),
                            rndd.nextInt(1, 5),
                            rndd.nextInt(1, 10),
                            BigInteger.valueOf(1000), domain.asDomain(), LEX, rnd));
            lMultivariatePolynomialZp parsed = parse(poly.toString(), domain, LEX, Arrays.copyOf(vars, poly.nVariables));
            assertEquals(poly, parsed);
        }
    }

    @Test
    public void testCoefficient1() throws Exception {
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp poly = parse0("5+6*b+7*b^2+3*a^2+15*a^2*b^2+a^3+11*a^3*b+6*a^3*b^2", vars);
        assertEquals(parse0("7+15*a^2+6*a^3", vars), poly.coefficientOf(1, 2));
        assertEquals(parse0("1+11*b+6*b^2", vars), poly.coefficientOf(0, 3));

        for (int v = 0; v < poly.nVariables; v++) {
            final int var = v;
            lMultivariatePolynomialZp r = IntStream.rangeClosed(0, poly.degree(var))
                    .mapToObj(i -> poly.coefficientOf(var, i).multiply(new lMonomialTerm(poly.nVariables, var, i, 1)))
                    .reduce(poly.createZero(), AMultivariatePolynomial::add);
            assertEquals(poly, r);
        }
    }

    @Test
    public void testRenameVariables1() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp poly = parse0("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
        for (int i = 1; i < poly.nVariables; i++)
            for (int j = 0; j < i; j++)
                assertEquals(poly, swapVariables(swapVariables(poly, i, j), i, j));
    }


    @Test
    public void testRenameVariables2() throws Exception {
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp poly = parse0("1 + a^4 + b^6 + c^5 + d^3 + e^2", vars);
        int[] variables = ArraysUtil.sequence(poly.nVariables);
        int[] degrees = poly.degrees();
        int lastVar = 2;

        ArraysUtil.quickSort(ArraysUtil.negate(degrees), variables);
        ArraysUtil.negate(degrees);

        assertEquals(0, variables[lastVar]);
        lMultivariatePolynomialZp renamed = renameVariables(poly, variables);
        assertDescending(renamed.degrees());
        assertEquals(poly, renameVariables(renamed, inverse(variables)));
    }

    @Test
    public void testSubstitute1() throws Exception {
        String[] vars = {"a", "b"};
        lIntegersModulo domain = new lIntegersModulo(17);
        lMultivariatePolynomialZp poly = parse("1 + a^2*b^2 + a^3*b^3 + a*b^3 + b^3 + a^2 + 2", domain, vars);
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
        lIntegersModulo domain = new lIntegersModulo(1321349);
        lMultivariatePolynomialZp
                a = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain);

        int[] vars = {1, 2, 3, 4};
        long[] shifts = {762555, 207901, 752954, 112652};
        long[] bShifts = ArraysUtil.negate(shifts.clone());

        assertEquals(a, a.shift(vars, shifts).shift(vars, bShifts));
    }

    @Test
    public void testSubstitute3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1321349);
        lMultivariatePolynomialZp
                poly = parse("7*a*b^3*c^5*d^5*e^4 + 2*a^2*c^2*d^3*e^4", domain),
                expected = parse("8*a^2*d^3*e^4 + 8*a^2*c*d^3*e^4 + 2*a^2*c^2*d^3*e^4 + 224*a*b^3*d^5*e^4 + 560*a*b^3*c*d^5*e^4 + 560*a*b^3*c^2*d^5*e^4 + 280*a*b^3*c^3*d^5*e^4 + 70*a*b^3*c^4*d^5*e^4 + 7*a*b^3*c^5*d^5*e^4", domain);

        assertEquals(expected, poly.shift(2, 2));
    }

    //    @Test
//    public void testMonomialContent1() throws Exception {
//        String[] vars = {"a", "b", "c", "d", "e"};
//        lMultivariatePolynomialZp poly = parse("5+6*b*e+7*b^2+3*a^2+d+15*a^2*b^2+a^3+11*a^3*b*e^9+6*a^3*b^2+c*e^3", vars);
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