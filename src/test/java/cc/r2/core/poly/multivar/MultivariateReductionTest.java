package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.lIntegersModulo;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;

import static cc.r2.core.poly.Integers.Integers;
import static cc.r2.core.poly.multivar.MonomialTerm.*;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;
import static org.junit.Assert.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateReductionTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        String[] vars = {"a", "b"};
        Comparator<DegreeVector> ordering = LEX;
        MultivariatePolynomial<BigInteger> dividend = parse("a*b^2 + 1", ordering, vars);
        MultivariatePolynomial<BigInteger> f1 = parse("a*b + 1", ordering, vars);
        MultivariatePolynomial<BigInteger> f2 = parse("b + 1", ordering, vars);


        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(dividend, f1, f2);
        assertQuotientReminder(qd, dividend, f1, f2);
        assertArrayEquals(new MultivariatePolynomial[]{
                parse("b", ordering, vars),
                parse("-1", ordering, vars),
                parse("2", ordering, vars)
        }, qd);

        qd = divideAndRemainder(dividend, f2, f1);
        assertQuotientReminder(qd, dividend, f2, f1);
        assertArrayEquals(new MultivariatePolynomial[]{
                parse("a*b - a", ordering, vars),
                parse("0", ordering, vars),
                parse("a+1", ordering, vars)
        }, qd);
    }

    @Test
    public void test2() throws Exception {
        String[] vars = {"a", "b"};
        Comparator<DegreeVector> ordering = LEX;
        MultivariatePolynomial<BigInteger> dividend = parse("a^2*b+a*b^2+b^2", ordering, vars);
        MultivariatePolynomial<BigInteger> f1 = parse("a*b - 1", ordering, vars);
        MultivariatePolynomial<BigInteger> f2 = parse("b^2 - 1", ordering, vars);


        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(dividend, f1, f2);
        assertQuotientReminder(qd, dividend, f1, f2);
        assertArrayEquals(new MultivariatePolynomial[]{
                parse("a+b", ordering, vars),
                parse("1", ordering, vars),
                parse("a+b+1", ordering, vars)
        }, qd);

        assertQuotientReminder(dividend, f2, f1);
    }

    @Test
    public void test3() throws Exception {
        MultivariatePolynomial<BigInteger> a = randomPolynomial(5, 10, 10, getRandom());
        MultivariatePolynomial<BigInteger> b = randomPolynomial(5, 10, 10, getRandom());
        for (Comparator<DegreeVector> order : Arrays.asList(LEX, GRLEX, GREVLEX)) {
            MultivariatePolynomial c = a.clone().multiply(b).setOrdering(order);
            assertTrue(divideAndRemainder(c, a.setOrdering(order))[1].isZero());
            assertTrue(divideAndRemainder(c, b.setOrdering(order))[1].isZero());
        }
    }

    @Test
    public void test4_random() throws Exception {
        testRandomReduce(its(1000, 10000), 5, 3, 1, 5, LEX);
        testRandomReduce(its(1000, 10000), 5, 1, 1, 15, LEX);
    }

    @Test
    public void test5_random() throws Exception {
        for (Comparator<DegreeVector> ord : Arrays.asList(LEX, GRLEX, GREVLEX)) {
            testRandomReduce(its(300, 3000), 5, 3, 1, 5, ord);
            testRandomReduce(its(300, 3000), 5, 1, 1, 15, ord);
        }
    }

    @Test
    public void test6() throws Exception {
        String[] vars = {"a", "b"};
        Domain<BigInteger> domain = new IntegersModulo(2);
        Comparator<DegreeVector> ordering = LEX;
        MultivariatePolynomial<BigInteger> dividend = parse("a^2*b+a*b^2+b^2", domain, ordering, vars);
        MultivariatePolynomial<BigInteger> f1 = parse("a*b - 1", domain, ordering, vars);
        MultivariatePolynomial<BigInteger> f2 = parse("b^2 - 1", domain, ordering, vars);

        MultivariatePolynomial<BigInteger>[] qd;
        qd = divideAndRemainder(dividend, f1, f2);
        assertQuotientReminder(qd, dividend, f1, f2);
        assertArrayEquals(new MultivariatePolynomial[]{
                parse("a+b", ordering, vars),
                parse("1", ordering, vars),
                parse("a+b+1", ordering, vars)
        }, qd);

        assertQuotientReminder(dividend, f2, f1);
    }

    @Test
    public void test7() throws Exception {
        String[] vars = {"a", "b"};
        Domain<BigInteger> domain = new IntegersModulo(2);
        Comparator<DegreeVector> ordering = LEX;
        MultivariatePolynomial<BigInteger> dividend = parse("a^2*b+a*b^2+b^2", domain, ordering, vars);
        MultivariatePolynomial<BigInteger> divider = parse("1", domain, ordering, vars);
        Assert.assertArrayEquals(new MultivariatePolynomial[]{dividend, dividend.createZero()},
                divideAndRemainder(dividend, divider));
    }

    @Test
    public void test8() throws Exception {
        MultivariatePolynomial<BigInteger> dividend = parse("b*c^2+6*b*c^9+2*b*c^16+4*b^15*c^9+3*b^15*c^16+b^15*c^23+4*b^29*c^16+3*b^29*c^23+b^29*c^30+5*a*c+2*a*c^8+3*a*c^15+6*a*b^14*c^8+a*b^14*c^15+5*a*b^14*c^22+6*a*b^28*c^15+a*b^28*c^22+5*a*b^28*c^29+4*a^7*b*c^2+3*a^7*b*c^9+6*a^7*b*c^16+a^7*b*c^23+a^7*b^15*c^9+6*a^7*b^15*c^16+a^7*b^15*c^23+4*a^7*b^15*c^30+6*a^7*b^29*c^30+4*a^7*b^29*c^37+6*a^8*c+a^8*c^8+2*a^8*c^15+5*a^8*c^22+5*a^8*b^14*c^8+2*a^8*b^14*c^15+5*a^8*b^14*c^22+6*a^8*b^14*c^29+2*a^8*b^28*c^29+6*a^8*b^28*c^36+4*a^11*b^11*c^8+3*a^11*b^11*c^15+a^11*b^11*c^22+2*a^11*b^25*c^15+5*a^11*b^25*c^22+4*a^11*b^25*c^29+2*a^11*b^39*c^22+5*a^11*b^39*c^29+4*a^11*b^39*c^36+2*a^12*b^10*c^7+5*a^12*b^10*c^14+4*a^12*b^10*c^21+a^12*b^24*c^14+6*a^12*b^24*c^21+2*a^12*b^24*c^28+a^12*b^38*c^21+6*a^12*b^38*c^28+2*a^12*b^38*c^35+4*a^14*b*c^2+3*a^14*b*c^9+4*a^14*b*c^23+a^14*b*c^30+6*a^14*b^8*c^2+2*a^14*b^8*c^9+a^14*b^8*c^16+5*a^14*b^15*c^23+a^14*b^15*c^30+4*a^14*b^15*c^37+2*a^14*b^22*c^9+2*a^14*b^22*c^16+2*a^14*b^22*c^23+4*a^14*b^29*c^44+a^14*b^36*c^16+3*a^14*b^36*c^23+6*a^15*c+a^15*c^8+6*a^15*c^22+5*a^15*c^29+2*a^15*b^7*c+3*a^15*b^7*c^8+5*a^15*b^7*c^15+4*a^15*b^14*c^22+5*a^15*b^14*c^29+6*a^15*b^14*c^36+3*a^15*b^21*c^8+3*a^15*b^21*c^15+3*a^15*b^21*c^22+6*a^15*b^28*c^43+5*a^15*b^35*c^15+a^15*b^35*c^22+2*a^18*b^11*c^8+5*a^18*b^11*c^15+3*a^18*b^11*c^22+4*a^18*b^11*c^29+4*a^18*b^25*c^15+3*a^18*b^25*c^22+4*a^18*b^25*c^29+2*a^18*b^25*c^36+3*a^18*b^39*c^36+2*a^18*b^39*c^43+a^19*b^10*c^7+6*a^19*b^10*c^14+5*a^19*b^10*c^21+2*a^19*b^10*c^28+2*a^19*b^24*c^14+5*a^19*b^24*c^21+2*a^19*b^24*c^28+a^19*b^24*c^35+5*a^19*b^38*c^35+a^19*b^38*c^42+6*a^21*b*c^16+4*a^21*b*c^23+4*a^21*b*c^30+2*a^21*b^8*c^2+2*a^21*b^8*c^9+6*a^21*b^8*c^16+4*a^21*b^8*c^23+a^21*b^15*c^37+2*a^21*b^22*c^9+6*a^21*b^22*c^16+4*a^21*b^22*c^23+a^21*b^22*c^30+6*a^21*b^36*c^30+2*a^22*c^15+6*a^22*c^22+6*a^22*c^29+3*a^22*b^7*c+3*a^22*b^7*c^8+2*a^22*b^7*c^15+6*a^22*b^7*c^22+5*a^22*b^14*c^36+3*a^22*b^21*c^8+2*a^22*b^21*c^14+2*a^22*b^21*c^15+5*a^22*b^21*c^21+6*a^22*b^21*c^22+4*a^22*b^21*c^28+5*a^22*b^21*c^29+a^22*b^35*c^21+6*a^22*b^35*c^28+2*a^22*b^35*c^29+2*a^22*b^35*c^35+a^22*b^49*c^28+6*a^22*b^49*c^35+2*a^22*b^49*c^42+2*a^25*b^11*c^8+5*a^25*b^11*c^15+2*a^25*b^11*c^29+4*a^25*b^11*c^36+3*a^25*b^18*c^8+a^25*b^18*c^15+4*a^25*b^18*c^22+6*a^25*b^25*c^29+4*a^25*b^25*c^36+2*a^25*b^25*c^43+a^25*b^32*c^15+a^25*b^32*c^22+a^25*b^32*c^29+2*a^25*b^39*c^50+4*a^25*b^46*c^22+5*a^25*b^46*c^29+a^26*b^10*c^7+6*a^26*b^10*c^14+a^26*b^10*c^28+2*a^26*b^10*c^35+5*a^26*b^17*c^7+4*a^26*b^17*c^14+2*a^26*b^17*c^21+3*a^26*b^24*c^28+2*a^26*b^24*c^35+a^26*b^24*c^42+4*a^26*b^31*c^14+4*a^26*b^31*c^21+4*a^26*b^31*c^28+a^26*b^38*c^49+2*a^26*b^45*c^21+6*a^26*b^45*c^28+4*a^28*b*c^30+a^28*b^8*c^2+3*a^28*b^8*c^9+4*a^28*b^8*c^16+a^28*b^8*c^23+4*a^28*b^8*c^30+6*a^28*b^15*c^2+6*a^28*b^15*c^9+a^28*b^15*c^16+5*a^28*b^22*c^23+a^28*b^22*c^37+6*a^28*b^29*c^9+6*a^28*b^29*c^16+4*a^28*b^43*c^16+6*a^29*c^29+5*a^29*b^7*c+a^29*b^7*c^8+6*a^29*b^7*c^15+5*a^29*b^7*c^22+6*a^29*b^7*c^29+2*a^29*b^14*c+2*a^29*b^14*c^8+5*a^29*b^14*c^15+a^29*b^21*c^14+6*a^29*b^21*c^21+4*a^29*b^21*c^22+5*a^29*b^21*c^28+2*a^29*b^21*c^35+5*a^29*b^21*c^36+2*a^29*b^28*c^8+2*a^29*b^28*c^15+2*a^29*b^35*c^21+5*a^29*b^35*c^28+2*a^29*b^35*c^35+a^29*b^35*c^42+6*a^29*b^42*c^15+5*a^29*b^49*c^42+a^29*b^49*c^49+3*a^32*b^11*c^22+2*a^32*b^11*c^29+2*a^32*b^11*c^36+a^32*b^18*c^8+a^32*b^18*c^15+3*a^32*b^18*c^22+2*a^32*b^18*c^29+4*a^32*b^25*c^43+a^32*b^32*c^15+3*a^32*b^32*c^22+2*a^32*b^32*c^29+4*a^32*b^32*c^36+3*a^32*b^46*c^36+5*a^33*b^10*c^21+a^33*b^10*c^28+a^33*b^10*c^35+4*a^33*b^17*c^7+4*a^33*b^17*c^14+5*a^33*b^17*c^21+a^33*b^17*c^28+2*a^33*b^24*c^42+4*a^33*b^31*c^14+5*a^33*b^31*c^21+a^33*b^31*c^28+2*a^33*b^31*c^35+5*a^33*b^45*c^35+6*a^35*b^8*c^16+a^35*b^8*c^30+6*a^35*b^15*c^2+6*a^35*b^15*c^9+5*a^35*b^15*c^16+4*a^35*b^15*c^23+a^35*b^29*c^9+5*a^35*b^29*c^23+2*a^36*b^7*c^15+5*a^36*b^7*c^29+2*a^36*b^14*c+2*a^36*b^14*c^8+4*a^36*b^14*c^15+6*a^36*b^14*c^22+a^36*b^21*c^14+6*a^36*b^21*c^21+a^36*b^21*c^35+2*a^36*b^21*c^42+5*a^36*b^28*c^8+5*a^36*b^28*c^14+4*a^36*b^28*c^21+4*a^36*b^28*c^22+2*a^36*b^28*c^28+3*a^36*b^35*c^35+2*a^36*b^35*c^42+a^36*b^35*c^49+4*a^36*b^42*c^21+4*a^36*b^42*c^28+4*a^36*b^42*c^35+a^36*b^49*c^56+2*a^36*b^56*c^28+6*a^36*b^56*c^35+2*a^39*b^11*c^36+4*a^39*b^18*c^8+5*a^39*b^18*c^15+2*a^39*b^18*c^22+4*a^39*b^18*c^29+2*a^39*b^18*c^36+3*a^39*b^25*c^8+3*a^39*b^25*c^15+4*a^39*b^25*c^22+6*a^39*b^32*c^29+4*a^39*b^32*c^43+3*a^39*b^39*c^15+3*a^39*b^39*c^22+2*a^39*b^53*c^22+a^40*b^10*c^35+2*a^40*b^17*c^7+6*a^40*b^17*c^14+a^40*b^17*c^21+2*a^40*b^17*c^28+a^40*b^17*c^35+5*a^40*b^24*c^7+5*a^40*b^24*c^14+2*a^40*b^24*c^21+3*a^40*b^31*c^28+2*a^40*b^31*c^42+5*a^40*b^38*c^14+5*a^40*b^38*c^21+a^40*b^52*c^21+4*a^42*b^15*c^2+5*a^42*b^15*c^16+4*a^42*b^15*c^30+5*a^42*b^22*c^2+3*a^42*b^22*c^9+a^42*b^36*c^9+6*a^43*b^14*c+4*a^43*b^14*c^15+6*a^43*b^14*c^29+4*a^43*b^21*c+a^43*b^21*c^8+5*a^43*b^21*c^28+a^43*b^21*c^35+a^43*b^21*c^42+4*a^43*b^28*c^14+4*a^43*b^28*c^21+5*a^43*b^28*c^28+a^43*b^28*c^35+5*a^43*b^35*c^8+2*a^43*b^35*c^49+4*a^43*b^42*c^21+5*a^43*b^42*c^28+a^43*b^42*c^35+2*a^43*b^42*c^42+5*a^43*b^56*c^42+3*a^46*b^18*c^22+4*a^46*b^18*c^36+3*a^46*b^25*c^8+3*a^46*b^25*c^15+6*a^46*b^25*c^22+2*a^46*b^25*c^29+4*a^46*b^39*c^15+6*a^46*b^39*c^29+5*a^47*b^17*c^21+2*a^47*b^17*c^35+5*a^47*b^24*c^7+5*a^47*b^24*c^14+3*a^47*b^24*c^21+a^47*b^24*c^28+2*a^47*b^38*c^14+3*a^47*b^38*c^28+a^49*b^22*c^2+6*a^49*b^22*c^16+5*a^50*b^21*c+2*a^50*b^21*c^15+a^50*b^21*c^42+2*a^50*b^28*c^14+6*a^50*b^28*c^21+a^50*b^28*c^28+2*a^50*b^28*c^35+a^50*b^28*c^42+5*a^50*b^35*c^14+5*a^50*b^35*c^21+2*a^50*b^35*c^28+3*a^50*b^42*c^35+2*a^50*b^42*c^49+5*a^50*b^49*c^21+5*a^50*b^49*c^28+a^50*b^63*c^28");
        MultivariatePolynomial<BigInteger> divider = parse("c", "a", "b", "c");

        IntegersModulo domain = new IntegersModulo(7);
        dividend = dividend.setDomain(domain);
        divider = divider.setDomain(domain);
        assertQuotientReminder(divideAndRemainder(dividend, divider), dividend, divider);
    }

    @Test
    public void test9() throws Exception {
        MultivariatePolynomial<BigInteger> dividend = parse("b*c^2+6*b*c^9+2*b*c^16+4*b^15*c^9+3*b^15*c^16+b^15*c^23+4*b^29*c^16+3*b^29*c^23+b^29*c^30+5*a*c+2*a*c^8+3*a*c^15+6*a*b^14*c^8+a*b^14*c^15+5*a*b^14*c^22+6*a*b^28*c^15+a*b^28*c^22+5*a*b^28*c^29+4*a^7*b*c^2+3*a^7*b*c^9+6*a^7*b*c^16+a^7*b*c^23+a^7*b^15*c^9+6*a^7*b^15*c^16+a^7*b^15*c^23+4*a^7*b^15*c^30+6*a^7*b^29*c^30+4*a^7*b^29*c^37+6*a^8*c+a^8*c^8+2*a^8*c^15+5*a^8*c^22+5*a^8*b^14*c^8+2*a^8*b^14*c^15+5*a^8*b^14*c^22+6*a^8*b^14*c^29+2*a^8*b^28*c^29+6*a^8*b^28*c^36+4*a^11*b^11*c^8+3*a^11*b^11*c^15+a^11*b^11*c^22+2*a^11*b^25*c^15+5*a^11*b^25*c^22+4*a^11*b^25*c^29+2*a^11*b^39*c^22+5*a^11*b^39*c^29+4*a^11*b^39*c^36+2*a^12*b^10*c^7+5*a^12*b^10*c^14+4*a^12*b^10*c^21+a^12*b^24*c^14+6*a^12*b^24*c^21+2*a^12*b^24*c^28+a^12*b^38*c^21+6*a^12*b^38*c^28+2*a^12*b^38*c^35+4*a^14*b*c^2+3*a^14*b*c^9+4*a^14*b*c^23+a^14*b*c^30+6*a^14*b^8*c^2+2*a^14*b^8*c^9+a^14*b^8*c^16+5*a^14*b^15*c^23+a^14*b^15*c^30+4*a^14*b^15*c^37+2*a^14*b^22*c^9+2*a^14*b^22*c^16+2*a^14*b^22*c^23+4*a^14*b^29*c^44+a^14*b^36*c^16+3*a^14*b^36*c^23+6*a^15*c+a^15*c^8+6*a^15*c^22+5*a^15*c^29+2*a^15*b^7*c+3*a^15*b^7*c^8+5*a^15*b^7*c^15+4*a^15*b^14*c^22+5*a^15*b^14*c^29+6*a^15*b^14*c^36+3*a^15*b^21*c^8+3*a^15*b^21*c^15+3*a^15*b^21*c^22+6*a^15*b^28*c^43+5*a^15*b^35*c^15+a^15*b^35*c^22+2*a^18*b^11*c^8+5*a^18*b^11*c^15+3*a^18*b^11*c^22+4*a^18*b^11*c^29+4*a^18*b^25*c^15+3*a^18*b^25*c^22+4*a^18*b^25*c^29+2*a^18*b^25*c^36+3*a^18*b^39*c^36+2*a^18*b^39*c^43+a^19*b^10*c^7+6*a^19*b^10*c^14+5*a^19*b^10*c^21+2*a^19*b^10*c^28+2*a^19*b^24*c^14+5*a^19*b^24*c^21+2*a^19*b^24*c^28+a^19*b^24*c^35+5*a^19*b^38*c^35+a^19*b^38*c^42+6*a^21*b*c^16+4*a^21*b*c^23+4*a^21*b*c^30+2*a^21*b^8*c^2+2*a^21*b^8*c^9+6*a^21*b^8*c^16+4*a^21*b^8*c^23+a^21*b^15*c^37+2*a^21*b^22*c^9+6*a^21*b^22*c^16+4*a^21*b^22*c^23+a^21*b^22*c^30+6*a^21*b^36*c^30+2*a^22*c^15+6*a^22*c^22+6*a^22*c^29+3*a^22*b^7*c+3*a^22*b^7*c^8+2*a^22*b^7*c^15+6*a^22*b^7*c^22+5*a^22*b^14*c^36+3*a^22*b^21*c^8+2*a^22*b^21*c^14+2*a^22*b^21*c^15+5*a^22*b^21*c^21+6*a^22*b^21*c^22+4*a^22*b^21*c^28+5*a^22*b^21*c^29+a^22*b^35*c^21+6*a^22*b^35*c^28+2*a^22*b^35*c^29+2*a^22*b^35*c^35+a^22*b^49*c^28+6*a^22*b^49*c^35+2*a^22*b^49*c^42+2*a^25*b^11*c^8+5*a^25*b^11*c^15+2*a^25*b^11*c^29+4*a^25*b^11*c^36+3*a^25*b^18*c^8+a^25*b^18*c^15+4*a^25*b^18*c^22+6*a^25*b^25*c^29+4*a^25*b^25*c^36+2*a^25*b^25*c^43+a^25*b^32*c^15+a^25*b^32*c^22+a^25*b^32*c^29+2*a^25*b^39*c^50+4*a^25*b^46*c^22+5*a^25*b^46*c^29+a^26*b^10*c^7+6*a^26*b^10*c^14+a^26*b^10*c^28+2*a^26*b^10*c^35+5*a^26*b^17*c^7+4*a^26*b^17*c^14+2*a^26*b^17*c^21+3*a^26*b^24*c^28+2*a^26*b^24*c^35+a^26*b^24*c^42+4*a^26*b^31*c^14+4*a^26*b^31*c^21+4*a^26*b^31*c^28+a^26*b^38*c^49+2*a^26*b^45*c^21+6*a^26*b^45*c^28+4*a^28*b*c^30+a^28*b^8*c^2+3*a^28*b^8*c^9+4*a^28*b^8*c^16+a^28*b^8*c^23+4*a^28*b^8*c^30+6*a^28*b^15*c^2+6*a^28*b^15*c^9+a^28*b^15*c^16+5*a^28*b^22*c^23+a^28*b^22*c^37+6*a^28*b^29*c^9+6*a^28*b^29*c^16+4*a^28*b^43*c^16+6*a^29*c^29+5*a^29*b^7*c+a^29*b^7*c^8+6*a^29*b^7*c^15+5*a^29*b^7*c^22+6*a^29*b^7*c^29+2*a^29*b^14*c+2*a^29*b^14*c^8+5*a^29*b^14*c^15+a^29*b^21*c^14+6*a^29*b^21*c^21+4*a^29*b^21*c^22+5*a^29*b^21*c^28+2*a^29*b^21*c^35+5*a^29*b^21*c^36+2*a^29*b^28*c^8+2*a^29*b^28*c^15+2*a^29*b^35*c^21+5*a^29*b^35*c^28+2*a^29*b^35*c^35+a^29*b^35*c^42+6*a^29*b^42*c^15+5*a^29*b^49*c^42+a^29*b^49*c^49+3*a^32*b^11*c^22+2*a^32*b^11*c^29+2*a^32*b^11*c^36+a^32*b^18*c^8+a^32*b^18*c^15+3*a^32*b^18*c^22+2*a^32*b^18*c^29+4*a^32*b^25*c^43+a^32*b^32*c^15+3*a^32*b^32*c^22+2*a^32*b^32*c^29+4*a^32*b^32*c^36+3*a^32*b^46*c^36+5*a^33*b^10*c^21+a^33*b^10*c^28+a^33*b^10*c^35+4*a^33*b^17*c^7+4*a^33*b^17*c^14+5*a^33*b^17*c^21+a^33*b^17*c^28+2*a^33*b^24*c^42+4*a^33*b^31*c^14+5*a^33*b^31*c^21+a^33*b^31*c^28+2*a^33*b^31*c^35+5*a^33*b^45*c^35+6*a^35*b^8*c^16+a^35*b^8*c^30+6*a^35*b^15*c^2+6*a^35*b^15*c^9+5*a^35*b^15*c^16+4*a^35*b^15*c^23+a^35*b^29*c^9+5*a^35*b^29*c^23+2*a^36*b^7*c^15+5*a^36*b^7*c^29+2*a^36*b^14*c+2*a^36*b^14*c^8+4*a^36*b^14*c^15+6*a^36*b^14*c^22+a^36*b^21*c^14+6*a^36*b^21*c^21+a^36*b^21*c^35+2*a^36*b^21*c^42+5*a^36*b^28*c^8+5*a^36*b^28*c^14+4*a^36*b^28*c^21+4*a^36*b^28*c^22+2*a^36*b^28*c^28+3*a^36*b^35*c^35+2*a^36*b^35*c^42+a^36*b^35*c^49+4*a^36*b^42*c^21+4*a^36*b^42*c^28+4*a^36*b^42*c^35+a^36*b^49*c^56+2*a^36*b^56*c^28+6*a^36*b^56*c^35+2*a^39*b^11*c^36+4*a^39*b^18*c^8+5*a^39*b^18*c^15+2*a^39*b^18*c^22+4*a^39*b^18*c^29+2*a^39*b^18*c^36+3*a^39*b^25*c^8+3*a^39*b^25*c^15+4*a^39*b^25*c^22+6*a^39*b^32*c^29+4*a^39*b^32*c^43+3*a^39*b^39*c^15+3*a^39*b^39*c^22+2*a^39*b^53*c^22+a^40*b^10*c^35+2*a^40*b^17*c^7+6*a^40*b^17*c^14+a^40*b^17*c^21+2*a^40*b^17*c^28+a^40*b^17*c^35+5*a^40*b^24*c^7+5*a^40*b^24*c^14+2*a^40*b^24*c^21+3*a^40*b^31*c^28+2*a^40*b^31*c^42+5*a^40*b^38*c^14+5*a^40*b^38*c^21+a^40*b^52*c^21+4*a^42*b^15*c^2+5*a^42*b^15*c^16+4*a^42*b^15*c^30+5*a^42*b^22*c^2+3*a^42*b^22*c^9+a^42*b^36*c^9+6*a^43*b^14*c+4*a^43*b^14*c^15+6*a^43*b^14*c^29+4*a^43*b^21*c+a^43*b^21*c^8+5*a^43*b^21*c^28+a^43*b^21*c^35+a^43*b^21*c^42+4*a^43*b^28*c^14+4*a^43*b^28*c^21+5*a^43*b^28*c^28+a^43*b^28*c^35+5*a^43*b^35*c^8+2*a^43*b^35*c^49+4*a^43*b^42*c^21+5*a^43*b^42*c^28+a^43*b^42*c^35+2*a^43*b^42*c^42+5*a^43*b^56*c^42+3*a^46*b^18*c^22+4*a^46*b^18*c^36+3*a^46*b^25*c^8+3*a^46*b^25*c^15+6*a^46*b^25*c^22+2*a^46*b^25*c^29+4*a^46*b^39*c^15+6*a^46*b^39*c^29+5*a^47*b^17*c^21+2*a^47*b^17*c^35+5*a^47*b^24*c^7+5*a^47*b^24*c^14+3*a^47*b^24*c^21+a^47*b^24*c^28+2*a^47*b^38*c^14+3*a^47*b^38*c^28+a^49*b^22*c^2+6*a^49*b^22*c^16+5*a^50*b^21*c+2*a^50*b^21*c^15+a^50*b^21*c^42+2*a^50*b^28*c^14+6*a^50*b^28*c^21+a^50*b^28*c^28+2*a^50*b^28*c^35+a^50*b^28*c^42+5*a^50*b^35*c^14+5*a^50*b^35*c^21+2*a^50*b^35*c^28+3*a^50*b^42*c^35+2*a^50*b^42*c^49+5*a^50*b^49*c^21+5*a^50*b^49*c^28+a^50*b^63*c^28+2*a^53*b^25*c^8+6*a^53*b^25*c^22+2*a^53*b^25*c^36+6*a^53*b^32*c^8+5*a^53*b^32*c^15+4*a^53*b^46*c^15+a^54*b^24*c^7+3*a^54*b^24*c^21+a^54*b^24*c^35+3*a^54*b^31*c^7+6*a^54*b^31*c^14+2*a^54*b^45*c^14+4*a^56*b^29*c^2+6*a^57*b^28*c+5*a^57*b^28*c^28+2*a^57*b^28*c^42+5*a^57*b^35*c^14+5*a^57*b^35*c^21+3*a^57*b^35*c^28+a^57*b^35*c^35+2*a^57*b^49*c^21+3*a^57*b^49*c^35+4*a^60*b^32*c^8+3*a^60*b^32*c^22+2*a^61*b^31*c^7+5*a^61*b^31*c^21+a^64*b^35*c^14+3*a^64*b^35*c^28+a^64*b^35*c^42+3*a^64*b^42*c^14+6*a^64*b^42*c^21+2*a^64*b^56*c^21+2*a^67*b^39*c^8+a^68*b^38*c^7+2*a^71*b^42*c^14+5*a^71*b^42*c^28+a^78*b^49*c^14");
        MultivariatePolynomial<BigInteger> divider = parse("c", "a", "b", "c");

        IntegersModulo domain = new IntegersModulo(7);
        dividend = dividend.setDomain(domain);
        divider = divider.setDomain(domain);
        assertQuotientReminder(divideAndRemainder(dividend, divider), dividend, divider);
    }

    @Test
    public void test10() throws Exception {
        MultivariatePolynomial<BigInteger> dividend = parse("738*b^18*c^20+28*b^26*c^14+393*a*b^8*c^13+32*a*b^13*c^16+32*a*b^16*c^7+377*a*b^18*c^17+24*a*b^21*c^10+20*a*b^26*c^11+722*a^2*b^9*c^18+42*a^2*b^11*c^25+16*a^2*b^17*c^12+49*a^2*b^24*c^23+48*a^3*b*c^18+36*a^3*b^6*c^21+30*a^3*b^11*c^22+56*a^3*b^14*c^16+361*a^3*b^18*c^14+42*a^3*b^19*c^19+35*a^3*b^24*c^20+8*a^3*b^26*c^8+24*a^4*b^2*c^23+28*a^4*b^15*c^21+12*a^5*b^11*c^19+14*a^5*b^24*c^17+393*a^6*b^8*c^15+623*a^6*b^8*c^17+32*a^6*b^16*c^9+730*a^6*b^16*c^11+262*a^7*b^10*c^17+393*a^7*b^12*c^22+722*a^7*b^18*c^11+32*a^7*b^20*c^16+42*a^7*b^27*c^25+48*a^8*b*c^20+44*a^8*b*c^22+377*a^8*b^14*c^14+56*a^8*b^14*c^18+752*a^8*b^14*c^20+48*a^8*b^17*c^18+20*a^8*b^22*c^8+36*a^8*b^22*c^21+30*a^8*b^27*c^22+32*a^9*b^3*c^22+48*a^9*b^5*c^27+738*a^9*b^16*c^20+968*a^9*b^17*c^18+24*a^9*b^18*c^23+56*a^9*b^18*c^25+726*a^9*b^25*c^12+14*a^9*b^29*c^8+30*a^10*b^7*c^19+7*a^10*b^13*c^24+16*a^10*b^19*c+35*a^10*b^20*c^17+12*a^10*b^24*c^4+12*a^10*b^27*c^19+10*a^10*b^29*c^5+8*a^11*b^3*c^17+706*a^11*b^7*c^22+6*a^11*b^8*c^20+38*a^11*b^10*c^23+5*a^11*b^13*c^21+4*a^11*b^15*c^16+8*a^11*b^20*c^6+745*a^11*b^23*c^21+4*a^12*b^4*c^22+4*a^12*b^29*c^2+6*a^13*c^27+2*a^13*b^13*c^18+7*a^13*b^13*c^25+48*a^13*b^17*c^20+44*a^13*b^17*c^22+32*a^14*b^19*c^22+48*a^14*b^21*c^27+35*a^15*b^18*c^11+16*a^15*b^19*c^3+365*a^15*b^19*c^5+30*a^15*b^23*c^19+8*a^16*b^3*c^19+708*a^16*b^3*c^21+40*a^16*b^8*c^4+30*a^16*b^13*c^7+25*a^16*b^18*c^8+361*a^16*b^21*c^5+16*a^16*b^23*c^10+38*a^16*b^26*c^23+706*a^17*b^5*c^21+8*a^17*b^7*c^26+20*a^17*b^9*c^9+10*a^17*b^25*c^2+5*a^18*b^9*c^18+6*a^18*b^16*c^27+10*a^18*b^18*c^5+49*a^18*b^23*c^20+363*a^18*b^28*c^6+707*a^19*b^12*c^22+56*a^19*b^13*c^13+42*a^19*b^18*c^16+35*a^19*b^23*c^17+28*a^20*b^14*c^18+2*a^20*b^18*c^10+a^21*b^2*c^26+40*a^21*b^8*c^6+387*a^21*b^8*c^8+14*a^21*b^23*c^14+377*a^22*b^10*c^8+40*a^22*b^12*c^13+25*a^23*b^14*c^5+56*a^24*b^13*c^15+752*a^24*b^13*c^17+382*a^24*b^17*c^9+738*a^25*b^15*c^17+56*a^25*b^17*c^22+5*a^26*b^7*c^13+35*a^26*b^19*c^14+745*a^27*b^22*c^18+7*a^29*b^12*c^22");
        MultivariatePolynomial<BigInteger> divider = parse("b", "a", "b", "c");
        assertFalse(dividesQ(dividend, divider));
    }

    @Test
    public void test11() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(1361);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                dividend = lMultivariatePolynomialZp.parse("a^2*b", domain, vars),
                divider = lMultivariatePolynomialZp.parse("b", domain, vars);
        System.out.println(MultivariateReduction.divideExact(dividend, divider));
    }

    static void testRandomReduce(int nIterations, int nVariables, int nDividers,
                                 int minSize, int maxDegree,
                                 Comparator<DegreeVector> ordering) {
        testRandomReduce(nIterations, nVariables, nDividers, minSize, maxDegree, ordering, Integers, getRandom());
    }

    @SuppressWarnings("unchecked")
    static void testRandomReduce(int nIterations, int nVariables, int nDividers,
                                 int minSize, int maxDegree,
                                 Comparator<DegreeVector> ordering,
                                 Domain<BigInteger> domain, RandomGenerator rnd) {
        RandomDataGenerator rndd = new RandomDataGenerator(rnd);
        int maxGeneratedDividendSize = 0;
        int maxFoundQuotientSize = 0;
        for (int n = 0; n < nIterations; n++) {

            MultivariatePolynomial<BigInteger>[]
                    dividers = new MultivariatePolynomial[nDividers],
                    quotients = new MultivariatePolynomial[dividers.length];

            MultivariatePolynomial<BigInteger> dividend = MultivariatePolynomial.zero(nVariables, domain, ordering);
            for (int j = 0; j < dividers.length; j++) {
                dividers[j] = randomPolynomial(nVariables, rndd.nextInt(1, maxDegree), rndd.nextInt(minSize, maxDegree), BigInteger.valueOf(100), domain, ordering, rnd);
                if (dividers[j].isZero()) {
                    --j;
                    continue;
                }
                quotients[j] = randomPolynomial(nVariables, rndd.nextInt(1, maxDegree), rndd.nextInt(minSize, maxDegree), BigInteger.valueOf(100), domain, ordering, rnd);
                dividend = dividend.add(dividers[j].clone().multiply(quotients[j]));
            }
            if (dividend.size() > maxGeneratedDividendSize)
                maxGeneratedDividendSize = dividend.size();

            MultivariatePolynomial<BigInteger>[] qd = divideAndRemainder(dividend, dividers);
            assertNotNull(qd);
            int qSize = Arrays.stream(qd).mapToInt(MultivariatePolynomial::size).max().getAsInt();
            if (qSize > maxFoundQuotientSize)
                maxFoundQuotientSize = qSize;
            assertQuotientReminder(qd, dividend, dividers);
            if (nDividers == 1)
                assertTrue(qd[1].isZero());
        }

        System.out.println("Maximal dividend size meet: " + maxGeneratedDividendSize);
        System.out.println("Maximal quotient size meet: " + maxFoundQuotientSize);
    }

    private static MultivariatePolynomial[] setOrdering(MultivariatePolynomial[] p, Comparator<DegreeVector> ordering) {
        return Arrays.stream(p).map(x -> x.setOrdering(ordering)).toArray(MultivariatePolynomial[]::new);
    }

    public static <E> void assertQuotientReminder(MultivariatePolynomial<E> dividend, MultivariatePolynomial<E>... dividers) {
        assertQuotientReminder(divideAndRemainder(dividend, dividers), dividend, dividers);
    }

    public static <E> void assertQuotientReminder(MultivariatePolynomial<E>[] quotRem, MultivariatePolynomial<E> dividend, MultivariatePolynomial<E>... dividers) {
        assertEquals(quotRem.length - 1, dividers.length);
        MultivariatePolynomial<E> r = dividend.createZero();
        for (int i = 0; i < dividers.length; i++)
            r.add(dividers[i].clone().multiply(quotRem[i]));
        r.add(quotRem[quotRem.length - 1]);
        assertEquals(dividend, r);
    }
}