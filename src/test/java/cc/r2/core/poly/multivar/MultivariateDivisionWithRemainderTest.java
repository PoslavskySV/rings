package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import cc.r2.core.test.AbstractTest;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;

import static cc.r2.core.poly.multivar.MultivariateDivisionWithRemainder.divideAndRemainder;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar.RandomMultivariatePolynomial.randomPolynomial;
import static org.junit.Assert.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateDivisionWithRemainderTest extends AbstractTest {
    @Test
    public void test1() throws Exception {
        String[] vars = {"a", "b"};
        Comparator<DegreeVector> ordering = LEX;
        MultivariatePolynomial dividend = parse("a*b^2 + 1", ordering, vars);
        MultivariatePolynomial f1 = parse("a*b + 1", ordering, vars);
        MultivariatePolynomial f2 = parse("b + 1", ordering, vars);


        MultivariatePolynomial[] qd;
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
        MultivariatePolynomial dividend = parse("a^2*b+a*b^2+b^2", ordering, vars);
        MultivariatePolynomial f1 = parse("a*b - 1", ordering, vars);
        MultivariatePolynomial f2 = parse("b^2 - 1", ordering, vars);


        MultivariatePolynomial[] qd;
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
        MultivariatePolynomial a = randomPolynomial(5, 10, 10, getRandom());
        MultivariatePolynomial b = randomPolynomial(5, 10, 10, getRandom());
        for (Comparator<DegreeVector> order : Arrays.asList(LEX, REVLEX, GRLEX, GREVLEX)) {
            MultivariatePolynomial c = a.clone().multiply(b).setOrdering(order);
            assertTrue(divideAndRemainder(c, a.setOrdering(order))[1].isZero());
            assertTrue(divideAndRemainder(c, b.setOrdering(order))[1].isZero());
        }
    }

    @Test
    public void test4_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = 100;
        int nVariables = 5;
        int nDividers = 5;
        Comparator<DegreeVector> ordering = LEX;
        for (int n = 0; n < nIterations; n++) {

            MultivariatePolynomial[]
                    dividers = new MultivariatePolynomial[nDividers],
                    quotients = new MultivariatePolynomial[dividers.length];

            MultivariatePolynomial dividend = MultivariatePolynomial.zero(nVariables, ordering);
            for (int j = 0; j < dividers.length; j++) {
                dividers[j] = randomPolynomial(nVariables, rndd.nextInt(1, 10), rndd.nextInt(5, 10), BigInteger.valueOf(100), ordering, rnd);
                quotients[j] = randomPolynomial(nVariables, rndd.nextInt(1, 10), rndd.nextInt(5, 10), BigInteger.valueOf(100), ordering, rnd);
                dividend = dividend.add(dividers[j].clone().multiply(quotients[j]));
            }

            MultivariatePolynomial[] qd = divideAndRemainder(dividend, dividers);
            assertQuotientReminder(qd, dividend, dividers);
        }
    }

    public static void assertQuotientReminder(MultivariatePolynomial dividend, MultivariatePolynomial... dividers) {
        assertQuotientReminder(divideAndRemainder(dividend, dividers), dividend, dividers);
    }

    public static void assertQuotientReminder(MultivariatePolynomial[] quotRem, MultivariatePolynomial dividend, MultivariatePolynomial... dividers) {
        assertEquals(quotRem.length - 1, dividers.length);
        MultivariatePolynomial r = dividend.createZero();
        for (int i = 0; i < dividers.length; i++)
            r.add(dividers[i].clone().multiply(quotRem[i]));
        r.add(quotRem[quotRem.length - 1]);
        assertEquals(dividend, r);
    }
}