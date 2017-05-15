package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersModulo;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Test;

import java.util.Arrays;
import java.util.Comparator;

import static cc.r2.core.poly.Integers.Integers;
import static cc.r2.core.poly.multivar.MonomialTerm.*;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideAndRemainder;
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
        MultivariatePolynomial dividend = parse("a^2*b+a*b^2+b^2", domain, ordering, vars);
        MultivariatePolynomial f1 = parse("a*b - 1", domain, ordering, vars);
        MultivariatePolynomial f2 = parse("b^2 - 1", domain, ordering, vars);

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