package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.multivar.MultivariatePolynomial.DegreeVector;
import cc.r2.core.test.AbstractTest;
import org.junit.Test;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.poly.generics.IntegersDomain.IntegersDomain;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.LEX;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariatePolynomialTest extends AbstractTest {
    @Test
    public void testArithmetics1() throws Exception {
        MultivariatePolynomial<BigInteger> a = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", LEX);
        assertEquals(ZERO, a.cc());
        assertEquals(BigInteger.ONE, a.lc());
        assertEquals(BigInteger.ONE, a.clone().increment().cc());
        assertEquals(BigInteger.NEGATIVE_ONE, a.clone().decrement().cc());
        MultivariatePolynomial b = MultivariatePolynomial.parse("a*b - a^2 + c^3*b^2", LEX);
        assertEquals(MultivariatePolynomial.parse("2*a^2", LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(MultivariatePolynomial.parse("2*a*b + 2*c^3*b^2", LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(MultivariatePolynomial.parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", LEX), a.multiply(b));
    }

    @Test
    public void testZero1() throws Exception {
        MultivariatePolynomial<BigInteger> p = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", LEX);

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
    public void testZero2() throws Exception {
        MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.create(
                IntegersDomain, LEX, new DegreeVector[]{
                        new DegreeVector(1, 2, 3),
                        new DegreeVector(0, 1, 2),
                        new DegreeVector(0, 1, 2),
                        new DegreeVector(3, 43, 1)
                }, new BigInteger[]{ZERO, FIVE, FIVE, TEN}
        );
        assertEquals(2, poly.size());
        assertEquals(parse("10*b*c^2 + 10*a^3*b^43*c", LEX), poly);
    }

    private static <E> void assertZero(MultivariatePolynomial<E> a) {
        assertEquals(a.size(), 0);
        assertTrue(a.isZero());
        assertTrue(a.isConstant());
        assertTrue(a.lt().getKey().isZeroVector());
        assertEquals(a.lc(), a.domain.getZero());
        assertEquals(a.cc(), a.domain.getZero());
        assertEquals(a.degree(), 0);
        assertTrue(a.multiDegree().isZeroVector());
    }

    @Test
    public void testArithmetics2() throws Exception {
        ModularDomain algebra = new ModularDomain(17);
        MultivariatePolynomial<BigInteger> a = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", algebra, LEX);
        assertEquals(ZERO, a.cc());
        assertEquals(BigInteger.ONE, a.lc());
        assertEquals(BigInteger.ONE, a.clone().increment().cc());
        assertEquals(algebra.getNegativeOne(), a.clone().decrement().cc());
        MultivariatePolynomial<BigInteger> b = MultivariatePolynomial.parse("a*b - a^2 + c^3*b^2", algebra, LEX);
        assertEquals(MultivariatePolynomial.parse("2*a^2", algebra, LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(MultivariatePolynomial.parse("2*a*b + 2*c^3*b^2", algebra, LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(MultivariatePolynomial.parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", algebra, LEX), a.multiply(b));
    }

    @Test
    public void testZeroVariables() throws Exception {
        MultivariatePolynomial<BigInteger> p0 = MultivariatePolynomial.parse("23", LEX);
        assertEquals(0, p0.nVariables);
        assertEquals(1, p0.size());
        assertEquals(0, p0.clone().subtract(p0).size());
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
}