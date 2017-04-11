package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.multivar.MultivariatePolynomial.DegreeVector;
import org.junit.Test;

import static cc.r2.core.number.BigInteger.*;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.LEX;
import static cc.r2.core.poly.multivar.MultivariatePolynomial.parse;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariatePolynomialTest {
    @Test
    public void testArith1() throws Exception {
        MultivariatePolynomial a = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", LEX);
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
        MultivariatePolynomial p = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", LEX);

        MultivariatePolynomial a = p.clone();
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
        MultivariatePolynomial poly = MultivariatePolynomial.create(
                new BigInteger[]{ZERO, FIVE, FIVE, TEN},
                new DegreeVector[]{
                        new DegreeVector(1, 2, 3),
                        new DegreeVector(0, 1, 2),
                        new DegreeVector(0, 1, 2),
                        new DegreeVector(3, 43, 1)
                }, LEX);
        assertEquals(2, poly.size());
        assertEquals(parse("10*b*c^2 + 10*a^3*b^43*c", LEX), poly);
    }

    private static void assertZero(MultivariatePolynomial a) {
        assertEquals(a.size(), 0);
        assertTrue(a.isZero());
        assertTrue(a.isConstant());
        assertTrue(a.lt().getKey().isZeroVector());
        assertEquals(a.lc(), ZERO);
        assertEquals(a.cc(), ZERO);
        assertEquals(a.degree(), 0);
        assertTrue(a.multiDegree().isZeroVector());
    }
}