package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariatePolynomialTest {
    @Test
    public void testArith1() throws Exception {
        MultivariatePolynomial a = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", MultivariatePolynomial.LEX);
        assertEquals(BigInteger.ZERO, a.cc());
        assertEquals(BigInteger.ONE, a.lc());
        assertEquals(BigInteger.ONE, a.clone().increment().cc());
        assertEquals(BigInteger.NEGATIVE_ONE, a.clone().decrement().cc());

        MultivariatePolynomial b = MultivariatePolynomial.parse("a*b - a^2 + c^3*b^2", MultivariatePolynomial.LEX);
        assertEquals(MultivariatePolynomial.parse("2*a^2", MultivariatePolynomial.LEX, "a", "b", "c"), a.clone().subtract(b));
        assertEquals(MultivariatePolynomial.parse("2*a*b + 2*c^3*b^2", MultivariatePolynomial.LEX, "a", "b", "c"), a.clone().add(b));
        assertEquals(MultivariatePolynomial.parse("-a^4 + a^2*b^2 + 2*a*b^3*c^3 + b^4*c^6", MultivariatePolynomial.LEX), a.multiply(b));
    }

    @Test
    public void testZero1() throws Exception {
        MultivariatePolynomial p = MultivariatePolynomial.parse("a*b + a^2 + c^3*b^2", MultivariatePolynomial.LEX);

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

    private static void assertZero(MultivariatePolynomial a) {
        assertEquals(a.size(), 0);
        assertTrue(a.isZero());
        assertTrue(a.isConstant());
        assertTrue(a.lt().getKey().isZeroVector());
        assertEquals(a.lc(), BigInteger.ZERO);
        assertEquals(a.cc(), BigInteger.ZERO);
        assertEquals(a.degree(), 0);
        assertTrue(a.multiDegree().isZeroVector());
    }
}