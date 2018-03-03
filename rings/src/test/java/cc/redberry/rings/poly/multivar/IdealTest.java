package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Rational;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import org.junit.Test;

import static cc.redberry.rings.Rings.Q;
import static cc.redberry.rings.Rings.Z;
import static cc.redberry.rings.poly.multivar.MonomialOrder.GREVLEX;
import static cc.redberry.rings.poly.multivar.MonomialOrder.LEX;
import static org.junit.Assert.assertEquals;

/**
 * @since 2.3
 */
public class IdealTest {
    @Test
    public void test1() throws Exception {
        String[] vars = {"x", "y", "z"};
        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                I1 = Ideal.parse(new String[]{"x", "y"}, Q, MonomialOrder.LEX, vars),
                I2 = Ideal.parse(new String[]{"y^2", "z"}, Q, MonomialOrder.LEX, vars),
                I3 = Ideal.parse(new String[]{"y^2", "x*z", "y*z"}, Q, MonomialOrder.LEX, vars);
        assertEquals(I3, I1.intersection(I2));
    }

    @Test
    public void test2() throws Exception {
        String[] vars = {"x1", "x2", "x3"};
        String[] gensI = {"x1^2*x2 - 1", "x1 + x2 - x3"};
        String[] gensJ = {"x1^2*x3^2 - x2^2", "x1^3 + x2 - x3^2 - 1"};

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                I = Ideal.parse(gensI, Q, vars),
                J = Ideal.parse(gensJ, Q, vars);

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedIJ =
                Ideal.parse(new String[]{
                        "x1^4+x1^3*x2-x1^3*x3-x1*x3^2-x2*x3^2+x3^3+x1*x2+x2^2-x2*x3-x1-x2+x3",
                        "x1*x3^4+x2*x3^4-x3^5-x1^2*x2^2-x1*x2^3+x1*x2^2*x3-x1*x2*x3^2-x2^2*x3^2+x2*x3^3+x1*x3^2+x2*x3^2-x3^3",
                        "x1^3*x3^2+x1^2*x2*x3^2-x1^2*x3^3-x1*x2^2-x2^3+x2^2*x3",
                        "x1^3*x2^3-2*x1^3*x2^2*x3-x1^2*x2^2*x3^2+x1^2*x2*x3^3-x2^3*x3^2+2*x2^2*x3^3-x2*x3^4+x1*x2^3+2*x2^4-3*x2^3*x3+x2^2*x3^2-x1^3-x2^3+2*x2^2*x3-x2*x3^2+x3^2-x2+1",
                        "x2^3*x3^4-2*x2^2*x3^5+x2*x3^6-x1*x2^5+2*x1*x2^4*x3-x1*x2^3*x3^2-x2^4*x3^2+2*x2^3*x3^3-x2^2*x3^4+x2^3*x3^2-2*x2^2*x3^3+x2*x3^4-x3^4+x1*x2^2+x2*x3^2-x3^2",
                        "x1^2*x2^3*x3^2-2*x1^2*x2^2*x3^3+x2^3*x3^4-2*x2^2*x3^5+x2*x3^6-x1*x2^5+2*x1^3*x2^2*x3+2*x1*x2^4*x3+2*x1^2*x2^2*x3^2-x1*x2^3*x3^2-x2^4*x3^2-x1^2*x2*x3^3+2*x2^3*x3^3-x2^2*x3^4-x2^5+2*x2^4*x3-x1^2*x2*x3^2+x2^3*x3^2-4*x2^2*x3^3+2*x2*x3^4-x1*x2^3-2*x2^4+3*x2^3*x3-x1^2*x3^2-x2^2*x3^2+x1^3+x2^3-2*x2^2*x3+x2*x3^2+x2^2-x3^2+x2-1"
                }, Q, vars);
        assertEquals(expectedIJ, I.multiply(J));


        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedIntersectIJ = Ideal.parse(new String[]{
                "x1^4+x1^3*x2-x1^3*x3-x1*x3^2-x2*x3^2+x3^3+x1*x2+x2^2-x2*x3-x1-x2+x3",
                "x1*x3^4+x2*x3^4-x3^5-x1^2*x2^2-x1*x2^3+x1*x2^2*x3-x1*x2*x3^2-x2^2*x3^2+x2*x3^3+x1*x3^2+x2*x3^2-x3^3",
                "x1^3*x3^2+x1^2*x2*x3^2-x1^2*x3^3-x1*x2^2-x2^3+x2^2*x3",
                "x1^3*x2^3-2*x1^3*x2^2*x3-x1^2*x2^2*x3^2+x1^2*x2*x3^3-x2^3*x3^2+2*x2^2*x3^3-x2*x3^4+x1*x2^3+2*x2^4-3*x2^3*x3+x2^2*x3^2-x1^3-x2^3+2*x2^2*x3-x2*x3^2+x3^2-x2+1",
                "x2^3*x3^4-2*x2^2*x3^5+x2*x3^6-x1*x2^5+2*x1*x2^4*x3-x1*x2^3*x3^2-x2^4*x3^2+2*x2^3*x3^3-x2^2*x3^4+x2^3*x3^2-2*x2^2*x3^3+x2*x3^4-x3^4+x1*x2^2+x2*x3^2-x3^2",
                "8*x1^2*x2^3*x3^2-16*x1^2*x2^2*x3^3+7*x2^3*x3^4-14*x2^2*x3^5+7*x2*x3^6-7*x1*x2^5+16*x1^3*x2^2*x3+14*x1*x2^4*x3+16*x1^2*x2^2*x3^2-7*x1*x2^3*x3^2-7*x2^4*x3^2-8*x1^2*x2*x3^3+14*x2^3*x3^3-7*x2^2*x3^4-8*x2^5+16*x2^4*x3-8*x1^2*x2*x3^2+7*x2^3*x3^2-30*x2^2*x3^3+15*x2*x3^4-8*x1*x2^3-16*x2^4+24*x2^3*x3-8*x1^2*x3^2-8*x2^2*x3^2+x3^4+8*x1^3-x1*x2^2+8*x2^3-16*x2^2*x3+7*x2*x3^2+8*x2^2-7*x3^2+8*x2-8",
        }, Q, vars);
        assertEquals(expectedIntersectIJ, I.intersection(J));

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedQuotIJ =
                Ideal.parse(new String[]{
                        "x1+x2-x3",
                        "x2^3-2*x2^2*x3+x2*x3^2-1",
                }, Q, vars);
        assertEquals(expectedQuotIJ, I.quotient(J));
        assertEquals(expectedQuotIJ, I.intersection(J).quotient(J));

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedQuotJI =
                Ideal.parse(new String[]{
                        "x1^3-x3^2+x2-1",
                        "x3^4-x1*x2^2-x2*x3^2+x3^2",
                        "x1^2*x3^2-x2^2"
                }, Q, vars);
        assertEquals(expectedQuotJI, J.quotient(I));
        assertEquals(expectedQuotJI, I.intersection(J).quotient(I));

        assertEquals(I, I.multiply(J).quotient(J));
        assertEquals(J, I.multiply(J).quotient(I));
    }

    @Test
    public void test3() throws Exception {
        String[] vars = {"u0", "u1", "u2", "u3"};

        String[] gensI0 = {"u3*u3 + u2*u2 + u1*u1 + u0*u0 + u1*u1 + u2*u2 + u3*u3 - u0", "u3*0 + u2*0 + u1*u3 + u0*u2 + u1*u1 + u2*u0 + u3*u1 - u2"};
        String[] gensJ0 = {"u3*0 + u2*u3 + u1*u2 + u0*u1 + u1*u0 + u2*u1 + u3*u2 - u1", "u3 + u2 + u1 + u0 + u1 + u2 + u3 - 1"};
        String[] gensK = {"u0*u1^2 - u2^3", "u3^2*u0^2 - u1^3*u2^3"};

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>>
                I0 = Ideal.parse(gensI0, Q, vars),
                J0 = Ideal.parse(gensJ0, Q, vars),
                K = Ideal.parse(gensK, Q, vars),
                I = I0.multiply(K),
                J = J0.multiply(K);

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedI0 = Ideal.parse(new String[]{
                "u1^2+2*u0*u2+2*u1*u3-u2",
                "u0^2-4*u0*u2+2*u2^2-4*u1*u3+2*u3^2-u0+2*u2"
        }, Q, vars);
        assertEquals(expectedI0, I0);

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedI = Ideal.parse(new String[]{
                "u0*u1^4+2*u0^2*u1^2*u2-u1^2*u2^3-2*u0*u2^4+2*u0*u1^3*u3-2*u1*u2^3*u3-u0*u1^2*u2+u2^4",
                "u0^3*u1^2-4*u0^2*u1^2*u2+2*u0*u1^2*u2^2-u0^2*u2^3+4*u0*u2^4-2*u2^5-4*u0*u1^3*u3+4*u1*u2^3*u3+2*u0*u1^2*u3^2-2*u2^3*u3^2-u0^2*u1^2+2*u0*u1^2*u2+u0*u2^3-2*u2^4",
                "u1^5*u2^3+2*u0*u1^3*u2^4+2*u1^4*u2^3*u3-u1^3*u2^4-u0^2*u1^2*u3^2-2*u0^3*u2*u3^2-2*u0^2*u1*u3^3+u0^2*u2*u3^2",
                "u0^2*u1^3*u2^3-4*u0*u1^3*u2^4+2*u1^3*u2^5-4*u1^4*u2^3*u3+2*u1^3*u2^3*u3^2-u0*u1^3*u2^3+2*u1^3*u2^4-u0^4*u3^2+4*u0^3*u2*u3^2-2*u0^2*u2^2*u3^2+4*u0^2*u1*u3^3-2*u0^2*u3^4+u0^3*u3^2-2*u0^2*u2*u3^2",
                "u1^3*u2^6+2*u0*u1*u2^7+2*u1^2*u2^6*u3-u1*u2^7-2*u0^4*u2*u3^2-4*u0^2*u1^2*u2*u3^2+2*u0*u1^2*u2^2*u3^2-u0^2*u2^3*u3^2+4*u0*u2^4*u3^2-2*u2^5*u3^2-2*u0^3*u1*u3^3-4*u0*u1^3*u3^3+4*u1*u2^3*u3^3+2*u0*u1^2*u3^4-2*u2^3*u3^4-u0^2*u1^2*u3^2+u0^3*u2*u3^2+2*u0*u1^2*u2*u3^2+u0*u2^3*u3^2-2*u2^4*u3^2",
                "u0^2*u1*u2^6-4*u0*u1*u2^7+2*u1*u2^8-4*u1^2*u2^6*u3+2*u1*u2^6*u3^2-u0*u1*u2^6+2*u1*u2^7-u0^5*u3^2+4*u0^4*u2*u3^2-2*u0^3*u2^2*u3^2+4*u0^3*u1*u3^3-2*u0^3*u3^4+u0^4*u3^2-2*u0^3*u2*u3^2",
                "2*u1^4*u2^5+5*u0*u1^2*u2^6-4*u1^2*u2^7-4*u2^9+4*u1^3*u2^5*u3+2*u0*u1*u2^6*u3+2*u1^4*u2^3*u3^2+4*u0*u1^2*u2^4*u3^2-4*u2^7*u3^2+4*u1^3*u2^3*u3^3+u0^2*u1^2*u2^4+2*u1^4*u2^4-u1^2*u2^6-u0*u2^7-2*u1*u2^6*u3-u0*u1^2*u2^4+u2^7-u0^4*u1*u3^2+4*u0^3*u1*u2*u3^2-2*u0^2*u1*u2^2*u3^2-2*u0^4*u3^3-4*u0^2*u2^2*u3^3-2*u0^2*u1*u3^4-4*u0^2*u3^5+u0^3*u1*u3^2-2*u0^2*u1*u2*u3^2+2*u0^3*u3^3",
                "u0*u1^2*u2^7-4*u1^2*u2^8-4*u2^10+2*u0*u1*u2^7*u3+2*u1^4*u2^4*u3^2+4*u0*u1^2*u2^5*u3^2-4*u2^8*u3^2+4*u1^3*u2^4*u3^3+u0^2*u1^2*u2^5-5*u0*u1^2*u2^6+5*u1^2*u2^7-u0*u2^8+4*u2^9-4*u1^3*u2^5*u3-2*u0*u1*u2^6*u3-2*u1*u2^7*u3-2*u1^4*u2^3*u3^2-4*u0*u1^2*u2^4*u3^2+4*u2^7*u3^2-4*u1^3*u2^3*u3^3-u0^2*u1^2*u2^4-2*u1^4*u2^4-u0*u1^2*u2^5+u1^2*u2^6+u0*u2^7+u2^8+2*u1*u2^6*u3+3*u0^4*u1*u2*u3^2+8*u0^2*u1^3*u2*u3^2+4*u0^3*u1*u2^2*u3^2-4*u0*u1^3*u2^2*u3^2-8*u0*u1*u2^4*u3^2+4*u1*u2^5*u3^2-2*u0^4*u2*u3^3-8*u0*u1^2*u2^2*u3^3+8*u2^5*u3^3-4*u0*u1^3*u3^4-2*u0^2*u1*u2*u3^4+4*u1*u2^3*u3^4-8*u0*u1^2*u3^5-4*u0^2*u2*u3^5+8*u2^3*u3^5+u0*u1^2*u2^4-u2^7+u0^4*u1*u3^2+2*u0^2*u1^3*u3^2-5*u0^3*u1*u2*u3^2-4*u0*u1^3*u2*u3^2-2*u0*u1*u2^3*u3^2+4*u1*u2^4*u3^2+2*u0^4*u3^3+4*u0^2*u1^2*u3^3+2*u0^3*u2*u3^3+4*u0^2*u2^2*u3^3-4*u0*u2^3*u3^3+2*u0^2*u1*u3^4+4*u0^2*u3^5-u0^3*u1*u3^2+2*u0^2*u1*u2*u3^2-2*u0^3*u3^3",
                "u1^2*u2^9+2*u0*u2^10+2*u1*u2^9*u3-u2^10-2*u0^5*u1*u2*u3^2-14*u0^2*u1^3*u2^2*u3^2-u0^3*u1*u2^3*u3^2+8*u0*u1^3*u2^3*u3^2+14*u0*u1*u2^5*u3^2-8*u1*u2^6*u3^2+36*u0^2*u1^2*u2^2*u3^3-2*u0^3*u2^3*u3^3-36*u0*u2^5*u3^3+2*u0^2*u1^3*u3^4+40*u0*u1^3*u2*u3^4-2*u0*u1*u2^3*u3^4-40*u1*u2^4*u3^4+4*u0^2*u1^2*u3^5-4*u0*u2^3*u3^5+u0^4*u1*u2*u3^2-6*u0^2*u1^3*u2*u3^2+10*u0*u1^3*u2^2*u3^2+6*u0*u1*u2^4*u3^2-10*u1*u2^5*u3^2-12*u0*u1^2*u2^2*u3^3+12*u2^5*u3^3+2*u0*u1^3*u3^4-2*u1*u2^3*u3^4+4*u0*u1^2*u3^5-4*u2^3*u3^5-u0^2*u1^3*u3^2+2*u0*u1^3*u2*u3^2+u0*u1*u2^3*u3^2-2*u1*u2^4*u3^2-2*u0^2*u1^2*u3^3+2*u0*u2^3*u3^3",
                "u0^2*u2^9-4*u0*u2^10+2*u2^11-4*u1*u2^9*u3+2*u2^9*u3^2-u0*u2^9+2*u2^10-u0^6*u1*u3^2+4*u0^5*u1*u2*u3^2-2*u0^4*u1*u2^2*u3^2+56*u0^2*u1^2*u2^2*u3^3+4*u0^3*u2^3*u3^3-32*u0*u1^2*u2^3*u3^3-56*u0*u2^5*u3^3+32*u2^6*u3^3-2*u0^4*u1*u3^4+16*u0^2*u1^3*u3^4+64*u0*u1^3*u2*u3^4-16*u0*u1*u2^3*u3^4-64*u1*u2^4*u3^4-8*u0^2*u1^2*u3^5-32*u0*u1^2*u2*u3^5+8*u0*u2^3*u3^5+32*u2^4*u3^5+u0^5*u1*u3^2-2*u0^4*u1*u2*u3^2+24*u0^2*u1^2*u2*u3^3-40*u0*u1^2*u2^2*u3^3-24*u0*u2^4*u3^3+40*u2^5*u3^3+16*u0*u1^3*u3^4-16*u1*u2^3*u3^4-8*u0*u1^2*u3^5+8*u2^3*u3^5+4*u0^2*u1^2*u3^3-8*u0*u1^2*u2*u3^3-4*u0*u2^3*u3^3+8*u2^4*u3^3"
        }, Q, vars);

        Ideal<Monomial<Rational<BigInteger>>, MultivariatePolynomial<Rational<BigInteger>>> expectedJ = Ideal.parse(new String[]{
                "u0^2*u1^2+2*u0*u1^3+2*u0*u1^2*u2-u0*u2^3-2*u1*u2^3-2*u2^4+2*u0*u1^2*u3-2*u2^3*u3-u0*u1^2+u2^3",
                "4*u0*u1^4+2*u0*u1^3*u2-4*u1^2*u2^3-2*u1*u2^4+4*u0*u1^3*u3-2*u0*u1^2*u2*u3-4*u1*u2^3*u3+2*u2^4*u3-u0*u1^3+u1*u2^3",
                "u0*u1^3*u2^3+2*u1^4*u2^3+2*u1^3*u2^4+2*u1^3*u2^3*u3-u1^3*u2^3-u0^3*u3^2-2*u0^2*u1*u3^2-2*u0^2*u2*u3^2-2*u0^2*u3^3+u0^2*u3^2",
                "u0*u1*u2^6+2*u1^2*u2^6+2*u1*u2^7+2*u1*u2^6*u3-u1*u2^6-u0^4*u3^2-2*u0^3*u1*u3^2-2*u0^3*u2*u3^2-2*u0^3*u3^3+u0^3*u3^2",
                "8*u1^4*u2^4+4*u1^3*u2^5+4*u1^2*u2^6+2*u1*u2^7+8*u1^4*u2^3*u3+2*u0*u1^2*u2^4*u3+16*u1^3*u2^4*u3+4*u1*u2^6*u3-2*u2^7*u3+8*u1^3*u2^3*u3^2-4*u1^4*u2^3-4*u1^3*u2^4-u1*u2^6-6*u1^3*u2^3*u3+u1^3*u2^3-4*u0^3*u1*u3^2-2*u0^3*u2*u3^2-8*u0^2*u1*u2*u3^2-4*u0^2*u2^2*u3^2-4*u0^3*u3^3-8*u0^2*u1*u3^3-16*u0^2*u2*u3^3-8*u0^2*u3^4+u0^3*u3^2+4*u0^2*u1*u3^2+4*u0^2*u2*u3^2+6*u0^2*u3^3-u0^2*u3^2",
                "16*u1^5*u2^3-4*u1^3*u2^5-4*u1^2*u2^6-2*u1*u2^7+8*u1^4*u2^3*u3-2*u0*u1^2*u2^4*u3-24*u1^3*u2^4*u3-4*u1*u2^6*u3+2*u2^7*u3-8*u1^3*u2^3*u3^2+4*u1^3*u2^4+u1*u2^6+6*u1^3*u2^3*u3-u1^3*u2^3+4*u0^3*u1*u3^2+32*u0*u1^3*u3^2+2*u0^3*u2*u3^2+32*u0*u1^2*u2*u3^2+4*u0^2*u2^2*u3^2-16*u0*u2^3*u3^2-32*u1*u2^3*u3^2-32*u2^4*u3^2+4*u0^3*u3^3-8*u0^2*u1*u3^3+32*u0*u1^2*u3^3+24*u0^2*u2*u3^3-32*u2^3*u3^3+8*u0^2*u3^4-u0^3*u3^2-16*u0*u1^2*u3^2-4*u0^2*u2*u3^2+16*u2^3*u3^2-6*u0^2*u3^3+u0^2*u3^2",
                "4*u1^3*u2^6+2*u1^2*u2^7+4*u1^2*u2^6*u3-2*u1*u2^7*u3-u1^2*u2^6-2*u0^3*u1*u2*u3^2-24*u0*u1^3*u2*u3^2-16*u0*u1^2*u2^2*u3^2-4*u0^2*u2^3*u3^2+24*u1*u2^4*u3^2+16*u2^5*u3^2-4*u0^3*u1*u3^3-16*u0*u1^3*u3^3+2*u0^3*u2*u3^3-40*u0*u1^2*u2*u3^3+16*u1*u2^3*u3^3+40*u2^4*u3^3-16*u0*u1^2*u3^4+16*u2^3*u3^4+u0^3*u1*u3^2+12*u0*u1^3*u3^2+16*u0*u1^2*u2*u3^2-12*u1*u2^3*u3^2-16*u2^4*u3^2+16*u0*u1^2*u3^3-16*u2^3*u3^3-4*u0*u1^2*u3^2+4*u2^3*u3^2",
                "u0*u2^9+2*u1*u2^9+2*u2^10+2*u2^9*u3-u2^9-u0^5*u1*u3^2-2*u0^4*u1*u2*u3^2+28*u0*u1^3*u2^2*u3^2-2*u0^3*u2^3*u3^2+16*u0*u1^2*u2^3*u3^2-28*u1*u2^5*u3^2-16*u2^6*u3^2-2*u0^4*u1*u3^3+48*u0*u1^3*u2*u3^3+68*u0*u1^2*u2^2*u3^3-48*u1*u2^4*u3^3-68*u2^5*u3^3+16*u0*u1^3*u3^4+64*u0*u1^2*u2*u3^4-16*u1*u2^3*u3^4-64*u2^4*u3^4+16*u0*u1^2*u3^5-16*u2^3*u3^5+u0^4*u1*u3^2-28*u0*u1^3*u2*u3^2-24*u0*u1^2*u2^2*u3^2+28*u1*u2^4*u3^2+24*u2^5*u3^2-20*u0*u1^3*u3^3-58*u0*u1^2*u2*u3^3+20*u1*u2^3*u3^3+58*u2^4*u3^3-24*u0*u1^2*u3^4+24*u2^3*u3^4+7*u0*u1^3*u3^2+12*u0*u1^2*u2*u3^2-7*u1*u2^3*u3^2-12*u2^4*u3^2+12*u0*u1^2*u3^3-12*u2^3*u3^3-2*u0*u1^2*u3^2+2*u2^3*u3^2",
                "16*u1^2*u2^8+8*u1*u2^9-16*u1^2*u2^7*u3+8*u1*u2^8*u3-8*u2^9*u3-8*u0*u1^2*u2^5*u3^2-16*u1^3*u2^5*u3^2-8*u1*u2^7*u3^2+8*u2^8*u3^2+32*u1^4*u2^3*u3^3+8*u0*u1^2*u2^4*u3^3+32*u1^3*u2^4*u3^3+16*u1*u2^6*u3^3-8*u2^7*u3^3+32*u1^3*u2^3*u3^4+8*u1^2*u2^7+4*u0*u1^2*u2^5*u3+8*u1^3*u2^5*u3-8*u1^2*u2^6*u3+4*u1*u2^7*u3-4*u2^8*u3-48*u1^4*u2^3*u3^2-8*u0*u1^2*u2^4*u3^2-40*u1^3*u2^4*u3^2-20*u1*u2^6*u3^2+8*u2^7*u3^2-56*u1^3*u2^3*u3^3+4*u1^2*u2^6+24*u1^4*u2^3*u3+2*u0*u1^2*u2^4*u3+16*u1^3*u2^4*u3+8*u1*u2^6*u3-2*u2^7*u3-16*u0^3*u1*u2^2*u3^2+128*u0*u1^3*u2^2*u3^2-8*u0^3*u2^3*u3^2+64*u0*u1^2*u2^3*u3^2+36*u1^3*u2^3*u3^2-128*u1*u2^5*u3^2-64*u2^6*u3^2+8*u0^4*u1*u3^3+16*u0^3*u1*u2*u3^3+416*u0*u1^3*u2*u3^3-8*u0^3*u2^2*u3^3+448*u0*u1^2*u2^2*u3^3-416*u1*u2^4*u3^3-448*u2^5*u3^3+192*u0*u1^3*u3^4+8*u0^3*u2*u3^4+608*u0*u1^2*u2*u3^4+16*u0^2*u2^2*u3^4-192*u1*u2^3*u3^4-608*u2^4*u3^4-16*u0^3*u3^5-32*u0^2*u1*u3^5+192*u0*u1^2*u3^5-32*u0^2*u2*u3^5-192*u2^3*u3^5-32*u0^2*u3^6-4*u1^4*u2^3-2*u1^3*u2^4-u1*u2^6-10*u1^3*u2^3*u3-8*u0^3*u1*u2*u3^2-128*u0*u1^3*u2*u3^2-96*u0*u1^2*u2^2*u3^2+128*u1*u2^4*u3^2+96*u2^5*u3^2+8*u0^3*u1*u3^3-176*u0*u1^3*u3^3-4*u0^3*u2*u3^3-384*u0*u1^2*u2*u3^3-8*u0^2*u2^2*u3^3+176*u1*u2^3*u3^3+384*u2^4*u3^3+20*u0^3*u3^4+48*u0^2*u1*u3^4-224*u0*u1^2*u3^4+40*u0^2*u2*u3^4+224*u2^3*u3^4+56*u0^2*u3^5+u1^3*u2^3-4*u0^3*u1*u3^2+32*u0*u1^3*u3^2+48*u0*u1^2*u2*u3^2-32*u1*u2^3*u3^2-48*u2^4*u3^2-8*u0^3*u3^3-24*u0^2*u1*u3^3+80*u0*u1^2*u3^3-16*u0^2*u2*u3^3-80*u2^3*u3^3-36*u0^2*u3^4+u0^3*u3^2+4*u0^2*u1*u3^2-8*u0*u1^2*u3^2+2*u0^2*u2*u3^2+8*u2^3*u3^2+10*u0^2*u3^3-u0^2*u3^2"
        }, Q, vars);

        assertEquals(expectedI, I);
        assertEquals(expectedJ, J);
        assertEquals(I0, I.quotient(K));
        assertEquals(J0, J.quotient(K));
    }

    @Test
    public void test4() throws Exception {
        String[] vars = {"x", "y", "z"};

        Ideal<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> ideal = Ideal.parse(new String[]{
                "x^2 - y*z*x^2 + 2",
                "x^2*y*z^3 - y*z^2 + 2*x^5",
                "x*y*z^3 - y*z^12 + 2*x*y*z^5"
        }, Z, GREVLEX, vars);

        System.out.println(ideal.dimension());
        System.out.println(ideal.degree());

        System.out.println(ideal);
        System.out.println(ideal.changeOrder(LEX));
    }

    @Test
    public void test5() throws Exception {
        String[] vars = {"x", "y", "z"};

        Ideal<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> lex = Ideal.parse(new String[]{
                "x^2 - y*z*x^2 + 2",
                "x^2*y*z^3 - y*z^2 + 2*x^5",
                "x*y*z^3 - y*z^12 + 2*x*y*z^5"
        }, Z, LEX, vars);
        assertEquals(lex, Ideal.create(lex.getOriginalGenerators(), GREVLEX).changeOrder(LEX));

        assertEquals(3, lex.nBasisGenerators());
        assertEquals(0, lex.dimension());
        assertEquals(62, lex.degree());
        assertEquals(UnivariatePolynomial.parse("1+3*x+6*x^2+10*x^3+14*x^4+16*x^5+12*x^6", Q), lex.hilbertSeries().numerator);


        Ideal<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> graded = lex.changeOrder(GREVLEX);
        assertEquals(graded, Ideal.create(lex.getOriginalGenerators(), GREVLEX));

        assertEquals(18, graded.nBasisGenerators());

        Ideal<Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> graded2 = graded.square();

        assertEquals(graded2, graded2.intersection(graded));
        assertEquals(graded, graded2.union(graded));

        assertEquals(0, graded2.dimension());
        assertEquals(248, graded2.degree());
    }
}