package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Coder;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.redberry.rings.Rings.*;

/**
 * @since 1.0
 */
public class RationalsTest extends AbstractTest {
    @Test
    public void test1() throws Exception {
        for (int i = 19; i < 1000; i++)
            Assert.assertNotNull(Frac(UnivariateRingZp64(17)).randomElement());
    }

    @Test
    public void test2() throws Exception {
        Rational<BigInteger> rat1 = new Rational<>(Z, Z.valueOf(2), Z.valueOf(3));
        Rational<BigInteger> rat2 = new Rational<>(Z, Z.valueOf(0), Z.valueOf(13));
        Assert.assertEquals(rat2, rat1.subtract(rat1));
        Rational<BigInteger> rat3 = new Rational<>(Z, Z.valueOf(2), Z.valueOf(2));
        Rational<BigInteger> rat4 = new Rational<>(Z, Z.valueOf(13), Z.valueOf(13));
        Assert.assertEquals(rat3, rat4);
    }

    @Test
    public void test3() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = MultivariateRing(7, Z);
        Rationals<MultivariatePolynomial<BigInteger>> fracRing = Frac(polyRing);
        String[] vars = {"x", "y", "z", "a", "b", "c", "d"};
        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?> coder = Coder.mkRationalsCoder(fracRing, Coder.mkMultivariateCoder(polyRing, vars));


        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);
        String[] els = {"x", "y", "z", "a", "b", "c", "d",
                "(a+b/c)", "(c-d/c)", "(a-x/y)", "(x -y^2)",
                "x^3", "y^3", "z^3", "a^3", "b^3", "c^3", "d^3",
                "x^4", "y^4", "z^4", "a^4", "b^4", "c^4", "d^4"};

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 100; ++i) {
            if (i != 0)
                sb.append(rnd.nextBoolean() ? "+" : "-");

            for (int j = 0; j < 1 + rnd.nextInt(15); ++j) {
                if (j != 0)
                    sb.append(rnd.nextBoolean() ? "*" : "/");
                sb.append(els[rnd.nextInt(els.length)]);
            }
        }

        String exprString = sb.toString();
        for (int i = 0; i < 10; ++i) {
            long start = System.nanoTime();
            Rational<MultivariatePolynomial<BigInteger>> r = coder.parse(exprString);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(2467 + 24, r.numerator().size() + r.denominator().size());
        }

        // 104ms
        // 92ms
        // 85ms
        // 91ms
        // 87ms
        // 84ms
    }

    @Test
    public void test4() {
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = MultivariateRing(7, Z);
        Rationals<MultivariatePolynomial<BigInteger>> fracRing = Frac(polyRing);
        String[] vars = {"x", "y", "z", "a", "b", "c", "d"};
        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?> coder = Coder.mkRationalsCoder(fracRing, Coder.mkMultivariateCoder(polyRing, vars));


        RandomGenerator rnd = getRandom();
        String[] els = {"x", "y", "z", "a", "b", "c", "d",
                "(a+b/c)", "(c-d/c)", "(a-x/y)", "(x -y^2)",
                "x^3", "y^3", "z^3", "a^3", "b^3", "c^3", "d^3",
                "x^4", "y^4", "z^4", "a^4", "b^4", "c^4", "d^4"};

        for (int n = 0; n < 10; ++n) {
            rnd.setSeed(n);
            System.out.println("=> " + n);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < 100; ++i) {
                if (i != 0)
                    sb.append(rnd.nextBoolean() ? "+" : "-");

                for (int j = 0; j < 1 + rnd.nextInt(15); ++j) {
                    if (j != 0)
                        sb.append(rnd.nextBoolean() ? "*" : "/");
                    sb.append(els[rnd.nextInt(els.length)]);
                }
            }

            String exprString = sb.toString();

            Rational<MultivariatePolynomial<BigInteger>> r = coder.parse(exprString);

            String numString = exprString;
            MultivariatePolynomial<BigInteger>
                    num = r.numerator(),
                    den = r.denominator();
            for (int i = 0; i < vars.length; i++) {
                numString = numString.replace(vars[i], Integer.toString(i + 1));
                num = num.evaluate(i, i + 1);
                den = den.evaluate(i, i + 1);
            }

            Assert.assertEquals(new Rational<>(polyRing, num, den), coder.parse(numString));
        }
    }
}