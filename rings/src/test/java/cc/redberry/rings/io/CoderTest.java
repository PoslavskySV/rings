package cc.redberry.rings.io;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateQuotientRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.Ideal;
import cc.redberry.rings.poly.multivar.Monomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.*;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class CoderTest extends AbstractTest {
    @Ignore
    @Test
    public void test0_performance() {
        String[] vars = {"x", "y", "z", "x1", "x2", "x3", "x4"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(7, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; ++i)
            pVars.put(vars[i], polyRing.variable(i));

        MultivariateRing<MultivariatePolynomial<BigInteger>> baseRing = polyRing;
        HashMap<String, MultivariatePolynomial<BigInteger>> eVars = pVars;

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, Function.identity());

        MultivariatePolynomial<BigInteger> p = polyRing.parse("x*x2*x4 + 2 + z - y*z*x^2 + x1*x2*x3*x4*x*y*z");
        String expression = polyRing.pow(p, 30).toString(vars);
        System.out.println(expression.length());


        for (int i = 0; i < 1000; ++i) {
            long start = System.nanoTime();
            MultivariatePolynomial<BigInteger> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(expression));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            MultivariatePolynomial<BigInteger> b = optimizedParser.parse(Tokenizer.mkTokenizer(expression));
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));

            start = System.nanoTime();
            MultivariatePolynomial<BigInteger> c = polyRing.parse(expression, vars);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println("\n\n");


            assert a.equals(b);
            assert a.equals(c);
        }
    }

    @Test
    public void test1() {
        String[] vars = new String[]{"x", "y", "z"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        MultivariateRing<MultivariatePolynomial<BigInteger>> baseRing = polyRing;
        HashMap<String, MultivariatePolynomial<BigInteger>> eVars = pVars;

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, Function.identity());

        String expressionStr = "x + y + z*x - z*y";

        MultivariatePolynomial<BigInteger> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        MultivariatePolynomial<BigInteger> b = optimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        assertEquals(a, b);
    }

    @Test
    public void test2() {
        String[] vars = new String[]{"x", "y", "z"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        MultivariateRing<MultivariatePolynomial<BigInteger>> baseRing = polyRing;
        HashMap<String, MultivariatePolynomial<BigInteger>> eVars = pVars;

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, Function.identity());

        String expressionStr = "-x^2*y^2*z^3 + x*y^2 + z^2*x - z*y*x";

        MultivariatePolynomial<BigInteger> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        MultivariatePolynomial<BigInteger> b = optimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        assertEquals(a, b);
    }

    @Test
    public void test3() {
        String[] vars = new String[]{"x", "y", "z"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        Rationals<MultivariatePolynomial<BigInteger>> baseRing = Rings.Frac(polyRing);
        Map<String, Rational<MultivariatePolynomial<BigInteger>>> eVars = pVars.entrySet()
                .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> new Rational<>(polyRing, e.getValue())));

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

        String expressionStr = "z - y/2 + z*y/(3 + 4)";


        Rational<MultivariatePolynomial<BigInteger>> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        Rational<MultivariatePolynomial<BigInteger>> b = optimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        assertEquals(a, b);
    }

    @Test
    public void test4() {
        String[] vars = new String[]{"x", "y", "z"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        Rationals<MultivariatePolynomial<BigInteger>> baseRing = Rings.Frac(polyRing);
        Map<String, Rational<MultivariatePolynomial<BigInteger>>> eVars = pVars.entrySet()
                .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> new Rational<>(polyRing, e.getValue())));

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

        String expressionStr = "(-5942283839318488889)*x^4*y^4";

        Rational<MultivariatePolynomial<BigInteger>> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        Rational<MultivariatePolynomial<BigInteger>> b = optimizedParser.parse(Tokenizer.mkTokenizer(expressionStr));
        assertEquals(a, b);
    }

    @Test
    public void test5() {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);
        String[] vars = new String[]{"x", "y", "z", "p", "q"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);

        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        Rationals<MultivariatePolynomial<BigInteger>> baseRing = Rings.Frac(polyRing);
        Map<String, Rational<MultivariatePolynomial<BigInteger>>> eVars = pVars.entrySet()
                .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> new Rational<>(polyRing, e.getValue())));

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

        DescriptiveStatistics
                notOpt = new DescriptiveStatistics(),
                opt = new DescriptiveStatistics();
        for (int i = 0; i < 100; ++i) {
            MultivariatePolynomial<BigInteger> poly = polyRing.randomElement(100, 100, rnd);

            if (rnd.nextBoolean())
                poly.increment();
            else if (rnd.nextBoolean())
                poly.decrement();
            String stringExpr = poly.toString(vars);

            long start = System.nanoTime();
            Rational<MultivariatePolynomial<BigInteger>> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(stringExpr));
            notOpt.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            Rational<MultivariatePolynomial<BigInteger>> b = optimizedParser.parse(Tokenizer.mkTokenizer(stringExpr));
            opt.addValue(System.nanoTime() - start);

            assertEquals(a, b);
            assertTrue(a.isIntegral());
            assertEquals(poly, a.numerator());
        }

        System.out.println(TimeUnits.statisticsNanotime(notOpt));
        System.out.println(TimeUnits.statisticsNanotime(opt));
    }

    @Test
    public void test6() {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);
        String[] vars = new String[]{"x", "y", "z", "p", "q"};
        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> polyRing = Rings.MultivariateRing(vars.length, Rings.Q);

        HashMap<String, MultivariatePolynomial<Rational<BigInteger>>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));

        Rationals<MultivariatePolynomial<Rational<BigInteger>>> baseRing = Rings.Frac(polyRing);
        Map<String, Rational<MultivariatePolynomial<Rational<BigInteger>>>> eVars = pVars.entrySet()
                .stream().collect(Collectors.toMap(Map.Entry::getKey, e -> new Rational<>(polyRing, e.getValue())));

        Coder<Rational<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?>
                notOptimizedParser = Coder.mkCoder(baseRing, eVars, null, null, null);

        Coder<Rational<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?>
                optimizedParser = Coder.mkCoder(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

        DescriptiveStatistics
                notOpt = new DescriptiveStatistics(),
                opt = new DescriptiveStatistics();
        for (int i = 0; i < 100; ++i) {
            MultivariatePolynomial<Rational<BigInteger>> poly = polyRing.randomElement(100, 100, rnd);

            if (rnd.nextBoolean())
                poly.increment();
            else if (rnd.nextBoolean())
                poly.decrement();
            String stringExpr = poly.toString(vars);

            long start = System.nanoTime();
            Rational<MultivariatePolynomial<Rational<BigInteger>>> a = notOptimizedParser.parse(Tokenizer.mkTokenizer(stringExpr));
            notOpt.addValue(System.nanoTime() - start);

            start = System.nanoTime();
            Rational<MultivariatePolynomial<Rational<BigInteger>>> b = optimizedParser.parse(Tokenizer.mkTokenizer(stringExpr));
            opt.addValue(System.nanoTime() - start);

            assertEquals(a, b);
            assertTrue(a.isIntegral());
            assertEquals(poly, a.numerator());
        }

        System.out.println(TimeUnits.statisticsNanotime(notOpt));
        System.out.println(TimeUnits.statisticsNanotime(opt));
    }

    @Test
    public void test7() {
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> parser = Coder.mkUnivariateCoder(UnivariateRing(Z), "x");
        UnivariatePolynomial<BigInteger> poly = parser.parse(Tokenizer.mkTokenizer("x + (x^2 + 123*x^3 - 1)^2"));
        assertEquals("1+x-2*x^2-246*x^3+x^4+246*x^5+15129*x^6", poly.toString());
    }

    @Test
    public void test8() {
        FiniteField<UnivariatePolynomialZp64> gf = Rings.GF(17, 3);
        Coder<UnivariatePolynomialZp64, ?, ?> gfParser = Coder.mkUnivariateCoder(gf, mkVars(gf, "t"));

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> polyRing = Rings.MultivariateRing(3, gf);
        Coder<MultivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> polyParser = Coder.mkMultivariateCoder(polyRing, gfParser, "x", "y", "z");

        Rationals<MultivariatePolynomial<UnivariatePolynomialZp64>> fracRing = Rings.Frac(polyRing);
        Coder<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> fracParser =
                Coder.mkNestedCoder(
                        fracRing,
                        new HashMap<>(),
                        polyParser,
                        p -> new Rational<>(polyRing, p));

        Rational<MultivariatePolynomial<UnivariatePolynomialZp64>> rat = fracParser.parse(Tokenizer.mkTokenizer("x*(t + 1 - 1/t^99) + t/y"));
        System.out.println(rat);
    }

    @Test
    public void test9() {
        String[] vars = new String[]{"a", "b", "c"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> polyRing = Rings.MultivariateRing(vars.length, Z);
        HashMap<String, MultivariatePolynomial<BigInteger>> pVars = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            pVars.put(vars[i], polyRing.variable(i));


        Coder<MultivariatePolynomial<BigInteger>, Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> parser = Coder.mkMultivariateCoder(polyRing, vars);
        System.out.println(parser.parse("-5*a^22*c*3^13 + 5*a^32*b^24*c*3 + a^31*c*3^42 + c^66").toString(vars));
    }

    @Test
    public void test10() {
        Coder<Rational<BigInteger>, ?, ?> p = Coder.mkCoder(Q);
        assertEquals(p.parse("1/2*3"), new Rational<>(Z, Z.valueOf(3), Z.valueOf(2)));
    }

    @Test
    public void test11() {
        Coder<Rational<BigInteger>, ?, ?> p = Coder.mkCoder(Q);
        assertEquals(p.parse("-1/2*3/5*9*9/3^2/4/6*7 + 1^2/3 - 3 - 2*3*5/8*9/6*9/3^2/4*7*6^2"), new Rational<>(Z, Z.valueOf(-85879), Z.valueOf(240)));
    }

    @Test
    public void test12_random() {
        // complicated ring

        // Gf(17, 3) as polys over "t"
        FiniteField<UnivariatePolynomialZp64> gf17p3 = Rings.GF(17, 3);
        Coder<UnivariatePolynomialZp64, ?, ?> gfCoder = Coder.mkPolynomialCoder(gf17p3, "t");

        // Gf(17, 3)[x, y, z]
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> mRing = Rings.MultivariateRing(3, gf17p3);
        Coder<MultivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> mCoder = Coder.mkMultivariateCoder(mRing, gfCoder, "x", "y", "z");
        //System.out.println(mRing.randomElementTree().toString(mCoder));

        // Frac(Gf(17, 3)[x, y, z])
        Rationals<MultivariatePolynomial<UnivariatePolynomialZp64>> fRing = Rings.Frac(mRing);
        Coder<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> fCoder = Coder.mkRationalsCoder(fRing, mCoder);
        //System.out.println(fRing.randomElementTree().toString(fCoder));

        // Frac(Gf(17, 3)[x, y, z])[W]
        UnivariateRing<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>> uRing = Rings.UnivariateRing(fRing);
        Coder<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>, ?, ?> wCoder = Coder.mkUnivariateCoder(uRing, fCoder, "W");
        //System.out.println(uRing.randomElementTree().toString(wCoder));

        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 100); ++i)
            assertDecode(wCoder, uRing.randomElementTree(rnd));
    }

    @Test
    public void test13_random() {
        // complicated ring

        MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>> mRing = Rings.MultivariateRing(3, Rings.Q);
        Coder<MultivariatePolynomial<Rational<BigInteger>>, ?, ?> mCoder = Coder.mkMultivariateCoder(mRing, "x1", "xx1", "xxx123");
        //System.out.println(mRing.randomElementTree().toString(mCoder));

        Rationals<MultivariatePolynomial<Rational<BigInteger>>> fRing = Rings.Frac(mRing);
        Coder<Rational<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?> fCoder = Coder.mkRationalsCoder(fRing, mCoder);
        //System.out.println(fRing.randomElementTree().toString(fCoder));

        UnivariateRing<UnivariatePolynomial<Rational<MultivariatePolynomial<Rational<BigInteger>>>>> uRing = Rings.UnivariateRing(fRing);
        Coder<UnivariatePolynomial<Rational<MultivariatePolynomial<Rational<BigInteger>>>>, ?, ?> wCoder = Coder.mkUnivariateCoder(uRing, fCoder, "W");
        //System.out.println(uRing.randomElementTree().toString(wCoder));

        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(10, 10); ++i)
            assertDecode(wCoder, uRing.randomElementTree(rnd));
    }

    @Test
    public void test14_random() {
        // complicated ring

        UnivariateRing<UnivariatePolynomial<BigInteger>> uRing = Rings.UnivariateRing(Rings.Z);
        Coder<UnivariatePolynomial<BigInteger>, ?, ?> uCoder = Coder.mkUnivariateCoder(uRing, "p");
        //System.out.println(uRing.randomElementTree().toString(uCoder));

        Rationals<UnivariatePolynomial<BigInteger>> fRing = Rings.Frac(uRing);
        Coder<Rational<UnivariatePolynomial<BigInteger>>, ?, ?> fCoder = Coder.mkRationalsCoder(fRing, uCoder);
        //System.out.println(fRing.randomElementTree().toString(fCoder));

        MultivariateRing<MultivariatePolynomial<Rational<UnivariatePolynomial<BigInteger>>>> mRing = Rings.MultivariateRing(2, fRing);
        Coder<MultivariatePolynomial<Rational<UnivariatePolynomial<BigInteger>>>, ?, ?> mCoder = Coder.mkMultivariateCoder(mRing, fCoder, "x1", "x2");
        //System.out.println(mRing.randomElementTree().toString(mCoder));

        RandomGenerator rnd = getRandom();
        for (int i = 0; i < its(100, 100); ++i)
            assertDecode(mCoder, mRing.randomElementTree(rnd));
    }

    @Test
    public void test15_entities() {
        // complicated ring

        UnivariatePolynomialZp64 gen = Coder
                .mkPolynomialCoder(Rings.UnivariateRingZp64(17), "t")
                .parse("6+7*t+14*t^2+t^3");

        // Gf(17, 3) as polys over "t"
        FiniteField<UnivariatePolynomialZp64> gf17p3 = Rings.GF(gen);
        Coder<UnivariatePolynomialZp64, ?, ?> gfCoder = Coder.mkPolynomialCoder(gf17p3, "t");
        assertEquals("(Z/17)[t]/<6+7*t+14*t^2+t^3>", gf17p3.toString(gfCoder));

        // Gf(17, 3)[x, y, z]
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> mRing = Rings.MultivariateRing(3, gf17p3);
        Coder<MultivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> mCoder = Coder.mkMultivariateCoder(mRing, gfCoder, "x", "y", "z");
        assertEquals("((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x, y, z]", mRing.toString(mCoder));

        // Frac(Gf(17, 3)[x, y, z])
        Rationals<MultivariatePolynomial<UnivariatePolynomialZp64>> fRing = Frac(mRing);
        Coder<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> fCoder = Coder.mkRationalsCoder(fRing, mCoder);
        assertEquals("Frac(((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x, y, z])", fRing.toString(fCoder));

        // Frac(Gf(17, 3)[x, y, z])[W]
        UnivariateRing<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>> wRing = Rings.UnivariateRing(fRing);
        Coder<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>, ?, ?> wCoder = Coder.mkUnivariateCoder(wRing, fCoder, "W");
        assertEquals("(Frac(((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x, y, z]))[W]", wRing.toString(wCoder));


        //x,y,z,t,W
        List<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>> wPolys =
                Arrays.asList(
                        wCoder.parse("x + y - z - t"),
                        wCoder.parse("(x + t * y / z - W)^2"),
                        wCoder.parse("z + y + W^2 + t^5 / x - x^4"),
                        wCoder.parse("1 + x / t + y + t * x^2 + t * y^2 + z*W^5"));

        PolynomialFactorDecomposition<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>
                wFactors = wRing.factor(wRing.multiply(wPolys));
        assertEquals(4, wFactors.sumExponents());
        assertEquals(wFactors.multiply(), wCoder.parse(wFactors.toString(wCoder)));


        // Frac(Gf(17, 3)[x, y, z])[W][u,v]
        MultivariateRing<MultivariatePolynomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>> uvRing = Rings.MultivariateRing(2, wRing);
        Coder<MultivariatePolynomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>, ?, ?> uvCoder = Coder.mkMultivariateCoder(uvRing, wCoder, "u", "v");
        assertEquals("((Frac(((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x, y, z]))[W])[u, v]", uvRing.toString(uvCoder));

        List<MultivariatePolynomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>> uvPolys = Arrays.asList(
                uvCoder.parse("u^2 + v^2 - 1/z - t"),
                uvCoder.parse("(u + t * v / z - W)^2"),
                uvCoder.parse("z + v + W^2 + t^5 / x - W^4"),
                uvCoder.parse("1 + x / t + u^3 + t * v^2 + t * u^2 + z*W^5"));

        PolynomialFactorDecomposition<MultivariatePolynomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>>
                uvFactors = uvRing.factor(uvRing.multiply(uvPolys));
        assertEquals(uvFactors.multiply(), uvCoder.parse(uvFactors.toString(uvCoder)));


        uvPolys = Arrays.asList(
                uvCoder.parse("t * u^2 + v^2"),
                uvCoder.parse("u^2 - t * v^2"));
        Ideal<Monomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>, MultivariatePolynomial<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>>
                ideal = Ideal.create(uvPolys);
        assertEquals("<v^2, u^2>", ideal.toString(uvCoder));

        UnivariateQuotientRing<UnivariatePolynomial<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>>>
                wQuot = Rings.UnivariateQuotientRing(wRing, wCoder.parse("1 + W^2 + W^4"));

        assertDecode(wCoder, wQuot.valueOf(wCoder.parse("1 + x / t + y + t * x^2 + t * y^2 + z*W^116")));
        assertDecode(wCoder, wQuot.valueOf(wCoder.parse("1 + x / t + y + t * x^2 + t * y^2 + z*W^1161 + W^22")));
        assertEquals("(Frac(((Z/17)[t]/<6+7*t+14*t^2+t^3>)[x, y, z]))[W]/<1+W^2+W^4>", wQuot.toString(wCoder));
    }

    private static <E> void assertDecode(Coder<E, ?, ?> coder, E element) {
        String encoded = coder.encode(element);
        E decoded = coder.decode(encoded);
        assertEquals(encoded, element, decoded);
    }

    private static <P extends IPolynomial<P>> Map<String, P> mkVars(IPolynomialRing<P> ring, String... vars) {
        Map<String, P> map = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            map.put(vars[i], ring.variable(i));
        return map;
    }
}