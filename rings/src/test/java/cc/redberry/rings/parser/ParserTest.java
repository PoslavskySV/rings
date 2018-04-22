package cc.redberry.rings.parser;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FiniteField;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.MultivariateRing;
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

import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.parser.Parser.mkMultivariateParser;
import static cc.redberry.rings.parser.Parser.mkParser;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class ParserTest extends AbstractTest {
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

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, Function.identity());

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

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, Function.identity());

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

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<MultivariatePolynomial<BigInteger>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, Function.identity());

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

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

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

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

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

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<Rational<MultivariatePolynomial<BigInteger>>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

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
            assertEquals(poly, a.numerator);
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

        Parser<Rational<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?>
                notOptimizedParser = mkParser(baseRing, eVars, null, null, null);

        Parser<Rational<MultivariatePolynomial<Rational<BigInteger>>>, ?, ?>
                optimizedParser = mkParser(baseRing, eVars, polyRing, pVars, p -> new Rational<>(polyRing, p));

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
            assertEquals(poly, a.numerator);
        }

        System.out.println(TimeUnits.statisticsNanotime(notOpt));
        System.out.println(TimeUnits.statisticsNanotime(opt));
    }

    @Test
    public void test7() {
        Parser<UnivariatePolynomial<BigInteger>, ?, ?> parser = Parser.mkUnivariateParser(UnivariateRing(Z), "x");
        UnivariatePolynomial<BigInteger> poly = parser.parse(Tokenizer.mkTokenizer("x + (x^2 + 123*x^3 - 1)^2"));
        assertEquals("1+x-2*x^2-246*x^3+x^4+246*x^5+15129*x^6", poly.toString());
    }

    @Test
    public void test8() {
        FiniteField<UnivariatePolynomialZp64> gf = Rings.GF(17, 3);
        Parser<UnivariatePolynomialZp64, ?, ?> gfParser = Parser.mkUnivariateParser(gf, mkVars(gf, "t"));

        MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> polyRing = Rings.MultivariateRing(3, gf);
        Parser<MultivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> polyParser = Parser.mkMultivariateParser(polyRing, gfParser, "x", "y", "z");

        Rationals<MultivariatePolynomial<UnivariatePolynomialZp64>> fracRing = Rings.Frac(polyRing);
        Parser<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> fracParser =
                Parser.mkNestedParser(
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


        Parser<MultivariatePolynomial<BigInteger>, Monomial<BigInteger>, MultivariatePolynomial<BigInteger>> parser = mkMultivariateParser(polyRing, vars);
        System.out.println(parser.parse("-5*a^22*c*3^13 + 5*a^32*b^24*c*3 + a^31*c*3^42 + c^66").toString(vars));
    }

    @Test
    public void test10() {
        Parser<Rational<BigInteger>, ?, ?> p = Parser.mkParser(Q);
        assertEquals(p.parse("1/2*3"), new Rational<>(Z, Z.valueOf(3), Z.valueOf(2)));
    }

    @Test
    public void test11() {
        Parser<Rational<BigInteger>, ?, ?> p = Parser.mkParser(Q);
        assertEquals(p.parse("-1/2*3/5*9*9/3^2/4/6*7 + 1^2/3 - 3 - 2*3*5/8*9/6*9/3^2/4*7*6^2"), new Rational<>(Z, Z.valueOf(-85879), Z.valueOf(240)));
    }

    private static <P extends IPolynomial<P>> Map<String, P> mkVars(IPolynomialRing<P> ring, String... vars) {
        Map<String, P> map = new HashMap<>();
        for (int i = 0; i < vars.length; i++)
            map.put(vars[i], ring.variable(i));
        return map;
    }
}