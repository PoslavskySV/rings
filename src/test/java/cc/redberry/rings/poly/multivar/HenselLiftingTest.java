package cc.redberry.rings.poly.multivar;

import cc.redberry.combinatorics.Combinatorics;
import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.HenselLifting.Evaluation;
import cc.redberry.rings.poly.multivar.HenselLifting.IEvaluation;
import cc.redberry.rings.poly.multivar.HenselLifting.lEvaluation;
import cc.redberry.rings.poly.multivar.MultivariateFactorization.IEvaluationLoop;
import cc.redberry.rings.poly.multivar.MultivariateFactorization.lEvaluationLoop;
import cc.redberry.rings.poly.test.FactorizationInput;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.test.AbstractTest;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.IntStream;
import java.util.stream.StreamSupport;

import static cc.redberry.rings.poly.multivar.MultivariateFactorizationTest.*;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64.parse;

/**
 * @since 1.0
 */
public class HenselLiftingTest {

    @Test
    public void test1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        parse("a^15*b + a^15 - 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                        parse("a^5*b^6 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars)
                },
                base = multiply(factors);

        assert MultivariateGCD.PolynomialGCD(factors).isConstant();

        lEvaluation evaluation = new lEvaluation(base.nVariables, new long[]{1}, domain, base.ordering);

        MultivariatePolynomialZp64[] uFactors = evaluation.evaluateFrom(factors, 1);
        MultivariatePolynomialZp64[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(MultivariatePolynomialZp64[]::new);
        HenselLifting.multivariateLift0(base,
                uFactors,
                factorsLC,
                evaluation,
                base.degrees());
        Assert.assertEquals(base, multiply(uFactors));
    }

    @Test
    public void test2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(43313);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        parse("36045*b^2+23621*a*b+21517*a^2", domain, vars),
                        parse("12894*b^3+22166*a+31033*a*b^2+25906*a^2*b^3", domain, vars),
                        parse("2387*b+11677*a+14775*a^2+25925*a^2*b", domain, vars),
                        parse("1+17708*a*b^2+7251*a*b^5+12898*a^2*b^2+12277*a^3*b^4+23269*a^5*b", domain, vars),
                        parse("27799+34918*a+25070*a^2+2145*a^2*b", domain, vars),
                },
                base = multiply(factors);

        assert MultivariateGCD.PolynomialGCD(factors).isConstant();

        lEvaluation evaluation = new lEvaluation(base.nVariables, new long[]{1146}, domain, base.ordering);

        MultivariatePolynomialZp64[] uFactors = evaluation.evaluateFrom(factors, 1);
        System.out.println(allCoprime(uFactors));
        System.out.println(IntStream.range(0, uFactors.length).allMatch(i -> factors[i].degree(0) == uFactors[i].degree(0)));

        MultivariatePolynomialZp64[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(MultivariatePolynomialZp64[]::new);
        HenselLifting.multivariateLift0(base,
                uFactors,
                factorsLC,
                evaluation,
                base.degrees());
        Assert.assertEquals(base, multiply(uFactors));
    }


    @Test
    public void testEvaluation1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(97);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                poly = parse("b*a^2 + b + a^2 + 2 + a^3*b^4 - a^62*b + 3*b^55*a^55 + b^66 + 3*c^55 + c*a*b + 3", domain, vars);
        lEvaluation evaluation = new lEvaluation(poly.nVariables, new long[]{2, 3}, domain, poly.ordering);

        Assert.assertEquals(evaluation.evaluate(poly, 1), evaluation.modImage(poly, 1, 1));
        Assert.assertEquals(parse("53 + a^2 + 49*a^3 + 22*a^55 + 2*b + a^2*b + 32*a^3*b + 84*a^55*b + 96*a^62*b + a*b*c + 3*c^55", domain, vars),
                evaluation.modImage(poly, 1, 2));
        Assert.assertEquals(parse("31 + a^2 + 26*a^55 + 77*b + a^2*b + 93*a^55*b + 96*a^62*b + 96*b^2 + 60*a^55*b^2 + 52*b^3 + 4*a^55*b^3 + 6*b^4 + a^3*b^4 + 27*a^55*b^4 + 15*b^5 + 88*a^55*b^5 + 25*b^6 + 75*a^55*b^6 + a*b*c + 3*c^55", domain, vars),
                evaluation.modImage(poly, 1, 7));
        Assert.assertEquals(parse("89 + a^2 + 22*a^55 + 33*b + a^2*b + 81*a^55*b + 96*a^62*b + 83*b^2 + 26*a^55*b^2 + 54*b^3 + 32*a^55*b^3 + 52*b^4 + a^3*b^4 + 89*a^55*b^4 + 14*b^5 + 56*a^55*b^5 + 19*b^6 + 8*a^55*b^6 + 58*b^7 + 9*a^55*b^7 + 32*b^8 + 23*a^55*b^8 + 57*b^9 + 42*a^55*b^9 + 58*b^10 + 90*a^55*b^10 + 56*b^11 + 53*a^55*b^11 + 24*b^12 + 4*a^55*b^12 + 86*b^13 + 92*a^55*b^13 + 33*b^14 + 64*a^55*b^14 + 71*b^15 + 71*a^55*b^15 + 44*b^16 + 42*a^55*b^16 + 49*b^17 + 92*a^55*b^17 + 26*b^18 + 94*a^55*b^18 + 77*b^19 + 78*a^55*b^19 + 76*b^20 + 24*a^55*b^20 + 37*b^21 + 76*a^55*b^21 + 3*b^22 + 64*a^55*b^22 + 80*b^23 + 71*a^55*b^23 + 77*b^24 + 89*a^55*b^24 + 71*b^25 + 7*a^55*b^25 + 92*b^26 + 87*a^55*b^26 + 65*b^27 + 49*a^55*b^27 + 25*b^28 + 60*a^55*b^28 + 30*b^29 + 19*a^55*b^29 + 95*b^30 + 16*a^55*b^30 + 32*b^31 + 3*a^55*b^31 + 85*b^32 + 62*a^55*b^32 + 24*b^33 + 91*a^55*b^33 + 32*b^34 + 55*a^55*b^34 + 79*b^35 + 58*a^55*b^35 + 57*b^36 + 82*a^55*b^36 + 30*b^37 + 62*a^55*b^37 + 44*b^38 + 64*a^55*b^38 + 61*b^39 + 95*a^55*b^39 + 68*b^40 + 17*a^55*b^40 + 57*b^41 + 35*a^55*b^41 + 76*b^42 + 96*a^55*b^42 + 2*b^43 + 29*a^55*b^43 + 83*b^44 + 37*a^55*b^44 + 42*b^45 + 92*a^55*b^45 + 3*b^46 + 42*a^55*b^46 + a*b*c + 3*c^55", domain, vars),
                evaluation.modImage(poly, 1, 47));


        Assert.assertEquals(evaluation.evaluate(poly, 2), evaluation.modImage(poly, 2, 1));
        Assert.assertEquals(parse("52 + a^2 + b + a^2*b + 96*a^62*b + a^3*b^4 + 3*a^55*b^55 + b^66 + 5*c + a*b*c", domain, vars),
                evaluation.modImage(poly, 2, 2));
        Assert.assertEquals(parse("87 + a^2 + b + a^2*b + 96*a^62*b + a^3*b^4 + 3*a^55*b^55 + b^66 + 9*c + a*b*c + 7*c^2 + 93*c^3 + 79*c^4 + 4*c^5 + 64*c^6", domain, vars),
                evaluation.modImage(poly, 2, 7));
        Assert.assertEquals(parse("52 + a^2 + b + a^2*b + 96*a^62*b + a^3*b^4 + 3*a^55*b^55 + b^66 + 36*c + a*b*c + 58*c^2 + 65*c^3 + 70*c^4 + 29*c^5 + 12*c^6 + 9*c^7 + 80*c^8 + 51*c^9 + 59*c^10 + 44*c^11 + 62*c^12 + 13*c^13 + 96*c^14 + 71*c^15 + 28*c^16 + 84*c^17 + 53*c^18 + 19*c^19 + 81*c^20 + 74*c^21 + 96*c^22 + 71*c^23 + 27*c^24 + 57*c^25 + 15*c^26 + 48*c^27 + 57*c^28 + 67*c^29 + 24*c^30 + 3*c^31 + 9*c^32 + 62*c^33 + 63*c^34 + 39*c^35 + 10*c^36 + 91*c^37 + 96*c^38 + 95*c^39 + 76*c^40 + 91*c^41 + 50*c^42 + 68*c^43 + 40*c^44 + 13*c^45 + 63*c^46", domain, vars),
                evaluation.modImage(poly, 2, 47));
    }

    @Test
    public void testEvaluation2() throws Exception {
        // small characteristic
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                poly = parse("b*a^2 + b + a^2 + 2 + a^3*b^4 - a^62*b + 3*b^55*a^55 + b^66 + 3*c^55 + c*a*b + 3", domain, vars);
        lEvaluation evaluation = new lEvaluation(poly.nVariables, new long[]{1, 3}, domain, poly.ordering);

        Assert.assertEquals(evaluation.evaluate(poly, 1), evaluation.modImage(poly, 1, 1));
        Assert.assertEquals(parse("a^2 + a^3 + b + a^2*b + a^55*b + a^62*b + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 2));
        Assert.assertEquals(parse("1 + a^2 + a^55 + b + a^2*b + a^55*b + a^62*b + b^2 + a^55*b^2 + a^55*b^3 + a^3*b^4 + a^55*b^4 + a^55*b^5 + a^55*b^6 + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 7));
        Assert.assertEquals(parse("1 + a^2 + b + a^2*b + a^62*b + b^2 + a^3*b^4 + a^55*b^7 + a^55*b^23 + a^55*b^39 + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 47));


        Assert.assertEquals(evaluation.evaluate(poly, 2), evaluation.modImage(poly, 2, 1));
        Assert.assertEquals(parse("1 + a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + c + a*b*c", domain, vars),
                evaluation.modImage(poly, 2, 2));
        Assert.assertEquals(parse("a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + c + a*b*c + c^2 + c^3 + c^4 + c^5 + c^6", domain, vars),
                evaluation.modImage(poly, 2, 7));
        Assert.assertEquals(parse("1 + a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + a*b*c + c^7 + c^23 + c^39", domain, vars),
                evaluation.modImage(poly, 2, 47));
    }

    @Test
    public void testEvaluation3() throws Exception {
        // small characteristic
        IntegersZp domain = new IntegersZp(2);
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>
                poly = MultivariatePolynomial.parse("b*a^2 + b + a^2 + 2 + a^3*b^4 - a^62*b + 3*b^55*a^55 + b^66 + 3*c^55 + c*a*b + 3", domain, vars);
        Evaluation<BigInteger> evaluation = new Evaluation<>(poly.nVariables, new BigInteger[]{BigInteger.ONE, BigInteger.ONE}, domain, poly.ordering);

        Assert.assertEquals(evaluation.evaluate(poly, 1), evaluation.modImage(poly, 1, 1));
        Assert.assertEquals(MultivariatePolynomial.parse("a^2 + a^3 + b + a^2*b + a^55*b + a^62*b + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 2));
        Assert.assertEquals(MultivariatePolynomial.parse("1 + a^2 + a^55 + b + a^2*b + a^55*b + a^62*b + b^2 + a^55*b^2 + a^55*b^3 + a^3*b^4 + a^55*b^4 + a^55*b^5 + a^55*b^6 + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 7));
        Assert.assertEquals(MultivariatePolynomial.parse("1 + a^2 + b + a^2*b + a^62*b + b^2 + a^3*b^4 + a^55*b^7 + a^55*b^23 + a^55*b^39 + a*b*c + c^55", domain, vars),
                evaluation.modImage(poly, 1, 47));


        Assert.assertEquals(evaluation.evaluate(poly, 2), evaluation.modImage(poly, 2, 1));
        Assert.assertEquals(MultivariatePolynomial.parse("1 + a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + c + a*b*c", domain, vars),
                evaluation.modImage(poly, 2, 2));
        Assert.assertEquals(MultivariatePolynomial.parse("a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + c + a*b*c + c^2 + c^3 + c^4 + c^5 + c^6", domain, vars),
                evaluation.modImage(poly, 2, 7));
        Assert.assertEquals(MultivariatePolynomial.parse("1 + a^2 + b + a^2*b + a^62*b + a^3*b^4 + a^55*b^55 + b^66 + a*b*c + c^7 + c^23 + c^39", domain, vars),
                evaluation.modImage(poly, 2, 47));
    }

    @Test
    public void testHenselLiftingRandom1() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource baseSource = new MultivariateFactorizationTest.lSampleDecompositionSource(
                3, 5,
                2, 4,
                2, 6,
                1, 4);
        baseSource.minModulusBits = 15;
        baseSource.maxModulusBits = 30;
        FactorizationInput.SampleDecompositionSource<MultivariatePolynomialZp64> source
                = orderVarsByDegree(
                filterMonomialContent(
                        filterNonSquareFree(baseSource)));
        testHenselLift(source, AbstractTest.its(100, 1000), lEvaluationLoop::new, HenselLifting::multivariateLift0, true, 2);
        testHenselLift(source, AbstractTest.its(100, 1000), lEvaluationLoop::new, HenselLifting::multivariateLift0, true, 1);
    }

    @Benchmark(runAnyway = true)
    @Test
    public void testBivariateLifting1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(67);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                base = parse("33*b+34*b^2+36*b^3+59*b^4+50*b^5+52*b^6+66*b^7+17*b^8+33*a^4+34*a^4*b+36*a^4*b^2+28*a^4*b^3+14*a^4*b^4+6*a^4*b^5+57*a^4*b^6+17*a^4*b^7+52*a^4*b^8+66*a^4*b^9+17*a^4*b^10+11*a^5*b+59*a^5*b^2+50*a^5*b^3+53*a^5*b^4+65*a^5*b^5+45*a^5*b^6+56*a^5*b^7+53*a^5*b^8+a^5*b^9+66*a^5*b^10+17*a^5*b^11+33*a^6*b^2+11*a^6*b^4+3*a^6*b^5+a^6*b^7+3*a^8*b^2+64*a^8*b^3+52*a^8*b^4+2*a^8*b^5+14*a^8*b^6+52*a^8*b^7+66*a^8*b^8+17*a^8*b^9+11*a^9+59*a^9*b+50*a^9*b^2+65*a^9*b^3+56*a^9*b^4+45*a^9*b^5+42*a^9*b^6+51*a^9*b^7+47*a^9*b^8+54*a^9*b^9+20*a^9*b^10+a^9*b^11+66*a^9*b^12+17*a^9*b^13+33*a^10*b+a^10*b^2+10*a^10*b^3+56*a^10*b^4+13*a^10*b^6+4*a^10*b^7+66*a^10*b^8+18*a^10*b^9+11*a^11*b^2+3*a^11*b^3+2*a^11*b^5+11*a^11*b^7+a^11*b^10+a^13*b^2+66*a^13*b^3+17*a^13*b^4+a^13*b^5+66*a^13*b^6+18*a^13*b^7+66*a^13*b^8+17*a^13*b^9+a^13*b^10+66*a^13*b^11+17*a^13*b^12+a^14*b+66*a^14*b^2+20*a^14*b^3+a^14*b^4+21*a^14*b^6+66*a^14*b^7+18*a^14*b^8+a^14*b^9+66*a^14*b^10+17*a^14*b^11+11*a^15*b+3*a^15*b^2+14*a^15*b^4+3*a^15*b^5+11*a^15*b^6+2*a^15*b^7+13*a^15*b^9+a^15*b^12+a^16*b^3+a^16*b^8+a^19*b^3+a^19*b^6+a^19*b^8+a^19*b^11+a^20*b^2+a^20*b^5+a^20*b^7+a^20*b^10", domain);

        MultivariatePolynomialZp64[] uFactors = {
                parse("b^2+b^5+b^7+b^10", domain, vars),
                parse("33+a", domain, vars),
                parse("22+a", domain, vars),
                parse("32+8*a+a^2", domain, vars),
                parse("32+59*a+a^2", domain, vars),
                parse("24+5*a+15*a^2+45*a^3+a^4", domain, vars),
                parse("28+56*a+45*a^2+23*a^3+a^4", domain, vars),
                parse("19+a^6", domain, vars)
        };

        lEvaluation evaluation = new lEvaluation(2, new long[]{56}, domain, base.ordering);
        int degree = base.degree(1) + 1;
        for (int i = 0; i < AbstractTest.its(10, 100); i++) {
            long start = System.nanoTime();
            MultivariatePolynomialZp64[] lifted = uFactors.clone();
            HenselLifting.bivariateLiftNoLCCorrection0(base, lifted, evaluation, degree);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(base, evaluation.modImage(base.createOne().multiply(lifted), 1, degree));
        }
    }

    @Test
    public void testBivariateLifting2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(62653);
        String[] vars = {"a", "b"};

        MultivariatePolynomialZp64[] factors = {
                MultivariatePolynomialZp64.parse("17096+6578*a*b^2+54905*a^3", domain, vars),
                MultivariatePolynomialZp64.parse("43370+32368*a^2*b^2+45712*a^2*b^4+52302*a^4+23776*a^4*b^2", domain, vars)
        };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        lEvaluation evaluation = new lEvaluation(2, new long[]{0}, domain, base.ordering);

        MultivariatePolynomialZp64[] uFactors = evaluation.evaluateFrom(factors, 1);

        int degree = base.degree(1) + 1;
        HenselLifting.bivariateLift0(base, uFactors, null, evaluation, degree);
        Assert.assertTrue(evaluation.modImage(base.clone().subtract(base.createOne().multiply(uFactors)), 1, degree).isZero());
    }

    @Test
    public void testBivariateLiftingRandom1() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource baseSource = new MultivariateFactorizationTest.lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 4);
        baseSource.minModulusBits = 15;
        baseSource.maxModulusBits = 30;
        FactorizationInput.SampleDecompositionSource<MultivariatePolynomialZp64> source
                = orderVarsByDegree(
                filterMonomialContent(
                        filterNonSquareFree(baseSource)));
        testHenselLift(source, AbstractTest.its(500, 1000), lEvaluationLoop::new,
                ((base, factors, factorsLC, evaluation, degreeBounds, from) ->
                        HenselLifting.bivariateLift0(base, factors, factorsLC, evaluation, degreeBounds[1])), true, 1);
    }

    @Test
    public void testBivariateLiftingRandom2() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource baseSource = new MultivariateFactorizationTest.lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 4);
        baseSource.minModulusBits = 15;
        baseSource.maxModulusBits = 30;
        FactorizationInput.SampleDecompositionSource<MultivariatePolynomialZp64> source
                = orderVarsByDegree(
                filterMonomialContent(
                        filterNonSquareFree(baseSource)));
        testHenselLift(source, AbstractTest.its(500, 1000), lEvaluationLoop::new,
                ((base, factors, factorsLC, evaluation, degreeBounds, from) ->
                        HenselLifting.bivariateLiftNoLCCorrection0(base, factors, evaluation, degreeBounds[1])), false, 1);
    }


    @Test
    public void test4() throws Exception {
        PrivateRandom.getRandom().setSeed(50);
        IntegersZp64 domain = new IntegersZp64(592346501);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("8864159 + 332216825*a + 307171438*a^2 + 574396609*a^3 + b", domain, vars),
                b = MultivariatePolynomialZp64.parse("364341910 + 56968290*a + 134477777*a^2 + 264733241*b + 223672725*a*b + 365910146*a^2*b + 448183856*b^2 + 56041492*a*b^2 + 1386*a^2*b^2", domain, vars),
                base = a.clone().multiply(b);

        lEvaluation evaluation = new lEvaluation(base.nVariables, new long[]{0}, domain, base.ordering);

        MultivariatePolynomialZp64[] uFactors = {
                evaluation.evaluateFrom(a, 1),
                evaluation.evaluateFrom(b, 1),
        };

        uFactors[1].multiplyByLC(uFactors[0]);
        uFactors[0].monic();

        System.out.println(uFactors[0]);
        System.out.println(uFactors[1]);

        MultivariatePolynomialZp64[] factorsLC = {
                null,
                base.lc(0)
        };

        HenselLifting.multivariateLift0(base,
                uFactors,
                factorsLC,
                evaluation,
                base.degrees());
        Assert.assertEquals(base, multiply(uFactors));
    }

    @Test
    public void test5() throws Exception {
        IntegersZp64 domain = new IntegersZp64(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        parse("a^15*b*c^2 + 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                        parse("a^5*b^6*c^3 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                        //parse("a^3*b*c + a*b^2 + 1", ring, vars)
                },
                base = multiply(factors);

        assert MultivariateGCD.PolynomialGCD(factors).isConstant();

        lEvaluation evaluation = new lEvaluation(base.nVariables, new long[]{1, 2}, domain, base.ordering);

        MultivariatePolynomialZp64[] uFactors = evaluation.evaluateFrom(factors, 1);
        HenselLifting.multivariateLiftAutomaticLC(base, uFactors, evaluation);
        Assert.assertEquals(base, multiply(uFactors));

        MultivariatePolynomialZp64[] biFactors = evaluation.evaluateFrom(factors, 2);
        HenselLifting.multivariateLiftAutomaticLC(base, biFactors, evaluation, 2);
        Assert.assertEquals(base, multiply(biFactors));
    }

    @Test
    public void test5a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        parse("a^15*b*c^2 + 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                        parse("a^5*b^6*c^3 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                        //parse("a^3*b*c + a*b^2 + 1", ring, vars)
                },
                base = multiply(factors);

        assert MultivariateGCD.PolynomialGCD(factors).isConstant();

        lEvaluation evaluation = new lEvaluation(base.nVariables, new long[]{1, 2}, domain, base.ordering);


        MultivariatePolynomialZp64[] biFactors = evaluation.evaluateFrom(factors, 2);
        MultivariatePolynomialZp64 lc = base.lc(0);

//         imposing leading coefficients
        MultivariatePolynomialZp64 lcCorrection = evaluation.evaluateFrom(lc, 2);
//        assert lcCorrection.isConstant();

        System.out.println(Arrays.stream(biFactors).map(p -> p.lc(0)).reduce(base.createOne(), (a, b) -> a.clone().multiply(b)));
        System.out.println(evaluation.evaluateFrom(base.lc(0), 2));

        for (MultivariatePolynomialZp64 factor : biFactors) {
            MultivariatePolynomialZp64 r = MultivariateDivision.divideExact(lcCorrection, factor.lc(0));
            factor.multiply(r);
//            assert factor.lt().exponents[0] == factor.degree(0);
//            factor.setLC(0, evaluation.evaluateFrom(factor.lc(0), 1));
//            factor.monicWithLC(lcCorrection.lcAsPoly());
        }

        MultivariatePolynomialZp64 tmp = base.clone().multiply(PolynomialMethods.polyPow(lc, biFactors.length - 1, true));

        HenselLifting.multivariateLift0(tmp, biFactors, ArraysUtil.arrayOf(lc, biFactors.length), evaluation, base.degrees(), 2);

        for (MultivariatePolynomialZp64 factor : biFactors)
            factor.set(HenselLifting.primitivePart(factor));

        System.out.println(Arrays.toString(factors));
        System.out.println(Arrays.toString(biFactors));
        Assert.assertEquals(base, multiply(biFactors));
    }

    @Ignore
    @Test
    public void test6() throws Exception {
        // todo: discover
        IntegersZp64 domain = new IntegersZp64(33554467);
        String[] vars = {"a", "b", "c", "d", "e", "f", "g"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("25078271*b^5*c*d^2*e^2*f^6*g^3+22985334*a*b*c^2*d*e^7*f^3*g^7+19249719*a*b^7*d^5*e^6*g^2+6865506*a^2*b^5*c^3*d^6*e^6*f^3*g^5+20943085*a^2*b^5*c^8*d^3*e^3*f^7+733087*a^3*c^3*d^4*f^4*g^2+24327652*a^3*b^2*c^2*d^2*e^2*f^3*g^5+2508535*a^3*b^3*c*d^3*e^5*f^2*g^2+9991244*a^3*b^4*c^5*e^5*f^5*g^3+22044750*a^3*b^7*c*d^8*e*f^6+8526153*a^4*c^8*d*e^8*f^4*g^6+15162335*a^4*b^8*c^3*d^4*f^4*g^6+21943911*a^5*b*c^3*d^2*e^5*g^2+7268253*a^5*b^8*c^4*d^4*e*f*g^5+11265450*a^6*b^3*c^5*d^5*e+1307471*a^6*b^5*c^4*d^3*e*f^7+27352310*a^7*b^2*c^2*d^6*e^3*f^3*g+18596343*a^8*b^3*e^4*f^2*g+477464*a^8*b^4*c^3*d*e^3*f^5*g^3+20723946*a^8*b^4*c^8*d^3*e^2*g^3", domain, vars),
                        MultivariatePolynomialZp64.parse("24999489*c^6*d^5*g^7+31605719*b^5*c^5*d^4*e^3*g^6+33475465*b^8*c^5*d^6*f^7+21150942*a*c^4*d^4*e^3*f^3*g+30544835*a*b^3*d^7*e*f^8*g^5+8725705*a*b^8*c^6*d^5*e^4*f*g+3830207*a^2*d^2*e*f^7*g^8+31725230*a^2*b^8*c*e*f^8*g^2+5924640*a^3*c*d^4*e^8+14191319*a^3*b*c^3*e^7*f^3*g^8+5482302*a^3*b^2*c^2*d^5*f^8*g^2+350050*a^4*b^3*c^6*d^6*e^7*f^6+246147*a^4*b^4*c^3*d^7*e^5*g^8+27052604*a^5*c^4*d^4*f+2073523*a^5*b^4*c^4*d^7*e^4*f^2*g^5+21322895*a^5*b^5*c*d^3*e^5*f^5*g^4+19375356*a^5*b^6*c^7*e^3*f^2*g^6+15676776*a^6*b^7*c^3*d^8*e^3*f^6*g^6+9971731*a^7*b^4*c^3*d*e^5*f*g^5+16734963*a^8*b^8*c^7*d^4*e*g^7", domain, vars),
                        MultivariatePolynomialZp64.parse("21113415*b^5*c^6*d^4*e^4*f^2*g^3+20864231*b^8*c*d^6*e^5*f^8*g^5+33448448*a*b^3*c^6*d*e^4*f^7*g^2+31133965*a*b^4*c^2*d^2*e^7*f^6*g^2+27612593*a*b^5*d^5*e^2*f^7*g^4+17128197*a*b^7*c^3*d^6*e^2+4469686*a^2*b^5*c^4*d^8*e^4*f^4*g^7+1374035*a^3*c^8*e^7*f*g^5+10414621*a^3*b^6*c^5*d^7*e^7*f^6*g^6+10872067*a^3*b^8*c^3*d*e^4*f^8*g^4+6381772*a^4*b^2*c^6*d^6*e^6*f^3*g^3+26978581*a^4*b^5*d^6*e^5*f^7+30602413*a^4*b^8*c^8*d^4*e^5*f^3*g^3+13372094*a^5*b^3*c^3*d^7*e^5*f^8*g^3+25263857*a^5*b^5*c*d^7*e^6*g^5+4204332*a^6*c^2*d^2*e*f^6*g^2+13228578*a^6*b^2*c^5*d^7*e^6*f^8*g^6+17934510*a^6*b^8*c^4*d^5*e^3*f^4+17371834*a^7*b^4*c^2*d^8*e^4*f^2*g+8745908*a^8*b*c^4*d^7*e^5*f*g^6", domain, vars)
                };
        for (int i = 0; i < factors.length; i++)
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 2, 3);

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        System.out.println(base.size());
        System.out.println();
//        if(true) return;
        long[] values = {32161817, 2822446, 31477240, 17500389, 4067697, 32147770};
        lEvaluation evaluation = new lEvaluation(base.nVariables, values, domain, base.ordering);

        MultivariatePolynomialZp64[] lcFactors = Arrays.stream(factors).map(f -> f.lc(0)).toArray(MultivariatePolynomialZp64[]::new);

        for (int i = 4; i <= base.nVariables; i++) {

            MultivariatePolynomialZp64 iBase = evaluation.evaluateFrom(base, i);
            MultivariatePolynomialZp64[] iLcFactors = lcFactors.clone();
            for (int j = 0; j < iLcFactors.length; j++)
                iLcFactors[j] = evaluation.evaluateFrom(iLcFactors[j], i);

            for (int j = 0; j < 1111; j++) {

                MultivariatePolynomialZp64[] biImages = Arrays.stream(factors).map(f -> evaluation.evaluateFrom(f, 2)).toArray(MultivariatePolynomialZp64[]::new);
                long start = System.nanoTime();
//                HenselLifting.multivariateLift0_old(iBase, biImages, iLcFactors, evaluation, base.degrees(), 2);
//                System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
//                assert multiply(biImages).equals(iBase);

                biImages = Arrays.stream(factors).map(f -> evaluation.evaluateFrom(f, 2)).toArray(MultivariatePolynomialZp64[]::new);
                start = System.nanoTime();
                HenselLifting.multivariateLift0(iBase, biImages, iLcFactors, evaluation, base.degrees(), 2);
                System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
                assert multiply(biImages).equals(iBase);
                System.out.println();
            }

        }

//        assert MultivariateSquareFreeFactorization.isSquareFree(base);
//        for (int i = 0; i < its(20, 20); i++) {
//            long start = System.nanoTime();
//            FactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
//            Assert.assertEquals(3, decomposition.size());
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
//        }
    }


    /* ==================================== Test data =============================================== */

    private interface Lifting<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        void lift(Poly base, Poly[] factors, Poly[] factorsLC, IEvaluation<Term, Poly> evaluation, int[] degreeBounds, int from);
    }

    final class BivariateLift<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
            implements Lifting<Term, Poly> {
        @Override
        public void lift(Poly base, Poly[] factors, Poly[] factorsLC, IEvaluation<Term, Poly> evaluation, int[] degreeBounds, int from) {
            HenselLifting.bivariateLift0(base, factors, factorsLC, evaluation, degreeBounds[1]);
        }
    }

    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void testHenselLift(FactorizationInput.SampleDecompositionSource<Poly> source, int nIterations,
                        Function<Poly, IEvaluationLoop<Term, Poly>> evalFactory,
                        Lifting<Term, Poly> algorithm, boolean correctLC, int from) {
        System.out.println("Testing Hensel lifting ");
        System.out.println("Input source: " + source);

        DescriptiveStatistics timing = new DescriptiveStatistics();

        int prevProgress = -1, currProgress;
        main:
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                timing.clear();

            if ((currProgress = (int) (100.0 * n / nIterations)) != prevProgress) {
                prevProgress = currProgress;
                System.out.print(">");
                System.out.flush();
            }
            FactorizationInput.SampleDecomposition<Poly> sample = source.next();
            if (!allCoprime(sample.factors)) {
                --n;
                continue;
            }
            if (Arrays.stream(sample.factors).anyMatch(p -> p.degree(0) == 0)) {
                --n;
                continue;
            }
            try {

                Poly factory = sample.poly;
                IEvaluationLoop<Term, Poly> evaluations = evalFactory.apply(factory);

                for (int nAttempt = 0; nAttempt < 64; nAttempt++) {
                    IEvaluation<Term, Poly> evaluation = evaluations.next();
                    Poly[] uFactors = Arrays.stream(sample.factors).map(p -> evaluation.evaluateFrom(p, from)).toArray(factory::createArray);
                    if (!allCoprime(uFactors))
                        continue;

                    if (!IntStream.range(0, uFactors.length).allMatch(i -> sample.factors[i].degree(0) == uFactors[i].degree(0)))
                        continue;

                    Poly[] factorsLC = Arrays.stream(sample.factors).map(p -> p.lc(0)).toArray(factory::createArray);

                    long start = System.nanoTime();
                    algorithm.lift(sample.poly, uFactors, correctLC ? factorsLC : null, evaluation, sample.poly.degrees(), from);
                    timing.addValue(System.nanoTime() - start);


                    if (correctLC || Arrays.stream(factorsLC).allMatch(Poly::isConstant))
                        Assert.assertArrayEquals(sample.factors, uFactors);
                    else
                        Assert.assertEquals(sample.poly, evaluation.modImage(multiply(uFactors), sample.poly.degrees()));
                    continue main;
                }
                throw new RuntimeException();
            } catch (Throwable throwable) {
                System.out.println("\n============ Error ============");
                System.out.println("Domain: " + sample.poly.coefficientRingToString());
                System.out.println("Polynomial: " + sample.poly);
                System.out.println("Expected factorization: " + Arrays.toString(sample.factors));
                throw throwable;
            }
        }

        System.out.println(source.statisticsToString());

        System.out.println("\n============ Timings ============");
        System.out.println("Stats: " + TimeUnits.statisticsNanotime(timing));
    }

    static <Poly extends IPolynomial<Poly>> Poly multiply(Poly... p) {
        return p[0].createOne().multiply(p);
    }

    private static <Poly extends AMultivariatePolynomial> boolean allCoprime(Poly... factors) {
        return StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(Combinatorics.combinations(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());
    }
}