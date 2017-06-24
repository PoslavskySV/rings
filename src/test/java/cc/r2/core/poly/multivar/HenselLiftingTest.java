package cc.r2.core.poly.multivar;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.multivar.HenselLifting.AllProductsCache;
import cc.r2.core.poly.multivar.HenselLifting.UMultiDiophantineSolver;
import cc.r2.core.poly.multivar.HenselLifting.lEvaluation;
import cc.r2.core.poly.univar.RandomPolynomials;
import cc.r2.core.poly.univar.UnivariateGCD;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;
import java.util.BitSet;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.stream.StreamSupport;

import static cc.r2.core.poly.multivar.HenselLifting.liftWang;
import static cc.r2.core.poly.multivar.HenselLifting.modImage;
import static cc.r2.core.poly.multivar.lMultivariatePolynomialZp.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class HenselLiftingTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*b + 2 + b^2*a - b^4", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^2*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("1 + x - x^2 + x^2*y", domain, vars),
                b = parse("2 + x + x^2 + 2*y*x^2 + x^2*y^5 + y", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = a.evaluate(1, 0),
                bF = b.evaluate(1, 0);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 - 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test5() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"a", "b", "c"};
        lMultivariatePolynomialZp
                a = parse("a^15 + a^15*b^2*c^3 + a^15*c + 7*a^15- 2*a*b^4 - 3*c*b + 2 + b^2*a - c*b^4 + c^3", domain, vars),
                b = parse("a^5  + a^5*b*c - a^5*b + 2*a^5 + a*b^2 - 3*b^2 + b + 2 - a^3*b^6*c^3", domain, vars),
                base = a.clone().multiply(b);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();


        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Test
    public void test6() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp
                b = parse("1227874+3587706*b+5373508*a+7197578*a^2+a^3", domain, vars),
                a = parse("9540707+24*a", domain, vars),
                base = parse("7717493+597721*b+6517458*b^2+361611*a+9241048*a*b+9607947*a*b^2+3165308*a^2+9338813*a^3+24*a^4", domain, vars);

        assert MultivariateGCD.PolynomialGCD(a, b).isConstant();

        lMultivariatePolynomialZp
                aF = modImage(a.clone(), 1),
                bF = modImage(b.clone(), 1);

        HenselLifting.liftPair(base, aF, bF);
        Assert.assertEquals(base, aF.clone().multiply(bF));
    }

    @Ignore
    @Test
    public void testHenselLiftingRandom() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = 1000;
        for (int i = 0; i < nIterations; i++) {
            System.out.println(i);
            long modulus = getModulusRandom(rndd.nextInt(5, 30));
            lIntegersModulo domain = new lIntegersModulo(modulus);

            int nVars = rndd.nextInt(2, 3);

            lMultivariatePolynomialZp[] polys = new lMultivariatePolynomialZp[2];
            for (int j = 0; j < polys.length; j++) {
                polys[j] = MultivariatePolynomial.asLongPolyZp(
                        RandomMultivariatePolynomial.randomPolynomial(nVars,
                                rndd.nextInt(5, 15), rndd.nextInt(5, 25), domain.asDomain(), DegreeVector.LEX, rnd));
            }

            lMultivariatePolynomialZp a = polys[0], b = polys[1];
//            //make monic
//            a = a.add(a.createUnivariateMonomial(0, a.degree(0) + 2));
//            b = b.add(b.createUnivariateMonomial(0, b.degree(0) + 2));

            System.out.println("gcd");
            if (!MultivariateGCD.PolynomialGCD(a, b).isOne()) {
                System.out.println("bad inp");
                --i;
                continue;
            }
            System.out.println("gcd...done");

            int nEvaluations = 10;
            evals:
            for (int n = 0; n < nEvaluations; n++) {
                long[] substitutions = new long[nVars - 1];
                for (int j = 0; j < substitutions.length; j++) {
                    do {
                        substitutions[j] = domain.randomElement(rnd);
                    } while (substitutions[j] == 0);
                }


                lEvaluation evaluation = new lEvaluation(nVars, substitutions, domain, DegreeVector.LEX);
                lMultivariatePolynomialZp
                        ua = evaluation.evaluateFrom(a, 1),
                        ub = evaluation.evaluateFrom(b, 1);
                if (!MultivariateGCD.PolynomialGCD(ua, ub).isOne()) {
                    System.out.println("bad p");
                    --n;
                    continue;
                }

                if (!ua.getSkeleton().equals(a.getSkeleton(0)) || !ub.getSkeleton().equals(b.getSkeleton(0))) {
                    System.out.println("bad p");
                    --n;
                    continue;
                }

                try {
                    //System.out.println(a);
                    //System.out.println(b);
                    liftWang(a.clone().multiply(b), ua, ub, a.lc(0), b.lc(0), evaluation);
                    Assert.assertEquals(a, ua);
                    Assert.assertEquals(b, ub);
                } catch (Throwable thr) {
                    System.out.println(domain);
                    System.out.println(a);
                    System.out.println(b);
                    System.out.println(Arrays.toString(evaluation.values));
                    throw thr;
                }
            }
        }
    }

    @Test
    public void testHenselLifting1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(9607987);
        String[] vars = {"a", "b", "c"};

        lMultivariatePolynomialZp
                a = parse("1 + c*b + b*a*c^5 + b*c*2*a^2", domain, vars),
                b = parse("1227874+3587706*b+5373508*a+7197578*c^2*a^2+a^3", domain, vars),
                base = a.clone().multiply(b);

        long[] subs = {11, 12};
        lEvaluation evaluation = new lEvaluation(vars.length, subs, domain, DegreeVector.LEX);

        lMultivariatePolynomialZp
                ua = evaluation.evaluateFrom(a, 1),
                ub = evaluation.evaluateFrom(b, 1);

        HenselLifting.liftWang(base, ua, ub, a.lc(0), b.lc(0), evaluation);

        Assert.assertEquals(a, ua);
        Assert.assertEquals(b, ub);
    }

    @Test
    public void testProducts1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        for (int n = 0; n < 10; n++) {

            for (int nFactors : new int[]{4, 5, 6, 7}) {
                lUnivariatePolynomialZp[] factors = new lUnivariatePolynomialZp[nFactors];
                for (int i = 0; i < factors.length; i++)
                    factors[i] = RandomPolynomials.randomMonicPoly(5, domain.modulus, getRandom());

                AllProductsCache<lUnivariatePolynomialZp> cache = new AllProductsCache<>(factors);
                for (int i = 1; i < factors.length; i++) {
                    IntCombinationsGenerator combinations = new IntCombinationsGenerator(factors.length, i);
                    for (int[] combination : combinations) {
                        lUnivariatePolynomialZp expected = factors[0].createOne();
                        BitSet key = new BitSet(factors.length);
                        for (int k : combination) {
                            key.set(k);
                            expected = expected.multiply(factors[k]);
                        }

                        Assert.assertEquals(expected, cache.multiply(key));
                    }
                }
            }

        }
    }

    @Test
    public void testDiophantineSolver1() throws Exception {
        RandomGenerator random = getRandom();
        lIntegersModulo domain = new lIntegersModulo(17);
        for (int n = 0; n < 10; n++) {

            for (int nFactors : new int[]{4, 5, 6, 7}) {
                lUnivariatePolynomialZp[] factors = new lUnivariatePolynomialZp[nFactors];

                out:
                for (int i = 0; i < factors.length; i++) {
                    factors[i] = RandomPolynomials.randomMonicPoly(5, domain.modulus, random);
                    for (int j = 0; j < i; j++) {
                        if (!UnivariateGCD.PolynomialGCD(factors[j], factors[i]).isConstant()) {
                            --i;
                            continue out;
                        }
                    }
                }

                assert UnivariateGCD.PolynomialGCD(factors).isConstant();

                AllProductsCache<lUnivariatePolynomialZp> cache = new AllProductsCache<>(factors);
                UMultiDiophantineSolver<lUnivariatePolynomialZp> solver = new UMultiDiophantineSolver<>(cache);
                lUnivariatePolynomialZp rhs = RandomPolynomials.randomMonicPoly(5, domain.modulus, random);
                solver.solve(rhs);

                lUnivariatePolynomialZp actual = rhs.createZero();
                for (int i = 0; i < factors.length; i++)
                    actual = actual.add(cache.except(i).clone().multiply(solver.solution[i]));

                Assert.assertEquals(rhs, actual);
            }
        }
    }


    @Test
    public void testDiophantineSolver2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        lUnivariatePolynomialZp[] factors = {
                lUnivariatePolynomialZ.create(2, 11, 30, 1).modulus(domain),
                lUnivariatePolynomialZ.create(54, 0, 24, 1).modulus(domain),
                lUnivariatePolynomialZ.create(31, 0, 60, 1).modulus(domain)
        };


        assert UnivariateGCD.PolynomialGCD(factors).isConstant();

        AllProductsCache<lUnivariatePolynomialZp> cache = new AllProductsCache<>(factors);
        UMultiDiophantineSolver<lUnivariatePolynomialZp> solver = new UMultiDiophantineSolver<>(cache);
        lUnivariatePolynomialZp rhs = lUnivariatePolynomialZ.one().modulus(domain);
        solver.solve(rhs);

        lUnivariatePolynomialZp actual = rhs.createZero();
        for (int i = 0; i < factors.length; i++)
            actual = actual.add(cache.except(i).clone().multiply(solver.solution[i]));

        Assert.assertEquals(rhs, actual);
    }


    @Test
    public void testMultiLift1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^5*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                c = parse("x^5 - y*x^2 - 3*y^2 + y + 2 - y^3*x^5", domain, vars),
                d = parse("x^5 - y*x^2 - 3*y + x + 2 - y*x", domain, vars);

        lMultivariatePolynomialZp[] factors = {a, b, c, d};
        assert StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(new IntCombinationsGenerator(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());

        lMultivariatePolynomialZp[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(lMultivariatePolynomialZp[]::new);
        lMultivariatePolynomialZp base = a.createOne().multiply(factors);

        long[] vals = {31};
        lEvaluation evaluation = new lEvaluation(a.nVariables, vals, a.domain, a.ordering);

        lMultivariatePolynomialZp[] uFactors = evaluation.evaluateFrom(factors, 1);
        liftWang(base, uFactors, factorsLC, evaluation);

        Assert.assertArrayEquals(factors, uFactors);
    }

    @Test
    public void testMultiLift2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y", "z"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^2*y - x^2", domain, vars),
                b = parse("x^5 + y*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                c = parse("x^5 - y*x^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                d = parse("x^5 - y*x^2 - 3*y + x + z - y*x", domain, vars);

        lMultivariatePolynomialZp[] factors = {c, d};
        assert StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(new IntCombinationsGenerator(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());

        lMultivariatePolynomialZp base = a.createOne().multiply(factors);

        long[] vals = {31, 32};
        lEvaluation evaluation = new lEvaluation(a.nVariables, vals, a.domain, a.ordering);

        lMultivariatePolynomialZp[] uFactors = evaluation.evaluateFrom(factors, 1);

        liftWang(base, uFactors, null, evaluation);

        Assert.assertArrayEquals(factors, uFactors);
    }

    @Test
    public void testMultiLift3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(SmallPrimes.nextPrime(66));
        String[] vars = {"x", "y", "z"};
        lMultivariatePolynomialZp
                a = parse("x^5 - 2*x*y^4 - 3*x + 2 + x^5*y - x^2 - z", domain, vars),
                b = parse("x^5 + z*y^2 - 3*y^2 + y + 2 - y^3*x^2", domain, vars),
                c = parse("x^5 - y*x^2 - 3*y^2 + y + 2 - z - y^3*x^5 + z^3*x^5", domain, vars),
                d = parse("x^5 - y*x^2 - 3*y + x + 2 - y*x + y*z*x^5", domain, vars);

        lMultivariatePolynomialZp[] factors = {a, b, c, d};
        assert StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(new IntCombinationsGenerator(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());


        lMultivariatePolynomialZp[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(lMultivariatePolynomialZp[]::new);
        lMultivariatePolynomialZp base = a.createOne().multiply(factors);

        assert MultivariateGCD.PolynomialGCD(a, b, c, d).isConstant();

        long[] vals = {31, 32};
        lEvaluation evaluation = new lEvaluation(a.nVariables, vals, a.domain, a.ordering);


        lMultivariatePolynomialZp[] uFactors = evaluation.evaluateFrom(factors, 1);
        liftWang(base, uFactors, factorsLC, evaluation);

        Assert.assertArrayEquals(factors, uFactors);
    }

    @Test
    public void testMultiLift4() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(17);
        String[] vars = {"a", "b", "c", "d", "e"};
        lMultivariatePolynomialZp
                a = parse("6*d+14*c^2*e^2+14*b^2*c*d^2*e+3*a*c*d^2*e+5*a*b*c^2*d*e^2+10*a^2*d^2*e^2+10*a^2*b*c*d*e", domain, vars),
                b = parse("2*b^2*c*d^3*e^3+3*a^2*b^4*c^2*d^3*e^6+9*a^2*b^6*c*d^2*e^3+4*a^3*b^5*c^2*e^3+2*a^3*b^5*c^6*d+11*a^5*b*c^6*d^3*e^5+11*a^6*b^2*c^6*d^4*e", domain, vars),
                c = parse("16*c^3*d^2*e^3+2*b^3*c^3+4*a^2*b^2*d*e^3+2*a^2*b^3*c^2*d^2*e^2", domain, vars);

        lMultivariatePolynomialZp[] factors = {a, b, c};
        assert StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(new IntCombinationsGenerator(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());


        lMultivariatePolynomialZp[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(lMultivariatePolynomialZp[]::new);
        lMultivariatePolynomialZp base = a.createOne().multiply(factors);

        long[] vals = {31, 32, 33, 32};
        lEvaluation evaluation = new lEvaluation(a.nVariables, vals, a.domain, a.ordering);


        lMultivariatePolynomialZp[] uFactors = evaluation.evaluateFrom(factors, 1);
        liftWang(base, uFactors, factorsLC, evaluation);

        Assert.assertArrayEquals(factors, uFactors);
    }

    @Test
    public void testMultiLiftRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(100, 500);

        int prevPercent = -1, currPercent;
        main:
        for (int n = 0; n < nIterations; n++) {

            lIntegersModulo domain = new lIntegersModulo(getModulusRandom(rndd.nextInt(4, 30)));

            int nVariables = 5;
            lMultivariatePolynomialZp[] factors = new lMultivariatePolynomialZp[rndd.nextInt(2, 5)];
            out:
            for (int i = 0; i < factors.length; i++) {
                factors[i] = RandomMultivariatePolynomial.randomPolynomial(nVariables,
                        rndd.nextInt(2, 5),
                        rndd.nextInt(3, 7),
                        domain, DegreeVector.LEX, rnd);
                for (int j = 0; j < i; j++) {
                    if (!MultivariateGCD.PolynomialGCD(factors[j], factors[i]).isConstant()) {
                        --i;
                        continue out;
                    }
                }
            }

            lMultivariatePolynomialZp[] factorsLC = Arrays.stream(factors).map(f -> f.lc(0)).toArray(lMultivariatePolynomialZp[]::new);
            lMultivariatePolynomialZp base = factors[0].createOne().multiply(factors);

            long[] substitutions = new long[nVariables - 1];
            for (int j = 0; j < substitutions.length; j++)
                substitutions[j] = domain.randomNonZeroElement(rnd);


            lEvaluation evaluation = new lEvaluation(nVariables, substitutions, domain, base.ordering);
            lMultivariatePolynomialZp[] uFactors = evaluation.evaluateFrom(factors, 1);

            if (!allCoprime(uFactors)) {
                --n;
                continue;
            }

            if (Arrays.stream(uFactors).anyMatch(lMultivariatePolynomialZp::isConstant)) {
                --n;
                continue;
            }

            for (int i = 0; i < factors.length; i++)
                if (uFactors[i].degree() != factors[i].degree(0)) {
                    --n;
                    continue main;
                }


            //boolean b = execute(() -> liftWang(base, uFactors, factorsLC, evaluation), 100_000);
            liftWang(base, uFactors, factorsLC, evaluation);
            try {
                Assert.assertArrayEquals(factors, uFactors);
            } catch (Throwable r) {
                System.out.println(domain);
                System.out.println(Arrays.toString(substitutions));
                System.out.println(base.size());
                for (lMultivariatePolynomialZp factor : factors)
                    System.out.println(factor);

                System.out.println(Arrays.equals(factors, uFactors));
                throw r;
            }
            if ((currPercent = (int) (100. * n / nIterations)) != prevPercent) {
                System.out.print(">");
                System.out.flush();
                prevPercent = currPercent;
            }
        }
    }

//    static boolean execute(Runnable runnable, long milliseconds) throws InterruptedException {
//        final AtomicBoolean mutex = new AtomicBoolean(false);
//        CountDownLatch latch = new CountDownLatch(1);
//        Thread base = new Thread(() -> {
//            try {
//                Thread.sleep(milliseconds);
//            } catch (InterruptedException e) {}
//            latch.countDown();
//        });
//
//        Thread exec = new Thread(() -> {runnable.run(); mutex.set(true); latch.countDown();});
//        base.start();
//        exec.start();
//
//        latch.await();
//
//        base.stop();
//        exec.stop();
//        base.join();
//        exec.join();
//
//        return mutex.get();
//    }

    private static <Poly extends AMultivariatePolynomial> boolean allCoprime(Poly... factors) {
        return StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(new IntCombinationsGenerator(factors.length, 2), Spliterator.ORDERED), false)
                .allMatch(arr -> MultivariateGCD.PolynomialGCD(factors[arr[0]], factors[arr[1]]).isOne());
    }
}