package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly2.HenselLifting.LiftableQuintet;
import cc.r2.core.test.Benchmark;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import static cc.r2.core.poly2.Factorization.factor;
import static cc.r2.core.poly2.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.poly2.SquareFreeFactorization.SquareFreePart;
import static cc.r2.core.poly2.SquareFreeFactorization.isSquareFree;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class HenselLiftingTest extends AbstractPolynomialTest {

    @SuppressWarnings("unchecked")
    private static <PolyZ extends IMutablePolynomialZ<PolyZ>, PolyZp extends IMutablePolynomialZp<PolyZp>>
    PolyZp modulus(PolyZ poly, long modulus) {
        if (poly instanceof bMutablePolynomialZ)
            return (PolyZp) ((bMutablePolynomialZ) poly).modulus(modulus);
        else
            return (PolyZp) ((MutablePolynomialZ) poly).modulus(modulus);
    }

    @SuppressWarnings("unchecked")
    private static <PolyZ extends IMutablePolynomialZ<PolyZ>, PolyZp extends IMutablePolynomialZp<PolyZp>>
    PolyZp modulus(PolyZ poly, BigInteger modulus) {
        if (poly instanceof bMutablePolynomialZ)
            return (PolyZp) ((bMutablePolynomialZ) poly).modulus(modulus);
        else
            return (PolyZp) ((MutablePolynomialZ) poly).modulus(modulus.longValueExact());
    }

    @SuppressWarnings("unchecked")
    private static <PolyZ extends IMutablePolynomialZ<PolyZ>, PolyZp extends IMutablePolynomialZp<PolyZp>>
    HenselLifting.QuadraticLiftAbstract<PolyZp> createQuadraticLift(long modulus, PolyZ poly, PolyZp aFactor, PolyZp bFactor) {
        if (poly instanceof bMutablePolynomialZ)
            return (HenselLifting.QuadraticLiftAbstract<PolyZp>) HenselLifting.createQuadraticLift(BigInteger.valueOf(modulus),
                    (bMutablePolynomialZ) poly, (bMutablePolynomialMod) aFactor, (bMutablePolynomialMod) bFactor);
        else
            return (HenselLifting.QuadraticLiftAbstract<PolyZp>) HenselLifting.createQuadraticLift(modulus,
                    (MutablePolynomialZ) poly, (MutablePolynomialMod) aFactor, (MutablePolynomialMod) bFactor);
    }

    @SuppressWarnings("unchecked")
    private static <PolyZ extends IMutablePolynomialZ<PolyZ>, PolyZp extends IMutablePolynomialZp<PolyZp>>
    LiftableQuintet<PolyZp> createLinearLift(long modulus, PolyZ poly, PolyZp aFactor, PolyZp bFactor) {
        if (poly instanceof bMutablePolynomialZ)
            return (LiftableQuintet<PolyZp>) HenselLifting.createLinearLift(modulus,
                    (bMutablePolynomialZ) poly, ((bMutablePolynomialMod) aFactor).toLong(), ((bMutablePolynomialMod) bFactor).toLong());
        else
            return (LiftableQuintet<PolyZp>) HenselLifting.createLinearLift(modulus,
                    (MutablePolynomialZ) poly, (MutablePolynomialMod) aFactor, (MutablePolynomialMod) bFactor);
    }

    <PolyZ extends IMutablePolynomialZ<PolyZ>, PolyZp extends IMutablePolynomialZp<PolyZp>>
    void testHenselStepRandomModulus(
            PolyZ poly,
            PolyZ aFactor,
            PolyZ bFactor,
            boolean quadratic,
            int nIterations,
            int nPrimes) {

        out:
        for (long modulus : getModulusArray(nPrimes, 0, 5, 0)) {
            PolyZp aMod = modulus(aFactor, modulus);
            PolyZp bMod = modulus(bFactor, modulus);
            if (aMod.degree() != aFactor.degree() || bMod.degree() != bFactor.degree())
                continue;

            if (!PolynomialGCD.PolynomialGCD(aMod, bMod).isConstant())
                continue;
            LiftableQuintet<PolyZp> hensel =
                    quadratic ? createQuadraticLift(modulus, poly, aMod, bMod) : createLinearLift(modulus, poly, aMod, bMod);
            assertHenselLift(hensel);

            for (int i = 0; i < nIterations - 1; i++) {
                try {
                    hensel.lift();
                } catch (ArithmeticException ex) {
                    if (!"long overflow".equals(ex.getMessage()))
                        throw ex;
                    continue out;
                }
                assertHenselLift(hensel);
            }
            hensel.liftLast();
            assertHenselLift(hensel);
        }
    }

    @Test
    public void testHensel1() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(-2, -1, 2, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(-2, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 10, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 10, 20);
    }

    @Test
    public void testHensel2() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(3, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(4, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 10, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 10, 20);
    }

    @Test
    public void testHensel3() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(12, 13112, 1112, 31, 112);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(22, 112311, 12);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);

        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 10, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 10, 20);
    }

    @Test
    public void testHensel4_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = (int) its(100, 1000);
        for (int i = 0; i < nIterations; i++) {
            bMutablePolynomialZ a = RandomPolynomials.randomPoly(rndd.nextInt(2, 5), BigInteger.INT_MAX_VALUE, rnd);
            bMutablePolynomialZ b = RandomPolynomials.randomPoly(rndd.nextInt(2, 5), BigInteger.INT_MAX_VALUE, rnd);
            bMutablePolynomialZ poly = a.clone().multiply(b);
            testHenselStepRandomModulus(poly, a, b, true, 5, 5);
            testHenselStepRandomModulus(poly, a, b, false, 5, 5);
        }
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus, int nIterations, boolean quadratic) {
        base = SquareFreePart(base);

        MutablePolynomialMod baseMod = base.modulus(modulus).toLong();
        if (baseMod.degree() != base.degree())
            return;

        if (!isSquareFree(baseMod))
            return;

        FactorDecomposition<MutablePolynomialMod> modularFactors = factor(baseMod);
        assertFactorization(baseMod, modularFactors);

        HenselLifting.LiftFactory<MutablePolynomialMod> factory = quadratic ? HenselLifting::createQuadraticLift : HenselLifting::createLinearLift;
        BigInteger newModulus = newModulus(modulus, nIterations, quadratic);
        List<bMutablePolynomialMod> tmp = HenselLifting.liftFactorization0(BigInteger.valueOf(modulus), newModulus, nIterations, base, modularFactors.factors, factory);
        bMutablePolynomialMod bm = base.modulus(newModulus);
        assertFactorization(bm, bm.lc(), tmp);
    }

    static BigInteger newModulus(long modulus, int nIterations, boolean quadratic) {
        BigInteger result, bModulus = result = BigInteger.valueOf(modulus);
        for (int i = 0; i < nIterations; i++)
            result = result.multiply(quadratic ? result : bModulus);
        return result;
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus, boolean quadratic) {
        testMultiFactorHenselLifting(base, modulus, HenselLifting.nIterations(BigInteger.valueOf(modulus), base.mignotteBound().shiftLeft(1), quadratic).nIterations, quadratic);
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus) {
        testMultiFactorHenselLifting(base, modulus, true);
        testMultiFactorHenselLifting(base, modulus, false);
    }

    @Test
    public void testHensel5() throws Exception {
        bMutablePolynomialZ base = bMutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 5, 5, true);
        testMultiFactorHenselLifting(base, 5, 5, false);
    }

    @Test
    public void testHensel6() throws Exception {
        bMutablePolynomialZ base = bMutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 11, 5, true);
        testMultiFactorHenselLifting(base, 11, 5, false);
    }

    @Test
    public void testHensel7_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(100, 1000); ++i) {
            bMutablePolynomialZ poly = RandomPolynomials.randomPoly(rndd.nextInt(5, 15), BigInteger.LONG_MAX_VALUE, rnd);
            for (long modulus : getModulusArray((int) its(5, 15), 0, 5, 0)) {
                MutablePolynomialMod polyMod = poly.modulus(modulus).toLong();
                if (polyMod.isConstant() || polyMod.isMonomial())
                    continue;

                while (poly.lc().mod(BigInteger.valueOf(modulus)).isZero())
                    poly.data[poly.degree] = BigInteger.valueOf(rnd.nextInt());

                try {
                    testMultiFactorHenselLifting(poly, modulus);
                } catch (AssertionError err) {
                    System.out.println(poly.toStringForCopy());
                    System.out.println(modulus);
                    throw err;
                }
            }
        }
    }

    @Test
    @Benchmark
    public void testHensel8_linear_vs_quadratic() throws Exception {
        long modulus = SmallPrimes.nextPrime(2132131231);
        RandomGenerator rnd = getRandom();
        bMutablePolynomialZ aFactor, bFactor;
        do {
            aFactor = RandomPolynomials.randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
            bFactor = RandomPolynomials.randomPoly(5, BigInteger.LONG_MAX_VALUE, rnd);
        } while (!PolynomialGCD.PolynomialGCD(aFactor.modulus(modulus), bFactor.modulus(modulus)).isConstant());

        bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);

        int nTrials = 20;
        int maxIterations;
        LiftableQuintet<bMutablePolynomialMod> linearLift, quadraticLift;
        DescriptiveStatistics acc = new DescriptiveStatistics();

        List<double[]> linear = new ArrayList<>(), quadratic = new ArrayList<>();

        //warm-up
        for (int j = 0; j < 64; j++) {
            linearLift = createLinearLift(modulus, poly, aFactor.modulus(modulus), bFactor.modulus(modulus));
            linearLift.lift(j);
        }

        maxIterations = 12 * 12;
        for (int nIts = 1; nIts < maxIterations; nIts++) {
            System.out.println(nIts);
            acc.clear();
            for (int j = 0; j < nTrials; j++) {
                long start = System.nanoTime();
                linearLift = createLinearLift(modulus, poly, aFactor.modulus(modulus), bFactor.modulus(modulus));
                linearLift.lift(nIts);
                acc.addValue(System.nanoTime() - start);
                assertHenselLift(linearLift);
            }
            linear.add(new double[]{nIts, acc.getPercentile(0.5), acc.getStandardDeviation()});
        }


        //warm-up
        for (int j = 0; j < 6; j++) {
            quadraticLift = createQuadraticLift(modulus, poly, aFactor.modulus(modulus), bFactor.modulus(modulus));
            quadraticLift.lift(j);
        }

        maxIterations = 12;
        for (int nIts = 1; nIts < maxIterations; nIts++) {
            System.out.println(nIts);
            acc.clear();
            for (int j = 0; j < nTrials; j++) {
                long start = System.nanoTime();
                quadraticLift = createQuadraticLift(modulus, poly, aFactor.modulus(modulus), bFactor.modulus(modulus));
                quadraticLift.lift(nIts);
                acc.addValue(System.nanoTime() - start);
                assertHenselLift(quadraticLift);
            }
            quadratic.add(new double[]{nIts, acc.getPercentile(0.5), acc.getStandardDeviation()});
        }

        System.out.println("linear = " + toString(linear) + ";");
        System.out.println("quadratic = " + toString(quadratic) + ";");
    }

    private static String toString(List<double[]> d) {
        StringBuilder sb = new StringBuilder().append("{");
        for (int i = 0; ; i++) {
            sb
                    .append("{")
                    .append(Arrays.stream(d.get(i)).mapToObj(Double::toString).collect(Collectors.joining(",")))
                    .append("}");
            if (i == d.size() - 1)
                break;
            sb.append(",");
        }
        sb.append("}");
        return sb.toString();
    }


    @Test
    @Benchmark(runAnyway = true)
    public void testHensel9_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                linear = new DescriptiveStatistics(),
                quadratic = new DescriptiveStatistics(),
                adaptive = new DescriptiveStatistics();

        int nFactors = 6;
        int nIterations = (int) its(30, 150);
        for (int n = 0; n < nIterations; n++) {
            if (nIterations / 10 == n) {
                linear.clear(); quadratic.clear(); adaptive.clear();
            }
            if (n % 10 == 0)
                System.out.println(n);

            bMutablePolynomialZ poly = bMutablePolynomialZ.one();
            for (int i = 0; i < nFactors; i++)
                poly = poly.multiply(RandomPolynomials.randomPoly(rndd.nextInt(5, 15), BigInteger.valueOf(1000), rnd));
            poly = SquareFreeFactorization.SquareFreePart(poly);

            bMutablePolynomialMod polyMod;
            long modulus;
            do {
                modulus = getModulusRandom(rndd.nextInt(5, 28));
                polyMod = poly.modulus(modulus);
            } while (!SquareFreeFactorization.isSquareFree(poly.modulus(modulus)) || polyMod.degree() != poly.degree());

            BigInteger desiredBound = poly.mignotteBound().shiftLeft(1).multiply(poly.lc());
            lFactorDecomposition<MutablePolynomialMod> modularFactors = factor(polyMod.toLong());
            BigInteger bModulus = BigInteger.valueOf(modulus);

            long start;
            try {
                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsQuad = HenselLifting.liftFactorization(bModulus, desiredBound, poly, modularFactors.factors, true);
                quadratic.addValue(System.nanoTime() - start);
                assertFactorization(poly.modulus(factorsQuad.get(0).modulus), poly.lc(), factorsQuad);

                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsLin = HenselLifting.liftFactorization(bModulus, desiredBound, poly, modularFactors.factors, false);
                linear.addValue(System.nanoTime() - start);
                assertFactorization(poly.modulus(factorsLin.get(0).modulus), poly.lc(), factorsLin);


                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsAdp = HenselLifting.liftFactorizationAdaptive(bModulus, desiredBound, poly, modularFactors.factors);
                adaptive.addValue(System.nanoTime() - start);
                assertFactorization(poly.modulus(factorsAdp.get(0).modulus), poly.lc(), factorsAdp);
            } catch (Throwable e) {
                System.out.println("====");
                System.out.println(modulus);
                System.out.println(poly);
                System.out.println(modularFactors);
                System.out.println(poly.toStringForCopy());
                System.out.println("====");
                throw e;
            }
        }

        System.out.println("Quadratic lifting : " + TimeUnits.statisticsNanotime(quadratic, true));
        System.out.println("Linear lifting   : " + TimeUnits.statisticsNanotime(linear, true));
        System.out.println("Adaptive lifting   : " + TimeUnits.statisticsNanotime(adaptive, true));
    }


    private static <T extends IMutablePolynomialZp<T>> void assertHenselLift(LiftableQuintet<T> lift) {
        Assert.assertEquals(lift.toString(), lift.polyMod(), lift.aFactorMod().clone().multiply(lift.bFactorMod()));
        Assert.assertTrue(lift.toString(), lift.aCoFactorMod() == null && lift.bCoFactorMod() == null ||
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod()))
                        .isOne());
    }
}