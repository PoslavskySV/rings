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
            int nPrimes) {

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

            for (int i = 0; i < 10; i++) {
                try {
                    hensel.lift();
                } catch (ArithmeticException ex) {
                    if (!"long overflow".equals(ex.getMessage()))
                        throw ex;
                    break;
                }
                assertHenselLift(hensel);
            }
        }
    }

    @Test
    public void testHensel1() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(-2, -1, 2, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(-2, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 20);
    }

    @Test
    public void testHensel2() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(3, 1);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(4, 1);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 20);
    }

    @Test
    public void testHensel3() throws Exception {
        MutablePolynomialZ aFactor = MutablePolynomialZ.create(12, 13112, 1112, 31, 112);
        MutablePolynomialZ bFactor = MutablePolynomialZ.create(22, 112311, 12);
        MutablePolynomialZ poly = aFactor.clone().multiply(bFactor);

        testHenselStepRandomModulus(poly, aFactor, bFactor, true, 20);
        testHenselStepRandomModulus(poly, aFactor, bFactor, false, 20);
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus, BigInteger newModulus, int nIterations, boolean quadratic) {
        base = SquareFreePart(base);

        MutablePolynomialMod baseMod = base.modulus(modulus).toLong();
        if (baseMod.degree() != base.degree())
            return;

        if (!isSquareFree(baseMod))
            return;

        FactorDecomposition<MutablePolynomialMod> modularFactors = factor(baseMod);
        assertFactorization(baseMod, modularFactors);

        HenselLifting.LiftFactory<MutablePolynomialMod> factory = quadratic ? HenselLifting::createQuadraticLift : HenselLifting::createLinearLift;
        List<bMutablePolynomialMod> tmp = HenselLifting.liftFactorization0(BigInteger.valueOf(modulus), newModulus, nIterations, base, modularFactors.factors, factory);
        bMutablePolynomialMod bm = base.modulus(newModulus);
        assertFactorization(bm, bm.lc(), tmp);
    }

    static BigInteger newModulus(long modulus, int nIterations, boolean quadratic) {
        BigInteger result, bModulus = result = BigInteger.valueOf(modulus);
        if (quadratic)
            for (int i = 0; i < nIterations; i++)
                result = result.multiply(result);
        else
            for (int i = 0; i < nIterations; i++)
                result = result.multiply(bModulus);
        return result;
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus, int nIterations, boolean quadratic) {
        testMultiFactorHenselLifting(base, modulus, newModulus(modulus, nIterations, quadratic), nIterations, quadratic);
    }

    static void testMultiFactorHenselLifting(bMutablePolynomialZ base, long modulus, int nIterations) {
        testMultiFactorHenselLifting(base, modulus, nIterations, true);
        testMultiFactorHenselLifting(base, modulus, nIterations, false);
    }

    @Test
    public void testHensel4() throws Exception {
        bMutablePolynomialZ base = bMutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 5, 5);
    }

    @Test
    public void testHensel5() throws Exception {
        bMutablePolynomialZ base = bMutablePolynomialZ.create(1, 2, 3, 4, 1, 6, 7, 1, 5, 4, 3, 2, 1);
        testMultiFactorHenselLifting(base, 11, 5);
    }

    @Test
    public void testHensel6_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int i = 0; i < its(1000, 10000); ++i) {
            bMutablePolynomialZ poly = RandomPolynomials.randomPoly(rndd.nextInt(5, 15), BigInteger.LONG_MAX_VALUE, rnd);
            for (long modulus : getModulusArray((int) its(5, 15), 0, 5, 0)) {
                MutablePolynomialMod polyMod = poly.modulus(modulus).toLong();
                if (polyMod.isConstant() || polyMod.isMonomial())
                    continue;

                while (poly.lc().mod(BigInteger.valueOf(modulus)).isZero())
                    poly.data[poly.degree] = BigInteger.valueOf(rnd.nextInt());

                try {
                    testMultiFactorHenselLifting(poly, modulus, 10);
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
    public void testStrategy_performance() throws Exception {
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


//    @Test
//    public void name() throws Exception {
//        bMutablePolynomialZ a = bMutablePolynomialZ.create(new BigInteger("12423527654545234"), new BigInteger("121422353452133312"), new BigInteger("-1231212412")).primitivePart();
//        bMutablePolynomialZ b = bMutablePolynomialZ.create(new BigInteger("-12312423527654545234"), new BigInteger("1214223534501333122"), new BigInteger("1223")).primitivePart();
//        bMutablePolynomialZ poly = a.clone().multiply(b);
//
//
//        BigInteger desiredBound = poly.mignotteBound().shiftLeft(1).multiply(poly.lc());
//        long modulus = 17;
//        BigInteger bModulus = BigInteger.valueOf(modulus);
//
//        MutablePolynomialMod aMod = a.modulus(modulus).monic().multiply(poly.lc()).toLong();
//        MutablePolynomialMod bMod = b.modulus(modulus).monic().toLong();
////        bMutablePolynomialMod[] lift = HenselLifting.liftQuadratic(modulus, poly, aMod, bMod,HenselLifting.nIterations(bModulus, desiredBound, true).nIterations);
//        bMutablePolynomialMod[] lift = lift(modulus, desiredBound, poly, aMod, bMod);
//        System.out.println(a);
//        System.out.println(lift[0].normalSymmetricForm().primitivePart());
//        System.out.println(b);
//        System.out.println(lift[1].normalSymmetricForm().primitivePart());
//
//    }

//    static final int always_linear = 4;
//    static final int linear_addon = 10;
//
//    private static bMutablePolynomialMod[] lift(
//            long modulus,
//            BigInteger desiredBound,
//            bMutablePolynomialZ poly,
//            MutablePolynomialMod aFactor,
//            MutablePolynomialMod bFactor) {
//
//        BigInteger bModulus = BigInteger.valueOf(modulus);
//
//        int nLinearIterations = HenselLifting.nIterations(bModulus, desiredBound, false).nIterations;
//        if (nLinearIterations <= always_linear)
//            // for small number of iterations we always use linear lift
//            return HenselLifting.liftLinear(modulus, poly, aFactor, bFactor, nLinearIterations);
//
//
//        BigInteger finalModulus;
//
//        // lift once with fast linear lift  => finalModulus = modulus^2
//        HenselLifting.bLinearLift linearLift = HenselLifting.createLinearLift(modulus, poly, aFactor, bFactor);
//        linearLift.lift();
//        finalModulus = bModulus.multiply(bModulus);
//
//        assertHenselLift(linearLift);
//
//        //calculate how much quadratic iterations are required:
//        BigInteger tmp = finalModulus, prevTmp = finalModulus;
//        int nQuadraticIterations = 0;
//        while (tmp.compareTo(desiredBound) < 0) {
//            prevTmp = tmp;
//            tmp = tmp.multiply(tmp);
//            ++nQuadraticIterations;
//        }
//
//        //now we need to decide whether to do all iterations quadratic or do last few iterations with linear lift
//
//        nLinearIterations = 0;
//        while (prevTmp.compareTo(desiredBound) < 0) {
//            prevTmp = prevTmp.multiply(bModulus);
//            ++nLinearIterations;
//            if (nLinearIterations > linear_addon)
//                break;
//        }
//
//        if (nLinearIterations <= linear_addon)
//            --nQuadraticIterations;
//
//        HenselLifting.bQuadraticLift quadraticLift =
//                new HenselLifting.bQuadraticLift(linearLift.modulus, poly, linearLift.aFactorMod(), linearLift.bFactorMod(), linearLift.aCoFactorMod(), linearLift.bCoFactorMod());
//        assertHenselLift(quadraticLift);
//        quadraticLift.liftWithCoFactors(nQuadraticIterations);
//        assertHenselLift(quadraticLift);
//
//        if (nLinearIterations > 0) {
//            HenselLifting.bLinearLift finalLift = new HenselLifting.bLinearLift(
//                    poly,
//                    quadraticLift.aFactorMod().normalForm(false),
//                    quadraticLift.bFactorMod().normalForm(false),
//                    quadraticLift.aCoFactorMod().normalForm(false),
//                    quadraticLift.bCoFactorMod().normalForm(false),
//                    aFactor,
//                    linearLift.aFactorModMonic,
//                    bFactor,
//                    linearLift.aCoFactorMod,
//                    linearLift.bCoFactorMod,
//                    bModulus,
//                    quadraticLift.aFactor.modulus);
//            finalLift.lift(nLinearIterations);
//            return new bMutablePolynomialMod[]{finalLift.aFactorMod(), finalLift.bFactorMod()};
//        } else
//            return new bMutablePolynomialMod[]{quadraticLift.aFactorMod(), quadraticLift.bFactorMod()};
//    }

    @Test
//    @Benchmark
    public void testLiftFactorization7_random() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        DescriptiveStatistics
                linear = new DescriptiveStatistics(),
                quadratic = new DescriptiveStatistics(),
                adaptive = new DescriptiveStatistics();

        int nFactors = 6;

        int nIterations = 500;
        for (int n = 0; n < nIterations; n++) {
            if (nIterations / 10 == n) {
                linear.clear(); quadratic.clear(); adaptive.clear();
            }
            if (n % 10 == 0)
                System.out.println(n);

            bMutablePolynomialZ poly = bMutablePolynomialZ.one();
            for (int i = 0; i < nFactors; i++)
                poly = poly.multiply(RandomPolynomials.randomPoly(15, BigInteger.valueOf(1000), rnd));
            poly = SquareFreeFactorization.SquareFreePart(poly);

            bMutablePolynomialMod polyMod;
            long modulus;
            do {
                modulus = getModulusRandom(rndd.nextInt(5, 28));
                polyMod = poly.modulus(modulus);
            } while (!SquareFreeFactorization.isSquareFree(poly.modulus(modulus)) || polyMod.degree() != poly.degree());

            BigInteger desiredBound = poly.mignotteBound().shiftLeft(1).multiply(poly.lc());


            lFactorDecomposition<MutablePolynomialMod> modularFactors = factor(polyMod.toLong());

            long start;

            try {
                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsQuad = HenselLifting.liftFactorization(BigInteger.valueOf(modulus), desiredBound, poly, modularFactors.factors, true);
                quadratic.addValue(System.nanoTime() - start);
                assertFactorization(poly.modulus(factorsQuad.get(0).modulus), poly.lc(), factorsQuad);

                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsLin = HenselLifting.liftFactorization(BigInteger.valueOf(modulus), desiredBound, poly, modularFactors.factors, false);
                linear.addValue(System.nanoTime() - start);
                assertFactorization(poly.modulus(factorsLin.get(0).modulus), poly.lc(), factorsLin);


                start = System.nanoTime();
                List<bMutablePolynomialMod> factorsAdp = HenselLifting.liftFactorizationAdaptive(modulus, desiredBound, poly, modularFactors.factors);
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

    @Test
    public void test8() throws Exception {
        bMutablePolynomialZ a = bMutablePolynomialZ.create(1, 2, 3, 5, 3, 2, 1),
                b = bMutablePolynomialZ.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = bMutablePolynomialZ.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = bMutablePolynomialZ.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = bMutablePolynomialZ.create(-11, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = bMutablePolynomialZ.create(-33, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        bMutablePolynomialZ poly = a.clone().multiply(b, c, d, e, f, g, g.clone().increment(), f.clone().increment());


        BigInteger bound2 = poly.mignotteBound().shiftLeft(1);
        long modulus = SmallPrimes.nextPrime((1 << 24) + getRandom().nextInt((1 << 30) - (1 << 24)));

        BigInteger bModulus = BigInteger.valueOf(modulus);
        System.out.println(HenselLifting.nIterations(bModulus, bound2, true).nIterations);
        System.out.println(HenselLifting.nIterations(bModulus, bound2, false).nIterations);
    }

    @Test
    public void test9() throws Exception {

        long modulus = 47;

        bMutablePolynomialZ poly = bMutablePolynomialZ.create(new BigInteger("-177177428820536591540197383566490586277126749363929600"), new BigInteger("-1242657761787104232614791039136879496461412157650740032"), new BigInteger("-3531511936294694064576189424935490231581922528167398384"), new BigInteger("-6084332472010238642690531118562682194933790282930088528"), new BigInteger("-8155005410062241493973240591280440926232487083775272016"), new BigInteger("-7839734274576181093562990906498360629649729979582690786"), new BigInteger("-4119348054365418839188718893550617808924749793249241838"), new BigInteger("-2712242727268458360861305710699004712722364939263631786"), new BigInteger("1227216583725964856823790760658316098620048897741354534"), new BigInteger("24216429013482071530263092295805264001700345143490431746"), new BigInteger("57671175987391189945931547837750570813205372792009964557"), new BigInteger("87200443377506038701107745628666899803969637223604694229"), new BigInteger("123440595280492568817512621770439276863460170374867019756"), new BigInteger("141770990209049511779896925094689465658428943289230136776"), new BigInteger("189550330769616265235959109762342671695099325855439757554"), new BigInteger("253066664492310181423363373361886052180627702748343049872"), new BigInteger("349690647617511993952396995696207383721081506571270980034"), new BigInteger("385456540302527465810061000786893447210367730308081138880"), new BigInteger("526866589182598293360556604021941129984912695269353483172"), new BigInteger("765797785491342931353682694929398362438502639720493473156"), new BigInteger("966971357443907258082284841995269669603555146493004574174"), new BigInteger("1168265250555763072891035530664055486734808802908000647754"), new BigInteger("1269209241381104745347036867940121193000545757835221565430"), new BigInteger("1523535283048282686017926990086174818918281441025286597837"), new BigInteger("1677472657492590273768473003309222712393329316924511265599"), new BigInteger("2083913978522957720466548156842310688231217834524452688610"), new BigInteger("2266127869204685322167064924772588655515387725483542308342"), new BigInteger("2784407394732921716498948785209554459890861766467556821452"), new BigInteger("3104179084632504150033843744335052026437222778870064551980"), new BigInteger("4235682671687461394918042331750328317791767115259774262131"), new BigInteger("4557025057925223504438748406962807333333149796916369188912"), new BigInteger("5092588002693402313118369660025982089982733429946756571159"), new BigInteger("5314768197913419577870355709508994301610513735390392713370"), new BigInteger("6252966863454040256113320182427111998278403516315588173142"), new BigInteger("6670913992742472496941111278925909264435344705043341681550"), new BigInteger("7884753946704403783653549706041937369552018964643863189666"), new BigInteger("9008837449392851768612381498404297215968396324085662967388"), new BigInteger("9336003989485608390135320355234083929457851356098142275465"), new BigInteger("10097566287661235636580479618188603302810021003071530833708"), new BigInteger("11022623760470402033164736077095766649180620012058536418387"), new BigInteger("11239780324726301886300638634942157048717648229470414508312"), new BigInteger("12289809156000700748415279405643438073739864230381362746881"), new BigInteger("13852149720504256772086228934587254962187955653551299958355"), new BigInteger("14725309233915842775834182407584085701881528687756330321758"), new BigInteger("15074080241605073700764800505544424302028253949756050974640"), new BigInteger("15804766063468982057041136185763564553319562465024476451909"), new BigInteger("16418308941221585820167221517027435991224496333709144296904"), new BigInteger("16413343207312873176572266765701880034124955472002133144905"), new BigInteger("18369813255638839310163461623813778338295129145744693186521"), new BigInteger("18949852721822855419925853256822131371750658301732527667720"), new BigInteger("18137458155980050467357864684759740609092038978136646628529"), new BigInteger("19507072053261047871615982347484111140161893573587177550359"), new BigInteger("18554703470549512001213120213000104847376351637194086031597"), new BigInteger("18624570538826514850806550075687383127155673342485166070664"), new BigInteger("20061831661729964166954857551615379954105548688030028460276"), new BigInteger("19239909520601616457455981279960446708581064930222996221881"), new BigInteger("19290037029536733425297368874881224166397218373218437576240"), new BigInteger("19538607683391486694419903458354668605396142434730251119248"), new BigInteger("18101586016275652703682443628809054015794575036970160804332"), new BigInteger("17906945055103474583220176287452413362258828449718169553332"), new BigInteger("17831172966650169243461404677746016583992656879256999955015"), new BigInteger("16675177599596883902630921509123643652600598728688912352347"), new BigInteger("15945051454879083236187822786350979743166673855126384118532"), new BigInteger("14970679657524904323134848413518297943976032930575660518079"), new BigInteger("13103369159971691529355326700363794023077614526225223906407"), new BigInteger("12032968169639164958030857999802636726396680665604981501439"), new BigInteger("11568391488664415241918649716346772138228730461790411362728"), new BigInteger("10052500280149582491392799958108504415447798547792066931648"), new BigInteger("9090125536215949820264014431857460526208269773417578203036"), new BigInteger("7797628304143296786827206593371923724914405441259925079394"), new BigInteger("6494128693586457608302291176996842792789538277952110578478"), new BigInteger("5969815195253556590182696813056320001710282017000357861979"), new BigInteger("5002554504357626647949739865343629345265769463671696897444"), new BigInteger("4255610300133212741425481905605300821899878576174890868044"), new BigInteger("3490544760353135530679689751314055469331465124914418050913"), new BigInteger("2706421822313658365510867773283920195964495673025481429946"), new BigInteger("2228343481180518545763889405070938452743092715783943938920"), new BigInteger("1692918918745349977967263518198756388947291164771144787124"), new BigInteger("1340650058000612912542776587313434050794685378401371384016"), new BigInteger("969056777548212700002722865690942291519058912131423810087"), new BigInteger("640072278163165090143492947374467346535008127480197450509"), new BigInteger("490000958445739055624285102462829799131708617146840664716"), new BigInteger("284674715458815006664121504826441837470031106451664063706"), new BigInteger("183154936426515710323970496340799967708193589538342385210"), new BigInteger("125563861024650235940594703844208936124009621318618008264"), new BigInteger("48882336833837794768720001359770640671506214813691426632"), new BigInteger("28452604269975079673827028056891217079753980379117303480"), new BigInteger("11529952705667314043872877825589099431569872299323184448"), new BigInteger("3524331934626291661437350174828103426365065122099201984"), new BigInteger("5052775290319813160076968738276160892860747696737155200"), new BigInteger("1518494176512432373504770375706925168293178376011458560"));
        lFactorDecomposition<MutablePolynomialMod> factors = factor(poly.modulus(modulus).toLong());

        System.out.println(poly.degree);
        System.out.println(poly.modulus(modulus).degree);
        List<bMutablePolynomialMod> lifted = HenselLifting.liftFactorization(BigInteger.valueOf(modulus), poly.mignotteBound(), poly, factors.factors, true);
    }

    //    @Test
//    public void adasdasd() throws Exception {
//        bMutablePolynomialZ aFactor = bMutablePolynomialZ.create(12, 13112, 1112, 31, 112);
//        bMutablePolynomialZ bFactor = bMutablePolynomialZ.create(22, 112311, 12);
//        bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
//
////        System.out.println(aFactor);
////        System.out.println(bFactor);
//
//        BigInteger modulus = BigInteger.valueOf(17);
//        aFactor = aFactor.modulus(modulus).monic().normalForm(true).multiply(poly.lc());
//        bFactor = bFactor.modulus(modulus).monic().normalForm(true);
//
//        HenselLifting.bLinearLift lift = createLinearLift(modulus.longValue(), poly,
//                aFactor.modulus(modulus).toLong(),
//                bFactor.modulus(modulus).toLong());
////        System.out.println(lift.poly);
//        assertHenselLift(lift);
//        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
//        System.out.println(lift.aFactor);
//        System.out.println(lift.bFactor);
//
//        BigInteger lmodulus = BigInteger.valueOf(17);
//
//        lift.lift();
//        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
//        System.out.println(lift.aFactor);
//        System.out.println(lift.bFactor);
//        lift.lift();
//        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
//        System.out.println(lift.aFactor);
//        System.out.println(lift.bFactor);
//
//        lift.lift();
//        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
//        System.out.println(lift.aFactor);
//        System.out.println(lift.bFactor);
//    }

    private static <T extends IMutablePolynomialZp<T>> void assertHenselLift(LiftableQuintet<T> lift) {
        Assert.assertEquals(lift.toString(), lift.polyMod(), lift.aFactorMod().clone().multiply(lift.bFactorMod()));
        Assert.assertTrue(lift.toString(), lift.aCoFactorMod() == null && lift.bCoFactorMod() == null ||
                lift.aFactorMod().clone().multiply(lift.aCoFactorMod())
                        .add(lift.bFactorMod().clone().multiply(lift.bCoFactorMod()))
                        .isOne());
    }
}