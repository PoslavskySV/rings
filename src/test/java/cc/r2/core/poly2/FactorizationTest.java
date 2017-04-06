package cc.r2.core.poly2;

import cc.r2.core.combinatorics.IntCombinationsGenerator;
import cc.r2.core.number.BigInteger;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.test.Benchmark;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.List;

import static cc.r2.core.poly2.Factorization.*;
import static cc.r2.core.poly2.FactorizationTestUtil.assertFactorization;
import static cc.r2.core.poly2.PolynomialArithmetics.polyMod;
import static cc.r2.core.poly2.SquareFreeFactorization.isSquareFree;
import static cc.r2.core.util.TimeUnits.nanosecondsToString;

/**
 * Created by poslavsky on 27/02/2017.
 */
public class FactorizationTest extends AbstractPolynomialTest {
    @Test
    public void test1() throws Exception {
        Assert.assertTrue(factor(MutablePolynomialZ.create(3, 7).modulus(19)).get(0).isMonic());
    }

    @Test
    public void test2() throws Exception {
        BigInteger modulus = BigInteger.LONG_MAX_VALUE;
        modulus = modulus.multiply(modulus).increment().nextProbablePrime();
        bMutablePolynomialMod poly = bMutablePolynomialZ.create(
                BigInteger.valueOf(Long.MAX_VALUE),
                BigInteger.valueOf(Long.MAX_VALUE - 1),
                BigInteger.valueOf(Long.MAX_VALUE - 2)).modulus(modulus);
        for (int i = 0; i < 5; i++)
            poly = poly.square().add(poly.derivative()).increment();
        bFactorDecomposition<bMutablePolynomialMod> fct = factor(poly);
        Assert.assertEquals(7, fct.size());
        assertFactorization(poly, fct);
    }

    @Test
    public void test3() throws Exception {
        long modulus = 13;
        MutablePolynomialMod poly = MutablePolynomialZ.create(5, 8, 1, 5, 7, 0, 0, 1, 5, 7, 0, 9, 3, 2).modulus(modulus);
        assertFactorization(poly, factor(poly));
    }

//    static <PolyZp extends IMutablePolynomialZp<PolyZp>> void assertHenselLift(QuadraticLiftAbstract<PolyZp> lift) {
//        Assert.assertEquals(lift.toString(), lift.polyMod(), lift.aFactor.clone().multiply(lift.bFactor));
//        Assert.assertEquals(lift.toString(), lift.polyMod().createOne(),
//                lift.aFactor.clone().multiply(lift.aCoFactor)
//                        .add(lift.bFactor.clone().multiply(lift.bCoFactor)));
//    }

    private static <T extends IMutablePolynomialZp<T>> void assertHenselLift(HenselLifting.QuadraticLiftAbstract<T> lift) {
        Assert.assertEquals(lift.toString(), lift.polyMod(), lift.aFactor.clone().multiply(lift.bFactor));
        Assert.assertTrue(lift.toString(), lift.aCoFactor == null && lift.bCoFactor == null
                || lift.aFactor.clone().multiply(lift.aCoFactor).add(lift.bFactor.clone().multiply(lift.bCoFactor)).isOne());
    }


    @Test
    public void testFactorization1() throws Exception {
        MutablePolynomialZ a = MutablePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
        MutablePolynomialZ poly = a.multiply(b).primitivePart();
        Assert.assertTrue(SquareFreeFactorization.isSquareFree(poly));

        GlobalRandom.getRandom().setSeed(685922130507849253L);
    }

    //
//    @Test
//    public void testHensel5() throws Exception {
//        MutablePolynomialZ base = MutablePolynomialZ.create(69, 30);
//        long modulus = 29;
//
//        MutablePolynomialMod baseMod = base.modulus(modulus);
//        System.out.println(factor(baseMod));
//        List<MutablePolynomialMod> fact = liftFactorization(modulus, 1, base, Collections.singletonList(baseMod.monic()));
//        MutablePolynomialMod sin = fact.get(0);
//        System.out.println(sin.multiply(base.lc()));
//
//
//        System.out.println(factor(base.modulus(modulus)));
//        testMultiFactorHenselLifting(base, modulus, 1);
//
//    }

//
//    @Test
//    public void test1() throws Exception {
//        long modulus = 127;
//        MutablePolynomialMod poly =
//                MutablePolynomialZ.create(1, 2, 3, 4).modulus(modulus)
//                        .multiply(MutablePolynomialZ.create(1, 2, 3, 4, 5).modulus(modulus))
//                        .square()
//                        .multiply(MutablePolynomialZ.create(2, 3, 4).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 34).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, 32, 42, 12, 23, 31).modulus(modulus))
//                        .multiply(MutablePolynomialZ.create(22, -32, 42, 12, 23, 31).modulus(modulus));
//        poly = poly.multiply(poly.shiftLeft(5)).multiply(poly.shiftRight(5));
//
//        FactorDecomposition<MutablePolynomialMod> factorization = Factorization.factor(poly);
//        System.out.println(factorization);
//        System.out.println(poly);
//        System.out.println(factorization.toPolynomial(poly));
//
//    }
//

    @Ignore
    @Test
    public void name() throws Exception {
//        System.out.println(factorBigPrime(MutablePolynomialZ.create(-1, 0, 1)));
//        System.out.println(factorBigPrime(MutablePolynomialZ.create(-4, 0, 1)));

        MutablePolynomialZ a = MutablePolynomialZ.create(1, 11, 121, 1, 1, 1, 1, 123);
        MutablePolynomialZ b = MutablePolynomialZ.create(1, 1, 1, 3, 1, 1, 2, 3, 4, 5, 6);
        MutablePolynomialZ poly = a.multiply(b).primitivePart();
        Assert.assertTrue(SquareFreeFactorization.isSquareFree(poly));

//        MutablePolynomialMod moduloImage = poly.modulus(6101).monic();
//        System.out.println(moduloImage);
////        System.out.println(poly);
//        System.out.println(factor(moduloImage));
//        FactorDecomposition<MutablePolynomialMod> modFact = factor(moduloImage.monic());
//        System.out.println(modFact);

//        System.out.println(modFact.get(0).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(0).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(2).clone().multiply(modFact.get(1)).multiply(poly.lc()).normalSymmetricForm());
//        System.out.println(modFact.get(2).clone().multiply(modFact.get(3)).multiply(poly.lc()).normalSymmetricForm());
//
//        System.out.println(modFact);

//        poly = MutablePolynomialZ.create(46225, 0, -5596840, 0, 13950764, 0, -7453176, 0, 1513334, 0, -141912, 0, 6476, 0, -136, 0, 1);
        System.out.println(poly);

        System.out.println(Factorization.factorSquareFree(poly));

//        if(true)
//            return;
        int k = 0;
        for (int i = 0; i < 10000; i++) {
            if (i % 100 == 0)
                System.out.println(i);
            long seed = GlobalRandom.getRandom().nextLong();
            GlobalRandom.getRandom().setSeed(seed);
            long start = System.nanoTime();
            Assert.assertEquals(2, Factorization.factorSquareFree(poly).size());
            double time = ((System.nanoTime() - start) / 1000. / 1000.);
            if (time > 1000) {
                System.out.println("======   " + time);
                System.out.println(seed);
                break;
            }
            if (i > 1000)
                System.out.println(time);
        }
        System.out.println(k);

//
//        0.31519400000000003
//        0.37378500000000003
//        0.419362
//        0.358527
//        0.415434
//        0.520177
//        0.534273
//        0.401569
//        0.376355

    }

    @Test
    @Benchmark
    public void tmpcompare() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        rnd.setSeed(123);
        long modulus = getModulusRandom(15);
        BigInteger bigModulus = BigInteger.LONG_MAX_VALUE;
        bigModulus = bigModulus.multiply(bigModulus).nextProbablePrime();

        DescriptiveStatistics modT = new DescriptiveStatistics(), bT = new DescriptiveStatistics();
        int nIterations = 10000;
        for (int i = 0; i < nIterations; i++) {
            if (i % 100 == 0)
                System.out.println(i);
            if (i == nIterations / 10) {
                modT.clear();
                bT.clear();
            }
            MutablePolynomialMod poly = RandomPolynomials.randomMonicPoly(rndd.nextInt(10, 15), modulus, rnd);
            long start = System.nanoTime();
            lFactorDecomposition<MutablePolynomialMod> fct = factor(poly);
            modT.addValue(System.nanoTime() - start);
            assertFactorization(poly, fct);


            if (i % 100 == 0) {
                bMutablePolynomialMod bPoly = poly.toBigPoly().setModulusUnsafe(bigModulus);
                start = System.nanoTime();
                bFactorDecomposition<bMutablePolynomialMod> bfct = factor(bPoly);
                bT.addValue(System.nanoTime() - start);
                assertFactorization(bPoly, bfct);
            }
        }

        System.out.println("Machine-size:");
        System.out.println(modT);
        System.out.println("BigInteger:");
        System.out.println(bT);
    }

    private static int randomModulusInf() {
        return (1 << 20) + GlobalRandom.getRandom().nextInt((1 << 30) - (1 << 20));
    }

    @Test
    public void test12() throws Exception {
        long modulus = 689177113;
        BigInteger bModulus = BigInteger.valueOf(modulus);

        bMutablePolynomialZ a = bMutablePolynomialZ.create(1, 2, 3, 5, 3, 2, 1),
                b = bMutablePolynomialZ.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = bMutablePolynomialZ.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = bMutablePolynomialZ.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = bMutablePolynomialZ.create(-11, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = bMutablePolynomialZ.create(-33, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        bMutablePolynomialZ poly = a.clone().multiply(b, c, d, e, f, g, g.clone().increment(), f.clone().increment());
        Assert.assertTrue(isSquareFree(poly));

        MutablePolynomialMod polyMod = poly.modulus(modulus).toLong();
        Assert.assertTrue(isSquareFree(polyMod));

        BigInteger bound2 = poly.mignotteBound().shiftLeft(1);
        lFactorDecomposition<MutablePolynomialMod> lModFactors = factor(polyMod);
        System.out.println(HenselLifting.nIterations(BigInteger.valueOf(modulus), bound2,false).nIterations);

        List<bMutablePolynomialMod> bModFactorsL = HenselLifting.liftFactorization(bModulus, bound2, poly, lModFactors.factors, false);
        List<bMutablePolynomialMod> bModFactorsA = HenselLifting.liftFactorizationAdaptive(modulus, bound2, poly, lModFactors.factors);

        System.out.println(bound2);
        System.out.println(bModFactorsL.get(0).modulus);
        System.out.println(bModFactorsA.get(0).modulus);

        System.out.println(Factorization.reconstructFactorsZ(poly, new bFactorDecomposition<>(bModFactorsL, BigInteger.ONE)));
        System.out.println(Factorization.reconstructFactorsZ(poly, new bFactorDecomposition<>(bModFactorsA, BigInteger.ONE)));
    }

    @Test
//    @Benchmark
    public void testZxPerformance() throws Exception {
        bMutablePolynomialZ a = bMutablePolynomialZ.create(1, 2, 3, 5, 3, 2, 1),
                b = bMutablePolynomialZ.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = bMutablePolynomialZ.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = bMutablePolynomialZ.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = bMutablePolynomialZ.create(-11, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = bMutablePolynomialZ.create(-33, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        bMutablePolynomialZ poly = a.clone().multiply(b, c, d, e, f, g, g.clone().increment(), f.clone().increment());
        System.out.println(poly);

//        poly.increment();
//        long[] primes = {3, 5, 7, 11, 13, 17, 19, 23, SmallPrimes.nextPrime(12343), SmallPrimes.nextPrime(123413)};
//        for (long prime : primes) {
//            System.out.println(factor(poly.modulus(prime).toLong()).size());
//        }
//
//        for (int i = 0; i < 10; i++) {
//            System.out.println(factor(poly.modulus(SmallPrimes.nextPrime(randomModulusInf())).toLong()).size());
//        }
//
//        if(true) return;
        BigInteger x = poly.mignotteBound();
        BigInteger r = BigInteger.ONE.shiftLeft(15);
        for (int i = 0; ; i++) {
            r = r.multiply(r);
            if (r.compareTo(x) >= 0)
                break;
            System.out.println(i);
        }
//        System.out.println(x);
//        System.out.println(poly.norm2Double() * Math.pow(2, poly.degree));

//        if (true) return;
        N_MODULAR_FACTORIZATION_TRIALS = 2;
        System.out.println(poly);
        for (int i = 0; i < 400; i++) {
            if (i == 10) {
                t_CZ = t_LIFT = t_RECOMB = 0;
            }
//            System.out.println(i);
            long start = System.nanoTime();
            bFactorDecomposition<bMutablePolynomialZ> fct = factorSquareFree(poly);
            if (i > 10) {
                System.out.println("Ti: " + nanosecondsToString(System.nanoTime() - start));
//                System.out.println("CZ: " + nanosecondsToString(t_CZ));
//                System.out.println("HE: " + nanosecondsToString(t_LIFT));
//                System.out.println("RE: " + nanosecondsToString(t_RECOMB));
//
//                System.out.println("HE_MULT: " + nanosecondsToString(t_HE_MUL));
//                System.out.println("HE_LIFT: " + nanosecondsToString(t_HE_LIFT));
//                System.out.println("HE_CREA: " + nanosecondsToString(t_HE_CREATE));
                t_CZ = t_LIFT = t_RECOMB = t_HE_LIFT = t_HE_CREATE = t_HE_MUL = 0;
            }
            Assert.assertEquals(9, fct.size());
            assertFactorization(poly, fct);
        }
    }

    @Test
    public void test4() throws Exception {

        bMutablePolynomialZ a = bMutablePolynomialZ.create(1, 2, 3, 5, 3, 2, 1),
                b = bMutablePolynomialZ.create(1, 2, -12443241213L, 412312, 3, 2, 123423554351L),
                c = bMutablePolynomialZ.create(-1, -2, -12443241213L, 412312, 3, 2, 123423554351L),
                d = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L),
                e = bMutablePolynomialZ.create(-11111111, -2, -12441213L, 412312, 3, 2, 1234235543L),
                f = bMutablePolynomialZ.create(-1, -2222121, -12441213L, 412312, 3, 2, 1234235543L),
                g = bMutablePolynomialZ.create(-1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L, -1, -2, -12441213L, 412312, 3, 2, 1234235543L, -12441213L, 412312, 3, 2, 1234235543L);
        bMutablePolynomialZ poly = a.clone().multiply(b, c, d, e, f, g);
        long modulus = 49566989;
        MutablePolynomialMod modPoly = poly.modulus(modulus).toLong();
        Assert.assertTrue(SquareFreeFactorization.isSquareFree(modPoly));

        bFactorDecomposition<bMutablePolynomialMod> modularFactors = lFactorDecomposition.convert(factor(modPoly.monic()));
        modularFactors.canonicalForm();
        System.out.println(modularFactors);
        bMutablePolynomialMod factory = modularFactors.get(0);
        bMutablePolynomialMod
                aFactor = factory.createConstant(poly.lc()),
                bFactor = factory.createOne();

        int nHalf = modularFactors.size() / 2, i = 0;
        for (; i < nHalf; ++i)
            aFactor = aFactor.multiply(modularFactors.get(i));
        for (; i < modularFactors.size(); ++i)
            bFactor = bFactor.multiply(modularFactors.get(i));


        HenselLifting.bQuadraticLift henselInput = HenselLifting.createQuadraticLift(BigInteger.valueOf(modulus), poly, aFactor, bFactor);
        assertHenselLift(henselInput);
        System.out.println(henselInput);
    }

    @Test
    public void azzzname() throws Exception {
        int[] set = {121, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
        TIntHashSet sums = new TIntHashSet();
        for (int i = 1; i < set.length; i++) {
            IntCombinationsGenerator gen = new IntCombinationsGenerator(set.length, i);
            int[] com;
            while ((com = gen.take()) != null) {
                int sum = 0;
                for (int j = 0; j < com.length; j++)
                    sum += set[com[j]];

                sums.add(sum);
            }
        }
        System.out.println(sums);
        System.out.println(sums.size());
        System.out.println(1 << set.length);
    }

    @Test
    public void adasdasd_QUAD() throws Exception {
        bMutablePolynomialZ aFactor = bMutablePolynomialZ.create(12, 13112, 1112, 31, 112);
        bMutablePolynomialZ bFactor = bMutablePolynomialZ.create(22, 112311, 1);
        bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);

//        System.out.println(aFactor);
//        System.out.println(bFactor);

        BigInteger modulus = BigInteger.valueOf(17);

        HenselLifting.bQuadraticLift he = HenselLifting.createQuadraticLift(modulus, poly, aFactor.modulus(modulus).monic().multiply(poly.lc()), bFactor.modulus(modulus).monic());
//        System.out.println(lift.poly);
        assertHenselLift(he);
        he.lift();
        System.out.println(he.aFactor);
        System.out.println(he.bFactor);
        assertHenselLift(he);
        he.lift();
        System.out.println(he.aFactor);
        System.out.println(he.bFactor);
    }

    @Test
    public void adasdasd() throws Exception {
        bMutablePolynomialZ aFactor = bMutablePolynomialZ.create(12, 13112, 1112, 31, 112);
        bMutablePolynomialZ bFactor = bMutablePolynomialZ.create(22, 112311, 12);
        bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);

//        System.out.println(aFactor);
//        System.out.println(bFactor);

        BigInteger modulus = BigInteger.valueOf(17);
        aFactor = aFactor.modulus(modulus).monic().normalForm(true).multiply(poly.lc());
        bFactor = bFactor.modulus(modulus).monic().normalForm(true);

        LinearLift lift = createLift(modulus, poly, aFactor, bFactor);
//        System.out.println(lift.poly);
        assertHenselLift(lift);
        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
        System.out.println(lift.aFactor);
        System.out.println(lift.bFactor);

        BigInteger lmodulus = BigInteger.valueOf(17);

        lift.lift();
        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
        System.out.println(lift.aFactor);
        System.out.println(lift.bFactor);
        lift.lift();
        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
        System.out.println(lift.aFactor);
        System.out.println(lift.bFactor);

        lift.lift();
        System.out.println(lift.poly.clone().subtract(lift.aFactor.clone().multiply(lift.bFactor)).modulus(lift.modulus));
        System.out.println(lift.aFactor);
        System.out.println(lift.bFactor);
    }


    @Test
    public void testLinearLiftRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (int n = 0; n < 1000; n++) {
            bMutablePolynomialZ aFactor = RandomPolynomials.randomPoly(rndd.nextInt(3, 5), BigInteger.valueOf(100), rnd);
            bMutablePolynomialZ bFactor = RandomPolynomials.randomPoly(rndd.nextInt(3, 5), BigInteger.valueOf(100), rnd);
            if (!PolynomialGCD.PolynomialGCD(aFactor, bFactor).isConstant()) {
                --n;
                continue;
            }
            bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
            BigInteger modulus = BigInteger.valueOf(getModulusRandom(10));
            aFactor = aFactor.modulus(modulus).monic().normalForm(true).multiply(poly.lc());
            bFactor = bFactor.modulus(modulus).monic().normalForm(true);
            if (!PolynomialGCD.PolynomialGCD(aFactor.modulus(modulus), bFactor.modulus(modulus)).isConstant()) {
                --n;
                continue;
            }

            LinearLift lift = createLift(modulus, poly, aFactor, bFactor);
            assertHenselLift(lift);
            for (int i = 0; i < 10; i++) {
                lift.lift();
                assertHenselLift(lift);
            }
        }
    }

    /** creates liftable quintet */
    static LinearLift createLift(BigInteger modulus,
                                 bMutablePolynomialZ poly,
                                 bMutablePolynomialZ aFactor,
                                 bMutablePolynomialZ bFactor) {
        bMutablePolynomialMod[] xgcd = HenselLifting.monicExtendedEuclid(aFactor.modulus(modulus), bFactor.modulus(modulus));
        return new LinearLift(modulus, poly, aFactor, bFactor, xgcd[1].normalForm(true), xgcd[2].normalForm(true));
    }


    private static void assertHenselLift(LinearLift lift) {
        Assert.assertEquals(lift.toString(), lift.baseModulus(), lift.aFactor.clone().multiply(lift.bFactor).modulus(lift.modulus));
        Assert.assertTrue(lift.toString(), lift.aCoFactor == null && lift.bCoFactor == null
                || lift.aFactor.clone().multiply(lift.aCoFactor).add(lift.bFactor.clone().multiply(lift.bCoFactor)).modulus(lift.modulus).isOne());
    }

    @Test
    public void testdivv() throws Exception {
        BigInteger bound2 = new BigInteger("528143776232283612903301650916953940598417708048841122384169806827815822163968");
        System.out.println(bound2.bitLength());
        int linearIts = 0, quadIts = 0;

        BigInteger bModulus = BigInteger.valueOf(SmallPrimes.nextPrime(419566989));
        BigInteger liftedModulus = bModulus;
        while (liftedModulus.compareTo(bound2) < 0) {
            liftedModulus = liftedModulus.multiply(liftedModulus);
            ++quadIts;
        }
        System.out.println("Quad modulus: " + liftedModulus);

        liftedModulus = bModulus;
        while (liftedModulus.compareTo(bound2) < 0) {
            liftedModulus = liftedModulus.multiply(bModulus);
            ++linearIts;
        }
        System.out.println("Lin modulus: " + liftedModulus);


        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        DescriptiveStatistics quad = new DescriptiveStatistics(), lin = new DescriptiveStatistics();
        for (int n = 0; n < 1000; n++) {
//            System.out.println(n);
            if (n / 10 == 0) {
                quad.clear(); lin.clear();
            }

            if (n % 100 == 0)
                System.out.println(n);
            bMutablePolynomialZ aFactor = RandomPolynomials.randomPoly(rndd.nextInt(10, 15), BigInteger.valueOf(100), rnd);
            bMutablePolynomialZ bFactor = RandomPolynomials.randomPoly(rndd.nextInt(10, 15), BigInteger.valueOf(100), rnd);
            if (!PolynomialGCD.PolynomialGCD(aFactor, bFactor).isConstant()) {
                --n;
                continue;
            }
            bMutablePolynomialZ poly = aFactor.clone().multiply(bFactor);
            BigInteger modulus = bModulus;
            aFactor = aFactor.modulus(modulus).monic().normalForm(true).multiply(poly.lc());
            bFactor = bFactor.modulus(modulus).monic().normalForm(true);
            if (!PolynomialGCD.PolynomialGCD(aFactor.modulus(modulus), bFactor.modulus(modulus)).isConstant()) {
                --n;
                continue;
            }

            long start = System.nanoTime();
            LinearLift linLift = createLift(modulus, poly, aFactor, bFactor);
//            assertHenselLift(linLift);
            for (int i = 0; i < linearIts; i++) {
                linLift.lift();
//                assertHenselLift(linLift);
            }
            lin.addValue(System.nanoTime() - start);

            start = System.nanoTime();
//            System.out.println(modulus);
//            System.out.println(poly.toStringForCopy());
//            System.out.println(aFactor.modulus(modulus).toStringForCopy());
//            System.out.println(bFactor.modulus(modulus).toStringForCopy());
            HenselLifting.bQuadraticLift quadLift = HenselLifting.createQuadraticLift(modulus, poly, aFactor.modulus(modulus), bFactor.modulus(modulus));
//            assertHenselLift(quadLift);
            for (int i = 0; i < quadIts; i++) {
                quadLift.lift();
//                assertHenselLift(quadLift);
            }
            quad.addValue(System.nanoTime() - start);

            Assert.assertTrue(linLift.modulus.compareTo(bound2) > 0);
            Assert.assertTrue(quadLift.modulus.compareTo(bound2) > 0);
        }

        System.out.println("Quadratic: " + nanosecondsToString((long) quad.getMean()) + " ± " + nanosecondsToString((long) quad.getStandardDeviation()));
        System.out.println("Linear: " + nanosecondsToString((long) lin.getMean()) + " ± " + nanosecondsToString((long) lin.getStandardDeviation()));
    }

    static final class LinearLift {
        final BigInteger initialModulus;
        final bMutablePolynomialZ poly;
        BigInteger modulus;
        bMutablePolynomialZ aFactor, bFactor, aCoFactor, bCoFactor;
        final MutablePolynomialMod aFactorMod, aFactorModMonic, bFactorMod, aCoFactorMod, bCoFactorMod;
        final DivisionWithRemainder.InverseModMonomial<MutablePolynomialMod> aFactorModMonicInv, bFactorModInv;

        public LinearLift(BigInteger modulus, bMutablePolynomialZ poly, bMutablePolynomialZ aFactor, bMutablePolynomialZ bFactor, bMutablePolynomialZ aCoFactor, bMutablePolynomialZ bCoFactor) {
            assert modulus.isLong();
            assert bFactor.isMonic();
            this.poly = poly;
            this.modulus = modulus;
            this.aFactor = aFactor;
            this.bFactor = bFactor;
            this.aCoFactor = aCoFactor;
            this.bCoFactor = bCoFactor;
            this.initialModulus = modulus;
            aFactorMod = aFactor.modulus(initialModulus).toLong();
            aFactorModMonic = aFactorMod.clone().monic();
            bFactorMod = bFactor.modulus(initialModulus).toLong();
            aCoFactorMod = aCoFactor.modulus(initialModulus).toLong();
            bCoFactorMod = bCoFactor.modulus(initialModulus).toLong();
            aFactorModMonicInv = DivisionWithRemainder.fastDivisionPreConditioning(aFactorModMonic);
            bFactorModInv = DivisionWithRemainder.fastDivisionPreConditioning(bFactorMod);
        }

        bMutablePolynomialMod baseModulus() {
            return poly.modulus(modulus);
        }

        void lift() {
            BigInteger newModulus = modulus.multiply(initialModulus);
            bMutablePolynomialZ c = poly.clone().subtract(aFactor.clone().multiply(bFactor)).divideOrNull(modulus);
            MutablePolynomialMod cc = c.modulus(initialModulus).toLong();


            MutablePolynomialMod aAdd = cc.clone();
            aAdd = polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);
            aAdd = aAdd.multiply(bCoFactorMod);
            aAdd = polyMod(aAdd, aFactorModMonic, aFactorModMonicInv, false);

            MutablePolynomialMod bAdd = cc.clone();
            bAdd = polyMod(bAdd, bFactorMod, bFactorModInv, false);
            bAdd = bAdd.multiply(aCoFactorMod);
            bAdd = polyMod(bAdd, bFactorMod, bFactorModInv, false);

            aFactor = aFactor.add(aAdd.normalForm(false).toBigPoly().multiply(modulus));
            bFactor = bFactor.add(bAdd.normalForm(false).toBigPoly().multiply(modulus));

            bMutablePolynomialZ r = aCoFactor.clone().multiply(aFactor).add(bCoFactor.clone().multiply(bFactor)).decrement().negate().divideOrNull(modulus);
            MutablePolynomialMod rr = r.modulus(initialModulus).toLong();

            aAdd = rr.clone();
            aAdd = polyMod(aAdd, bFactorMod, bFactorModInv, false);
            aAdd = aAdd.multiply(aCoFactorMod);
            aAdd = polyMod(aAdd, bFactorMod, bFactorModInv, false);

            bAdd = rr.clone();
            bAdd = polyMod(bAdd, aFactorModMonic, aFactorModMonicInv, false);
            bAdd = bAdd.multiply(bCoFactorMod);
            bAdd = polyMod(bAdd, aFactorModMonic, aFactorModMonicInv, false);

            aCoFactor = aCoFactor.add(aAdd.normalForm(false).toBigPoly().multiply(modulus));
            bCoFactor = bCoFactor.add(bAdd.normalForm(false).toBigPoly().multiply(modulus));

            modulus = newModulus;
        }

        @Override
        public String toString() {
            return "LinearLift{" +
                    "poly=" + poly +
                    ", modulus=" + modulus +
                    ", aFactor=" + aFactor +
                    ", bFactor=" + bFactor +
                    ", aCoFactor=" + aCoFactor +
                    ", bCoFactor=" + bCoFactor +
                    '}';
        }
    }

}