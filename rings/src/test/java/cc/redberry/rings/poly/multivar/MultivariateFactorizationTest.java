package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.test.FactorizationInput;
import cc.redberry.rings.poly.univar.*;
import cc.redberry.rings.test.Benchmark;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.RandomDataGenerator;
import cc.redberry.rings.util.RandomUtil;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.text.SimpleDateFormat;
import java.util.*;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static cc.redberry.rings.Rings.*;
import static cc.redberry.rings.poly.PolynomialMethods.Factor;
import static cc.redberry.rings.poly.PolynomialMethods.polyPow;
import static cc.redberry.rings.poly.multivar.MultivariateFactorization.bivariateDenseFactorSquareFreeInGF;
import static cc.redberry.rings.poly.multivar.MultivariateFactorization.factorPrimitiveInGF;

/**
 * @since 1.0
 */
public class MultivariateFactorizationTest extends AMultivariateTest {
    @Ignore
    @Test
    public void testBivariate1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(67);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("b*a^5 + a^4*b^2 + 11 + b^3", domain),
                b = MultivariatePolynomialZp64.parse("a^6*b + 66*b + 17*b^2 + 1", domain),
                c = MultivariatePolynomialZp64.parse("b^3*a^4 + a^4 + b", domain),
                d = MultivariatePolynomialZp64.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);

        System.out.println(base);
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            Assert.assertEquals(4, bivariateDenseFactorSquareFreeInGF(base).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testBivariate2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(62653);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64[] factors = {
                MultivariatePolynomialZp64.parse("17096+6578*a*b^2+54905*a^3", domain, vars),
                MultivariatePolynomialZp64.parse("43370+32368*a^2*b^2+45712*a^2*b^4+52302*a^4+23776*a^4*b^2", domain, vars)
        };
        MultivariatePolynomialZp64 poly = factors[0].createOne().multiply(factors);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factorization = bivariateDenseFactorSquareFreeInGF(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }

    @Test
    public void testBivariateRandom3() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 5);
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(100, 1000),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"));
    }

    @Test
    public void testBivariateRandom4() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                5, 10,
                3, 6);
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(10, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"));
    }

    @Test
    public void testBivariate5() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1336151);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64[] factors = {
                MultivariatePolynomialZp64.parse("319792+402081*a^3", domain, vars),
                MultivariatePolynomialZp64.parse("685686+694157*a", domain, vars),
                MultivariatePolynomialZp64.parse("616781+1057293*b^2+158725*a+730076*a*b^2", domain, vars)
        };
        MultivariatePolynomialZp64 poly = factors[0].createOne().multiply(factors);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factorization = bivariateDenseFactorSquareFreeInGF(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }

    @Test
    public void testBivariate6() throws Exception {
        IntegersZp64 domain = new IntegersZp64(57352861);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64[] factors = {
                MultivariatePolynomialZp64.parse("15042434+15817122*b", domain, vars),
                MultivariatePolynomialZp64.parse("39330400+51579304*a^2", domain, vars)
        };
        MultivariatePolynomialZp64 poly = factors[0].createOne().multiply(factors);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factorization = bivariateDenseFactorSquareFreeInGF(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }


    @Test
    public void testBivaraiteSmallDomain7() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1 + b*a^5 + a^4*b^2 + b^3", domain, vars),
                b = MultivariatePolynomialZp64.parse("a + a^2 + b + a^6*b + 66*b + 17*b^2 + 1", domain, vars),
                c = MultivariatePolynomialZp64.parse("b^3*a^4 + a^4 + b", domain, vars),
                d = MultivariatePolynomialZp64.parse("a^5 + b^5*a^5 + b^2 + 3", domain, vars),
                base = a.clone().multiply(b, c, d);
        //System.out.println(base);
        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = bivariateDenseFactorSquareFreeInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(5, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }

        //49ms 35ms 30ms 28ms 38ms 38ms 28ms 28ms 33ms 33ms 33ms 34ms
    }

    @Test
    public void testBivaraiteSmallDomain5Random8() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 5);
        source.minModulusBits = 2;
        source.maxModulusBits = 3;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(50, 500),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization"));
    }

    @Test
    public void testBivaraiteSmallDomain5Random8a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(5);
        String[] vars = {"x", "y"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+x^3*y+x^6*y^4+2*x^6*y^5+2*x^6*y^6", domain, vars),
                b = MultivariatePolynomialZp64.parse("x^5+4*y^6+2*x^5*y^2", domain, vars),
                c = MultivariatePolynomialZp64.parse("1+x^2+4*x^3+3*x^6+4*x^3*y^4+x^4*y^4", domain, vars),
                d = MultivariatePolynomialZp64.parse("1+2*x^4*y^2+x^3*y^4+2*x^6*y^6", domain, vars),
                e = MultivariatePolynomialZp64.parse("1+3*x^4*y+3*x^3*y^4+4*x^4*y^5", domain, vars),
                base = a.clone().multiply(b, c, d, e);

        MultivariatePolynomial<BigInteger> bBase = base.toBigPoly();
        //System.out.println(base);
        for (int i = 0; i < its(40, 40); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> decomposition = bivariateDenseFactorSquareFreeInGF(bBase);
            System.out.println("" + i + "  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(5, decomposition.size());
            FactorDecompositionTest.assertFactorization(bBase, decomposition);
        }
    }

    @Test
    public void testBivaraiteSmallDomain5Random9() throws Exception {
        lSampleDecompositionSource source = new lSampleDecompositionSource(
                3, 5,
                2, 2,
                5, 10,
                3, 6);
        source.minModulusBits = 2;
        source.maxModulusBits = 3;
        testFactorizationAlgorithm(filterNonSquareFree(filterMonomialContent(source)), its(10, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization (small ring)"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInGF, "Bivariate dense factorization (small ring)"));
    }

    @Test
    public void testBivariateZ10() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("a^5 + a^4*b^2 + 11 + b^3"),
                b = MultivariatePolynomial.parse("a^6 + 66*b + 17*b^2 + 1"),
                c = MultivariatePolynomial.parse("a^4 + a + b*a^3 + 1"),
                d = MultivariatePolynomial.parse("a^5 + b^5*a^2 + b^2 + 3"),
                base = a.clone().multiply(b, c, d);

        //System.out.println(base);
        for (int i = 0; i < its(10, 10); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(4, factors.size());
        }
    }

    @Test
    public void testBivariateZ11() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("a^5*b + a^4*b^2 + 11 + b^3"),
                b = MultivariatePolynomial.parse("a^6 + 66*b + 17*b^2 + 1"),
                c = MultivariatePolynomial.parse("a^4 + a + b*a^3 + 1"),
                d = MultivariatePolynomial.parse("2*a^5*b^2 + a^5 + b^5*a^2 + b^2 + 3"),
                base = a.clone().multiply(b, c, d);

        //System.out.println(base);
        for (int i = 0; i < its(10, 10); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(4, factors.size());
        }
    }

    @Test
    public void testBivariateZRandom12() throws Exception {
        FactorizationInput.SampleDecompositionSource<MultivariatePolynomial<BigInteger>> source = new SampleDecompositionSourceZ(new lSampleDecompositionSource(
                3, 5,
                2, 2,
                2, 6,
                1, 5));

        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(100, 500),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInZ, "Bivariate dense factorization over Z"));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testBivariateZ12() throws Exception {
        MultivariatePolynomial<BigInteger> arr[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("177522327160552394341592645703385805765313887009291257319200617192383987428268805512690222419723419649 - 3920053743382812791836333767690842621932892411052334923146683840*a*b^2 + 7033850605913160576188071720524348507893675322801660689547349600868470384122764842150502780707785405529567188627563575595567387688045087798899184801064943616*a^5*b^2"),
                MultivariatePolynomial.parse("-1 - 659270292030815756363107128406733811547883294165085698816218454359750870703729555961546094888322347526082389105851760640*a^4 - 46529261163122475264752208961046337576273813308064326806338441809831450970874815332813704557963583488*a*b^2 + 60352646264726941172425539163284931317263872840808599407691511712170692052757127936*a^3*b^2 + 32892567946099972602587490607087056551535159861075199340856544477051037413463904644350103239197101283178014229071574120289916042759634944*a^5*b^2 + 1082364104922388503041167216*a^3*b^3"),
                MultivariatePolynomial.parse("1 - 396600905*a^4 + 18362289248902900227792935099132168196824077316766057448880233870344222135220897354236903227821776949155243835602763776*a^5*b + 54170215828526082015490079767653396191028300490271559531210492014654817736810561536*a*b^2 + 22296468805300143007110377483353074672219983699113622589773332810841273786269609285556407919332417536*a^3*b^3 - 4757394916849468388862409325495215636800300577101565824900427776*a^3*b^4"),
                MultivariatePolynomial.parse("-74245850854698227527322345985623966516292508703727252494450071089311844567647393096878536833433600*a + 25571940100073228627212382584833434076020757173045774321723224487468121499656011190161153995868899392979405069323465864668564085877571041780306601408943225575357784608735232*a^5*b^2 + 39224050727535633086898269296188590065607395610645363228145583907532677889236133260211280579276217051952112979635702373291087160595880341698628015404428693615610586950795264*b^5 - 74103934761117655866120723755998673350750638260236056304102159227089816262619525*a^4*b^5"),
                MultivariatePolynomial.parse("-1 - 106855402320422968799577513764486364196813269278986918073314897723185477854153087041104406250905236480*a^2*b^2 + 8068843832123230886122491909026167511493726107894055575979325820379758367061684714670315223007513134082333814492609438418585790965664643241520071083651235840*a^2*b^5 + 35672870377471658349867981049557362295743394300549931316635340036585883874243293966739476828669154066698149127319931379517583333329668421045053945301708161614217556217781288960*a^4*b^5"),},
                base = multiply(arr).multiply(-512);

        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(5, factors.size());
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testBivariateZ13() throws Exception {
        MultivariatePolynomial<BigInteger> arr[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1 + 1404971693217937738271243675250700150144919317697793951880016259765182277604375*a^3*b - 12133268343571183672190*a^2*b^2"),
                MultivariatePolynomial.parse("-1 - 312795970327582098045170081592831115925175*a^3 - 77950950035775481760560*a^3*b + 5401*a^3*b^2 + 40302*a^3*b^4"),
                MultivariatePolynomial.parse("1 - 5782055583736088058003458647322869334017577347037817545999736578862858345204522030604840232262500*a^6 - 43979083297211167074370*a^3*b^5 + 111820046899539828771390*a^6*b^6")},
                base = multiply(arr).multiply(-1);

        for (int i = 0; i < its(2, 12); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(3, factors.size());
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testBivariateZ14() throws Exception {
        MultivariatePolynomial<BigInteger> arr[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("(-3110505)*b-4683135*a-9030239*a^2"),
                MultivariatePolynomial.parse("1+7266251*a*b^5+11178392*a^2*b^2+2162182*a^2*b^6+5303702*a^6"),
                MultivariatePolynomial.parse("4789608+3904604*b+3917626*a*b+3416219*a^2*b")},
                base = multiply(arr);


        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            PrivateRandom.getRandom().setSeed(i);
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(3, factors.size());
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testBivariateZ15() throws Exception {
        MultivariatePolynomial<BigInteger> arr[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("(-93)-59*b+164*a^5+168*a^7*b+226*a^8*b"),
                MultivariatePolynomial.parse("1+245*b^3-211*a*b^2+110*a^5*b^2"),
                MultivariatePolynomial.parse("5-2*a*b-245*a^2+54*a^3*b^4-92*a^4*b^4")},
                base = multiply(arr);


        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            PrivateRandom.getRandom().setSeed(i);
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(3, factors.size());
        }
    }

    @Test
    public void testBivariateZRandom16() throws Exception {
        FactorizationInput.SampleDecompositionSource<MultivariatePolynomial<BigInteger>> source = new SampleDecompositionSourceZ(new lSampleDecompositionSource(
                3, 3,
                2, 2,
                3, 5,
                1, 6));
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(50, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::bivariateDenseFactorSquareFreeInZ, "Bivariate dense factorization over Z"));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testBivariateZ17() throws Exception {
        MultivariatePolynomial<BigInteger> arr[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1-679*a^5*b^5-552*a^5*b^7+675*a^7*b^5"),
                MultivariatePolynomial.parse("1-813*b-997*a^2*b^2-230*a^4*b^4+1294*a^5*b"),
                MultivariatePolynomial.parse("(-1)-8*a*b^3")},
                base = multiply(arr);


        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            PrivateRandom.getRandom().setSeed(i);
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.bivariateDenseFactorSquareFreeInZ(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            FactorDecompositionTest.assertFactorization(base, factors);
            Assert.assertEquals(3, factors.size());
        }
    }

    @Test
    public void testBivariateBenchmarkSingular() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        MultivariatePolynomialZp64 poly = MultivariatePolynomialZp64.parse("x^4120 + x^4118*y^2 + x^3708*y^400 + x^3706*y^402 + x^2781*y^1300 + x^2779*y^1302 + x^1339*y^2700 + x^927*y^3100 + y^4000 + x^7172*y^4167 + x^8349*y^4432 + x^8347*y^4434 + x^6760*y^4567 + x^5833*y^5467 + x^5568*y^7132 + x^11401*y^8599", domain);

        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            Assert.assertEquals(2, bivariateDenseFactorSquareFreeInGF(poly).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testBivariateBenchmarkSingular2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        MultivariatePolynomialZp64 poly = MultivariatePolynomialZp64.parse("b^1300+a^927*b^400+a^1339+a^5568*b^4432", domain);

        for (int i = 0; i < its(2, 2); i++) {
            long start = System.nanoTime();
            Assert.assertEquals(1, bivariateDenseFactorSquareFreeInGF(poly).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    static <Poly extends IPolynomial<Poly>>
    Poly multiply(Poly... factors) {
        return factors[0].createOne().multiply(factors);
    }

    @Test
    public void testGCDFreeBasis1() throws Exception {
        long modulus = 17;
        UnivariatePolynomialZp64
                a = UnivariatePolynomialZ64.create(1, 2, 3).modulus(modulus),
                b = UnivariatePolynomialZ64.create(3, 2, 1, 2).modulus(modulus),
                c = UnivariatePolynomialZ64.create(1, 0, 0, 2, 3).modulus(modulus),
                d = UnivariatePolynomialZ64.create(1, 11, 0, 12, 4).modulus(modulus);

        PolynomialFactorDecomposition<UnivariatePolynomialZp64>
                d1 = PolynomialFactorDecomposition.of(Arrays.asList(multiply(a, b, b, b), multiply(a, c, c), multiply(a, a, d))),
                d2 = PolynomialFactorDecomposition.of(Arrays.asList(multiply(a, c, b), multiply(b, d, c), multiply(d, c, d))),
                d3 = PolynomialFactorDecomposition.of(Arrays.asList(multiply(c, c, d), multiply(b, b, c), multiply(a, b, c, d)));

        @SuppressWarnings("unchecked")
        PolynomialFactorDecomposition<UnivariatePolynomialZp64>[] decomps = new PolynomialFactorDecomposition[]{
                d1.clone(), d2.clone(), d3.clone()
        };

        MultivariateFactorization.GCDFreeBasis(decomps);

        Assert.assertEquals(d1.multiply(), decomps[0].multiply());
        Assert.assertEquals(d2.multiply(), decomps[1].multiply());
        Assert.assertEquals(d3.multiply(), decomps[2].multiply());

        System.out.println(d1.multiply().equals(decomps[0].multiply()));
        System.out.println(Arrays.toString(UnivariateDivision.divideAndRemainder(d1.multiply(), decomps[0].multiply(), true)));
        for (PolynomialFactorDecomposition<UnivariatePolynomialZp64> decomp : decomps)
            System.out.println(decomp.size() + " => " + decomp);
    }

    @Test
    public void testGCDFreeBasis2() throws Exception {
        UnivariatePolynomial<BigInteger>
                a = UnivariatePolynomial.create(0, 1),
                b = UnivariatePolynomial.create(4, 8),
                c = UnivariatePolynomial.create(2);

        PolynomialFactorDecomposition<UnivariatePolynomial<BigInteger>>
                d1 = PolynomialFactorDecomposition.of(Arrays.asList(a)),
                d2 = PolynomialFactorDecomposition.of(Arrays.asList(b.primitivePart())).addFactor(b.contentAsPoly(), 1),
                d3 = PolynomialFactorDecomposition.of(Arrays.asList(c));

        @SuppressWarnings("unchecked")
        PolynomialFactorDecomposition<UnivariatePolynomialZp64>[] decomps = new PolynomialFactorDecomposition[]{
                d1.clone(), d2.clone(), d3.clone()
        };

        System.out.println(Arrays.asList(decomps));

        MultivariateFactorization.GCDFreeBasis(decomps);

        System.out.println(Arrays.asList(decomps));
    }

    @Test
    public void testGCDFreeBasis3() throws Exception {
        UnivariatePolynomialZp64
                a = UnivariatePolynomialZ64.create(9, 7, 10, 1).modulus(17),
                b = UnivariatePolynomialZ64.create(16, 2, 13, 1).modulus(17),
                c = UnivariatePolynomialZ64.create(3, 1).modulus(17),
                d = UnivariatePolynomialZ64.create(13, 14, 8, 1).modulus(17);

        PolynomialFactorDecomposition<UnivariatePolynomialZp64>
                d1 = PolynomialFactorDecomposition.of(a),
                d2 = PolynomialFactorDecomposition.of(b),
                d3 = PolynomialFactorDecomposition.of(c),
                d4 = PolynomialFactorDecomposition.of(d);

        @SuppressWarnings("unchecked")
        PolynomialFactorDecomposition<UnivariatePolynomialZp64>[] decomps = new PolynomialFactorDecomposition[]{
                d1.clone(), d2.clone(), d3.clone(), d4.clone(),
        };
//        System.out.println(Arrays.toString(Stream.of(a, b, c, d).map(p -> Factor(p).toString()).toArray()));
        System.out.println(Arrays.asList(decomps));
        MultivariateFactorization.GCDFreeBasis(decomps);
        System.out.println(Arrays.asList(decomps));

        Assert.assertEquals(d1.multiply(), decomps[0].multiply());
        Assert.assertEquals(d2.multiply(), decomps[1].multiply());
        Assert.assertEquals(d3.multiply(), decomps[2].multiply());
        Assert.assertEquals(d4.multiply(), decomps[3].multiply());
    }

    @Test
    public void testMultivariateFactorization1() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization1a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(42);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization1b() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("c*b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6*b + 66*b*c + 17*b^2 + c", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(328);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization2() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization3() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization4() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 100); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization4a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(31);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization4b() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(88);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization4c() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a^5*c + 2*a^5*b^2*c + 3*a^5*c^2 + a^5*b^2*c^2 + a^5*c^3 + 11*c + b^3 + a^5", domain),
                b = MultivariatePolynomialZp64.parse("a^6 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(531);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization5() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1361);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("2*a^5*c^2*b^2 + 11*c + b^3 + 1", domain),
                b = MultivariatePolynomialZp64.parse("a^6*b*c^3 + 66*b*c + 17*b^2 + c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization6() throws Exception {
        IntegersZp64 domain = new IntegersZp64(27239);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+10435*c^2+21950*a*b^2*c^2+17887*a*b^3*c+4648*a^2*c+862*a^2*b*c", domain),
                b = MultivariatePolynomialZp64.parse("1+21170*b^2*c+7162*b^3+18183*a^2*b^2*c+16794*a^3*b+3096*a^3*b^3*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization7() throws Exception {
        IntegersZp64 domain = new IntegersZp64(63185123);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("7962775*c^3+54287330*b^3+48565396*a^2+26248711*a^3*b^3+10971203*a^3*b^3*c", domain),
                b = MultivariatePolynomialZp64.parse("1+48442198*b^2+36965231*b^3+35212338*a*b^2*c^3+62918195*a^2*b^2*c+47759030*a^3*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization8() throws Exception {
        IntegersZp64 domain = new IntegersZp64(829657);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+81633*a+270565*a^2*b*c+799187*a^2*b^2+159093*a^3*b+562717*a^3*b^3*c", domain),
                b = MultivariatePolynomialZp64.parse("1+73615*a*b^2+92694*a*b^2*c^3+582676*a^3*b*c^3+144867*a^3*b^2*c^2+132332*a^3*b^2*c^3", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization9() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1734917);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+1179031*b^3+548360*a*b*c^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b", domain),
                b = MultivariatePolynomialZp64.parse("439433*b*c+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b", domain),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Benchmark
    @Test
    public void testMultivariateFactorization10() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+1179031*b^3+548360*a*b*c^2*e^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b*d^3+a^15", domain, vars),
                        MultivariatePolynomialZp64.parse("439433*b*c*d*e+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b+a^15", domain, vars),
                        MultivariatePolynomialZp64.parse("439433*d*c+197065*d*e+264505*a*c^3*d+1075508*a*d^3*e^3+1338483*a^15*e +a^15 + b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("433*d^2*c+165*d*e+265*a^4*c^3*d+107*a*d^3*b+1338*a^15*e +a^15*b +a^15 + b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("433*d^2*e+165*d*e+265*b^4*c^3*d+107*c*d^3*a+1338*a^15*e +a^15*e + c^2 + a^15*b + a^15 + a^15*d + a^15*e", domain, vars),
                },
                base = factors[0].createOne().multiply(factors);
//        System.out.println(PolynomialFactorDecomposition.create(Arrays.asList(factors)));
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(100, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base, false);
            Assert.assertEquals(5, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
        //1391 ms
        //1529 ms
        //1524 ms
        //1351 ms

        // update 14/11/17
        //705ms
        //644ms
        //681ms
        //688ms
    }

    @Test
    public void testMultivariateFactorization10a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+1179031*b^3+548360*a*b*c^2*e^2+18887*a*b^3*c+295179*a*b^3*c^3+175792*a^2*b*d^3+a^15", domain, vars),
                        MultivariatePolynomialZp64.parse("439433*b*c*d*e+197065*a+264505*a*b^3*c+1075508*a*b^3*c^3+1338483*a^2*b+a^15", domain, vars),
                        MultivariatePolynomialZp64.parse("439433*d*c+197065*d*e+264505*a*c^3*d+1075508*a*d^3*e^3+1338483*a^15*e +a^15 + b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("433*d^2*c+165*d*e+265*a^4*c^3*d+107*a*d^3*b+1338*a^15*e +a^15*b +a^15 + b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("433*d^2*e+165*d*e+265*b^4*c^3*d+107*c*d^3*a+1338*a^15*e +a^15*e + c^2 + a^15*b + a^15 + a^15*d + a^15*e", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(5, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization11() throws Exception {
        IntegersZp64 domain = new IntegersZp64(1734917);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("a + b + c + d", domain, vars),
                b = MultivariatePolynomialZp64.parse("b + c + d + e", domain, vars),
                c = MultivariatePolynomialZp64.parse("c + d + e + a", domain, vars),
                d = MultivariatePolynomialZp64.parse("2*d + 2*e + 2*a + 2*b", domain, vars),
                e = MultivariatePolynomialZp64.parse("e + a + b + 2*c", domain, vars),
                base = a.clone().multiply(b, c, d, e);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorToPrimitive(base);
        Assert.assertEquals(5, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization12() throws Exception {
        IntegersZp64 domain = new IntegersZp64(74017943);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("18370804*b^3*c^3+40724543*a+25831118*a*b+28120978*a*b*c^3+49314822*a*b^3*c^3", domain),
                b = MultivariatePolynomialZp64.parse("58629076*b*c^3+37897966*a*b^3*c^2+55834047*a^2*c^3+18265939*a^3*c^3+43535405*a^3*b", domain);

        a = AMultivariatePolynomial.renameVariables(a, new int[]{1, 2, 0});
        b = AMultivariatePolynomial.renameVariables(b, new int[]{1, 2, 0});
        MultivariatePolynomialZp64 base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);

        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        System.out.println(factors);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization13() throws Exception {
        IntegersZp64 domain = new IntegersZp64(192149);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+81770*b+19081*b*c+23953*a*b^3*c+7807*a^3*b+14026*a^3*b^2*c^3", domain),
                b = MultivariatePolynomialZp64.parse("1+105163*a^2+81015*a^2*c+166076*a^3*c^3+106464*a^3*b^2*c^2+43621*a^3*b^3", domain);
        MultivariatePolynomialZp64 base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization14() throws Exception {

        IntegersZp64 domain = new IntegersZp64(386039);
        MultivariatePolynomialZp64
                a = MultivariatePolynomialZp64.parse("1+377446*b*c+302126*a*b^2+97219*a*b^2*c^2+92497*a*b^2*c^3+84001*a^3*b^2", domain),
                b = MultivariatePolynomialZp64.parse("1+248663*a*b*c^2+10589*a*b^3*c^3+62097*a^2*c+81842*a^2*b^2*c+51504*a^3", domain);
        a = AMultivariatePolynomial.renameVariables(a, new int[]{0, 2, 1});
        b = AMultivariatePolynomial.renameVariables(b, new int[]{0, 2, 1});
        MultivariatePolynomialZp64 base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> factors = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(2, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    public void testMultivariateFactorization15() throws Exception {
        IntegersZp64 domain = new IntegersZp64(34957081);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("2999550*b*c^3*d^2+14700809*a*c^2+13282494*a*b^3*c*d^3+30075047*a^2*c^2*d+2736476*a^3*d^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+7381919*c^2+33667094*b*c^2*d^2+26355114*b^3*c*d^3+30536438*b^3*c^2*d^3+10734561*a*b*c", domain, vars),
                        MultivariatePolynomialZp64.parse("1+24559556*b^2*c^2*d^3+13753085*a*b*c^2*d+22081133*a*b*c^3*d^3+30781594*a^3*b*c^2+27334226*a^3*b^3*d^3", domain, vars),
                };

        for (int i = 0; i < factors.length; i++) {
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 0, 1);
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 1, 3);
        }

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(3, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization16() throws Exception {
        IntegersZp64 domain = new IntegersZp64(316797977);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("157129769*b*c*d^3+234842760*b*c^3*d^3+105538252*a+55274980*a*b^2*c^2+89647854*a^3*b^3*c^2*d^2", domain, vars),
                        MultivariatePolynomialZp64.parse("241626121*d^2+47627151*b^2*c*d+150262012*a^2*b*c^2+299159387*a^2*b^3*c^2+53788517*a^3*b*c^3*d^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+76011411*b*d^3+189305430*b*c^3*d+218732499*a^2*b^2*c*d^2+125990992*a^2*b^3*c^2*d+36953173*a^3*b*c^2*d", domain, vars),
                        MultivariatePolynomialZp64.parse("1+299415864*a*b^2*c^3*d^3+154985851*a^2*c*d+157246866*a^2*b^2*c^3*d^3+32838497*a^3*b^3*c*d+41239905*a^3*b^3*c*d^2", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base, false);
        Assert.assertEquals(4, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization17() throws Exception {
        IntegersZp64 domain = new IntegersZp64(125617);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+71543*b*c+89032*b*c*d+101233*b^2*c+69912*b^2*c^2*d^2+122146*a*c^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+62395*b*c^2*d^2+111331*a*b*c*d^2+13129*a^3*b^2*c^2*d^2+54277*a^3*b^3*c^2+36488*a^3*b^3*c^2*d^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+64768*a*b^2*d+66817*a*b^3*c*d+19563*a^2+13861*a^3*b*c^3+76958*a^3*b^3*d", domain, vars),
                };

        for (int i = 0; i < factors.length; i++)
            factors[i] = AMultivariatePolynomial.swapVariables(factors[i], 0, 1);

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        for (int i = 0; i < its(10, 10); i++) {
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base, false);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorizationRandom1() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 2,
                3, 3,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(100, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"));
    }

    @Test
    public void testMultivariateFactorizationRandom2() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 4,
                3, 4,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(100, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"));
    }

    @Test
    public void testMultivariateFactorizationRandom3() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 4,
                5, 6,
                5, 5,
                3, 3);
        source.minModulusBits = 15;
        source.maxModulusBits = 30;
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(100, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields"));
    }

    @Test
    public void testMultivariateFactorizationRandom3a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(337965077);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+128786660*b*c*f^2+38739797*a*b^3*c^2*d*e*f+159449306*a^2*c^2*d^2*e^3*f^3+298952491*a^2*b*c^3*d*e*f^3+263798205*a^3*c^2*e*f^3", domain, vars),
                        MultivariatePolynomialZp64.parse("69412172*c*f+175964784*c^3*f+319203880*a^3*c*d^3*f+154158650*a^3*b^2*c^3*d^3*f^2+309716243*a^3*b^3*d^3*e^3", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        for (int i = 0; i < its(10, 10); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(2, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }

        //7s
        //6s
        //7s
        //8s
        //5s
        //10s
        //6s
    }


    @Test
    public void testMultivariateFactorization18_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+71543*b*c+89032*b*c*d+101233*b^2*c+69912*b^2*c^2*d^2+122146*a*c^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+62395*b*c^2*d^2+111331*a*b*c*d^2+13129*a^3*b^2*c^2*d^2+54277*a^3*b^3*c^2+36488*a^3*b^3*c^2*d^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+64768*a*b^2*d+66817*a*b^3*c*d+19563*a^2+13861*a^3*b*c^3+76958*a^3*b^3*d", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        for (int i = 0; i < its(10, 30); i++) {
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization19_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+2*b*c^3*d*f^3+2*b^2*c^3*d^2*e^2*f^3+2*a*b^3*d*e^2*f^2+2*a^2*b*c^3*d^3*e*f^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+2*a^2*b^2*d*e^2*f^2+a^3*b^3*c*d^3*f", domain, vars),
                        MultivariatePolynomialZp64.parse("1+a*b^2*c^2*e*f^2+2*a^3*e^2+2*a^3*b*c*d^2*e^2+2*a^3*b^2*c*e", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*c*d^3*e^3+a*b*c^2*d*e^3*f^3+2*a^2*e^2*f", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(0);
        for (int i = 1; i < its(10, 10); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization20_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("b*d^3*e^2+a*c^3*d^2+a*b^3*c^3*e^3+a^3*b^2*d^3+a^3*b^3*c^2*d*e^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^2*c^3*d^3*e^2+a^3*c*e^2+a^3*b^2*c^3*d", domain, vars),
                MultivariatePolynomialZp64.parse("1+a*b*c*e^2+a^3*b^3*c*d^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c*d^2*e^2+a^3*c*e", domain, vars),},
                base = multiply(arr);
        //System.out.println(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization21_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(11);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+5*a*c*e+6*a^2*c*d+9*a^2*b^3*c^3*e+2*a^3*c^3*d^3*e^2+9*a^3*b^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*b^2*d*e+7*a*c^2*d*e^2+5*a^2*d*e^3+a^2*b^2*d^3*e^3+8*a^3*c^2*d*e^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+9*b^2*d*e+10*b^2*d^3*e+8*b^2*c^2*d+3*a*b^2*c^2*d^3+8*a*b^3*c^2*d^3*e^2", domain, vars),
        },
                base = multiply(arr);
        //System.out.println(base);
        for (int i = 11; i < its(12, 12); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization22_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(11);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("12*b*c^2*d+3*a*f+8*a^2*b*d^3*f^3+13*a^2*b*c*d*e*f^3+16*a^3*c^3*d^2*e^3", domain, vars),
                MultivariatePolynomialZp64.parse("16*c*d^2*f+9*c*d^2*e*f^2+17*c*d^3*e^3+14*a^3*b^3*d*f^2+8*a^3*b^3*c^2*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+8*a*b^2*c*d*f+18*a*b^3*d^3*e^2*f+22*a^2*c*d*e*f+a^2*b^2*c*d^2*e*f^2+12*a^3*b^3*d^2*f^2", domain, vars),

        },
                base = multiply(arr);
        //System.out.println(base);
        for (int i = 0; i < its(2, 2); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization23_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+a^2*b^2*c*e^3*f^2+a^3*b^2*d^2*e^2*f^2+a^3*b^2*c^3*e*f^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+a^2*b*c*d^3*f+a^3*c^3*d*e^2*f+a^3*b^3*c^3*d*e^2*f^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^2*c^2*d^3+a^3*b*c^2*e^3*f^3", domain, vars),

        },
                base = multiply(arr);
        //System.out.println(base);
        for (int i = 11; i < its(12, 12); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization24_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(29);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+8*b+11*b^2*c*d^3+4*a*b^2*c^3*d+23*a^3*b*c^2*d^2*e^2+25*a^3*b^2*d^2*e", domain, vars),
                MultivariatePolynomialZp64.parse("24*b^2*c^2*e^2+22*a*d*e^2+18*a*b*c^3*d^3+22*a^2*b*c*d^2*e^3+19*a^2*b^3*c*d^2*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+9*a*b*c^3*e+20*a*b^3*c*d^3*e^2+25*a^3*b*c*d+17*a^3*b^2*c^3+16*a^3*b^3*c*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+11*b*c^2*d*e^2+18*a*c^3*d^3*e^3+8*a*b*d^2*e^3+11*a*b^3*c^2*d^2*e+5*a^3*b*c^3*d^3*e^2", domain, vars),
        }, base = multiply(arr);
        System.out.println(MultivariateSquareFreeFactorization.isSquareFree(base));
        System.out.println(MultivariateFactorization.factorToPrimitive(base).size());
        //System.out.println(base);
        for (int i = 0; i < its(12, 12); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization25_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+b*c*d^3*e^2*f^2+b^2*e+a^2*b*d^2*f", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c^3*d*e^2+a^2*b^3*c*e^3*f^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c*d^2+a*b*c^2*d*e^2+a^2*b^3*c^2*f^3+a^3*b^2*c^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+a*b*c*d*e*f^3+a*b^3*c*d*e^3*f^2+a^2*c*d^2*e*f^3+a^3*b^3*d^3*e", domain, vars),
        }, base = multiply(arr);
        for (int i = 0; i < its(1, 1); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }

        //80s
    }

    @Test
    public void testMultivariateFactorization26_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(17);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("15*b*c^3*d^3*e^2+b^2*d^2*e^3+15*a*b^2*c^3*e+3*a^2+16*a^2*c*d^3*e^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^2*c*d^3*e+3*a*b^2+8*a^3*b^2*c*d^2*e^2+2*a^3*b^2*c^2*d*e^2+12*a^3*b^2*c^3*d^2*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+14*b^3*d^3+10*a*b*c^3*d^2+12*a^2*b^2*d^2*e^2+6*a^2*b^3*c^3*d*e^3+11*a^3*b^2*c^3*d", domain, vars),
                MultivariatePolynomialZp64.parse("15*c^3*e^3+a^2*b^3*d^3+9*a^2*b^3*c^3+15*a^3*b*d^3", domain, vars),
        }, base = multiply(arr);
        for (int i = 0; i < its(1, 1); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization27_SmallDomain() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+b^3*c^3*d^2+a^2*b*c^3*e*f^2+2*a^3*b^2*c*d*e^3*f", domain, vars),
                MultivariatePolynomialZp64.parse("1+b*c^3*d*e*f^2+b^2*c*e+2*b^3*c^3*d^2+a^3*b*c*d^3*e*f^3", domain, vars),
        }, base = multiply(arr);
        for (int i = 0; i < its(1, 1); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(2, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorizationRandom4_SmallDomain() throws Exception {
        MultivariateFactorizationTest.lSampleDecompositionSource source = new MultivariateFactorizationTest.lSampleDecompositionSource(
                2, 4,
                5, 6,
                5, 5,
                3, 3);
        source.minModulusBits = 2;
        source.maxModulusBits = 5;
        //DETALIZATION_PERCENT = 1000;
        //PRINT_FACTORS = true;
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(100, 100),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields (small ring)"),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInGF, "Multivariate factorization in finite fields (small ring)"));
    }

    @Test
    public void testMultivariateFactorizationRandom4_SmallDomain_a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(11);
        String[] vars = {"a", "b", "c", "d", "e", "f"};
        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+10*b^2*c^3*d^3*e^2*f^2+10*b^3*d^2*e*f^3+7*a*b*c^2*d*e^2*f^3+9*a*b^3*e^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*c^3*d^2*e*f^2+10*c^3*d^2*e^3+6*b*d^2*e^3*f^3+9*a^3*b*c^3*d*e*f", domain, vars),
        }, base = multiply(arr);
        for (int i = 43; i < its(50, 50); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition decomposition = MultivariateFactorization.FactorInGF(base.toBigPoly());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(2, decomposition.size());
            FactorDecompositionTest.assertFactorization(base.toBigPoly(), decomposition);
        }
    }

    @Test(timeout = 200000)
    public void testMultivariateFactorizationRandom4_SmallDomain_b() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};

        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+2*c^3*d^2+2*b^3*c^3*d^3*e+a*c^3*d*e+2*a^2*b^3*c^2*d^2*e^3+a^2*b^3*c^3*e^2", domain, vars),
                MultivariatePolynomialZp64.parse("1+b^3*c^2*d^3*e^3+a*c^3*d*e^2+2*a^3*e^3+2*a^3*b^3*d*e^3+2*a^3*b^3*c*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*a*b^3*c+a^2*d^3*e", domain, vars),
                MultivariatePolynomialZp64.parse("1+2*b^3*c^3*d^3*e+2*a*b^2*c*d^2*e^3+a*b^3*c^2*d*e^2+a^3*b^2*c^3*d^2", domain, vars),
        }, base = multiply(arr);
        System.out.println(base.sparsity());
        for (int i = 0; i < its(5, 5); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i + 1);
            timestamp();
            PolynomialFactorDecomposition decomposition = MultivariateFactorization.FactorInGF(base);
            timeElapsed();
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
        //12s
        //3s
        //15s
        //18s
    }

    @Test
    public void testMultivariateFactorizationRandom4_SmallDomain_c() throws Exception {
        IntegersZp64 domain = new IntegersZp64(5);
        String[] vars = {"a", "b", "c", "d", "e", "f"};

        MultivariatePolynomialZp64 arr[] = {
                MultivariatePolynomialZp64.parse("1+2*c^3*d^2*e*f^2+3*a*b*c*e^3*f^2+4*a^2*d*e+4*a^3*c^3*d^2*e^3", domain, vars),
                MultivariatePolynomialZp64.parse("1+4*b^2*c^2*e^3*f^2+b^3*c^2*d*e^3+a*b^3*c*d^3*f^3+a*b^3*c^2*e*f^2+3*a*b^3*c^2*e^3*f", domain, vars),
        }, base = multiply(arr);
        for (int i = 0; i < its(5, 5); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i + 1);
            timestamp();
            PolynomialFactorDecomposition decomposition = MultivariateFactorization.FactorInGF(base);
            timeElapsed();
            Assert.assertEquals(2, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
    }

    @Test
    public void testMultivariateFactorization18() throws Exception {
        IntegersZp64 domain = new IntegersZp64(11);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+3*c^3*d^2+9*a*c*e^3+a^2*b*c^3*d^2*e^2+2*a^2*b^2*c^3*d*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+5*a*b*c*d^3*e^3+6*a^3*d^2*e^2+5*a^3*b*c^2*d*e+9*a^3*b^2*c*d^2*e^2+2*a^3*b^3*c*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+8*b^3*c^2*d+2*a^2*c*e^3+6*a^2*c^2*d^3*e+10*a^3*c^3*d^2+7*a^3*b*c^3*d^3*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+4*a*b^2*c^2*d^3*e^3+4*a^2*b*c*d^2+3*a^3*b*c^3*e^3+3*a^3*b^2*c^2*d^3*e^2", domain, vars),
                };

        //System.out.println(PolynomialFactorDecomposition.create(Arrays.asList(factors)));

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            //PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(4, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }
        // 3s
        // 1529ms
        // 1494ms
        // 1465ms
        // 622ms
        // 1464ms
        // 1243ms
        // 1107ms
        // 1092ms
        // 1037ms

        // update 14/11/17

        //204ms
        //609ms
        //142ms
        //175ms
        //541ms
    }

    @Test
    public void testMultivariateFactorization18a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(11);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+3*c^3*d^2+9*a*c*e^3+a^2*b*c^3*d^2*e^2+2*a^2*b^2*c^3*d*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+5*a*b*c*d^3*e^3+6*a^3*d^2*e^2+5*a^3*b*c^2*d*e+9*a^3*b^2*c*d^2*e^2+2*a^3*b^3*c*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+8*b^3*c^2*d+2*a^2*c*e^3+6*a^2*c^2*d^3*e+10*a^3*c^3*d^2+7*a^3*b*c^3*d^3*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+4*a*b^2*c^2*d^3*e^3+4*a^2*b*c*d^2+3*a^3*b*c^3*e^3+3*a^3*b^2*c^2*d^3*e^2", domain, vars),
                };

        //System.out.println(PolynomialFactorDecomposition.create(Arrays.asList(factors)));

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(9);

        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
        Assert.assertEquals(4, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }


    @Test
    public void testMultivariateFactorization19() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d^2*e+a*b^3*c^2*d^2+a^3*b^2*c*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*d^3*e^2+a*d^3*e^3+a^2*b*c^3*d^2*e^3+a^3*b^3*c^2*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c*d^3+b*c^2*d^2*e^2+a^3*b^2*c*d^3*e^2", domain, vars),
                };

        //System.out.println(PolynomialFactorDecomposition.create(Arrays.asList(factors)));

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);

        // sho
        for (int i = 0; i < its(10, 10); i++) {
            //cc.r2.core.poly.univar.PrivateRandom.getRandom().setSeed(i);
            //PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(3, decomposition.size());
            FactorDecompositionTest.assertFactorization(base, decomposition);
        }

        // Timings:

        // 127ms
        // 124ms
        // 126ms
        // 244ms
        // 88ms
        // 64ms
        // 152ms
        // 174ms
        // 285ms
    }

    @Test
    public void testMultivariateFactorization19a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d^2*e+a*b^3*c^2*d^2+a^3*b^2*c*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*d^3*e^2+a*d^3*e^3+a^2*b*c^3*d^2*e^3+a^3*b^3*c^2*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c*d^3+b*c^2*d^2*e^2+a^3*b^2*c*d^3*e^2", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(5);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
        Assert.assertEquals(3, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization19b() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d^2*e+a*b^3*c^2*d^2+a^3*b^2*c*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*d^3*e^2+a*d^3*e^3+a^2*b*c^3*d^2*e^3+a^3*b^3*c^2*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c*d^3+b*c^2*d^2*e^2+a^3*b^2*c*d^3*e^2", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(2);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
        Assert.assertEquals(3, decomposition.size());
        FactorDecompositionTest.assertFactorization(base, decomposition);
    }

    @Test
    public void testMultivariateFactorization19c() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d^2*e+a*b^3*c^2*d^2+a^3*b^2*c*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*d^3*e^2+a*d^3*e^3+a^2*b*c^3*d^2*e^3+a^3*b^3*c^2*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c*d^3+b*c^2*d^2*e^2+a^3*b^2*c*d^3*e^2", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(7);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
        Assert.assertEquals(3, decomposition.size());
    }

    @Test//(timeout = 10000L)
    public void testMultivariateFactorization19d() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d^2*e+a*b^3*c^2*d^2+a^3*b^2*c*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^3*d^3*e^2+a*d^3*e^3+a^2*b*c^3*d^2*e^3+a^3*b^3*c^2*d^2*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c*d^3+b*c^2*d^2*e^2+a^3*b^2*c*d^3*e^2", domain, vars),
                };

        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        //cc.r2.core.poly.univar.PrivateRandom.getRandom().setSeed(10);
        PrivateRandom.getRandom().setSeed(10);
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
        Assert.assertEquals(3, decomposition.size());
    }

    @Test
    public void testMultivariateFactorization20() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+b^3*e^2+2*a*b*c^3*d*e^2+a^2*b^2*c^2*e+a^3*b^3*c^3*d^3*e", domain, vars),
                        MultivariatePolynomialZp64.parse("1+2*b^3*c*e^2+2*b^3*c^2*d*e^2+a^2*b*d^2*e^3+2*a^3*b^2*c^2*d^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+2*b*c^3*d^2*e+a^2*b*c*e^3+a^3*b*c^3*d^3+a^3*b^3*c^3*d^2*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+c^2*e+b^3*c^2*d^2*e+2*a*b^3*c^3*d^2*e+2*a^2*c^2*d*e^3+a^3*c^3*d^2", domain, vars),
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(2, 2); i++) {
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            Assert.assertEquals(4, decomposition.size());
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization21() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+2*a*c*d^2+a^2*b*c^2*d^3+2*a^2*b^3*d*e", domain, vars),
                        MultivariatePolynomialZp64.parse("b^3*c*e+2*a*d^3+2*a^3*b^3*c*e^2", domain, vars),
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(20, 20); i++) {
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            Assert.assertEquals(2, decomposition.size());
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization22() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^2*c^2*d+a*b^3*c^3+a^2*b^2*c*d^2+2*a^2*b^3*c*d^2*e", domain, vars),
                        MultivariatePolynomialZp64.parse("c^2*d+2*a^3*b^3*e", domain, vars),
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        //System.out.println(factorToPrimitive(base));
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(20, 20); i++) {
            //long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.factorPrimitiveInGF(base);
            Assert.assertEquals(2, decomposition.size());
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization23() throws Exception {
        IntegersZp64 domain = new IntegersZp64(3);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+2*a*b*c*e^2+a^2*b^2*c^3*e+a^3*d^3*e^3+a^3*b^2*c^2*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+d*e^3+a^3*b*e^2+a^3*b*d^3*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b^2*c^2*d*e+2*b^2*c^3*d*e^2+a^3*c^3*d*e^2+a^3*b^3*c^3*d*e^3", domain, vars),
                        MultivariatePolynomialZp64.parse("2*b^3*c*d^3*e+2*a*b^2*c^2*e^3+2*a^2*c^3*e+2*a^2*b*c^2*d^3*e^2+a^3*b^3*d^2", domain, vars),
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(20, 20); i++) {
            PrivateRandom.getRandom().setSeed(3);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            Assert.assertEquals(4, decomposition.size());
            System.out.println("=> " + i + "  " + TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization24() throws Exception {
        IntegersZp64 domain = new IntegersZp64(2);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("1+a*b^3*d+a^3*c*d*e^3+a^3*b*c^2*e^3+a^3*b^2*c^3*d*e", domain, vars),
                        MultivariatePolynomialZp64.parse("1+a^2*b*c*d*e^3+a^3*b^2*c*d", domain, vars),
                        MultivariatePolynomialZp64.parse("1+b*c^2*e^2+a^3*b*c^2*e+a^3*b^3*c^3*d*e^2", domain, vars)
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(20, 20); i++) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            Assert.assertEquals(3, decomposition.size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }

        // 577ms
        // 534ms
        // 177ms
        // 220ms
        // 414ms
        // 365ms
        // 202ms
        // 194ms
        // 322ms
    }

    @Ignore
    @Test
    public void testMultivariateFactorization26() throws Exception {
        IntegersZp64 domain = new IntegersZp64(7);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("b*c*d*e+5*a^2*c*d*e^2+5*a^3*b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+2*a*d^2+4*a*b*c^2*d*e+5*a^2*b^2*c*d*e^3+a^2*b^2*c^2*e^2+3*a^2*b^2*c^3*d*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("3*c^3*d+5*a*b*d^3*e+2*a^2*b^3*c*d+3*a^3*d*e^2+3*a^3*b^2*c^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+6*a*b*d^2+6*a*b^2*c^2*d^2*e^2+a^2*b^2*c^2*d^3+4*a^3*c^3*d^2", domain, vars)
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(20, 20); i++) {
            System.out.println(i + 20);
            PrivateRandom.getRandom().setSeed(i + 20);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
            Assert.assertEquals(4, decomposition.size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization26a() throws Exception {
        IntegersZp64 domain = new IntegersZp64(7);
        String[] vars = {"a", "b", "c", "d", "e"};
        MultivariatePolynomialZp64
                factors[] =
                {
                        MultivariatePolynomialZp64.parse("b*c*d*e+5*a^2*c*d*e^2+5*a^3*b^2", domain, vars),
                        MultivariatePolynomialZp64.parse("1+2*a*d^2+4*a*b*c^2*d*e+5*a^2*b^2*c*d*e^3+a^2*b^2*c^2*e^2+3*a^2*b^2*c^3*d*e^2", domain, vars),
                        MultivariatePolynomialZp64.parse("3*c^3*d+5*a*b*d^3*e+2*a^2*b^3*c*d+3*a^3*d*e^2+3*a^3*b^2*c^3", domain, vars),
                        MultivariatePolynomialZp64.parse("1+6*a*b*d^2+6*a*b^2*c^2*d^2*e^2+a^2*b^2*c^2*d^3+4*a^3*c^3*d^2", domain, vars)
                };
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(27);
        long start = System.nanoTime();
        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
        Assert.assertEquals(4, decomposition.size());
        System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
    }

    @Ignore
    @Test
    public void testMultivariateFactorization25() throws Exception {
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
        MultivariatePolynomialZp64 base = factors[0].createOne().multiply(factors);

        PolynomialFactorDecomposition<MultivariatePolynomialZp64> decomposition = MultivariateFactorization.FactorInGF(base);
        Assert.assertEquals(3, decomposition.size());
    }

    @Test
    public void testMultivariateFactorization27() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("a^5 + a^4*b^2*c^2 + 11*c + b^3 + 2"),
                b = MultivariatePolynomial.parse("2*a^6 - 66*b*c + 17*b^2 + c"),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization28() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + 2"),
                b = MultivariatePolynomial.parse("2*a^6 - 66*b*c + 17*b^2 + c"),
                base = a.clone().multiply(b);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(2, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization29() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("b*a^5 + a^4*b^2*c^2 + 11*c + b^3 + 2"),
                b = MultivariatePolynomial.parse("2*a^6 - 66*b*c + 17*b^2 + c"),
                c = MultivariatePolynomial.parse("2*a^6 + 4*b^2*a^6 - 66*b*c + 17*b^2 + b + 1", "a", "b", "c"),
                base = a.clone().multiply(b, c);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(3, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization30() throws Exception {
        MultivariatePolynomial<BigInteger>
                a = MultivariatePolynomial.parse("1+28136331*b*c^2-19241375*a*b^2*c+31673065*a*b^2*c^2-21766903*a^2*b^2*c^2"),
                b = MultivariatePolynomial.parse("1-13098013*a^2*b^2*c-32210690*a^2*b^3*c^2+20904590*a^3*b^2*c^3"),
                c = MultivariatePolynomial.parse("1+18527870*a*b-17012638*a*b*c+10557740*a^2*c-9843426*a^2*b^2"),
                d = MultivariatePolynomial.parse("(-4727463)-17684585*b+11364285*a*b^2+15438263*a^2*c^2-14902664*a^2*b*c^2"),
                e = MultivariatePolynomial.parse("(-4994407)+32080518*a*b^2*c"),
                base = a.clone().multiply(b, c, d, e);
        //System.out.println(base);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(5, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization30a() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+28136331*a*c^2-19241375*a^2*b*c+31673065*a^2*b*c^2-21766903*a^2*b^2*c^2"),
                MultivariatePolynomial.parse("1-13098013*a^2*b^2*c+20904590*a^2*b^3*c^3-32210690*a^3*b^2*c^2"),
                MultivariatePolynomial.parse("1+10557740*b^2*c+18527870*a*b-17012638*a*b*c-9843426*a^2*b^2"),
                MultivariatePolynomial.parse("(-4727463)+15438263*b^2*c^2-17684585*a-14902664*a*b^2*c^2+11364285*a^2*b"),
                MultivariatePolynomial.parse("(-4994407)+32080518*a^2*b*c"),
        };

        MultivariatePolynomial<BigInteger> base = multiply(pp);

        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(1);
        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
        Assert.assertEquals(5, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization30b() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+28136331*a*c^2-19241375*a^2*b*c+31673065*a^2*b*c^2-21766903*a^2*b^2*c^2"),
                MultivariatePolynomial.parse("1-13098013*a^2*b^2*c+20904590*a^2*b^3*c^3-32210690*a^3*b^2*c^2"),
                MultivariatePolynomial.parse("1+10557740*b^2*c+18527870*a*b-17012638*a*b*c-9843426*a^2*b^2"),
                MultivariatePolynomial.parse("(-4727463)+15438263*b^2*c^2-17684585*a-14902664*a*b^2*c^2+11364285*a^2*b"),
                MultivariatePolynomial.parse("(-4994407)+32080518*a^2*b*c"),
        };

        MultivariatePolynomial<BigInteger> base = multiply(pp);

        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(3);
        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
        Assert.assertEquals(5, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization30c() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+28136331*a*c^2-19241375*a^2*b*c+31673065*a^2*b*c^2-21766903*a^2*b^2*c^2"),
                MultivariatePolynomial.parse("1-13098013*a^2*b^2*c+20904590*a^2*b^3*c^3-32210690*a^3*b^2*c^2"),
                MultivariatePolynomial.parse("1+10557740*b^2*c+18527870*a*b-17012638*a*b*c-9843426*a^2*b^2"),
                MultivariatePolynomial.parse("(-4727463)+15438263*b^2*c^2-17684585*a-14902664*a*b^2*c^2+11364285*a^2*b"),
                MultivariatePolynomial.parse("(-4994407)+32080518*a^2*b*c"),
        };

        MultivariatePolynomial<BigInteger> base = multiply(pp);

        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        PrivateRandom.getRandom().setSeed(9);
        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
        Assert.assertEquals(5, factors.size());
        FactorDecompositionTest.assertFactorization(base, factors);
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization31() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+47250*b^5+9661*a^2*b-118880*a^5*b*c^3"),
                MultivariatePolynomial.parse("1-10526*a*b*c^2+11461*a*b^2*c^3+107742*a^3-113101*a^3*c+89716*a^3*b*c^3"),
                MultivariatePolynomial.parse("1-77553*b^3+60044*a^3*b^2*c^2"),
                MultivariatePolynomial.parse("1+102974*b*c-19396*a^2*b^2*c"),
                MultivariatePolynomial.parse("1+81674*b^5*c^6+107381*a*b^6*c^3+5984*a^5*b^5*c^6-96446*a^6*c^6"),
        };

//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 1, 2)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);

        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(5, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization32() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1-9*b*c+435*b^2*c-460*a*b-410*a*b^2+4*a*b^2*c"),
                MultivariatePolynomial.parse("1-157*a^2*c^2+340*a^2*b*c"),
                MultivariatePolynomial.parse("(-457)*c+120*b-141*a^2*b^2*c^2"),
                MultivariatePolynomial.parse("1-458*a^2-118*a^2*b^2*c^3+318*a^2*b^3-256*a^2*b^4*c^2+238*a^2*b^4*c^4+192*a^5*b^3"),
        };
        MultivariatePolynomial<BigInteger> base = multiply(pp);

        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 0; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(4, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization33() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+106564457*a*b^4*c*d^2+325561484*a^3*b*d^3-494009011*a^3*b^4*d^2-340064640*a^4*b*c^5+5486146*a^5*b^3*c^4*d^2"),
                MultivariatePolynomial.parse("1-263084145*a^2*b^5*d^2+435163540*a^2*b^5*c^3*d^5+538303140*a^4*b^5*c^4*d^5"),
                MultivariatePolynomial.parse("303047887*b^5*c^4*d^5+65189925*a*c^5+310563050*a^4*b^4*d^5+197570681*a^5*b^5*c^4*d^2"),
                MultivariatePolynomial.parse("1+332010051*a*c^4*d^2+112027379*a^3*b*c^2*d^4+51064879*a^3*b^4*c*d^3-396600905*a^4-10945831*a^5*b^3*c^3*d"),
                MultivariatePolynomial.parse("1-327267807*b*c^4+106677035*a*b^3*d^2-379344971*a^5*b^5*c^3*d^2"),
        };
//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 0, 3)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 1; i < its(1, 1); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(5, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization34() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1-3778*a^2*b^2*c+13207*a^3*b*c^4"),
                MultivariatePolynomial.parse("1-13694*a^3*b^5*c-16924*a^6*c^5+34818*a^6*b^6*c"),
                MultivariatePolynomial.parse("1+30327*a^3*c^2+24272*a^3*b*c-5401*a^3*b^2-40302*a^3*b^4"),
        };
//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 1, 2)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);
        //System.out.println(base);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 3; i < its(10, 10); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(3, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization35() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("57335917*d-53377618*b*c^2+20557031*a*b^2*c*d+28290557*a*b^2*c^2*d^2+219143656*a^2*b*d+80576655*a^2*b^2*c^2*d"),
                MultivariatePolynomial.parse("77388697*c^2-29509948*a^2*c^3*d-15585632*a^2*b^3*c^2*d+193714032*a^3*b^2*d^3-129223469*a^3*b^3*c^3*d"),
                MultivariatePolynomial.parse("1-207609883*a*b*c*d^2+165695566*a*b^2*c^2*d^2+116765978*a^2*c*d^2+89347609*a^2*c^2*d^2+3700023*a^2*b*c^2*d^2"),
                MultivariatePolynomial.parse("1-20532829*a^4*b^4*c^3*d^3+181984105*a^5*b^3*c^5*d^4+142483093*a^5*b^5*c^2"),
        };
//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 1, 2)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);
        //System.out.println(base);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 7; i < its(10, 10); i++) {
            System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(4, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization36() throws Exception {
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("1+324*c^2+119*a*b^3*c+461*a^2*c^3+597*a^3*b*c^3"),
                MultivariatePolynomial.parse("1+141*b^3*c^5-134*a*b^4*c"),
                MultivariatePolynomial.parse("1-519*a*b^2+98*a^2*b*c+362*a^2*b^2*c"),
        };
//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 1, 2)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);
        //System.out.println(base);
        assert MultivariateSquareFreeFactorization.isSquareFree(base);
        for (int i = 7; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorPrimitiveInZ(base);
            Assert.assertEquals(3, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization37() throws Exception {
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger>[] pp = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("(-1787086)*c^2+135144*b^2*c^2-459684*a*c^2+1379391*a^2*b*c-1357466*a^2*b*c^2", vars),
                MultivariatePolynomial.parse("1295770*a^2*b*c^2-418269*a^4*b*c+414642*a^6*b^6*c^3", vars),
                MultivariatePolynomial.parse("1274683*a+56290*a*c", vars),
                MultivariatePolynomial.parse("(-493283)*a^4*b^2*c^3-1192675*a^4*b^5*c^2", vars),
                MultivariatePolynomial.parse("476108*a*c^2-1504550*a*b*c-1519177*a*b*c^2+1007772*a*b^2*c", vars),
        };
//        for (MultivariatePolynomial<BigInteger> f : pp)
//            System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", AMultivariatePolynomial.swapVariables(f, 1, 2)));

        MultivariatePolynomial<BigInteger> base = multiply(pp);
        //System.out.println(base);
        for (int i = 1; i < its(10, 10); i++) {
            //System.out.println(i);
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.factorToPrimitive(base);
            Assert.assertEquals(3, factors.size());
            FactorDecompositionTest.assertFactorization(base, factors);
            //System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization38RandomZ() throws Exception {
        //DETALIZATION_PERCENT = 5000;
        //PRINT_FACTORS = true;

        FactorizationInput.SampleDecompositionSource<MultivariatePolynomial<BigInteger>> source = new SampleDecompositionSourceZ(new lSampleDecompositionSource(
                3, 5,
                3, 4,
                2, 6,
                1, 5));
        testFactorizationAlgorithm(filterNonPrimitive(filterNonSquareFree(filterMonomialContent(source))), its(25, 75),
                FactorizationInput.FactorizationAlgorithm.named(MultivariateFactorization::factorPrimitiveInZ, "Multivariate factorization over Z"));
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization39() throws Exception {
        String[] vars = {"a", "b", "c"};
        MultivariatePolynomial<BigInteger> pp[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("a^2 + b^2 + c^2", vars),
                MultivariatePolynomial.parse("a*b - b*c - c*a", vars).square(),
                MultivariatePolynomial.parse("3", vars),
                MultivariatePolynomial.parse("a - b - c", vars).square(),
                MultivariatePolynomial.parse("a - b ", vars).square(),
                MultivariatePolynomial.parse("b + c", vars).square(),
                MultivariatePolynomial.parse("a^6 + c^2 + a*c^2 + 1", vars).square()
        }, poly = multiply(pp);

        PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.Factor(poly);
        Assert.assertEquals(6, factors.size());
        Assert.assertEquals(poly, factors.multiply());

        poly = poly.setRing(new IntegersZp(Integer.MAX_VALUE));
        factors = MultivariateFactorization.Factor(poly);
        Assert.assertEquals(6, factors.size());
        Assert.assertEquals(poly, factors.multiply());
    }


    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization40() throws Exception {
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomial<BigInteger> pp[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("31694*c^2*d^2-9343*b^2-12027*b^2*d^2+40289*a^2*b*c*d", vars),
                MultivariatePolynomial.parse("(-114834)*b^2*d+119706*a*c^2-105913*a*b^2+129857*a*b^2*d", vars),
                MultivariatePolynomial.parse("(-27997)*b^4*d^5-120509*b^5*c^3*d^5-103731*a*c^3*d^2-71946*a*b^4*c^5*d^5+113296*a^3*b^4*c", vars),
                MultivariatePolynomial.parse("1-16532*b*c*d^2-101779*a*b^2*c^4*d^4-74636*a^2*b^2*c^4-110612*a^3*c^4*d^4-97709*a^3*b^2*c+120010*a^5*b^3*c*d^2", vars),
        }, poly = multiply(pp);

        for (int i = 21; i < 22; i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factors = MultivariateFactorization.Factor(poly);
            Assert.assertEquals(4, factors.size());
            Assert.assertEquals(poly, factors.multiply());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testMultivariateFactorization41() throws Exception {
        String[] vars = {"a", "b", "c", "d"};
        MultivariatePolynomial<Rational<BigInteger>> pp[] = new MultivariatePolynomial[]{
                MultivariatePolynomial.parse("(31694/3)*c^2*d^2-(9343/1111)*b^2-12027*b^2*d^2+40289*a^2*b*c*d", Rings.Q, vars),
                MultivariatePolynomial.parse("(-114834)*b^2*d+119706*a*c^2-105913*a*b^2+129857*a*b^2*d", Rings.Q, vars),
                MultivariatePolynomial.parse("(-27997)*b^4*d^5-(120509/123)*b^5*c^3*d^5-103731*a*c^3*d^2-71946*a*b^4*c^5*d^5+113296*a^3*b^4*c", Rings.Q, vars),
                MultivariatePolynomial.parse("1-16532*b*c*d^2-101779*a*b^2*c^4*d^4-(74636/1234342)*a^2*b^2*c^4-110612*a^3*c^4*d^4-97709*a^3*b^2*c+120010*a^5*b^3*c*d^2", Rings.Q, vars),
        }, poly = multiply(pp);

        for (int i = 21; i < 22; i++) {
            PrivateRandom.getRandom().setSeed(i);
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<Rational<BigInteger>>> factors = MultivariateFactorization.Factor(poly);
            Assert.assertEquals(4, factors.size());
            Assert.assertEquals(poly, factors.multiply());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testMultivariateFactorization42() throws Exception {
        Ring<BigInteger> ring = Rings.Zp(524287);
        MultivariatePolynomial<BigInteger>
                p1 = polyPow(MultivariatePolynomial.parse("1 + 3*a*b + 5*b*c + 7*c*d + 9*d*e + 11*e*f + 13*f*g + 15*g*a", ring), 3),
                p2 = polyPow(MultivariatePolynomial.parse("1 + 3*a*c + 5*b*d + 7*c*e + 9*f*e + 11*g*f + 13*f*a + 15*g*b", ring), 3),
                p3 = polyPow(MultivariatePolynomial.parse("1 + 3*a*d + 5*b*e + 7*c*f + 9*f*g + 11*g*a + 13*f*b + 15*g*c", ring), 3),
                poly = p1.multiply(p2, p3);
        poly.decrement();
        Assert.assertEquals(3, factorPrimitiveInGF(poly).size());
    }

    @Test(timeout = 600_000)
    public void testMultivariateFactorization43() throws Exception {
        Ring<BigInteger> ring = Rings.Z;
        MultivariatePolynomial<BigInteger>
                p1 = polyPow(MultivariatePolynomial.parse("1 + 3*a*b + 5*b*c + 7*c*d + 9*d*e + 11*e*f + 13*f*g + 15*g*a", ring), 3),
                p2 = polyPow(MultivariatePolynomial.parse("1 + 3*a*c + 5*b*d + 7*c*e + 9*f*e + 11*g*f + 13*f*a + 15*g*b", ring), 3),
                p3 = polyPow(MultivariatePolynomial.parse("1 + 3*a*d + 5*b*e + 7*c*f + 9*f*g + 11*g*a + 13*f*b + 15*g*c", ring), 3),
                poly = p1.multiply(p2, p3);
        poly.decrement();
        List<MultivariatePolynomial<BigInteger>> l = new ArrayList<>(Arrays.asList(poly.derivative()));
        l.add(poly);
        Assert.assertEquals(2, Factor(poly).size());
    }

    @Test(timeout = 10000)
    public void testMultivariateFactorization44() throws Exception {
        for (int i = 0; i < 50; i++) {
            PrivateRandom.getRandom().setSeed(i);
            MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.parse("25104542697922606374180382796611903776*x^5*y^7*z^4+20534636095763889368932021297717635120*x^6*y^7*z^5+76297700033840376318392368872790943968*x^7*y^4*z^7-26574785740999558717642330744030288047*x^7*y^6*z^4+17324112943714387939253008546475979840*x^8*y^3*z^3+62408844645806356433165500595194302160*x^8*y^4*z^8+14170517235134724116322359611873960800*x^9*y^3*z^4-80766061159845213267443892965900078721*x^9*y^3*z^7-18338696512889824461965469610690129230*x^10*y^2*z^3");
            Assert.assertEquals(5, Factor(poly).size());
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testNestedRing1() throws Exception {
        for (Ring<?> inner1 : Arrays.asList(
                Rings.GF(7, 3),
                Rings.UnivariateRing(Rings.Z),
                Rings.UnivariateRingZp64(17))) {
            MultivariateRing<?> inner2 = MultivariateRing(2, inner1);
            UnivariateRing<?> inner3 = Rings.UnivariateRing(inner2);
            MultivariateRing<?> ring = MultivariateRing(2, inner3);

            for (int i = 0; i < its(4, 8); ++i) {
                AMultivariatePolynomial
                        a = ring.randomElement(3, 5),
                        b = ring.randomElement(3, 5);

                AMultivariatePolynomial poly = multiply(a, b);
                PolynomialFactorDecomposition<?> factors = Factor(poly);
                Assert.assertTrue(factors.sumExponents() >= 2);
                Assert.assertEquals(poly, factors.multiply());
            }
        }
    }

    @Test
    @SuppressWarnings("unchecked")
    public void testNestedRing2() throws Exception {

        for (UnivariateRing inner2 : Arrays.asList(
                Rings.UnivariateRingZp64(17),
                Rings.UnivariateRing(Rings.Z)
        )) {

            for (int i = 0; i < 8; i++) {
                System.out.println(i);
                Rationals<IUnivariatePolynomial> inner3 = Frac(inner2);
                MultivariateRing<MultivariatePolynomial<Rational>> inner4 = (MultivariateRing) MultivariateRing(2, inner3);
                UnivariateRing<UnivariatePolynomial<MultivariatePolynomial<Rational>>> ring = UnivariateRing(inner4);

                RandomGenerator rnd = PrivateRandom.getRandom();
                UnivariatePolynomial<MultivariatePolynomial<Rational<UnivariatePolynomialZp64>>> poly =
                        IntStream.range(0, 4)
                                .mapToObj(__ -> {
                                    IUnivariatePolynomial num = inner2.randomElement(3, rnd);
                                    IUnivariatePolynomial den = inner2.randomElement(3, rnd);
                                    if (num instanceof UnivariatePolynomial) {
                                        num = ((UnivariatePolynomial<BigInteger>) num).setRing(Zp(10)).setRing(Z);
                                        den = ((UnivariatePolynomial<BigInteger>) den).setRing(Zp(10)).setRing(Z);
                                    }
                                    Rational rat = new Rational(inner2, num, den);
                                    UnivariatePolynomial r = ring.randomElement(3, rnd).add(inner4.randomElement(3, 3).add(rat));
                                    return r;
                                }).reduce(ring.getOne(), (a, b) -> ring.multiply((UnivariatePolynomial) a, (UnivariatePolynomial) b));

                PolynomialFactorDecomposition<?> fs = Factor(poly);

                Assert.assertTrue(4 <= fs.sumExponents());
                Assert.assertEquals(poly, fs.multiply());

            }
        }
    }

    @Test
    public void testMultivariateFactorization45() throws Exception {
        MultivariatePolynomialZp64 poly = MultivariatePolynomialZp64.parse(
                "7+4*x4^2+8*x4^4+16*x4^6+6*x3+9*x3*x4+x3*x4^2+12*x3*x4^3+2*x3*x4^4+8*x3*x4^5+4*x3*x4^6+8*x3*x4^7+16*x3*x4^9+2*x3^2+15*x3^2*x4+13*x3^2*x4^2+3*x3^2*x4^3+7*x3^2*x4^4+2*x3^2*x4^5+5*x3^2*x4^6+2*x3^2*x4^7+4*x3^2*x4^8+4*x3^2*x4^9+9*x3^3+11*x3^3*x4+16*x3^3*x4^2+2*x3^3*x4^3+6*x3^3*x4^4+10*x3^3*x4^5+14*x3^3*x4^6+15*x3^3*x4^7+x3^3*x4^8+15*x3^4+7*x3^4*x4+14*x3^4*x4^2+9*x3^4*x4^3+11*x3^4*x4^4+11*x3^4*x4^5+8*x3^4*x4^7+14*x3^4*x4^8+8*x3^5+4*x3^5*x4+12*x3^5*x4^2+9*x3^5*x4^3+7*x3^5*x4^4+14*x3^5*x4^5+12*x3^5*x4^8+x3^6*x4+15*x3^6*x4^3+4*x3^6*x4^4+12*x3^6*x4^5+11*x3^6*x4^6+8*x3^7*x4+x3^7*x4^4+7*x3^7*x4^6+2*x3^8*x4+8*x3^8*x4^4+2*x3^9*x4^4+11*x2+15*x2*x4^2+11*x2*x4^4+x2*x4^6+15*x2*x3+4*x2*x3*x4+15*x2*x3*x4^2+x2*x3*x4^3+4*x2*x3*x4^4+5*x2*x3*x4^5+7*x2*x3*x4^6+16*x2*x3*x4^7+11*x2*x3*x4^9+3*x2*x3^2+4*x2*x3^2*x4+9*x2*x3^2*x4^2+9*x2*x3^2*x4^3+x2*x3^2*x4^4+4*x2*x3^2*x4^5+10*x2*x3^2*x4^6+4*x2*x3^2*x4^7+7*x2*x3^2*x4^8+7*x2*x3^2*x4^9+16*x2*x3^3*x4+8*x2*x3^3*x4^2+15*x2*x3^3*x4^3+x2*x3^3*x4^4+x2*x3^3*x4^5+11*x2*x3^3*x4^6+5*x2*x3^3*x4^7+6*x2*x3^3*x4^8+10*x2*x3^4+2*x2*x3^4*x4+3*x2*x3^4*x4^2+8*x2*x3^4*x4^3+8*x2*x3^4*x4^4+12*x2*x3^4*x4^5+15*x2*x3^4*x4^6+14*x2*x3^4*x4^7+16*x2*x3^4*x4^8+16*x2*x3^5+12*x2*x3^5*x4+4*x2*x3^5*x4^2+2*x2*x3^5*x4^3+2*x2*x3^5*x4^4+16*x2*x3^5*x4^5+8*x2*x3^5*x4^6+4*x2*x3^5*x4^8+10*x2*x3^6*x4+7*x2*x3^6*x4^3+6*x2*x3^6*x4^4+4*x2*x3^6*x4^5+15*x2*x3^6*x4^6+9*x2*x3^7*x4+10*x2*x3^7*x4^4+8*x2*x3^7*x4^6+12*x2*x3^8*x4+14*x2*x3^8*x4^4+12*x2*x3^9*x4^4+9*x2^2+8*x2^2*x4^2+9*x2^2*x4^4+9*x2^2*x4^6+8*x2^2*x3+14*x2^2*x3*x4+9*x2^2*x3*x4^2+9*x2^2*x3*x4^3+12*x2^2*x3*x4^6+14*x2^2*x3*x4^9+6*x2^2*x3^2+15*x2^2*x3^2*x4+11*x2^2*x3^2*x4^2+12*x2^2*x3^2*x4^3+10*x2^2*x3^2*x4^4+12*x2^2*x3^2*x4^5+16*x2^2*x3^2*x4^6+12*x2^2*x3^2*x4^8+12*x2^2*x3^2*x4^9+12*x2^2*x3^3+14*x2^2*x3^3*x4+8*x2^2*x3^3*x4^2+6*x2^2*x3^3*x4^3+5*x2^2*x3^3*x4^4+8*x2^2*x3^3*x4^5+4*x2^2*x3^3*x4^6+11*x2^2*x3^3*x4^7+3*x2^2*x3^3*x4^8+14*x2^2*x3^4+3*x2^2*x3^4*x4+10*x2^2*x3^4*x4^2+15*x2^2*x3^4*x4^3+12*x2^2*x3^4*x4^4+10*x2^2*x3^4*x4^5+12*x2^2*x3^4*x4^6+7*x2^2*x3^4*x4^7+8*x2^2*x3^4*x4^8+3*x2^2*x3^5+16*x2^2*x3^5*x4+2*x2^2*x3^5*x4^2+15*x2^2*x3^5*x4^3+3*x2^2*x3^5*x4^4+8*x2^2*x3^5*x4^5+3*x2^2*x3^5*x4^6+2*x2^2*x3^5*x4^8+10*x2^2*x3^6*x4+7*x2^2*x3^6*x4^3+6*x2^2*x3^6*x4^4+2*x2^2*x3^6*x4^5+16*x2^2*x3^6*x4^6+13*x2^2*x3^7*x4+10*x2^2*x3^7*x4^4+4*x2^2*x3^7*x4^6+6*x2^2*x3^8*x4+7*x2^2*x3^8*x4^4+6*x2^2*x3^9*x4^4+2*x2^3+15*x2^3*x4^2+x2^3*x4^4+8*x2^3*x4^6+9*x2^3*x3+16*x2^3*x3*x4+4*x2^3*x3*x4^2+15*x2^3*x3*x4^3+13*x2^3*x3*x4^4+x2^3*x3*x4^5+14*x2^3*x3*x4^6+x2^3*x3*x4^7+5*x2^3*x3*x4^9+3*x2^3*x3^2+5*x2^3*x3^2*x4+x2^3*x3^2*x4^2+2*x2^3*x3^2*x4^3+2*x2^3*x3^2*x4^4+16*x2^3*x3^2*x4^5+6*x2^3*x3^2*x4^6+13*x2^3*x3^2*x4^7+14*x2^3*x3^2*x4^8+14*x2^3*x3^2*x4^9+x2^3*x3^3*x4+7*x2^3*x3^3*x4^2+7*x2^3*x3^3*x4^3+16*x2^3*x3^3*x4^4+x2^3*x3^3*x4^5+10*x2^3*x3^3*x4^6+10*x2^3*x3^3*x4^7+12*x2^3*x3^3*x4^8+3*x2^3*x3^4*x4+7*x2^3*x3^4*x4^2+8*x2^3*x3^4*x4^3+12*x2^3*x3^4*x4^4+15*x2^3*x3^4*x4^5+15*x2^3*x3^4*x4^6+11*x2^3*x3^4*x4^7+15*x2^3*x3^4*x4^8+4*x2^3*x3^5+9*x2^3*x3^5*x4+8*x2^3*x3^5*x4^2+13*x2^3*x3^5*x4^3+3*x2^3*x3^5*x4^4+15*x2^3*x3^5*x4^5+8*x2^3*x3^5*x4^6+8*x2^3*x3^5*x4^8+15*x2^3*x3^6*x4+3*x2^3*x3^6*x4^3+9*x2^3*x3^6*x4^4+8*x2^3*x3^6*x4^5+13*x2^3*x3^6*x4^6+4*x2^3*x3^7*x4+15*x2^3*x3^7*x4^4+16*x2^3*x3^7*x4^6+7*x2^3*x3^8*x4+11*x2^3*x3^8*x4^4+7*x2^3*x3^9*x4^4+12*x2^4+14*x2^4*x4^2+5*x2^4*x4^4+9*x2^4*x4^6+15*x2^4*x3+8*x2^4*x3*x4+16*x2^4*x3*x4^2+8*x2^4*x3*x4^3+6*x2^4*x3*x4^4+11*x2^4*x3*x4^5+5*x2^4*x3*x4^6+7*x2^4*x3*x4^7+3*x2^4*x3*x4^9+4*x2^4*x3^2+13*x2^4*x3^2*x4+12*x2^4*x3^2*x4^2+7*x2^4*x3^2*x4^3+2*x2^4*x3^2*x4^4+13*x2^4*x3^2*x4^5+2*x2^4*x3^2*x4^6+6*x2^4*x3^2*x4^7+5*x2^4*x3^2*x4^8+5*x2^4*x3^2*x4^9+11*x2^4*x3^3+x2^4*x3^3*x4+11*x2^4*x3^3*x4^2+11*x2^4*x3^3*x4^3+6*x2^4*x3^3*x4^4+3*x2^4*x3^3*x4^5+9*x2^4*x3^3*x4^6+6*x2^4*x3^3*x4^7+14*x2^4*x3^3*x4^8+10*x2^4*x3^4+8*x2^4*x3^4*x4+10*x2^4*x3^4*x4^2+11*x2^4*x3^4*x4^3+15*x2^4*x3^4*x4^4+9*x2^4*x3^4*x4^5+10*x2^4*x3^4*x4^7+9*x2^4*x3^4*x4^8+12*x2^4*x3^5+10*x2^4*x3^5*x4+15*x2^4*x3^5*x4^2+8*x2^4*x3^5*x4^4+9*x2^4*x3^5*x4^5+15*x2^4*x3^5*x4^8+8*x2^4*x3^6*x4^3+15*x2^4*x3^6*x4^5+x2^4*x3^6*x4^6+13*x2^4*x3^7*x4+13*x2^4*x3^7*x4^6+11*x2^4*x3^8*x4+10*x2^4*x3^8*x4^4+11*x2^4*x3^9*x4^4+2*x2^5+16*x2^5*x4^2+8*x2^5*x4^4+5*x2^5*x4^6+x2^5*x3+10*x2^5*x3*x4+11*x2^5*x3*x4^2+x2^5*x3*x4^3+2*x2^5*x3*x4^4+7*x2^5*x3*x4^5+16*x2^5*x3*x4^6+8*x2^5*x3*x4^7+13*x2^5*x3*x4^9+4*x2^5*x3^2+15*x2^5*x3^2*x4+13*x2^5*x3^2*x4^2+4*x2^5*x3^2*x4^3+2*x2^5*x3^2*x4^4+15*x2^5*x3^2*x4^5+12*x2^5*x3^2*x4^6+2*x2^5*x3^2*x4^7+16*x2^5*x3^2*x4^8+16*x2^5*x3^2*x4^9+9*x2^5*x3^3+7*x2^5*x3^3*x4+7*x2^5*x3^3*x4^2+11*x2^5*x3^3*x4^3+13*x2^5*x3^3*x4^4+9*x2^5*x3^3*x4^5+3*x2^5*x3^3*x4^6+9*x2^5*x3^3*x4^7+4*x2^5*x3^3*x4^8+9*x2^5*x3^4+4*x2^5*x3^4*x4+15*x2^5*x3^4*x4^2+2*x2^5*x3^4*x4^3+16*x2^5*x3^4*x4^4+4*x2^5*x3^4*x4^5+6*x2^5*x3^4*x4^6+15*x2^5*x3^4*x4^7+5*x2^5*x3^4*x4^8+4*x2^5*x3^5+12*x2^5*x3^5*x4+14*x2^5*x3^5*x4^2+12*x2^5*x3^5*x4^3+4*x2^5*x3^5*x4^4+5*x2^5*x3^5*x4^5+10*x2^5*x3^5*x4^6+14*x2^5*x3^5*x4^8+8*x2^5*x3^6*x4+15*x2^5*x3^6*x4^3+15*x2^5*x3^6*x4^4+14*x2^5*x3^6*x4^5+10*x2^5*x3^6*x4^6+11*x2^5*x3^7*x4+8*x2^5*x3^7*x4^4+11*x2^5*x3^7*x4^6+8*x2^5*x3^8*x4+15*x2^5*x3^8*x4^4+8*x2^5*x3^9*x4^4+4*x2^6+x2^6*x4^2+6*x2^6*x4^4+2*x2^6*x4^6+4*x2^6*x3+10*x2^6*x3*x4+6*x2^6*x3*x4^2+8*x2^6*x3*x4^3+16*x2^6*x3*x4^4+16*x2^6*x3*x4^5+6*x2^6*x3*x4^6+13*x2^6*x3*x4^7+7*x2^6*x3*x4^9+13*x2^6*x3^2+11*x2^6*x3^2*x4+2*x2^6*x3^2*x4^2+9*x2^6*x3^2*x4^3+14*x2^6*x3^2*x4^4+16*x2^6*x3^2*x4^5+3*x2^6*x3^2*x4^6+16*x2^6*x3^2*x4^7+6*x2^6*x3^2*x4^8+6*x2^6*x3^2*x4^9+8*x2^6*x3^3+8*x2^6*x3^3*x4+10*x2^6*x3^3*x4^2+5*x2^6*x3^3*x4^3+6*x2^6*x3^3*x4^4+12*x2^6*x3^3*x4^5+5*x2^6*x3^3*x4^6+14*x2^6*x3^3*x4^7+10*x2^6*x3^3*x4^8+12*x2^6*x3^4+9*x2^6*x3^4*x4+6*x2^6*x3^4*x4^2+3*x2^6*x3^4*x4^3+2*x2^6*x3^4*x4^4+11*x2^6*x3^4*x4^5+14*x2^6*x3^4*x4^6+12*x2^6*x3^4*x4^7+4*x2^6*x3^4*x4^8+5*x2^6*x3^5+9*x2^6*x3^5*x4+x2^6*x3^5*x4^2+15*x2^6*x3^5*x4^3+9*x2^6*x3^5*x4^4+4*x2^6*x3^5*x4^5+12*x2^6*x3^5*x4^6+x2^6*x3^5*x4^8+16*x2^6*x3^6*x4+7*x2^6*x3^6*x4^3+13*x2^6*x3^6*x4^4+x2^6*x3^6*x4^5+8*x2^6*x3^6*x4^6+x2^6*x3^7*x4+16*x2^6*x3^7*x4^4+2*x2^6*x3^7*x4^6+3*x2^6*x3^8*x4+12*x2^6*x3^8*x4^4+3*x2^6*x3^9*x4^4+10*x2^7+4*x2^7*x4^2+10*x2^7*x4^4+9*x2^7*x4^6+4*x2^7*x3+10*x2^7*x3*x4+2*x2^7*x3*x4^2+8*x2^7*x3*x4^3+9*x2^7*x3*x4^4+6*x2^7*x3*x4^5+x2^7*x3*x4^6+2*x2^7*x3*x4^7+4*x2^7*x3*x4^9+12*x2^7*x3^2+4*x2^7*x3^2*x4+15*x2^7*x3^2*x4^5+x2^7*x3^2*x4^6+9*x2^7*x3^2*x4^7+x2^7*x3^2*x4^8+x2^7*x3^2*x4^9+12*x2^7*x3^3+14*x2^7*x3^3*x4+8*x2^7*x3^3*x4^2+16*x2^7*x3^3*x4^3+6*x2^7*x3^3*x4^4+8*x2^7*x3^3*x4^5+13*x2^7*x3^3*x4^6+8*x2^7*x3^3*x4^7+13*x2^7*x3^3*x4^8+7*x2^7*x3^4+12*x2^7*x3^4*x4+10*x2^7*x3^4*x4^2+7*x2^7*x3^4*x4^3+14*x2^7*x3^4*x4^4+11*x2^7*x3^4*x4^5+14*x2^7*x3^4*x4^6+2*x2^7*x3^4*x4^7+12*x2^7*x3^4*x4^8+7*x2^7*x3^5+14*x2^7*x3^5*x4+3*x2^7*x3^5*x4^2+14*x2^7*x3^5*x4^3+12*x2^7*x3^5*x4^4+12*x2^7*x3^5*x4^5+12*x2^7*x3^5*x4^6+3*x2^7*x3^5*x4^8+15*x2^7*x3^6*x4+13*x2^7*x3^6*x4^3+9*x2^7*x3^6*x4^4+3*x2^7*x3^6*x4^5+7*x2^7*x3^6*x4^6+13*x2^7*x3^7*x4+15*x2^7*x3^7*x4^4+6*x2^7*x3^7*x4^6+9*x2^7*x3^8*x4+2*x2^7*x3^8*x4^4+9*x2^7*x3^9*x4^4+13*x2^8+x2^8*x4^2+2*x2^8*x4^4+12*x2^8*x4^6+16*x2^8*x3+14*x2^8*x3*x4+2*x2^8*x3*x4^3+10*x2^8*x3*x4^4+3*x2^8*x3*x4^5+9*x2^8*x3*x4^6+6*x2^8*x3*x4^7+2*x2^8*x3*x4^9+15*x2^8*x3^2*x4+15*x2^8*x3^2*x4^2+11*x2^8*x3^2*x4^3+16*x2^8*x3^2*x4^4+15*x2^8*x3^2*x4^5+14*x2^8*x3^2*x4^6+10*x2^8*x3^2*x4^7+9*x2^8*x3^2*x4^8+9*x2^8*x3^2*x4^9+10*x2^8*x3^3*x4^3+16*x2^8*x3^3*x4^4+2*x2^8*x3^3*x4^5+12*x2^8*x3^3*x4^6+4*x2^8*x3^3*x4^7+15*x2^8*x3^3*x4^8+2*x2^8*x3^4+12*x2^8*x3^4*x4+2*x2^8*x3^4*x4^2+4*x2^8*x3^4*x4^3+14*x2^8*x3^4*x4^4+10*x2^8*x3^4*x4^5+16*x2^8*x3^4*x4^6+x2^8*x3^4*x4^7+6*x2^8*x3^4*x4^8+5*x2^8*x3^5+16*x2^8*x3^5*x4+10*x2^8*x3^5*x4^2+7*x2^8*x3^5*x4^3+12*x2^8*x3^5*x4^4+6*x2^8*x3^5*x4^5+4*x2^8*x3^5*x4^6+10*x2^8*x3^5*x4^8+x2^8*x3^6*x4+8*x2^8*x3^6*x4^3+4*x2^8*x3^6*x4^4+10*x2^8*x3^6*x4^5+12*x2^8*x3^6*x4^6+6*x2^8*x3^7*x4+x2^8*x3^7*x4^4+3*x2^8*x3^7*x4^6+13*x2^8*x3^8*x4+x2^8*x3^8*x4^4+13*x2^8*x3^9*x4^4+3*x2^9+6*x2^9*x4^2+2*x2^9*x4^4+5*x2^9*x4^6+3*x2^9*x3+7*x2^9*x3*x4+8*x2^9*x3*x4^2+8*x2^9*x3*x4^3+10*x2^9*x3*x4^4+12*x2^9*x3*x4^5+6*x2^9*x3*x4^7+16*x2^9*x3^2*x4+2*x2^9*x3^2*x4^2+x2^9*x3^2*x4^3+6*x2^9*x3^2*x4^4+8*x2^9*x3^2*x4^5+9*x2^9*x3^2*x4^6+10*x2^9*x3^2*x4^7+16*x2^9*x3^3*x4+2*x2^9*x3^3*x4^2+14*x2^9*x3^3*x4^3+16*x2^9*x3^3*x4^4+6*x2^9*x3^3*x4^5+15*x2^9*x3^3*x4^6+3*x2^9*x3^4+5*x2^9*x3^4*x4+15*x2^9*x3^4*x4^2+10*x2^9*x3^4*x4^3+3*x2^9*x3^4*x4^4+2*x2^9*x3^4*x4^5+6*x2^9*x3^4*x4^6+3*x2^9*x3^5+8*x2^9*x3^5*x4+8*x2^9*x3^5*x4^3+5*x2^9*x3^5*x4^4+10*x2^9*x3^5*x4^6+3*x2^9*x3^6*x4^3+11*x2^9*x3^7*x4+2*x2^10*x4^2+8*x2^10*x4^4+12*x2^10*x4^6+9*x2^10*x3*x4+2*x2^10*x3*x4^2+2*x2^10*x3*x4^3+8*x2^10*x3*x4^4+11*x2^10*x3*x4^5+12*x2^10*x3*x4^6+15*x2^10*x3*x4^7+14*x2^10*x3*x4^9+9*x2^10*x3^2*x4+2*x2^10*x3^2*x4^3+9*x2^10*x3^2*x4^4+5*x2^10*x3^2*x4^5+8*x2^10*x3^2*x4^6+8*x2^10*x3^2*x4^7+12*x2^10*x3^2*x4^8+12*x2^10*x3^2*x4^9+6*x2^10*x3^3*x4+13*x2^10*x3^3*x4^3+16*x2^10*x3^3*x4^4+2*x2^10*x3^3*x4^5+2*x2^10*x3^3*x4^6+11*x2^10*x3^3*x4^7+3*x2^10*x3^3*x4^8+6*x2^10*x3^4*x4+2*x2^10*x3^4*x4^2+13*x2^10*x3^4*x4^3+7*x2^10*x3^4*x4^4+2*x2^10*x3^4*x4^5+x2^10*x3^4*x4^6+7*x2^10*x3^4*x4^7+8*x2^10*x3^4*x4^8+9*x2^10*x3^5*x4+2*x2^10*x3^5*x4^2+4*x2^10*x3^5*x4^3+6*x2^10*x3^5*x4^4+8*x2^10*x3^5*x4^5+13*x2^10*x3^5*x4^6+2*x2^10*x3^5*x4^8+9*x2^10*x3^6*x4+4*x2^10*x3^6*x4^3+2*x2^10*x3^6*x4^4+2*x2^10*x3^6*x4^5+16*x2^10*x3^6*x4^6+6*x2^10*x3^7*x4+9*x2^10*x3^7*x4^4+4*x2^10*x3^7*x4^6+6*x2^10*x3^8*x4+7*x2^10*x3^8*x4^4+6*x2^10*x3^9*x4^4+7*x1+2*x1*x4^2+16*x1*x4^4+4*x1*x4^6+13*x1*x3+2*x1*x3*x4+13*x1*x3*x4^2+x1*x3*x4^3+12*x1*x3*x4^4+2*x1*x3*x4^5+14*x1*x3*x4^7+7*x1*x3^2+x1*x3^2*x4+9*x1*x3^2*x4^2+x1*x3^2*x4^3+12*x1*x3^2*x4^4+13*x1*x3^2*x4^5+3*x1*x3^2*x4^6+12*x1*x3^2*x4^7+8*x1*x3^3+14*x1*x3^3*x4+11*x1*x3^3*x4^2+x1*x3^3*x4^4+5*x1*x3^3*x4^5+5*x1*x3^3*x4^6+14*x1*x3^4+6*x1*x3^4*x4+12*x1*x3^4*x4^2+7*x1*x3^4*x4^4+11*x1*x3^4*x4^5+2*x1*x3^4*x4^6+10*x1*x3^5+2*x1*x3^5*x4+13*x1*x3^5*x4^3+6*x1*x3^5*x4^4+9*x1*x3^5*x4^6+13*x1*x3^6*x4+10*x1*x3^6*x4^3+x1*x3^6*x4^4+2*x1*x3^7*x4+13*x1*x3^7*x4^4+12*x1*x2+2*x1*x2*x4^2+7*x1*x2*x4^4+7*x1*x2*x4^6+x1*x2*x3+16*x1*x2*x3*x4+7*x1*x2*x3*x4^2+12*x1*x2*x3*x4^3+4*x1*x2*x3*x4^4+16*x1*x2*x3*x4^7+15*x1*x2*x3^2+10*x1*x2*x3^2*x4+15*x1*x2*x3^2*x4^2+14*x1*x2*x3^2*x4^3+3*x1*x2*x3^2*x4^4+7*x1*x2*x3^2*x4^5+x1*x2*x3^2*x4^6+4*x1*x2*x3^2*x4^7+13*x1*x2*x3^3+15*x1*x2*x3^3*x4^2+7*x1*x2*x3^3*x4^3+10*x1*x2*x3^3*x4^4+13*x1*x2*x3^3*x4^5+13*x1*x2*x3^3*x4^6+10*x1*x2*x3^4+16*x1*x2*x3^4*x4+4*x1*x2*x3^4*x4^2+16*x1*x2*x3^4*x4^3+13*x1*x2*x3^4*x4^4+15*x1*x2*x3^4*x4^5+12*x1*x2*x3^4*x4^6+9*x1*x2*x3^5+9*x1*x2*x3^5*x4+10*x1*x2*x3^5*x4^3+16*x1*x2*x3^5*x4^4+3*x1*x2*x3^5*x4^6+10*x1*x2*x3^6*x4+9*x1*x2*x3^6*x4^3+6*x1*x2*x3^6*x4^4+12*x1*x2*x3^7*x4+10*x1*x2*x3^7*x4^4+10*x1*x2^2+5*x1*x2^2*x4^2+10*x1*x2^2*x4^4+12*x1*x2^2*x4^6+16*x1*x2^2*x3*x4+15*x1*x2^2*x3*x4^2+11*x1*x2^2*x3*x4^3+2*x1*x2^2*x3*x4^4+12*x1*x2^2*x3*x4^5+8*x1*x2^2*x3*x4^7+9*x1*x2^2*x3^2+2*x1*x2^2*x3^2*x4+3*x1*x2^2*x3^2*x4^2+15*x1*x2^2*x3^2*x4^3+15*x1*x2^2*x3^2*x4^4+15*x1*x2^2*x3^2*x4^5+9*x1*x2^2*x3^2*x4^6+2*x1*x2^2*x3^2*x4^7+16*x1*x2^2*x3^3+10*x1*x2^2*x3^3*x4+16*x1*x2^2*x3^3*x4^2+15*x1*x2^2*x3^3*x4^3+2*x1*x2^2*x3^3*x4^4+15*x1*x2^2*x3^3*x4^5+15*x1*x2^2*x3^3*x4^6+6*x1*x2^2*x3^4*x4+2*x1*x2^2*x3^4*x4^2+9*x1*x2^2*x3^4*x4^3+7*x1*x2^2*x3^4*x4^4+16*x1*x2^2*x3^4*x4^5+6*x1*x2^2*x3^4*x4^6+13*x1*x2^2*x3^5+x1*x2^2*x3^5*x4+5*x1*x2^2*x3^5*x4^3+6*x1*x2^2*x3^5*x4^4+10*x1*x2^2*x3^5*x4^6+5*x1*x2^2*x3^6*x4+13*x1*x2^2*x3^6*x4^3+3*x1*x2^2*x3^6*x4^4+6*x1*x2^2*x3^7*x4+5*x1*x2^2*x3^7*x4^4+15*x1*x2^3+13*x1*x2^3*x4^2+3*x1*x2^3*x4^4+14*x1*x2^3*x4^6+12*x1*x2^3*x3+12*x1*x2^3*x3*x4+10*x1*x2^3*x3*x4^2+8*x1*x2^3*x3*x4^4+x1*x2^3*x3*x4^5+15*x1*x2^3*x3*x4^7+15*x1*x2^3*x3^2+x1*x2^3*x3^2*x4+12*x1*x2^3*x3^2*x4^2+4*x1*x2^3*x3^2*x4^3+15*x1*x2^3*x3^2*x4^4+10*x1*x2^3*x3^2*x4^5+2*x1*x2^3*x3^2*x4^6+8*x1*x2^3*x3^2*x4^7+6*x1*x2^3*x3^3+7*x1*x2^3*x3^3*x4+13*x1*x2^3*x3^3*x4^2+16*x1*x2^3*x3^3*x4^3+x1*x2^3*x3^3*x4^4+9*x1*x2^3*x3^3*x4^5+9*x1*x2^3*x3^3*x4^6+7*x1*x2^3*x3^4+8*x1*x2^3*x3^4*x4+8*x1*x2^3*x3^4*x4^2+12*x1*x2^3*x3^4*x4^3+15*x1*x2^3*x3^4*x4^4+13*x1*x2^3*x3^4*x4^5+7*x1*x2^3*x3^4*x4^6+x1*x2^3*x3^5+7*x1*x2^3*x3^5*x4+3*x1*x2^3*x3^5*x4^3+8*x1*x2^3*x3^5*x4^4+6*x1*x2^3*x3^5*x4^6+3*x1*x2^3*x3^6*x4+x1*x2^3*x3^6*x4^3+12*x1*x2^3*x3^6*x4^4+7*x1*x2^3*x3^7*x4+3*x1*x2^3*x3^7*x4^4+7*x1*x2^4+5*x1*x2^4*x4^2+16*x1*x2^4*x4^4+5*x1*x2^4*x4^6+4*x1*x2^4*x3+14*x1*x2^4*x3*x4+14*x1*x2^4*x3*x4^2+15*x1*x2^4*x3*x4^3+15*x1*x2^4*x3*x4^4+2*x1*x2^4*x3*x4^5+9*x1*x2^4*x3*x4^7+3*x1*x2^4*x3^2+8*x1*x2^4*x3^2*x4+6*x1*x2^4*x3^2*x4^2+6*x1*x2^4*x3^2*x4^3+8*x1*x2^4*x3^2*x4^4+14*x1*x2^4*x3^2*x4^5+8*x1*x2^4*x3^2*x4^6+15*x1*x2^4*x3^2*x4^7+4*x1*x2^4*x3^3+15*x1*x2^4*x3^3*x4+x1*x2^4*x3^3*x4^2+15*x1*x2^4*x3^3*x4^3+8*x1*x2^4*x3^3*x4^4+2*x1*x2^4*x3^3*x4^5+2*x1*x2^4*x3^3*x4^6+9*x1*x2^4*x3^4+12*x1*x2^4*x3^4*x4+15*x1*x2^4*x3^4*x4^2+11*x1*x2^4*x3^4*x4^3+14*x1*x2^4*x3^4*x4^4+x1*x2^4*x3^4*x4^5+11*x1*x2^4*x3^4*x4^6+4*x1*x2^4*x3^5+8*x1*x2^4*x3^5*x4+12*x1*x2^4*x3^5*x4^3+12*x1*x2^4*x3^5*x4^4+7*x1*x2^4*x3^5*x4^6+12*x1*x2^4*x3^6*x4+4*x1*x2^4*x3^6*x4^3+14*x1*x2^4*x3^6*x4^4+11*x1*x2^4*x3^7*x4+12*x1*x2^4*x3^7*x4^4+5*x1*x2^5+6*x1*x2^5*x4^2+16*x1*x2^5*x4^6+6*x1*x2^5*x3+4*x1*x2^5*x3*x4+4*x1*x2^5*x3*x4^2+12*x1*x2^5*x3*x4^3+14*x1*x2^5*x3*x4^4+3*x1*x2^5*x3*x4^5+5*x1*x2^5*x3*x4^7+11*x1*x2^5*x3^2+11*x1*x2^5*x3^2*x4+8*x1*x2^5*x3^2*x4^2+9*x1*x2^5*x3^2*x4^3+8*x1*x2^5*x3^2*x4^4+4*x1*x2^5*x3^2*x4^5+12*x1*x2^5*x3^2*x4^6+14*x1*x2^5*x3^2*x4^7+6*x1*x2^5*x3^3+9*x1*x2^5*x3^3*x4+10*x1*x2^5*x3^3*x4^2+7*x1*x2^5*x3^3*x4^3+11*x1*x2^5*x3^3*x4^4+3*x1*x2^5*x3^3*x4^5+3*x1*x2^5*x3^3*x4^6+8*x1*x2^5*x3^4+6*x1*x2^5*x3^4*x4+14*x1*x2^5*x3^4*x4^2+8*x1*x2^5*x3^4*x4^3+7*x1*x2^5*x3^4*x4^4+10*x1*x2^5*x3^4*x4^5+8*x1*x2^5*x3^4*x4^6+6*x1*x2^5*x3^5+3*x1*x2^5*x3^5*x4+x1*x2^5*x3^5*x4^3+6*x1*x2^5*x3^5*x4^4+2*x1*x2^5*x3^5*x4^6+x1*x2^5*x3^6*x4+6*x1*x2^5*x3^6*x4^3+4*x1*x2^5*x3^6*x4^4+8*x1*x2^5*x3^7*x4+x1*x2^5*x3^7*x4^4+15*x1*x2^6+15*x1*x2^6*x4^2+5*x1*x2^6*x4^4+6*x1*x2^6*x4^6+5*x1*x2^6*x3+5*x1*x2^6*x3*x4+4*x1*x2^6*x3*x4^2+2*x1*x2^6*x3*x4^3+x1*x2^6*x3*x4^4+9*x1*x2^6*x3*x4^5+4*x1*x2^6*x3*x4^7+13*x1*x2^6*x3^2+16*x1*x2^6*x3^2*x4+7*x1*x2^6*x3^2*x4^2+4*x1*x2^6*x3^2*x4^3+8*x1*x2^6*x3^2*x4^4+4*x1*x2^6*x3^2*x4^5+13*x1*x2^6*x3^2*x4^6+x1*x2^6*x3^2*x4^7+5*x1*x2^6*x3^3+2*x1*x2^6*x3^3*x4+8*x1*x2^6*x3^3*x4^2+11*x1*x2^6*x3^3*x4^3+16*x1*x2^6*x3^3*x4^4+16*x1*x2^6*x3^3*x4^5+16*x1*x2^6*x3^3*x4^6+10*x1*x2^6*x3^4+13*x1*x2^6*x3^4*x4+x1*x2^6*x3^4*x4^2+10*x1*x2^6*x3^4*x4^3+x1*x2^6*x3^4*x4^4+8*x1*x2^6*x3^4*x4^5+3*x1*x2^6*x3^4*x4^6+15*x1*x2^6*x3^5+14*x1*x2^6*x3^5*x4+11*x1*x2^6*x3^5*x4^3+13*x1*x2^6*x3^5*x4^4+5*x1*x2^6*x3^5*x4^6+11*x1*x2^6*x3^6*x4+15*x1*x2^6*x3^6*x4^3+10*x1*x2^6*x3^6*x4^4+3*x1*x2^6*x3^7*x4+11*x1*x2^6*x3^7*x4^4+8*x1*x2^7+16*x1*x2^7*x4^2+2*x1*x2^7*x4^4+x1*x2^7*x4^6+6*x1*x2^7*x3+7*x1*x2^7*x3*x4+10*x1*x2^7*x3*x4^2+10*x1*x2^7*x3*x4^3+3*x1*x2^7*x3*x4^4+2*x1*x2^7*x3*x4^5+12*x1*x2^7*x3*x4^7+13*x1*x2^7*x3^2+10*x1*x2^7*x3^2*x4+3*x1*x2^7*x3^2*x4^2+3*x1*x2^7*x3^2*x4^3+8*x1*x2^7*x3^2*x4^4+10*x1*x2^7*x3^2*x4^5+5*x1*x2^7*x3^2*x4^6+3*x1*x2^7*x3^2*x4^7+2*x1*x2^7*x3^3+9*x1*x2^7*x3^3*x4+7*x1*x2^7*x3^3*x4^2+2*x1*x2^7*x3^3*x4^3+10*x1*x2^7*x3^3*x4^4+14*x1*x2^7*x3^3*x4^5+14*x1*x2^7*x3^3*x4^6+4*x1*x2^7*x3^4+8*x1*x2^7*x3^4*x4+3*x1*x2^7*x3^4*x4^2+15*x1*x2^7*x3^4*x4^4+7*x1*x2^7*x3^4*x4^5+9*x1*x2^7*x3^4*x4^6+11*x1*x2^7*x3^5+6*x1*x2^7*x3^5*x4+16*x1*x2^7*x3^5*x4^3+8*x1*x2^7*x3^5*x4^4+15*x1*x2^7*x3^5*x4^6+16*x1*x2^7*x3^6*x4+11*x1*x2^7*x3^6*x4^3+13*x1*x2^7*x3^6*x4^4+9*x1*x2^7*x3^7*x4+16*x1*x2^7*x3^7*x4^4+2*x1*x2^8+14*x1*x2^8*x4^2+12*x1*x2^8*x4^4+9*x1*x2^8*x4^6+7*x1*x2^8*x3+7*x1*x2^8*x3*x4+12*x1*x2^8*x3*x4^2+4*x1*x2^8*x3*x4^3+10*x1*x2^8*x3*x4^4+12*x1*x2^8*x3*x4^5+6*x1*x2^8*x3*x4^7+14*x1*x2^8*x3^2+15*x1*x2^8*x3^2*x4+16*x1*x2^8*x3^2*x4^2+14*x1*x2^8*x3^2*x4^3+10*x1*x2^8*x3^2*x4^4+12*x1*x2^8*x3^2*x4^5+11*x1*x2^8*x3^2*x4^6+10*x1*x2^8*x3^2*x4^7+8*x1*x2^8*x3^3+12*x1*x2^8*x3^3*x4^2+12*x1*x2^8*x3^3*x4^3+15*x1*x2^8*x3^3*x4^4+7*x1*x2^8*x3^3*x4^5+7*x1*x2^8*x3^3*x4^6+x1*x2^8*x3^4+5*x1*x2^8*x3^4*x4+10*x1*x2^8*x3^4*x4^2+7*x1*x2^8*x3^4*x4^3+3*x1*x2^8*x3^4*x4^4+12*x1*x2^8*x3^4*x4^5+13*x1*x2^8*x3^4*x4^6+14*x1*x2^8*x3^5+6*x1*x2^8*x3^5*x4+8*x1*x2^8*x3^5*x4^3+5*x1*x2^8*x3^5*x4^4+16*x1*x2^8*x3^5*x4^6+8*x1*x2^8*x3^6*x4+14*x1*x2^8*x3^6*x4^3+15*x1*x2^8*x3^6*x4^4+13*x1*x2^8*x3^7*x4+8*x1*x2^8*x3^7*x4^4+10*x1*x2^9+12*x1*x2^9*x4^2+8*x1*x2^9*x4^4+3*x1*x2^9*x3+12*x1*x2^9*x3*x4+13*x1*x2^9*x3*x4^2+12*x1*x2^9*x3*x4^3+x1*x2^9*x3*x4^5+11*x1*x2^9*x3^2+15*x1*x2^9*x3^2*x4+3*x1*x2^9*x3^2*x4^2+3*x1*x2^9*x3^2*x4^3+9*x1*x2^9*x3^2*x4^4+13*x1*x2^9*x3^2*x4^5+11*x1*x2^9*x3^3+9*x1*x2^9*x3^3*x4+10*x1*x2^9*x3^3*x4^3+15*x1*x2^9*x3^3*x4^4+7*x1*x2^9*x3^4+10*x1*x2^9*x3^4*x4+11*x1*x2^9*x3^4*x4^3+6*x1*x2^9*x3^4*x4^4+12*x1*x2^9*x3^5*x4+10*x1*x2^9*x3^5*x4^4+13*x1*x2^10+11*x1*x2^10*x4^2+10*x1*x2^10*x4^4+12*x1*x2^10*x4^6+13*x1*x2^10*x3+11*x1*x2^10*x3*x4+9*x1*x2^10*x3*x4^2+x1*x2^10*x3*x4^3+2*x1*x2^10*x3*x4^4+5*x1*x2^10*x3*x4^5+8*x1*x2^10*x3*x4^7+2*x1*x2^10*x3^2*x4+16*x1*x2^10*x3^2*x4^2+11*x1*x2^10*x3^2*x4^3+15*x1*x2^10*x3^2*x4^4+9*x1*x2^10*x3^2*x4^5+9*x1*x2^10*x3^2*x4^6+2*x1*x2^10*x3^2*x4^7+12*x1*x2^10*x3^3*x4+16*x1*x2^10*x3^3*x4^2+6*x1*x2^10*x3^3*x4^3+2*x1*x2^10*x3^3*x4^4+15*x1*x2^10*x3^3*x4^5+15*x1*x2^10*x3^3*x4^6+13*x1*x2^10*x3^4+6*x1*x2^10*x3^4*x4+2*x1*x2^10*x3^4*x4^2+10*x1*x2^10*x3^4*x4^3+7*x1*x2^10*x3^4*x4^4+16*x1*x2^10*x3^4*x4^5+6*x1*x2^10*x3^4*x4^6+13*x1*x2^10*x3^5+14*x1*x2^10*x3^5*x4+5*x1*x2^10*x3^5*x4^3+6*x1*x2^10*x3^5*x4^4+10*x1*x2^10*x3^5*x4^6+5*x1*x2^10*x3^6*x4+13*x1*x2^10*x3^6*x4^3+3*x1*x2^10*x3^6*x4^4+6*x1*x2^10*x3^7*x4+5*x1*x2^10*x3^7*x4^4+14*x1^2+10*x1^2*x4^2+12*x1^2*x4^6+6*x1^2*x3*x4+7*x1^2*x3*x4^2+6*x1^2*x3*x4^3+4*x1^2*x3*x4^4+14*x1^2*x3*x4^5+16*x1^2*x3*x4^7+13*x1^2*x3^2+7*x1^2*x3^2*x4+7*x1^2*x3^2*x4^2+11*x1^2*x3^2*x4^3+x1^2*x3^2*x4^4+7*x1^2*x3^2*x4^5+10*x1^2*x3^2*x4^6+4*x1^2*x3^2*x4^7+3*x1^2*x3^3+5*x1^2*x3^3*x4+4*x1^2*x3^3*x4^2+5*x1^2*x3^3*x4^3+7*x1^2*x3^3*x4^4+x1^2*x3^3*x4^5+11*x1^2*x3^3*x4^6+13*x1^2*x3^4+3*x1^2*x3^4*x4+2*x1^2*x3^4*x4^2+16*x1^2*x3^4*x4^3+12*x1^2*x3^4*x4^4+4*x1^2*x3^4*x4^5+x1^2*x3^4*x4^6+16*x1^2*x3^5+4*x1^2*x3^5*x4+3*x1^2*x3^5*x4^4+13*x1^2*x3^5*x4^6+14*x1^2*x3^6*x4+16*x1^2*x3^6*x4^3+5*x1^2*x3^6*x4^4+6*x1^2*x3^7*x4+14*x1^2*x3^7*x4^4+6*x1^2*x2+9*x1^2*x2*x4^2+4*x1^2*x2*x4^6+7*x1^2*x2*x3+11*x1^2*x2*x3*x4+13*x1^2*x2*x3*x4^3+7*x1^2*x2*x3*x4^4+x1^2*x2*x3*x4^5+11*x1^2*x2*x3*x4^7+15*x1^2*x2*x3^2+8*x1^2*x2*x3^2*x4+6*x1^2*x2*x3^2*x4^2+5*x1^2*x2*x3^2*x4^3+6*x1^2*x2*x3^2*x4^4+9*x1^2*x2*x3^2*x4^6+7*x1^2*x2*x3^2*x4^7+13*x1^2*x2*x3^3+10*x1^2*x2*x3^3*x4+7*x1^2*x2*x3^3*x4^2+10*x1^2*x2*x3^3*x4^3+8*x1^2*x2*x3^3*x4^4+6*x1^2*x2*x3^3*x4^5+15*x1^2*x2*x3^3*x4^6+10*x1^2*x2*x3^4+x1^2*x2*x3^4*x4+12*x1^2*x2*x3^4*x4^2+6*x1^2*x2*x3^4*x4^3+4*x1^2*x2*x3^4*x4^4+7*x1^2*x2*x3^4*x4^5+6*x1^2*x2*x3^4*x4^6+11*x1^2*x2*x3^5+x1^2*x2*x3^5*x4+x1^2*x2*x3^5*x4^4+10*x1^2*x2*x3^5*x4^6+16*x1^2*x2*x3^6*x4+11*x1^2*x2*x3^6*x4^3+13*x1^2*x2*x3^6*x4^4+2*x1^2*x2*x3^7*x4+16*x1^2*x2*x3^7*x4^4+9*x1^2*x2^2+13*x1^2*x2^2*x4^2+11*x1^2*x2^2*x4^4+2*x1^2*x2^2*x4^6+7*x1^2*x2^2*x3+5*x1^2*x2^2*x3*x4+8*x1^2*x2^2*x3*x4^2+16*x1^2*x2^2*x3*x4^3+12*x1^2*x2^2*x3*x4^4+7*x1^2*x2^2*x3*x4^5+14*x1^2*x2^2*x3*x4^7+9*x1^2*x2^2*x3^2+14*x1^2*x2^2*x3^2*x4+15*x1^2*x2^2*x3^2*x4^2+6*x1^2*x2^2*x3^2*x4^3+9*x1^2*x2^2*x3^2*x4^4+8*x1^2*x2^2*x3^2*x4^5+13*x1^2*x2^2*x3^2*x4^6+12*x1^2*x2^2*x3^2*x4^7+3*x1^2*x2^2*x3^3+13*x1^2*x2^2*x3^3*x4+12*x1^2*x2^2*x3^3*x4^2+5*x1^2*x2^2*x3^3*x4^3+14*x1^2*x2^2*x3^3*x4^4+3*x1^2*x2^2*x3^3*x4^5+16*x1^2*x2^2*x3^3*x4^6+7*x1^2*x2^2*x3^4+10*x1^2*x2^2*x3^4*x4+6*x1^2*x2^2*x3^4*x4^2+8*x1^2*x2^2*x3^4*x4^3+6*x1^2*x2^2*x3^4*x4^4+12*x1^2*x2^2*x3^4*x4^5+3*x1^2*x2^2*x3^4*x4^6+14*x1^2*x2^2*x3^5+7*x1^2*x2^2*x3^5*x4+10*x1^2*x2^2*x3^5*x4^4+5*x1^2*x2^2*x3^5*x4^6+8*x1^2*x2^2*x3^6*x4+14*x1^2*x2^2*x3^6*x4^3+15*x1^2*x2^2*x3^6*x4^4+x1^2*x2^2*x3^7*x4+8*x1^2*x2^2*x3^7*x4^4+11*x1^2*x2^3+15*x1^2*x2^3*x4^2+4*x1^2*x2^3*x4^4+8*x1^2*x2^3*x4^6+2*x1^2*x2^3*x3+7*x1^2*x2^3*x3*x4+11*x1^2*x2^3*x3*x4^2+x1^2*x2^3*x3*x4^3+14*x1^2*x2^3*x3*x4^4+12*x1^2*x2^3*x3*x4^5+5*x1^2*x2^3*x3*x4^7+7*x1^2*x2^3*x3^2+8*x1^2*x2^3*x3^2*x4+8*x1^2*x2^3*x3^2*x4^2+15*x1^2*x2^3*x3^2*x4^3+14*x1^2*x2^3*x3^2*x4^4+11*x1^2*x2^3*x3^2*x4^5+x1^2*x2^3*x3^2*x4^6+14*x1^2*x2^3*x3^2*x4^7+11*x1^2*x2^3*x3^3+2*x1^2*x2^3*x3^3*x4+14*x1^2*x2^3*x3^3*x4^2+15*x1^2*x2^3*x3^3*x4^3+8*x1^2*x2^3*x3^3*x4^4+12*x1^2*x2^3*x3^3*x4^5+13*x1^2*x2^3*x3^3*x4^6+11*x1^2*x2^3*x3^4+8*x1^2*x2^3*x3^4*x4+7*x1^2*x2^3*x3^4*x4^2+14*x1^2*x2^3*x3^4*x4^3+15*x1^2*x2^3*x3^4*x4^4+14*x1^2*x2^3*x3^4*x4^5+12*x1^2*x2^3*x3^4*x4^6+5*x1^2*x2^3*x3^5+8*x1^2*x2^3*x3^5*x4+8*x1^2*x2^3*x3^5*x4^4+3*x1^2*x2^3*x3^5*x4^6+15*x1^2*x2^3*x3^6*x4+5*x1^2*x2^3*x3^6*x4^3+9*x1^2*x2^3*x3^6*x4^4+4*x1^2*x2^3*x3^7*x4+15*x1^2*x2^3*x3^7*x4^4+16*x1^2*x2^4+15*x1^2*x2^4*x4^2+8*x1^2*x2^4*x4^4+15*x1^2*x2^4*x4^6+16*x1^2*x2^4*x3+6*x1^2*x2^4*x3*x4+3*x1^2*x2^4*x3*x4^2+3*x1^2*x2^4*x3*x4^3+5*x1^2*x2^4*x3*x4^4+3*x1^2*x2^4*x3*x4^5+3*x1^2*x2^4*x3*x4^7+3*x1^2*x2^4*x3^2+3*x1^2*x2^4*x3^2*x4+9*x1^2*x2^4*x3^2*x4^2+8*x1^2*x2^4*x3^2*x4^4+3*x1^2*x2^4*x3^2*x4^5+4*x1^2*x2^4*x3^2*x4^6+5*x1^2*x2^4*x3^2*x4^7+12*x1^2*x2^4*x3^3+12*x1^2*x2^4*x3^3*x4+5*x1^2*x2^4*x3^3*x4^2+12*x1^2*x2^4*x3^3*x4^3+3*x1^2*x2^4*x3^3*x4^4+14*x1^2*x2^4*x3^3*x4^5+x1^2*x2^4*x3^3*x4^6+8*x1^2*x2^4*x3^4+7*x1^2*x2^4*x3^4*x4+11*x1^2*x2^4*x3^4*x4^2+7*x1^2*x2^4*x3^4*x4^3+11*x1^2*x2^4*x3^4*x4^4+5*x1^2*x2^4*x3^4*x4^5+14*x1^2*x2^4*x3^4*x4^6+3*x1^2*x2^4*x3^5+x1^2*x2^4*x3^5*x4+7*x1^2*x2^4*x3^5*x4^4+12*x1^2*x2^4*x3^5*x4^6+9*x1^2*x2^4*x3^6*x4+3*x1^2*x2^4*x3^6*x4^3+2*x1^2*x2^4*x3^6*x4^4+16*x1^2*x2^4*x3^7*x4+9*x1^2*x2^4*x3^7*x4^4+3*x1^2*x2^5+x1^2*x2^5*x4^2+8*x1^2*x2^5*x4^4+14*x1^2*x2^5*x4^6+15*x1^2*x2^5*x3+13*x1^2*x2^5*x3*x4^2+5*x1^2*x2^5*x3*x4^3+16*x1^2*x2^5*x3*x4^4+13*x1^2*x2^5*x3*x4^5+13*x1^2*x2^5*x3*x4^7+x1^2*x2^5*x3^2+3*x1^2*x2^5*x3^2*x4+2*x1^2*x2^5*x3^2*x4^2+8*x1^2*x2^5*x3^2*x4^3+6*x1^2*x2^5*x3^2*x4^4+13*x1^2*x2^5*x3^2*x4^5+6*x1^2*x2^5*x3^2*x4^6+16*x1^2*x2^5*x3^2*x4^7+x1^2*x2^5*x3^3+14*x1^2*x2^5*x3^3*x4+16*x1^2*x2^5*x3^3*x4^2+14*x1^2*x2^5*x3^3*x4^3+3*x1^2*x2^5*x3^3*x4^4+4*x1^2*x2^5*x3^3*x4^5+10*x1^2*x2^5*x3^3*x4^6+4*x1^2*x2^5*x3^4+x1^2*x2^5*x3^4*x4+8*x1^2*x2^5*x3^4*x4^2+2*x1^2*x2^5*x3^4*x4^3+4*x1^2*x2^5*x3^4*x4^4+16*x1^2*x2^5*x3^4*x4^5+4*x1^2*x2^5*x3^4*x4^6+13*x1^2*x2^5*x3^5+x1^2*x2^5*x3^5*x4^4+x1^2*x2^5*x3^5*x4^6+5*x1^2*x2^5*x3^6*x4+13*x1^2*x2^5*x3^6*x4^3+3*x1^2*x2^5*x3^6*x4^4+7*x1^2*x2^5*x3^7*x4+5*x1^2*x2^5*x3^7*x4^4+6*x1^2*x2^6+13*x1^2*x2^6*x4^2+x1^2*x2^6*x4^6+9*x1^2*x2^6*x3+16*x1^2*x2^6*x3*x4+10*x1^2*x2^6*x3*x4^2+13*x1^2*x2^6*x3*x4^3+6*x1^2*x2^6*x3*x4^4+2*x1^2*x2^6*x3*x4^5+7*x1^2*x2^6*x3*x4^7+10*x1^2*x2^6*x3^2+2*x1^2*x2^6*x3^2*x4+6*x1^2*x2^6*x3^2*x4^2+10*x1^2*x2^6*x3^2*x4^4+10*x1^2*x2^6*x3^2*x4^5+15*x1^2*x2^6*x3^2*x4^6+6*x1^2*x2^6*x3^2*x4^7+12*x1^2*x2^6*x3^3+3*x1^2*x2^6*x3^3*x4+6*x1^2*x2^6*x3^3*x4^2+2*x1^2*x2^6*x3^3*x4^3+2*x1^2*x2^6*x3^3*x4^4+10*x1^2*x2^6*x3^3*x4^5+8*x1^2*x2^6*x3^3*x4^6+4*x1^2*x2^6*x3^4+13*x1^2*x2^6*x3^4*x4+3*x1^2*x2^6*x3^4*x4^2+6*x1^2*x2^6*x3^4*x4^3+x1^2*x2^6*x3^4*x4^4+6*x1^2*x2^6*x3^4*x4^5+10*x1^2*x2^6*x3^4*x4^6+7*x1^2*x2^6*x3^5+15*x1^2*x2^6*x3^5*x4+13*x1^2*x2^6*x3^5*x4^4+11*x1^2*x2^6*x3^5*x4^6+4*x1^2*x2^6*x3^6*x4+7*x1^2*x2^6*x3^6*x4^3+16*x1^2*x2^6*x3^6*x4^4+9*x1^2*x2^6*x3^7*x4+4*x1^2*x2^6*x3^7*x4^4+10*x1^2*x2^7+10*x1^2*x2^7*x4^2+5*x1^2*x2^7*x4^4+3*x1^2*x2^7*x4^6+4*x1^2*x2^7*x3+9*x1^2*x2^7*x3*x4+15*x1^2*x2^7*x3*x4^2+13*x1^2*x2^7*x3*x4^3+x1^2*x2^7*x3*x4^4+14*x1^2*x2^7*x3*x4^5+4*x1^2*x2^7*x3*x4^7+15*x1^2*x2^7*x3^2+8*x1^2*x2^7*x3^2*x4+6*x1^2*x2^7*x3^2*x4^2+11*x1^2*x2^7*x3^2*x4^3+4*x1^2*x2^7*x3^2*x4^4+15*x1^2*x2^7*x3^2*x4^5+11*x1^2*x2^7*x3^2*x4^6+x1^2*x2^7*x3^2*x4^7+5*x1^2*x2^7*x3^3+16*x1^2*x2^7*x3^3*x4+x1^2*x2^7*x3^3*x4^2+11*x1^2*x2^7*x3^3*x4^3+8*x1^2*x2^7*x3^3*x4^4+13*x1^2*x2^7*x3^3*x4^5+7*x1^2*x2^7*x3^3*x4^6+12*x1^2*x2^7*x3^4*x4+9*x1^2*x2^7*x3^4*x4^2+4*x1^2*x2^7*x3^4*x4^3+14*x1^2*x2^7*x3^4*x4^4+x1^2*x2^7*x3^4*x4^5+13*x1^2*x2^7*x3^4*x4^6+4*x1^2*x2^7*x3^5+16*x1^2*x2^7*x3^5*x4+12*x1^2*x2^7*x3^5*x4^4+16*x1^2*x2^7*x3^5*x4^6+12*x1^2*x2^7*x3^6*x4+4*x1^2*x2^7*x3^6*x4^3+14*x1^2*x2^7*x3^6*x4^4+10*x1^2*x2^7*x3^7*x4+12*x1^2*x2^7*x3^7*x4^4+16*x1^2*x2^8+5*x1^2*x2^8*x4^2+x1^2*x2^8*x4^4+10*x1^2*x2^8*x4^6+9*x1^2*x2^8*x3+9*x1^2*x2^8*x3*x4+8*x1^2*x2^8*x3*x4^2+10*x1^2*x2^8*x3*x4^3+9*x1^2*x2^8*x3*x4^4+9*x1^2*x2^8*x3*x4^5+2*x1^2*x2^8*x3*x4^7+x1^2*x2^8*x3^2+11*x1^2*x2^8*x3^2*x4+7*x1^2*x2^8*x3^2*x4^2+4*x1^2*x2^8*x3^2*x4^3+13*x1^2*x2^8*x3^2*x4^4+8*x1^2*x2^8*x3^2*x4^5+14*x1^2*x2^8*x3^2*x4^6+9*x1^2*x2^8*x3^2*x4^7+12*x1^2*x2^8*x3^3+2*x1^2*x2^8*x3^3*x4+9*x1^2*x2^8*x3^3*x4^2+13*x1^2*x2^8*x3^3*x4^3+11*x1^2*x2^8*x3^3*x4^4+15*x1^2*x2^8*x3^3*x4^5+12*x1^2*x2^8*x3^3*x4^6+9*x1^2*x2^8*x3^4+5*x1^2*x2^8*x3^4*x4+13*x1^2*x2^8*x3^4*x4^2+3*x1^2*x2^8*x3^4*x4^3+3*x1^2*x2^8*x3^4*x4^4+9*x1^2*x2^8*x3^4*x4^5+15*x1^2*x2^8*x3^4*x4^6+2*x1^2*x2^8*x3^5+2*x1^2*x2^8*x3^5*x4+5*x1^2*x2^8*x3^5*x4^4+8*x1^2*x2^8*x3^5*x4^6+6*x1^2*x2^8*x3^6*x4+2*x1^2*x2^8*x3^6*x4^3+7*x1^2*x2^8*x3^6*x4^4+5*x1^2*x2^8*x3^7*x4+6*x1^2*x2^8*x3^7*x4^4+2*x1^2*x2^9+15*x1^2*x2^9*x4^2+x1^2*x2^9*x4^4+5*x1^2*x2^9*x3+6*x1^2*x2^9*x3*x4+16*x1^2*x2^9*x3*x4^2+15*x1^2*x2^9*x3*x4^3+13*x1^2*x2^9*x3*x4^5+7*x1^2*x2^9*x3^2+11*x1^2*x2^9*x3^2*x4+11*x1^2*x2^9*x3^2*x4^2+5*x1^2*x2^9*x3^2*x4^3+10*x1^2*x2^9*x3^2*x4^4+16*x1^2*x2^9*x3^2*x4^5+4*x1^2*x2^9*x3^3+9*x1^2*x2^9*x3^3*x4+7*x1^2*x2^9*x3^3*x4^3+11*x1^2*x2^9*x3^3*x4^4+12*x1^2*x2^9*x3^4+13*x1^2*x2^9*x3^4*x4+4*x1^2*x2^9*x3^4*x4^3+x1^2*x2^9*x3^4*x4^4+9*x1^2*x2^9*x3^5*x4+13*x1^2*x2^9*x3^5*x4^4+11*x1^2*x2^10+2*x1^2*x2^10*x4^2+4*x1^2*x2^10*x4^4+2*x1^2*x2^10*x4^6+15*x1^2*x2^10*x3+12*x1^2*x2^10*x3*x4+4*x1^2*x2^10*x3*x4^2+12*x1^2*x2^10*x3*x4^3+12*x1^2*x2^10*x3*x4^4+8*x1^2*x2^10*x3*x4^5+14*x1^2*x2^10*x3*x4^7+8*x1^2*x2^10*x3^2+11*x1^2*x2^10*x3^2*x4^2+14*x1^2*x2^10*x3^2*x4^3+4*x1^2*x2^10*x3^2*x4^4+4*x1^2*x2^10*x3^2*x4^5+13*x1^2*x2^10*x3^2*x4^6+12*x1^2*x2^10*x3^2*x4^7+8*x1^2*x2^10*x3^3+2*x1^2*x2^10*x3^3*x4+12*x1^2*x2^10*x3^3*x4^2+x1^2*x2^10*x3^3*x4^3+3*x1^2*x2^10*x3^3*x4^5+16*x1^2*x2^10*x3^3*x4^6+10*x1^2*x2^10*x3^4+12*x1^2*x2^10*x3^4*x4+6*x1^2*x2^10*x3^4*x4^2+13*x1^2*x2^10*x3^4*x4^3+14*x1^2*x2^10*x3^4*x4^4+12*x1^2*x2^10*x3^4*x4^5+3*x1^2*x2^10*x3^4*x4^6+14*x1^2*x2^10*x3^5+6*x1^2*x2^10*x3^5*x4+12*x1^2*x2^10*x3^5*x4^4+5*x1^2*x2^10*x3^5*x4^6+8*x1^2*x2^10*x3^6*x4+14*x1^2*x2^10*x3^6*x4^3+15*x1^2*x2^10*x3^6*x4^4+x1^2*x2^10*x3^7*x4+8*x1^2*x2^10*x3^7*x4^4+16*x1^3+3*x1^3*x4^2+15*x1^3*x4^4+3*x1^3*x4^6+8*x1^3*x3+2*x1^3*x3*x4+5*x1^3*x3*x4^2+2*x1^3*x3*x4^3+4*x1^3*x3*x4^4+8*x1^3*x3*x4^5+16*x1^3*x3*x4^7+3*x1^3*x3^2+9*x1^3*x3^2*x4+12*x1^3*x3^2*x4^2+2*x1^3*x3^2*x4^3+8*x1^3*x3^2*x4^4+5*x1^3*x3^2*x4^5+10*x1^3*x3^2*x4^6+4*x1^3*x3^2*x4^7+4*x1^3*x3^3+2*x1^3*x3^3*x4+12*x1^3*x3^3*x4^2+6*x1^3*x3^3*x4^3+9*x1^3*x3^3*x4^4+6*x1^3*x3^3*x4^5+11*x1^3*x3^3*x4^6+14*x1^3*x3^4+12*x1^3*x3^4*x4+9*x1^3*x3^4*x4^2+14*x1^3*x3^4*x4^4+12*x1^3*x3^4*x4^5+x1^3*x3^4*x4^6+16*x1^3*x3^5+6*x1^3*x3^5*x4+14*x1^3*x3^5*x4^3+12*x1^3*x3^5*x4^4+13*x1^3*x3^5*x4^6+4*x1^3*x3^6*x4+16*x1^3*x3^6*x4^3+16*x1^3*x3^6*x4^4+10*x1^3*x3^7*x4+4*x1^3*x3^7*x4^4+7*x1^3*x2+13*x1^3*x2*x4^2+6*x1^3*x2*x4^4+x1^3*x2*x4^6+11*x1^3*x2*x3+10*x1^3*x2*x3*x4+10*x1^3*x2*x3*x4^2+7*x1^3*x2*x3*x4^3+7*x1^3*x2*x3*x4^4+2*x1^3*x2*x3*x4^5+11*x1^3*x2*x3*x4^7+10*x1^3*x2*x3^2+10*x1^3*x2*x3^2*x4^2+9*x1^3*x2*x3^2*x4^3+2*x1^3*x2*x3^2*x4^4+10*x1^3*x2*x3^2*x4^5+9*x1^3*x2*x3^2*x4^6+7*x1^3*x2*x3^2*x4^7+9*x1^3*x2*x3^3+3*x1^3*x2*x3^3*x4+4*x1^3*x2*x3^3*x4^2+9*x1^3*x2*x3^3*x4^3+2*x1^3*x2*x3^3*x4^5+15*x1^3*x2*x3^3*x4^6+3*x1^3*x2*x3^4+2*x1^3*x2*x3^4*x4+3*x1^3*x2*x3^4*x4^2+2*x1^3*x2*x3^4*x4^3+8*x1^3*x2*x3^4*x4^4+4*x1^3*x2*x3^4*x4^5+6*x1^3*x2*x3^4*x4^6+11*x1^3*x2*x3^5+12*x1^3*x2*x3^5*x4+16*x1^3*x2*x3^5*x4^3+2*x1^3*x2*x3^5*x4^4+10*x1^3*x2*x3^5*x4^6+7*x1^3*x2*x3^6*x4+11*x1^3*x2*x3^6*x4^3+11*x1^3*x2*x3^6*x4^4+9*x1^3*x2*x3^7*x4+7*x1^3*x2*x3^7*x4^4+4*x1^3*x2^2+15*x1^3*x2^2*x4^2+10*x1^3*x2^2*x4^4+9*x1^3*x2^2*x4^6+9*x1^3*x2^2*x3+9*x1^3*x2^2*x3*x4+14*x1^3*x2^2*x3*x4^2+10*x1^3*x2^2*x3*x4^3+12*x1^3*x2^2*x3*x4^4+3*x1^3*x2^2*x3*x4^5+14*x1^3*x2^2*x3*x4^7+x1^3*x2^2*x3^2+6*x1^3*x2^2*x3^2*x4+8*x1^3*x2^2*x3^2*x4^2+8*x1^3*x2^2*x3^2*x4^3+8*x1^3*x2^2*x3^2*x4^4+14*x1^3*x2^2*x3^2*x4^5+13*x1^3*x2^2*x3^2*x4^6+12*x1^3*x2^2*x3^2*x4^7+9*x1^3*x2^2*x3^3+2*x1^3*x2^2*x3^3*x4+2*x1^3*x2^2*x3^3*x4^2+9*x1^3*x2^2*x3^3*x4^3+6*x1^3*x2^2*x3^3*x4^4+x1^3*x2^2*x3^3*x4^5+16*x1^3*x2^2*x3^3*x4^6+2*x1^3*x2^2*x3^4+5*x1^3*x2^2*x3^4*x4+10*x1^3*x2^2*x3^4*x4^2+14*x1^3*x2^2*x3^4*x4^3+3*x1^3*x2^2*x3^4*x4^4+2*x1^3*x2^2*x3^4*x4^5+3*x1^3*x2^2*x3^4*x4^6+14*x1^3*x2^2*x3^5+14*x1^3*x2^2*x3^5*x4+8*x1^3*x2^2*x3^5*x4^3+5*x1^3*x2^2*x3^5*x4^4+5*x1^3*x2^2*x3^5*x4^6+12*x1^3*x2^2*x3^6*x4+14*x1^3*x2^2*x3^6*x4^3+14*x1^3*x2^2*x3^6*x4^4+13*x1^3*x2^2*x3^7*x4+12*x1^3*x2^2*x3^7*x4^4+14*x1^3*x2^3+16*x1^3*x2^3*x4^2+9*x1^3*x2^3*x4^4+2*x1^3*x2^3*x4^6+13*x1^3*x2^3*x3+14*x1^3*x2^3*x3*x4+6*x1^3*x2^3*x3*x4^2+11*x1^3*x2^3*x3*x4^3+14*x1^3*x2^3*x3*x4^4+16*x1^3*x2^3*x3*x4^5+5*x1^3*x2^3*x3*x4^7+10*x1^3*x2^3*x3^2+7*x1^3*x2^3*x3^2*x4+10*x1^3*x2^3*x3^2*x4^2+9*x1^3*x2^3*x3^2*x4^3+15*x1^3*x2^3*x3^2*x4^4+6*x1^3*x2^3*x3^2*x4^5+x1^3*x2^3*x3^2*x4^6+14*x1^3*x2^3*x3^2*x4^7+7*x1^3*x2^3*x3^3+16*x1^3*x2^3*x3^3*x4+8*x1^3*x2^3*x3^3*x4^2+13*x1^3*x2^3*x3^3*x4^3+7*x1^3*x2^3*x3^3*x4^4+4*x1^3*x2^3*x3^3*x4^5+13*x1^3*x2^3*x3^3*x4^6+9*x1^3*x2^3*x3^4+3*x1^3*x2^3*x3^4*x4+6*x1^3*x2^3*x3^4*x4^2+10*x1^3*x2^3*x3^4*x4^3+12*x1^3*x2^3*x3^4*x4^4+8*x1^3*x2^3*x3^4*x4^5+12*x1^3*x2^3*x3^4*x4^6+5*x1^3*x2^3*x3^5+5*x1^3*x2^3*x3^5*x4+15*x1^3*x2^3*x3^5*x4^3+3*x1^3*x2^3*x3^5*x4^4+3*x1^3*x2^3*x3^5*x4^6+14*x1^3*x2^3*x3^6*x4+5*x1^3*x2^3*x3^6*x4^3+5*x1^3*x2^3*x3^6*x4^4+x1^3*x2^3*x3^7*x4+14*x1^3*x2^3*x3^7*x4^4+4*x1^3*x2^4+13*x1^3*x2^4*x4^2+4*x1^3*x2^4*x4^4+8*x1^3*x2^4*x4^6+14*x1^3*x2^4*x3+x1^3*x2^4*x3*x4^2+x1^3*x2^4*x3*x4^3+5*x1^3*x2^4*x3*x4^4+6*x1^3*x2^4*x3*x4^5+3*x1^3*x2^4*x3*x4^7+15*x1^3*x2^4*x3^2+12*x1^3*x2^4*x3^2*x4+5*x1^3*x2^4*x3^2*x4^2+15*x1^3*x2^4*x3^2*x4^3+13*x1^3*x2^4*x3^2*x4^4+x1^3*x2^4*x3^2*x4^5+4*x1^3*x2^4*x3^2*x4^6+5*x1^3*x2^4*x3^2*x4^7+16*x1^3*x2^4*x3^3+4*x1^3*x2^4*x3^3*x4+15*x1^3*x2^4*x3^3*x4^2+3*x1^3*x2^4*x3^3*x4^3+12*x1^3*x2^4*x3^3*x4^4+16*x1^3*x2^4*x3^3*x4^5+x1^3*x2^4*x3^3*x4^6+9*x1^3*x2^4*x3^4+7*x1^3*x2^4*x3^4*x4+7*x1^3*x2^4*x3^4*x4^2+11*x1^3*x2^4*x3^4*x4^3+11*x1^3*x2^4*x3^4*x4^4+15*x1^3*x2^4*x3^4*x4^5+14*x1^3*x2^4*x3^4*x4^6+3*x1^3*x2^4*x3^5+3*x1^3*x2^4*x3^5*x4+9*x1^3*x2^4*x3^5*x4^3+7*x1^3*x2^4*x3^5*x4^4+12*x1^3*x2^4*x3^5*x4^6+5*x1^3*x2^4*x3^6*x4+3*x1^3*x2^4*x3^6*x4^3+3*x1^3*x2^4*x3^6*x4^4+4*x1^3*x2^4*x3^7*x4+5*x1^3*x2^4*x3^7*x4^4+2*x1^3*x2^5+5*x1^3*x2^5*x4^4+12*x1^3*x2^5*x4^6+5*x1^3*x2^5*x3+11*x1^3*x2^5*x3*x4+4*x1^3*x2^5*x3*x4^2+4*x1^3*x2^5*x3*x4^3+16*x1^3*x2^5*x3*x4^4+2*x1^3*x2^5*x3*x4^5+13*x1^3*x2^5*x3*x4^7+15*x1^3*x2^5*x3^2*x4+11*x1^3*x2^5*x3^2*x4^2+15*x1^3*x2^5*x3^2*x4^3+16*x1^3*x2^5*x3^2*x4^4+4*x1^3*x2^5*x3^2*x4^5+6*x1^3*x2^5*x3^2*x4^6+16*x1^3*x2^5*x3^2*x4^7+9*x1^3*x2^5*x3^3+11*x1^3*x2^5*x3^3*x4+14*x1^3*x2^5*x3^3*x4^2+3*x1^3*x2^5*x3^3*x4^3+15*x1^3*x2^5*x3^3*x4^4+7*x1^3*x2^5*x3^3*x4^5+10*x1^3*x2^5*x3^3*x4^6+3*x1^3*x2^5*x3^4+2*x1^3*x2^5*x3^4*x4^2+10*x1^3*x2^5*x3^4*x4^3+14*x1^3*x2^5*x3^4*x4^5+4*x1^3*x2^5*x3^4*x4^6+13*x1^3*x2^5*x3^5+2*x1^3*x2^5*x3^5*x4+5*x1^3*x2^5*x3^5*x4^3+x1^3*x2^5*x3^5*x4^6+16*x1^3*x2^5*x3^6*x4+13*x1^3*x2^5*x3^6*x4^3+13*x1^3*x2^5*x3^6*x4^4+6*x1^3*x2^5*x3^7*x4+16*x1^3*x2^5*x3^7*x4^4+13*x1^3*x2^6+4*x1^3*x2^6*x4^2+6*x1^3*x2^6*x4^4+13*x1^3*x2^6*x4^6+6*x1^3*x2^6*x3+4*x1^3*x2^6*x3*x4+6*x1^3*x2^6*x3*x4^2+6*x1^3*x2^6*x3*x4^4+6*x1^3*x2^6*x3*x4^5+7*x1^3*x2^6*x3*x4^7+9*x1^3*x2^6*x3^2*x4+x1^3*x2^6*x3^2*x4^2+14*x1^3*x2^6*x3^2*x4^3+11*x1^3*x2^6*x3^2*x4^4+6*x1^3*x2^6*x3^2*x4^5+15*x1^3*x2^6*x3^2*x4^6+6*x1^3*x2^6*x3^2*x4^7+6*x1^3*x2^6*x3^3+10*x1^3*x2^6*x3^3*x4+x1^3*x2^6*x3^3*x4^2+6*x1^3*x2^6*x3^3*x4^3+9*x1^3*x2^6*x3^3*x4^4+9*x1^3*x2^6*x3^3*x4^5+8*x1^3*x2^6*x3^3*x4^6+15*x1^3*x2^6*x3^4*x4+5*x1^3*x2^6*x3^4*x4^2+9*x1^3*x2^6*x3^4*x4^4+x1^3*x2^6*x3^4*x4^5+10*x1^3*x2^6*x3^4*x4^6+7*x1^3*x2^6*x3^5+4*x1^3*x2^6*x3^5*x4^3+15*x1^3*x2^6*x3^5*x4^4+11*x1^3*x2^6*x3^5*x4^6+6*x1^3*x2^6*x3^6*x4+7*x1^3*x2^6*x3^6*x4^3+7*x1^3*x2^6*x3^6*x4^4+15*x1^3*x2^6*x3^7*x4+6*x1^3*x2^6*x3^7*x4^4+14*x1^3*x2^7+10*x1^3*x2^7*x4^2+13*x1^3*x2^7*x4^4+5*x1^3*x2^7*x4^6+6*x1^3*x2^7*x3+9*x1^3*x2^7*x3*x4+10*x1^3*x2^7*x3*x4^2+14*x1^3*x2^7*x3*x4^3+x1^3*x2^7*x3*x4^4+3*x1^3*x2^7*x3*x4^5+4*x1^3*x2^7*x3*x4^7+16*x1^3*x2^7*x3^2+x1^3*x2^7*x3^2*x4+4*x1^3*x2^7*x3^2*x4^2+13*x1^3*x2^7*x3^2*x4^3+14*x1^3*x2^7*x3^2*x4^4+10*x1^3*x2^7*x3^2*x4^5+11*x1^3*x2^7*x3^2*x4^6+x1^3*x2^7*x3^2*x4^7+10*x1^3*x2^7*x3^3+6*x1^3*x2^7*x3^3*x4+3*x1^3*x2^7*x3^3*x4^2+16*x1^3*x2^7*x3^3*x4^3+x1^3*x2^7*x3^3*x4^4+10*x1^3*x2^7*x3^3*x4^5+7*x1^3*x2^7*x3^3*x4^6+6*x1^3*x2^7*x3^4+5*x1^3*x2^7*x3^4*x4+15*x1^3*x2^7*x3^4*x4^2+9*x1^3*x2^7*x3^4*x4^3+3*x1^3*x2^7*x3^4*x4^4+3*x1^3*x2^7*x3^4*x4^5+13*x1^3*x2^7*x3^4*x4^6+4*x1^3*x2^7*x3^5+8*x1^3*x2^7*x3^5*x4+12*x1^3*x2^7*x3^5*x4^3+5*x1^3*x2^7*x3^5*x4^4+16*x1^3*x2^7*x3^5*x4^6+x1^3*x2^7*x3^6*x4+4*x1^3*x2^7*x3^6*x4^3+4*x1^3*x2^7*x3^6*x4^4+11*x1^3*x2^7*x3^7*x4+x1^3*x2^7*x3^7*x4^4+13*x1^3*x2^8+6*x1^3*x2^8*x4^2+16*x1^3*x2^8*x4^4+11*x1^3*x2^8*x4^6+7*x1^3*x2^8*x3+16*x1^3*x2^8*x3*x4+9*x1^3*x2^8*x3*x4^2+x1^3*x2^8*x3*x4^3+9*x1^3*x2^8*x3*x4^4+9*x1^3*x2^8*x3*x4^5+2*x1^3*x2^8*x3*x4^7+6*x1^3*x2^8*x3^2+2*x1^3*x2^8*x3^2*x4+13*x1^3*x2^8*x3^2*x4^2+2*x1^3*x2^8*x3^2*x4^3+13*x1^3*x2^8*x3^2*x4^4+9*x1^3*x2^8*x3^2*x4^5+14*x1^3*x2^8*x3^2*x4^6+9*x1^3*x2^8*x3^2*x4^7+7*x1^3*x2^8*x3^3+10*x1^3*x2^8*x3^3*x4+10*x1^3*x2^8*x3^3*x4^2+7*x1^3*x2^8*x3^3*x4^3+2*x1^3*x2^8*x3^3*x4^4+5*x1^3*x2^8*x3^3*x4^5+12*x1^3*x2^8*x3^3*x4^6+15*x1^3*x2^8*x3^4+12*x1^3*x2^8*x3^4*x4+16*x1^3*x2^8*x3^4*x4^2+15*x1^3*x2^8*x3^4*x4^3+14*x1^3*x2^8*x3^4*x4^4+10*x1^3*x2^8*x3^4*x4^5+15*x1^3*x2^8*x3^4*x4^6+2*x1^3*x2^8*x3^5+13*x1^3*x2^8*x3^5*x4+6*x1^3*x2^8*x3^5*x4^3+12*x1^3*x2^8*x3^5*x4^4+8*x1^3*x2^8*x3^5*x4^6+9*x1^3*x2^8*x3^6*x4+2*x1^3*x2^8*x3^6*x4^3+2*x1^3*x2^8*x3^6*x4^4+14*x1^3*x2^8*x3^7*x4+9*x1^3*x2^8*x3^7*x4^4+4*x1^3*x2^9+2*x1^3*x2^9*x4^2+4*x1^3*x2^9*x4^4+11*x1^3*x2^9*x3+7*x1^3*x2^9*x3*x4+3*x1^3*x2^9*x3*x4^2+x1^3*x2^9*x3*x4^3+12*x1^3*x2^9*x3*x4^5+14*x1^3*x2^9*x3^2+12*x1^3*x2^9*x3^2*x4+8*x1^3*x2^9*x3^2*x4^2+11*x1^3*x2^9*x3^2*x4^3+14*x1^3*x2^9*x3^2*x4^4+3*x1^3*x2^9*x3^2*x4^5+6*x1^3*x2^9*x3^3+15*x1^3*x2^9*x3^3*x4+x1^3*x2^9*x3^3*x4^3+12*x1^3*x2^9*x3^3*x4^4+x1^3*x2^9*x3^4+8*x1^3*x2^9*x3^4*x4+6*x1^3*x2^9*x3^4*x4^3+15*x1^3*x2^9*x3^4*x4^4+5*x1^3*x2^9*x3^5*x4+8*x1^3*x2^9*x3^5*x4^4+5*x1^3*x2^10+x1^3*x2^10*x4^2+2*x1^3*x2^10*x4^4+9*x1^3*x2^10*x4^6+2*x1^3*x2^10*x3+5*x1^3*x2^10*x3*x4+11*x1^3*x2^10*x3*x4^2+10*x1^3*x2^10*x3*x4^3+12*x1^3*x2^10*x3*x4^4+8*x1^3*x2^10*x3*x4^5+14*x1^3*x2^10*x3*x4^7+3*x1^3*x2^10*x3^2+5*x1^3*x2^10*x3^2*x4+11*x1^3*x2^10*x3^2*x4^2+x1^3*x2^10*x3^2*x4^3+4*x1^3*x2^10*x3^2*x4^4+11*x1^3*x2^10*x3^2*x4^5+13*x1^3*x2^10*x3^2*x4^6+12*x1^3*x2^10*x3^2*x4^7+12*x1^3*x2^10*x3^3+2*x1^3*x2^10*x3^3*x4+2*x1^3*x2^10*x3^3*x4^2+5*x1^3*x2^10*x3^3*x4^4+x1^3*x2^10*x3^3*x4^5+16*x1^3*x2^10*x3^3*x4^6+16*x1^3*x2^10*x3^4+10*x1^3*x2^10*x3^4*x4+10*x1^3*x2^10*x3^4*x4^2+6*x1^3*x2^10*x3^4*x4^4+2*x1^3*x2^10*x3^4*x4^5+3*x1^3*x2^10*x3^4*x4^6+14*x1^3*x2^10*x3^5+12*x1^3*x2^10*x3^5*x4+8*x1^3*x2^10*x3^5*x4^3+10*x1^3*x2^10*x3^5*x4^4+5*x1^3*x2^10*x3^5*x4^6+12*x1^3*x2^10*x3^6*x4+14*x1^3*x2^10*x3^6*x4^3+14*x1^3*x2^10*x3^6*x4^4+13*x1^3*x2^10*x3^7*x4+12*x1^3*x2^10*x3^7*x4^4+13*x1^4+2*x1^4*x4^2+8*x1^4*x4^4+2*x1^4*x3+14*x1^4*x3*x4+16*x1^4*x3*x4^2+9*x1^4*x3*x4^3+13*x1^4*x3*x4^5+11*x1^4*x3^2+16*x1^4*x3^2*x4+11*x1^4*x3^2*x4^2+2*x1^4*x3^2*x4^3+13*x1^4*x3^2*x4^4+16*x1^4*x3^2*x4^5+15*x1^4*x3^3+3*x1^4*x3^3*x4+4*x1^4*x3^3*x4^3+16*x1^4*x3^3*x4^4+12*x1^4*x3^4+5*x1^4*x3^4*x4+15*x1^4*x3^4*x4^3+3*x1^4*x3^4*x4^4+9*x1^4*x3^5*x4+5*x1^4*x3^5*x4^4+15*x1^4*x2+15*x1^4*x2*x4^2+14*x1^4*x2*x4^4+10*x1^4*x2*x3*x4+11*x1^4*x2*x3*x4^2+6*x1^4*x2*x3*x4^3+10*x1^4*x2*x3*x4^5+2*x1^4*x2*x3^2+11*x1^4*x2*x3^2*x4+15*x1^4*x2*x3^2*x4^2+10*x1^4*x2*x3^2*x4^4+11*x1^4*x2*x3^2*x4^5+5*x1^4*x2*x3^3+14*x1^4*x2*x3^3*x4+7*x1^4*x2*x3^3*x4^3+11*x1^4*x2*x3^3*x4^4+4*x1^4*x2*x3^4+13*x1^4*x2*x3^4*x4+5*x1^4*x2*x3^4*x4^3+x1^4*x2*x3^4*x4^4+3*x1^4*x2*x3^5*x4+13*x1^4*x2*x3^5*x4^4+13*x1^4*x2^2+13*x1^4*x2^2*x4^2+7*x1^4*x2^2*x4^4+8*x1^4*x2^2*x3+9*x1^4*x2^2*x3*x4+14*x1^4*x2^2*x3*x4^2+x1^4*x2^2*x3*x4^3+5*x1^4*x2^2*x3*x4^5+14*x1^4*x2^2*x3^2*x4+16*x1^4*x2^2*x3^2*x4^2+8*x1^4*x2^2*x3^2*x4^3+5*x1^4*x2^2*x3^2*x4^4+14*x1^4*x2^2*x3^2*x4^5+11*x1^4*x2^2*x3^3+4*x1^4*x2^2*x3^3*x4+12*x1^4*x2^2*x3^3*x4^3+14*x1^4*x2^2*x3^3*x4^4+2*x1^4*x2^2*x3^4+15*x1^4*x2^2*x3^4*x4+11*x1^4*x2^2*x3^4*x4^3+9*x1^4*x2^2*x3^4*x4^4+10*x1^4*x2^2*x3^5*x4+15*x1^4*x2^2*x3^5*x4^4+12*x1^4*x2^3+2*x1^4*x2^3*x4^2+11*x1^4*x2^3*x4^4+2*x1^4*x2^3*x3+5*x1^4*x2^3*x3*x4^2+3*x1^4*x2^3*x3*x4^3+3*x1^4*x2^3*x3*x4^5+15*x1^4*x2^3*x3^2+5*x1^4*x2^3*x3^2*x4+13*x1^4*x2^3*x3^2*x4^2+2*x1^4*x2^3*x3^2*x4^3+3*x1^4*x2^3*x3^2*x4^4+5*x1^4*x2^3*x3^2*x4^5+10*x1^4*x2^3*x3^3+9*x1^4*x2^3*x3^3*x4+14*x1^4*x2^3*x3^3*x4^3+5*x1^4*x2^3*x3^3*x4^4+8*x1^4*x2^3*x3^4+9*x1^4*x2^3*x3^4*x4+10*x1^4*x2^3*x3^4*x4^3+2*x1^4*x2^3*x3^4*x4^4+6*x1^4*x2^3*x3^5*x4+9*x1^4*x2^3*x3^5*x4^4+12*x1^4*x2^4+14*x1^4*x2^4*x4^2+10*x1^4*x2^4*x4^4+4*x1^4*x2^4*x3+12*x1^4*x2^4*x3*x4+3*x1^4*x2^4*x3*x4^2+13*x1^4*x2^4*x3*x4^3+12*x1^4*x2^4*x3*x4^5+8*x1^4*x2^4*x3^2+3*x1^4*x2^4*x3^2*x4+x1^4*x2^4*x3^2*x4^2+4*x1^4*x2^4*x3^2*x4^3+12*x1^4*x2^4*x3^2*x4^4+3*x1^4*x2^4*x3^2*x4^5+6*x1^4*x2^4*x3^3+10*x1^4*x2^4*x3^3*x4+5*x1^4*x2^4*x3^3*x4^3+3*x1^4*x2^4*x3^3*x4^4+15*x1^4*x2^4*x3^4+2*x1^4*x2^4*x3^4*x4+6*x1^4*x2^4*x3^4*x4^3+8*x1^4*x2^4*x3^4*x4^4+7*x1^4*x2^4*x3^5*x4+2*x1^4*x2^4*x3^5*x4^4+4*x1^4*x2^5+15*x1^4*x2^5*x4^2+15*x1^4*x2^5*x4^4+10*x1^4*x2^5*x3+14*x1^4*x2^5*x3*x4+13*x1^4*x2^5*x3*x4^2+10*x1^4*x2^5*x3*x4^3+x1^4*x2^5*x3*x4^5+14*x1^4*x2^5*x3^2+13*x1^4*x2^5*x3^2*x4+10*x1^4*x2^5*x3^2*x4^2+10*x1^4*x2^5*x3^2*x4^3+x1^4*x2^5*x3^2*x4^4+13*x1^4*x2^5*x3^2*x4^5+9*x1^4*x2^5*x3^3+x1^4*x2^5*x3^3*x4+16*x1^4*x2^5*x3^3*x4^3+13*x1^4*x2^5*x3^3*x4^4+14*x1^4*x2^5*x3^4+3*x1^4*x2^5*x3^4*x4+9*x1^4*x2^5*x3^4*x4^3+12*x1^4*x2^5*x3^4*x4^4+2*x1^4*x2^5*x3^5*x4+3*x1^4*x2^5*x3^5*x4^4+8*x1^4*x2^6+13*x1^4*x2^6*x4^2+12*x1^4*x2^6*x4^4+11*x1^4*x2^6*x3+x1^4*x2^6*x3*x4+7*x1^4*x2^6*x3*x4^2+3*x1^4*x2^6*x3*x4^3+11*x1^4*x2^6*x3*x4^5+15*x1^4*x2^6*x3^2+7*x1^4*x2^6*x3^2*x4+8*x1^4*x2^6*x3^2*x4^2+11*x1^4*x2^6*x3^2*x4^3+11*x1^4*x2^6*x3^2*x4^4+7*x1^4*x2^6*x3^2*x4^5+14*x1^4*x2^6*x3^3+11*x1^4*x2^6*x3^3*x4+6*x1^4*x2^6*x3^3*x4^3+7*x1^4*x2^6*x3^3*x4^4+x1^4*x2^6*x3^4+16*x1^4*x2^6*x3^4*x4+14*x1^4*x2^6*x3^4*x4^3+13*x1^4*x2^6*x3^4*x4^4+5*x1^4*x2^6*x3^5*x4+16*x1^4*x2^6*x3^5*x4^4+16*x1^4*x2^7+5*x1^4*x2^7*x4^2+2*x1^4*x2^7*x4^4+12*x1^4*x2^7*x3+9*x1^4*x2^7*x3*x4+4*x1^4*x2^7*x3*x4^2+10*x1^4*x2^7*x3*x4^3+16*x1^4*x2^7*x3*x4^5+5*x1^4*x2^7*x3^2+4*x1^4*x2^7*x3^2*x4+7*x1^4*x2^7*x3^2*x4^2+12*x1^4*x2^7*x3^2*x4^3+16*x1^4*x2^7*x3^2*x4^4+4*x1^4*x2^7*x3^2*x4^5+8*x1^4*x2^7*x3^3+3*x1^4*x2^7*x3^3*x4+x1^4*x2^7*x3^3*x4^3+4*x1^4*x2^7*x3^3*x4^4+3*x1^4*x2^7*x3^4+14*x1^4*x2^7*x3^4*x4+8*x1^4*x2^7*x3^4*x4^3+5*x1^4*x2^7*x3^4*x4^4+15*x1^4*x2^7*x3^5*x4+14*x1^4*x2^7*x3^5*x4^4+5*x1^4*x2^8+5*x1^4*x2^8*x4^2+x1^4*x2^8*x4^4+2*x1^4*x2^8*x3+5*x1^4*x2^8*x3*x4+2*x1^4*x2^8*x3*x4^2+6*x1^4*x2^8*x3*x4^3+8*x1^4*x2^8*x3*x4^5+9*x1^4*x2^8*x3^2+2*x1^4*x2^8*x3^2*x4+12*x1^4*x2^8*x3^2*x4^2+2*x1^4*x2^8*x3^2*x4^3+8*x1^4*x2^8*x3^2*x4^4+2*x1^4*x2^8*x3^2*x4^5+4*x1^4*x2^8*x3^3+16*x1^4*x2^8*x3^3*x4+9*x1^4*x2^8*x3^3*x4^3+2*x1^4*x2^8*x3^3*x4^4+10*x1^4*x2^8*x3^4+7*x1^4*x2^8*x3^4*x4+4*x1^4*x2^8*x3^4*x4^3+11*x1^4*x2^8*x3^4*x4^4+16*x1^4*x2^8*x3^5*x4+7*x1^4*x2^8*x3^5*x4^4+10*x1^4*x2^9+2*x1^4*x2^9*x3+4*x1^4*x2^9*x3*x4+8*x1^4*x2^9*x3*x4^3+11*x1^4*x2^9*x3^2+2*x1^4*x2^9*x3^2*x4^3+14*x1^4*x2^9*x3^3*x4+9*x1^4*x2^10+14*x1^4*x2^10*x4^2+7*x1^4*x2^10*x4^4+16*x1^4*x2^10*x3+12*x1^4*x2^10*x3*x4+14*x1^4*x2^10*x3*x4^2+16*x1^4*x2^10*x3*x4^3+5*x1^4*x2^10*x3*x4^5+13*x1^4*x2^10*x3^2+14*x1^4*x2^10*x3^2*x4+16*x1^4*x2^10*x3^2*x4^2+16*x1^4*x2^10*x3^2*x4^3+5*x1^4*x2^10*x3^2*x4^4+14*x1^4*x2^10*x3^2*x4^5+11*x1^4*x2^10*x3^3+6*x1^4*x2^10*x3^3*x4+12*x1^4*x2^10*x3^3*x4^3+14*x1^4*x2^10*x3^3*x4^4+2*x1^4*x2^10*x3^4+15*x1^4*x2^10*x3^4*x4+11*x1^4*x2^10*x3^4*x4^3+9*x1^4*x2^10*x3^4*x4^4+10*x1^4*x2^10*x3^5*x4+15*x1^4*x2^10*x3^5*x4^4+16*x1^5+15*x1^5*x4^4+6*x1^5*x3+15*x1^5*x3*x4+9*x1^5*x3*x4^2+10*x1^5*x3*x4^3+2*x1^5*x3*x4^5+6*x1^5*x3^2+14*x1^5*x3^2*x4+5*x1^5*x3^2*x4^2+6*x1^5*x3^2*x4^3+5*x1^5*x3^2*x4^4+9*x1^5*x3^2*x4^5+8*x1^5*x3^3+14*x1^5*x3^3*x4+14*x1^5*x3^3*x4^4+9*x1^5*x3^4+15*x1^5*x3^4*x4+8*x1^5*x3^4*x4^3+9*x1^5*x3^4*x4^4+14*x1^5*x3^5*x4+15*x1^5*x3^5*x4^4+6*x1^5*x2+5*x1^5*x2*x4^2+5*x1^5*x2*x4^4+x1^5*x2*x3+12*x1^5*x2*x3*x4+3*x1^5*x2*x3*x4^2+5*x1^5*x2*x3*x4^3+12*x1^5*x2*x3*x4^5+14*x1^5*x2*x3^2+16*x1^5*x2*x3^2*x4+13*x1^5*x2*x3^2*x4^2+x1^5*x2*x3^2*x4^3+13*x1^5*x2*x3^2*x4^4+3*x1^5*x2*x3^2*x4^5+14*x1^5*x2*x3^3+15*x1^5*x2*x3^3*x4+16*x1^5*x2*x3^3*x4^4+3*x1^5*x2*x3^4+5*x1^5*x2*x3^4*x4+14*x1^5*x2*x3^4*x4^3+3*x1^5*x2*x3^4*x4^4+16*x1^5*x2*x3^5*x4+5*x1^5*x2*x3^5*x4^4+9*x1^5*x2^2+10*x1^5*x2^2*x4^2+11*x1^5*x2^2*x4^4+3*x1^5*x2^2*x3+6*x1^5*x2^2*x3*x4+10*x1^5*x2^2*x3*x4^2+4*x1^5*x2^2*x3*x4^3+6*x1^5*x2^2*x3*x4^5+3*x1^5*x2^2*x3^2+8*x1^5*x2^2*x3^2*x4+15*x1^5*x2^2*x3^2*x4^2+3*x1^5*x2^2*x3^2*x4^3+15*x1^5*x2^2*x3^2*x4^4+10*x1^5*x2^2*x3^2*x4^5+7*x1^5*x2^2*x3^3+16*x1^5*x2^2*x3^3*x4+8*x1^5*x2^2*x3^3*x4^4+10*x1^5*x2^2*x3^4+11*x1^5*x2^2*x3^4*x4+7*x1^5*x2^2*x3^4*x4^3+10*x1^5*x2^2*x3^4*x4^4+8*x1^5*x2^2*x3^5*x4+11*x1^5*x2^2*x3^5*x4^4+15*x1^5*x2^3+11*x1^5*x2^3*x4^2+10*x1^5*x2^3*x4^4+2*x1^5*x2^3*x3+9*x1^5*x2^3*x3*x4+6*x1^5*x2^3*x3*x4^2+10*x1^5*x2^3*x3*x4^3+7*x1^5*x2^3*x3*x4^5+15*x1^5*x2^3*x3^2*x4+9*x1^5*x2^3*x3^2*x4^2+2*x1^5*x2^3*x3^2*x4^3+9*x1^5*x2^3*x3^2*x4^4+6*x1^5*x2^3*x3^2*x4^5+11*x1^5*x2^3*x3^3+3*x1^5*x2^3*x3^3*x4+15*x1^5*x2^3*x3^3*x4^4+6*x1^5*x2^3*x3^4+10*x1^5*x2^3*x3^4*x4+11*x1^5*x2^3*x3^4*x4^3+6*x1^5*x2^3*x3^4*x4^4+15*x1^5*x2^3*x3^5*x4+10*x1^5*x2^3*x3^5*x4^4+14*x1^5*x2^4+3*x1^5*x2^4*x4^2+6*x1^5*x2^4*x4^4+10*x1^5*x2^4*x3+9*x1^5*x2^4*x3*x4+7*x1^5*x2^4*x3*x4^2+14*x1^5*x2^4*x3*x4^3+11*x1^5*x2^4*x3*x4^5+13*x1^5*x2^4*x3^2+9*x1^5*x2^4*x3^2*x4+2*x1^5*x2^4*x3^2*x4^2+10*x1^5*x2^4*x3^2*x4^3+2*x1^5*x2^4*x3^2*x4^4+7*x1^5*x2^4*x3^2*x4^5+10*x1^5*x2^4*x3^3+11*x1^5*x2^4*x3^3*x4+9*x1^5*x2^4*x3^3*x4^4+7*x1^5*x2^4*x3^4+6*x1^5*x2^4*x3^4*x4+10*x1^5*x2^4*x3^4*x4^3+7*x1^5*x2^4*x3^4*x4^4+9*x1^5*x2^4*x3^5*x4+6*x1^5*x2^4*x3^5*x4^4+12*x1^5*x2^5*x4^2+9*x1^5*x2^5*x4^4+5*x1^5*x2^5*x3+14*x1^5*x2^5*x3*x4+2*x1^5*x2^5*x3*x4^2+15*x1^5*x2^5*x3*x4^3+8*x1^5*x2^5*x3*x4^5+5*x1^5*x2^5*x3^2*x4+3*x1^5*x2^5*x3^2*x4^2+5*x1^5*x2^5*x3^2*x4^3+3*x1^5*x2^5*x3^2*x4^4+2*x1^5*x2^5*x3^2*x4^5+15*x1^5*x2^5*x3^3+14*x1^5*x2^5*x3^3*x4+5*x1^5*x2^5*x3^3*x4^4+2*x1^5*x2^5*x3^4+9*x1^5*x2^5*x3^4*x4+15*x1^5*x2^5*x3^4*x4^3+2*x1^5*x2^5*x3^4*x4^4+5*x1^5*x2^5*x3^5*x4+9*x1^5*x2^5*x3^5*x4^4+3*x1^5*x2^6+3*x1^5*x2^6*x4^2+14*x1^5*x2^6*x4^4+9*x1^5*x2^6*x3+5*x1^5*x2^6*x3*x4^2+15*x1^5*x2^6*x3*x4^3+3*x1^5*x2^6*x3*x4^5+11*x1^5*x2^6*x3^2+4*x1^5*x2^6*x3^2*x4+16*x1^5*x2^6*x3^2*x4^2+9*x1^5*x2^6*x3^2*x4^3+16*x1^5*x2^6*x3^2*x4^4+5*x1^5*x2^6*x3^2*x4^5+12*x1^5*x2^6*x3^3+6*x1^5*x2^6*x3^3*x4+4*x1^5*x2^6*x3^3*x4^4+5*x1^5*x2^6*x3^4+14*x1^5*x2^6*x3^4*x4+12*x1^5*x2^6*x3^4*x4^3+5*x1^5*x2^6*x3^4*x4^4+4*x1^5*x2^6*x3^5*x4+14*x1^5*x2^6*x3^5*x4^4+8*x1^5*x2^7+6*x1^5*x2^7*x4^2+8*x1^5*x2^7*x4^4+9*x1^5*x2^7*x3+2*x1^5*x2^7*x3*x4+15*x1^5*x2^7*x3*x4^2+7*x1^5*x2^7*x3*x4^3+9*x1^5*x2^7*x3*x4^5+13*x1^5*x2^7*x3^2+12*x1^5*x2^7*x3^2*x4+14*x1^5*x2^7*x3^2*x4^2+9*x1^5*x2^7*x3^2*x4^3+14*x1^5*x2^7*x3^2*x4^4+15*x1^5*x2^7*x3^2*x4^5+2*x1^5*x2^7*x3^3+8*x1^5*x2^7*x3^3*x4+12*x1^5*x2^7*x3^3*x4^4+15*x1^5*x2^7*x3^4+8*x1^5*x2^7*x3^4*x4+2*x1^5*x2^7*x3^4*x4^3+15*x1^5*x2^7*x3^4*x4^4+12*x1^5*x2^7*x3^5*x4+8*x1^5*x2^7*x3^5*x4^4+15*x1^5*x2^8+11*x1^5*x2^8*x4^2+4*x1^5*x2^8*x4^4+3*x1^5*x2^8*x3+3*x1^5*x2^8*x3*x4+16*x1^5*x2^8*x3*x4^2+6*x1^5*x2^8*x3*x4^3+13*x1^5*x2^8*x3*x4^5+7*x1^5*x2^8*x3^2+6*x1^5*x2^8*x3^2*x4+7*x1^5*x2^8*x3^2*x4^2+3*x1^5*x2^8*x3^2*x4^3+7*x1^5*x2^8*x3^2*x4^4+16*x1^5*x2^8*x3^2*x4^5+x1^5*x2^8*x3^3+11*x1^5*x2^8*x3^3*x4+6*x1^5*x2^8*x3^3*x4^4+16*x1^5*x2^8*x3^4+4*x1^5*x2^8*x3^4*x4+x1^5*x2^8*x3^4*x4^3+16*x1^5*x2^8*x3^4*x4^4+6*x1^5*x2^8*x3^5*x4+4*x1^5*x2^8*x3^5*x4^4+10*x1^5*x2^9+11*x1^5*x2^9*x4^2+12*x1^5*x2^9*x3+15*x1^5*x2^9*x3*x4+14*x1^5*x2^9*x3*x4^3+14*x1^5*x2^9*x3^2+12*x1^5*x2^9*x3^2*x4^3+10*x1^5*x2^9*x3^3*x4+6*x1^5*x2^10+9*x1^5*x2^10*x4^2+11*x1^5*x2^10*x4^4+11*x1^5*x2^10*x3+3*x1^5*x2^10*x3*x4+10*x1^5*x2^10*x3*x4^2+2*x1^5*x2^10*x3*x4^3+6*x1^5*x2^10*x3*x4^5+9*x1^5*x2^10*x3^2+8*x1^5*x2^10*x3^2*x4+15*x1^5*x2^10*x3^2*x4^2+11*x1^5*x2^10*x3^2*x4^3+15*x1^5*x2^10*x3^2*x4^4+10*x1^5*x2^10*x3^2*x4^5+7*x1^5*x2^10*x3^3+14*x1^5*x2^10*x3^3*x4+8*x1^5*x2^10*x3^3*x4^4+10*x1^5*x2^10*x3^4+11*x1^5*x2^10*x3^4*x4+7*x1^5*x2^10*x3^4*x4^3+10*x1^5*x2^10*x3^4*x4^4+8*x1^5*x2^10*x3^5*x4+11*x1^5*x2^10*x3^5*x4^4+6*x1^6+12*x1^6*x4^2+3*x1^6*x4^4+14*x1^6*x3+11*x1^6*x3*x4+8*x1^6*x3*x4^2+9*x1^6*x3*x4^3+15*x1^6*x3*x4^5+14*x1^6*x3^2+10*x1^6*x3^2*x4+9*x1^6*x3^2*x4^2+14*x1^6*x3^2*x4^3+6*x1^6*x3^2*x4^4+8*x1^6*x3^2*x4^5+15*x1^6*x3^3+9*x1^6*x3^3*x4+6*x1^6*x3^3*x4^3+10*x1^6*x3^3*x4^4+12*x1^6*x3^4+x1^6*x3^4*x4+15*x1^6*x3^4*x4^3+4*x1^6*x3^4*x4^4+3*x1^6*x3^5*x4+x1^6*x3^5*x4^4+8*x1^6*x2+3*x1^6*x2*x4^2+x1^6*x2*x4^4+x1^6*x2*x3+9*x1^6*x2*x3*x4+14*x1^6*x2*x3*x4^2+11*x1^6*x2*x3*x4^3+5*x1^6*x2*x3*x4^5+14*x1^6*x2*x3^2+9*x1^6*x2*x3^2*x4+3*x1^6*x2*x3^2*x4^2+x1^6*x2*x3^2*x4^3+2*x1^6*x2*x3^2*x4^4+14*x1^6*x2*x3^2*x4^5+5*x1^6*x2*x3^3+16*x1^6*x2*x3^3*x4+2*x1^6*x2*x3^3*x4^3+9*x1^6*x2*x3^3*x4^4+4*x1^6*x2*x3^4+6*x1^6*x2*x3^4*x4+5*x1^6*x2*x3^4*x4^3+7*x1^6*x2*x3^4*x4^4+x1^6*x2*x3^5*x4+6*x1^6*x2*x3^5*x4^4+8*x1^6*x2^2+4*x1^6*x2^2*x4^2+9*x1^6*x2^2*x4^4+13*x1^6*x2^2*x3+9*x1^6*x2^2*x3*x4+7*x1^6*x2^2*x3*x4^2+13*x1^6*x2^2*x3*x4^3+11*x1^6*x2^2*x3*x4^5+4*x1^6*x2^2*x3^2+13*x1^6*x2^2*x3^2*x4+10*x1^6*x2^2*x3^2*x4^2+13*x1^6*x2^2*x3^2*x4^3+x1^6*x2^2*x3^2*x4^4+7*x1^6*x2^2*x3^2*x4^5+11*x1^6*x2^2*x3^3+11*x1^6*x2^2*x3^3*x4+x1^6*x2^2*x3^3*x4^3+13*x1^6*x2^2*x3^3*x4^4+2*x1^6*x2^2*x3^4+3*x1^6*x2^2*x3^4*x4+11*x1^6*x2^2*x3^4*x4^3+12*x1^6*x2^2*x3^4*x4^4+9*x1^6*x2^2*x3^5*x4+3*x1^6*x2^2*x3^5*x4^4+10*x1^6*x2^3+16*x1^6*x2^3*x4^2+2*x1^6*x2^3*x4^4+14*x1^6*x2^3*x3+11*x1^6*x2^3*x3*x4^2+2*x1^6*x2^3*x3*x4^3+10*x1^6*x2^3*x3*x4^5+4*x1^6*x2^3*x3^2+x1^6*x2^3*x3^2*x4+6*x1^6*x2^3*x3^2*x4^2+14*x1^6*x2^3*x3^2*x4^3+4*x1^6*x2^3*x3^2*x4^4+11*x1^6*x2^3*x3^2*x4^5+10*x1^6*x2^3*x3^3+3*x1^6*x2^3*x3^3*x4+4*x1^6*x2^3*x3^3*x4^3+x1^6*x2^3*x3^3*x4^4+8*x1^6*x2^3*x3^4+12*x1^6*x2^3*x3^4*x4+10*x1^6*x2^3*x3^4*x4^3+14*x1^6*x2^3*x3^4*x4^4+2*x1^6*x2^3*x3^5*x4+12*x1^6*x2^3*x3^5*x4^4+11*x1^6*x2^4+5*x1^6*x2^4*x4^2+8*x1^6*x2^4*x4^4+14*x1^6*x2^4*x3+7*x1^6*x2^4*x3*x4+10*x1^6*x2^4*x3*x4^2+10*x1^6*x2^4*x3*x4^3+6*x1^6*x2^4*x3*x4^5+10*x1^6*x2^4*x3^2+4*x1^6*x2^4*x3^2*x4+7*x1^6*x2^4*x3^2*x4^2+14*x1^6*x2^4*x3^2*x4^3+16*x1^6*x2^4*x3^2*x4^4+10*x1^6*x2^4*x3^2*x4^5+6*x1^6*x2^4*x3^3+11*x1^6*x2^4*x3^3*x4+16*x1^6*x2^4*x3^3*x4^3+4*x1^6*x2^4*x3^3*x4^4+15*x1^6*x2^4*x3^4+14*x1^6*x2^4*x3^4*x4+6*x1^6*x2^4*x3^4*x4^3+5*x1^6*x2^4*x3^4*x4^4+8*x1^6*x2^4*x3^5*x4+14*x1^6*x2^4*x3^5*x4^4+5*x1^6*x2^5+3*x1^6*x2^5*x4^2+12*x1^6*x2^5*x4^4+15*x1^6*x2^5*x3+4*x1^6*x2^5*x3*x4+15*x1^6*x2^5*x3*x4^2+8*x1^6*x2^5*x3*x4^3+9*x1^6*x2^5*x3*x4^5+8*x1^6*x2^5*x3^2+6*x1^6*x2^5*x3^2*x4+2*x1^6*x2^5*x3^2*x4^2+15*x1^6*x2^5*x3^2*x4^3+7*x1^6*x2^5*x3^2*x4^4+15*x1^6*x2^5*x3^2*x4^5+9*x1^6*x2^5*x3^3+15*x1^6*x2^5*x3^3*x4+7*x1^6*x2^5*x3^3*x4^3+6*x1^6*x2^5*x3^3*x4^4+14*x1^6*x2^5*x3^4+4*x1^6*x2^5*x3^4*x4+9*x1^6*x2^5*x3^4*x4^3+16*x1^6*x2^5*x3^4*x4^4+12*x1^6*x2^5*x3^5*x4+4*x1^6*x2^5*x3^5*x4^4+2*x1^6*x2^6+9*x1^6*x2^6*x4^2+13*x1^6*x2^6*x4^4+2*x1^6*x2^6*x3+7*x1^6*x2^6*x3*x4+12*x1^6*x2^6*x3*x4^2+14*x1^6*x2^6*x3*x4^3+14*x1^6*x2^6*x3*x4^5+15*x1^6*x2^6*x3^2+15*x1^6*x2^6*x3^2*x4+5*x1^6*x2^6*x3^2*x4^2+2*x1^6*x2^6*x3^2*x4^3+9*x1^6*x2^6*x3^2*x4^4+12*x1^6*x2^6*x3^2*x4^5+14*x1^6*x2^6*x3^3+10*x1^6*x2^6*x3^3*x4+9*x1^6*x2^6*x3^3*x4^3+15*x1^6*x2^6*x3^3*x4^4+x1^6*x2^6*x3^4+10*x1^6*x2^6*x3^4*x4+14*x1^6*x2^6*x3^4*x4^3+6*x1^6*x2^6*x3^4*x4^4+13*x1^6*x2^6*x3^5*x4+10*x1^6*x2^6*x3^5*x4^4+15*x1^6*x2^7+10*x1^6*x2^7*x4^2+5*x1^6*x2^7*x4^4+5*x1^6*x2^7*x3+9*x1^6*x2^7*x3*x4+2*x1^6*x2^7*x3*x4^2+4*x1^6*x2^7*x3*x4^3+8*x1^6*x2^7*x3*x4^5+12*x1^6*x2^7*x3^2+11*x1^6*x2^7*x3^2*x4+15*x1^6*x2^7*x3^2*x4^2+5*x1^6*x2^7*x3^2*x4^3+10*x1^6*x2^7*x3^2*x4^4+2*x1^6*x2^7*x3^2*x4^5+8*x1^6*x2^7*x3^3+5*x1^6*x2^7*x3^3*x4+10*x1^6*x2^7*x3^3*x4^3+11*x1^6*x2^7*x3^3*x4^4+3*x1^6*x2^7*x3^4+13*x1^6*x2^7*x3^4*x4+8*x1^6*x2^7*x3^4*x4^3+x1^6*x2^7*x3^4*x4^4+5*x1^6*x2^7*x3^5*x4+13*x1^6*x2^7*x3^5*x4^4+9*x1^6*x2^8+6*x1^6*x2^8*x4^2+11*x1^6*x2^8*x4^4+7*x1^6*x2^8*x3+x1^6*x2^8*x3*x4+x1^6*x2^8*x3*x4^2+3*x1^6*x2^8*x3*x4^3+4*x1^6*x2^8*x3*x4^5+8*x1^6*x2^8*x3^2+14*x1^6*x2^8*x3^2*x4+16*x1^6*x2^8*x3^2*x4^2+7*x1^6*x2^8*x3^2*x4^3+5*x1^6*x2^8*x3^2*x4^4+x1^6*x2^8*x3^2*x4^5+4*x1^6*x2^8*x3^3+3*x1^6*x2^8*x3^3*x4+5*x1^6*x2^8*x3^3*x4^3+14*x1^6*x2^8*x3^3*x4^4+10*x1^6*x2^8*x3^4+15*x1^6*x2^8*x3^4*x4+4*x1^6*x2^8*x3^4*x4^3+9*x1^6*x2^8*x3^4*x4^4+11*x1^6*x2^8*x3^5*x4+15*x1^6*x2^8*x3^5*x4^4+2*x1^6*x2^9+8*x1^6*x2^9*x4^2+6*x1^6*x2^9*x3+13*x1^6*x2^9*x3*x4+7*x1^6*x2^9*x3*x4^3+2*x1^6*x2^9*x3^2+6*x1^6*x2^9*x3^2*x4^3+3*x1^6*x2^9*x3^3*x4+10*x1^6*x2^10+12*x1^6*x2^10*x4^2+9*x1^6*x2^10*x4^4+9*x1^6*x2^10*x3+3*x1^6*x2^10*x3*x4+7*x1^6*x2^10*x3*x4^2+14*x1^6*x2^10*x3*x4^3+11*x1^6*x2^10*x3*x4^5+9*x1^6*x2^10*x3^2+13*x1^6*x2^10*x3^2*x4+10*x1^6*x2^10*x3^2*x4^2+9*x1^6*x2^10*x3^2*x4^3+x1^6*x2^10*x3^2*x4^4+7*x1^6*x2^10*x3^2*x4^5+11*x1^6*x2^10*x3^3+7*x1^6*x2^10*x3^3*x4+x1^6*x2^10*x3^3*x4^3+13*x1^6*x2^10*x3^3*x4^4+2*x1^6*x2^10*x3^4+3*x1^6*x2^10*x3^4*x4+11*x1^6*x2^10*x3^4*x4^3+12*x1^6*x2^10*x3^4*x4^4+9*x1^6*x2^10*x3^5*x4+3*x1^6*x2^10*x3^5*x4^4+5*x1^7+13*x1^7*x4^2+12*x1^7*x3+13*x1^7*x3*x4+14*x1^7*x3*x4^3+12*x1^7*x3^2+12*x1^7*x3^2*x4^3+3*x1^7*x3^3*x4+7*x1^7*x2+10*x1^7*x2*x4^2+4*x1^7*x2*x3+10*x1^7*x2*x3*x4+16*x1^7*x2*x3*x4^3+4*x1^7*x2*x3^2+4*x1^7*x2*x3^2*x4^3+x1^7*x2*x3^3*x4+4*x1^7*x2^2+5*x1^7*x2^2*x4^2+2*x1^7*x2^2*x3+5*x1^7*x2^2*x3*x4+8*x1^7*x2^2*x3*x4^3+2*x1^7*x2^2*x3^2+2*x1^7*x2^2*x3^2*x4^3+9*x1^7*x2^2*x3^3*x4+3*x1^7*x2^3+3*x1^7*x2^3*x4^2+8*x1^7*x2^3*x3+3*x1^7*x2^3*x3*x4+15*x1^7*x2^3*x3*x4^3+8*x1^7*x2^3*x3^2+8*x1^7*x2^3*x3^2*x4^3+2*x1^7*x2^3*x3^3*x4+6*x1^7*x2^4+12*x1^7*x2^4*x4^2+15*x1^7*x2^4*x3+12*x1^7*x2^4*x3*x4+9*x1^7*x2^4*x3*x4^3+15*x1^7*x2^4*x3^2+15*x1^7*x2^4*x3^2*x4^3+8*x1^7*x2^4*x3^3*x4+10*x1^7*x2^5+x1^7*x2^5*x4^2+14*x1^7*x2^5*x3+x1^7*x2^5*x3*x4+5*x1^7*x2^5*x3*x4^3+14*x1^7*x2^5*x3^2+14*x1^7*x2^5*x3^2*x4^3+12*x1^7*x2^5*x3^3*x4+5*x1^7*x2^6+11*x1^7*x2^6*x4^2+x1^7*x2^6*x3+11*x1^7*x2^6*x3*x4+4*x1^7*x2^6*x3*x4^3+x1^7*x2^6*x3^2+x1^7*x2^6*x3^2*x4^3+13*x1^7*x2^6*x3^3*x4+16*x1^7*x2^7*x4^2+3*x1^7*x2^7*x3+16*x1^7*x2^7*x3*x4+12*x1^7*x2^7*x3*x4^3+3*x1^7*x2^7*x3^2+3*x1^7*x2^7*x3^2*x4^3+5*x1^7*x2^7*x3^3*x4+11*x1^7*x2^8+8*x1^7*x2^8*x4^2+10*x1^7*x2^8*x3+8*x1^7*x2^8*x3*x4+6*x1^7*x2^8*x3*x4^3+10*x1^7*x2^8*x3^2+10*x1^7*x2^8*x3^2*x4^3+11*x1^7*x2^8*x3^3*x4+16*x1^7*x2^9+5*x1^7*x2^10+5*x1^7*x2^10*x4^2+2*x1^7*x2^10*x3+5*x1^7*x2^10*x3*x4+8*x1^7*x2^10*x3*x4^3+2*x1^7*x2^10*x3^2+2*x1^7*x2^10*x3^2*x4^3+9*x1^7*x2^10*x3^3*x4+x1^8+x1^8*x4^2+9*x1^8*x3+15*x1^8*x3*x4+2*x1^8*x3*x4^3+9*x1^8*x3^2*x4^3+10*x1^8*x3^3*x4+10*x1^8*x2+6*x1^8*x2*x4^2+3*x1^8*x2*x3+5*x1^8*x2*x3*x4+12*x1^8*x2*x3*x4^3+3*x1^8*x2*x3^2*x4^3+9*x1^8*x2*x3^3*x4+4*x1^8*x2^2+3*x1^8*x2^2*x4^2+10*x1^8*x2^2*x3+11*x1^8*x2^2*x3*x4+6*x1^8*x2^2*x3*x4^3+10*x1^8*x2^2*x3^2*x4^3+13*x1^8*x2^2*x3^3*x4+10*x1^8*x2^3+12*x1^8*x2^3*x4^2+6*x1^8*x2^3*x3+10*x1^8*x2^3*x3*x4+7*x1^8*x2^3*x3*x4^3+6*x1^8*x2^3*x3^2*x4^3+x1^8*x2^3*x3^3*x4+4*x1^8*x2^4+14*x1^8*x2^4*x4^2+7*x1^8*x2^4*x3+6*x1^8*x2^4*x3*x4+11*x1^8*x2^4*x3*x4^3+7*x1^8*x2^4*x3^2*x4^3+4*x1^8*x2^4*x3^3*x4+6*x1^8*x2^5+4*x1^8*x2^5*x4^2+2*x1^8*x2^5*x3+9*x1^8*x2^5*x3*x4+8*x1^8*x2^5*x3*x4^3+2*x1^8*x2^5*x3^2*x4^3+6*x1^8*x2^5*x3^3*x4+14*x1^8*x2^6+10*x1^8*x2^6*x4^2+5*x1^8*x2^6*x3+14*x1^8*x2^6*x3*x4+3*x1^8*x2^6*x3*x4^3+5*x1^8*x2^6*x3^2*x4^3+15*x1^8*x2^6*x3^3*x4+11*x1^8*x2^7+13*x1^8*x2^7*x4^2+15*x1^8*x2^7*x3+8*x1^8*x2^7*x3*x4+9*x1^8*x2^7*x3*x4^3+15*x1^8*x2^7*x3^2*x4^3+11*x1^8*x2^7*x3^3*x4+5*x1^8*x2^8+15*x1^8*x2^8*x4^2+16*x1^8*x2^8*x3+4*x1^8*x2^8*x3*x4+13*x1^8*x2^8*x3*x4^3+16*x1^8*x2^8*x3^2*x4^3+14*x1^8*x2^8*x3^3*x4+9*x1^8*x2^10+3*x1^8*x2^10*x4^2+10*x1^8*x2^10*x3+11*x1^8*x2^10*x3*x4+6*x1^8*x2^10*x3*x4^3+10*x1^8*x2^10*x3^2*x4^3+13*x1^8*x2^10*x3^3*x4+11*x1^9+6*x1^9*x4^2+8*x1^9*x3+16*x1^9*x3*x4+15*x1^9*x3*x4^3+7*x1^9*x3^2+8*x1^9*x3^2*x4^3+5*x1^9*x3^3*x4+5*x1^9*x2+2*x1^9*x2*x4^2+14*x1^9*x2*x3+11*x1^9*x2*x3*x4+5*x1^9*x2*x3*x4^3+8*x1^9*x2*x3^2+14*x1^9*x2*x3^2*x4^3+13*x1^9*x2*x3^3*x4+14*x1^9*x2^2+x1^9*x2^2*x4^2+7*x1^9*x2^2*x3+14*x1^9*x2^2*x3*x4+11*x1^9*x2^2*x3*x4^3+4*x1^9*x2^2*x3^2+7*x1^9*x2^2*x3^2*x4^3+15*x1^9*x2^2*x3^3*x4+14*x1^9*x2^3+4*x1^9*x2^3*x4^2+11*x1^9*x2^3*x3+5*x1^9*x2^3*x3*x4+10*x1^9*x2^3*x3*x4^3+16*x1^9*x2^3*x3^2+11*x1^9*x2^3*x3^2*x4^3+9*x1^9*x2^3*x3^3*x4+4*x1^9*x2^4+16*x1^9*x2^4*x4^2+10*x1^9*x2^4*x3+3*x1^9*x2^4*x3*x4+6*x1^9*x2^4*x3*x4^3+13*x1^9*x2^4*x3^2+10*x1^9*x2^4*x3^2*x4^3+2*x1^9*x2^4*x3^3*x4+10*x1^9*x2^5+7*x1^9*x2^5*x4^2+15*x1^9*x2^5*x3+13*x1^9*x2^5*x3*x4+9*x1^9*x2^5*x3*x4^3+11*x1^9*x2^5*x3^2+15*x1^9*x2^5*x3^2*x4^3+3*x1^9*x2^5*x3^3*x4+5*x1^9*x2^6+9*x1^9*x2^6*x4^2+12*x1^9*x2^6*x3+7*x1^9*x2^6*x3*x4+14*x1^9*x2^6*x3*x4^3+2*x1^9*x2^6*x3^2+12*x1^9*x2^6*x3^2*x4^3+16*x1^9*x2^6*x3^3*x4+16*x1^9*x2^7+10*x1^9*x2^7*x4^2+2*x1^9*x2^7*x3+4*x1^9*x2^7*x3*x4+8*x1^9*x2^7*x3*x4^3+6*x1^9*x2^7*x3^2+2*x1^9*x2^7*x3^2*x4^3+14*x1^9*x2^7*x3^3*x4+3*x1^9*x2^8+5*x1^9*x2^8*x4^2+x1^9*x2^8*x3+2*x1^9*x2^8*x3*x4+4*x1^9*x2^8*x3*x4^3+3*x1^9*x2^8*x3^2+x1^9*x2^8*x3^2*x4^3+7*x1^9*x2^8*x3^3*x4+6*x1^9*x2^9+x1^9*x2^10*x4^2+7*x1^9*x2^10*x3+14*x1^9*x2^10*x3*x4+11*x1^9*x2^10*x3*x4^3+4*x1^9*x2^10*x3^2+7*x1^9*x2^10*x3^2*x4^3+15*x1^9*x2^10*x3^3*x4+10*x1^10+9*x1^10*x2+13*x1^10*x2^2+x1^10*x2^3+4*x1^10*x2^4+6*x1^10*x2^5+15*x1^10*x2^6+11*x1^10*x2^7+14*x1^10*x2^8+13*x1^10*x2^10+x1^11+6*x1^11*x2+3*x1^11*x2^2+12*x1^11*x2^3+14*x1^11*x2^4+4*x1^11*x2^5+10*x1^11*x2^6+13*x1^11*x2^7+15*x1^11*x2^8+3*x1^11*x2^10+6*x1^12+2*x1^12*x2+x1^12*x2^2+4*x1^12*x2^3+16*x1^12*x2^4+7*x1^12*x2^5+9*x1^12*x2^6+10*x1^12*x2^7+5*x1^12*x2^8+x1^12*x2^10",
                Zp64(17));
        Assert.assertEquals(4, Factor(poly).size());
    }

    @Test(timeout = 200_000)
    public void testMultivariateFactorization46() throws Exception {
        MultivariatePolynomial<BigInteger> poly = PolynomialMethods.polyPow(MultivariatePolynomial.parse("z+3*y^2+x^3*y"), 50).decrement();
        Assert.assertEquals(6, PolynomialMethods.Factor(poly).size());
    }

    @Test
    public void testMultivariateFactorization47() throws Exception {
        MultivariatePolynomialZp64 poly = PolynomialMethods.polyPow(
                MultivariatePolynomialZp64.parse("x^2*y - z^2*x^3 - y^2*z - 1", Rings.Zp64(2)), 20).decrement();

        FactorDecompositionTest.assertFactorization(poly, Factor(poly));
    }

    @Test(timeout = 200_000L)
    public void testMultivariateFactorization48() throws Exception {
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(5, Z);

        MultivariatePolynomial<BigInteger>
                p1 = ring.parse("52037*x1^4*x2^3*x3^2*x4*x5+18900*x1^3*x2^10*x3*x4+42072*x1^4*x2*x4^7*x5^4+25153*x1^3*x2^10*x4^5*x5+40049*x2*x3^9*x4^3*x5^10+21479*x1^4*x3^5*x4^14*x5^2+46360*x1^2*x2^5*x3^8*x4^9*x5^2+1426*x1^5*x2^2*x3^2*x4^7*x5^12+53121*x1^3*x2^6*x3^11*x5^9+55693*x1^13*x2^2*x3^10*x5^4+47321*x1^10*x2*x3^7*x4^5*x5^7+53426*x2^11*x3^13*x4^4*x5^2+41556*x2^11*x3*x4^7*x5^12+51863*x1^2*x2^8*x3^10*x4^12*x5+15656*x1^5*x2^8*x3^6*x4^7*x5^9+12962*x1^8*x2^6*x3^5*x4^7*x5^9+51360*x1^5*x2^13*x3^8*x4^5*x5^8+63981*x1^6*x2^14*x3^5*x4^13*x5^3+36882*x1^13*x2^6*x3^8*x4^6*x5^14+3510*x1^5*x2^5*x3^12*x4^12*x5^13+21301*x1^10*x2^14*x3^11*x4^9*x5^9"),
                p2 = ring.parse("26692*x1^2*x2*x3^2*x5^13+62118*x1^2*x2^6*x3^5*x4^2*x5^7+828*x1^8*x2^6*x3^3*x4^4*x5+63108*x1^14*x2^2*x3^3*x4*x5^4+35392*x1^2*x2^2*x3^14*x4^6+18686*x1^6*x2^9*x3^2*x4^3*x5^5+55647*x2^4*x3^13*x4^6*x5^2+52635*x1^14*x2^5*x4^5*x5+191*x1^2*x3^8*x4^10*x5^6+65201*x1^12*x2*x3^8*x4*x5^5+3551*x1^9*x2*x3^12*x4^4*x5+23009*x1^6*x2^9*x3^6*x4^6+21159*x1^2*x2^9*x3^3*x4^13*x5^2+38666*x1^9*x2^9*x3^3*x5^10+25692*x1^12*x2^9*x3^5*x4^3*x5^5+21256*x1*x2^4*x3^12*x4^12*x5^6+32978*x1^6*x3^14*x4^6*x5^11+35661*x1^12*x2^3*x3^3*x4^14*x5^9+10212*x1^13*x2^14*x3^6*x4^8*x5+65089*x1^9*x2^13*x3^10*x4*x5^10+31555*x1^9*x2^8*x3^6*x4^9*x5^12"),
                p3 = ring.parse("51940*x3^12*x4*x5^9+28440*x1^7*x2^2*x3^6*x4^4*x5^6+64657*x1^4*x2^9*x3*x4^12+49150*x1*x2^7*x3^6*x4^3*x5^10+39701*x1^7*x2^4*x3^8*x5^8+55736*x1^11*x2^4*x3^7*x5^6+7442*x1^12*x2*x3^11*x4^5+20639*x1^7*x2^3*x4^13*x5^8+45382*x1^4*x2^2*x3^14*x4^4*x5^8+34045*x1^6*x2*x3^10*x4^9*x5^8+57079*x1^11*x2^9*x3^4*x4^3*x5^8+15678*x1^14*x2^5*x3^11*x4^2*x5^3+35697*x1^6*x2^11*x3^4*x4^6*x5^9+20453*x1^11*x2^12*x3^6*x4^7*x5+20337*x1^14*x2^9*x3^5*x4^7*x5^6+32687*x1^13*x2^8*x3^12*x4^3*x5^5+18470*x1^10*x2^6*x3^12*x4^10*x5^3+26468*x1^13*x2^3*x4^13*x5^13+53498*x1^13*x2^4*x3^13*x4^9*x5^6+51463*x1^9*x2^10*x3^10*x4^13*x5^5+17452*x1^11*x2^12*x3^12*x4^14*x5^11");
        List<MultivariatePolynomial<BigInteger>> o = Arrays.asList(p1, p2, p3);
        List<MultivariatePolynomial<BigInteger>> m = Stream.of(p1, p2, p3).map(p -> p.setRing(Zp(17)).setRing(Z)).collect(Collectors.toList());

        MultivariatePolynomial<BigInteger> poly = p1.createOne().multiply(m);
        Assert.assertEquals(3, Factor(poly).size());
    }


    @Ignore
    @Test
    public void testMultivariateFactorization49() throws Exception {
        RandomGenerator rnd = getRandom();
        rnd.setSeed(1);

        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = MultivariateRing(7, Z);
        Supplier<MultivariatePolynomial<BigInteger>> random
                = () -> randomPolynomial(ring.nVariables(), 0, 15,
                20, Z, MonomialOrder.GREVLEX, r -> new BigInteger(16, r), rnd);

        MultivariatePolynomial<BigInteger> poly = ring.multiply(random.get(), random.get(), random.get());
        for (int i = 0; i < 1000; ++i) {
            long start = System.nanoTime();
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> factor = Factor(poly);
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            Assert.assertEquals(poly, factor.multiply());
        }

        // 101s  <- new
        // 98s  <- new

        // 116s  <- old
        // 103s  <- old
    }

    static <E> MultivariatePolynomial<E> randomPolynomial(int nVars, int minDegree, int maxDegree, int size, Ring<E> ring,
                                                          Comparator<DegreeVector> ordering,
                                                          Function<RandomGenerator, E> method, RandomGenerator rnd) {
        @SuppressWarnings("unchecked")
        Monomial<E>[] terms = new Monomial[size];
        for (int i = 0; i < size; i++) {
            E cfx = method.apply(rnd);
            int[] exponents = RandomUtil.randomIntArray(nVars, minDegree, maxDegree, rnd);
            terms[i] = new Monomial<>(exponents, cfx);
        }
        return MultivariatePolynomial.create(nVars, ring, ordering, terms);
    }

    /* ==================================== Test data =============================================== */

    static double DETALIZATION_PERCENT = 100;
    static boolean PRINT_FACTORS = false;
    static boolean DO_ASSERTION = true;

    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void testFactorizationAlgorithm(FactorizationInput.SampleDecompositionSource<Poly> source,
                                    int nIterations,
                                    FactorizationInput.FactorizationAlgorithm<Poly> algorithm) {
        System.out.println("Testing factorization algorithm " + algorithm.name);
        System.out.println("Input source: " + source);

        DescriptiveStatistics
                lTiming = new DescriptiveStatistics(),
                bTiming = new DescriptiveStatistics();

        int prevProgress = -1, currProgress;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.asList(lTiming, bTiming).forEach(DescriptiveStatistics::clear);

            if ((currProgress = (int) (DETALIZATION_PERCENT * n / nIterations)) != prevProgress) {
                prevProgress = currProgress;
                System.out.print(">");
                System.out.flush();
            }
            FactorizationInput.SampleDecomposition<Poly> sample = source.next();
            PolynomialFactorDecomposition<Poly> decomposition = null;
            try {
                if (PRINT_FACTORS) {
                    System.out.println("\n");
                    System.out.println("#N = " + n + "   " + new SimpleDateFormat("HH:mm:ss").format(new Date()));
                    for (Poly factor : sample.factors)
                        System.out.println(String.format("MultivariatePolynomial.parse(\"%s\"),", factor));
                    System.out.println();
                }
                long start = System.nanoTime();
                decomposition = algorithm.algorithm.apply(sample.poly);
                lTiming.addValue(System.nanoTime() - start);
                if (DO_ASSERTION)
                    sample.assertFactorization(decomposition);
            } catch (Throwable throwable) {
                System.out.println("============ Error ============");
                System.out.println("Ring: " + sample.poly.coefficientRingToString());
                System.out.println("Polynomial: " + sample.poly);
                System.out.println("Expected factorization: " + Arrays.toString(sample.factors));
                System.out.println("Actual decomposition: " + decomposition);
                throw throwable;
            }
        }

        System.out.println(source.statisticsToString());

        System.out.println("\n============ Timings ============");
        System.out.println("Statistics : " + TimeUnits.statisticsNanotime(bTiming));
    }

    public static void
    testFactorizationAlgorithm(FactorizationInput.SampleDecompositionSource<MultivariatePolynomialZp64> source,
                               int nIterations,
                               FactorizationInput.FactorizationAlgorithm<MultivariatePolynomialZp64> lAlgorithm,
                               FactorizationInput.FactorizationAlgorithm<MultivariatePolynomial<BigInteger>> bAlgorithm) {
        System.out.println("Testing factorization algorithm " + lAlgorithm.name);
        System.out.println("Input source: " + source);

        DescriptiveStatistics
                lTiming = new DescriptiveStatistics(),
                bTiming = new DescriptiveStatistics();

        int prevProgress = -1, currProgress;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10)
                Arrays.asList(lTiming, bTiming).forEach(DescriptiveStatistics::clear);

            if ((currProgress = (int) (DETALIZATION_PERCENT * n / nIterations)) != prevProgress) {
                prevProgress = currProgress;
                System.out.print(">");
                System.out.flush();
            }
            FactorizationInput.SampleDecomposition<MultivariatePolynomialZp64> sample = source.next();
//            System.out.println(sample.poly.ring);
//            System.out.println(PolynomialFactorDecomposition.create(Arrays.asList(sample.factors)));
            PolynomialFactorDecomposition<MultivariatePolynomialZp64> lDecomposition = null;
            PolynomialFactorDecomposition<MultivariatePolynomial<BigInteger>> bDecomposition = null;
            try {
                if (PRINT_FACTORS) {
                    System.out.println("\n");
                    System.out.println("#N = " + n + "   " + new SimpleDateFormat("HH:mm:ss").format(new Date()));
                    System.out.println(sample.poly.ring);
                    for (MultivariatePolynomialZp64 factor : sample.factors)
                        System.out.println(String.format("MultivariatePolynomialZp64.parse(\"%s\", ring, vars),", factor));
                    System.out.println();
                }

                long start;
                if (lAlgorithm != null) {
                    start = System.nanoTime();
                    lDecomposition = lAlgorithm.algorithm.apply(sample.poly);
                    lTiming.addValue(System.nanoTime() - start);
                    if (DO_ASSERTION)
                        sample.assertFactorization(lDecomposition);
                }

                if (bAlgorithm != null) {
                    FactorizationInput.SampleDecomposition<MultivariatePolynomial<BigInteger>>
                            bSample = toBigPoly(sample);
                    start = System.nanoTime();
                    bDecomposition = bAlgorithm.algorithm.apply(bSample.poly);
                    bTiming.addValue(System.nanoTime() - start);
                    if (DO_ASSERTION)
                        bSample.assertFactorization(bDecomposition);
                }
            } catch (Throwable throwable) {
                System.out.println("============ Error ============");
                System.out.println("Ring: " + sample.poly.ring);
                System.out.println("Polynomial: " + sample.poly);
                System.out.println("Expected factorization: " + Arrays.toString(sample.factors));
                System.out.println("Actual decomposition (longs): " + lDecomposition);
                System.out.println("Actual decomposition (BigInts): " + bDecomposition);
                throw throwable;
            }
        }

        System.out.println(source.statisticsToString());

        System.out.println("\n============ Timings ============");
        System.out.println("Machine integers: " + TimeUnits.statisticsNanotime(lTiming));
        System.out.println("Big integers    : " + TimeUnits.statisticsNanotime(bTiming));
    }

    @SuppressWarnings("unchecked")
    static FactorizationInput.SampleDecomposition<MultivariatePolynomial<BigInteger>>
    toBigPoly(FactorizationInput.SampleDecomposition<MultivariatePolynomialZp64> decomposition) {
        return new FactorizationInput.SampleDecomposition<>(
                Arrays.stream(decomposition.factors)
                        .map(MultivariatePolynomialZp64::toBigPoly)
                        .toArray(MultivariatePolynomial[]::new));
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorizationInput.SampleDecompositionSource<Poly> orderVarsByDegree(FactorizationInput.SampleDecompositionSource<Poly> source) {
        return new FactorizationInput.SampleDecompositionSource<Poly>() {
            @Override
            public FactorizationInput.SampleDecomposition<Poly> next0() {
                FactorizationInput.SampleDecomposition<Poly> sample = source.next0();
                int[] degrees = sample.poly.degrees();
                int[] variables = ArraysUtil.sequence(0, degrees.length);
                ArraysUtil.quickSort(ArraysUtil.negate(degrees), variables);
                return new FactorizationInput.SampleDecomposition<>(
                        Arrays.stream(sample.factors)
                                .map(p -> AMultivariatePolynomial.renameVariables(p, variables))
                                .toArray(sample.poly::createArray));
            }

            @Override
            public String toString() {
                return super.toString() + " (vars ordered)";
            }
        };
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorizationInput.SampleDecompositionSource<Poly> filterMonomialContent(FactorizationInput.SampleDecompositionSource<Poly> source) {
        return new FactorizationInput.SampleDecompositionSource<Poly>() {
            @Override
            public FactorizationInput.SampleDecomposition<Poly> next0() {
                FactorizationInput.SampleDecomposition<Poly> sample = source.next0();
                for (Poly factor : sample.factors) {
                    if (!factor.monomialContent().isZeroVector())
                        factor.increment();
                    assert factor.monomialContent().isZeroVector();
                }
                return new FactorizationInput.SampleDecomposition<>(sample.factors);
            }

            @Override
            public String toString() {
                return super.toString() + " (content filtered)";
            }
        };
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorizationInput.SampleDecompositionSource<Poly> filterNonPrimitive(FactorizationInput.SampleDecompositionSource<Poly> source) {
        return new FactorizationInput.SampleDecompositionSource<Poly>() {
            @Override
            public FactorizationInput.SampleDecomposition<Poly> next0() {
                FactorizationInput.SampleDecomposition<Poly> sample = source.next0();
                if (MultivariateFactorization.factorToPrimitive(sample.poly).isTrivial()) {
                    for (Poly factor : sample.factors)
                        factor.primitivePart();
                    return new FactorizationInput.SampleDecomposition<>(sample.factors);
                } else return next0();
            }

            @Override
            public String toString() {
                return super.toString() + " (square-free)";
            }
        };
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    FactorizationInput.SampleDecompositionSource<Poly> filterNonSquareFree(FactorizationInput.SampleDecompositionSource<Poly> source) {
        return new FactorizationInput.SampleDecompositionSource<Poly>() {
            @Override
            public FactorizationInput.SampleDecomposition<Poly> next0() {
                FactorizationInput.SampleDecomposition<Poly> sample = source.next0();
                if (MultivariateSquareFreeFactorization.isSquareFree(sample.poly))
                    return sample;
                else return next0();
            }

            @Override
            public String toString() {
                return super.toString() + " (square-free)";
            }
        };
    }

    public static final class lSampleDecompositionSource
            extends FactorizationInput.SampleDecompositionSource<MultivariatePolynomialZp64> {
        final int
                nFactorsMin, nFactorsMax,
                nVarsMin, nVarsMax,
                minSize, maxSize,
                minDegree, maxDegree;

        int minModulusBits = 10, maxModulusBits = 32;

        public lSampleDecompositionSource(int nFactorsMin, int nFactorsMax,
                                          int nVarsMin, int nVarsMax,
                                          int minSize, int maxSize,
                                          int minDegree, int maxDegree) {
            this.nFactorsMin = nFactorsMin;
            this.nFactorsMax = nFactorsMax;
            this.nVarsMin = nVarsMin;
            this.nVarsMax = nVarsMax;
            this.minSize = minSize;
            this.maxSize = maxSize;
            this.minDegree = minDegree;
            this.maxDegree = maxDegree;
        }

        final RandomGenerator rnd = getRandom();
        final RandomDataGenerator rndd = getRandomData();

        @Override
        public FactorizationInput.SampleDecomposition<MultivariatePolynomialZp64> next0() {
            IntegersZp64 domain = new IntegersZp64(
                    getModulusRandom(rndd.nextInt(minModulusBits, maxModulusBits)));

            int nVariables = rndd.nextInt(nVarsMin, nVarsMax);
            int nFactors = rndd.nextInt(nFactorsMin, nFactorsMax);

            MultivariatePolynomialZp64[] factors = new MultivariatePolynomialZp64[nFactors];
            while (nFactors > 0) {
                MultivariatePolynomialZp64 sample =
                        RandomMultivariatePolynomials.randomPolynomial(nVariables,
                                rndd.nextInt(minDegree, maxDegree),
                                rndd.nextInt(minSize, maxSize),
                                domain, MonomialOrder.DEFAULT, rnd);

                if (sample.isConstant() || sample.isMonomial())
                    continue;
                factors[--nFactors] = sample;
            }
            return new FactorizationInput.SampleDecomposition<>(factors);
        }

        private static String range(int from, int to) {
            return "[" + from + ", " + to + "]";
        }

        @Override
        public String toString() {
            return "Sample data: " +
                    " #factors ∈ " + range(nFactorsMin, nFactorsMax) +
                    " #variables ∈ " + range(nVarsMin, nVarsMax) +
                    ", deg ∈ " + range(minDegree, maxDegree) +
                    ", size ∈ " + range(minSize, maxSize);
        }
    }

    public static final class SampleDecompositionSourceZ
            extends FactorizationInput.SampleDecompositionSource<MultivariatePolynomial<BigInteger>> {
        final lSampleDecompositionSource lSource;

        public SampleDecompositionSourceZ(lSampleDecompositionSource lSource) {
            this.lSource = lSource;
        }

        @Override
        public FactorizationInput.SampleDecomposition<MultivariatePolynomial<BigInteger>> next0() {
            return new FactorizationInput.SampleDecomposition<>(Arrays.stream(lSource.next().factors)
                    .map(MultivariatePolynomialZp64::asPolyZSymmetric).toArray(MultivariatePolynomial[]::new));
        }
    }
}