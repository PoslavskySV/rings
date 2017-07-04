package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.univar.IrreduciblePolynomials;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MultivariateFactorizationTest extends AbstractPolynomialTest {
    @Ignore
    @Test
    public void test2() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(67);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("b*a^5 + a^4*b^2 + 11 + b^3", domain),
                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b + 17*b^2 + 1", domain),
                c = lMultivariatePolynomialZp.parse("b^3*a^4 + a^4 + b", domain),
                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b, c, d);

        System.out.println(base);
        for (int i = 0; i < 1000; i++) {
            long start = System.nanoTime();
            Assert.assertEquals(4, MultivariateFactorization.bivariateDenseFactorSquareFree(base).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
        }
    }

    @Test
    public void testBivariateFactorizationRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        DescriptiveStatistics
                machine = new DescriptiveStatistics(),
                big = new DescriptiveStatistics();
        int nIterations = its(1000, 10000);
        int prevPercent = -1, currPercent;
        for (int n = 0; n < nIterations; n++) {
            if (n == nIterations / 10) {

            }
            if ((currPercent = (int) (100.0 * n / nIterations)) != prevPercent) {
                prevPercent = currPercent;
                System.out.print(">");
                System.out.flush();
            }
            int nFactors = rndd.nextInt(2, 4);
            lIntegersModulo domain = new lIntegersModulo(getModulusRandom(rndd.nextInt(10, 30)));
            lMultivariatePolynomialZp[] factors = new lMultivariatePolynomialZp[nFactors];
            for (int i = 0; i < nFactors; i++)
                do {
                    factors[i] = RandomMultivariatePolynomial.randomPolynomial(2, rndd.nextInt(1, 4), rndd.nextInt(2, 5), domain, DegreeVector.LEX, rnd);
                    factors[i] = factors[i].divideDegreeVectorOrNull(factors[i].monomialContent());
                } while (factors[i].isConstant() || factors[i].nUsedVariables() == 1);

            lMultivariatePolynomialZp poly = factors[0].createOne().multiply(factors);
            if (!MultivariateSquareFreeFactorization.isSquareFree(poly)) {
                --n;
                continue;
            }

            poly = poly.divideDegreeVectorOrNull(poly.monomialContent());
            try {
                long start = System.nanoTime();
                FactorDecomposition<lMultivariatePolynomialZp> lFactorization = MultivariateFactorization.bivariateDenseFactorSquareFree(poly);
                machine.addValue(System.nanoTime() - start);
                FactorDecompositionTest.assertFactorization(poly, lFactorization);
                Assert.assertTrue(lFactorization.size() >= nFactors);

                start = System.nanoTime();
                FactorDecomposition<MultivariatePolynomial<BigInteger>> bFactorization = MultivariateFactorization.bivariateDenseFactorSquareFree(poly.toBigPoly());
                big.addValue(System.nanoTime() - start);
                FactorDecompositionTest.assertFactorization(poly.toBigPoly(), bFactorization);
                Assert.assertTrue(bFactorization.size() >= nFactors);
            } catch (Throwable t) {
                System.out.println(domain);
                System.out.println(Arrays.toString(factors));
                throw t;
            }
        }

        System.out.println();
        System.out.println("Machine integers: " + TimeUnits.statisticsNanotime(machine));
        System.out.println("Big integers    : " + TimeUnits.statisticsNanotime(big));
    }

    @Test
    public void test1() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(62653);
        String[] vars = {"a", "b"};
        lMultivariatePolynomialZp[] factors = {
                lMultivariatePolynomialZp.parse("17096+6578*a*b^2+54905*a^3", domain, vars),
                lMultivariatePolynomialZp.parse("43370+32368*a^2*b^2+45712*a^2*b^4+52302*a^4+23776*a^4*b^2", domain, vars)
        };
        lMultivariatePolynomialZp poly = factors[0].createOne().multiply(factors);
        FactorDecomposition<lMultivariatePolynomialZp> factorization = MultivariateFactorization.bivariateDenseFactorSquareFree(poly);
        FactorDecompositionTest.assertFactorization(poly, factorization);
        Assert.assertTrue(factorization.size() >= factors.length);
    }

//    @Test
//    public void test3() throws Exception {
//        lIntegersModulo domain = new lIntegersModulo(2);
//        lMultivariatePolynomialZp
//                a = lMultivariatePolynomialZp.parse("b*a^5 + a^4*b^2 + 11 + b^3", domain),
//                b = lMultivariatePolynomialZp.parse("a^6*b + 66*b + 17*b^2 + 1", domain),
//                c = lMultivariatePolynomialZp.parse("b^3*a^4 + a^4 + b", domain),
//                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
//                base = a.clone().multiply(b, c, d);
//        System.out.println(base);
//
////        System.out.println(base.asOverUnivariate(0));
//        UnivariatePolynomial<lUnivariatePolynomialZp> uBase =
//                base.asOverUnivariate(1).asUnivariate();
//
//        FiniteField<lUnivariatePolynomialZp> ff = new FiniteField<>(
//                IrreduciblePolynomials.randomIrreduciblePolynomial(domain.modulus, base.degreeSum() + 111, getRandom()));
//
//        System.out.println(ff.irreducible.toStringForCopy());
//        uBase = uBase.setDomain(ff);
//        FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>>
//                uFactors = Factorization.factor(uBase);
//
//        List<lMultivariatePolynomialZp> factors = new ArrayList<>();
//        for (int i = 0; i < uFactors.size(); i++) {
//            UnivariatePolynomial<lUnivariatePolynomialZp> uFactor = uFactors.get(i);
//            UnivariatePolynomial<lMultivariatePolynomialZp> factor =
//                    uFactor.mapElements(new MultivariatePolynomials<>(base),
//                            p -> AMultivariatePolynomial.asMultivariate(p, 2, 1, base.ordering));
//            factors.add(lMultivariatePolynomialZp.asMultivariate(factor, 0));
//        }
//
//
//        System.out.println(factors.size());
////        for (int i = 0; i < 1000; i++) {
////            long start = System.nanoTime();
////            Assert.assertEquals(4, MultivariateFactorization.bivariateDenseFactorSquareFree(base).size());
////            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
////        }
//    }

    @Test
    public void test3() throws Exception {
        lIntegersModulo domain = new lIntegersModulo(2);
        lMultivariatePolynomialZp
                a = lMultivariatePolynomialZp.parse("1 + b*a^5 + a^4*b^2 + 11 + b^3", domain),
                b = lMultivariatePolynomialZp.parse("a + a^2 + b + a^6*b + 66*b + 17*b^2 + 1", domain),
//                c = lMultivariatePolynomialZp.parse("b^3*a^4 + a^4 + b", domain),
//                d = lMultivariatePolynomialZp.parse("a^5 + b^5*a^5 + b^2 + 3", domain),
                base = a.clone().multiply(b);
//        System.out.println(base);

        System.out.println(MultivariateGCD.PolynomialGCD(base, base.derivative(0)));

        System.out.println(MultivariateSquareFreeFactorization.isSquareFree(base));

        FiniteField<lUnivariatePolynomialZp> ff = new FiniteField<>(
                IrreduciblePolynomials.randomIrreduciblePolynomial(domain.modulus, 4, getRandom()));


        System.out.println(MultivariateFactorization.bivariateDenseFactorSquareFree(base.mapCoefficients(ff, ff::valueOf)));
//        for (int i = 0; i < 1000; i++) {
//            long start = System.nanoTime();
//            Assert.assertEquals(4, MultivariateFactorization.bivariateDenseFactorSquareFree(base).size());
//            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
//        }
    }
}