package cc.r2.core.poly.multivar;

import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.FactorDecompositionTest;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.util.RandomDataGenerator;
import cc.r2.core.util.TimeUnits;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import static cc.r2.core.poly.multivar.MultivariateFactorization.LIFT;
import static cc.r2.core.poly.multivar.MultivariateFactorization.RECO;
import static cc.r2.core.poly.multivar.MultivariateFactorization.UNI;

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
//        base = base.shift(1, 1);
//        System.out.println(MultivariateFactorization.bivariateFactorSquareFree(base));
        for (int i = 0; i < 1000; i++) {
            UNI = LIFT = RECO = 0;
            long start = System.nanoTime();
            Assert.assertEquals(4, MultivariateFactorization.bivariateFactorSquareFree(base).size());
            System.out.println(TimeUnits.nanosecondsToString(System.nanoTime() - start));
            System.out.println(TimeUnits.nanosecondsToString(UNI));
            System.out.println(TimeUnits.nanosecondsToString(LIFT));
            System.out.println(TimeUnits.nanosecondsToString(RECO));
            System.out.println();
        }
    }

    @Test
    public void testBivariateFactorizationRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();

        int nIterations = its(100, 1000);
        int prevPercent = -1, currPercent;
        for (int n = 0; n < nIterations; n++) {
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
            FactorDecomposition<lMultivariatePolynomialZp> factorization = MultivariateFactorization.bivariateFactorSquareFree(poly);
            FactorDecompositionTest.assertFactorization(poly, factorization);
            Assert.assertTrue(factorization.size() >= nFactors);
        }
    }
}