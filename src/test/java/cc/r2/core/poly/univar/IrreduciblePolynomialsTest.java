package cc.r2.core.poly.univar;

import cc.r2.core.bigint.BigInteger;
import cc.r2.core.poly.Domains;
import cc.r2.core.poly.test.APolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.univar.UnivariateDivision.InverseModMonomial;
import cc.r2.core.util.RandomDataGenerator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.poly.univar.IrreduciblePolynomials.irreducibleQ;
import static cc.r2.core.poly.univar.UnivariatePolynomialArithmetic.createMonomialMod;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class IrreduciblePolynomialsTest extends APolynomialTest {
    @Test
    public void testIrreducibleRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(5000, 1000);
        for (int i = 0; i < nIterations; i++) {
            long modulus = getModulusRandom(rndd.nextInt(5, 15));
            UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(5, 15), modulus, rnd);
            FactorDecomposition<UnivariatePolynomialZp64> factors = UnivariateFactorization.Factor(poly);
            try {
                Assert.assertEquals(factors.size() == 1, irreducibleQ(poly));
            } catch (Throwable e) {
                System.out.println("Modulus: " + poly.domain.modulus);
                System.out.println("Poly: " + poly.toStringForCopy());
                System.out.println("Factors: " + factors);
                System.out.println("Irred: " + irreducibleQ(poly));
                throw e;
            }
        }
    }

    @Test
    public void test1() throws Exception {
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(18, 92, 51, 36, 61, 93, 14, 13, 45, 11, 21, 79, 61, 1).modulus(97);
        Assert.assertFalse(irreducibleQ(poly));
    }

    @Test
    public void test2() throws Exception {
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(42, 73, 0, 79, 47, 1).modulus(89);
        Assert.assertTrue(UnivariateFactorization.Factor(poly).toString(), irreducibleQ(poly));
    }

    @Test
    public void test3() throws Exception {
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(952, 1768, 349, 1839, 1538, 1851, 941, 167, 1).modulus(1861);
        int exponent = 7;

        InverseModMonomial<UnivariatePolynomialZp64> invMod = UnivariateDivision.fastDivisionPreConditioning(poly);
        UnivariatePolynomialZp64 xq = createMonomialMod(poly.modulus(), poly, invMod);
        TIntObjectMap<UnivariatePolynomialZp64> cache = new TIntObjectHashMap<>();

        UnivariatePolynomialZp64 actual = IrreduciblePolynomials.composition(xq.clone(), exponent, poly, invMod, cache);
        UnivariatePolynomialZp64 expected = composition(xq.clone(), exponent, poly, invMod);
        UnivariatePolynomialZp64 expected0 = createMonomialMod(Domains.Z.pow(BigInteger.valueOf(poly.modulus()), exponent), poly, invMod);
        Assert.assertEquals(expected, actual);
        Assert.assertEquals(expected0, actual);
    }

    private static UnivariatePolynomialZp64 composition(UnivariatePolynomialZp64 xq, int n, UnivariatePolynomialZp64 poly, InverseModMonomial<UnivariatePolynomialZp64> invMod) {
        UnivariatePolynomialZp64 composition = xq.clone();
        for (int i = 1; i < n; i++)
            composition = ModularComposition.composition(composition, xq, poly, invMod);
        return composition;
    }

    @Test
    public void testCachedComposition1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(100, 1000);
        for (int i = 0; i < nIterations; i++) {
            if (i % 10 == 0)
                System.out.println(i);
            long modulus = getModulusRandom(rndd.nextInt(5, 15));
            UnivariatePolynomialZp64 poly = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(5, 10), modulus, rnd);

            InverseModMonomial<UnivariatePolynomialZp64> invMod = UnivariateDivision.fastDivisionPreConditioning(poly);
            UnivariatePolynomialZp64 xq = createMonomialMod(poly.modulus(), poly, invMod);
            TIntObjectMap<UnivariatePolynomialZp64> cache = new TIntObjectHashMap<>();

            int exponent = rndd.nextInt(1, 1024);
            UnivariatePolynomialZp64 actual = IrreduciblePolynomials.composition(xq.clone(), exponent, poly, invMod, cache);
            UnivariatePolynomialZp64 expected = composition(xq.clone(), exponent, poly, invMod);
            String msg = new StringBuilder()
                    .append("\npoly: " + poly.toStringForCopy())
                    .append("\nmodulus: " + modulus)
                    .append("\nexponent: " + exponent)
                    .toString();
            Assert.assertEquals(msg, expected, actual);
            if (exponent <= 64)
                Assert.assertEquals(msg, createMonomialMod(Domains.Z.pow(BigInteger.valueOf(poly.modulus()), exponent), poly, invMod), actual);

            for (int ci : cache.keys()) {
                msg = msg + "\nexponent: " + ci;
                Assert.assertEquals(msg, composition(xq.clone(), ci, poly, invMod), cache.get(ci));
                if (ci <= 64)
                    Assert.assertEquals(msg, createMonomialMod(Domains.Z.pow(poly.coefficientDomainCardinality(), ci), poly, invMod), cache.get(ci));
            }
        }
    }
}