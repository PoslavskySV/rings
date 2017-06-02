package cc.r2.core.poly.univar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.AbstractPolynomialTest;
import cc.r2.core.poly.FactorDecomposition;
import cc.r2.core.poly.univar.DivisionWithRemainder.InverseModMonomial;
import cc.r2.core.util.RandomDataGenerator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import static cc.r2.core.number.BigIntegerArithmetics.safePow;
import static cc.r2.core.poly.univar.IrreduciblePolynomials.irreducibleQ;
import static cc.r2.core.poly.univar.PolynomialArithmetics.createMonomialMod;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class IrreduciblePolynomialsTest extends AbstractPolynomialTest {
    @Test
    public void testIrreducibleRandom1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        int nIterations = its(5000, 1000);
        for (int i = 0; i < nIterations; i++) {
            long modulus = getModulusRandom(rndd.nextInt(5, 15));
            lUnivariatePolynomialZp poly = RandomPolynomials.randomMonicPoly(rndd.nextInt(5, 15), modulus, rnd);
            FactorDecomposition<lUnivariatePolynomialZp> factors = Factorization.factor(poly);
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
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(18, 92, 51, 36, 61, 93, 14, 13, 45, 11, 21, 79, 61, 1).modulus(97);
        Assert.assertFalse(irreducibleQ(poly));
    }

    @Test
    public void test2() throws Exception {
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(42, 73, 0, 79, 47, 1).modulus(89);
        Assert.assertTrue(Factorization.factor(poly).toString(), irreducibleQ(poly));
    }

    @Test
    public void test3() throws Exception {
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(952, 1768, 349, 1839, 1538, 1851, 941, 167, 1).modulus(1861);
        int exponent = 7;

        InverseModMonomial<lUnivariatePolynomialZp> invMod = DivisionWithRemainder.fastDivisionPreConditioning(poly);
        lUnivariatePolynomialZp xq = createMonomialMod(poly.modulus(), poly, invMod);
        TIntObjectMap<lUnivariatePolynomialZp> cache = new TIntObjectHashMap<>();

        lUnivariatePolynomialZp actual = IrreduciblePolynomials.composition(xq.clone(), exponent, poly, invMod, cache);
        lUnivariatePolynomialZp expected = composition(xq.clone(), exponent, poly, invMod);
        lUnivariatePolynomialZp expected0 = createMonomialMod(safePow(BigInteger.valueOf(poly.modulus()), exponent), poly, invMod);
        Assert.assertEquals(expected, actual);
        Assert.assertEquals(expected0, actual);
    }

    private static lUnivariatePolynomialZp composition(lUnivariatePolynomialZp xq, int n, lUnivariatePolynomialZp poly, InverseModMonomial<lUnivariatePolynomialZp> invMod) {
        lUnivariatePolynomialZp composition = xq.clone();
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
            lUnivariatePolynomialZp poly = RandomPolynomials.randomMonicPoly(rndd.nextInt(5, 10), modulus, rnd);

            InverseModMonomial<lUnivariatePolynomialZp> invMod = DivisionWithRemainder.fastDivisionPreConditioning(poly);
            lUnivariatePolynomialZp xq = createMonomialMod(poly.modulus(), poly, invMod);
            TIntObjectMap<lUnivariatePolynomialZp> cache = new TIntObjectHashMap<>();

            int exponent = rndd.nextInt(1, 1024);
            lUnivariatePolynomialZp actual = IrreduciblePolynomials.composition(xq.clone(), exponent, poly, invMod, cache);
            lUnivariatePolynomialZp expected = composition(xq.clone(), exponent, poly, invMod);
            String msg = new StringBuilder()
                    .append("\npoly: " + poly.toStringForCopy())
                    .append("\nmodulus: " + modulus)
                    .append("\nexponent: " + exponent)
                    .toString();
            Assert.assertEquals(msg, expected, actual);
            if (exponent <= 64)
                Assert.assertEquals(msg, createMonomialMod(safePow(BigInteger.valueOf(poly.modulus()), exponent), poly, invMod), actual);

            for (int ci : cache.keys()) {
                msg = msg + "\nexponent: " + ci;
                Assert.assertEquals(msg, composition(xq.clone(), ci, poly, invMod), cache.get(ci));
                if (ci <= 64)
                    Assert.assertEquals(msg, createMonomialMod(safePow(poly.coefficientDomainCardinality(), ci), poly, invMod), cache.get(ci));
            }
        }
    }
}