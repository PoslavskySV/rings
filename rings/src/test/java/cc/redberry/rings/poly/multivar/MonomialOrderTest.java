package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MonomialOrder.*;
import cc.redberry.rings.util.ArraysUtil;
import cc.redberry.rings.util.ZipUtil;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

import static cc.redberry.rings.poly.multivar.MonomialOrder.*;

/**
 *
 */
public class MonomialOrderTest extends AMultivariateTest {
    @Test
    public void test1() throws Exception {
        RandomGenerator random = getRandom();

        int[] permutation = ArraysUtil.sequence(0, 7);
        for (int i = 0; i < 100; ++i) {
            ArraysUtil.shuffle(permutation, random);
            int[] inversePermutation = MultivariateGCD.inversePermutation(permutation);

            MultivariatePolynomial<BigInteger> poly = RandomMultivariatePolynomials.randomPolynomial(permutation.length, 10, 20, random);
            poly = poly.setOrdering(MonomialOrder.GREVLEX);

            MultivariatePolynomial<BigInteger> renamed = AMultivariatePolynomial.renameVariables(poly, permutation);
            Monomial<BigInteger>[] rArray = renamed.toArray();
            for (int j = 0; j < rArray.length; j++)
                rArray[j] = AMultivariatePolynomial.renameVariables(rArray[j], inversePermutation);

            Monomial<BigInteger>[] pArray = poly.setOrdering(new GrevLexWithPermutation(permutation)).toArray();
            Assert.assertArrayEquals(pArray, rArray);
        }
    }

    @Test
    public void test2() {
        for (Object o : Arrays.asList(LEX, GRLEX, ALEX, GREVLEX))
            Assert.assertEquals(o, ZipUtil.uncompress(ZipUtil.compress(o)));
    }
}