package cc.r2.core.poly;

import cc.r2.core.bigint.BigInteger;
import cc.r2.core.poly.univar.UnivariateFactorization;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariatePolynomialZ64;
import cc.r2.core.poly.univar.UnivariatePolynomialZp64;
import org.junit.Assert;
import org.junit.Test;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class FiniteFieldTest {
    @Test
    public void test1() throws Exception {
        UnivariatePolynomialZp64 irreducible = UnivariatePolynomialZ64.create(-1, -1, 0, 1).modulus(3);
        FiniteField<UnivariatePolynomialZp64> domain = new FiniteField<>(irreducible);
        UnivariatePolynomialZp64 poly = UnivariatePolynomialZ64.create(1, 1, 1, 1).modulus(3);
        UnivariatePolynomialZp64 el = domain.valueOf(poly);
        Assert.assertTrue(domain.isOne(domain.multiply(el, domain.reciprocal(el))));
    }

    @Test
    public void test2() throws Exception {
        FiniteField<UnivariatePolynomialZp64> domain = FiniteField.GF27;
        UnivariatePolynomialZp64
                c0 = domain.valueOf(UnivariatePolynomialZ64.create(1, 2, 3, 4, 5).modulus(3)),
                c1 = domain.valueOf(UnivariatePolynomialZ64.create(1, -2, 3, -4, -5).modulus(3)),
                c2 = domain.valueOf(UnivariatePolynomialZ64.create(11, 12, 13, 14, 15).modulus(3)),
                c3 = domain.add(c0, c1),
                c4 = domain.subtract(c1, c2),
                c5 = domain.multiply(c0, c1);

        UnivariatePolynomial<UnivariatePolynomialZp64>
                poly1 = UnivariatePolynomial.create(domain, c0, c1, c2, c3, c4, c5),
                poly2 = UnivariatePolynomial.create(domain, c5, c4, c3, c2, c1, c0),
                poly = poly1.clone().multiply(poly2).multiply(poly1.clone().add(poly2));

//        System.out.println(poly1.monic());
        FactorDecomposition<UnivariatePolynomial<UnivariatePolynomialZp64>> factors = UnivariateFactorization.FactorInGF(poly);
        System.out.println(factors);
        System.out.println(poly);
        System.out.println(factors.toPolynomial());
    }

    @Test
    public void test3() throws Exception {
        UnivariatePolynomialZp64 irreducible = UnivariatePolynomialZ64.create(1, 1, 1, 1, 1).modulus(2);
        FiniteField<UnivariatePolynomialZp64> field = new FiniteField<>(irreducible);
        Assert.assertEquals(16, field.cardinality().intValue());
        Assert.assertEquals(4, field.perfectPowerExponent().intValue());
    }

    @Test
    public void test4() throws Exception {
        IntegersZp domain = new IntegersZp(5);
        FiniteField<UnivariatePolynomial<BigInteger>> ff = new FiniteField<>(UnivariatePolynomial.create(domain, 1, 1, 1, 1));
        Set<UnivariatePolynomial<BigInteger>> set = new HashSet<>();
        for (UnivariatePolynomial<BigInteger> polynomial : ff)
            set.add(polynomial);

        Assert.assertEquals(ff.cardinality().intValue(), set.size());
    }

    @Test
    public void test5() throws Exception {
        FiniteField<UnivariatePolynomialZp64> ff = new FiniteField<>(UnivariatePolynomialZp64.create(7, new long[]{1, 1, 1, 1}));
        Set<UnivariatePolynomialZp64> set = new HashSet<>();
        for (UnivariatePolynomialZp64 polynomial : ff)
            set.add(polynomial);

        Assert.assertEquals(ff.cardinality().intValue(), set.size());
    }
}