package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.univar.UnivariateFactorization;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
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
        lUnivariatePolynomialZp irreducible = lUnivariatePolynomialZ.create(-1, -1, 0, 1).modulus(3);
        FiniteField<lUnivariatePolynomialZp> domain = new FiniteField<>(irreducible);
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(1, 1, 1, 1).modulus(3);
        lUnivariatePolynomialZp el = domain.valueOf(poly);
        Assert.assertTrue(domain.isOne(domain.multiply(el, domain.reciprocal(el))));
    }

    @Test
    public void test2() throws Exception {
        FiniteField<lUnivariatePolynomialZp> domain = FiniteField.GF27;
        lUnivariatePolynomialZp
                c0 = domain.valueOf(lUnivariatePolynomialZ.create(1, 2, 3, 4, 5).modulus(3)),
                c1 = domain.valueOf(lUnivariatePolynomialZ.create(1, -2, 3, -4, -5).modulus(3)),
                c2 = domain.valueOf(lUnivariatePolynomialZ.create(11, 12, 13, 14, 15).modulus(3)),
                c3 = domain.add(c0, c1),
                c4 = domain.subtract(c1, c2),
                c5 = domain.multiply(c0, c1);

        UnivariatePolynomial<lUnivariatePolynomialZp>
                poly1 = UnivariatePolynomial.create(domain, c0, c1, c2, c3, c4, c5),
                poly2 = UnivariatePolynomial.create(domain, c5, c4, c3, c2, c1, c0),
                poly = poly1.clone().multiply(poly2).multiply(poly1.clone().add(poly2));

//        System.out.println(poly1.monic());
        FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>> factors = UnivariateFactorization.factorInGF(poly);
        System.out.println(factors);
        System.out.println(poly);
        System.out.println(factors.toPolynomial());
    }

    @Test
    public void test3() throws Exception {
        lUnivariatePolynomialZp irreducible = lUnivariatePolynomialZ.create(1, 1, 1, 1, 1).modulus(2);
        FiniteField<lUnivariatePolynomialZp> field = new FiniteField<>(irreducible);
        Assert.assertEquals(16, field.cardinality().intValue());
        Assert.assertEquals(4, field.perfectPowerExponent().intValue());
    }

    @Test
    public void test4() throws Exception {
        IntegersModulo domain = new IntegersModulo(5);
        FiniteField<UnivariatePolynomial<BigInteger>> ff = new FiniteField<>(UnivariatePolynomial.create(domain, 1, 1, 1, 1));
        Set<UnivariatePolynomial<BigInteger>> set = new HashSet<>();
        for (UnivariatePolynomial<BigInteger> polynomial : ff)
            set.add(polynomial);

        Assert.assertEquals(ff.cardinality().intValue(), set.size());
    }

    @Test
    public void test5() throws Exception {
        FiniteField<lUnivariatePolynomialZp> ff = new FiniteField<>(lUnivariatePolynomialZp.create(7, new long[]{1, 1, 1, 1}));
        Set<lUnivariatePolynomialZp> set = new HashSet<>();
        for (lUnivariatePolynomialZp polynomial : ff)
            set.add(polynomial);

        Assert.assertEquals(ff.cardinality().intValue(), set.size());
    }
}