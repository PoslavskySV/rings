package cc.r2.core.poly;

import cc.r2.core.poly.univar2.*;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class FiniteFieldTest {
    @Test
    public void test1() throws Exception {
        lUnivariatePolynomialZp irreducible = lUnivariatePolynomialZ.create(-1, -1, 0, 1).modulus(3);
        FiniteField domain = new FiniteField(irreducible);
        lUnivariatePolynomialZp poly = lUnivariatePolynomialZ.create(1, 1, 1, 1).modulus(3);
        lUnivariatePolynomialZp el = domain.valueOf(poly);
        Assert.assertTrue(domain.isOne(domain.multiply(el, domain.reciprocal(el))));
    }

    @Test
    public void test2() throws Exception {
        FiniteField domain = FiniteField.GF27;
        lUnivariatePolynomialZp
                c0 = domain.valueOf(lUnivariatePolynomialZ.create(1, 2, 3, 4, 5)),
                c1 = domain.valueOf(lUnivariatePolynomialZ.create(1, -2, 3, -4, -5)),
                c2 = domain.valueOf(lUnivariatePolynomialZ.create(11, 12, 13, 14, 15)),
                c3 = domain.add(c0, c1),
                c4 = domain.subtract(c1, c2),
                c5 = domain.multiply(c0, c1);

        UnivariatePolynomial<lUnivariatePolynomialZp>
                poly1 = UnivariatePolynomial.create(domain, c0, c1, c2, c3, c4, c5),
                poly2 = UnivariatePolynomial.create(domain, c5, c4, c3, c2, c1, c0),
                poly = poly1.clone().multiply(poly2).multiply(poly1.clone().add(poly2));

//        System.out.println(poly1.monic());
        FactorDecomposition<UnivariatePolynomial<lUnivariatePolynomialZp>> factors = Factorization.factorZp(poly);
        System.out.println(factors);
        System.out.println(poly);
        System.out.println(factors.toPolynomial());
    }
}