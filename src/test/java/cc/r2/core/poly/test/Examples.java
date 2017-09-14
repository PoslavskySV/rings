package cc.r2.core.poly.test;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domains;
import cc.r2.core.poly.univar.UnivariateGCD;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import org.junit.Test;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class Examples {
    @Test
    public void example1() throws Exception {

        UnivariatePolynomial<BigInteger> poly1 = UnivariatePolynomial.create(Domains.Z, 1, 3, 5, 7, 9, 11, 6);
        UnivariatePolynomial<BigInteger> poly2 = UnivariatePolynomial.create(Domains.Z, 6, 11, 9, 7, 5, 3, 1);
        UnivariatePolynomial<BigInteger> gcd = UnivariateGCD.PolynomialGCD(poly1, poly2);





    }
}
