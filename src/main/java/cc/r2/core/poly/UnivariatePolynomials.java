package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.univar.DivisionWithRemainder;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.RandomPolynomials;
import cc.r2.core.poly.univar.UnivariateGCD;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariatePolynomials<Poly extends IUnivariatePolynomial<Poly>> extends APolynomialsDomain<Poly> {
    public UnivariatePolynomials(Poly factory) {
        super(factory);
    }

    @Override
    public Poly remainder(Poly a, Poly b) {return DivisionWithRemainder.remainder(a, b, true);}

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {return DivisionWithRemainder.divideAndRemainder(a, b, true);}

    @Override
    public Poly gcd(Poly a, Poly b) {return UnivariateGCD.PolynomialGCD(a, b);}

//    @Override
//    public Poly randomElement(RandomGenerator rnd) {
//        return RandomPolynomials.randomPoly(32, factory.);
//    }
}
