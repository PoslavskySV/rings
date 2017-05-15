package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.univar2.DivisionWithRemainder;
import cc.r2.core.poly.univar2.IUnivariatePolynomial;
import cc.r2.core.poly.univar2.UnivariateGCD;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariatePolynomials<Poly extends IUnivariatePolynomial<Poly>> implements Domain<Poly> {
    public final Poly factory;

    public UnivariatePolynomials(Poly factory) {
        this.factory = factory;
    }

    @Override
    public boolean isField() {return false;}

    @Override
    public BigInteger cardinality() {return null;}

    @Override
    public BigInteger characteristics() {return factory.coefficientDomainCharacteristics();}

    @Override
    public Poly add(Poly a, Poly b) {return a.clone().add(b);}

    @Override
    public Poly subtract(Poly a, Poly b) {return a.clone().subtract(b);}

    @Override
    public Poly multiply(Poly a, Poly b) {return a.clone().multiply(b);}

    @Override
    public Poly negate(Poly val) {return val.clone().negate();}

    @Override
    public int signum(Poly a) {return a.signum();}

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {return DivisionWithRemainder.divideAndRemainder(a, b, true);}

    @Override
    public Poly reciprocal(Poly a) {throw new UnsupportedOperationException();}

    @Override
    public Poly gcd(Poly a, Poly b) {return UnivariateGCD.PolynomialGCD(a, b);}

    @Override
    public Poly getZero() {return factory.createZero();}

    @Override
    public Poly getOne() {return factory.createOne();}

    @Override
    public boolean isZero(Poly poly) {return poly.isZero();}

    @Override
    public boolean isOne(Poly poly) {return poly.isOne();}

    @Override
    public Poly valueOf(long val) {return factory.createOne().multiply(val);}

    @Override
    public Poly valueOf(Poly val) {return val;}

    @Override
    public int compare(Poly o1, Poly o2) {return o1.compareTo(o2);}

    @Override
    public Poly randomElement(RandomGenerator rnd) {
        throw new UnsupportedOperationException();
    }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[][] createArray2d(int length) {
        Poly[] array = createArray(0);
        return (Poly[][]) Array.newInstance(array.getClass(), length);
    }

    @Override
    public Poly[][] createArray2d(int m, int n) {
        Poly[][] arr = createArray2d(m);
        for (int i = 0; i < arr.length; i++)
            arr[i] = createArray(n);
        return arr;
    }

    @Override
    @SuppressWarnings("unchecked")
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        UnivariatePolynomials<?> that = (UnivariatePolynomials<?>) o;
        return that.factory.getClass().equals(factory.getClass()) && factory.sameDomainWith((Poly) that.factory);
    }

    @Override
    public int hashCode() {
        return factory.hashCode();
    }
}
