package cc.r2.core.poly.generics;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.univar.DivisionWithRemainder;
import cc.r2.core.poly.univar.IMutablePolynomial;
import cc.r2.core.poly.univar.PolynomialGCD;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class UnivariatePolynomialDomain<T extends IMutablePolynomial<T>> implements Domain<T> {
    public final T factory;

    public UnivariatePolynomialDomain(T factory) {
        this.factory = factory;
    }

    @Override
    public boolean isField() {return false;}

    @Override
    public BigInteger size() {return null;}

    @Override
    public T add(T a, T b) {
        return a.clone().add(b);
    }

    @Override
    public T subtract(T a, T b) {
        return a.clone().subtract(b);
    }

    @Override
    public T multiply(T a, T b) {
        return a.clone().multiply(b);
    }

    @Override
    public T negate(T val) {
        return val.clone().negate();
    }

    @Override
    public int signum(T a) {
        return a.signum();
    }

    @Override
    public T[] divideAndRemainder(T a, T b) {
        return DivisionWithRemainder.divideAndRemainder(a, b, true);
    }

    @Override
    public T reciprocal(T a) {
        throw new UnsupportedOperationException();
    }

    @Override
    public T gcd(T a, T b) {
        return PolynomialGCD.PolynomialGCD(a, b);
    }

    @Override
    public T getZero() {
        return factory.createZero();
    }

    @Override
    public T getOne() {
        return factory.createOne();
    }

    @Override
    public boolean isZero(T t) {
        return t.isZero();
    }

    @Override
    public boolean isOne(T t) {
        return t.isOne();
    }

    @Override
    public T valueOf(long val) {
        return factory.createOne().multiply(val);
    }

    @Override
    public T valueOf(T val) {
        return val;
    }

    @Override
    public int compare(T o1, T o2) {
        return o1.compareTo(o2);
    }

    @Override
    public Domain<T> getExtension() {return null;}
}
