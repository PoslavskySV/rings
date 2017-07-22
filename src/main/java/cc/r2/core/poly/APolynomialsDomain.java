package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;

import java.lang.reflect.Array;
import java.util.Iterator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
abstract class APolynomialsDomain<Poly extends IGeneralPolynomial<Poly>> extends ADomain<Poly> {
    public final Poly factory;

    APolynomialsDomain(Poly factory) {
        this.factory = factory.createZero();
    }

    @Override
    public final boolean isField() {return false;}

    @Override
    public final BigInteger cardinality() {return null;}

    @Override
    public final BigInteger characteristics() {return factory.coefficientDomainCharacteristics();}

    @Override
    public final Poly add(Poly a, Poly b) {return a.clone().add(b);}

    @Override
    public final Poly subtract(Poly a, Poly b) {return a.clone().subtract(b);}

    @Override
    public final Poly multiply(Poly a, Poly b) {return a.clone().multiply(b);}

    @Override
    public final Poly negate(Poly val) {return val.clone().negate();}

    @Override
    public Poly addMutable(Poly a, Poly b) {
        return a.add(b);
    }

    @Override
    public Poly subtractMutable(Poly a, Poly b) {
        return a.subtract(b);
    }

    @Override
    public Poly multiplyMutable(Poly a, Poly b) {
        return a.multiply(b);
    }

    @Override
    public Poly negateMutable(Poly val) {
        return val.negate();
    }

    @Override
    public final int signum(Poly a) {return a.signum();}

    @Override
    public final Poly reciprocal(Poly a) {
        if (a.isConstant())
            return divideExact(getOne(), a);
        throw new ArithmeticException("not divisible: 1 / " + a);
    }

    @Override
    public final Poly getZero() {return factory.createZero();}

    @Override
    public final Poly getOne() {return factory.createOne();}

    @Override
    public final boolean isZero(Poly poly) {return poly.isZero();}

    @Override
    public final boolean isOne(Poly poly) {return poly.isOne();}

    @Override
    public boolean isUnit(Poly poly) {
        return poly.isOverField() ? poly.isConstant() : isOne(poly);
    }

    @Override
    public final Poly valueOf(long val) {return factory.createOne().multiply(val);}

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return factory.createOne().multiplyByBigInteger(val);
    }

    @Override
    public final Poly valueOf(Poly val) {return val;}

    @Override
    public Poly copy(Poly element) {
        return element.clone();
    }

    @Override
    public final int compare(Poly o1, Poly o2) {return o1.compareTo(o2);}

    @Override
    @SuppressWarnings("unchecked")
    public final Poly[][] createArray2d(int length) {
        Poly[] array = createArray(0);
        return (Poly[][]) Array.newInstance(array.getClass(), length);
    }

    @Override
    public final Poly[][] createArray2d(int m, int n) {
        Poly[][] arr = createArray2d(m);
        for (int i = 0; i < arr.length; i++)
            arr[i] = createArray(n);
        return arr;
    }

    @Override
    @SuppressWarnings("unchecked")
    public final boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        APolynomialsDomain<?> that = (APolynomialsDomain<?>) o;
        return that.factory.getClass().equals(factory.getClass()) && factory.sameDomainWith((Poly) that.factory);
    }

    @Override
    public final int hashCode() {
        return factory.hashCode();
    }

    @Override
    public Poly parse(String string) {
        return factory.parsePoly(string);
    }

    @Override
    public Iterator<Poly> iterator() {
        throw new UnsupportedOperationException("Domain of infinite cardinality.");
    }
}
