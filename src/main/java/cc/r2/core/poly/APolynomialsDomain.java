package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.multivar.AMultivariatePolynomial;
import cc.r2.core.poly.univar.IUnivariatePolynomial;

import java.util.Arrays;
import java.util.Iterator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
abstract class APolynomialsDomain<Poly extends IPolynomial<Poly>> extends ADomain<Poly> implements PolynomialDomain<Poly> {
    private static final long serialVersionUID = 1L;

    public final Poly factory;

    APolynomialsDomain(Poly factory) {
        this.factory = factory.createZero();
    }

    @Override
    public final boolean isField() {return false;}

    @Override
    public final BigInteger cardinality() {return null;}

    @Override
    public final BigInteger characteristic() {return factory.coefficientDomainCharacteristics();}

    @Override
    public final Poly add(Poly a, Poly b) {return a.clone().add(b);}

    @Override
    public final Poly subtract(Poly a, Poly b) {return a.clone().subtract(b);}

    @Override
    public final Poly multiply(Poly a, Poly b) {return a.clone().multiply(b);}

    @Override
    public final Poly negate(Poly element) {return element.clone().negate();}

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
    public Poly negateMutable(Poly element) {
        return element.negate();
    }

    @Override
    public final Poly reciprocal(Poly element) {
        if (element.isConstant())
            return divideExact(getOne(), element);
        throw new ArithmeticException("not divisible: 1 / " + element);
    }

    @Override
    public final Poly getZero() {return factory.createZero();}

    @Override
    public final Poly getOne() {return factory.createOne();}

    @Override
    public final boolean isZero(Poly element) {return element.isZero();}

    @Override
    public final boolean isOne(Poly element) {return element.isOne();}

    @Override
    public boolean isUnit(Poly element) {
        return element.isOverField() ? element.isConstant() : isOne(element);
    }

    @Override
    public final Poly valueOf(long val) {return factory.createOne().multiply(val);}

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return factory.createOne().multiplyByBigInteger(val);
    }

    @Override
    public final Poly valueOf(Poly val) {
        if (factory.sameDomainWith(val))
            return val;
        else
            return val.setDomainFrom(factory);
    }

    @Override
    public Poly copy(Poly element) {
        return element.clone();
    }

    @Override
    public final int compare(Poly o1, Poly o2) {return o1.compareTo(o2);}

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
    public Poly parse(String string, String[] variables) {
        return factory.parsePoly(string, variables);
    }

    @Override
    public Iterator<Poly> iterator() {
        throw new UnsupportedOperationException("Domain of infinite cardinality.");
    }

    private int nVariables() {
        return factory instanceof IUnivariatePolynomial ? 1 : ((AMultivariatePolynomial) factory).nVariables;
    }

    @Override
    public String toString(String[] variables) {
        return toString(factory.coefficientDomainToString(), variables);
    }

    @Override
    public String toString(String coefficientDomain, String[] variables) {
        if (factory.isOverFiniteField())
            coefficientDomain = "(" + coefficientDomain + ")";
        return coefficientDomain + Arrays.toString(Arrays.copyOf(variables, nVariables()));
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(nVariables()));
    }
}
