package cc.redberry.rings.poly;

import cc.redberry.rings.ARing;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.IStringifier;

import java.util.Iterator;

/**
 * @since 1.0
 */
abstract class APolynomialRing<Poly extends IPolynomial<Poly>> extends ARing<Poly> implements IPolynomialRing<Poly> {
    private static final long serialVersionUID = 1L;

    /** the factory polynomial */
    final Poly factory;

    APolynomialRing(Poly factory) {
        this.factory = factory.createZero();
    }

    @Override
    public Poly factory() {return factory;}

    @Override
    public boolean isEuclideanRing() {return factory.isOverField();}

    @Override
    public final boolean isField() {return false;}

    @Override
    public final BigInteger cardinality() {return null;}

    @Override
    public final BigInteger characteristic() {return factory.coefficientRingCharacteristic();}

    @Override
    public final Poly add(Poly a, Poly b) {return a.clone().add(b);}

    @Override
    public final Poly subtract(Poly a, Poly b) {return a.clone().subtract(b);}

    @Override
    public final Poly multiply(Poly a, Poly b) {return a.clone().multiply(b);}

    @Override
    public final Poly negate(Poly element) {return element.clone().negate();}

    @Override
    public Poly pow(Poly base, BigInteger exponent) {
        return PolynomialMethods.polyPow(base, exponent, true);
    }

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
        return element.isOverField() ? element.isConstant() : (isOne(element) || isMinusOne(element));
    }

    @Override
    public final Poly valueOf(long val) {return factory.createOne().multiply(val);}

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return factory.createOne().multiplyByBigInteger(val);
    }

    @Override
    public final Poly valueOf(Poly val) {
        if (factory.sameCoefficientRingWith(val))
            return val;
        else
            return val.setCoefficientRingFrom(factory);
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

        APolynomialRing<?> that = (APolynomialRing<?>) o;
        return that.factory.getClass().equals(factory.getClass()) && factory.sameCoefficientRingWith((Poly) that.factory);
    }

    @Override
    public final int hashCode() {
        return getOne().hashCode();
    }

    @Override
    public Iterator<Poly> iterator() {
        throw new UnsupportedOperationException("Ring of infinite cardinality.");
    }

    @Override
    public Poly parse(String string) {
        return factory.parsePoly(string);
    }

    @Override
    public String toString(IStringifier<Poly> stringifier) {
        String cfRing = factory.coefficientRingToString(stringifier);
        if (cfRing.length() > 2)
            cfRing = "(" + cfRing + ")";
        StringBuilder sb = new StringBuilder();
        sb.append(cfRing);
        sb.append("[");
        int nVars = nVariables();
        for (int i = 0; ; ++i) {
            sb.append(stringifier.getBinding(variable(i), IStringifier.defaultVar(i, nVars)));
            if (i == nVars - 1)
                break;
            sb.append(", ");
        }
        sb.append("]");
        return sb.toString();
    }

    public String toString(String... variables) {
        return toString(IStringifier.mkPolyStringifier(factory, variables));
    }

    @Override
    public String toString() {
        return toString(IStringifier.defaultVars(nVariables()));
    }
}
