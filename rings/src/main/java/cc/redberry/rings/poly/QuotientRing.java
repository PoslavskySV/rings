package cc.redberry.rings.poly;

import cc.redberry.rings.ARing;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.IStringifier;
import cc.redberry.rings.poly.multivar.AMonomial;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.Ideal;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * Multivariate quotient ring
 */
public class QuotientRing<Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        extends ARing<Poly> implements IPolynomialRing<Poly> {
    private static final long serialVersionUID = 1L;
    /** the base ring */
    public final MultivariateRing<Poly> baseRing;
    /** the ideal */
    public final Ideal<Term, Poly> ideal;
    /** factory element */
    private final Poly factory;

    public QuotientRing(MultivariateRing<Poly> baseRing, Ideal<Term, Poly> ideal) {
        this.baseRing = baseRing;
        this.ideal = ideal;
        this.factory = ideal.getBasisGenerator(0).createZero();
    }

    @Override
    public int nVariables() {
        return factory.nVariables;
    }

    @Override
    public Poly factory() {
        return factory;
    }

    @Override
    public Poly variable(int variable) {
        return factory.createMonomial(variable, 1);
    }

    @Override
    public boolean isField() {
        return factory.isOverField() && ideal.isMaximal();
    }

    @Override
    public boolean isEuclideanRing() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public BigInteger cardinality() {
        return factory.coefficientRingCardinality().isZero()
                ? BigInteger.ZERO
                : ideal.dimension() != 0 ? null : BigInteger.valueOf(ideal.degree()).multiply(factory.coefficientRingCardinality());
    }

    @Override
    public BigInteger characteristic() {
        return factory.coefficientRingCharacteristic();
    }

    public Poly mod(Poly el) {
        return ideal.normalForm(el);
    }

    public Poly normalForm(Poly el) {
        return mod(el);
    }

    @Override
    public Poly add(Poly a, Poly b) {
        return mod(baseRing.add(a, b));
    }

    @Override
    public Poly subtract(Poly a, Poly b) {
        return mod(baseRing.subtract(a, b));
    }

    @Override
    public Poly multiply(Poly a, Poly b) {
        return mod(baseRing.multiply(a, b));
    }

    @Override
    public Poly negate(Poly element) {
        return mod(baseRing.negate(element));
    }

    @Override
    public Poly copy(Poly element) {
        return baseRing.copy(element);
    }

    @Override
    public Poly[] divideAndRemainder(Poly dividend, Poly divider) {
        if (baseRing.isUnit(divider))
            return createArray(multiply(dividend, baseRing.reciprocal(divider)), getZero());
        if (isField())
            return createArray(multiply(dividend, reciprocal(divider)), getZero());
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public Poly reciprocal(Poly element) {
        if (isOne(element))
            return element;
        if (isMinusOne(element))
            return element;
        if (baseRing.isUnit(element))
            return valueOf(baseRing.reciprocal(element));
        if (isField()) {
            if (!element.isConstant())
                element = mod(element);
            if (!element.isConstant())
                throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
            return baseRing.getOne().divideByLC(element);
        }
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public Poly getZero() {
        return baseRing.getZero();
    }

    @Override
    public Poly getOne() {
        return baseRing.getOne();
    }

    @Override
    public boolean isZero(Poly element) {
        return baseRing.isZero(element);
    }

    @Override
    public boolean isOne(Poly element) {
        return baseRing.isOne(element);
    }

    @Override
    public boolean isUnit(Poly element) {
        return baseRing.isUnit(element);
    }

    @Override
    public Poly valueOf(long val) {
        return mod(baseRing.valueOf(val));
    }

    @Override
    public Poly valueOfBigInteger(BigInteger val) {
        return mod(baseRing.valueOfBigInteger(val));
    }

    @Override
    public Poly valueOf(Poly val) {
        return mod(baseRing.valueOf(val));
    }

    @Override
    public Iterator<Poly> iterator() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public int compare(Poly o1, Poly o2) {
        return baseRing.compare(o1, o2);
    }

    @Override
    public Poly parse(String string) {
        return valueOf(baseRing.parse(string));
    }

    @Override
    public Poly randomElement(RandomGenerator rnd) {
        return valueOf(baseRing.randomElement(rnd));
    }

    @Override
    public Poly randomElementTree(RandomGenerator rnd) {
        return valueOf(baseRing.randomElementTree(rnd));
    }

    @Override
    public String toString(IStringifier<Poly> stringifier) {
        return baseRing.toString(stringifier) + "/<" + ideal.toString(stringifier) + ">";
    }

    @Override
    public String toString() {
        return toString(IStringifier.dummy());
    }
}
