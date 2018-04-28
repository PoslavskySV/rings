package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * Parent class for quotient rings
 */
public abstract class AQuotientRing<E> implements Ring<E> {
    /** the base ring */
    public final Ring<E> baseRing;

    protected AQuotientRing(Ring<E> baseRing) {
        this.baseRing = baseRing;
    }

    /** modulo operation */
    public abstract E mod(E el);

    @Override
    public boolean isField() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public boolean isEuclideanRing() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public BigInteger cardinality() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public BigInteger characteristic() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public boolean isPerfectPower() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public BigInteger perfectPowerBase() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public BigInteger perfectPowerExponent() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public E add(E a, E b) {
        return mod(baseRing.add(a, b));
    }

    @Override
    public E subtract(E a, E b) {
        return mod(baseRing.subtract(a, b));
    }

    @Override
    public E multiply(E a, E b) {
        return mod(baseRing.multiply(a, b));
    }

    @Override
    public E negate(E element) {
        return mod(baseRing.negate(element));
    }

    @Override
    public E copy(E element) {
        return baseRing.copy(element);
    }

    @Override
    public E[] divideAndRemainder(E dividend, E divider) {
        if (baseRing.isUnit(divider))
            return createArray(multiply(dividend, baseRing.reciprocal(divider)), getZero());
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public E reciprocal(E element) {
        if (isOne(element))
            return element;
        if (isMinusOne(element))
            return element;
        if (baseRing.isUnit(element))
            return valueOf(baseRing.reciprocal(element));
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public E getZero() {
        return baseRing.getZero();
    }

    @Override
    public E getOne() {
        return baseRing.getOne();
    }

    @Override
    public boolean isZero(E element) {
        return baseRing.isZero(element);
    }

    @Override
    public boolean isOne(E element) {
        return baseRing.isOne(element);
    }

    @Override
    public boolean isUnit(E element) {
        return baseRing.isUnit(element);
    }

    @Override
    public E valueOf(long val) {
        return mod(baseRing.valueOf(val));
    }

    @Override
    public E valueOfBigInteger(BigInteger val) {
        return mod(baseRing.valueOfBigInteger(val));
    }

    @Override
    public E valueOf(E val) {
        return mod(baseRing.valueOf(val));
    }

    @Override
    public Iterator<E> iterator() {
        throw new UnsupportedOperationException("Algebraic structure of ring is unknown");
    }

    @Override
    public int compare(E o1, E o2) {
        return baseRing.compare(o1, o2);
    }

    @Override
    public E parse(String string) {
        return valueOf(baseRing.parse(string));
    }

    @Override
    public E randomElement(RandomGenerator rnd) {
        return valueOf(baseRing.randomElement(rnd));
    }

    @Override
    public E randomElementTree(RandomGenerator rnd) {
        return valueOf(baseRing.randomElementTree(rnd));
    }
}
