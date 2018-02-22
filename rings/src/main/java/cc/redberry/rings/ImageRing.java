package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.StreamSupport;

/**
 *
 */
public class ImageRing<F, I> implements Ring<I> {
    /** the ring */
    public final Ring<F> ring;
    public final Function<F, I> imageFunc;
    public final Function<I, F> inverseFunc;

    public ImageRing(Ring<F> ring, Function<I, F> inverseFunc, Function<F, I> imageFunc) {
        this.ring = ring;
        this.inverseFunc = inverseFunc;
        this.imageFunc = imageFunc;
    }

    public I image(F el) {return imageFunc.apply(el);}

    public I[] image(F[] el) {
        I[] array = createArray(el.length);
        for (int i = 0; i < array.length; i++)
            array[i] = image(el[i]);
        return array;
    }

    public F inverse(I el) {return inverseFunc.apply(el);}

    public F[] inverse(I[] el) {
        F[] array = ring.createArray(el.length);
        for (int i = 0; i < array.length; i++)
            array[i] = inverse(el[i]);
        return array;
    }

    @Override
    public boolean isField() {
        return ring.isField();
    }

    @Override
    public boolean isEuclideanRing() {
        return ring.isEuclideanRing();
    }

    @Override
    public BigInteger cardinality() {
        return ring.cardinality();
    }

    @Override
    public BigInteger characteristic() {
        return ring.characteristic();
    }

    @Override
    public boolean isPerfectPower() {
        return ring.isPerfectPower();
    }

    @Override
    public BigInteger perfectPowerBase() {
        return ring.perfectPowerBase();
    }

    @Override
    public BigInteger perfectPowerExponent() {
        return ring.perfectPowerExponent();
    }

    @Override
    public I add(I a, I b) {
        return image(ring.add(inverse(a), inverse(b)));
    }

    @Override
    public I subtract(I a, I b) {
        return image(ring.subtract(inverse(a), inverse(b)));
    }

    @Override
    public I multiply(I a, I b) {
        return image(ring.multiply(inverse(a), inverse(b)));
    }

    @Override
    public I negate(I element) {
        return image(ring.negate(inverse(element)));
    }

    @Override
    public I increment(I element) {
        return image(ring.increment(inverse(element)));
    }

    @Override
    public I decrement(I element) {
        return image(ring.decrement(inverse(element)));
    }

    @Override
    public I add(I... elements) {
        return image(ring.add(inverse(elements)));
    }

    @Override
    public I multiply(I... elements) {
        return image(ring.multiply(inverse(elements)));
    }

    @Override
    public I abs(I el) {
        return image(ring.abs(inverse(el)));
    }

    @Override
    public I copy(I element) {
        return element;
    }

    @Override
    public I[] divideAndRemainder(I dividend, I divider) {
        return image(ring.divideAndRemainder(inverse(dividend), inverse(divider)));
    }

    @Override
    public I quotient(I dividend, I divider) {
        return image(ring.quotient(inverse(dividend), inverse(divider)));
    }

    @Override
    public I remainder(I dividend, I divider) {
        return image(ring.remainder(inverse(dividend), inverse(divider)));
    }

    @Override
    public I reciprocal(I element) {
        return image(ring.reciprocal(inverse(element)));
    }

    @Override
    public I getZero() {
        return image(ring.getZero());
    }

    @Override
    public I getOne() {
        return image(ring.getOne());
    }

    @Override
    public boolean isZero(I element) {
        return ring.isZero(inverse(element));
    }

    @Override
    public boolean isOne(I element) {
        return ring.isOne(inverse(element));
    }

    @Override
    public boolean isUnit(I element) {
        return ring.isUnit(inverse(element));
    }

    @Override
    public I valueOf(long val) {
        return image(ring.valueOf(val));
    }

    @Override
    public I valueOfBigInteger(BigInteger val) {
        return image(ring.valueOfBigInteger(val));
    }

    @Override
    public I valueOf(I val) {
        return image(ring.valueOf(inverse(val)));
    }

    @Override
    public Iterator<I> iterator() {
        return StreamSupport
                .stream(Spliterators.spliteratorUnknownSize(ring.iterator(), Spliterator.ORDERED), false)
                .map(this::image).iterator();
    }

    @Override
    public I gcd(I a, I b) {
        return image(ring.gcd(inverse(a), inverse(b)));
    }

    @Override
    public I[] extendedGCD(I a, I b) {
        return image(ring.extendedGCD(inverse(a), inverse(b)));
    }

    @Override
    public I lcm(I a, I b) {
        return image(ring.lcm(inverse(a), inverse(b)));
    }

    @Override
    public I gcd(I... elements) {
        return image(ring.gcd(inverse(elements)));
    }

    @Override
    public I gcd(Iterable<I> elements) {
        return image(ring.gcd(() -> StreamSupport.stream(elements.spliterator(), false).map(this::inverse).iterator()));
    }

    @Override
    public int signum(I element) {
        return ring.signum(inverse(element));
    }


    @Override
    public FactorDecomposition<I> factorSquareFree(I element) {
        return ring.factorSquareFree(inverse(element)).mapTo(this, this::image);
    }

    @Override
    public FactorDecomposition<I> factor(I element) {
        return ring.factor(inverse(element)).mapTo(this, this::image);
    }

    @Override
    public I parse(String string) {
        return image(ring.parse(string));
    }

    @Override
    public I pow(I base, int exponent) {
        return image(ring.pow(inverse(base), exponent));
    }

    @Override
    public I pow(I base, long exponent) {
        return image(ring.pow(inverse(base), exponent));
    }

    @Override
    public I pow(I base, BigInteger exponent) {
        return null;
    }

    @Override
    public I factorial(long num) {
        return image(ring.factorial(num));
    }

    @Override
    public I randomElement(RandomGenerator rnd) {
        return image(ring.randomElement(rnd));
    }

    @Override
    public int compare(I o1, I o2) {
        return ring.compare(inverse(o1), inverse(o2));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        ImageRing<?, ?> that = (ImageRing<?, ?>) o;

        if (!ring.equals(that.ring)) return false;
        if (!inverseFunc.equals(that.inverseFunc)) return false;
        return imageFunc.equals(that.imageFunc);
    }

    @Override
    public int hashCode() {
        int result = ring.hashCode();
        result = 31 * result + inverseFunc.hashCode();
        result = 31 * result + imageFunc.hashCode();
        return result;
    }
}
