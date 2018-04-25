package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.IParser;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * The ring of rationals (Q).
 *
 * @since 1.0
 */
public final class Rationals<E> implements Ring<Rational<E>> {
    private static final long serialVersionUID = 1L;
    /** Ring that numerator and denominator belongs to */
    public final Ring<E> ring;

    public Rationals(Ring<E> ring) {
        this.ring = ring;
    }

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public boolean isEuclideanRing() { return true; }

    @Override
    public BigInteger cardinality() {
        return null;
    }

    @Override
    public BigInteger characteristic() {
        return BigInteger.ZERO;
    }

    @Override
    public boolean isPerfectPower() {
        return false;
    }

    @Override
    public BigInteger perfectPowerBase() {
        return null;
    }

    @Override
    public BigInteger perfectPowerExponent() {
        return null;
    }

    @Override
    public Rational<E> add(Rational<E> a, Rational<E> b) {
        return a.add(b);
    }

    @Override
    public Rational<E> addMutable(Rational<E> a, Rational<E> b) {
        return a.addMutable(b);
    }

    @Override
    public Rational<E> subtract(Rational<E> a, Rational<E> b) {
        return a.subtract(b);
    }

    @Override
    public Rational<E> subtractMutable(Rational<E> a, Rational<E> b) {
        return a.subtractMutable(b);
    }

    @Override
    public Rational<E> multiply(Rational<E> a, Rational<E> b) {
        return a.multiply(b);
    }

    @Override
    public Rational<E> multiplyMutable(Rational<E> a, Rational<E> b) {
        return a.multiplyMutable(b);
    }

    @Override
    public Rational<E> negate(Rational<E> element) {
        return element.negate();
    }

    @Override
    public Rational<E> negateMutable(Rational<E> element) {
        return element.negateMutable();
    }

    @Override
    public int signum(Rational<E> element) {
        return element.signum();
    }

    @Override
    @SuppressWarnings("unchecked")
    public Rational<E>[] divideAndRemainder(Rational<E> dividend, Rational<E> divider) {
        return new Rational[]{dividend.divide(divider), Rational.zero(ring)};
    }

    @Override
    public Rational<E> divideExactMutable(Rational<E> dividend, Rational<E> divider) {
        return dividend.divideMutable(divider);
    }

    @Override
    public Rational<E> reciprocal(Rational<E> element) {
        return element.reciprocal();
    }

    @Override
    public Rational<E> gcd(Rational<E> a, Rational<E> b) {
        return Rational.one(ring);
    }

    @Override
    public FactorDecomposition<Rational<E>> factorSquareFree(Rational<E> element) {
        FactorDecomposition<Rational<E>> numFactors = ring.factorSquareFree(element.numerator)
                .mapTo(this, p -> new Rational<>(ring, p));
        FactorDecomposition<Rational<E>> denFactors = ring.factorSquareFree(element.denominator)
                .mapTo(this, p -> new Rational<>(ring, ring.getOne(), p));

        return numFactors.addAll(denFactors);
    }

    @Override
    public FactorDecomposition<Rational<E>> factor(Rational<E> element) {
        FactorDecomposition<Rational<E>> numFactors = ring.factor(element.numerator)
                .mapTo(this, p -> new Rational<>(ring, p));
        FactorDecomposition<Rational<E>> denFactors = ring.factor(element.denominator)
                .mapTo(this, p -> new Rational<>(ring, ring.getOne(), p));

        return numFactors.addAll(denFactors);
    }

    @Override
    public Rational<E> getZero() {
        return Rational.zero(ring);
    }

    @Override
    public Rational<E> getOne() {
        return Rational.one(ring);
    }

    @Override
    public boolean isZero(Rational element) {
        return element.isZero();
    }

    @Override
    public boolean isOne(Rational element) {
        return element.isOne();
    }

    @Override
    public boolean isUnit(Rational element) {
        return !isZero(element);
    }

    @Override
    public Rational<E> valueOf(long val) {
        return new Rational<>(ring, ring.valueOf(val));
    }

    @Override
    public Rational<E> valueOfBigInteger(BigInteger val) {
        return new Rational<>(ring, ring.valueOfBigInteger(val));
    }

    @Override
    public Rational<E> copy(Rational<E> element) {
        return new Rational<>(true, ring, ring.copy(element.numerator), ring.copy(element.denominator));
    }

    @Override
    public Rational<E> valueOf(Rational<E> val) {
        if (val.ring.equals(ring))
            return val;
        else return new Rational<>(ring, ring.valueOf(val.numerator), ring.valueOf(val.denominator));
    }

    @Override
    @SuppressWarnings("unchecked")
    public Rational<E>[][] createArray2d(int length) {
        return new Rational[length][];
    }

    @Override
    @SuppressWarnings("unchecked")
    public Rational<E>[][] createArray2d(int m, int n) {
        return new Rational[m][n];
    }

    @Override
    public int compare(Rational<E> o1, Rational<E> o2) {
        return o1.compareTo(o2);
    }

    @Override
    public Rational<E> getNegativeOne() {
        return Rational.one(ring).negateMutable();
    }

    @Override
    @SuppressWarnings("unchecked")
    public Rational<E>[] createArray(int length) {
        return new Rational[length];
    }

    @Override
    public Rational<E> randomElement(RandomGenerator rnd) {
        long den;
        E eden;
        do {den = rnd.nextInt();} while (ring.isZero(eden = ring.valueOf(den)));
        return new Rational<>(ring, ring.valueOf(rnd.nextInt()), eden);
    }

    public Rational<E> randomNonTrivialElement(RandomGenerator rnd) {
        E den;
        do {den = ring.randomElement(rnd);} while (ring.isZero(den));
        return new Rational<>(ring, ring.randomElement(rnd), den);
    }

    public Rational<E> parse(IParser<E> parser, String string) {
        int level = 0;
        int indexOfDiv = -1;
        for (int i = 0; i < string.length(); i++) {
            char c = string.charAt(i);
            if (c == '(')
                ++level;
            else if (c == ')')
                --level;
            if (level == 0 && c == '/')
                if (indexOfDiv != -1)
                    throw new IllegalArgumentException("not a valid rational");
                else indexOfDiv = i;
        }
        if (indexOfDiv == -1)
            return new Rational<>(ring, parser.parse(removeParenthesis(string)));
        return new Rational<>(ring,
                parser.parse(removeParenthesis(string.substring(0, indexOfDiv)).trim()),
                parser.parse(removeParenthesis(string.substring(indexOfDiv + 1)).trim()));
    }

    @Override
    public Rational<E> parse(String string) {
        return parse(ring, string);
    }

    private static String removeParenthesis(String string) {
        string = string.trim();
        if (!string.startsWith("(") || !string.endsWith(")"))
            return string;

        int level = 0;
        for (int i = 0; i < (string.length() - 1); i++) {
            char c = string.charAt(i);
            if (c == '(')
                ++level;
            else if (c == ')')
                --level;
            if (level == 0)
                return string;
        }
        return string.substring(1, string.length() - 1);
    }

    @Override
    public Iterator<Rational<E>> iterator() {
        throw new UnsupportedOperationException("Ring of infinite cardinality.");
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Rationals<?> rationals = (Rationals<?>) o;

        return ring.equals(rationals.ring);
    }

    @Override
    public int hashCode() {
        return ring.hashCode();
    }

    @Override
    public String toString() {
        return ring.equals(Rings.Z) ? "Q" : "Q[" + ring + "]";
    }
}
