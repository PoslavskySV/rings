package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.Rational;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Iterator;

/**
 * The domain of rationals (Q).
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Rationals implements Domain<Rational> {
    private static final long serialVersionUID = 1L;

    /** The domain of rationals (Q) */
    public static final Rationals Rationals = new Rationals();

    private Rationals() {}

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public BigInteger cardinality() {
        return null;
    }

    @Override
    public BigInteger characteristics() {
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
    public Rational add(Rational a, Rational b) {
        return a.add(b);
    }

    @Override
    public Rational subtract(Rational a, Rational b) {
        return a.subtract(b);
    }

    @Override
    public Rational multiply(Rational a, Rational b) {
        return a.multiply(b);
    }

    @Override
    public Rational negate(Rational element) {
        return element.negate();
    }

    @Override
    public int signum(Rational element) {
        return element.signum();
    }

    @Override
    public Rational[] divideAndRemainder(Rational dividend, Rational divider) {
        return new Rational[]{dividend.divide(divider), Rational.ZERO};
    }

    @Override
    public Rational reciprocal(Rational element) {
        return element.reciprocal();
    }

    @Override
    public Rational gcd(Rational a, Rational b) {
        return Rational.ONE;
    }

    @Override
    public Rational getZero() {
        return Rational.ZERO;
    }

    @Override
    public Rational getOne() {
        return Rational.ONE;
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
    public Rational valueOf(long val) {
        return new Rational(val);
    }

    @Override
    public Rational valueOfBigInteger(BigInteger val) {
        return new Rational(val);
    }

    @Override
    public Rational copy(Rational element) {
        return element;
    }

    @Override
    public Rational valueOf(Rational val) {
        return val;
    }

    @Override
    public Rational[][] createArray2d(int length) {
        return new Rational[length][];
    }

    @Override
    public Rational[][] createArray2d(int m, int n) {
        return new Rational[m][n];
    }

    @Override
    public int compare(Rational o1, Rational o2) {
        return o1.compareTo(o2);
    }

    @Override
    public Rational getNegativeOne() {
        return Rational.MINUS_ONE;
    }

    @Override
    public Rational[] createArray(int length) {
        return new Rational[length];
    }

    @Override
    public Rational randomElement(RandomGenerator rnd) {
        long den;
        do {den = rnd.nextInt();} while (den == 0);
        return new Rational(rnd.nextInt(), den);
    }

    @Override
    public Rational parse(String string) {
        string = string.trim();
        if (string.startsWith("(") && string.endsWith(")"))
            string = string.substring(1, string.length() - 1);
        String[] nd = string.split("\\/");
        if (nd.length == 1)
            return new Rational(new BigInteger(nd[0]));
        return new Rational(new BigInteger(nd[0]), new BigInteger(nd[1]));
    }

    @Override
    public Iterator<Rational> iterator() {
        throw new UnsupportedOperationException("Domain of infinite cardinality.");
    }
}
