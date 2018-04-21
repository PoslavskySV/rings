package cc.redberry.rings.poly.univar;

import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.poly.multivar.MonomialOrder;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.ArraysUtil;

import java.util.*;
import java.util.function.*;
import java.util.stream.Collector;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import static cc.redberry.rings.bigint.BigInteger.ONE;
import static cc.redberry.rings.bigint.BigInteger.ZERO;
import static cc.redberry.rings.bigint.BigIntegerUtil.abs;

/**
 * Univariate polynomial over generic ring.
 *
 * @since 1.0
 */
public final class UnivariatePolynomial<E> implements IUnivariatePolynomial<UnivariatePolynomial<E>>, Iterable<E> {
    private static final long serialVersionUID = 1L;
    /** The coefficient ring */
    public final Ring<E> ring;
    /** list of coefficients { x^0, x^1, ... , x^degree } */
    E[] data;
    /** points to the last non zero element in the data array */
    int degree;

    private UnivariatePolynomial(Ring<E> ring, E[] data, int degree) {
        this.ring = ring;
        this.data = data;
        this.degree = degree;
        assert data.length > 0;
    }

    private UnivariatePolynomial(Ring<E> ring, E[] data) {
        this(ring, data, data.length - 1);
        fixDegree();
    }

    /**
     * Parse string into polynomial
     */
    public static <E> UnivariatePolynomial<E> parse(String string, Ring<E> ring) {
        return Parser.parse(ring, ring, string);
    }

    /**
     * Creates new univariate polynomial over specified ring with the specified coefficients. Note: the array {@code
     * data} will not be copied.
     *
     * @param ring the ring
     * @param data the coefficients
     * @return new univariate polynomial over specified ring with specified coefficients
     */
    @SuppressWarnings("unchecked")
    public static <E> UnivariatePolynomial<E> create(Ring<E> ring, E... data) {
        ring.setToValueOf(data);
        return new UnivariatePolynomial<>(ring, data);
    }

    /** skips {@code ring.setToValueOf(data)} */
    public static <E> UnivariatePolynomial<E> createUnsafe(Ring<E> ring, E[] data) {
        return new UnivariatePolynomial<>(ring, data);
    }

    /**
     * Creates univariate polynomial over specified ring (with integer elements) with the specified coefficients
     *
     * @param ring the ring
     * @param data the coefficients
     * @return new univariate polynomial over specified ring with specified coefficients
     */
    public static UnivariatePolynomial<BigInteger> create(Ring<BigInteger> ring, long... data) {
        return create(ring, ring.valueOf(data));
    }

    /**
     * Creates new univariate Z[x] polynomial
     *
     * @param data the data
     * @return new univariate Z[x] polynomial with specified coefficients
     */
    public static UnivariatePolynomial<BigInteger> create(long... data) {
        return create(Rings.Z, data);
    }

    /**
     * Creates constant polynomial over specified ring
     *
     * @param ring     the ring
     * @param constant the value
     * @return constant polynomial over specified ring
     */
    @SuppressWarnings("unchecked")
    public static <E> UnivariatePolynomial<E> constant(Ring<E> ring, E constant) {
        return create(ring, ring.createArray(constant));
    }

    /**
     * Creates zero polynomial over specified ring
     *
     * @param ring the ring
     * @return zero polynomial over specified ring
     */
    public static <E> UnivariatePolynomial<E> zero(Ring<E> ring) {
        return constant(ring, ring.getZero());
    }

    /**
     * Creates unit polynomial over specified ring
     *
     * @param ring the ring
     * @return unit polynomial over specified ring
     */
    public static <E> UnivariatePolynomial<E> one(Ring<E> ring) {
        return constant(ring, ring.getOne());
    }

    /**
     * Converts poly over BigIntegers to machine-sized polynomial in Z
     *
     * @param poly the polynomial over BigIntegers
     * @return machine-sized polynomial in Z
     * @throws ArithmeticException if some of coefficients is out of long range
     */
    public static UnivariatePolynomialZ64 asOverZ64(UnivariatePolynomial<BigInteger> poly) {
        long[] data = new long[poly.degree + 1];
        for (int i = 0; i < data.length; i++)
            data[i] = poly.data[i].longValueExact();
        return UnivariatePolynomialZ64.create(data);
    }

    /**
     * Converts Zp[x] poly over BigIntegers to machine-sized polynomial in Zp
     *
     * @param poly the Z/p polynomial over BigIntegers
     * @return machine-sized polynomial in Z/p
     * @throws IllegalArgumentException if {@code poly.ring} is not {@link IntegersZp}
     * @throws ArithmeticException      if some of {@code poly} elements is out of long range
     */
    public static UnivariatePolynomialZp64 asOverZp64(UnivariatePolynomial<BigInteger> poly) {
        if (!(poly.ring instanceof IntegersZp))
            throw new IllegalArgumentException("Not a modular ring: " + poly.ring);
        long[] data = new long[poly.degree + 1];
        for (int i = 0; i < data.length; i++)
            data[i] = ((BigInteger) poly.data[i]).longValueExact();
        return UnivariatePolynomialZp64.create(((IntegersZp) poly.ring).modulus.longValueExact(), data);
    }

    /**
     * Converts Zp[x] polynomial to Z[x] polynomial formed from the coefficients of this represented in symmetric
     * modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @param poly Zp polynomial
     * @return Z[x] version of the poly with coefficients represented in symmetric modular form ({@code -modulus/2 <=
     * cfx <= modulus/2}).
     * @throws IllegalArgumentException is {@code poly.ring} is not a {@link IntegersZp}
     */
    public static UnivariatePolynomial<BigInteger> asPolyZSymmetric(UnivariatePolynomial<BigInteger> poly) {
        if (!(poly.ring instanceof IntegersZp))
            throw new IllegalArgumentException("Not a modular ring: " + poly.ring);
        IntegersZp ring = (IntegersZp) poly.ring;
        BigInteger[] newData = new BigInteger[poly.degree + 1];
        for (int i = poly.degree; i >= 0; --i)
            newData[i] = ring.symmetricForm(poly.data[i]);
        return UnivariatePolynomial.createUnsafe(Rings.Z, newData);
    }

    @Override
    public int degree() {return degree;}

    /**
     * Returns i-th coefficient of this poly
     */
    public E get(int i) { return i > degree ? ring.getZero() : data[i];}

    /**
     * Sets i-th coefficient of this poly with specified value
     */
    public UnivariatePolynomial<E> set(int i, E el) {
        el = ring.valueOf(el);
        if (ring.isZero(el)) {
            if (i > degree)
                return this;
            data[i] = el;
            fixDegree();
            return this;
        }
        ensureCapacity(i);
        data[i] = el;
        fixDegree();
        return this;
    }

    /**
     * Sets the leading coefficient of this poly
     */
    public UnivariatePolynomial<E> setLC(E lc) {
        return set(degree, lc);
    }

    @Override
    public int firstNonZeroCoefficientPosition() {
        if (isZero()) return -1;
        int i = 0;
        while (ring.isZero(data[i])) ++i;
        return i;
    }

    /**
     * Returns a copy of this with elements reduced to a new coefficient ring
     *
     * @param newRing the new ring
     * @return a copy of this with elements reduced to a new coefficient ring
     */
    public UnivariatePolynomial<E> setRing(Ring<E> newRing) {
        if (ring == newRing)
            return clone();
        E[] newData = Arrays.copyOf(data, degree + 1);
        newRing.setToValueOf(newData);
        return new UnivariatePolynomial<>(newRing, newData);
    }

    /** internal API */
    public UnivariatePolynomial<E> setRingUnsafe(Ring<E> newRing) {
        return new UnivariatePolynomial<>(newRing, data, degree);
    }

    /**
     * Returns the leading coefficient
     *
     * @return leading coefficient
     */
    public E lc() {return data[degree];}

    @Override
    public UnivariatePolynomial<E> lcAsPoly() {return createConstant(lc());}

    @Override
    public UnivariatePolynomial<E> ccAsPoly() {return createConstant(cc());}

    @Override
    public UnivariatePolynomial<E> getAsPoly(int i) {return createConstant(get(i));}

    /**
     * Returns the constant coefficient
     *
     * @return constant coefficient
     */
    public E cc() {return data[0];}

    @Override
    public void ensureInternalCapacity(int desiredCapacity) {
        if (data.length < desiredCapacity) {
            int oldLength = data.length;
            data = Arrays.copyOf(data, desiredCapacity);
            fillZeroes(data, oldLength, data.length);
        }
    }

    /**
     * Ensures that the capacity of internal storage is enough for storing polynomial of the {@code desiredDegree}. The
     * degree of {@code this} is set to {@code desiredDegree} if the latter is greater than the former.
     *
     * @param desiredDegree desired degree
     */
    final void ensureCapacity(int desiredDegree) {
        if (degree < desiredDegree)
            degree = desiredDegree;

        if (data.length < (desiredDegree + 1)) {
            int oldLen = data.length;
            data = Arrays.copyOf(data, desiredDegree + 1);
            fillZeroes(data, oldLen, data.length);
        }
    }

    /**
     * Removes zeroes from the end of {@code data} and adjusts the degree
     */
    final void fixDegree() {
        int i = degree;
        while (i >= 0 && ring.isZero(data[i])) --i;
        if (i < 0) i = 0;

        if (i != degree) {
            degree = i;
            // not necessary to fillZeroes here!
            // fillZeroes(data, degree + 1, data.length);
        }
    }

    @Override
    public UnivariatePolynomial<E> getRange(int from, int to) {
        return new UnivariatePolynomial<>(ring, Arrays.copyOfRange(data, from, to));
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E>[] createArray(int length) {
        return new UnivariatePolynomial[length];
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E>[] createArray(UnivariatePolynomial<E> a, UnivariatePolynomial<E> b) {
        return new UnivariatePolynomial[]{a, b};
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E>[][] createArray2d(int length) {
        return new UnivariatePolynomial[length][];
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E>[][] createArray2d(int length1, int length2) {
        return new UnivariatePolynomial[length1][length2];
    }

    @Override
    public boolean sameCoefficientRingWith(UnivariatePolynomial<E> oth) {
        return ring.equals(oth.ring);
    }

    @Override
    public UnivariatePolynomial<E> setCoefficientRingFrom(UnivariatePolynomial<E> poly) {
        return setRing(poly.ring);
    }

    /**
     * Creates new poly with the specified coefficients (over the same ring)
     *
     * @param data the data
     * @return polynomial
     */
    public UnivariatePolynomial<E> createFromArray(E[] data) {
        ring.setToValueOf(data);
        return new UnivariatePolynomial<>(ring, data);
    }

    @Override
    public UnivariatePolynomial<E> createMonomial(int degree) {return createMonomial(ring.getOne(), degree);}

    /**
     * Creates linear polynomial of form {@code cc + x * lc} (over the same ring)
     *
     * @param cc the  constant coefficient
     * @param lc the  leading coefficient
     * @return {@code cc + x * lc}
     */
    public UnivariatePolynomial<E> createLinear(E cc, E lc) {
        return createFromArray(ring.createArray(cc, lc));
    }

    /**
     * Creates monomial {@code coefficient * x^degree} (over the same ring)
     *
     * @param coefficient monomial coefficient
     * @param degree      monomial degree
     * @return {@code coefficient * x^degree}
     */
    public UnivariatePolynomial<E> createMonomial(E coefficient, int degree) {
        coefficient = ring.valueOf(coefficient);
        E[] data = ring.createZeroesArray(degree + 1);
        data[degree] = coefficient;
        return new UnivariatePolynomial<>(ring, data);
    }

    /**
     * Creates constant polynomial with specified value (over the same ring)
     *
     * @param val the value
     * @return constant polynomial with specified value
     */
    public UnivariatePolynomial<E> createConstant(E val) {
        E[] array = ring.createArray(1);
        array[0] = val;
        return createFromArray(array);
    }

    @Override
    public UnivariatePolynomial<E> createZero() {return createConstant(ring.getZero());}

    @Override
    public UnivariatePolynomial<E> createOne() {return createConstant(ring.getOne());}

    @Override
    public boolean isZeroAt(int i) {return i >= data.length || ring.isZero(data[i]);}

    @Override
    public final UnivariatePolynomial<E> setZero(int i) {
        if (i < data.length)
            data[i] = ring.getZero();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> setFrom(int indexInThis, UnivariatePolynomial<E> poly, int indexInPoly) {
        ensureCapacity(indexInThis);
        data[indexInThis] = poly.get(indexInPoly);
        fixDegree();
        return this;
    }

    @Override
    public boolean isZero() {return ring.isZero(data[degree]);}

    @Override
    public boolean isOne() {return degree == 0 && ring.isOne(data[0]);}

    @Override
    public boolean isMonic() {return ring.isOne(lc());}

    @Override
    public boolean isUnitCC() {return ring.isOne(cc());}

    @Override
    public boolean isConstant() {return degree == 0;}

    @Override
    public boolean isMonomial() {
        for (int i = degree - 1; i >= 0; --i)
            if (!ring.isZero(data[i]))
                return false;
        return true;
    }

    @Override
    public int signumOfLC() {
        return ring.signum(lc());
    }

    @Override
    public boolean isOverField() {
        return ring.isField();
    }

    @Override
    public boolean isOverFiniteField() {
        return ring.isFinite();
    }

    @Override
    public boolean isOverZ() {return ring.equals(Rings.Z);}

    @Override
    public BigInteger coefficientRingCardinality() {
        return ring.cardinality();
    }

    @Override
    public BigInteger coefficientRingCharacteristic() {
        return ring.characteristic();
    }

    @Override
    public boolean isOverPerfectPower() {
        return ring.isPerfectPower();
    }

    @Override
    public BigInteger coefficientRingPerfectPowerBase() {
        return ring.perfectPowerBase();
    }

    @Override
    public BigInteger coefficientRingPerfectPowerExponent() {
        return ring.perfectPowerExponent();
    }

    /**
     * Returns Mignotte's bound (sqrt(n+1) * 2^n max |this|) of the poly
     */
    public static BigInteger mignotteBound(UnivariatePolynomial<BigInteger> poly) {
        return ONE.shiftLeft(poly.degree).multiply(norm2(poly));
    }

    /**
     * Returns L1 norm of the polynomial, i.e. sum of abs coefficients
     */
    public static BigInteger norm1(UnivariatePolynomial<BigInteger> poly) {
        BigInteger norm = ZERO;
        for (int i = poly.degree; i >= 0; --i)
            norm = norm.add(abs(poly.data[i]));
        return norm;
    }

    /**
     * Returns L2 norm of the polynomial, i.e. a square root of a sum of coefficient squares.
     */
    public static BigInteger norm2(UnivariatePolynomial<BigInteger> poly) {
        BigInteger norm = ZERO;
        for (int i = poly.degree; i >= 0; --i)
            norm = norm.add(poly.data[i].multiply(poly.data[i]));
        return BigIntegerUtil.sqrtCeil(norm);
    }

    /**
     * Returns L2 norm of the poly, i.e. a square root of a sum of coefficient squares.
     */
    public static double norm2Double(UnivariatePolynomial<BigInteger> poly) {
        double norm = 0;
        for (int i = poly.degree; i >= 0; --i) {
            double d = poly.data[i].doubleValue();
            norm += d * d;
        }
        return Math.sqrt(norm);
    }

    /**
     * Returns max abs coefficient of the poly
     */
    public E maxAbsCoefficient() {
        E el = ring.abs(data[0]);
        for (int i = 1; i <= degree; i++)
            el = ring.max(el, ring.abs(data[i]));
        return el;
    }

    /**
     */
    public E normMax() {
        return maxAbsCoefficient();
    }

    private void fillZeroes(E[] data, int from, int to) {
        for (int i = from; i < to; ++i)
            data[i] = ring.getZero(); //invoke getZero() at each cycle
    }

    @Override
    public UnivariatePolynomial<E> toZero() {
        fillZeroes(data, 0, degree + 1);
        degree = 0;
        return this;
    }

    @Override
    public UnivariatePolynomial<E> set(UnivariatePolynomial<E> oth) {
        if (oth == this)
            return this;
        this.data = oth.data.clone();
        this.degree = oth.degree;
        return this;
    }

    @Override
    public final UnivariatePolynomial<E> setAndDestroy(UnivariatePolynomial<E> oth) {
        this.data = oth.data;
        oth.data = null; // destroy
        this.degree = oth.degree;
        assert data.length > 0;
        return this;
    }

    @Override
    public UnivariatePolynomial<E> shiftLeft(int offset) {
        if (offset == 0)
            return this;
        if (offset > degree)
            return toZero();

        System.arraycopy(data, offset, data, 0, degree - offset + 1);
        fillZeroes(data, degree - offset + 1, degree + 1);
        degree = degree - offset;
        return this;
    }

    @Override
    public UnivariatePolynomial<E> shiftRight(int offset) {
        if (offset == 0)
            return this;
        int degree = this.degree;
        ensureCapacity(offset + degree);
        System.arraycopy(data, 0, data, offset, degree + 1);
        fillZeroes(data, 0, offset);
        return this;
    }

    @Override
    public UnivariatePolynomial<E> truncate(int newDegree) {
        if (newDegree >= degree)
            return this;
        fillZeroes(data, newDegree + 1, degree + 1);
        degree = newDegree;
        fixDegree();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> reverse() {
        ArraysUtil.reverse(data, 0, degree + 1);
        fixDegree();
        return this;
    }

    /**
     * Returns the content of the poly
     *
     * @return polynomial content
     */
    public E content() {
        if (degree == 0)
            return data[0];
        return ring.gcd(this);
//        E gcd = data[degree];
//        for (int i = degree - 1; i >= 0; --i)
//            gcd = ring.gcd(gcd, data[i]);
//        return gcd;
    }

    @Override
    public UnivariatePolynomial<E> contentAsPoly() {
        return createConstant(content());
    }

    @Override
    public UnivariatePolynomial<E> primitivePart() {
        E content = content();
        if (signumOfLC() < 0 && ring.signum(content) > 0)
            content = ring.negate(content);
        if (ring.isMinusOne(content))
            return negate();
        return primitivePart0(content);
    }

    @Override
    public UnivariatePolynomial<E> primitivePartSameSign() {
        return primitivePart0(content());
    }

    private UnivariatePolynomial<E> primitivePart0(E content) {
        if (isZero())
            return this;
        if (ring.isOne(content))
            return this;
        for (int i = degree; i >= 0; --i) {
            data[i] = ring.divideOrNull(data[i], content);
            if (data[i] == null)
                return null;
        }
        return this;
    }

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    public E evaluate(long point) {
        return evaluate(ring.valueOf(point));
    }

    /**
     * Evaluates this poly at a given {@code point} (via Horner method).
     *
     * @param point {@code point}
     * @return value at {@code point}
     */
    public E evaluate(E point) {
        if (ring.isZero(point))
            return cc();

        point = ring.valueOf(point);
        E res = ring.getZero();
        for (int i = degree; i >= 0; --i)
            res = ring.addMutable(ring.multiplyMutable(res, point), data[i]);
        return res;
    }

    @Override
    public UnivariatePolynomial<E> composition(UnivariatePolynomial<E> value) {
        if (value.isOne())
            return this.clone();
        if (value.isZero())
            return ccAsPoly();

        UnivariatePolynomial<E> result = createZero();
        for (int i = degree; i >= 0; --i)
            result = result.multiply(value).add(data[i]);
        return result;
    }

    /**
     * Shifts variable x -> x + value and returns the result (new instance)
     *
     * @param value shift amount
     * @return a copy of this with x -> x + value
     */
    public UnivariatePolynomial<E> shift(E value) {
        return composition(createLinear(value, ring.getOne()));
    }

    /**
     * Add constant to this.
     *
     * @param val some number
     * @return this + val
     */
    public UnivariatePolynomial<E> add(E val) {
        data[0] = ring.add(data[0], ring.valueOf(val));
        fixDegree();
        return this;
    }

    /**
     * Subtract constant from this.
     *
     * @param val some number
     * @return this + val
     */
    public UnivariatePolynomial<E> subtract(E val) {
        data[0] = ring.subtract(data[0], ring.valueOf(val));
        fixDegree();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> decrement() {
        return subtract(createOne());
    }

    @Override
    public UnivariatePolynomial<E> increment() {
        return add(createOne());
    }

    @Override
    public UnivariatePolynomial<E> add(UnivariatePolynomial<E> oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth);

        assertSameCoefficientRingWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = ring.add(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /**
     * Adds {@code coefficient*x^exponent} to {@code this}
     *
     * @param coefficient monomial coefficient
     * @param exponent    monomial exponent
     * @return {@code this + coefficient*x^exponent}
     */
    public UnivariatePolynomial<E> addMonomial(E coefficient, int exponent) {
        if (ring.isZero(coefficient))
            return this;

        ensureCapacity(exponent);
        data[exponent] = ring.add(data[exponent], ring.valueOf(coefficient));
        fixDegree();
        return this;
    }

    /**
     * Adds {@code oth * factor} to {@code this}
     *
     * @param oth    the polynomial
     * @param factor the factor
     * @return {@code this + oth * factor} modulo {@code modulus}
     */
    public UnivariatePolynomial<E> addMul(UnivariatePolynomial<E> oth, E factor) {
        if (oth.isZero())
            return this;

        factor = ring.valueOf(factor);
        if (ring.isZero(factor))
            return this;

        assertSameCoefficientRingWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = ring.add(data[i], ring.multiply(factor, oth.data[i]));
        fixDegree();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> subtract(UnivariatePolynomial<E> oth) {
        if (oth.isZero())
            return this;
        if (isZero())
            return set(oth).negate();

        assertSameCoefficientRingWith(oth);
        ensureCapacity(oth.degree);
        for (int i = oth.degree; i >= 0; --i)
            data[i] = ring.subtract(data[i], oth.data[i]);
        fixDegree();
        return this;
    }

    /**
     * Subtracts {@code factor * x^exponent * oth} from {@code this}
     *
     * @param oth      the polynomial
     * @param factor   the factor
     * @param exponent the exponent
     * @return {@code this - factor * x^exponent * oth}
     */
    public UnivariatePolynomial<E> subtract(UnivariatePolynomial<E> oth, E factor, int exponent) {
        if (oth.isZero())
            return this;

        factor = ring.valueOf(factor);
        if (ring.isZero(factor))
            return this;

        assertSameCoefficientRingWith(oth);
        for (int i = oth.degree + exponent; i >= exponent; --i)
            data[i] = ring.subtract(data[i], ring.multiply(factor, oth.data[i - exponent]));

        fixDegree();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> negate() {
        for (int i = degree; i >= 0; --i)
            if (!ring.isZero(data[i]))
                data[i] = ring.negate(data[i]);
        return this;
    }

    /**
     * Multiplies {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public UnivariatePolynomial<E> multiply(E factor) {
        factor = ring.valueOf(factor);
        if (ring.isOne(factor))
            return this;

        if (ring.isZero(factor))
            return toZero();

        for (int i = degree; i >= 0; --i)
            data[i] = ring.multiply(data[i], factor);
        return this;
    }

    @Override
    public UnivariatePolynomial<E> multiplyByLC(UnivariatePolynomial<E> other) {
        return multiply(other.lc());
    }

    @Override
    public UnivariatePolynomial<E> multiply(long factor) {
        return multiply(ring.valueOf(factor));
    }

    @Override
    public UnivariatePolynomial<E> divideByLC(UnivariatePolynomial<E> other) {
        return divideOrNull(other.lc());
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of
     * the elements can't be exactly divided by the {@code factor}. NOTE: if {@code null} is returned, the content of
     * {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public UnivariatePolynomial<E> divideOrNull(E factor) {
        factor = ring.valueOf(factor);
        if (ring.isZero(factor))
            throw new ArithmeticException("Divide by zero");
        if (ring.isOne(factor))
            return this;
        if (ring.isMinusOne(factor))
            return negate();
        if (ring.isField()) // this is typically much faster
            return multiply(ring.reciprocal(factor));

        for (int i = degree; i >= 0; --i) {
            E l = ring.divideOrNull(data[i], factor);
            if (l == null)
                return null;
            data[i] = l;
        }
        return this;
    }

    /**
     * Divides this polynomial by a {@code factor} or throws exception if exact division is not possible
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor}
     * @throws ArithmeticException if exact division is not possible
     */
    public UnivariatePolynomial<E> divideExact(E factor) {
        UnivariatePolynomial<E> r = divideOrNull(factor);
        if (r == null)
            throw new ArithmeticException("not divisible " + this + " / " + factor);
        return r;
    }

    @Override
    public UnivariatePolynomial<E> monic() {
        if (isZero())
            return this;
        return divideOrNull(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor}.
     *
     * @param factor the factor
     * @return {@code this}
     */
    public UnivariatePolynomial<E> monic(E factor) {
        E lc = lc();
        return multiply(factor).divideOrNull(lc);
    }

    @Override
    public UnivariatePolynomial<E> monicWithLC(UnivariatePolynomial<E> other) {
        if (lc().equals(other.lc()))
            return this;
        return monic(other.lc());
    }

    @Override
    public UnivariatePolynomial<E> multiplyByBigInteger(BigInteger factor) {
        return multiply(ring.valueOfBigInteger(factor));
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E> multiply(UnivariatePolynomial<E> oth) {
        if (isZero())
            return this;
        if (oth.isZero())
            return toZero();
        if (this == oth)
            return square();

        if (oth.degree == 0)
            return multiply(oth.data[0]);
        if (degree == 0) {
            E factor = data[0];
            data = oth.data.clone();
            degree = oth.degree;
            return multiply(factor);
        }

        if (ring instanceof IntegersZp) {
            // faster method with exact operations
            UnivariatePolynomial<E>
                    iThis = setRingUnsafe((Ring<E>) Rings.Z),
                    iOth = oth.setRingUnsafe((Ring<E>) Rings.Z);
            data = iThis.multiply(iOth).data;
            ring.setToValueOf(data);
        } else
            data = multiplySafe0(oth);
        degree += oth.degree;
        fixDegree();
        return this;
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E> square() {
        if (isZero())
            return this;
        if (degree == 0)
            return multiply(data[0]);

        if (ring instanceof IntegersZp) {
            // faster method with exact operations
            UnivariatePolynomial<E> iThis = setRingUnsafe((Ring<E>) Rings.Z);
            data = iThis.square().data;
            ring.setToValueOf(data);
        } else
            data = squareSafe0();
        degree += degree;
        fixDegree();
        return this;
    }

    @Override
    public UnivariatePolynomial<E> derivative() {
        if (isConstant())
            return createZero();
        E[] newData = ring.createArray(degree);
        for (int i = degree; i > 0; --i)
            newData[i - 1] = ring.multiply(data[i], ring.valueOf(i));
        return createFromArray(newData);
    }

    @Override
    public UnivariatePolynomial<E> clone() {
        return new UnivariatePolynomial<>(ring, data.clone(), degree);
    }

    @Override
    public UnivariatePolynomial<E> parsePoly(String string) {
        return parse(string, ring);
    }

    public UnivariatePolynomial<E> parsePoly(String string, ElementParser<E> eParser, String variable) {
        return Parser.parse(ring, eParser, string);
    }

    @Override
    public Iterator<E> iterator() {
        return new It();
    }

    /**
     * Returns a sequential {@code Stream} with coefficients of this as its source.
     *
     * @return a sequential {@code Stream} over the coefficients in this polynomial
     */
    public Stream<E> stream() {
        return StreamSupport.stream(spliterator(), false);
    }

    @Override
    public Spliterator<E> spliterator() {
        return Arrays.spliterator(data, 0, degree + 1);
    }

    /**
     * Applies transformation function to this and returns the result. This method is equivalent of {@code
     * stream().map(mapper).collect(new PolynomialCollector<>(ring))}.
     *
     * @param ring   ring of the new polynomial
     * @param mapper function that maps coefficients of this to coefficients of the result
     * @param <T>    result elements type
     * @return a new polynomial with the coefficients obtained from this by applying {@code mapper}
     */
    public <T> UnivariatePolynomial<T> mapCoefficients(Ring<T> ring, Function<E, T> mapper) {
        return stream().map(mapper).collect(new PolynomialCollector<>(ring));
    }

    /**
     * Applies transformation function to this and returns the result. This method is equivalent of {@code
     * stream().map(mapper).collect(new PolynomialCollector<>(ring))}.
     *
     * @param ring   ring of the new polynomial
     * @param mapper function that maps coefficients of this to coefficients of the result
     * @return a new polynomial with the coefficients obtained from this by applying {@code mapper}
     */
    public UnivariatePolynomialZp64 mapCoefficients(IntegersZp64 ring, ToLongFunction<E> mapper) {
        return UnivariatePolynomialZp64.create(ring, stream().mapToLong(mapper).toArray());
    }

    private static final class ListToPoly<E> implements Function<List<E>, UnivariatePolynomial<E>> {
        final Ring<E> ring;

        public ListToPoly(Ring<E> ring) {
            this.ring = ring;
        }

        @Override
        public UnivariatePolynomial<E> apply(List<E> es) {
            return UnivariatePolynomial.create(ring, es.toArray(ring.createArray(es.size())));
        }
    }

    /**
     * Collector which collects stream of element to a UnivariatePolynomial
     *
     * @param <E> element type
     */
    public static final class PolynomialCollector<E>
            implements Collector<E, List<E>, UnivariatePolynomial<E>> {
        final Supplier<List<E>> supplier = ArrayList::new;
        final BiConsumer<List<E>, E> accumulator = List::add;
        final BinaryOperator<List<E>> combiner = (l, r) -> {
            l.addAll(r);
            return l;
        };
        final Function<List<E>, UnivariatePolynomial<E>> finisher;
        final Ring<E> ring;

        public PolynomialCollector(Ring<E> ring) {
            this.ring = ring;
            this.finisher = new ListToPoly<>(ring);
        }

        @Override
        public Supplier<List<E>> supplier() {
            return supplier;
        }

        @Override
        public BiConsumer<List<E>, E> accumulator() {
            return accumulator;
        }

        @Override
        public BinaryOperator<List<E>> combiner() {
            return combiner;
        }

        @Override
        public Function<List<E>, UnivariatePolynomial<E>> finisher() {
            return finisher;
        }

        @Override
        public Set<Characteristics> characteristics() {
            return Collections.emptySet();
        }
    }

    private final class It implements Iterator<E> {
        int index = 0;

        @Override
        public boolean hasNext() {
            synchronized (UnivariatePolynomial.this) {
                return index <= degree;
            }
        }

        @Override
        public E next() {
            synchronized (UnivariatePolynomial.this) {
                return data[index++];
            }
        }
    }

    /** internal API */
    public E[] getDataReferenceUnsafe() {return data;}

    @Override
    public MultivariatePolynomial<E> asMultivariate() {
        return MultivariatePolynomial.asMultivariate(this, 1, 0, MonomialOrder.DEFAULT);
    }

    @Override
    public int compareTo(UnivariatePolynomial<E> o) {
        int c = Integer.compare(degree, o.degree);
        if (c != 0)
            return c;
        for (int i = degree; i >= 0; --i) {
            c = ring.compare(data[i], o.data[i]);
            if (c != 0)
                return c;
        }
        return 0;
    }

    @Override
    public String coefficientRingToString() {
        return ring.toString();
    }

    @Override
    public String toString() {
        return toString(WithVariables.defaultVars(1));
    }

    @Override
    public String toString(String[] variables) {
        return toString(ring, variables[0]);
    }

    public String toString(ToStringSupport<E> cfxToString, String variable) {
        if (isZero())
            return "0";
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < data.length; i++) {
            String str = termToString(cfxToString, i, variable);
            if (sb.length() == 0 && str.startsWith("+"))
                str = str.substring(1);
            sb.append(str);
        }
        return sb.toString();
    }

    private String termToString(ToStringSupport<E> cfxToString, int i, String var) {
        E el = data[i];
        if (ring.isZero(el))
            return "";
        String coefficient;
        boolean needSeparator;
        if (ring.isOne(el)) {
            coefficient = "";
            needSeparator = false;
        } else if (!(ring instanceof IntegersZp) && ring.isMinusOne(el)) {
            coefficient = "-";
            needSeparator = false;
        } else {
            coefficient = cfxToString.toString(el);
            if (!coefficient.matches("[+\\-]?\\d+"))
                coefficient = "(" + coefficient + ")";
            needSeparator = true;
        }

        if (!coefficient.startsWith("-") && !coefficient.startsWith("+"))
            coefficient = "+" + coefficient;

        String m;
        if (i == 0)
            m = "";
        else
            m = ((needSeparator ? "*" : "") + var + (i == 1 ? "" : "^" + i));
        if (m.isEmpty())
            if (ring.isOne(el) || (!(ring instanceof IntegersZp) && ring.isMinusOne(el)))
                coefficient = coefficient + "1";

        return coefficient + m;
    }

    String toStringForCopy() {
        String s = ArraysUtil.toString(data, 0, degree + 1, x -> "new BigInteger(\"" + x + "\")");
        return "of(" + s.substring(1, s.length() - 1) + ")";
    }

    @Override
    public boolean equals(Object obj) {
        if (obj.getClass() != this.getClass())
            return false;
        @SuppressWarnings("unchecked")
        UnivariatePolynomial<E> oth = (UnivariatePolynomial<E>) obj;
        if (degree != oth.degree)
            return false;
        for (int i = 0; i <= degree; ++i)
            if (!(data[i].equals(oth.data[i])))
                return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = 1;
        for (int i = degree; i >= 0; --i)
            result = 31 * result + data[i].hashCode();
        return result;
    }

    /* =========================== Exact multiplication with safe arithmetics =========================== */


    /** switch to classical multiplication */
    static final long KARATSUBA_THRESHOLD = 1024L;
    /** when use Karatsuba fast multiplication */
    static final long
            MUL_CLASSICAL_THRESHOLD = 256L * 256L,
            MUL_MOD_CLASSICAL_THRESHOLD = 128L * 128L;

    /** switch algorithms */
    final E[] multiplySafe0(UnivariatePolynomial<E> oth) {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return multiplyClassicalSafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
        else
            return multiplyKaratsubaSafe(data, 0, degree + 1, oth.data, 0, oth.degree + 1);
    }

    /** switch algorithms */
    final E[] squareSafe0() {
        if (1L * (degree + 1) * (degree + 1) <= MUL_CLASSICAL_THRESHOLD)
            return squareClassicalSafe(data, 0, degree + 1);
        else
            return squareKaratsubaSafe(data, 0, degree + 1);
    }

    /**
     * Classical n*m multiplication algorithm
     *
     * @param a     the first multiplier
     * @param aFrom begin in a
     * @param aTo   end in a
     * @param b     the second multiplier
     * @param bFrom begin in b
     * @param bTo   end in b
     * @return the result
     */
    final E[] multiplyClassicalSafe(final E[] a, final int aFrom, final int aTo, final E[] b, final int bFrom, final int bTo) {
        E[] result = ring.createZeroesArray(aTo - aFrom + bTo - bFrom - 1);
        multiplyClassicalSafe(result, a, aFrom, aTo, b, bFrom, bTo);
        return result;
    }

    /**
     * Classical n*m multiplication algorithm
     *
     * @param result where to write the result
     * @param a      the first multiplier
     * @param aFrom  begin in a
     * @param aTo    end in a
     * @param b      the second multiplier
     * @param bFrom  begin in b
     * @param bTo    end in b
     */
    final void multiplyClassicalSafe(final E[] result, final E[] a, final int aFrom, final int aTo, final E[] b, final int bFrom, final int bTo) {
        if (aTo - aFrom > bTo - bFrom) {
            multiplyClassicalSafe(result, b, bFrom, bTo, a, aFrom, aTo);
            return;
        }
        for (int i = 0; i < aTo - aFrom; ++i) {
            E c = a[aFrom + i];
            if (!ring.isZero(c))
                for (int j = 0; j < bTo - bFrom; ++j)
                    result[i + j] = ring.addMutable(result[i + j], ring.multiply(c, b[bFrom + j]));
        }
    }

    /**
     * Karatsuba multiplication
     *
     * @param f     the first multiplier
     * @param g     the second multiplier
     * @param fFrom begin in f
     * @param fTo   end in f
     * @param gFrom begin in g
     * @param gTo   end in g
     * @return the result
     */
    E[] multiplyKaratsubaSafe(
            final E[] f, final int fFrom, final int fTo,
            final E[] g, final int gFrom, final int gTo) {
        // return zero
        if (fFrom >= fTo || gFrom >= gTo)
            return ring.createArray(0);

        // single element in f
        if (fTo - fFrom == 1) {
            E[] result = ring.createArray(gTo - gFrom);
            for (int i = gFrom; i < gTo; ++i)
                result[i - gFrom] = ring.multiply(f[fFrom], g[i]);
            return result;
        }
        // single element in g
        if (gTo - gFrom == 1) {
            E[] result = ring.createArray(fTo - fFrom);
            //single element in b
            for (int i = fFrom; i < fTo; ++i)
                result[i - fFrom] = ring.multiply(g[gFrom], f[i]);
            return result;
        }
        // linear factors
        if (fTo - fFrom == 2 && gTo - gFrom == 2) {
            E[] result = ring.createArray(3);
            //both a and b are linear
            result[0] = ring.multiply(f[fFrom], g[gFrom]);
            result[1] = ring.addMutable(ring.multiply(f[fFrom], g[gFrom + 1]), ring.multiply(f[fFrom + 1], g[gFrom]));
            result[2] = ring.multiply(f[fFrom + 1], g[gFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (gTo - gFrom) < KARATSUBA_THRESHOLD)
            return multiplyClassicalSafe(g, gFrom, gTo, f, fFrom, fTo);

        if (fTo - fFrom < gTo - gFrom)
            return multiplyKaratsubaSafe(g, gFrom, gTo, f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        //if we can't split b
        if (gFrom + split >= gTo) {
            E[] f0g = multiplyKaratsubaSafe(f, fFrom, fFrom + split, g, gFrom, gTo);
            E[] f1g = multiplyKaratsubaSafe(f, fFrom + split, fTo, g, gFrom, gTo);

            int oldLen = f0g.length, newLen = fTo - fFrom + gTo - gFrom - 1;
            E[] result = Arrays.copyOf(f0g, newLen);
            fillZeroes(result, oldLen, newLen);
            for (int i = 0; i < f1g.length; i++)
                result[i + split] = ring.addMutable(result[i + split], f1g[i]);
            return result;
        }

        int fMid = fFrom + split, gMid = gFrom + split;
        E[] f0g0 = multiplyKaratsubaSafe(f, fFrom, fMid, g, gFrom, gMid);
        E[] f1g1 = multiplyKaratsubaSafe(f, fMid, fTo, g, gMid, gTo);

        // f0 + f1
        E[] f0_plus_f1 = ring.createArray(Math.max(fMid - fFrom, fTo - fMid));
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        fillZeroes(f0_plus_f1, fMid - fFrom, f0_plus_f1.length);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = ring.add(f0_plus_f1[i - fMid], f[i]);

        //g0 + g1
        E[] g0_plus_g1 = ring.createArray(Math.max(gMid - gFrom, gTo - gMid));
        System.arraycopy(g, gFrom, g0_plus_g1, 0, gMid - gFrom);
        fillZeroes(g0_plus_g1, gMid - gFrom, g0_plus_g1.length);
        for (int i = gMid; i < gTo; ++i)
            g0_plus_g1[i - gMid] = ring.add(g0_plus_g1[i - gMid], g[i]);

        E[] mid = multiplyKaratsubaSafe(f0_plus_f1, 0, f0_plus_f1.length, g0_plus_g1, 0, g0_plus_g1.length);

        if (mid.length < f0g0.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f0g0.length);
            fillZeroes(mid, oldLen, mid.length);
        }
        if (mid.length < f1g1.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f1g1.length);
            fillZeroes(mid, oldLen, mid.length);
        }

        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = ring.subtractMutable(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = ring.subtractMutable(mid[i], f1g1[i]);


        int oldLen = f0g0.length;
        E[] result = Arrays.copyOf(f0g0, (fTo - fFrom) + (gTo - gFrom) - 1);
        fillZeroes(result, oldLen, result.length);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = ring.addMutable(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = ring.addMutable(result[i + 2 * split], f1g1[i]);

        return result;
    }

    E[] squareClassicalSafe(E[] a, int from, int to) {
        E[] x = ring.createZeroesArray((to - from) * 2 - 1);
        squareClassicalSafe(x, a, from, to);
        return x;
    }


    /**
     * Square the poly {@code data} using classical algorithm
     *
     * @param result result destination
     * @param data   the data
     * @param from   data from
     * @param to     end point in the {@code data}
     */
    void squareClassicalSafe(final E[] result, E[] data, int from, int to) {
        int len = to - from;
        for (int i = 0; i < len; ++i) {
            E c = data[from + i];
            if (!ring.isZero(c))
                for (int j = 0; j < len; ++j)
                    result[i + j] = ring.addMutable(result[i + j], ring.multiply(c, data[from + j]));
        }
    }

    /**
     * Karatsuba squaring
     *
     * @param f     the data
     * @param fFrom begin in f
     * @param fTo   end in f
     * @return the result
     */
    E[] squareKaratsubaSafe(final E[] f, final int fFrom, final int fTo) {
        if (fFrom >= fTo)
            return ring.createArray(0);
        if (fTo - fFrom == 1) {
            E[] r = ring.createArray(1);
            r[0] = ring.multiply(f[fFrom], f[fFrom]);
            return r;
        }
        if (fTo - fFrom == 2) {
            E[] result = ring.createArray(3);
            result[0] = ring.multiply(f[fFrom], f[fFrom]);
            result[1] = ring.multiplyMutable(ring.multiply(f[fFrom], f[fFrom + 1]), ring.valueOf(2));
            result[2] = ring.multiply(f[fFrom + 1], f[fFrom + 1]);
            return result;
        }
        //switch to classical
        if (1L * (fTo - fFrom) * (fTo - fFrom) < KARATSUBA_THRESHOLD)
            return squareClassicalSafe(f, fFrom, fTo);


        //we now split a and b into 2 parts:
        int split = (fTo - fFrom + 1) / 2;
        int fMid = fFrom + split;
        E[] f0g0 = squareKaratsubaSafe(f, fFrom, fMid);
        E[] f1g1 = squareKaratsubaSafe(f, fMid, fTo);

        // f0 + f1
        E[] f0_plus_f1 = ring.createArray(Math.max(fMid - fFrom, fTo - fMid));
        System.arraycopy(f, fFrom, f0_plus_f1, 0, fMid - fFrom);
        fillZeroes(f0_plus_f1, fMid - fFrom, f0_plus_f1.length);
        for (int i = fMid; i < fTo; ++i)
            f0_plus_f1[i - fMid] = ring.add(f0_plus_f1[i - fMid], f[i]);

        E[] mid = squareKaratsubaSafe(f0_plus_f1, 0, f0_plus_f1.length);

        if (mid.length < f0g0.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f0g0.length);
            fillZeroes(mid, oldLen, mid.length);
        }
        if (mid.length < f1g1.length) {
            int oldLen = mid.length;
            mid = Arrays.copyOf(mid, f1g1.length);
            fillZeroes(mid, oldLen, mid.length);
        }


        //subtract f0g0, f1g1
        for (int i = 0; i < f0g0.length; ++i)
            mid[i] = ring.subtractMutable(mid[i], f0g0[i]);
        for (int i = 0; i < f1g1.length; ++i)
            mid[i] = ring.subtractMutable(mid[i], f1g1[i]);


        int oldLen = f0g0.length;
        E[] result = Arrays.copyOf(f0g0, 2 * (fTo - fFrom) - 1);
        fillZeroes(result, oldLen, result.length);
        for (int i = 0; i < mid.length; ++i)
            result[i + split] = ring.addMutable(result[i + split], mid[i]);
        for (int i = 0; i < f1g1.length; ++i)
            result[i + 2 * split] = ring.addMutable(result[i + 2 * split], f1g1[i]);

        return result;
    }
}
