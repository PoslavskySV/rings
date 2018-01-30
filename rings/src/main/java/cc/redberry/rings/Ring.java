package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;

/**
 * Ring of elements. Mathematical operations defined in {@code Ring} interface include all <i>field</i> operations,
 * though the particular implementations may represent a more restricted rings (general rings, Euclidean rings etc.), in
 * which case some field operations (e.g. reciprocal) are not applicable (will throw exception).
 *
 * @param <E> the type of objects that may be operated by this ring
 * @since 1.0
 */
public interface Ring<E> extends
                         Comparator<E>,
                         Iterable<E>,
                         ToStringSupport<E>,
                         ElementParser<E>,
                         java.io.Serializable {
    /**
     * Returns whether this ring is a field
     *
     * @return whether this ring is a field
     */
    boolean isField();

    /**
     * Returns whether this ring is a Euclidean ring
     *
     * @return whether this ring is a Euclidean ring
     */
    boolean isEuclideanRing();

    /**
     * Returns whether this ring is finite
     *
     * @return whether this ring is finite
     */
    default boolean isFinite() {
        return cardinality() != null;
    }

    /**
     * Returns whether this ring is a finite field
     *
     * @return whether this ring is a finite field
     */
    default boolean isFiniteField() {
        return isField() && isFinite();
    }

    /**
     * Returns the number of elements in this ring (cardinality) or null if ring is infinite
     *
     * @return the number of elements in this ring (cardinality) or null if ring is infinite
     */
    BigInteger cardinality();

    /**
     * Returns characteristic of this ring
     *
     * @return characteristic of this ring
     */
    BigInteger characteristic();

    /**
     * Returns whether the cardinality is a perfect power (p^k with k > 1)
     *
     * @return whether the cardinality is a perfect power (p^k with k > 1)
     */
    boolean isPerfectPower();

    /**
     * Returns {@code base} so that {@code cardinality == base^exponent} or null if cardinality is not finite
     *
     * @return {@code base} so that {@code cardinality == base^exponent} or null if cardinality is not finite
     */
    BigInteger perfectPowerBase();

    /**
     * Returns {@code exponent} so that {@code cardinality == base^exponent} or null if cardinality is not finite
     *
     * @return {@code exponent} so that {@code cardinality == base^exponent} or null if cardinality is not finite
     */
    BigInteger perfectPowerExponent();

    /**
     * Add two elements
     *
     * @param a the first element
     * @param b the second element
     * @return a + b
     */
    E add(E a, E b);

    /**
     * Total of the array of elements
     *
     * @param elements elements to sum
     * @return sum of the array
     */
    @SuppressWarnings("unchecked")
    default E add(E... elements) {
        E r = elements[0];
        for (int i = 1; i < elements.length; i++)
            r = add(r, elements[i]);
        return r;
    }

    /**
     * Returns {@code element + 1}
     *
     * @param element the element
     * @return {@code element + 1}
     */
    default E increment(E element) {
        return add(element, getOne());
    }

    /**
     * Returns {@code element - 1}
     *
     * @param element the element
     * @return {@code element - 1}
     */
    default E decrement(E element) {
        return subtract(element, getOne());
    }

    /**
     * Subtracts {@code b} from {@code a}
     *
     * @param a the first element
     * @param b the second element
     * @return a - b
     */
    E subtract(E a, E b);

    /**
     * Multiplies two elements
     *
     * @param a the first element
     * @param b the second element
     * @return a * b
     */
    E multiply(E a, E b);

    /**
     * Multiplies the array of elements
     *
     * @param elements the elements
     * @return product of the array
     */
    @SuppressWarnings("unchecked")
    default E multiply(E... elements) {
        E r = elements[0];
        for (int i = 1; i < elements.length; i++)
            r = multiply(r, elements[i]);
        return r;
    }

    /**
     * Negates the given element
     *
     * @param element the ring element
     * @return -val
     */
    E negate(E element);

    /**
     * Adds two elements and destroys the initial content of {@code a}.
     *
     * @param a the first element (may be destroyed)
     * @param b the second element
     * @return a + b
     */
    default E addMutable(E a, E b) {return add(a, b);}

    /**
     * Subtracts {@code b} from {@code a} and destroys the initial content of {@code a}
     *
     * @param a the first element (may be destroyed)
     * @param b the second element
     * @return a - b
     */
    default E subtractMutable(E a, E b) {return subtract(a, b);}

    /**
     * Multiplies two elements and destroys the initial content of {@code a}
     *
     * @param a the first element (may be destroyed)
     * @param b the second element
     * @return a * b
     */
    default E multiplyMutable(E a, E b) { return multiply(a, b);}

    /**
     * Negates the given element and destroys the initial content of {@code element}
     *
     * @param element the ring element (may be destroyed)
     * @return -element
     */
    default E negateMutable(E element) { return negate(element);}

    /**
     * Makes a deep copy of the specified element (for immutable instances the same reference returned).
     *
     * @param element the element
     * @return deep copy of specified element
     */
    E copy(E element);

    /**
     * Returns -1 if {@code element < 0}, 0 if {@code element == 0} and 1 if {@code element > 0}, where comparison is
     * specified by {@link #compare(Object, Object)}
     *
     * @param element the element
     * @return -1 if {@code element < 0}, 0 if {@code element == 0} and 1 otherwise
     */
    default int signum(E element) {
        return Integer.compare(compare(element, getZero()), 0);
    }

    /**
     * Returns the abs value of element (no copy)
     */
    default E abs(E el) {
        return signum(el) < 0 ? negate(el) : el;
    }

    /**
     * Returns the max value (no copy)
     */
    default E max(E a, E b) {
        return compare(a, b) < 0 ? b : a;
    }

    /**
     * Returns the min value (no copy)
     */
    default E min(E a, E b) {
        return compare(a, b) > 0 ? b : a;
    }

    /**
     * Returns quotient and remainder of {@code dividend / divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code {quotient, remainder}}
     */
    E[] divideAndRemainder(E dividend, E divider);

    /**
     * Returns the quotient of {@code dividend / divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return the quotient of {@code dividend / divider}
     */
    default E quotient(E dividend, E divider) {
        return divideAndRemainder(dividend, divider)[0];
    }

    /**
     * Returns the remainder of {@code dividend / divider}
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return the remainder of {@code dividend / divider}
     */
    default E remainder(E dividend, E divider) {
        return divideAndRemainder(dividend, divider)[1];
    }

    /**
     * Divides {@code dividend} by {@code divider} or returns {@code null} if exact division is not possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider} or {@code null} if exact division is not possible
     */
    default E divideOrNull(E dividend, E divider) {
        if (isOne(divider))
            return dividend;
        E[] qd = divideAndRemainder(dividend, divider);
        if (qd == null)
            return null;
        if (!isZero(qd[1]))
            return null;
        return qd[0];
    }

    /**
     * Divides {@code dividend} by {@code divider} or throws {@code ArithmeticException} if exact division is not
     * possible
     *
     * @param dividend the dividend
     * @param divider  the divider
     * @return {@code dividend / divider}
     * @throws ArithmeticException if exact division is not possible
     */
    default E divideExact(E dividend, E divider) {
        E result = divideOrNull(dividend, divider);
        if (result == null)
            throw new ArithmeticException("not divisible: " + dividend + " / " + divider);
        return result;
    }

    /**
     * Gives the inverse element {@code element ^ (-1) }
     *
     * @param element the element
     * @return {@code element ^ (-1)}
     */
    E reciprocal(E element);

    /**
     * Returns the greatest common divisor of two elements
     *
     * @param a the first element
     * @param b the second element
     * @return gcd
     */
    default E gcd(E a, E b) {
        if (isField())
            return a;
        if (!isEuclideanRing())
            throw new UnsupportedOperationException("GCD is not supported in this ring");
        if (isZero(a)) return b;
        if (isZero(b)) return a;

        // run Euclidean algorithm by default
        E x = a, y = b, r;
        while (true) {
            r = remainder(x, y);
            if (r == null)
                throw new ArithmeticException("Not divisible with remainder: (" + x + ") / (" + y + ")");

            if (isZero(r))
                break;
            x = y;
            y = r;
        }
        return y;
    }

    /**
     * Returns array of {@code [gcd(a,b), s, t]} such that {@code s * a + t * b = gcd(a, b)}
     *
     * @throws UnsupportedOperationException if this is not the Euclidean ring and there is no special implementation
     *                                       provided by particular subtype
     */
    default E[] extendedGCD(E a, E b) {
        if (!isEuclideanRing())
            throw new UnsupportedOperationException("Extended GCD is not supported in this ring");

        if (isZero(a)) return createArray(b, getOne(), getOne());
        if (isZero(b)) return createArray(a, getOne(), getOne());

        if (isField()) {
            E[] r = createArray(3);
            r[0] = getOne();
            r[1] = divideExact(reciprocal(a), valueOf(2));
            r[2] = divideExact(reciprocal(b), valueOf(2));
            return r;
        }

        E s = getZero(), old_s = getOne();
        E t = getOne(), old_t = getZero();
        E r = b, old_r = a;

        E q;
        E tmp;
        while (!isZero(r)) {
            q = quotient(old_r, r);
            if (q == null)
                throw new ArithmeticException("Not divisible with remainder: (" + old_r + ") / (" + r + ")");

            tmp = old_r;
            old_r = r;
            r = subtract(tmp, multiply(q, r));

            tmp = old_s;
            old_s = s;
            s = subtract(tmp, multiply(q, s));

            tmp = old_t;
            old_t = t;
            t = subtract(tmp, multiply(q, t));
        }

        E[] result = createArray(3);
        result[0] = old_r;
        result[1] = old_s;
        result[2] = old_t;
        return result;
    }

    /**
     * Returns the least common multiple of two elements
     *
     * @param a the first element
     * @param b the second element
     * @return lcm
     */
    default E lcm(E a, E b) {
        if (isZero(a) || isZero(b))
            return getZero();
        return multiply(divideExact(a, gcd(a, b)), b);
    }

    /**
     * Returns greatest common divisor of specified elements
     *
     * @param elements the elements
     * @return gcd
     */
    @SuppressWarnings("unchecked")
    default E gcd(E... elements) {
        return gcd(Arrays.asList(elements));
    }

    /**
     * Returns greatest common divisor of specified elements
     *
     * @param elements the elements
     * @return gcd
     */
    default E gcd(Iterable<E> elements) {
        E gcd = null;
        for (E e : elements) {
            if (gcd == null)
                gcd = e;
            else
                gcd = gcd(gcd, e);
        }
        return gcd;
    }

    /**
     * Square-free factorization of specified element
     */
    default FactorDecomposition<E> factorSquareFree(E element) {
        throw new UnsupportedOperationException();
    }

    /**
     * Factor specified element
     */
    default FactorDecomposition<E> factor(E element) {
        throw new UnsupportedOperationException();
    }

    /**
     * Returns zero element of this ring
     *
     * @return 0
     */
    E getZero();

    /**
     * Returns unit element of this ring (one)
     *
     * @return 1
     */
    E getOne();

    /**
     * Returns negative unit element of this ring (minus one)
     *
     * @return -1
     */
    default E getNegativeOne() {
        return negate(getOne());
    }

    /**
     * Tests whether specified element is zero
     *
     * @param element the ring element
     * @return whether specified element is zero
     */
    boolean isZero(E element);

    /**
     * Tests whether specified element is one (exactly)
     *
     * @param element the ring element
     * @return whether specified element is exactly one
     * @see #isUnit(Object)
     */
    boolean isOne(E element);

    /**
     * Tests whether specified element is a ring unit
     *
     * @param element the ring element
     * @return whether specified element is a ring unit
     * @see #isOne(Object)
     */
    boolean isUnit(E element);

    /**
     * Tests whether specified element is a ring unit or zero
     *
     * @param element the ring element
     * @return whether specified element is a ring unit or zero
     */
    default boolean isUnitOrZero(E element) {
        return isUnit(element) || isZero(element);
    }

    /**
     * Tests whether specified element is minus one
     *
     * @param e the ring element
     * @return whether specified element is minus one
     */
    default boolean isMinusOne(E e) {
        return getNegativeOne().equals(e);
    }

    /**
     * Returns ring element associated with specified {@code long}
     *
     * @param val machine integer
     * @return ring element associated with specified {@code long}
     */
    E valueOf(long val);

    /**
     * Returns ring element associated with specified integer
     *
     * @param val integer
     * @return ring element associated with specified integer
     */
    E valueOfBigInteger(BigInteger val);

    /**
     * Converts array of machine integers to ring elements via {@link #valueOf(long)}
     *
     * @param elements array of machine integers
     * @return array of ring elements
     */
    default E[] valueOf(long[] elements) {
        E[] array = createArray(elements.length);
        for (int i = 0; i < elements.length; i++)
            array[i] = valueOf(elements[i]);
        return array;
    }

    /**
     * Converts a value from other ring to this ring.
     *
     * @param val some element from any ring
     * @return this ring element associated with specified {@code val}
     */
    E valueOf(E val);

    /**
     * Applies {@link #valueOf(Object)} inplace to the specified array
     *
     * @param elements the array
     */
    default void setToValueOf(E[] elements) {
        for (int i = 0; i < elements.length; i++)
            elements[i] = valueOf(elements[i]);
    }

    /**
     * Parse string into ring element
     *
     * @param string string
     * @return ring element
     */
    @Override
    default E parse(String string) {
        throw new UnsupportedOperationException();
    }

    /**
     * Creates generic array of ring elements of specified length
     *
     * @param length array length
     * @return array of ring elements of specified {@code length}
     */
    @SuppressWarnings("unchecked")
    default E[] createArray(int length) {
        return (E[]) Array.newInstance(getOne().getClass(), length);
    }

    /**
     * Creates 2d array of ring elements of specified length
     *
     * @param length array length
     * @return 2d array of ring elements of specified {@code length}
     */
    @SuppressWarnings("unchecked")
    default E[][] createArray2d(int length) {
        return (E[][]) Array.newInstance(createArray(0).getClass(), length);
    }

    /**
     * Creates 2d array of ring elements of specified shape
     *
     * @param m result length
     * @param n length of each array in the result
     * @return 2d array E[m][n]
     */
    @SuppressWarnings("unchecked")
    default E[][] createArray2d(int m, int n) {
        return (E[][]) Array.newInstance(getOne().getClass(), m, n);
    }

    /**
     * Creates array filled with zero elements
     *
     * @param length array length
     * @return array filled with zero elements of specified {@code length}
     */
    @SuppressWarnings("unchecked")
    default E[] createZeroesArray(int length) {
        E[] array = createArray(length);
        for (int i = 0; i < array.length; i++)
            // NOTE: getZero() is invoked each time in a loop in order to fill array with unique elements
            array[i] = getZero();
        return array;
    }

    /**
     * Creates 2d array of ring elements of specified shape filled with zero elements
     *
     * @param m result length
     * @param n length of each array in the result
     * @return 2d array E[m][n] filled with zero elements
     */
    @SuppressWarnings("unchecked")
    default E[][] createZeroesArray2d(int m, int n) {
        E[][] arr = createArray2d(m, n);
        for (E[] a : arr)
            for (int i = 0; i < a.length; ++i)
                a[i] = getZero();
        return arr;
    }

    /**
     * Creates generic array of {@code {a, b}}
     *
     * @param a the first element of array
     * @param b the second element of array
     * @return array {@code {a,b}}
     */
    @SuppressWarnings("unchecked")
    default E[] createArray(E a, E b) {
        E[] array = createArray(2);
        array[0] = a;
        array[1] = b;
        return array;
    }

    /**
     * Creates generic array of {@code {a, b, c}}
     */
    @SuppressWarnings("unchecked")
    default E[] createArray(E a, E b, E c) {
        E[] array = createArray(3);
        array[0] = a;
        array[1] = b;
        array[2] = c;
        return array;
    }

    /**
     * Creates generic array with single element
     *
     * @param element the element
     * @return array with single specified element
     */
    @SuppressWarnings("unchecked")
    default E[] createArray(E element) {
        E[] array = createArray(1);
        array[0] = element;
        return array;
    }

    /**
     * Returns {@code base} in a power of {@code exponent} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code exponent}
     */
    default E pow(E base, int exponent) {
        return pow(base, BigInteger.valueOf(exponent));
    }

    /**
     * Returns {@code base} in a power of {@code exponent} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code exponent}
     */
    default E pow(E base, long exponent) {
        return pow(base, BigInteger.valueOf(exponent));
    }

    /**
     * Returns {@code base} in a power of {@code exponent} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code exponent}
     */
    default E pow(E base, BigInteger exponent) {
        if (exponent.signum() < 0)
            throw new IllegalArgumentException();

        if (exponent.isOne())
            return base;

        E result = getOne();
        E k2p = copy(base); // <= copy the base (mutable operations are used below)
        for (; ; ) {
            if ((exponent.testBit(0)))
                result = multiplyMutable(result, k2p);
            exponent = exponent.shiftRight(1);
            if (exponent.isZero())
                return result;
            k2p = multiplyMutable(k2p, k2p);
        }
    }

    /**
     * Gives a product of {@code valueOf(1) * valueOf(2) * .... * valueOf(num) }
     *
     * @param num the number
     * @return {@code valueOf(1) * valueOf(2) * .... * valueOf(num) }
     */
    default E factorial(long num) {
        E result = getOne();
        for (int i = 2; i <= num; ++i)
            result = multiplyMutable(result, valueOf(i));
        return result;
    }

//    /**
//     * Returns the element which is next to the specified {@code element} (according to {@link #compare(Object, Object)})
//     * or {@code null} in the case of infinite cardinality
//     *
//     * @param element the element
//     * @return next element
//     */
//    E nextElement(E element);

    /**
     * Returns iterator over ring elements (for finite rings, otherwise throws exception)
     */
    @Override
    Iterator<E> iterator();

    /**
     * Returns a random element from this ring
     *
     * @return random element from this ring
     */
    default E randomElement() { return randomElement(Rings.privateRandom);}

    /**
     * Returns a random element from this ring
     *
     * @param rnd the source of randomness
     * @return random element from this ring
     */
    default E randomElement(RandomGenerator rnd) { return valueOf(rnd.nextLong());}

    /**
     * Returns a random non zero element from this ring
     *
     * @param rnd the source of randomness
     * @return random non zero element from this ring
     */
    default E randomNonZeroElement(RandomGenerator rnd) {
        E el;
        do {
            el = randomElement(rnd);
        } while (isZero(el));
        return el;
    }

    @Override
    default String toString(E element) {
        return element.toString();
    }

//    /**
//     * Returns ring with larger cardinality that contains all elements of this or null if there is no such ring.
//     *
//     * @return ring with larger cardinality that contains all elements of this or null if there is no such ring
//     */
//    Ring<E> getExtension();
}
