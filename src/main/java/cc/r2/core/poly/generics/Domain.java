package cc.r2.core.poly.generics;

import java.lang.reflect.Array;
import java.util.Comparator;

/**
 * Domain specification.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface Domain<E> extends Comparator<E> {
    /**
     * Adds two elements
     *
     * @param a the first element
     * @param b the second element
     * @return a + b
     */
    E add(E a, E b);

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
     * Negates a given element
     *
     * @param val the domain element
     * @return -val
     */
    E negate(E val);

    /**
     * Returns -1 if {@code a < 0}, 0 if {@code a == 0} and 1 if {@code a > 0}, where comparison is defined in terms
     * of {@link #compare(Object, Object)}
     *
     * @param a the element
     * @return -1 if {@code a < 0}, 0 if {@code a == 0} and 1 otherwise
     */
    int signum(E a);

    /**
     * Returns quotient and remainder of {@code a/b}
     *
     * @param a the dividend
     * @param b the divider
     * @return {@code {quotient, remainder}}
     */
    E[] divideAndRemainder(E a, E b);

    /**
     * Returns greatest common divisor of two elements
     *
     * @param a the first element
     * @param b the second element
     * @return gcd
     */
    E gcd(E a, E b);

    /**
     * Domain element representing zero
     *
     * @return 0
     */
    E getZero();

    /**
     * Domain element representing one
     *
     * @return 1
     */
    E getOne();

    /**
     * Domain element representing minus one
     *
     * @return -1
     */
    default E getNegativeOne() {
        return negate(getOne());
    }

    /**
     * Tests whether specified element is zero
     *
     * @param e the domain element
     * @return whether specified element is zero
     */
    boolean isZero(E e);

    /**
     * Tests whether specified element is one
     *
     * @param e the domain element
     * @return whether specified element is one
     */
    boolean isOne(E e);

    /**
     * Tests whether specified element is minus one
     *
     * @param e the domain element
     * @return whether specified element is minus one
     */
    default boolean isMinusOne(E e) {
        return negate(getOne()).equals(e);
    }

    /**
     * Value of machine integer (if supported)
     *
     * @param val machine integer
     * @return domain element associated with specified {@code long}
     */
    E valueOf(long val);

    /**
     * Converts a value from possibly other domain to this domain.
     *
     * @param val some element from any domain
     * @return this domain element associated with specified {@code val}
     */
    E valueOf(E val);

    /**
     * Parse string into domain element
     *
     * @param string string
     * @return domain element
     */
    default E parse(String string) {
        throw new UnsupportedOperationException();
    }

    /**
     * Creates generic array of domain elements of specified length
     *
     * @param length array length
     * @return empty array of domain elements of specified {@code length}
     */
    @SuppressWarnings("unchecked")
    default E[] createArray(int length) {
        return (E[]) Array.newInstance(getOne().getClass(), length);
    }

    /**
     * Returns {@code base} in a power of {@code exponent} (non negative)
     *
     * @param base     base
     * @param exponent exponent (non negative)
     * @return {@code base} in a power of {@code exponent}
     * @throws ArithmeticException if the result overflows a long
     */
    default E pow(E base, int exponent) {
        if (exponent < 0)
            throw new IllegalArgumentException();

        E result = getOne();
        E k2p = base;
        for (; ; ) {
            if ((exponent&1) != 0)
                result = multiply(result, k2p);
            exponent = exponent >> 1;
            if (exponent == 0)
                return result;
            k2p = multiply(k2p, k2p);
        }
    }
}
