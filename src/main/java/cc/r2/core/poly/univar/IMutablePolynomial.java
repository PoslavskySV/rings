package cc.r2.core.poly.univar;

import cc.r2.core.poly.IGeneralPolynomial;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
interface IMutablePolynomial<T extends IMutablePolynomial> extends IGeneralPolynomial<T>, Comparable<T> {

    /**
     * Returns {@code true} if this polynomial is monic
     *
     * @return whether {@code this} is monic
     */
    boolean isMonic();

    /**
     * Gives signum of leading coefficient
     *
     * @return signum of leading coefficient
     */
    int signum();

    /**
     * Returns position of the first non-zero coefficient, that is common monomial exponent (e.g. 2 for x^2 + x^3 + ...)
     *
     * @return position of the first non-zero coefficient
     */
    int firstNonZeroCoefficientPosition();

    /**
     * Returns the quotient {@code this / x^offset}, it is polynomial with coefficient list formed by shifting coefficients
     * of {@code this} to the left by {@code offset}.
     *
     * @param offset shift amount
     * @return the quotient {@code this / x^offset}
     */
    T shiftLeft(int offset);

    /**
     * Multiplies {@code this} by the {@code x^offset}.
     *
     * @param offset monomial exponent
     * @return {@code this * x^offset}
     */
    T shiftRight(int offset);

    /**
     * Returns the remainder {@code this rem x^(newDegree + 1)}, it is polynomial with the coefficient list truncated
     * to the new degree {@code newDegree}.
     *
     * @param newDegree new degree
     * @return remainder {@code this rem x^(newDegree + 1)}
     */
    T truncate(int newDegree);

    /**
     * Creates polynomial formed from the coefficients of this starting from {@code from} (inclusive) to {@code to} (exclusive)
     *
     * @param from the initial index of the range to be copied, inclusive
     * @param to   the final index of the range to be copied, exclusive.
     * @return polynomial formed from the range of coefficients
     */
    T getRange(int from, int to);

    /**
     * Reverses the coefficients of this
     *
     * @return reversed polynomial
     */
    T reverse();

    /**
     * Creates monomial {@code x^degree}
     *
     * @param degree monomial degree
     * @return {@code coefficient * x^degree}
     */
    T createMonomial(int degree);

    /** overcome Java generics... */
    T[] arrayNewInstance(int length);

    /** overcome Java generics... */
    default T[] arrayNewInstance(T a, T b) {
        T[] r = arrayNewInstance(2);
        r[0] = a; r[1] = b;
        return r;
    }

    /**
     * Returns the formal derivative of this poly
     *
     * @return the formal derivative
     */
    T derivative();

    /**
     * Deep copy of this
     *
     * @return deep copy of this
     */
    T clone();
}
