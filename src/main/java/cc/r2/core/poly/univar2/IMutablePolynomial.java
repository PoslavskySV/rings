package cc.r2.core.poly.univar2;

import cc.r2.core.poly.IGeneralPolynomial;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface IMutablePolynomial<Poly extends IMutablePolynomial<Poly>> extends IGeneralPolynomial<Poly>, Comparable<Poly> {
    /**
     * Check whether {@code oth} and {@code this} belongs to the same Zp domain
     *
     * @param oth other polynomial
     */
    void checkCompatible(Poly oth);

    /**
     * Divides this polynomial by the leading coefficient of {@code other} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param other other polynomial
     * @return {@code this} divided by the {@code other.lc()} or {@code null}
     */
    Poly divideByLC(Poly other);

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
    Poly shiftLeft(int offset);

    /**
     * Multiplies {@code this} by the {@code x^offset}.
     *
     * @param offset monomial exponent
     * @return {@code this * x^offset}
     */
    Poly shiftRight(int offset);

    /**
     * Returns the remainder {@code this rem x^(newDegree + 1)}, it is polynomial with the coefficient list truncated
     * to the new degree {@code newDegree}.
     *
     * @param newDegree new degree
     * @return remainder {@code this rem x^(newDegree + 1)}
     */
    Poly truncate(int newDegree);

    /**
     * Creates polynomial formed from the coefficients of this starting from {@code from} (inclusive) to {@code to} (exclusive)
     *
     * @param from the initial index of the range to be copied, inclusive
     * @param to   the final index of the range to be copied, exclusive.
     * @return polynomial formed from the range of coefficients
     */
    Poly getRange(int from, int to);

    /**
     * Reverses the coefficients of this
     *
     * @return reversed polynomial
     */
    Poly reverse();

    /**
     * Creates monomial {@code x^degree}
     *
     * @param degree monomial degree
     * @return {@code coefficient * x^degree}
     */
    Poly createMonomial(int degree);

    /** overcome Java generics... */
    Poly[] arrayNewInstance(int length);

    /** overcome Java generics... */
    default Poly[] arrayNewInstance(Poly a, Poly b) {
        Poly[] r = arrayNewInstance(2);
        r[0] = a; r[1] = b;
        return r;
    }

    /**
     * Returns the formal derivative of this poly
     *
     * @return the formal derivative
     */
    Poly derivative();

    /**
     * Deep copy of this
     *
     * @return deep copy of this
     */
    Poly clone();
}
