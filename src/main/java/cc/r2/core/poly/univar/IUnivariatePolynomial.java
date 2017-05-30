package cc.r2.core.poly.univar;

import cc.r2.core.poly.IGeneralPolynomial;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface IUnivariatePolynomial<Poly extends IUnivariatePolynomial<Poly>> extends IGeneralPolynomial<Poly> {
    /**
     * Returns whether i-th coefficient of this is zero
     *
     * @param i the position
     * @return whether i-th coefficient of this is zero
     */
    boolean isZeroAt(int i);

    /**
     * Returns a set of exponents of non-zero terms
     *
     * @return a set of exponents of non-zero terms
     */
    default TIntSet exponents() {
        TIntHashSet degrees = new TIntHashSet();
        for (int i = degree(); i >= 0; --i)
            if (!isZeroAt(i))
                degrees.add(i);
        return degrees;
    }

    /**
     * Divides this polynomial by the leading coefficient of {@code other} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code other.lc()}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param other the polynomial
     * @return {@code this} divided by the {@code other.lc()} or {@code null} if exact division is not possible
     */
    Poly divideByLC(Poly other);

    /**
     * Returns position of the first non-zero coefficient, that is common monomial exponent (e.g. 2 for x^2 + x^3 + ...). In the case of zero polynomial, -1 returned
     *
     * @return position of the first non-zero coefficient or -1 if this is zero
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

    /**
     * Returns the formal derivative of this poly
     *
     * @return the formal derivative
     */
    Poly derivative();

    /** content as a constant poly */
    Poly contentAsPoly();

    /** leading coefficient as a constant poly */
    Poly lcAsPoly();

    /**
     * Deep copy of this
     *
     * @return deep copy of this
     */
    @Override
    Poly clone();
}
