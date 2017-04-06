package cc.r2.core.poly2;

import cc.r2.core.number.BigInteger;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
interface IMutablePolynomialZp<T extends IMutablePolynomialZp> extends IMutablePolynomial<T> {
    void checkCompatibleModulus(T oth);

    /**
     * Sets {@code this} to its monic part (that is {@code this} multiplied by its inversed leading coefficient).
     *
     * @return {@code this}
     */
    T monic();

    /**
     * Multiplies {@code this} by {@code modInverse(other.lc())}
     *
     * @param other other polynomial
     * @return {@code this} multiplied by {@code modInverse(other.lc())}
     */
    T divideByLC(T other);

    /**
     * Returns BigInteger value of the modulus.
     *
     * @return BigInteger value of the modulus
     */
    BigInteger modulusAsBigInt();

    <PolyZ extends IMutablePolynomialZ<PolyZ>> PolyZ normalSymmetricForm();

    <PolyZ extends IMutablePolynomialZ<PolyZ>> PolyZ normalForm(boolean copy);
}
