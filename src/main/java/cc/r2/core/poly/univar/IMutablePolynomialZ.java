package cc.r2.core.poly.univar;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
interface IMutablePolynomialZ<T extends IMutablePolynomialZ<T>> extends IMutablePolynomial<T> {
    /**
     * Divides this polynomial by the leading coefficient of {@code other} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param other other polynomial
     * @return {@code this} divided by the {@code other.lc()} or {@code null}
     */
    T divideOrNullByLC(T other);

    @Override
    T clone();
}
