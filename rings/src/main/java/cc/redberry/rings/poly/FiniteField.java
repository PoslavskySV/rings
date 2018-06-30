package cc.redberry.rings.poly;

import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;

/**
 * Galois field.
 *
 * @param <Poly> type of polynomials representing elements of this Galois field
 * @see AlgebraicNumberField
 * @since 1.0
 */
public final class FiniteField<Poly extends IUnivariatePolynomial<Poly>>
        extends UnivariateQuotientRing<Poly> {
    private static final long serialVersionUID = 1L;
    /** GF(3^3) */
    public static final FiniteField<UnivariatePolynomialZp64> GF27 = new FiniteField<>(UnivariatePolynomialZ64.create(-1, -1, 0, 1).modulus(3));
    /** GF(17^5) */
    public static final FiniteField<UnivariatePolynomialZp64> GF17p5 = new FiniteField<>(UnivariatePolynomialZ64.create(11, 11, 0, 3, 9, 9).modulus(17).monic());

    /**
     * Constructs finite field from the specified irreducible polynomial. NOTE: irreducibility test for the input
     * polynomial is not performed here, use {@link IrreduciblePolynomials#irreducibleQ(IUnivariatePolynomial)} to test
     * irreducibility.
     *
     * @param irreducible irreducible polynomial
     * @see IrreduciblePolynomials
     */
    public FiniteField(Poly irreducible) {
        super(irreducible);
        if (!irreducible.isOverFiniteField())
            throw new IllegalArgumentException("Irreducible poly must be over finite field.");
    }

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public boolean isUnit(Poly element) {
        return !element.isZero();
    }

    @Override
    public Poly gcd(Poly a, Poly b) {
        return a;
    }

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {
        return a.createArray(multiply(a, reciprocal(b)), getZero());
    }

    @Override
    public Poly remainder(Poly dividend, Poly divider) {
        return getZero();
    }
}
