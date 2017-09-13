package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.Rational;
import cc.r2.core.poly.multivar.MonomialOrder;
import cc.r2.core.poly.multivar.MultivariatePolynomial;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.IrreduciblePolynomials;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import org.apache.commons.math3.random.Well19937c;

/**
 * Shortcut methods for basic algebraic domains.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class Domains {
    private Domains() {}

    /**
     * Domain of integers (Z)
     */
    public static final Integers Z = Integers.Integers;

    /**
     * Domain of rationals (Q)
     */
    public static final Rationals Q = Rationals.Rationals;

    /**
     * Domain of integers modulo {@code modulus} (with modulus < 2^63)
     *
     * @param modulus the modulus
     */
    public static lIntegersModulo Zp(long modulus) {return new lIntegersModulo(modulus);}

    /**
     * Domain of integers modulo {@code modulus} (arbitrary large modulus)
     *
     * @param modulus the modulus (arbitrary large)
     */
    public static IntegersModulo Zp(BigInteger modulus) {return new IntegersModulo(modulus);}

    /**
     * Galois field with the cardinality {@code prime ^ exponent} (with prime < 2^63)
     *
     * @param prime    the integer prime modulus
     * @param exponent the exponent (degree of modulo polynomial)
     */
    public static FiniteField<lUnivariatePolynomialZp> GF(long prime, int exponent) {
        if (exponent <= 0)
            throw new IllegalArgumentException("Exponent must be positive");
        // provide random generator with fixed seed to make the behavior predictable
        return new FiniteField<>(IrreduciblePolynomials.randomIrreduciblePolynomial(prime, exponent, new Well19937c(0x77f3dfae)));
    }

    /**
     * Galois field with the cardinality {@code prime ^ exponent} for arbitrary large {@code prime}
     *
     * @param prime    the integer (arbitrary large) prime modulus
     * @param exponent the exponent (degree of modulo polynomial)
     */
    public static FiniteField<UnivariatePolynomial<BigInteger>> GF(BigInteger prime, int exponent) {
        if (exponent <= 0)
            throw new IllegalArgumentException("Exponent must be positive");
        // provide random generator with fixed seed to make the behavior predictable
        return new FiniteField<>(IrreduciblePolynomials.randomIrreduciblePolynomial(Zp(prime), exponent, new Well19937c(0x77f3dfae)));
    }

    /**
     * Galois field with the specified irreducible generator. Note: there is no explicit check that input is irreducible
     *
     * @param irreducible irreducible univariate polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> FiniteField<Poly> GF(Poly irreducible) {
        return new FiniteField<>(irreducible);
    }

    /**
     * Domain of univariate polynomials over specified coefficient domain
     *
     * @param coefficientDomain the coefficient domain
     */
    public static <E> UnivariatePolynomials<UnivariatePolynomial<E>> Polynomials(Domain<E> coefficientDomain) {
        return new UnivariatePolynomials<>(UnivariatePolynomial.zero(coefficientDomain));
    }

    /**
     * Domain of univariate polynomials over integers (Z[x])
     */
    public static final UnivariatePolynomials<UnivariatePolynomial<BigInteger>> PolynomialsZ = Polynomials(Z);

    /**
     * Domain of univariate polynomials over rationals (Q[x])
     */
    public static final UnivariatePolynomials<UnivariatePolynomial<Rational>> PolynomialsQ = Polynomials(Q);

    /**
     * Domain of univariate polynomials over Zp integers (Zp[x])
     *
     * @param modulus the modulus
     */
    public static UnivariatePolynomials<lUnivariatePolynomialZp> PolynomialsZp(long modulus) {
        return new UnivariatePolynomials<>(lUnivariatePolynomialZp.zero(modulus));
    }

    /**
     * Domain of univariate polynomials over Zp integers (Zp[x])
     *
     * @param modulus the modulus
     */
    public static UnivariatePolynomials<lUnivariatePolynomialZp> PolynomialsZp(lIntegersModulo modulus) {
        return new UnivariatePolynomials<>(lUnivariatePolynomialZp.zero(modulus));
    }

    /**
     * Domain of univariate polynomials over Zp integers (Zp[x]) with arbitrary large modulus
     *
     * @param modulus the modulus (arbitrary large)
     */
    public static UnivariatePolynomials<UnivariatePolynomial<BigInteger>> PolynomialsZp(BigInteger modulus) {
        return Polynomials(Zp(modulus));
    }

    /**
     * Domain of univariate polynomials over Galois Field (GF[x])
     *
     * @param gf the finite field
     */
    public static <uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariatePolynomials<UnivariatePolynomial<uPoly>> PolynomialsGF(FiniteField<uPoly> gf) { return Polynomials(gf); }

    /**
     * Domain of multivariate polynomials with specified number of variables over specified coefficient domain
     *
     * @param nVariables        the number of variables
     * @param coefficientDomain the coefficient domain
     */
    public static <E> MultivariatePolynomials<MultivariatePolynomial<E>>
    MultivariatePolynomials(int nVariables, Domain<E> coefficientDomain) {
        return new MultivariatePolynomials<>(MultivariatePolynomial.zero(nVariables, coefficientDomain, MonomialOrder.LEX));
    }

    /**
     * Domain of multivariate polynomials over integers (Z[x1, x2, ...])
     *
     * @param nVariables the number of variables
     */
    public static MultivariatePolynomials<MultivariatePolynomial<BigInteger>>
    MultivariatePolynomialsZ(int nVariables) {
        return MultivariatePolynomials(nVariables, Z);
    }

    /**
     * Domain of multivariate polynomials over rationals (Q[x1, x2, ...])
     *
     * @param nVariables the number of variables
     */
    public static MultivariatePolynomials<MultivariatePolynomial<Rational>>
    MultivariatePolynomialsQ(int nVariables) {
        return MultivariatePolynomials(nVariables, Q);
    }

    /**
     * Domain of multivariate polynomials over Zp integers (Zp[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus
     */
    public static MultivariatePolynomials<lMultivariatePolynomialZp>
    MultivariatePolynomialsZp(int nVariables, long modulus) {
        return new MultivariatePolynomials<>(lMultivariatePolynomialZp.zero(nVariables, Zp(modulus), MonomialOrder.LEX));
    }

    /**
     * Domain of multivariate polynomials over Zp integers (Zp[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus
     */
    public static MultivariatePolynomials<lMultivariatePolynomialZp>
    MultivariatePolynomialsZp(int nVariables, lIntegersModulo modulus) {
        return new MultivariatePolynomials<>(lMultivariatePolynomialZp.zero(nVariables, modulus, MonomialOrder.LEX));
    }

    /**
     * Domain of multivariate polynomials over Zp integers (Zp[x1, x2, ...]) with arbitrary large modulus
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus (arbitrary large)
     */
    public static MultivariatePolynomials<MultivariatePolynomial<BigInteger>>
    MultivariatePolynomialsZp(int nVariables, BigInteger modulus) {
        return MultivariatePolynomials(nVariables, Zp(modulus));
    }

    /**
     * Domain of multivariate polynomials over Galois Field (GF[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param gf         the finite field
     */
    public static <uPoly extends IUnivariatePolynomial<uPoly>>
    MultivariatePolynomials<MultivariatePolynomial<uPoly>>
    MultivariatePolynomialsGF(int nVariables, FiniteField<uPoly> gf) {
        return MultivariatePolynomials(nVariables, gf);
    }
}
