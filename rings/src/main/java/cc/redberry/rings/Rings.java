package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.*;
import cc.redberry.rings.poly.multivar.*;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.IrreduciblePolynomials;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.random.Well44497b;

import java.util.Comparator;

/**
 * Common rings.
 *
 * @since 1.0
 */
public final class Rings {
    private Rings() {}

    static RandomGenerator privateRandom = new Well44497b(System.nanoTime());

    /**
     * Ring of integers (Z)
     */
    public static final Integers Z = Integers.Integers;

    /**
     * Ring of rationals (Q)
     */
    public static final Rationals<BigInteger> Q = new Rationals<>(Z);

    /**
     * Ring of rational functions over specified ring
     *
     * @param ring the ring that numerators and denominators belong to
     */
    public static <E> Rationals<E> Frac(Ring<E> ring) {
        return new Rationals<>(ring);
    }

    /**
     * Ring of integers modulo {@code modulus} (with modulus < 2^63)
     *
     * @param modulus the modulus
     */
    public static IntegersZp64 Zp64(long modulus) {return new IntegersZp64(modulus);}

    /**
     * Ring of integers modulo {@code modulus} (arbitrary large modulus)
     *
     * @param modulus the modulus (arbitrary large)
     */
    public static IntegersZp Zp(long modulus) {return new IntegersZp(modulus);}

    /**
     * Ring of integers modulo {@code modulus} (arbitrary large modulus)
     *
     * @param modulus the modulus (arbitrary large)
     */
    public static IntegersZp Zp(BigInteger modulus) {return new IntegersZp(modulus);}

    /**
     * Galois field with the cardinality {@code prime ^ exponent} (with prime < 2^63).
     *
     * @param prime    the integer prime modulus
     * @param exponent the exponent (degree of modulo polynomial)
     */
    public static FiniteField<UnivariatePolynomialZp64> GF(long prime, int exponent) {
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
     * Galois field with the specified irreducible generator. Note: there is no explicit check that input is
     * irreducible
     *
     * @param irreducible irreducible univariate polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> FiniteField<Poly> GF(Poly irreducible) {
        return new FiniteField<>(irreducible);
    }

    /**
     * Ring of univariate polynomials over specified coefficient ring
     *
     * @param coefficientRing the coefficient ring
     */
    public static <E> UnivariateRing<UnivariatePolynomial<E>> UnivariateRing(Ring<E> coefficientRing) {
        return new UnivariateRing<>(UnivariatePolynomial.zero(coefficientRing));
    }

    /**
     * Ring of univariate polynomials with specified factory
     *
     * @param factory factory
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> PolynomialRing<Poly> UnivariateRing(Poly factory) {
        return new UnivariateRing<>(factory);
    }

    /**
     * Ring of univariate polynomials over integers (Z[x])
     */
    public static final UnivariateRing<UnivariatePolynomial<BigInteger>> UnivariateRingZ = UnivariateRing(Z);

    /**
     * Ring of univariate polynomials over rationals (Q[x])
     */
    public static final UnivariateRing<UnivariatePolynomial<Rational<BigInteger>>> UnivariateRingQ = UnivariateRing(Q);

    /**
     * Ring of univariate polynomials over Zp integers (Zp[x])
     *
     * @param modulus the modulus
     */
    public static UnivariateRing<UnivariatePolynomialZp64> UnivariateRingZp64(long modulus) {
        return new UnivariateRing<>(UnivariatePolynomialZp64.zero(modulus));
    }

    /**
     * Ring of univariate polynomials over Zp integers (Zp[x])
     *
     * @param modulus the modulus
     */
    public static UnivariateRing<UnivariatePolynomialZp64> UnivariateRingZp64(IntegersZp64 modulus) {
        return new UnivariateRing<>(UnivariatePolynomialZp64.zero(modulus));
    }

    /**
     * Ring of univariate polynomials over Zp integers (Zp[x]) with arbitrary large modulus
     *
     * @param modulus the modulus (arbitrary large)
     */
    public static UnivariateRing<UnivariatePolynomial<BigInteger>> UnivariateRingZp(BigInteger modulus) {
        return UnivariateRing(Zp(modulus));
    }

    /**
     * Ring of univariate polynomials over Galois Field (GF[x])
     *
     * @param gf the finite field
     */
    public static <uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariateRing<UnivariatePolynomial<uPoly>> UnivariateRingGF(FiniteField<uPoly> gf) { return UnivariateRing(gf); }

    /**
     * Ring of multivariate polynomials with specified number of variables over specified coefficient ring
     *
     * @param nVariables      the number of variables
     * @param coefficientRing the coefficient ring
     * @param monomialOrder   the monomial order
     */
    public static <E> MultivariateRing<MultivariatePolynomial<E>>
    MultivariateRing(int nVariables, Ring<E> coefficientRing, Comparator<DegreeVector> monomialOrder) {
        return new MultivariateRing<>(MultivariatePolynomial.zero(nVariables, coefficientRing, monomialOrder));
    }

    /**
     * Ring of multivariate polynomials with specified number of variables over specified coefficient ring
     *
     * @param nVariables      the number of variables
     * @param coefficientRing the coefficient ring
     */
    public static <E> MultivariateRing<MultivariatePolynomial<E>>
    MultivariateRing(int nVariables, Ring<E> coefficientRing) {
        return MultivariateRing(nVariables, coefficientRing, MonomialOrder.LEX);
    }

    /**
     * Ring of multivariate polynomials with specified factory
     *
     * @param factory factory
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    PolynomialRing<Poly> MultivariateRing(Poly factory) {
        return new MultivariateRing<>(factory);
    }

    /**
     * Ring of multivariate polynomials over integers (Z[x1, x2, ...])
     *
     * @param nVariables the number of variables
     */
    public static MultivariateRing<MultivariatePolynomial<BigInteger>>
    MultivariateRingZ(int nVariables) {
        return MultivariateRing(nVariables, Z);
    }

    /**
     * Ring of multivariate polynomials over rationals (Q[x1, x2, ...])
     *
     * @param nVariables the number of variables
     */
    public static MultivariateRing<MultivariatePolynomial<Rational<BigInteger>>>
    MultivariateRingQ(int nVariables) {
        return MultivariateRing(nVariables, Q);
    }

    /**
     * Ring of multivariate polynomials over Zp integers (Zp[x1, x2, ...])
     *
     * @param nVariables    the number of variables
     * @param modulus       the modulus
     * @param monomialOrder the monomial order
     */
    public static MultivariateRing<MultivariatePolynomialZp64>
    MultivariateRingZp64(int nVariables, long modulus, Comparator<DegreeVector> monomialOrder) {
        return new MultivariateRing<>(MultivariatePolynomialZp64.zero(nVariables, Zp64(modulus), monomialOrder));
    }

    /**
     * Ring of multivariate polynomials over Zp machine integers (Zp[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus
     */
    public static MultivariateRing<MultivariatePolynomialZp64>
    MultivariateRingZp64(int nVariables, long modulus) {
        return MultivariateRingZp64(nVariables, modulus, MonomialOrder.LEX);
    }

    /**
     * Ring of multivariate polynomials over Zp integers (Zp[x1, x2, ...])
     *
     * @param nVariables    the number of variables
     * @param modulus       the modulus
     * @param monomialOrder monomial order
     */
    public static MultivariateRing<MultivariatePolynomialZp64>
    MultivariateRingZp64(int nVariables, IntegersZp64 modulus, Comparator<DegreeVector> monomialOrder) {
        return new MultivariateRing<>(MultivariatePolynomialZp64.zero(nVariables, modulus, monomialOrder));
    }

    /**
     * Ring of multivariate polynomials over Zp integers (Zp[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus
     */
    public static MultivariateRing<MultivariatePolynomialZp64>
    MultivariateRingZp64(int nVariables, IntegersZp64 modulus) {
        return MultivariateRingZp64(nVariables, modulus, MonomialOrder.LEX);
    }

    /**
     * Ring of multivariate polynomials over Zp integers (Zp[x1, x2, ...]) with arbitrary large modulus
     *
     * @param nVariables the number of variables
     * @param modulus    the modulus (arbitrary large)
     */
    public static MultivariateRing<MultivariatePolynomial<BigInteger>>
    MultivariateRingZp(int nVariables, BigInteger modulus) {
        return MultivariateRing(nVariables, Zp(modulus));
    }

    /**
     * Ring of multivariate polynomials over Galois Field (GF[x1, x2, ...])
     *
     * @param nVariables the number of variables
     * @param gf         the finite field
     */
    public static <uPoly extends IUnivariatePolynomial<uPoly>>
    MultivariateRing<MultivariatePolynomial<uPoly>>
    MultivariateRingGF(int nVariables, FiniteField<uPoly> gf) {
        return MultivariateRing(nVariables, gf);
    }


    /**
     * Generic factory for polynomial ring
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IPolynomial<Poly>> PolynomialRing<Poly> PolynomialRing(Poly factory) {
        if (factory instanceof IUnivariatePolynomial)
            return (PolynomialRing<Poly>) UnivariateRing((IUnivariatePolynomial) factory);
        else
            return (PolynomialRing<Poly>) MultivariateRing((AMultivariatePolynomial) factory);
    }
}
