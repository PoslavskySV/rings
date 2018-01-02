package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.util.ArraysUtil;

import java.io.Serializable;

/**
 * Algebraic operations (multiplication, division) and utility methods for monomials.
 *
 * @since 2.3
 */
public interface IMonomialAlgebra<Term extends AMonomial<Term>> extends Serializable {
    /**
     * Multiplies two terms
     */
    Term multiply(Term a, Term b);

    /**
     * Gives quotient {@code dividend / divider } or null is exact division is not possible
     */
    Term divideOrNull(Term dividend, Term divider);

    /**
     * Negates term
     */
    Term negate(Term term);

    /**
     * Whether term is zero
     */
    boolean isZero(Term term);

    /**
     * Whether term is one
     */
    boolean isOne(Term term);

    /**
     * Whether term is unit
     */
    boolean isUnit(Term term);

    /**
     * Whether term is constant
     */
    default boolean isConstant(Term term) {
        return term.isZeroVector();
    }

    /** creates term with specified exponents and unit coefficient */
    Term createTermWithUnitCoefficient(int[] exponents);

    /** creates generic array of specified length */
    Term[] createArray(int length);

    /** creates a unit term */
    Term getUnitTerm(int nVariables);

    /** creates a zero term */
    Term getZeroTerm(int nVariables);

    /**
     * Term algebra for terms over Zp
     */
    class MonomialAlgebraZp64 implements IMonomialAlgebra<MonomialZp64> {
        public final IntegersZp64 ring;

        public MonomialAlgebraZp64(IntegersZp64 ring) {
            this.ring = ring;
        }

        @Override
        public MonomialZp64 multiply(MonomialZp64 a, MonomialZp64 b) {
            DegreeVector dv = a.dvMultiply(b);
            return new MonomialZp64(dv, ring.multiply(a.coefficient, b.coefficient));
        }

        @Override
        public MonomialZp64 divideOrNull(MonomialZp64 dividend, MonomialZp64 divider) {
            DegreeVector dv = dividend.dvDivideOrNull(divider);
            if (dv == null)
                return null;
            return new MonomialZp64(dv, ring.divide(dividend.coefficient, divider.coefficient));
        }

        @Override
        public MonomialZp64 negate(MonomialZp64 term) {
            return term.setCoefficient(ring.negate(term.coefficient));
        }

        @Override
        public boolean isZero(MonomialZp64 term) {
            return term.coefficient == 0;
        }

        @Override
        public boolean isOne(MonomialZp64 term) {
            return isConstant(term) && term.coefficient == 1;
        }

        @Override
        public boolean isUnit(MonomialZp64 term) {
            return isConstant(term);
        }

        @Override
        public MonomialZp64 createTermWithUnitCoefficient(int[] exponents) {
            return new MonomialZp64(exponents, ArraysUtil.sum(exponents), 1L);
        }

        @Override
        public MonomialZp64[] createArray(int length) {
            return new MonomialZp64[length];
        }

        @Override
        public MonomialZp64 getUnitTerm(int nVariables) {
            return new MonomialZp64(nVariables, 1L);
        }

        @Override
        public MonomialZp64 getZeroTerm(int nVariables) {
            return new MonomialZp64(nVariables, 0L);
        }
    }

    /**
     * Generic term algebra
     */
    class MonomialAlgebra<E> implements IMonomialAlgebra<Monomial<E>> {
        public final Ring<E> ring;

        public MonomialAlgebra(Ring<E> ring) {
            this.ring = ring;
        }

        @Override
        public Monomial<E> multiply(Monomial<E> a, Monomial<E> b) {
            DegreeVector dv = a.dvMultiply(b);
            return new Monomial<>(dv, ring.multiply(a.coefficient, b.coefficient));
        }

        @Override
        public Monomial<E> divideOrNull(Monomial<E> dividend, Monomial<E> divider) {
            DegreeVector dv = dividend.dvDivideOrNull(divider);
            if (dv == null)
                return null;
            E div = ring.divideOrNull(dividend.coefficient, divider.coefficient);
            if (div == null)
                return null;
            return new Monomial<>(dv, div);
        }

        @Override
        public Monomial<E> negate(Monomial<E> term) {
            return term.setCoefficient(ring.negate(term.coefficient));
        }

        @Override
        public boolean isZero(Monomial<E> term) {
            return ring.isZero(term.coefficient);
        }

        @Override
        public boolean isOne(Monomial<E> term) {
            return isConstant(term) && ring.isOne(term.coefficient);
        }

        @Override
        public boolean isUnit(Monomial<E> term) {
            return isConstant(term) && ring.isUnit(term.coefficient);
        }

        @Override
        public Monomial<E> createTermWithUnitCoefficient(int[] exponents) {
            return new Monomial<>(exponents, ArraysUtil.sum(exponents), ring.getOne());
        }

        @Override
        @SuppressWarnings("unchecked")
        public Monomial<E>[] createArray(int length) {
            return new Monomial[length];
        }

        @Override
        public Monomial<E> getUnitTerm(int nVariables) {
            return new Monomial<>(nVariables, ring.getOne());
        }

        @Override
        public Monomial<E> getZeroTerm(int nVariables) {
            return new Monomial<>(nVariables, ring.getZero());
        }
    }
}
