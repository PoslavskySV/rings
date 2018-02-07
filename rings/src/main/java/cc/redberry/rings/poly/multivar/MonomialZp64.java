package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.bigint.BigInteger;

/**
 * Monomial with coefficient from Zp with p < 2^64
 *
 * @see cc.redberry.rings.poly.multivar.IMonomialAlgebra.MonomialAlgebraZp64
 * @since 1.0
 */
public final class MonomialZp64 extends AMonomial<MonomialZp64> {
    /** the coefficient */
    public final long coefficient;

    /**
     * @param degreeVector DegreeVector
     * @param coefficient  the coefficient
     */
    public MonomialZp64(DegreeVector degreeVector, long coefficient) {
        super(degreeVector);
        this.coefficient = coefficient;
    }

    /**
     * @param exponents   exponents
     * @param totalDegree total degree (sum of exponents)
     * @param coefficient the coefficient
     */
    public MonomialZp64(int[] exponents, int totalDegree, long coefficient) {
        super(exponents, totalDegree);
        this.coefficient = coefficient;
    }

    /**
     * @param exponents   exponents
     * @param coefficient the coefficient
     */
    public MonomialZp64(int[] exponents, long coefficient) {
        super(exponents);
        this.coefficient = coefficient;
    }

    public MonomialZp64(int nVariables, long coefficient) {
        this(new int[nVariables], 0, coefficient);
    }

    @Override
    public MonomialZp64 setCoefficientFrom(MonomialZp64 oth) {
        if (coefficient == oth.coefficient)
            return this;
        return new MonomialZp64(this, oth.coefficient);
    }

    @Override
    public MonomialZp64 setDegreeVector(DegreeVector oth) {
        if (this == oth)
            return this;
        if (oth == null)
            return null;
        if (oth.exponents == exponents)
            return this;
        return new MonomialZp64(oth, coefficient);
    }

    @Override
    public MonomialZp64 setDegreeVector(int[] exponents, int totalDegree) {
        if (this.exponents == exponents)
            return this;
        return new MonomialZp64(exponents, totalDegree, coefficient);
    }

    public MonomialZp64 setCoefficient(long c) {
        if (coefficient == c)
            return this;
        return new MonomialZp64(exponents, totalDegree, c);
    }

    public Monomial<BigInteger> toBigMonomial() {
        return new Monomial<>(this, BigInteger.valueOf(coefficient));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        MonomialZp64 that = (MonomialZp64) o;

        return coefficient == that.coefficient;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (int) (coefficient ^ (coefficient >>> 32));
        return result;
    }

    @Override
    public String toString() {
        String dvString = super.toString();
        String cfString = Long.toString(coefficient);
        if (dvString.isEmpty())
            return cfString;
        if (coefficient == 1)
            return dvString;
        return coefficient + "*" + dvString;
    }
}
