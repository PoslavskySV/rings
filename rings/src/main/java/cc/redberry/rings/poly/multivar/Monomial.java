package cc.redberry.rings.poly.multivar;

/**
 * Monomial with coefficient from generic ring
 *
 * @see cc.redberry.rings.poly.multivar.IMonomialAlgebra.MonomialAlgebra
 * @since 1.0
 */
public class Monomial<E> extends AMonomial<Monomial<E>> {
    /** the coefficient */
    public final E coefficient;

    /**
     * @param degreeVector DegreeVector
     * @param coefficient  the coefficient
     */
    public Monomial(DegreeVector degreeVector, E coefficient) {
        super(degreeVector);
        this.coefficient = coefficient;
    }

    /**
     * @param exponents   exponents
     * @param totalDegree total degree (sum of exponents)
     * @param coefficient the coefficient
     */
    public Monomial(int[] exponents, int totalDegree, E coefficient) {
        super(exponents, totalDegree);
        this.coefficient = coefficient;
    }

    /**
     * @param exponents   exponents
     * @param coefficient the coefficient
     */
    public Monomial(int[] exponents, E coefficient) {
        super(exponents);
        this.coefficient = coefficient;
    }

    public Monomial(int nVariables, E coefficient) {
        this(new int[nVariables], 0, coefficient);
    }

    @Override
    public Monomial<E> setCoefficientFrom(Monomial<E> oth) {
        return new Monomial<>(this, oth.coefficient);
    }

    @Override
    public Monomial<E> setDegreeVector(DegreeVector oth) {
        if (this == oth)
            return this;
        if (oth == null)
            return null;
        if (oth.exponents == exponents)
            return this;
        return new Monomial<>(oth, coefficient);
    }

    @Override
    public Monomial<E> setDegreeVector(int[] exponents, int totalDegree) {
        if (this.exponents == exponents)
            return this;
        return new Monomial<>(exponents, totalDegree, coefficient);
    }

    @Override
    public Monomial<E> forceSetDegreeVector(int[] exponents, int totalDegree) {
        return new Monomial<>(exponents, totalDegree, coefficient);
    }

    public Monomial<E> setCoefficient(E c) {
        return new Monomial<>(exponents, totalDegree, c);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        Monomial<?> monomial = (Monomial<?>) o;

        return coefficient.equals(monomial.coefficient);
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + coefficient.hashCode();
        return result;
    }

    @Override
    public String toString() {
        String dvString = super.toString();
        String cfString = coefficient.toString();
        if (dvString.isEmpty())
            return cfString;
        if (cfString.equals("1"))
            return dvString;
        return coefficient + "*" + dvString;
    }
}
