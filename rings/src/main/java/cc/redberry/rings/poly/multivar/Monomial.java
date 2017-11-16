package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.util.ArraysUtil;

/**
 * Monomial with coefficient from generic ring.
 *
 * @since 1.0
 */
public final class Monomial<E> extends DegreeVector<Monomial<E>> {
    private static final long serialVersionUID = 1L;
    /** The coefficient */
    public final E coefficient;

    /**
     * Creates monomial with specified degree vector and coefficient
     *
     * @param exponents   the degree vector
     * @param totalDegree sum of exponents
     * @param coefficient the coefficient
     */
    public Monomial(int[] exponents, int totalDegree, E coefficient) {
        super(exponents, totalDegree);
        this.coefficient = coefficient;
    }

    /**
     * Creates monomial with specified degree vector and coefficient
     *
     * @param exponents   the degree vector
     * @param coefficient the coefficient
     */
    public Monomial(int[] exponents, E coefficient) {
        this(exponents, ArraysUtil.sum(exponents), coefficient);
    }

    /**
     * Creates monomial with specified number of variables and one single variable with non zero exponent
     *
     * @param nVariables  the number of variables
     * @param variable    the only one variable with non zero exponent
     * @param exponent    the exponent
     * @param coefficient the coefficient
     */
    public Monomial(int nVariables, int variable, int exponent, E coefficient) {
        super(nVariables, variable, exponent);
        this.coefficient = coefficient;
    }

    @Override
    public Monomial<E> setDegreeVector(int[] newDegree, int newTotalDegree) {
        return new Monomial<>(newDegree, newTotalDegree, coefficient);
    }

    @Override
    public Monomial<E> setCoefficientFrom(Monomial<E> oth) {
        return new Monomial<>(exponents, totalDegree, oth.coefficient);
    }

    /** Set's the coefficient to {@code newRing.valueOf(coefficient) } */
    public Monomial<E> setRing(Ring<E> newRing) {
        E e = newRing.valueOf(coefficient);
        return coefficient == e ? this : new Monomial<>(exponents, totalDegree, e);
    }

    /** Set's monomial coefficient to a specified value */
    public Monomial<E> setCoefficient(E value) {
        return coefficient == value ? this : new Monomial<>(exponents, totalDegree, value);
    }

    /** Negates the coefficient */
    public Monomial<E> negate(Ring<E> ring) {
        return setCoefficient(ring.negate(coefficient));
    }

    /** Multiplies this by {@code oth} and sets the resulting coefficient to the specified value */
    public Monomial<E> multiply(DegreeVector oth, E coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++)
            newExponents[i] = exponents[i] + oth.exponents[i];
        return new Monomial<>(newExponents, totalDegree + oth.totalDegree, coefficient);
    }

    /**
     * Divides this by {@code oth} and sets the resulting coefficient to the specified value or returns null if exact
     * division is not possible
     */
    public Monomial<E> divide(DegreeVector oth, E coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return new Monomial<>(newExponents, totalDegree - oth.totalDegree, coefficient);
    }

    /**
     * Removes specified variable from degree vector  (number of variables will be reduced)
     *
     * @param variable    the variable
     * @param coefficient the new coefficient
     */
    public Monomial<E> without(int variable, E coefficient) {
        if (exponents.length == 1) {
            assert variable == 0;
            return new Monomial<>(new int[0], 0, coefficient);
        }
        return new Monomial<>(ArraysUtil.remove(exponents, variable), totalDegree - exponents[variable], coefficient);
    }

    /**
     * Set exponent of specified variable to zero and the coefficient to a new value
     *
     * @param variable    the variable
     * @param coefficient the new coefficient
     */
    public Monomial<E> setZero(int variable, E coefficient) {
        if (exponents.length == 1) {
            assert variable == 0;
            return new Monomial<>(new int[1], 0, coefficient);
        }
        int[] newExponents = exponents.clone();
        newExponents[variable] = 0;
        return new Monomial<>(newExponents, totalDegree - exponents[variable], coefficient);
    }

    /**
     * Set exponents of specified variables to zero and the coefficient to a new value
     *
     * @param variables   the variables
     * @param coefficient the new coefficient
     */
    public Monomial<E> setZero(int[] variables, E coefficient) {
        if (variables.length == 0)
            return setCoefficient(coefficient);
        int[] newExponents = exponents.clone();
        int totalDeg = totalDegree;
        for (int var : variables) {
            totalDeg -= exponents[var];
            newExponents[var] = 0;
        }
        return new Monomial<>(newExponents, totalDeg, coefficient);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        Monomial<?> that = (Monomial<?>) o;

        return coefficient != null ? coefficient.equals(that.coefficient) : that.coefficient == null;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + (coefficient != null ? coefficient.hashCode() : 0);
        return result;
    }

    /**
     * Creates monomial with all zero exponents and specified coefficient (constant monomial)
     *
     * @param nVariables  the number of variables
     * @param coefficient the coefficient
     * @return constant monomial with specified coefficient
     */
    public static <E> Monomial<E> withZeroExponents(int nVariables, E coefficient) {
        int[] exponents = nVariables < zeroDegreeVectors.length ? zeroDegreeVectors[nVariables] : new int[nVariables];
        return new Monomial<>(exponents, 0, coefficient);
    }

    static <E> Monomial<E> create(E coeff, int... exponents) {
        return new Monomial<>(exponents, coeff);
    }
}
