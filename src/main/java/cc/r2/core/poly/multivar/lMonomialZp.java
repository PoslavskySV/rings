package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.util.ArraysUtil;

/**
 * Monomial with machine integer coefficient.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lMonomialZp extends DegreeVector<lMonomialZp> {
    private static final long serialVersionUID = 1L;
    /** The coefficient */
    public final long coefficient;

    /**
     * Creates monomial with specified degree vector and coefficient
     *
     * @param exponents   the degree vector
     * @param totalDegree sum of exponents
     * @param coefficient the coefficient
     */
    public lMonomialZp(int[] exponents, int totalDegree, long coefficient) {
        super(exponents, totalDegree);
        this.coefficient = coefficient;
    }

    /**
     * Creates monomial with specified degree vector and coefficient
     *
     * @param exponents   the degree vector
     * @param coefficient the coefficient
     */
    public lMonomialZp(int[] exponents, long coefficient) {
        super(exponents);
        this.coefficient = coefficient;
    }

    /**
     * Creates monomial with specified number of variables and one single variable with non zero exponent
     *
     * @param nVariables  the number of variables
     * @param variable    the only one variable with non zero exponent
     * @param exponent    the exponent
     * @param coefficient the coefficient
     */
    public lMonomialZp(int nVariables, int variable, int exponent, long coefficient) {
        super(nVariables, variable, exponent);
        this.coefficient = coefficient;
    }

    /** Reduces the coefficient modulo new modulus */
    public lMonomialZp setDomain(lIntegersModulo newDomain) {
        long e = newDomain.modulus(coefficient);
        return coefficient == e ? this : new lMonomialZp(exponents, totalDegree, e);
    }

    /** Set's the coefficient to {@code newDomain.valueOf(coefficient) } */
    public <E> Monomial<E> setDomain(Domain<E> newDomain) {
        E e = newDomain.valueOf(coefficient);
        return new Monomial<>(exponents, totalDegree, e);
    }

    @Override
    public lMonomialZp setDegreeVector(int[] newDegree, int newTotalDegree) {
        return new lMonomialZp(newDegree, newTotalDegree, coefficient);
    }

    /** Set's monomial coefficient to a specified value */
    public lMonomialZp setCoefficient(long value) {
        return value == coefficient ? this : new lMonomialZp(exponents, totalDegree, value);
    }

    /** Negates the coefficient */
    public lMonomialZp negate(lIntegersModulo domain) {
        return setCoefficient(domain.negate(coefficient));
    }

    /** Multiplies this monomial by {@code oth} and sets the resulting coefficient to a specified value */
    public lMonomialZp multiply(DegreeVector oth, long coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++)
            newExponents[i] = exponents[i] + oth.exponents[i];
        return new lMonomialZp(newExponents, totalDegree + oth.totalDegree, coefficient);
    }

    /**
     * Divides this by {@code oth} and sets the resulting coefficient to a specified value or returns
     * null if exact division is not possible
     */
    public lMonomialZp divide(DegreeVector oth, long coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return new lMonomialZp(newExponents, totalDegree - oth.totalDegree, coefficient);
    }

    /**
     * Removes specified variable from degree vector  (number of variables will be reduced)
     *
     * @param variable    the variable
     * @param coefficient the new coefficient
     */
    public lMonomialZp without(int variable, long coefficient) {
        if (exponents.length == 1) {
            assert variable == 0;
            return new lMonomialZp(new int[0], 0, coefficient);
        }
        return new lMonomialZp(ArraysUtil.remove(exponents, variable), totalDegree - exponents[variable], coefficient);
    }

    /**
     * Set exponent of specified variable to zero and the coefficient to a new value
     *
     * @param variable    the variable
     * @param coefficient the new coefficient
     */
    public lMonomialZp setZero(int variable, long coefficient) {
        if (exponents.length == 1) {
            assert variable == 0;
            return new lMonomialZp(new int[1], 0, coefficient);
        }
        int[] newExponents = exponents.clone();
        newExponents[variable] = 0;
        return new lMonomialZp(newExponents, totalDegree - exponents[variable], coefficient);
    }

    /**
     * Set exponents of specified variables to zero and the coefficient to a new value
     *
     * @param variables   the variables
     * @param coefficient the new coefficient
     */
    public lMonomialZp setZero(int[] variables, long coefficient) {
        if (variables.length == 0)
            return setCoefficient(coefficient);
        int[] newExponents = exponents.clone();
        int totalDeg = totalDegree;
        for (int var : variables) {
            totalDeg -= exponents[var];
            newExponents[var] = 0;
        }
        return new lMonomialZp(newExponents, totalDeg, coefficient);
    }

    public Monomial<BigInteger> toBigMonomial() {
        return new Monomial<>(exponents, totalDegree, BigInteger.valueOf(coefficient));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        lMonomialZp that = (lMonomialZp) o;

        return coefficient == that.coefficient;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + Long.hashCode(coefficient);
        return result;
    }

    /**
     * Creates monomial with all zero exponents and specified coefficient (constant monomial)
     *
     * @param nVariables  the number of variables
     * @param coefficient the coefficient
     * @return constant monomial with specified coefficient
     */
    public static lMonomialZp withZeroExponents(int nVariables, long coefficient) {
        int[] exponents = nVariables < zeroDegreeVectors.length ? zeroDegreeVectors[nVariables] : new int[nVariables];
        return new lMonomialZp(exponents, 0, coefficient);
    }
}
