package cc.r2.core.poly.multivar2;

import cc.r2.core.util.ArraysUtil;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lMonomialTerm extends DegreeVector<lMonomialTerm> {
    final long coefficient;

    public lMonomialTerm(int[] exponents, int totalDegree, long coefficient) {
        super(exponents, totalDegree);
        this.coefficient = coefficient;
    }

    public lMonomialTerm(int nVariables, int position, int exponent, long coefficient) {
        super(nVariables, position, exponent);
        this.coefficient = coefficient;
    }

    public lMonomialTerm(int[] exponents, long coefficient) {
        super(exponents);
        this.coefficient = coefficient;
    }

    @Override
    protected lMonomialTerm setDegreeVector(int[] newDegree, int newTotalDegree) {
        return new lMonomialTerm(newDegree, newTotalDegree, coefficient);
    }

    /** Set's monomial coefficient to a specified value */
    public lMonomialTerm setCoefficient(long value) {
        return new lMonomialTerm(exponents, totalDegree, value);
    }

    /** Multiplies this by {@code oth} and sets the resulting coefficient to a specified value */
    public lMonomialTerm multiply(DegreeVector oth, long coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++)
            newExponents[i] = exponents[i] + oth.exponents[i];
        return new lMonomialTerm(newExponents, totalDegree + oth.totalDegree, coefficient);
    }

    /**
     * Divides this by {@code oth} and sets the resulting coefficient to a specified value or returns
     * null if exact division is not possible
     */
    public lMonomialTerm divide(DegreeVector oth, long coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return new lMonomialTerm(newExponents, totalDegree - oth.totalDegree, coefficient);
    }

    /** removes i-th variable from degree vector */
    lMonomialTerm without(int i) {
        return without(i, coefficient);
    }

    /** removes i-th variable from degree vector */
    lMonomialTerm without(int i, long coefficient) {
        if (exponents.length == 1) {
            assert i == 0;
            return new lMonomialTerm(new int[0], 0, coefficient);
        }
        return new lMonomialTerm(ArraysUtil.remove(exponents, i), totalDegree - exponents[i], coefficient);
    }

    @Override
    lMonomialTerm setZero(int i) {
        return setZero(i, coefficient);
    }

    /** set i-th exponent to zero and the coefficient to a new value */
    lMonomialTerm setZero(int i, long coefficient) {
        if (exponents.length == 1) {
            assert i == 0;
            return new lMonomialTerm(new int[0], 0, coefficient);
        }
        int[] newExponents = exponents.clone();
        newExponents[i] = 0;
        return new lMonomialTerm(newExponents, totalDegree - exponents[i], coefficient);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;

        lMonomialTerm that = (lMonomialTerm) o;

        return coefficient == that.coefficient;
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + Long.hashCode(coefficient);
        return result;
    }

    public static lMonomialTerm withZeroExponents(int nVariables, long coefficient) {
        int[] exponents = nVariables < zeroDegreeVectors.length ? zeroDegreeVectors[nVariables] : new int[nVariables];
        return new lMonomialTerm(exponents, 0, coefficient);
    }
}
