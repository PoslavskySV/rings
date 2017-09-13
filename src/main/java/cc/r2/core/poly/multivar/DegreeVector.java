package cc.r2.core.poly.multivar;

import cc.r2.core.util.ArraysUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Parent class for monomials which holds the degree vector of monomial.
 * Instances of this class are immutable (all structural operations return new instances).
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class DegreeVector<MonomialTerm extends DegreeVector<MonomialTerm>> implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    /** exponents */
    final int[] exponents;
    /** Sum of all exponents (total degree) */
    public final int totalDegree;

    protected DegreeVector(int[] exponents, int totalDegree) {
        this.exponents = exponents;
        this.totalDegree = totalDegree;
        //assertions();
    }

    protected DegreeVector(int nVariables, int position, int exponent) {
        this.exponents = new int[nVariables];
        this.exponents[position] = exponent;
        this.totalDegree = exponent;
        //assertions();
    }

    protected DegreeVector(int[] exponents) {
        this(exponents, ArraysUtil.sum(exponents));
    }

    private void assertions() {
        assert ArraysUtil.sum(exponents) == totalDegree;
    }

    /**
     * Set the degree vector to a new value
     *
     * @param newDegree      new degree vector
     * @param newTotalDegree sum of newDegree
     */
    public abstract MonomialTerm setDegreeVector(int[] newDegree, int newTotalDegree);

    /**
     * Set the degree vector to a new value
     *
     * @param newDegree new degree vector
     * @return new monomial
     */
    public final MonomialTerm setDegreeVector(int[] newDegree) {
        return setDegreeVector(newDegree, ArraysUtil.sum(newDegree));
    }

    /**
     * Set the degree vector to a new value
     *
     * @param dv new degree vector
     * @return new monomial
     */
    public final MonomialTerm setDegreeVector(DegreeVector dv) {
        return setDegreeVector(dv.exponents, dv.totalDegree);
    }

    /**
     * Returns  exponents of all variables except the specified variable to zero and exponent of specified variable to
     * one (so the result is just "var")
     *
     * @param var the variable
     * @return the new monomial
     */
    public final MonomialTerm singleVar(int var) {
        int[] newExponents = new int[exponents.length];
        newExponents[var] = 1;
        return setDegreeVector(newExponents, 1);
    }

    /**
     * Divides this monomial by {@code oth} not taking into account the coefficient of {@code oth} (division of
     * degree vectors), returns null if exact division is not possible (e.g. a^2*b^3 / a^3*b^5)
     *
     * @param oth the degree vector
     * @return {@code this / oth} or null if exact division is not possible
     */
    public final MonomialTerm divideOrNull(DegreeVector oth) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return setDegreeVector(newExponents, totalDegree - oth.totalDegree);
    }

    /** Tests whether this can be divided by {@code oth} degree vector */
    public final boolean dividesQ(DegreeVector oth) {
        for (int i = 0; i < exponents.length; i++) {
            if (exponents[i] < oth.exponents[i])
                return false;
        }
        return true;
    }

    /** Returns whether all exponents are zero */
    public final boolean isZeroVector() {
        return totalDegree == 0;
    }

    /** Joins new variable (with zero exponent) to this monomial */
    public final MonomialTerm joinNewVariable() {
        return setDegreeVector(Arrays.copyOf(exponents, exponents.length + 1), totalDegree);
    }

    /** internal API */
    final MonomialTerm joinNewVariables(int newNVariables, int[] mapping) {
        int[] newExponents = new int[newNVariables];
        int c = 0;
        for (int i : mapping)
            newExponents[i] = exponents[c++];

        return setDegreeVector(newExponents, totalDegree);
    }

    /**
     * Sets exponents of all variables except the specified variable to zero
     *
     * @param var the variable
     * @return new monomial
     */
    public final MonomialTerm select(int var) {
        int[] newExponents = new int[exponents.length];
        newExponents[var] = exponents[var];
        return setDegreeVector(newExponents, exponents[var]);
    }

    /**
     * Set's exponents of all variables except specified variables to zero
     *
     * @param variables the variables
     * @return new monomial
     */
    public final MonomialTerm select(int[] variables) {
        int[] exs = new int[exponents.length];
        int totalDegree = 0;
        for (int i : variables) {
            exs[i] = exponents[i];
            totalDegree += exs[i];
        }
        return setDegreeVector(exs, totalDegree);
    }

    /**
     * Set exponent of specified {@code var} to zero
     *
     * @param var the variable
     * @return new monomial
     */
    public final MonomialTerm setZero(int var) {
        int[] exs = exponents.clone();
        int totalDegree = this.totalDegree;
        exs[var] = 0;
        totalDegree -= exponents[var];
        return setDegreeVector(exs, totalDegree);
    }

    /**
     * Set exponents of specified variables to zero and return new Monomial
     *
     * @param variables the array of variables
     * @return new monomial
     */
    public final MonomialTerm setZero(int[] variables) {
        int[] exs = exponents.clone();
        int totalDegree = this.totalDegree;
        for (int i : variables) {
            exs[i] = 0;
            totalDegree -= exponents[i];
        }
        return setDegreeVector(exs, totalDegree);
    }

    /**
     * Drops specified variable (number of variables will be reduced)
     *
     * @param variable the variable
     * @return new monomial
     */
    public final MonomialTerm without(int variable) {
        return setDegreeVector(ArraysUtil.remove(exponents, variable), totalDegree - exponents[variable]);
    }

    /**
     * Inserts new variable
     *
     * @param variable the variable
     * @return new monomial
     */
    public final MonomialTerm insert(int variable) {
        return setDegreeVector(ArraysUtil.insert(exponents, variable, 0), totalDegree);
    }

    /**
     * Set's exponent of specified variable to specified value
     *
     * @param variable the variable
     * @param exponent new exponent
     * @return new monomial
     */
    public final MonomialTerm set(int variable, int exponent) {
        int[] newExponents = exponents.clone();
        newExponents[variable] = exponent;
        return setDegreeVector(newExponents, totalDegree - exponents[variable] + exponent);
    }

    private static String toString0(String var, int exp) {
        return exp == 0 ? "" : var + (exp == 1 ? "" : "^" + exp);
    }

    /**
     * String representation of this monomial with specified string names for variables
     *
     * @param vars string names of variables
     */
    public final String toString(String[] vars) {
        List<String> result = new ArrayList<>();
        for (int i = 0; i < exponents.length; i++)
            result.add(toString0(vars[i], exponents[i]));
        return result.stream().filter(s -> !s.isEmpty()).collect(Collectors.joining("*"));
    }

    @Override
    public final String toString() {
        return toString(MultivariatePolynomial.defaultVars(exponents.length));
    }

    public final String toStringArray() {
        return Arrays.toString(exponents);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MonomialTerm that = (MonomialTerm) o;

        if (totalDegree != that.totalDegree) return false;
        if (!Arrays.equals(exponents, that.exponents)) return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(exponents);
        result = 31 * result + totalDegree;
        return result;
    }

    /** cached zero degree vectors */
    static int[][] zeroDegreeVectors = new int[32][];

    static {
        for (int i = 0; i < zeroDegreeVectors.length; i++)
            zeroDegreeVectors[i] = new int[i];
    }

}