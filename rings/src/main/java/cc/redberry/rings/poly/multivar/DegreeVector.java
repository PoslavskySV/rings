package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.io.IStringifier;
import cc.redberry.rings.util.ArraysUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Degree vector. This is parent class for all monomials. Instances are immutable. All {@code DegreeVector} methods are
 * prefixed with "dv" (which expands to "degree vector"), which means that they affect only exponents (not the
 * coefficients).
 *
 * @see AMonomial
 * @since 1.0
 */
public class DegreeVector implements java.io.Serializable {
    private static final long serialVersionUID = 1L;
    /** exponents */
    public final int[] exponents;
    /** Sum of all exponents (total degree) */
    public final int totalDegree;

    /**
     * @param exponents   exponents
     * @param totalDegree total degree (sum of exponents)
     */
    public DegreeVector(int[] exponents, int totalDegree) {
        this.exponents = exponents;
        this.totalDegree = totalDegree;
        // assert ArraysUtil.sum(exponents) == totalDegree;
    }

    /**
     * @param exponents exponents
     */
    public DegreeVector(int[] exponents) {
        this(exponents, ArraysUtil.sum(exponents));
    }

    /** Returns number of variables */
    public final int nVariables() { return exponents.length;}

    /** Returns whether all exponents are zero */
    public final boolean isZeroVector() {
        return totalDegree == 0;
    }

    /** Returns the total degree in specified variables */
    public final int dvTotalDegree(int... variables) {
        int d = 0;
        for (int v : variables)
            d += exponents[v];
        return d;
    }

    /** Multiplies this by oth */
    public final DegreeVector dvMultiply(DegreeVector oth) {
        if (oth.isZeroVector())
            return this;
        int[] res = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++)
            res[i] = exponents[i] + oth.exponents[i];
        return new DegreeVector(res, totalDegree + oth.totalDegree);
    }

    /** Multiplies this by oth */
    public final DegreeVector dvMultiply(int[] oth) {
        int deg = totalDegree;
        int[] res = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            res[i] = exponents[i] + oth[i];
            deg += oth[i];
        }
        if (deg == 0)
            return this; // avoid copying
        return new DegreeVector(res, deg);
    }

    /** Multiplies this by variable^exponent */
    public final DegreeVector dvMultiply(int variable, int exponent) {
        int[] res = exponents.clone();
        res[variable] += exponent;
        if (res[variable] < 0)
            return null;
        return new DegreeVector(res, totalDegree + exponent);
    }

    /** Divides this by variable^exponent */
    public final DegreeVector dvDivideOrNull(int variable, int exponent) {
        return dvMultiply(variable, -exponent);
    }

    /** Gives quotient {@code this / oth } or null if exact division is not possible (e.g. a^2*b^3 / a^3*b^5) */
    public final DegreeVector dvDivideOrNull(DegreeVector divider) {
        if (divider.isZeroVector())
            return this;
        int[] res = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            res[i] = exponents[i] - divider.exponents[i];
            if (res[i] < 0)
                return null;
        }
        return new DegreeVector(res, totalDegree - divider.totalDegree);
    }

    /** Gives quotient {@code this / oth } or null if exact division is not possible (e.g. a^2*b^3 / a^3*b^5) */
    public final DegreeVector dvDivideOrNull(int[] divider) {
        int deg = totalDegree;
        int[] res = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            res[i] = exponents[i] - divider[i];
            if (res[i] < 0)
                return null;
            deg -= divider[i];
        }
        if (deg == 0)
            return this; // avoid copying
        return new DegreeVector(res, deg);
    }

    /**
     * Gives quotient {@code this / oth } or throws {@code ArithmeticException} if exact division is not possible (e.g.
     * a^2*b^3 / a^3*b^5)
     */
    public final DegreeVector dvDivideExact(DegreeVector divider) {
        DegreeVector quot = dvDivideOrNull(divider);
        if (quot == null)
            throw new ArithmeticException("not divisible");
        return quot;
    }

    /**
     * Gives quotient {@code this / oth } or throws {@code ArithmeticException} if exact division is not possible (e.g.
     * a^2*b^3 / a^3*b^5)
     */
    public final DegreeVector dvDivideExact(int[] divider) {
        DegreeVector quot = dvDivideOrNull(divider);
        if (quot == null)
            throw new ArithmeticException("not divisible");
        return quot;
    }

    /** Tests whether this can be divided by {@code oth} degree vector */
    public final boolean dvDivisibleBy(int[] oth) {
        for (int i = 0; i < exponents.length; i++)
            if (exponents[i] < oth[i])
                return false;
        return true;
    }

    /** Tests whether this can be divided by {@code oth} degree vector */
    public final boolean dvDivisibleBy(DegreeVector oth) {
        return dvDivisibleBy(oth.exponents);
    }

    /** Joins new variable (with zero exponent) to degree vector */
    public final DegreeVector dvJoinNewVariable() {
        return dvJoinNewVariables(1);
    }

    /** Joins new variables (with zero exponents) to degree vector */
    public final DegreeVector dvJoinNewVariables(int n) {
        return new DegreeVector(Arrays.copyOf(exponents, exponents.length + n), totalDegree);
    }

    /** internal API */
    public final DegreeVector dvJoinNewVariables(int newNVariables, int[] mapping) {
        int[] res = new int[newNVariables];
        int c = 0;
        for (int i : mapping)
            res[i] = exponents[c++];
        return new DegreeVector(res, totalDegree);
    }

    /** Sets the number of variables */
    public final DegreeVector dvSetNVariables(int n) {
        if (n == exponents.length)
            return this;
        if (n > exponents.length)
            return new DegreeVector(Arrays.copyOf(exponents, n), totalDegree);
        else
            return new DegreeVector(Arrays.copyOf(exponents, n));
    }

    /** Sets exponents of all variables except the specified variable to zero */
    public final DegreeVector dvSelect(int var) {
        int[] res = new int[exponents.length];
        res[var] = exponents[var];
        return new DegreeVector(res, exponents[var]);
    }

    /** Set's exponents of all variables except specified variables to zero */
    public final DegreeVector dvSelect(int[] variables) {
        int[] res = new int[exponents.length];
        int deg = 0;
        for (int i : variables) {
            res[i] = exponents[i];
            deg += exponents[i];
        }
        return new DegreeVector(res, deg);
    }

    /** Picks only specified exponents */
    public final DegreeVector dvDropSelect(int[] variables) {
        int[] res = new int[variables.length];
        int deg = 0;
        int c = 0;
        for (int i : variables) {
            res[c++] = exponents[i];
            deg += exponents[i];
        }
        return new DegreeVector(res, deg);
    }

    /**
     * Selects range from this
     *
     * @param from from inclusive
     * @param to   to exclusive
     */
    public final DegreeVector dvRange(int from, int to) {
        if (from == 0 && to == exponents.length) {
            return this;
        }
        return new DegreeVector(Arrays.copyOfRange(exponents, from, to));
    }

    /** Set exponent of specified {@code var} to zero */
    public final DegreeVector dvSetZero(int var) {
        int[] res = exponents.clone();
        res[var] = 0;
        return new DegreeVector(res, totalDegree - exponents[var]);
    }

    /** Set exponents of specified variables to zero */
    public final DegreeVector dvSetZero(int[] variables) {
        int[] res = exponents.clone();
        int deg = totalDegree;
        for (int i : variables) {
            deg -= exponents[i];
            res[i] = 0;
        }
        return new DegreeVector(res, deg);
    }

    /** Drops specified variable (number of variables will be reduced) */
    public final DegreeVector dvWithout(int variable) {
        return new DegreeVector(ArraysUtil.remove(exponents, variable), totalDegree - exponents[variable]);
    }

    /** Drops specified variables (number of variables will be reduced) */
    public final DegreeVector dvWithout(int[] variables) {
        return new DegreeVector(ArraysUtil.remove(exponents, variables));
    }

    /** Inserts new variable */
    public final DegreeVector dvInsert(int variable) {
        return new DegreeVector(ArraysUtil.insert(exponents, variable, 0), totalDegree);
    }

    /** Inserts new variables */
    public final DegreeVector dvInsert(int variable, int count) {
        return new DegreeVector(ArraysUtil.insert(exponents, variable, 0, count), totalDegree);
    }

    /**
     * Set's exponent of specified variable to specified value
     *
     * @param variable the variable
     * @param exponent new exponent
     */
    public final DegreeVector dvSet(int variable, int exponent) {
        if (exponents[variable] == exponent)
            return this;
        int deg = totalDegree - exponents[variable] + exponent;
        int[] res = exponents.clone();
        res[variable] = exponent;
        return new DegreeVector(res, deg);
    }

    final int firstNonZeroVariable() {
        for (int i = 0; i < exponents.length; ++i)
            if (exponents[i] != 0)
                return i;
        return -1;
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
    public String toString() {
        return toString(IStringifier.defaultVars(exponents.length));
    }

    public final String toStringArray() {
        return Arrays.toString(exponents);
    }

    public final boolean dvEquals(DegreeVector dVector) {
        return totalDegree == dVector.totalDegree && Arrays.equals(exponents, dVector.exponents);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        DegreeVector dVector = (DegreeVector) o;
        return dvEquals(dVector);
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(exponents);
        result = 31 * result + totalDegree;
        return result;
    }

}