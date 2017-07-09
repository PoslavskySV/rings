package cc.r2.core.poly.multivar;

import cc.r2.core.util.ArraysUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class DegreeVector<MonomialTerm extends DegreeVector> {
    /** exponents */
    final int[] exponents;
    /** sum of all exponents */
    final int totalDegree;

    DegreeVector(int[] exponents, int totalDegree) {
        this.exponents = exponents;
        this.totalDegree = totalDegree;
    }

    DegreeVector(int nVariables, int position, int exponent) {
        this.exponents = new int[nVariables];
        this.exponents[position] = exponent;
        this.totalDegree = exponent;
    }

    DegreeVector(int[] exponents) {
        this(exponents, ArraysUtil.sum(exponents));
    }

    /** internal method */
    abstract MonomialTerm setDegreeVector(int[] newDegree, int newTotalDegree);

    final MonomialTerm setDegreeVector(int[] newDegree) {
        return setDegreeVector(newDegree, ArraysUtil.sum(newDegree));
    }

    /** internal method */
    final MonomialTerm setDegreeVector(DegreeVector dv) {
        return setDegreeVector(dv.exponents, dv.totalDegree);
    }

    /** set i-th exponent to zero and return new Monomial */
    abstract MonomialTerm setZero(int var);

    /** set i-th exponent to zero and return new Monomial */
    abstract MonomialTerm setZero(int[] vars);

    /** Divide degree vector */
    final MonomialTerm divide(DegreeVector oth) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return setDegreeVector(newExponents, totalDegree - oth.totalDegree);
    }

    /** Whether divides degree vector */
    final boolean dividesQ(DegreeVector oth) {
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

    /** Adjoins new variable (with zero exponent) to this monomial */
    public final MonomialTerm joinNewVariable() {
        return setDegreeVector(Arrays.copyOf(exponents, exponents.length + 1), totalDegree);
    }

    public final MonomialTerm joinNewVariables(int newNVariables, int[] mapping) {
        int[] newExponents = new int[newNVariables];
        int c = 0;
        for (int i : mapping)
            newExponents[i] = exponents[c++];

        return setDegreeVector(newExponents, totalDegree);
    }

    /** Set's exponents of all variables except specified ones to zero */
    public final MonomialTerm of(int[] variables) {
        int[] exs = new int[exponents.length];
        int totalDegree = 0;
        for (int i : variables) {
            exs[i] = exponents[i];
            totalDegree += exs[i];
        }
        return setDegreeVector(exs, totalDegree);
    }

    /** Set's exponents of specified variables to zero */
    public final MonomialTerm except(int[] variables) {
        int[] exs = exponents.clone();
        int totalDegree = this.totalDegree;
        for (int i : variables) {
            exs[i] = 0;
            totalDegree -= exponents[i];
        }
        return setDegreeVector(exs, totalDegree);
    }

    public final MonomialTerm without(int variable) {
        return setDegreeVector(ArraysUtil.remove(exponents, variable), totalDegree - exponents[variable]);
    }

    public final MonomialTerm insert(int variable) {
        return setDegreeVector(ArraysUtil.insert(exponents, variable, 0), totalDegree);
    }

    /** Set's exponent of i-th variable to specified value */
    public final MonomialTerm set(int i, int exponent) {
        int[] newExponents = exponents.clone();
        newExponents[i] = exponent;
        return setDegreeVector(newExponents, totalDegree - exponents[i] + exponent);
    }

    private static String toString0(String var, int exp) {
        return exp == 0 ? "" : var + (exp == 1 ? "" : "^" + exp);
    }

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

    /**
     * Lexicographic monomial order
     */
    public static final Comparator<DegreeVector> LEX = (DegreeVector a, DegreeVector b) -> {
        for (int i = 0; i < a.exponents.length; ++i) {
            int c = Integer.compare(a.exponents[i], b.exponents[i]);
            if (c != 0)
                return c;
        }
        return 0;
    };

    /**
     * Antilexicographic monomial order
     */
    public static final Comparator<DegreeVector> ALEX = (DegreeVector a, DegreeVector b) -> LEX.compare(b, a);

    /**
     * Graded lexicographic monomial order
     */
    public static final Comparator<DegreeVector> GRLEX = (DegreeVector a, DegreeVector b) -> {
        int c = Integer.compare(a.totalDegree, b.totalDegree);
        return c != 0 ? c : LEX.compare(a, b);
    };

    /**
     * Graded reverse lexicographic monomial order
     */
    public static final Comparator<DegreeVector> GREVLEX = (DegreeVector a, DegreeVector b) -> {
        int c = Integer.compare(a.totalDegree, b.totalDegree);
        if (c != 0)
            return c;
        for (int i = a.exponents.length - 1; i >= 0; --i) {
            c = Integer.compare(b.exponents[i], a.exponents[i]);
            if (c != 0)
                return c;
        }
        return 0;
    };
}