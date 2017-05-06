package cc.r2.core.poly.multivar2;

import cc.r2.core.poly.generics.Domain;
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
public final class MonomialTerm<E> implements Cloneable {
    /** exponents */
    final int[] exponents;
    /** sum of all exponents */
    final int totalDegree;
    /** the factor */
    final E coefficient;

    MonomialTerm(int[] exponents, int totalDegree, E coefficient) {
        this.exponents = exponents;
        this.totalDegree = totalDegree;
        this.coefficient = coefficient;
    }

    MonomialTerm(int nVariables, int position, int exponent, E coefficient) {
        this.exponents = new int[nVariables];
        this.exponents[position] = exponent;
        this.totalDegree = exponent;
        this.coefficient = coefficient;
    }

    MonomialTerm(int[] exponents, E coefficient) {
        this(exponents, ArraysUtil.sum(exponents), coefficient);
    }

    MonomialTerm<E> setDomain(Domain<E> newDomain) {
        E e = newDomain.valueOf(coefficient);
        return coefficient == e ? this : new MonomialTerm<>(exponents, totalDegree, e);
    }

    public MonomialTerm<E> setCoefficient(E value) {
        return new MonomialTerm<>(exponents, totalDegree, value);
    }

    public MonomialTerm<E> joinNewVariable() {
        return new MonomialTerm<>(Arrays.copyOf(exponents, exponents.length + 1), totalDegree, coefficient);
    }


    public MonomialTerm<E> of(int[] variables) {
        int[] exs = new int[exponents.length];
        int totalDegree = 0;
        for (int i : variables) {
            exs[i] = exponents[i];
            totalDegree += exs[i];
        }
        return new MonomialTerm<E>(exs, totalDegree, coefficient);
    }

    public MonomialTerm<E> except(int[] variables) {
        int[] exs = exponents.clone();
        int totalDegree = this.totalDegree;
        for (int i : variables) {
            exs[i] = 0;
            totalDegree -= exponents[i];
        }
        return new MonomialTerm<E>(exs, totalDegree, coefficient);
    }

    public boolean isZeroVector() {
        return totalDegree == 0;
    }

    private static String toString0(String var, int exp) {
        return exp == 0 ? "" : var + (exp == 1 ? "" : "^" + exp);
    }

    public MonomialTerm<E> multiply(MonomialTerm<E> oth, E coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++)
            newExponents[i] = exponents[i] + oth.exponents[i];
        return new MonomialTerm<E>(newExponents, totalDegree + oth.totalDegree, coefficient);
    }

    public MonomialTerm<E> divide(MonomialTerm<E> oth, E coefficient) {
        int[] newExponents = new int[exponents.length];
        for (int i = 0; i < exponents.length; i++) {
            newExponents[i] = exponents[i] - oth.exponents[i];
            if (newExponents[i] < 0)
                return null;
        }
        return new MonomialTerm<>(newExponents, totalDegree - oth.totalDegree, coefficient);
    }

    MonomialTerm<E> without(int i) {
        return without(i, coefficient);
    }

    MonomialTerm<E> without(int i, E coefficient) {
        if (exponents.length == 1) {
            assert i == 0;
            return new MonomialTerm<>(new int[0], 0, coefficient);
        }
        return new MonomialTerm<>(ArraysUtil.remove(exponents, i), totalDegree - exponents[i], coefficient);
    }

    MonomialTerm<E> setZero(int i) {
        return setZero(i, coefficient);
    }

    MonomialTerm<E> setZero(int i, E coefficient) {
        if (exponents.length == 1) {
            assert i == 0;
            return new MonomialTerm<>(new int[0], 0, coefficient);
        }
        int[] newExponents = exponents.clone();
        newExponents[i] = 0;
        return new MonomialTerm<>(newExponents, totalDegree - exponents[i], coefficient);
    }

    MonomialTerm<E> set(int i, int exponent) {
        int[] newExponents = exponents.clone();
        newExponents[i] = exponent;
        return new MonomialTerm<>(newExponents, totalDegree - exponents[i] + exponent, coefficient);
    }

    public String toString(String[] vars) {
        List<String> result = new ArrayList<>();
        for (int i = 0; i < exponents.length; i++)
            result.add(toString0(vars[i], exponents[i]));
        return result.stream().filter(s -> !s.isEmpty()).collect(Collectors.joining("*"));
    }

    @Override
    public String toString() {
        return toString(MultivariatePolynomial.defaultVars(exponents.length));
    }

    public String toStringArray() {
        return Arrays.toString(exponents);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MonomialTerm that = (MonomialTerm) o;

        if (totalDegree != that.totalDegree) return false;
        if (!Arrays.equals(exponents, that.exponents)) return false;
        return coefficient != null ? coefficient.equals(that.coefficient) : that.coefficient == null;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(exponents);
        result = 31 * result + totalDegree;
        result = 31 * result + (coefficient != null ? coefficient.hashCode() : 0);
        return result;
    }

    /** cached zero degree vectors */
    private static int[][] zeroDegreeVectors = new int[32][];

    static {
        for (int i = 0; i < zeroDegreeVectors.length; i++)
            zeroDegreeVectors[i] = new int[i];
    }

    public static <E> MonomialTerm<E> withZeroExponents(int nVariables, E coefficient) {
        int[] exponents = nVariables < zeroDegreeVectors.length ? zeroDegreeVectors[nVariables] : new int[nVariables];
        return new MonomialTerm<>(exponents, 0, coefficient);
    }
    
    /**
     * Lexicographic monomial order
     */
    public static final Comparator<MonomialTerm> LEX = (MonomialTerm a, MonomialTerm b) -> {
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
    public static final Comparator<MonomialTerm> ALEX = (MonomialTerm a, MonomialTerm b) -> LEX.compare(b, a);

    /**
     * Graded lexicographic monomial order
     */
    public static final Comparator<MonomialTerm> GRLEX = (MonomialTerm a, MonomialTerm b) -> {
        int c = Integer.compare(a.totalDegree, b.totalDegree);
        return c != 0 ? c : LEX.compare(a, b);
    };

    /**
     * Graded reverse lexicographic monomial order
     */
    public static final Comparator<MonomialTerm> GREVLEX = (MonomialTerm a, MonomialTerm b) -> {
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
