package cc.r2.core.poly2.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.util.ArraysUtil;

import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomial {
    final int nVariables;
    final Comparator<DegreeVector> ordering;
    final TreeMap<DegreeVector, BigInteger> data;

    private MultivariatePolynomial(int nVariables, Comparator<DegreeVector> ordering, TreeMap<DegreeVector, BigInteger> data) {
        this.nVariables = nVariables;
        this.ordering = ordering;
        this.data = data;
    }

    /**
     * Creates multivariate polynomial from a list of coefficients and corresponding degree vectors
     *
     * @param factors  coefficients
     * @param vectors  degree vectors
     * @param ordering term ordering
     * @return multivariate polynomial
     */
    public static MultivariatePolynomial create(BigInteger[] factors, DegreeVector[] vectors, Comparator<DegreeVector> ordering) {
        if (factors.length != vectors.length)
            throw new IllegalArgumentException();
        TreeMap<DegreeVector, BigInteger> map = new TreeMap<>(ordering);
        for (int i = 0; i < factors.length; i++) {
            BigInteger f = factors[i];
            map.compute(vectors[i],
                    (thisVector, thisValue) -> thisValue == null ? f : thisValue.add(f));
        }
        return new MultivariatePolynomial(vectors[0].exponents.length, ordering, map);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    string polynomials
     * @param ordering  term ordering
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static MultivariatePolynomial parse(String string, Comparator<DegreeVector> ordering, String... variables) {
        return Parser.parse(string, ordering, variables);
    }

    private void ensureCompatible(MultivariatePolynomial oth) {
        if (nVariables != oth.nVariables)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /**
     * Returns the underlying term ordering.
     *
     * @return underlying term ordering
     */
    public Comparator<DegreeVector> getOrdering() {
        return ordering;
    }

    /**
     * Returns the leading coefficient of this polynomial.
     *
     * @return leading coefficient of this polynomial
     */
    public BigInteger lc() {
        return data.firstEntry().getValue();
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public BigInteger cc() {
        return data.getOrDefault(zeroDegreeVector(nVariables), BigInteger.ZERO);
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public BigInteger content() {
        BigInteger gcd = null;
        for (BigInteger cf : data.values()) {
            if (gcd == null)
                gcd = cf;
            else
                gcd = gcd.gcd(cf);
            if (gcd.isOne())
                break;
        }
        return gcd;
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public MultivariatePolynomial divideOrNull(BigInteger factor) {
        if (factor.isOne())
            return this;
        for (Map.Entry<DegreeVector, BigInteger> entry : data.entrySet()) {
            BigInteger[] qd = entry.getValue().divideAndRemainder(factor);
            if (!qd[1].isZero())
                return null;
            entry.setValue(qd[0]);
        }
        return this;
    }

    /**
     * Negates this and returns
     *
     * @return this negated
     */
    public MultivariatePolynomial negate() {
        for (Map.Entry<DegreeVector, BigInteger> entry : data.entrySet())
            entry.setValue(entry.getValue().negate());
        return this;
    }

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public MultivariatePolynomial multiply(BigInteger factor) {
        if (factor.isOne())
            return this;
        for (Map.Entry<DegreeVector, BigInteger> entry : data.entrySet())
            entry.setValue(entry.getValue().multiply(factor));
        return this;
    }


    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial add(MultivariatePolynomial oth) {
        ensureCompatible(oth);
        oth.data.entrySet().forEach(othElement -> add(data, othElement.getKey(), othElement.getValue()));
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial subtract(MultivariatePolynomial oth) {
        ensureCompatible(oth);
        oth.data.entrySet().forEach(othElement -> add(data, othElement.getKey(), othElement.getValue().negate()));
        return this;
    }

    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial add(BigInteger oth) {
        add(data, zeroDegreeVector(nVariables), oth);
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial subtract(BigInteger oth) {
        add(data, zeroDegreeVector(nVariables), oth.negate());
        return this;
    }

    /**
     * Adds 1 to this
     *
     * @return {@code this + 1}
     */
    public MultivariatePolynomial increment() {
        return add(BigInteger.ONE);
    }

    /**
     * Subtracts 1 from this
     *
     * @return {@code this - 1}
     */
    public MultivariatePolynomial decrement() {
        return subtract(BigInteger.ONE);
    }

    private static void add(TreeMap<DegreeVector, BigInteger> map, DegreeVector deg, BigInteger val) {
        map.compute(deg, (thisVector, thisValue) -> {
            if (thisValue == null)
                return val;
            BigInteger r = thisValue.add(val);
            return r.isZero() ? null : r;
        });
    }

    /**
     * Multiplies this by {@code oth} and returns this
     *
     * @param oth other polynomial
     * @return {@code this * oth}
     */
    public MultivariatePolynomial multiply(MultivariatePolynomial oth) {
        ensureCompatible(oth);
        TreeMap<DegreeVector, BigInteger> newMap = new TreeMap<>(ordering);
        for (Map.Entry<DegreeVector, BigInteger> othElement : oth.data.entrySet()) {
            for (Map.Entry<DegreeVector, BigInteger> thisElement : data.entrySet()) {
                add(newMap, thisElement.getKey().multiply(othElement.getKey()), thisElement.getValue().multiply(othElement.getValue()));
            }
        }
        data.clear();
        data.putAll(newMap);
        return this;
    }

    @SuppressWarnings("unchecked")
    public MultivariatePolynomial clone() {
        return new MultivariatePolynomial(nVariables, ordering, (TreeMap) data.clone());
    }


    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        MultivariatePolynomial that = (MultivariatePolynomial) o;

        if (nVariables != that.nVariables)
            return false;
        return data.equals(that.data);
    }

    @Override
    public int hashCode() {
        int result = nVariables;
        result = 31 * result + data.hashCode();
        return result;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Map.Entry<DegreeVector, BigInteger> term : data.entrySet()) {
            DegreeVector monomial = term.getKey();
            BigInteger coeff = term.getValue();
            if (coeff.isZero())
                continue;
            String monomialString = monomial.toString();
            if (first) {
                if (!coeff.isOne() || monomialString.isEmpty()) {
                    sb.append(coeff);
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
                first = false;
            } else {
                if (coeff.signum() > 0)
                    sb.append("+");
                if (coeff.isMinusOne())
                    sb.append("-");
                else if (!coeff.isOne() || monomialString.isEmpty()) {
                    sb.append(coeff);
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
            }
        }
        return sb.length() == 0 ? "0" : sb.toString();
    }


    /* ******************************************* Static things ******************************************* */

    /** cached zero degree vectors */
    private static DegreeVector[] zeroDegreeVectors = new DegreeVector[32];

    static {
        for (int i = 0; i < zeroDegreeVectors.length; i++)
            zeroDegreeVectors[i] = new DegreeVector(new int[i], 0);
    }

    private static DegreeVector zeroDegreeVector(int nVariables) {
        if (nVariables < zeroDegreeVectors.length)
            return zeroDegreeVectors[nVariables];
        return new DegreeVector(new int[nVariables], 0);
    }


    /**
     * Degree vector.
     */
    public static final class DegreeVector {
        private final int[] exponents;
        private final int totalDegree;

        private DegreeVector(int[] exponents, int totalDegree) {
            this.exponents = exponents;
            this.totalDegree = totalDegree;
        }

        public DegreeVector(int... exponents) {
            this.exponents = exponents;
            this.totalDegree = ArraysUtil.sum(exponents);
        }

        private static String toString0(char var, int exp) {
            return exp == 0 ? "" : Character.toString(var) + (exp == 1 ? "" : "^" + exp);
        }

        private DegreeVector multiply(DegreeVector dv) {
            int[] newExponents = new int[exponents.length];
            for (int i = 0; i < exponents.length; i++)
                newExponents[i] = exponents[i] + dv.exponents[i];
            return new DegreeVector(newExponents, totalDegree + dv.totalDegree);
        }

        @Override
        public String toString() {
            List<String> result = new ArrayList<>();
            char x = 'a';
            for (int exponent : exponents)
                result.add(toString0(x++, exponent));
            return result.stream().filter(s -> !s.isEmpty()).collect(Collectors.joining("*"));
        }

        public String toStringArray() {
            return Arrays.toString(exponents);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            DegreeVector that = (DegreeVector) o;

            if (totalDegree != that.totalDegree) return false;
            return Arrays.equals(exponents, that.exponents);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(exponents);
            result = 31 * result + totalDegree;
            return result;
        }
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
     * Reverse lexicographic monomial order
     */
    public static final Comparator<DegreeVector> REVLEX = (DegreeVector a, DegreeVector b) -> LEX.compare(b, a);

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
