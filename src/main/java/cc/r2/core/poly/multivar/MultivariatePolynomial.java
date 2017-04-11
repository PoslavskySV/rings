package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.util.ArraysUtil;

import java.util.*;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import static cc.r2.core.number.BigInteger.ZERO;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomial implements IGeneralPolynomial<MultivariatePolynomial> {
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
        if (factors.length == 0)
            throw new IllegalArgumentException("empty");
        TreeMap<DegreeVector, BigInteger> map = new TreeMap<>(ordering);
        for (int i = 0; i < factors.length; i++) {
            BigInteger f = factors[i];
            if (!f.isZero())
                add(map, vectors[i], f);
        }
        return new MultivariatePolynomial(vectors[0].exponents.length, ordering, map);
    }

    /**
     * Creates zero
     *
     * @param nVariable number of variables
     * @param ordering  the orderging
     * @return zero
     */
    public static MultivariatePolynomial zero(int nVariable, Comparator<DegreeVector> ordering) {
        return new MultivariatePolynomial(nVariable, ordering, new TreeMap<>(ordering));
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

    /** check whether number of variables is the same */
    private void ensureCompatible(MultivariatePolynomial oth) {
        if (nVariables != oth.nVariables)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /** check whether number of variables is the same */
    private void ensureCompatible(Entry<DegreeVector, BigInteger> oth) {
        if (nVariables != oth.getKey().exponents.length)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    private void release() {}

    /**
     * Returns the underlying term ordering.
     *
     * @return underlying term ordering
     */
    public Comparator<DegreeVector> getOrdering() {
        return ordering;
    }

    /**
     * Copies this with the new ordering {@code newOrdering}
     *
     * @param newOrdering the new ordering
     * @return a copy of this with a new ordering
     */
    public MultivariatePolynomial setOrdering(Comparator<DegreeVector> newOrdering) {
        TreeMap<DegreeVector, BigInteger> newData = new TreeMap<>(newOrdering);
        newData.putAll(data);
        return new MultivariatePolynomial(nVariables, newOrdering, newData);
    }

    /**
     * Returns the number of variables in this
     *
     * @return the number of variables
     */
    public int nVariables() {
        return nVariables;
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public MultivariatePolynomial createConstant(BigInteger val) {
        TreeMap<DegreeVector, BigInteger> data = new TreeMap<>(ordering);
        if (!val.isZero())
            data.put(zeroDegreeVector(nVariables), val);
        return new MultivariatePolynomial(nVariables, ordering, data);
    }

    @Override
    public MultivariatePolynomial createZero() {
        return createConstant(ZERO);
    }

    @Override
    public MultivariatePolynomial createOne() {
        return createConstant(BigInteger.ONE);
    }

    @Override
    public MultivariatePolynomial toZero() {
        data.clear();
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial set(MultivariatePolynomial oth) {
        ensureCompatible(oth);
        data.clear();
        data.putAll(oth.data);
        release();
        return null;
    }

    /**
     * Returns the number of terms in this polynomial
     *
     * @return number of terms
     */
    public int size() {return data.size();}

    @Override
    public boolean isZero() {
        return data.size() == 0;
    }

    @Override
    public boolean isOne() {
        return size() == 1 && data.firstEntry().getValue().isOne();
    }

    @Override
    public boolean isUnitCC() {
        return cc().isOne();
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && data.firstEntry().getKey().isZeroVector());
    }

    @Override
    public boolean isMonomial() {
        return size() <= 1;
    }

    /**
     * Returns the total degree of this polynomial, that is the maximal degree among all terms
     *
     * @return the total degree of this polynomial, that is the maximal degree among all terms
     */
    @Override
    public int degree() {
        int max = 0;
        for (DegreeVector db : data.keySet())
            max = Math.max(max, db.totalDegree);
        return max;
    }

    /**
     * Returns the leading term in this polynomial according to ordering
     *
     * @return the leading term in this polynomial according to ordering
     */
    public Entry<DegreeVector, BigInteger> lt() {
        return size() == 0 ? new EntryImpl(zeroDegreeVector(nVariables), ZERO) : data.lastEntry();
    }

    /**
     * Returns the largest degree vector with respect to this ordering
     *
     * @return the largest degree vector
     */
    public DegreeVector multiDegree() {
        return lt().getKey();
    }

    /**
     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the ordering.
     *
     * @return leading coefficient of this polynomial
     */
    public BigInteger lc() {
        return lt().getValue();
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public BigInteger cc() {
        return data.getOrDefault(zeroDegreeVector(nVariables), ZERO);
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

    @Override
    public MultivariatePolynomial primitivePart() {
        MultivariatePolynomial r = divideOrNull(content());
        assert r != null;
        release();
        return r;
    }

    @Override
    public MultivariatePolynomial primitivePartSameSign() {
        BigInteger c = content();
        if (c.signum() < 0)
            c = c.negate();
        MultivariatePolynomial r = divideOrNull(c);
        assert r != null;
        release();
        return r;
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
        for (Entry<DegreeVector, BigInteger> entry : data.entrySet()) {
            BigInteger[] qd = entry.getValue().divideAndRemainder(factor);
            if (!qd[1].isZero())
                return null;
            entry.setValue(qd[0]);
        }
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial negate() {
        for (Entry<DegreeVector, BigInteger> entry : data.entrySet())
            entry.setValue(entry.getValue().negate());
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial add(MultivariatePolynomial oth) {
        if (data == oth.data)
            return multiply(2);
        ensureCompatible(oth);
        oth.data.entrySet().forEach(othElement -> add(data, othElement.getKey(), othElement.getValue()));
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial subtract(MultivariatePolynomial oth) {
        if (data == oth.data)
            return toZero();
        ensureCompatible(oth);
        oth.data.entrySet().forEach(othElement -> add(data, othElement.getKey(), othElement.getValue().negate()));
        release();
        return this;
    }

    /**
     * Adds {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this + oth}
     */
    MultivariatePolynomial add(Entry<DegreeVector, BigInteger> term) {
        ensureCompatible(term);
        add(data, term.getKey(), term.getValue());
        release();
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
        release();
        return this;
    }

    /**
     * Subtracts {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this - oth}
     */
    MultivariatePolynomial subtract(Entry<DegreeVector, BigInteger> term) {
        ensureCompatible(term);
        add(data, term.getKey(), term.getValue().negate());
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public MultivariatePolynomial subtract(BigInteger oth) {
        add(data, zeroDegreeVector(nVariables), oth.negate());
        release();
        return this;
    }

    /**
     * Removes the leading term from this polynomial
     *
     * @return this - this.lt()
     */
    public MultivariatePolynomial subtractLt() {
        data.pollLastEntry();
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial increment() {
        return add(BigInteger.ONE);
    }

    @Override
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
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public MultivariatePolynomial multiply(BigInteger factor) {
        if (factor.isOne())
            return this;
        if (factor.isZero())
            return toZero();
        for (Entry<DegreeVector, BigInteger> entry : data.entrySet())
            entry.setValue(entry.getValue().multiply(factor));
        release();
        return this;
    }

    private MultivariatePolynomial loadFrom(TreeMap<DegreeVector, BigInteger> map) {
        data.clear();
        data.putAll(map);
        release();
        return this;
    }

    /**
     * Multiplies {@code this} by the {@code term}
     *
     * @param term the factor
     * @return {@code} this multiplied by the {@code term}
     */
    public MultivariatePolynomial multiply(Entry<DegreeVector, BigInteger> term) {
        ensureCompatible(term);
        if (term.getValue().isZero())
            return toZero();

        TreeMap<DegreeVector, BigInteger> newMap = new TreeMap<>(ordering);
        for (Entry<DegreeVector, BigInteger> thisElement : data.entrySet())
            newMap.put(thisElement.getKey().multiply(term.getKey()), thisElement.getValue().multiply(term.getValue()));

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial multiply(long factor) {
        return multiply(BigInteger.valueOf(factor));
    }

    @Override
    public MultivariatePolynomial multiply(MultivariatePolynomial oth) {
        ensureCompatible(oth);
        TreeMap<DegreeVector, BigInteger> newMap = new TreeMap<>(ordering);
        for (Entry<DegreeVector, BigInteger> othElement : oth.data.entrySet())
            for (Entry<DegreeVector, BigInteger> thisElement : data.entrySet())
                add(newMap, thisElement.getKey().multiply(othElement.getKey()), thisElement.getValue().multiply(othElement.getValue()));

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial square() {
        return multiply(this);
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

    public String toString(String... vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Entry<DegreeVector, BigInteger> term : data.entrySet()) {
            DegreeVector monomial = term.getKey();
            BigInteger coeff = term.getValue();
            if (coeff.isZero())
                continue;
            String monomialString = monomial.toString(vars);
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

    @Override
    public String toString() {
        return toString(defaultVars(nVariables));
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

    private static String[] defaultVars(int nVars) {
        char v = 'a';
        String[] vars = new String[nVars];
        for (int i = 0; i < nVars; i++)
            vars[i] = Character.toString(v++);
        return vars;
    }


    /**
     * Degree vector.
     */
    public static final class DegreeVector {
        final int[] exponents;
        final int totalDegree;

        private DegreeVector(int[] exponents, int totalDegree) {
            this.exponents = exponents;
            this.totalDegree = totalDegree;
        }

        public DegreeVector(int... exponents) {
            this.exponents = exponents;
            this.totalDegree = ArraysUtil.sum(exponents);
        }

        boolean isZeroVector() {
            return totalDegree == 0;
        }

        private static String toString0(String var, int exp) {
            return exp == 0 ? "" : var + (exp == 1 ? "" : "^" + exp);
        }

        private DegreeVector multiply(DegreeVector dv) {
            int[] newExponents = new int[exponents.length];
            for (int i = 0; i < exponents.length; i++)
                newExponents[i] = exponents[i] + dv.exponents[i];
            return new DegreeVector(newExponents, totalDegree + dv.totalDegree);
        }

        public String toString(String[] vars) {
            List<String> result = new ArrayList<>();
            for (int i = 0; i < exponents.length; i++)
                result.add(toString0(vars[i], exponents[i]));
            return result.stream().filter(s -> !s.isEmpty()).collect(Collectors.joining("*"));
        }

        @Override
        public String toString() {
            return toString(defaultVars(exponents.length));
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

    static final class EntryImpl implements Entry<DegreeVector, BigInteger> {
        final DegreeVector degreeVector;
        final BigInteger coefficient;

        public EntryImpl(DegreeVector degreeVector, BigInteger coefficient) {
            this.degreeVector = degreeVector;
            this.coefficient = coefficient;
        }

        @Override
        public DegreeVector getKey() {
            return degreeVector;
        }

        @Override
        public BigInteger getValue() {
            return coefficient;
        }

        @Override
        public BigInteger setValue(BigInteger value) {
            throw new IllegalStateException();
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            EntryImpl entry = (EntryImpl) o;

            if (!degreeVector.equals(entry.degreeVector)) return false;
            return coefficient.equals(entry.coefficient);
        }

        @Override
        public int hashCode() {
            int result = degreeVector.hashCode();
            result = 31 * result + coefficient.hashCode();
            return result;
        }
    }
}
