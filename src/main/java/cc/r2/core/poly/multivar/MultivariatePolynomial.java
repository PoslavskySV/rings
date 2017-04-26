package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.poly.generics.Domain;
import cc.r2.core.poly.generics.ModularDomain;
import cc.r2.core.poly.generics.UnivariatePolynomialDomain;
import cc.r2.core.poly.univar.bMutablePolynomialZ;
import cc.r2.core.poly.univar.bMutablePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.set.hash.TIntHashSet;

import java.util.*;
import java.util.Map.Entry;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static cc.r2.core.poly.generics.IntegersDomain.IntegersDomain;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomial<E> implements IGeneralPolynomial<MultivariatePolynomial<E>> {
    public final Domain<E> domain;
    public final Comparator<DegreeVector> ordering;
    public final int nVariables;
    final TreeMap<DegreeVector, E> data;

    private MultivariatePolynomial(Domain<E> domain, Comparator<DegreeVector> ordering, int nVariables, TreeMap<DegreeVector, E> data) {
        this.nVariables = nVariables;
        this.ordering = ordering;
        this.data = data;
        this.domain = domain;
    }

    /**
     * Creates multivariate polynomial from a list of coefficients and corresponding degree vectors
     *
     * @param domain   domain
     * @param ordering term ordering
     * @param vectors  degree vectors
     * @param factors  coefficients
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> create(Domain<E> domain, Comparator<DegreeVector> ordering, DegreeVector[] vectors, E[] factors) {
        if (factors.length != vectors.length)
            throw new IllegalArgumentException();
        if (factors.length == 0)
            throw new IllegalArgumentException("empty");
        TreeMap<DegreeVector, E> map = new TreeMap<>(ordering);
        for (int i = 0; i < factors.length; i++) {
            E f = domain.valueOf(factors[i]);
            add(map, vectors[i], f, domain);
        }
        return new MultivariatePolynomial<>(domain, ordering, vectors[0].exponents.length, map);
    }

    /**
     * Creates zero
     *
     * @param domain     domain
     * @param ordering   the orderging
     * @param nVariables number of variables
     * @return zero
     */
    public static <E> MultivariatePolynomial<E> zero(Domain<E> domain, Comparator<DegreeVector> ordering, int nVariables) {
        return new MultivariatePolynomial<>(domain, ordering, nVariables, new TreeMap<>(ordering));
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    string polynomials
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static MultivariatePolynomial<BigInteger> parse(String string, String... variables) {
        return Parser.parse(string, LEX, variables);
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
    public static MultivariatePolynomial<BigInteger> parse(String string, Comparator<DegreeVector> ordering, String... variables) {
        return Parser.parse(string, ordering, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    string polynomials
     * @param domain    domain
     * @param ordering  term ordering
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Domain<E> domain, Comparator<DegreeVector> ordering, String... variables) {
        return Parser.parse(string, domain, ordering, variables);
    }

    /**
     * Converts multivariate polynomial to univariate if it is effectively univariate
     *
     * @param poly the multivariate polynomial
     * @return the univariate polynomial
     * @throws IllegalArgumentException if {@code poly} is not actually a univariate polynomial
     */
    public static bMutablePolynomialZ asUnivariateZ(MultivariatePolynomial<BigInteger> poly) {
        return bMutablePolynomialZ.create(asUnivariate(poly));
    }

    /**
     * Converts multivariate polynomial to univariate if it is effectively univariate
     *
     * @param poly the multivariate polynomial
     * @return the univariate polynomial
     * @throws IllegalArgumentException if {@code poly} is not actually a univariate polynomial
     */
    public static bMutablePolynomialZp asUnivariateZp(MultivariatePolynomial<BigInteger> poly) {
        if (!(poly.domain instanceof ModularDomain))
            throw new IllegalArgumentException("multivariate poly is not Zp[x]");
        ModularDomain domain = (ModularDomain) poly.domain;
        BigInteger[] data = asUnivariate(poly);
        return bMutablePolynomialZp.createUnsafe(domain.modulus, data);
    }

    private static BigInteger[] asUnivariate(MultivariatePolynomial<BigInteger> poly) {
        int[] degrees = poly.degrees();
        int theVar = -1;
        for (int i = 0; i < degrees.length; i++) {
            if (degrees[i] != 0) {
                if (theVar != -1)
                    throw new IllegalArgumentException("not a univariate polynomial: " + poly);
                theVar = i;
            }
        }

        BigInteger[] data = new BigInteger[degrees[theVar] + 1];
        Arrays.fill(data, BigInteger.ZERO);
        for (Entry<DegreeVector, BigInteger> e : poly.data.entrySet()) {
            assert e.getValue() != null;
            data[e.getKey().exponents[theVar]] = e.getValue();
        }
        return data;
    }

    /**
     * Converts univariate polynomial to multivariate over integers domain.
     *
     * @param poly       univariate polynomial
     * @param nVariables number of variables in the result
     * @param variable   variable that will be used as a primary variable
     * @param ordering   ordering
     * @return multivariate polynomial
     */
    public static MultivariatePolynomial<BigInteger> asMultivariate(bMutablePolynomialZ poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        return asMultivariate(poly.getDataReferenceUnsafe(), poly.degree(), nVariables, variable, IntegersDomain, ordering);
    }

    /**
     * Converts univariate polynomial to multivariate over Zp domain.
     *
     * @param poly       univariate polynomial
     * @param nVariables number of variables in the result
     * @param variable   variable that will be used as a primary variable
     * @param ordering   ordering
     * @return multivariate polynomial
     */
    public static MultivariatePolynomial<BigInteger> asMultivariate(bMutablePolynomialZp poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        return asMultivariate(poly.getDataReferenceUnsafe(), poly.degree(), nVariables, variable, new ModularDomain(poly.modulusAsBigInt()), ordering);
    }

    private static MultivariatePolynomial<BigInteger> asMultivariate(BigInteger[] data, int degree, int nVariables, int variable, Domain<BigInteger> domain, Comparator<DegreeVector> ordering) {
        TreeMap<DegreeVector, BigInteger> map = new TreeMap<>(ordering);
        for (int i = 0; i <= degree; i++) {
            if (domain.isZero(data[i]))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;

            map.put(new DegreeVector(degreeVector, i), data[i]);
        }
        return new MultivariatePolynomial<>(domain, ordering, nVariables, map);
    }

    public static MultivariatePolynomial<bMutablePolynomialZ> convertZ(MultivariatePolynomial<BigInteger> poly, int variable) {
        TreeMap<DegreeVector, bMutablePolynomialZ> map = new TreeMap<>(poly.ordering);
        UnivariatePolynomialDomain<bMutablePolynomialZ> domain = new UnivariatePolynomialDomain<>(bMutablePolynomialZ.zero());
        for (Entry<DegreeVector, BigInteger> e : poly.data.entrySet()) {
            DegreeVector oldDV = e.getKey();
            DegreeVector newDV = oldDV.without(variable);
            add(map, newDV, bMutablePolynomialZ.monomial(e.getValue(), oldDV.exponents[variable]), domain);
        }
        return new MultivariatePolynomial<>(domain, poly.ordering, poly.nVariables - 1, map);
    }

    public static MultivariatePolynomial<bMutablePolynomialZp> convertZp(MultivariatePolynomial<BigInteger> poly, int variable) {
        if (!(poly.domain instanceof ModularDomain))
            throw new IllegalArgumentException("not a Zp[x] poly");
        BigInteger modulus = ((ModularDomain) poly.domain).modulus;
        TreeMap<DegreeVector, bMutablePolynomialZp> map = new TreeMap<>(poly.ordering);
        UnivariatePolynomialDomain<bMutablePolynomialZp> domain = new UnivariatePolynomialDomain<>(bMutablePolynomialZp.zero(modulus));
        for (Entry<DegreeVector, BigInteger> e : poly.data.entrySet()) {
            DegreeVector oldDV = e.getKey();
            DegreeVector newDV = oldDV.without(variable);
            add(map, newDV, bMutablePolynomialZp.monomial(modulus, e.getValue(), oldDV.exponents[variable]), domain);
        }
        return new MultivariatePolynomial<>(domain, poly.ordering, poly.nVariables - 1, map);
    }

    public static MultivariatePolynomial<BigInteger> fromZp(MultivariatePolynomial<bMutablePolynomialZp> poly, Domain<BigInteger> domain, int variable) {
        int nVariables = poly.nVariables + 1;
        MultivariatePolynomial<BigInteger> result = zero(domain, poly.ordering, nVariables);
        for (Entry<DegreeVector, bMutablePolynomialZp> entry : poly.data.entrySet()) {
            bMutablePolynomialZp upoly = entry.getValue();
            int[] dv = ArraysUtil.insert(entry.getKey().exponents, variable, 0);
            for (int i = 0; i <= upoly.degree(); ++i) {
                if (upoly.get(i).isZero())
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new EntryImpl<>(new DegreeVector(cdv), upoly.get(i)));
            }
        }
        return result;
    }

    /**
     * Renames variable {@code i} to {@code j} and {@code j} to {@code i}
     *
     * @param poly the polynomial
     * @param i    the first variable
     * @param j    the second variable
     * @return polynomial with variable {@code i} renamed to {@code j} and {@code j} renamed to {@code i}
     */
    public static <E> MultivariatePolynomial<E> swapVariables(MultivariatePolynomial<E> poly, int i, int j) {
        int[] newVariables = ArraysUtil.sequence(poly.nVariables);
        newVariables[i] = j;
        newVariables[j] = i;
        return renameVariables(poly, newVariables, poly.ordering);
    }

    /**
     * Rename variables from [0,1,...N] to [newVariables[0], newVariables[1], ..., newVariables[N]]
     *
     * @param poly         the polynomial
     * @param newVariables the new variables
     * @return renamed polynomial
     */
    public static <E> MultivariatePolynomial<E> renameVariables(MultivariatePolynomial<E> poly, int[] newVariables) {
        return renameVariables(poly, newVariables, poly.ordering);
    }

    /**
     * Rename variables from [0,1,...N] to [newVariables[0], newVariables[1], ..., newVariables[N]]
     *
     * @param poly         the polynomial
     * @param newVariables the new variables
     * @param newOrdering  the new ordering
     * @return renamed polynomial
     */
    public static <E> MultivariatePolynomial<E> renameVariables(MultivariatePolynomial<E> poly, int[] newVariables, Comparator<DegreeVector> newOrdering) {
        // NOTE: always return a copy of poly, even if order of variables is unchanged
        TreeMap<DegreeVector, E> data = new TreeMap<>(newOrdering);
        for (Entry<DegreeVector, E> e : poly.data.entrySet()) {
            DegreeVector dv = e.getKey();
            data.put(new DegreeVector(map(dv.exponents, newVariables), dv.totalDegree), e.getValue());
        }
        return new MultivariatePolynomial<>(poly.domain, newOrdering, poly.nVariables, data);
    }

    private static int[] map(int[] degrees, int[] mapping) {
        int[] newDegrees = new int[degrees.length];
        for (int i = 0; i < degrees.length; i++)
            newDegrees[i] = degrees[mapping[i]];
        return newDegrees;
    }

    /** check whether number of variables is the same */
    private void ensureCompatible(MultivariatePolynomial<E> oth) {
        if (nVariables != oth.nVariables)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /** check whether number of variables is the same */
    private void ensureCompatible(Entry<DegreeVector, E> oth) {
        ensureCompatible(oth.getKey());
    }

    /** check whether number of variables is the same */
    private void ensureCompatible(DegreeVector oth) {
        if (nVariables != oth.exponents.length)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /** release caches */
    private void release() {}

    /**
     * Copies this with the new ordering {@code newOrdering}
     *
     * @param newOrdering the new ordering
     * @return a copy of this with a new ordering
     */
    public MultivariatePolynomial<E> setOrdering(Comparator<DegreeVector> newOrdering) {
        TreeMap<DegreeVector, E> newData = new TreeMap<>(newOrdering);
        newData.putAll(data);
        return new MultivariatePolynomial<>(domain, newOrdering, nVariables, newData);
    }

    /**
     * Switches to another domain soecified by {@code newDomain}
     *
     * @param newDomain the new domain
     * @return a copy of this reduced to the domain specified by {@code newDomain}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> setDomain(Domain<E> newDomain) {
        TreeMap<DegreeVector, E> newData = (TreeMap) data.clone();
        for (Entry<DegreeVector, E> e : newData.entrySet())
            e.setValue(newDomain.valueOf(e.getValue()));
        return new MultivariatePolynomial<>(newDomain, ordering, nVariables, newData);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public MultivariatePolynomial<E> createConstant(E val) {
        TreeMap<DegreeVector, E> data = new TreeMap<>(ordering);
        if (!domain.isZero(val))
            data.put(zeroDegreeVector(nVariables), val);
        return new MultivariatePolynomial<>(domain, ordering, nVariables, data);
    }

    /**
     * Creates multivariate polynomial with one single specified term
     *
     * @param dv  term degree vector
     * @param val term coefficient
     * @return multivariate polynomial with one single specified term
     */
    public MultivariatePolynomial<E> create(DegreeVector dv, E val) {
        TreeMap<DegreeVector, E> data = new TreeMap<>(ordering);
        val = domain.valueOf(val);
        if (!domain.isZero(val))
            data.put(dv, val);
        return new MultivariatePolynomial<>(domain, ordering, nVariables, data);
    }

    /**
     * Creates multivariate polynomial from a list of coefficients and corresponding degree vectors
     *
     * @param vectors degree vectors
     * @param factors coefficients
     * @return multivariate polynomial
     */
    public MultivariatePolynomial<E> create(DegreeVector[] vectors, E[] factors) {
        if (factors.length != vectors.length)
            throw new IllegalArgumentException();
        if (factors.length == 0)
            throw new IllegalArgumentException("empty");
        TreeMap<DegreeVector, E> map = new TreeMap<>(ordering);
        for (int i = 0; i < factors.length; i++) {
            E f = factors[i];
            add(map, vectors[i], f, domain);
        }
        return new MultivariatePolynomial<>(domain, ordering, vectors[0].exponents.length, map);
    }

    @Override
    public MultivariatePolynomial<E> createZero() {
        return createConstant(domain.getZero());
    }

    @Override
    public MultivariatePolynomial<E> createOne() {
        return createConstant(domain.getOne());
    }

    /**
     * Creates linear polynomial of the form {@code cc + lc * variable}
     *
     * @param variable the variable
     * @param cc       the constant coefficient
     * @param lc       the leadingcoefficient
     * @return linear polynomial {@code cc + lc * variable}
     */
    public MultivariatePolynomial<E> createLinear(int variable, E cc, E lc) {
        int[] ccDegreeVector = new int[nVariables], lcDegreeVector = new int[nVariables];
        lcDegreeVector[variable] = 1;
        TreeMap<DegreeVector, E> data = new TreeMap<>(ordering);
        if (!domain.isZero(cc))
            data.put(new DegreeVector(ccDegreeVector, 0), cc);
        if (!domain.isZero(lc))
            data.put(new DegreeVector(lcDegreeVector, 1), lc);
        return new MultivariatePolynomial<E>(domain, ordering, nVariables, data);
    }

    @Override
    public MultivariatePolynomial<E> toZero() {
        data.clear();
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> set(MultivariatePolynomial<E> oth) {
        ensureCompatible(oth);
        data.clear();
        data.putAll(oth.data);
        release();
        return null;
    }

    /**
     * Returns a copy of this with {@code nVariables + 1}
     *
     * @return a copy of this with with one additional (last) variable added
     */
    public MultivariatePolynomial<E> joinNewVariable() {
        TreeMap<DegreeVector, E> map = new TreeMap<>(ordering);
        for (Entry<DegreeVector, E> e : data.entrySet())
            map.put(e.getKey().joinNewVariable(), e.getValue());
        return new MultivariatePolynomial<>(domain, ordering, nVariables + 1, map);
    }

    /**
     * Returns the number of really used variables
     *
     * @return the number of presenting variables
     */
    public int usedVariables() {
        int[] degrees = degrees();
        int r = 0;
        for (int i = 0; i < degrees.length; i++)
            if (degrees[i] != 0)
                ++r;
        return r;
    }

    /**
     * Returns a coefficient before {@code variable^exponent} as a multivariate polynomial
     *
     * @param variable the variable
     * @param exponent the exponent
     * @return coefficient before {@code variable^exponent} as a multivariate polynomial
     */
    public MultivariatePolynomial<E> coefficientOf(int variable, int exponent) {
        MultivariatePolynomial<E> result = createZero();
        for (Entry<DegreeVector, E> e : data.entrySet()) {
            DegreeVector dv = e.getKey();
            if (dv.exponents[variable] != exponent)
                continue;
            result.add(dv.setZero(variable), e.getValue());
        }
        return result;
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
        return size() == 1 && domain.isOne(data.firstEntry().getValue());
    }

    @Override
    public boolean isUnitCC() {
        return domain.isOne(cc());
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
     * Returns whether this poly is effectively univariate
     *
     * @return whether this poly is effectively univariate
     */
    public boolean isEffectiveUnivariate() {
        int[] degrees = degrees();
        boolean b = false;
        for (int i = 0; i < nVariables; i++) {
            if (degrees[i] != 0) {
                if (b)
                    return false;
                else
                    b = true;
            }
        }
        return true;
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
     * Returns the degree of this polynomial with respect to i-th variable
     *
     * @return the degree of this polynomial with respect to i-th variable
     */
    public int degree(int i) {
        int max = 0;
        for (DegreeVector db : data.keySet())
            max = Math.max(max, db.exponents[i]);
        return max;
    }

    /**
     * Returns an array of degrees of all variables, so that is i-th element of the result is the polynomial degree
     * with respect to i-th variable
     *
     * @return array of degrees
     */
    public int[] degrees() {
        int[] degrees = new int[nVariables];
        for (DegreeVector db : data.keySet())
            for (int i = 0; i < nVariables; i++)
                if (db.exponents[i] > degrees[i])
                    degrees[i] = db.exponents[i];
        return degrees;
    }

    /**
     * Returns the degrees in which {@code variable} occurs in this polynomial
     *
     * @return the degrees in which {@code variable} occurs in this polynomial
     */
    public int[] degrees(int variable) {
        TIntHashSet degrees = new TIntHashSet();
        for (DegreeVector db : data.keySet())
            degrees.add(db.exponents[variable]);
        return degrees.toArray();
    }

    /**
     * Returns the product of {@link #degrees()}
     *
     * @return product of {@link #degrees()}
     */
    public int degreeProduct() {
        int r = 1;
        for (int d : degrees())
            r *= d;
        return r;
    }

    /**
     * Returns the leading term in this polynomial according to ordering
     *
     * @return the leading term in this polynomial according to ordering
     */
    public Entry<DegreeVector, E> lt() {
        return size() == 0 ? new EntryImpl<>(zeroDegreeVector(nVariables), domain.getZero()) : data.lastEntry();
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
    public E lc() {
        return lt().getValue();
    }

    /**
     * Sets the leading coefficient to the specified value
     *
     * @param val new value for the lc
     * @return the leading coefficient to the specified value
     */
    public MultivariatePolynomial<E> setLC(E val) {
        if (isZero())
            return add(val);
        val = domain.valueOf(val);
        Entry<DegreeVector, E> lt = lt();
        data.put(lt.getKey(), val);
        return this;
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public E cc() {
        return data.getOrDefault(zeroDegreeVector(nVariables), domain.getZero());
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public E content() {
        E gcd = null;
        for (E cf : data.values()) {
            if (gcd == null)
                gcd = cf;
            else
                gcd = domain.gcd(gcd, cf);
            if (domain.isOne(gcd))
                break;
        }
        return gcd == null ? domain.getOne() : gcd;
    }

    /**
     * Returns the monomial content of this polynomial
     *
     * @return the monomial content of this polynomial
     */
    public DegreeVector monomialContent() {
        return commonContent(null);
    }

    /**
     * Returns common content of {@code this} and {@code monomial}
     *
     * @param monomial the monomial
     * @return common monomial factor of {@code this} and {@code monomial}
     */
    DegreeVector commonContent(DegreeVector monomial) {
        int[] exponents = monomial == null ? null : monomial.exponents.clone();
        for (DegreeVector degreeVector : data.keySet())
            if (exponents == null)
                exponents = degreeVector.exponents.clone();
            else
                setMin(degreeVector, exponents);
        if (exponents == null)
            return zeroDegreeVector(nVariables);
        return new DegreeVector(exponents);
    }

    static void setMin(DegreeVector degreeVector, int[] exponents) {
        int[] dv = degreeVector.exponents;
        for (int i = 0; i < exponents.length; ++i)
            if (dv[i] < exponents[i])
                exponents[i] = dv[i];
    }

    @Override
    public MultivariatePolynomial<E> primitivePart() {
        MultivariatePolynomial<E> r = divideOrNull(content());
        assert r != null;
        release();
        return r;
    }

    @Override
    public MultivariatePolynomial<E> primitivePartSameSign() {
        E c = content();
        if (domain.signum(c) < 0)
            c = domain.negate(c);
        MultivariatePolynomial<E> r = divideOrNull(c);
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
    public MultivariatePolynomial<E> divideOrNull(E factor) {
        if (domain.isOne(factor))
            return this;
        for (Entry<DegreeVector, E> entry : data.entrySet()) {
            E[] qd = domain.divideAndRemainder(entry.getValue(), factor);
            if (!domain.isZero(qd[1]))
                return null;
            entry.setValue(qd[0]);
        }
        release();
        return this;
    }

    /**
     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code monomial}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param monomial monomial degrees
     * @param factor   monomial factor
     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
     */
    public MultivariatePolynomial<E> divideOrNull(DegreeVector monomial, E factor) {
        if (monomial.isZeroVector())
            return divideOrNull(factor);
        TreeMap<DegreeVector, E> map = new TreeMap<>(ordering);
        for (Entry<DegreeVector, E> entry : data.entrySet()) {
            DegreeVector dv = entry.getKey().divide(monomial);
            if (dv == null)
                return null;
            E[] qd = domain.divideAndRemainder(entry.getValue(), factor);
            if (!domain.isZero(qd[1]))
                return null;
            map.put(dv, qd[0]);
        }
        loadFrom(map);
        release();
        return this;
    }

    /**
     * Makes this polynomial monic if possible, if not -- destroys this and returns null
     *
     * @return monic this or null if the domain does not support exact division by lc
     */
    public MultivariatePolynomial<E> monic() {
        if (domain.isField())
            return multiply(domain.reciprocal(lc()));
        else
            return divideOrNull(lc());
    }

    /**
     * Substitutes {@code value} for {@code variable}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, Object)})
     * @see #eliminate(int, Object)
     */
    public MultivariatePolynomial<E> evaluate(int variable, E value) {
        value = domain.valueOf(value);
        TreeMap<DegreeVector, E> newData = new TreeMap<>(ordering);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, domain);
        for (Entry<DegreeVector, E> el : data.entrySet()) {
            DegreeVector dv = el.getKey();
            add(newData, dv.setZero(variable), domain.multiply(el.getValue(), powers.pow(dv.exponents[variable])));
        }
        return new MultivariatePolynomial<>(domain, ordering, nVariables, newData);
    }

    /**
     * Substitutes {@code values} for {@code variables}.
     *
     * @param variables the variables
     * @param values    the values
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, Object)})
     * @see #eliminate(int, Object)
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> evaluate(int[] variables, E[] values) {
        return evaluate(new PrecomputedPowersHolder<>(values, domain), ones == null ? ones = ArraysUtil.arrayOf(1, nVariables) : ones, variables);
    }

    private int[] ones = null;

    /** substitutes {@code values} for {@code variables} */
    @SuppressWarnings("unchecked")
    MultivariatePolynomial<E> evaluate(PrecomputedPowersHolder<E> powers, int[] variables, int[] raiseFactors) {
        TreeMap<DegreeVector, E> newData = new TreeMap<>(ordering);
        for (Entry<DegreeVector, E> el : data.entrySet()) {
            DegreeVector dv = el.getKey(), newDv = dv;
            E value = el.getValue();
            for (int i = 0; i < variables.length; ++i) {
                value = domain.multiply(value, powers.pow(i, raiseFactors[i] * dv.exponents[variables[i]]));
                newDv = newDv.setZero(variables[i]);
            }

            add(newData, newDv, value);
        }
        return new MultivariatePolynomial<>(domain, ordering, nVariables, newData);
    }

    /**
     * Evaluates this polynomial at specified points
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[] evaluate(int variable, E... values) {
        return Arrays.stream(values).map(p -> evaluate(variable, p)).toArray(MultivariatePolynomial[]::new);
    }

    /**
     * Substitutes {@code value} for {@code variable}. NOTE: the resulting polynomial will
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     * @see #eliminate(int, long)
     */
    public MultivariatePolynomial<E> evaluate(int variable, long value) {
        return evaluate(variable, domain.valueOf(value));
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
     * @see #evaluate(int, Object)
     */
    public MultivariatePolynomial<E> eliminate(int variable, E value) {
        value = domain.valueOf(value);
        TreeMap<DegreeVector, E> newData = new TreeMap<>(ordering);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, domain);
        for (Entry<DegreeVector, E> el : data.entrySet()) {
            DegreeVector dv = el.getKey();
            add(newData, dv.without(variable), domain.multiply(el.getValue(), powers.pow(dv.exponents[variable])));
        }
        return new MultivariatePolynomial<>(domain, ordering, nVariables - 1, newData);
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
     * @see #evaluate(int, long)
     */
    public MultivariatePolynomial<E> eliminate(int variable, long value) {
        return eliminate(variable, domain.valueOf(value));
    }

    private static final int SIZE_OF_POWERS_CACHE = 32;

    /** cached powers used to save some time */
    static final class PrecomputedPowers<E> {
        private final E value;
        private final Domain<E> domain;
        private final E[] precomputedPowers;

        PrecomputedPowers(E value, Domain<E> domain) {
            this(SIZE_OF_POWERS_CACHE, value, domain);
        }

        PrecomputedPowers(int cacheSize, E value, Domain<E> domain) {
            this.value = domain.valueOf(value);
            this.domain = domain;
            this.precomputedPowers = domain.createArray(cacheSize);
        }

        E pow(int exponent) {
            if (exponent >= SIZE_OF_POWERS_CACHE)
                return domain.pow(value, exponent);

            if (precomputedPowers[exponent] != null)
                return precomputedPowers[exponent];

            E result = domain.getOne();
            E k2p = value;
            int rExp = 0, kExp = 1;
            for (; ; ) {
                if ((exponent&1) != 0)
                    precomputedPowers[rExp += kExp] = result = domain.multiply(result, k2p);
                exponent = exponent >> 1;
                if (exponent == 0)
                    return precomputedPowers[rExp] = result;
                precomputedPowers[kExp *= 2] = k2p = domain.multiply(k2p, k2p);
            }
        }
    }

    /** holds an array of precomputed powers */
    static final class PrecomputedPowersHolder<E> {
        private final int cacheSize;
        private final Domain<E> domain;
        private final PrecomputedPowers<E>[] powers;

        PrecomputedPowersHolder(E[] points, Domain<E> domain) {
            this(SIZE_OF_POWERS_CACHE, points, domain);
        }

        @SuppressWarnings("unchecked")
        PrecomputedPowersHolder(int cacheSize, E[] points, Domain<E> domain) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = new PrecomputedPowers[points.length];
            for (int i = 0; i < points.length; i++)
                powers[i] = points[i] == null ? null : new PrecomputedPowers<E>(cacheSize, points[i], domain);
        }

        PrecomputedPowersHolder(int cacheSize, Domain<E> domain, PrecomputedPowers<E>[] powers) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = powers;
        }

        void set(int i, E point) {
            if (powers[i] == null || !powers[i].value.equals(point))
                powers[i] = new PrecomputedPowers<E>(cacheSize, point, domain);
        }

        E pow(int i, int exponent) {
            return powers[i].pow(exponent);
        }
    }

    @Override
    public MultivariatePolynomial<E> negate() {
        for (Entry<DegreeVector, E> entry : data.entrySet())
            entry.setValue(domain.negate(entry.getValue()));
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> add(MultivariatePolynomial<E> oth) {
        if (data == oth.data)
            return multiply(2);
        ensureCompatible(oth);
        if (oth.isZero())
            return this;
        oth.data.entrySet().forEach(othElement -> add(data, othElement.getKey(), othElement.getValue()));
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> subtract(MultivariatePolynomial<E> oth) {
        if (data == oth.data)
            return toZero();
        ensureCompatible(oth);
        oth.data.entrySet().forEach(othElement -> subtract(data, othElement.getKey(), othElement.getValue()));
        release();
        return this;
    }

    /**
     * Adds {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this + oth}
     */
    MultivariatePolynomial<E> add(DegreeVector dv, E term) {
        ensureCompatible(dv);
        term = domain.valueOf(term);
        add(data, dv, term);
        release();
        return this;
    }

    /**
     * Adds terms to this polynomial and returns it
     *
     * @param vectors degree vectors of the terms
     * @param factors terms coefficients
     * @return {@code this + terms}
     */
    MultivariatePolynomial<E> add(DegreeVector[] vectors, E[] factors) {
        if (factors.length != vectors.length)
            throw new IllegalArgumentException();
        if (factors.length == 0)
            throw new IllegalArgumentException("empty");
        for (int i = 0; i < factors.length; i++) {
            E f = factors[i];
            add(data, vectors[i], f);
        }
        return this;
    }

    /**
     * Adds {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this + oth}
     */
    MultivariatePolynomial<E> add(Entry<DegreeVector, E> term) {
        return add(term.getKey(), term.getValue());
    }

    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial<E> add(E oth) {
        oth = domain.valueOf(oth);
        if (domain.isZero(oth))
            return this;
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
    MultivariatePolynomial<E> subtract(Entry<DegreeVector, E> term) {
        ensureCompatible(term);
        subtract(data, term.getKey(), term.getValue());
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public MultivariatePolynomial<E> subtract(E oth) {
        oth = domain.valueOf(oth);
        if (domain.isZero(oth))
            return this;
        subtract(data, zeroDegreeVector(nVariables), oth);
        release();
        return this;
    }

    /**
     * Removes the leading term from this polynomial
     *
     * @return this - this.lt()
     */
    public MultivariatePolynomial<E> subtractLt() {
        data.pollLastEntry();
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> increment() {
        return add(domain.getOne());
    }

    @Override
    public MultivariatePolynomial<E> decrement() {
        return subtract(domain.getOne());
    }

    private void add(TreeMap<DegreeVector, E> map, DegreeVector deg, E val) {
        add(map, deg, val, domain);
    }

    private static <E> void add(TreeMap<DegreeVector, E> map, DegreeVector deg, E val, Domain<E> domain) {
        if (domain.isZero(val))
            return;
        map.compute(deg, (thisVector, thisValue) -> {
            if (thisValue == null)
                return val;
            E r = domain.add(thisValue, val);
            return domain.isZero(r) ? null : r;
        });
    }

    private void subtract(TreeMap<DegreeVector, E> map, DegreeVector deg, E val) {
        subtract(map, deg, val, domain);
    }

    private static <E> void subtract(TreeMap<DegreeVector, E> map, DegreeVector deg, E val, Domain<E> domain) {
        map.compute(deg, (thisVector, thisValue) -> {
            if (thisValue == null)
                return domain.negate(val);
            E r = domain.subtract(thisValue, val);
            return domain.isZero(r) ? null : r;
        });
    }

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public MultivariatePolynomial<E> multiply(E factor) {
        factor = domain.valueOf(factor);
        if (domain.isOne(factor))
            return this;
        if (domain.isZero(factor))
            return toZero();
        for (Entry<DegreeVector, E> entry : data.entrySet())
            entry.setValue(domain.multiply(entry.getValue(), factor));
        release();
        return this;
    }

    private MultivariatePolynomial<E> loadFrom(TreeMap<DegreeVector, E> map) {
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
    public MultivariatePolynomial<E> multiply(Entry<DegreeVector, E> term) {
        return multiply(term.getKey(), term.getValue());
    }

    /**
     * Multiplies {@code this} by the {@code term}
     *
     * @param dv   degree vector
     * @param term the factor
     * @return {@code} this multiplied by the {@code term}
     */
    public MultivariatePolynomial<E> multiply(DegreeVector dv, E term) {
        ensureCompatible(dv);
        if (dv.isZeroVector())
            return multiply(term);
        if (domain.isZero(term))
            return toZero();

        TreeMap<DegreeVector, E> newMap = new TreeMap<>(ordering);
        for (Entry<DegreeVector, E> thisElement : data.entrySet()) {
            E m = domain.multiply(thisElement.getValue(), term);
            if (!domain.isZero(m))
                newMap.put(thisElement.getKey().multiply(dv), m);
        }

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> multiply(long factor) {
        return multiply(domain.valueOf(factor));
    }

    @Override
    public MultivariatePolynomial<E> multiply(MultivariatePolynomial<E> oth) {
        ensureCompatible(oth);
        TreeMap<DegreeVector, E> newMap = new TreeMap<>(ordering);
        for (Entry<DegreeVector, E> othElement : oth.data.entrySet())
            for (Entry<DegreeVector, E> thisElement : data.entrySet())
                add(newMap, thisElement.getKey().multiply(othElement.getKey()), domain.multiply(thisElement.getValue(), othElement.getValue()));

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> square() {
        return multiply(this);
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> clone() {
        return new MultivariatePolynomial<>(domain, ordering, nVariables, (TreeMap) data.clone());
    }

    /**
     * Returns skeleton of this poly
     *
     * @return skeleton of this poly
     */
    public Set<DegreeVector> getSkeleton() {
        return data.keySet();
    }

    /**
     * Returns skeleton of this poly with respect to specified {@code variables}
     *
     * @param variables the variables
     * @return skeleton of this poly with respect to specified {@code variables}
     */
    public Set<DegreeVector> getSkeleton(int... variables) {
        return data.keySet().stream().map(dv -> dv.of(variables)).collect(Collectors.toSet());
    }

    /**
     * Returns skeleton of this poly with respect to all except specified {@code variables}
     *
     * @param variables the variables to exclude
     * @return skeleton of this poly with respect to all except specified {@code variables}
     */
    public Set<DegreeVector> getSkeletonExcept(int... variables) {
        return data.keySet().stream().map(dv -> dv.except(variables)).collect(Collectors.toSet());
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton
     *
     * @param oth other multivariate polynomial
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton and {@code false} otherwise
     */
    public boolean sameSkeleton(MultivariatePolynomial<E> oth) {
        return getSkeleton().equals(oth.getSkeleton());
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to test
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables} and {@code false} otherwise
     */
    public boolean sameSkeleton(MultivariatePolynomial<E> oth, int... variables) {
        return getSkeleton(variables).equals(oth.getSkeleton(variables));
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect all except specified {@code variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to exclude
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to all except specified  {@code variables} and {@code false} otherwise
     */
    public boolean sameSkeletonExcept(MultivariatePolynomial<E> oth, int... variables) {
        return getSkeletonExcept(variables).equals(oth.getSkeletonExcept(variables));
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

    private static final Pattern nonTrivialCoefficientString = Pattern.compile("[\\+\\-\\*]");

    private static <E> String coeffToString(E coeff) {
        String cfs = coeff.toString();
        if (coeff instanceof BigInteger)
            return cfs;
        if (nonTrivialCoefficientString.matcher(cfs).find())
            return "(" + cfs + ")";
        else return cfs;
    }

    public String toString(String... vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Entry<DegreeVector, E> term : data.entrySet()) {
            DegreeVector monomial = term.getKey();
            E coeff = term.getValue();
            if (domain.isZero(coeff))
                continue;
            String monomialString = monomial.toString(vars);
            if (first) {
                if (!domain.isOne(coeff) || monomialString.isEmpty()) {
                    sb.append(coeffToString(coeff));
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
                first = false;
            } else {
                if (domain.signum(coeff) > 0)
                    sb.append("+");
                else {
                    sb.append("-");
                    coeff = domain.negate(coeff);
                }

                if (!domain.isOne(coeff) || monomialString.isEmpty()) {
                    sb.append(coeffToString(coeff));
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
        private static final DegreeVector EMPTY = new DegreeVector(new int[0]);
        final int[] exponents;
        final int totalDegree;

        public DegreeVector(int nVariables, int position, int exponent) {
            this.exponents = new int[nVariables];
            this.exponents[position] = exponent;
            this.totalDegree = exponent;
        }

        private DegreeVector(int[] exponents, int totalDegree) {
            this.exponents = exponents;
            this.totalDegree = totalDegree;
        }

        public DegreeVector(int[] exponents) {
            this.exponents = exponents;
            this.totalDegree = ArraysUtil.sum(exponents);
        }

        public DegreeVector joinNewVariable() {
            return new DegreeVector(Arrays.copyOf(exponents, exponents.length + 1), totalDegree);
        }

        public DegreeVector of(int[] variables) {
            int[] exs = new int[exponents.length];
            int totalDegree = 0;
            for (int i : variables) {
                exs[i] = exponents[i];
                totalDegree += exs[i];
            }
            return new DegreeVector(exs, totalDegree);
        }

        public DegreeVector except(int[] variables) {
            int[] exs = exponents.clone();
            int totalDegree = this.totalDegree;
            for (int i : variables) {
                exs[i] = 0;
                totalDegree -= exponents[i];
            }
            return new DegreeVector(exs, totalDegree);
        }

        public boolean isZeroVector() {
            return totalDegree == 0;
        }

        private static String toString0(String var, int exp) {
            return exp == 0 ? "" : var + (exp == 1 ? "" : "^" + exp);
        }

        public DegreeVector multiply(DegreeVector dv) {
            int[] newExponents = new int[exponents.length];
            for (int i = 0; i < exponents.length; i++)
                newExponents[i] = exponents[i] + dv.exponents[i];
            return new DegreeVector(newExponents, totalDegree + dv.totalDegree);
        }

        public DegreeVector divide(DegreeVector dv) {
            int[] newExponents = new int[exponents.length];
            for (int i = 0; i < exponents.length; i++) {
                newExponents[i] = exponents[i] - dv.exponents[i];
                if (newExponents[i] < 0)
                    return null;
            }
            return new DegreeVector(newExponents, totalDegree - dv.totalDegree);
        }

        DegreeVector without(int i) {
            if (exponents.length == 1) {
                assert i == 0;
                return EMPTY;
            }
            return new DegreeVector(ArraysUtil.remove(exponents, i), totalDegree - exponents[i]);
        }

        DegreeVector setZero(int i) {
            if (exponents.length == 1) {
                assert i == 0;
                return EMPTY;
            }
            int[] newExponents = exponents.clone();
            newExponents[i] = 0;
            return new DegreeVector(newExponents, totalDegree - exponents[i]);
        }

        DegreeVector set(int i, int exponent) {
            int[] newExponents = exponents.clone();
            newExponents[i] = exponent;
            return new DegreeVector(newExponents, totalDegree - exponents[i] + exponent);
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

    static final class EntryImpl<V> implements Entry<DegreeVector, V> {
        final DegreeVector degreeVector;
        final V coefficient;

        public EntryImpl(DegreeVector degreeVector, V coefficient) {
            this.degreeVector = degreeVector;
            this.coefficient = coefficient;
        }

        @Override
        public DegreeVector getKey() {
            return degreeVector;
        }

        @Override
        public V getValue() {
            return coefficient;
        }

        @Override
        public V setValue(V value) {
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
