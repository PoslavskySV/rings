package cc.r2.core.poly.multivar;

import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Collections;
import java.util.Comparator;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
abstract class AMultivariatePolynomial<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements IGeneralPolynomial<Poly> {
    /** number of variables */
    final int nVariables;
    /** the ordering */
    final Comparator<DegreeVector> ordering;
    /** the actual data */
    final MonomialsSet<Term> terms;
    @SuppressWarnings("unchecked")
    final Poly self = (Poly) this;

    AMultivariatePolynomial(int nVariables, Comparator<DegreeVector> ordering, MonomialsSet<Term> terms) {
        this.nVariables = nVariables;
        this.ordering = ordering;
        this.terms = terms;
    }

    /**
     * Renames variable {@code i} to {@code j} and {@code j} to {@code i}
     *
     * @param poly the polynomial
     * @param i    the first variable
     * @param j    the second variable
     * @return polynomial with variable {@code i} renamed to {@code j} and {@code j} renamed to {@code i}
     */
    static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>> P
    swapVariables(P poly, int i, int j) {
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
    static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>> P
    renameVariables(P poly, int[] newVariables) {
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
    public static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>> P
    renameVariables(P poly, int[] newVariables, Comparator<DegreeVector> newOrdering) {
        // NOTE: always return a copy of poly, even if order of variables is unchanged
        MonomialsSet<T> data = new MonomialsSet<>(newOrdering);
        for (T e : poly.terms)
            data.add(e.setDegreeVector(map(e.exponents, newVariables), e.totalDegree));
        return poly.create(data);
    }

    private static int[] map(int[] degrees, int[] mapping) {
        int[] newDegrees = new int[degrees.length];
        for (int i = 0; i < degrees.length; i++)
            newDegrees[i] = degrees[mapping[i]];
        return newDegrees;
    }

    abstract IUnivariatePolynomial asUnivariate();

    /* private factory */
    final Poly create(MonomialsSet<Term> terms) {
        return create(nVariables, terms);
    }

    /* private factory */
    abstract Poly create(int nVariables, MonomialsSet<Term> terms);

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param terms the monomial terms
     * @return multivariate polynomial
     */
    public final Poly create(Term... terms) {
        MonomialsSet<Term> monomials = new MonomialsSet<>(ordering);
        for (Term term : terms) {
            if (term.exponents.length != nVariables)
                throw new IllegalArgumentException();
            add(this.terms, term);
        }
        return create(monomials);
    }

    /**
     * Creates multivariate polynomial from with a single monomial terms
     *
     * @param term the monomial terms
     * @return multivariate polynomial
     */
    public final Poly create(Term term) {
        if (term.exponents.length != nVariables)
            throw new IllegalArgumentException();
        MonomialsSet<Term> monomials = new MonomialsSet<>(ordering);
        monomials.add(term);
        return create(monomials);
    }

    /**
     * Copies this with the new ordering {@code newOrdering}
     *
     * @param newOrdering the new ordering
     * @return a copy of this with a new ordering
     */
    public final Poly setOrdering(Comparator<DegreeVector> newOrdering) {
        MonomialsSet<Term> newData = new MonomialsSet<>(newOrdering);
        newData.putAll(terms);
        return create(newData);
    }

    /** release caches */
    protected void release() { /* add cache in the future */ }

    /**
     * Returns the number of terms in this polynomial
     *
     * @return number of terms
     */
    public final int size() {return terms.size();}

    @Override
    public final boolean isMonomial() {
        return size() <= 1;
    }

    @Override
    public final Poly toZero() {
        terms.clear();
        release();
        return self;
    }

    @Override
    public final Poly set(Poly oth) {
        checkSameDomainWith(oth);
        return loadFrom(oth.terms);
    }

    final Poly loadFrom(MonomialsSet<Term> map) {
        terms.clear();
        terms.putAll(map);
        release();
        return self;
    }

    /**
     * Returns a copy of this with {@code nVariables = nVariables + 1}
     *
     * @return a copy of this with one additional (last) variable added
     */
    public final Poly joinNewVariable() {
        MonomialsSet<Term> newData = new MonomialsSet<>(ordering);
        for (Term term : terms)
            newData.add(term.joinNewVariable());
        return create(nVariables + 1, newData);
    }

    /**
     * Returns the number of really used variables
     *
     * @return the number of presenting variables
     */
    public final int nUsedVariables() {
        int[] degrees = degrees();
        int r = 0;
        for (int d : degrees)
            if (d != 0)
                ++r;
        return r;
    }

    /**
     * Returns the total degree of this polynomial, that is the maximal total degree among all terms
     *
     * @return the total degree of this polynomial, that is the maximal total degree among all terms
     */
    @Override
    public int degree() {
        int max = 0;
        for (Term db : terms)
            max = Math.max(max, db.totalDegree);
        return max;
    }

    /**
     * Returns the degree of this polynomial with respect to i-th variable
     *
     * @return the degree of this polynomial with respect to i-th variable
     */
    public final int degree(int i) {
        int max = 0;
        for (Term db : terms)
            max = Math.max(max, db.exponents[i]);
        return max;
    }

    /**
     * Returns an array of degrees of all variables, so that is i-th element of the result is the polynomial degree
     * with respect to i-th variable
     *
     * @return array of degrees
     */
    public final int[] degrees() {
        int[] degrees = new int[nVariables];
        for (Term db : terms)
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
    public final int[] degrees(int variable) {
        TIntHashSet degrees = new TIntHashSet();
        for (Term db : terms)
            degrees.add(db.exponents[variable]);
        return degrees.toArray();
    }

    /**
     * Returns the sum of {@link #degrees()}
     *
     * @return sum of {@link #degrees()}
     */
    public final int degreeSum() {
        int r = 0;
        for (int d : degrees())
            r += d;
        return r;
    }

    /**
     * Returns whether this poly is effectively univariate
     *
     * @return whether this poly is effectively univariate
     */
    public final boolean isEffectiveUnivariate() {
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
     * Returns a coefficient before {@code variable^exponent} as a multivariate polynomial
     *
     * @param variable the variable
     * @param exponent the exponent
     * @return coefficient before {@code variable^exponent} as a multivariate polynomial
     */
    public final Poly coefficientOf(int variable, int exponent) {
        Poly result = createZero();
        for (Term e : terms) {
            if (e.exponents[variable] != exponent)
                continue;
            result.add(e.setZero(variable));
        }
        return result;
    }

    /**
     * Returns the leading term in this polynomial according to ordering
     *
     * @return the leading term in this polynomial according to ordering
     */
    public abstract Term lt();

    /**
     * Returns the monomial content of this polynomial
     *
     * @return the monomial content of this polynomial
     */
    //todo rename!!
    public final Term monomialContent() {
        return commonContent(null);
    }

    /**
     * Returns common content of {@code this} and {@code monomial}
     *
     * @param monomial the monomial
     * @return common monomial factor of {@code this} and {@code monomial}
     */
    final Term commonContent(Term monomial) {
        int[] exponents = monomial == null ? null : monomial.exponents.clone();
        for (Term degreeVector : terms)
            if (exponents == null)
                exponents = degreeVector.exponents.clone();
            else
                setMin(degreeVector.exponents, exponents);
        if (exponents == null)
            return getUnitTerm();
        return createTermWithUnitCoefficient(exponents);
    }

    static void setMin(int[] dv, int[] exponents) {
        for (int i = 0; i < exponents.length; ++i)
            if (dv[i] < exponents[i])
                exponents[i] = dv[i];
    }

    /** private term factory */
    abstract Term createTermWithUnitCoefficient(int[] exponents);

    /** private term factory */
    abstract Term getUnitTerm();

    /**
     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code monomial}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param monomial monomial degrees
     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
     */
    public final Poly divideDegreeVectorOrNull(DegreeVector monomial) {
        if (monomial.isZeroVector())
            return self;
        MonomialsSet<Term> map = new MonomialsSet<>(ordering);
        for (Term term : terms) {
            Term dv = term.divide(monomial);
            if (dv == null)
                return null;
            map.add(dv);
        }
        return loadFrom(map);
    }

    /** check whether number of variables is the same */
    final void checkSameDomainWith(Term oth) {
        if (nVariables != oth.exponents.length)
            throw new IllegalArgumentException("Combining multivariate polynomials from different fields");
    }

    /**
     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code monomial}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param monomial monomial degrees
     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
     */
    public abstract Poly divideOrNull(Term monomial);

    /** add term to polynomial represented as terms */
    abstract void add(MonomialsSet<Term> terms, Term term);

    /** subtract term from polynomial represented as terms */
    abstract void subtract(MonomialsSet<Term> terms, Term term);

    @Override
    public final Poly add(Poly oth) {
        if (terms == oth.terms)
            return multiply(2);
        checkSameDomainWith(oth);
        if (oth.isZero())
            return self;
        for (Term term : oth.terms)
            add(terms, term);
        release();
        return self;
    }

    @Override
    public final Poly subtract(Poly oth) {
        if (terms == oth.terms)
            return toZero();
        checkSameDomainWith(oth);
        if (oth.isZero())
            return self;
        for (Term term : oth.terms)
            subtract(terms, term);
        release();
        return self;
    }

    /**
     * Adds {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this + oth}
     */
    public final Poly add(Term term) {
        checkSameDomainWith(term);
        add(terms, term);
        release();
        return self;
    }

    /**
     * Adds {@code term} to this polynomial and returns it
     *
     * @param term some term
     * @return {@code this + oth}
     */
    public final Poly subtract(Term term) {
        checkSameDomainWith(term);
        subtract(terms, term);
        release();
        return self;
    }

    /**
     * Adds terms to this polynomial and returns it
     *
     * @param terms terms
     * @return {@code this + terms}
     */
    public final Poly add(Term... terms) {
        if (terms.length == 0)
            throw new IllegalArgumentException("empty");
        for (Term term : terms)
            add(term);
        return self;
    }

    /**
     * Removes the leading term from this polynomial
     *
     * @return this - this.lt()
     */
    public final Poly subtractLt() {
        terms.pollLastEntry();
        release();
        return self;
    }

    /**
     * Multiplies {@code this} by the {@code term}
     *
     * @param term the factor
     * @return {@code} this multiplied by the {@code term}
     */
    public abstract Poly multiply(Term term);

    /**
     * Returns skeleton of this poly
     *
     * @return skeleton of this poly
     */
    public final Set<DegreeVector> getSkeleton() {
        return Collections.unmodifiableSet(terms.keySet());
    }

    /**
     * Returns skeleton of this poly with respect to specified {@code variables}
     *
     * @param variables the variables
     * @return skeleton of this poly with respect to specified {@code variables}
     */
    public final Set<DegreeVector> getSkeleton(int... variables) {
        return terms.keySet().stream().map(dv -> dv.of(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
    }

    /**
     * Returns skeleton of this poly with respect to all except specified {@code variables}
     *
     * @param variables the variables to exclude
     * @return skeleton of this poly with respect to all except specified {@code variables}
     */
    public final Set<DegreeVector> getSkeletonExcept(int... variables) {
        return terms.keySet().stream().map(dv -> dv.except(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton
     *
     * @param oth other multivariate polynomial
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton and {@code false} otherwise
     */
    public final boolean sameSkeleton(Poly oth) {
        return getSkeleton().equals(oth.getSkeleton());
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to test
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables} and {@code false} otherwise
     */
    public final boolean sameSkeleton(Poly oth, int... variables) {
        return getSkeleton(variables).equals(oth.getSkeleton(variables));
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect all except specified {@code variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to exclude
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to all except specified  {@code variables} and {@code false} otherwise
     */
    public final boolean sameSkeletonExcept(Poly oth, int... variables) {
        return getSkeletonExcept(variables).equals(oth.getSkeletonExcept(variables));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AMultivariatePolynomial<?, ?> that = (AMultivariatePolynomial<?, ?>) o;

        if (nVariables != that.nVariables) return false;
        return terms.equals(that.terms);
    }

    @Override
    public int hashCode() {
        int result = nVariables;
        result = 31 * result + terms.hashCode();
        return result;
    }

    @Override
    public abstract Poly clone();

    /**
     * Evaluates {@code poly} at random point chosen in such way that the skeleton of evaluated version is the same as of the
     * original {@code poly} with respect to all except {@code variable} variables
     */
    abstract Poly evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd);
}
