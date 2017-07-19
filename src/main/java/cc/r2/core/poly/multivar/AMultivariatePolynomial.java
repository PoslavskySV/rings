package cc.r2.core.poly.multivar;

import cc.r2.core.poly.IGeneralPolynomial;
import cc.r2.core.poly.MultivariatePolynomials;
import cc.r2.core.poly.univar.IUnivariatePolynomial;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;
import java.util.stream.Collectors;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public abstract class AMultivariatePolynomial<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements IGeneralPolynomial<Poly>, Iterable<Term> {
    /** number of variables */
    final int nVariables;
    /** the ordering */
    final Comparator<DegreeVector> ordering;
    /** the actual data */
    final MonomialsSet<Term> terms;
    @SuppressWarnings("unchecked")
    private final Poly self = (Poly) this;

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
        if (i == j)
            return poly.clone();
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
    static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>>
    P renameVariables(P poly, int[] newVariables) {
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

    @SuppressWarnings("unchecked")
    public static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly asMultivariate(IUnivariatePolynomial poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        if (poly instanceof UnivariatePolynomial)
            return (Poly) MultivariatePolynomial.asMultivariate((UnivariatePolynomial) poly, nVariables, variable, ordering);
        else if (poly instanceof lUnivariatePolynomialZp)
            return (Poly) lMultivariatePolynomialZp.asMultivariate((lUnivariatePolynomialZp) poly, nVariables, variable, ordering);
        else
            throw new RuntimeException();
    }

    @SuppressWarnings("unchecked")
    public static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] asMultivariate(IUnivariatePolynomial[] polys, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        Poly p = asMultivariate(polys[0], nVariables, variable, ordering);
        Poly[] r = p.arrayNewInstance(polys.length);
        r[0] = p;
        for (int i = 1; i < polys.length; i++)
            r[i] = asMultivariate(polys[i], nVariables, variable, ordering);
        return r;
    }

    @SuppressWarnings("unchecked")
    public static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> asMultivariate(List<IUnivariatePolynomial> polys, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        List<Poly> r = polys.stream().map(p -> AMultivariatePolynomial.<Term, Poly>asMultivariate(p, nVariables, variable, ordering)).collect(Collectors.toList());
        return r;
    }

    public abstract IUnivariatePolynomial asUnivariate();

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
        add(monomials, term);
        return create(monomials);
    }

    /**
     * Creates variable^degree monomial
     *
     * @param variable the variable
     * @param degree   the monomial degree
     * @return variable^degree
     */
    public final Poly createUnivariateMonomial(int variable, int degree) {
        int[] degreeVector = new int[nVariables];
        degreeVector[variable] = degree;
        return create(createTermWithUnitCoefficient(degreeVector));
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
     * @return the number of terms
     */
    @Override
    public final int size() {return terms.size();}

    @Override
    public final Iterator<Term> iterator() {
        return terms.iterator();
    }

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

    /** Remove specified variable */
    public final Poly dropVariable(int variable, boolean eliminate) {
        MonomialsSet<Term> newData = new MonomialsSet<>(ordering);
        for (Term term : terms)
            newData.add(term.without(variable));
        return create(eliminate ? nVariables - 1 : nVariables, newData);
    }

    /** Insert variable */
    public final Poly insertVariable(int variable) {
        MonomialsSet<Term> newData = new MonomialsSet<>(ordering);
        for (Term term : terms)
            newData.add(term.insert(variable));
        return create(nVariables + 1, newData);
    }

    /** auxiliary method */
    final Poly setNVariables(int newNVariables) {
        if (newNVariables == nVariables)
            return self;
        MonomialsSet<Term> newData = new MonomialsSet<>(ordering);
        for (Term term : terms)
            newData.add(term.setDegreeVector(Arrays.copyOf(term.exponents, newNVariables)));
        return create(newNVariables, newData);
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

    public final Poly joinNewVariables(int newNVariables, int[] mapping) {
        MonomialsSet<Term> newData = new MonomialsSet<>(ordering);
        for (Term term : terms)
            newData.add(term.joinNewVariables(newNVariables, mapping));
        return create(newNVariables, newData);
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
     * Returns the multidegree of this polynomial i.e. exponents of the leading term (without copying)
     *
     * @return the multidegree of this polynomial i.e. exponents of the leading term (without copying)
     */
    public final int[] multidegree() {
        return lt().exponents;
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
        return ArraysUtil.sum(degrees());
    }

    /**
     * Returns whether this poly is effectively univariate
     *
     * @return whether this poly is effectively univariate
     */
    public final boolean isEffectiveUnivariate() {
        return univariateVariable() != -1;
    }

    /**
     * Returns -1 if this poly is not effectively univariate or variable in which it is univariate
     *
     * @return -1 if this poly is not effectively univariate or variable in which it is univariate
     */
    public final int univariateVariable() {
        if (isConstant())
            return 0;
        if (nVariables == 1)
            return 0;
        int[] degrees = degrees();
        int var = -1;
        for (int i = 0; i < nVariables; i++) {
            if (degrees[i] != 0) {
                if (var != -1)
                    return -1;
                else
                    var = i;
            }
        }
        return var;
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
     * Converts this polynomial to a univariate polynomial over specified variable with the multivariate coefficient
     * domain.
     *
     * @param variable variable which will be treated as univariate variable
     * @return univariate polynomial
     * @throws IllegalArgumentException if this is not effectively a univariate polynomial
     */
    public final UnivariatePolynomial<Poly> asUnivariate(int variable) {
        MultivariatePolynomials<Poly> domain = new MultivariatePolynomials<>(self);
        Poly[] univarData = domain.createZeroesArray(degree(variable) + 1);
        for (Term e : terms)
            univarData[e.exponents[variable]].add(e.set(variable, 0));
        return UnivariatePolynomial.createUnsafe(domain, univarData);
    }

    /**
     * Converts this to a multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     *
     * @param variable variable
     * @return multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     */
    public abstract MultivariatePolynomial<? extends IUnivariatePolynomial> asOverUnivariate(int variable);

    /**
     * Converts this to a multivariate polynomial with coefficients being univariate polynomials over {@code variable},
     * the resulting polynomial have (nVariable - 1) multivariate variables
     *
     * @param variable variable
     * @return multivariate polynomial with coefficients being univariate polynomials over {@code variable}, the
     * resulting polynomial have (nVariable - 1) multivariate variables
     */
    public abstract MultivariatePolynomial<? extends IUnivariatePolynomial> asOverUnivariateEliminate(int variable);

    /**
     * Converts this to a multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     *
     * @param variables the variables to separate
     * @return multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     */
    public abstract MultivariatePolynomial<Poly> asOverMultivariate(int... variables);

    /**
     * Converts this to a multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     *
     * @param variables the variables to separate
     * @return multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     */
    public abstract MultivariatePolynomial<Poly> asOverMultivariateEliminate(int... variables);

    /**
     * Convert univariate polynomial over multivariate polynomials to a normal multivariate poly
     *
     * @param uPoly    univariate polynomial over multivariate polynomials
     * @param variable the univariate variable
     * @return multivariate poly
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly asMultivariate(UnivariatePolynomial<Poly> uPoly, int variable) {
        Poly result = uPoly.domain.getZero();
        for (int i = uPoly.degree(); i >= 0; --i) {
            if (uPoly.isZeroAt(i))
                continue;
            result.add(result.createUnivariateMonomial(variable, i).multiply(uPoly.get(i)));
        }
        return result;
    }

    /**
     * Gives primitive part of this considered over univariate polynomial domain R[variable]
     *
     * @param variable the variable
     * @return primitive part of this considered over univariate polynomial domain R[variable]
     */
    public abstract Poly primitivePart(int variable);

    /**
     * Gives the content of this considered over univariate polynomial domain R[variable]
     *
     * @param variable the variable
     * @return the content of this considered over univariate polynomial domain R[variable]
     */
    public abstract IUnivariatePolynomial contentUnivariate(int variable);

    /**
     * Gives the content of this considered over univariate polynomial domain R[variable]
     *
     * @param variable the variable
     * @return the content of this considered over univariate polynomial domain R[variable]
     */
    @SuppressWarnings("unchecked")
    public final Poly content(int variable) {
        return asMultivariate(contentUnivariate(variable), nVariables, variable, ordering);
    }

    /**
     * Gives the content of this considered as R[x1, ... (except variable) ..., xN][variable]
     *
     * @param variable the variable
     * @return the content of this considered as R[x1, ... (except variable) ..., xN][variable]
     */
    @SuppressWarnings("unchecked")
    public final Poly contentExcept(int variable) {
        return asUnivariate(variable).content();
    }

    /**
     * Multiplies this by variable^exponent
     *
     * @param variable the variable
     * @param exponent the exponent
     * @return this multiplied by variable^exponent
     */
    public final Poly multiplyByMonomial(int variable, int exponent) {
        if (exponent == 0)
            return self;
        Collection<Term> oldData = new ArrayList<>(terms.values());

        terms.clear();
        for (Term term : oldData)
            terms.add(term.set(variable, term.exponents[variable] + exponent));

        return self;
    }

    /**
     * Returns the multivariate leading coefficient of this poly seen as univariate poly over specified variable.
     *
     * @param variable the variable
     * @return multivariate leading coefficient of this poly treated as univariate poly over specified variable
     */
    public final Poly lc(int variable) {
        int degree = degree(variable);
        Poly result = createZero();
        for (Term term : this)
            if (term.exponents[variable] == degree)
                result.add(term.set(variable, 0));

        return result;
    }

    /**
     * Set the leading coefficient of specified variable to a specified value
     *
     * @param variable the variable
     * @param lc       the leading coefficient
     */
    public final Poly setLC(int variable, Poly lc) {
        int degree = degree(variable);

        lc = lc.clone().multiplyByMonomial(variable, degree);
        Iterator<Map.Entry<DegreeVector, Term>> it = terms.entrySet().iterator();
        while (it.hasNext()) {
            Term term = it.next().getValue();
            if (term.exponents[variable] == degree)
                it.remove();
        }
        terms.putAll(lc.terms);
        return self;
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

    /** private term factory */
    abstract Term getZeroTerm();

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
     * Multiplies {@code this} by the {@code term}
     *
     * @param term the factor
     * @return {@code} this multiplied by the {@code term}
     */
    public final Poly multiplyByDegreeVector(DegreeVector term) {
        return multiply(getUnitTerm().setDegreeVector(term));
    }

    /** Whether monomial is zero */
    abstract boolean isZeroMonomial(Term a);

    /** Multiply two terms */
    abstract Term multiply(Term a, Term b);

    /** Divide a/b or return null if exact division is not possible */
    abstract Term divideOrNull(Term a, Term b);

    /**
     * Returns skeleton of this poly
     *
     * @return skeleton of this poly
     */
    public final Set<DegreeVector> getSkeleton() {
        return Collections.unmodifiableSet(terms.keySet());
    }

    public final Poly setAllCoefficientsToUnit() {
        Term unit = getUnitTerm();
        for (Map.Entry<DegreeVector, Term> entry : terms.entrySet())
            entry.setValue(unit.setDegreeVector(entry.getKey()));
        return self;
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

    /**
     * Gives partial derivative with respect to specified variable
     *
     * @param variable the variable
     * @return partial derivative with respect to specified variable
     */
    public final Poly derivative(int variable) {return derivative(variable, 1);}

    /**
     * Gives partial derivative of specified {@code order} with respect to specified variable
     *
     * @param variable the variable
     * @param order    derivative order
     * @return partial derivative of specified {@code order} with respect to specified variable
     */
    public abstract Poly derivative(int variable, int order);

    /**
     * Gives (unevaluated) coefficient of Taylor series expansion for specified variable that
     * is {@code derivative(poly, variable, order) / order! }, where the derivative is formal derivative and
     * calculated with arithmetic performed in Z domain (to overcome possible zeros in Zp).
     *
     * @param variable the variable
     * @param order    derivative order
     * @return {@code derivative(poly, variable, order) / order! }, where the derivative is formal derivative and
     * calculated with arithmetic performed in Z domain (to overcome possible zeros in Zp)
     */
    public abstract Poly seriesCoefficient(int variable, int order);

    /**
     * Gives the derivative vector
     *
     * @return derivative vector
     */
    public final Poly[] derivative() {
        Poly[] result = arrayNewInstance(nVariables);
        for (int i = 0; i < nVariables; ++i)
            result[i] = derivative(i);
        return result;
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
     * Evaluates {@code poly} at random point
     */
    abstract Poly evaluateAtRandom(int variable, RandomGenerator rnd);

    /**
     * Evaluates {@code poly} at random point chosen in such way that the skeleton of evaluated version is the same as of the
     * original {@code poly} with respect to all except {@code variable} variables
     */
    abstract Poly evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd);


    /**
     * Collector which collects stream of element to a UnivariatePolynomial
     */
    public static final class PolynomialCollector
            <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
            implements Collector<Term, Poly, Poly> {
        final Supplier<Poly> supplier;
        final BiConsumer<Poly, Term> accumulator = Poly::add;
        final BinaryOperator<Poly> combiner = (l, r) -> {l.add(r); return l;};
        final Function<Poly, Poly> finisher = Function.identity();

        public PolynomialCollector(Supplier<Poly> supplier) {
            this.supplier = supplier;
        }

        @Override
        public Supplier<Poly> supplier() {
            return supplier;
        }

        @Override
        public BiConsumer<Poly, Term> accumulator() {
            return accumulator;
        }

        @Override
        public BinaryOperator<Poly> combiner() {
            return combiner;
        }

        @Override
        public Function<Poly, Poly> finisher() {
            return finisher;
        }

        @Override
        public Set<Characteristics> characteristics() {
            return EnumSet.of(Characteristics.IDENTITY_FINISH);
        }
    }
}
