package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.WithVariables;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.util.ArraysUtil;
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
 * Parent class for multivariate polynomials.
 *
 * <p> <i><b>Variables:</b></i> <p> The number of variables is invariant, which means that binary arithmetic operations
 * on polynomials with different number of variables (see {@link #nVariables}) are prohibited. Of course all exponents
 * of particular variable may be zero, so e.g.
 * <pre>
 * <code>MultivariatePolynomial.parse("x^2 + 2*x*y + y^3", "x", "y", "z")
 * </code></pre>
 * will have nVariables == 3 while "z" is actually absent in the poly.
 *
 * <p> Particular string names of variables are not stored in the polynomial data structure, instead the variables are
 * treated as consequent integer numbers (0, 1, 2,...), where 0 states for the first variable, 1 for the second etc.
 * Information about variables is accessible by the integer index of the variable. The mapping between the variable
 * index and its string representation should be stored separately outside this class. For example:
 * <pre>
 * <code>// x is the first variable, y is the second
 * String[] variables = {"x", "y"};
 * MultivariatePolynomial&lt;BigInteger&gt; poly =
 *                  MultivariatePolynomial.parse("x^2 + 2*x*y + y^3", variables);
 *
 * // degree in x
 * int xDegree = poly.degree(0);
 * assert xDegree == 2;
 * // degree in y
 * int yDegree = poly.degree(1);
 * assert yDegree == 3;
 *
 * // will use the specified mapping and print x^2 + 2*x*y + y^3
 * System.out.println(poly.toString(variables));
 *
 * // will use the default mapping and print a^2 + 2*a*b + b^3
 * System.out.println(poly.toString());
 * </code>
 * </pre>
 *
 * <p> <i><b>Terms storage and ordering:</b></i>
 *
 * <p> Terms of multivariate polynomial are stored in a sorted map {@code DegreeVector -> Monomial} (see {@link
 * MonomialSet}). The order of monomials is defined by the {@code Comparator<DegreeVector>} which possible values are
 * {@link MonomialOrder#LEX}, {@link MonomialOrder#ALEX}, {@link MonomialOrder#GREVLEX} and {@link MonomialOrder#GRLEX}.
 * All operations on the instances of this will preserve the ordering of this. The leading term of the poly is defined
 * with respect to this ordering.
 *
 * @param <Term> type of monomials
 * @param <Poly> type of this (self-type)
 * @see IPolynomial
 * @see MultivariatePolynomialZp64
 * @see MultivariatePolynomial
 * @since 1.0
 */
public abstract class AMultivariatePolynomial<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
        implements IPolynomial<Poly>, Iterable<Term> {
    private static final long serialVersionUID = 1L;
    /** The number of variables */
    public final int nVariables;
    /** The ordering */
    public final Comparator<DegreeVector> ordering;
    /** the actual data */
    final MonomialSet<Term> terms;
    @SuppressWarnings("unchecked")
    private final Poly self = (Poly) this;

    AMultivariatePolynomial(int nVariables, Comparator<DegreeVector> ordering, MonomialSet<Term> terms) {
        this.nVariables = nVariables;
        this.ordering = ordering;
        this.terms = terms;
    }

    /**
     * Renames variable {@code i} to {@code j} and {@code j} to {@code i} (new instance created)
     *
     * @param poly the polynomial
     * @param i    the first variable
     * @param j    the second variable
     * @return polynomial with variable {@code i} renamed to {@code j} and {@code j} renamed to {@code i}
     */
    public static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>> P
    swapVariables(P poly, int i, int j) {
        if (i == j)
            return poly.clone();
        int[] newVariables = ArraysUtil.sequence(poly.nVariables);
        newVariables[i] = j;
        newVariables[j] = i;
        return renameVariables(poly, newVariables, poly.ordering);
    }

    /**
     * Rename variables from [0,1,...N] to [newVariables[0], newVariables[1], ..., newVariables[N]] (new instance
     * created)
     *
     * @param poly         the polynomial
     * @param newVariables the new variables
     * @return renamed polynomial
     */
    public static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>>
    P renameVariables(P poly, int[] newVariables) {
        return renameVariables(poly, newVariables, poly.ordering);
    }

    /**
     * Rename variables from [0,1,...N] to [newVariables[0], newVariables[1], ..., newVariables[N]] (new instance
     * created)
     *
     * @param poly         the polynomial
     * @param newVariables the new variables
     * @param newOrdering  the new ordering
     * @return renamed polynomial
     */
    public static <T extends DegreeVector<T>, P extends AMultivariatePolynomial<T, P>>
    P renameVariables(P poly, int[] newVariables, Comparator<DegreeVector> newOrdering) {
        // NOTE: always return a copy of poly, even if order of variables is unchanged
        MonomialSet<T> data = new MonomialSet<>(newOrdering);
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

    /**
     * Converts univariate polynomial to multivariate. Example:
     * <pre>
     * <code>//convert (x^2 + 1) in Z[x] to multivariate polynomial (c^2 + 1) in Z[a,b,c]
     * multivarPoly = asMultivariate(univarPoly, 3, 2, MonomialOrder.LEX)
     * </code>
     * </pre>
     *
     * @param poly       the univariate polynomial
     * @param nVariables the total number of variables in the result
     * @param variable   the univariate variable
     * @param ordering   the term order
     * @param <Term>     desired terms type
     * @param <Poly>     desired polynomial type
     * @return the multivariate polynomial
     */
    @SuppressWarnings("unchecked")
    public static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly asMultivariate(IUnivariatePolynomial poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        if (poly instanceof UnivariatePolynomial)
            return (Poly) MultivariatePolynomial.asMultivariate((UnivariatePolynomial) poly, nVariables, variable, ordering);
        else if (poly instanceof UnivariatePolynomialZp64)
            return (Poly) MultivariatePolynomialZp64.asMultivariate((UnivariatePolynomialZp64) poly, nVariables, variable, ordering);
        else
            throw new RuntimeException();
    }

    /**
     * Converts this to univariate polynomial or throws exception if conversion is impossible (more than one variable
     * have non zero exponents)
     *
     * @return univariate polynomial
     * @throws IllegalArgumentException if this is not effectively a univariate polynomial
     */
    public abstract IUnivariatePolynomial asUnivariate();

    /* private factory */
    final Poly create(MonomialSet<Term> terms) {
        return create(nVariables, ordering, terms);
    }

    /* private factory */
    final Poly create(int nVariables, MonomialSet<Term> terms) {
        return create(nVariables, ordering, terms);
    }

    /* private factory */
    abstract Poly create(int nVariables, Comparator<DegreeVector> ordering, MonomialSet<Term> terms);

    /**
     * Creates multivariate polynomial over the same ring as this from the list of monomials
     *
     * @param terms the monomials
     * @return multivariate polynomial
     */
    public final Poly create(Term... terms) {
        return create(Arrays.asList(terms));
    }

    /**
     * Creates multivariate polynomial over the same ring as this from the list of monomials
     *
     * @param terms the monomials
     * @return multivariate polynomial
     */
    public final Poly create(Iterable<Term> terms) {
        MonomialSet<Term> monomials = new MonomialSet<>(ordering);
        for (Term term : terms) {
            if (term.exponents.length != nVariables)
                throw new IllegalArgumentException();
            add(this.terms, term);
        }
        return create(monomials);
    }

    /**
     * Creates multivariate polynomial over the same ring as this with the single monomial
     *
     * @param term the monomial
     * @return multivariate polynomial
     */
    public final Poly create(Term term) {
        if (term.exponents.length != nVariables)
            throw new IllegalArgumentException();
        MonomialSet<Term> monomials = new MonomialSet<>(ordering);
        add(monomials, term);
        return create(monomials);
    }

    /**
     * Creates monomial over the same ring as this of the form {@code variable ^ degree}
     *
     * @param variable the variable
     * @param degree   the monomial degree
     * @return monomial {@code variable ^ degree}
     */
    public final Poly createMonomial(int variable, int degree) {
        int[] degreeVector = new int[nVariables];
        degreeVector[variable] = degree;
        return create(createTermWithUnitCoefficient(degreeVector));
    }

    /**
     * Makes a copy of this with the new ordering {@code newOrdering}
     *
     * @param newOrdering the new ordering
     * @return a copy of this with the new ordering
     */
    public final Poly setOrdering(Comparator<DegreeVector> newOrdering) {
        MonomialSet<Term> newData = new MonomialSet<>(newOrdering);
        newData.putAll(terms);
        return create(nVariables, newOrdering, newData);
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
        if (oth == this)
            return self;
        assertSameCoefficientRingWith(oth);
        return loadFrom(oth.terms);
    }

    final Poly loadFrom(MonomialSet<Term> map) {
        terms.clear();
        terms.putAll(map);
        release();
        return self;
    }

    /**
     * Makes a copy of this with the specified variable replaced with the unit
     *
     * @param variable  the variable
     * @param eliminate whether to reduces the resulting {@code nVariables} by one (eliminate variable from degree
     *                  vectors) or not
     */
    public final Poly dropVariable(int variable, boolean eliminate) {
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term term : terms)
            newData.add(term.without(variable));
        return create(eliminate ? nVariables - 1 : nVariables, newData);
    }

    /**
     * Makes a copy of this by inserting new variable (the indexes will be shifted)
     *
     * @param variable the variable
     */
    public final Poly insertVariable(int variable) {
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term term : terms)
            newData.add(term.insert(variable));
        return create(nVariables + 1, newData);
    }

    /** auxiliary method */
    final Poly setNVariables(int newNVariables) {
        if (newNVariables == nVariables)
            return self;
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term term : terms)
            newData.add(term.setDegreeVector(Arrays.copyOf(term.exponents, newNVariables)));
        return create(newNVariables, newData);
    }

    /**
     * Returns a copy of this with {@code nVariables = nVariables + 1}
     *
     * @return a copy of this with one additional (last) variable added
     * @see #insertVariable(int)
     */
    public final Poly joinNewVariable() {
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term term : terms)
            newData.add(term.joinNewVariable());
        return create(nVariables + 1, newData);
    }

    /** internal API */
    final Poly joinNewVariables(int newNVariables, int[] mapping) {
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term term : terms)
            newData.add(term.joinNewVariables(newNVariables, mapping));
        return create(newNVariables, newData);
    }

    /**
     * Returns the number of really used variables (those which are not units)
     *
     * @return the number of variables (those which are not units)
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
     * Returns the degree of this polynomial with respect to specified variable
     *
     * @param variable the variable
     * @return the degree of this polynomial with respect to specified variable
     */
    public final int degree(int variable) {
        int max = 0;
        for (Term db : terms)
            max = Math.max(max, db.exponents[variable]);
        return max;
    }

    /**
     * Returns an array of degrees of all variables, so that is i-th element of the result is the polynomial degree with
     * respect to i-th variable
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
     * Returns the array of exponents in which {@code variable} occurs in this polynomial
     *
     * @return the array of exponents in which {@code variable} occurs in this polynomial
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
     * Returns whether this poly is effectively univariate (not more than one variable is non-unit)
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
     * ring.
     *
     * @param variable variable which will be treated as univariate variable
     * @return univariate polynomial over the ring of multivariate coefficients
     * @throws IllegalArgumentException if this is not effectively a univariate polynomial
     */
    public final UnivariatePolynomial<Poly> asUnivariate(int variable) {
        MultivariateRing<Poly> ring = new MultivariateRing<>(self);
        Poly[] univarData = ring.createZeroesArray(degree(variable) + 1);
        for (Term e : terms)
            univarData[e.exponents[variable]].add(e.set(variable, 0));
        return UnivariatePolynomial.createUnsafe(ring, univarData);
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
     * the resulting polynomial have (nVariable - 1) multivariate variables (specified {@code variable} is eliminated)
     *
     * @param variable the variable
     * @return multivariate polynomial with coefficients being univariate polynomials over {@code variable}, the
     * resulting polynomial have (nVariable - 1) multivariate variables
     */
    public abstract MultivariatePolynomial<? extends IUnivariatePolynomial> asOverUnivariateEliminate(int variable);

    /**
     * Converts this to a multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     *
     * @param variables the variables to separate
     * @return multivariate polynomial with coefficients being multivariate polynomials polynomials over {@code
     * variables} that is polynomial in R[variables][other_variables]
     */
    public abstract MultivariatePolynomial<Poly> asOverMultivariate(int... variables);

    /**
     * Converts this to a multivariate polynomial with coefficients being multivariate polynomials polynomials over
     * {@code variables} that is polynomial in R[variables][other_variables]
     *
     * @param variables the variables to separate
     * @return multivariate polynomial with coefficients being multivariate polynomials polynomials over {@code
     * variables} that is polynomial in R[variables][other_variables]
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
        Poly result = uPoly.ring.getZero();
        for (int i = uPoly.degree(); i >= 0; --i) {
            if (uPoly.isZeroAt(i))
                continue;
            result.add(result.createMonomial(variable, i).multiply(uPoly.get(i)));
        }
        return result;
    }

    /**
     * Gives primitive part of this considered as R[variable][other_variables]
     *
     * @param variable the variable
     * @return primitive part of this considered as R[variable][other_variables]
     */
    public abstract Poly primitivePart(int variable);

    /**
     * Gives the content of this considered as R[variable][other_variables]
     *
     * @param variable the variable
     * @return the content of this considered as R[variable][other_variables]
     */
    public abstract IUnivariatePolynomial contentUnivariate(int variable);

    /**
     * Gives the content of this considered as R[variable][other_variables]
     *
     * @param variable the variable
     * @return the content of this considered as R[variable][other_variables]
     */
    @SuppressWarnings("unchecked")
    public final Poly content(int variable) {
        return asMultivariate(contentUnivariate(variable), nVariables, variable, ordering);
    }

    /**
     * Gives the content of this considered as R[other_variables][variable]
     *
     * @param variable the variable
     * @return the content of this considered as R[other_variables][variable]
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
     * Returns the leading coefficient of this viewed as R[other_variables][variable]
     *
     * @param variable the variable
     * @return multivariate leading coefficient of this viewed as R[other_variables][variable]
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
     * Set the leading coefficient of specified variable to a specified value (this is considered as
     * R[other_variables][variable])
     *
     * @param variable the variable
     * @param lc       the leading coefficient of this viewed as R[other_variables][variable]
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
            return createUnitTerm();
        return createTermWithUnitCoefficient(exponents);
    }

    static void setMin(int[] dv, int[] exponents) {
        for (int i = 0; i < exponents.length; ++i)
            if (dv[i] < exponents[i])
                exponents[i] = dv[i];
    }

    /** creates term with specified exponents and unit coefficient */
    public abstract Term createTermWithUnitCoefficient(int[] exponents);

    /** creates a unit term */
    public abstract Term createUnitTerm();

    /** creates a zero term */
    public abstract Term createZeroTerm();

    /**
     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of
     * the elements can't be exactly divided by the {@code monomial}. NOTE: if {@code null} is returned, the content of
     * {@code this} is destroyed.
     *
     * @param monomial monomial
     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
     */
    public final Poly divideDegreeVectorOrNull(DegreeVector monomial) {
        if (monomial.isZeroVector())
            return self;
        MonomialSet<Term> map = new MonomialSet<>(ordering);
        for (Term term : terms) {
            Term dv = term.divideOrNull(monomial);
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
     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of
     * the elements can't be exactly divided by the {@code monomial}. NOTE: if {@code null} is returned, the content of
     * {@code this} is destroyed.
     *
     * @param monomial monomial degrees
     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
     */
    public abstract Poly divideOrNull(Term monomial);

    /** add term to polynomial represented as terms */
    abstract void add(MonomialSet<Term> terms, Term term);

    /** subtract term from polynomial represented as terms */
    abstract void subtract(MonomialSet<Term> terms, Term term);

    @Override
    public final Poly add(Poly oth) {
        if (terms == oth.terms)
            return multiply(2);
        assertSameCoefficientRingWith(oth);
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
        assertSameCoefficientRingWith(oth);
        if (oth.isZero())
            return self;
        for (Term term : oth.terms)
            subtract(terms, term);
        release();
        return self;
    }

    /**
     * Adds {@code monomial} to this polynomial
     *
     * @param monomial the monomial
     * @return {@code this + monomial}
     */
    public final Poly add(Term monomial) {
        checkSameDomainWith(monomial);
        add(terms, monomial);
        release();
        return self;
    }

    /**
     * Subtracts {@code monomial} from this polynomial
     *
     * @param monomial the monomial
     * @return {@code this - monomial}
     */
    public final Poly subtract(Term monomial) {
        checkSameDomainWith(monomial);
        subtract(terms, monomial);
        release();
        return self;
    }

    /**
     * Adds monomials to this polynomial
     *
     * @param monomials terms
     * @return {@code this + monomials}
     */
    public final Poly add(Iterable<Term> monomials) {
        for (Term term : monomials)
            add(term);
        return self;
    }

    /**
     * Adds monomials to this polynomial
     *
     * @param monomials terms
     * @return {@code this + monomials}
     */
    public final Poly add(Term... monomials) {
        return add(Arrays.asList(monomials));
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
     * Multiplies {@code this} by the {@code monomial}
     *
     * @param monomial the monomial
     * @return {@code} this multiplied by the {@code monomial}
     */
    public abstract Poly multiply(Term monomial);

    /**
     * Multiplies {@code this} by the degree vector
     *
     * @param dv the degree vector
     * @return {@code} this multiplied by the degree vector
     */
    public final Poly multiplyByDegreeVector(DegreeVector dv) {
        return multiply(createUnitTerm().setDegreeVector(dv));
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

    /**
     * Set all coefficients to units
     */
    public final Poly setAllCoefficientsToUnit() {
        Term unit = createUnitTerm();
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
        return terms.keySet().stream().map(dv -> dv.select(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
    }

    /**
     * Returns skeleton of this poly with respect to all except specified {@code variables}
     *
     * @param variables the variables to exclude
     * @return skeleton of this poly with respect to all except specified {@code variables}
     */
    public final Set<DegreeVector> getSkeletonExcept(int... variables) {
        return terms.keySet().stream().map(dv -> dv.setZero(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton
     *
     * @param oth other multivariate polynomial
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton and {@code false} otherwise
     */
    public final boolean sameSkeletonQ(Poly oth) {
        return getSkeleton().equals(oth.getSkeleton());
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to test
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to specified {@code
     * variables} and {@code false} otherwise
     */
    public final boolean sameSkeletonQ(Poly oth, int... variables) {
        return getSkeleton(variables).equals(oth.getSkeleton(variables));
    }

    /**
     * Tests whether {@code this} and {@code oth} have the same skeleton with respect all except specified {@code
     * variables}
     *
     * @param oth       other multivariate polynomial
     * @param variables variables to exclude
     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to all except specified
     * {@code variables} and {@code false} otherwise
     */
    public final boolean sameSkeletonExceptQ(Poly oth, int... variables) {
        return getSkeletonExcept(variables).equals(oth.getSkeletonExcept(variables));
    }

    /**
     * Gives partial derivative with respect to specified variable (new instance created)
     *
     * @param variable the variable
     * @return partial derivative with respect to specified variable
     */
    public final Poly derivative(int variable) {return derivative(variable, 1);}

    /**
     * Gives partial derivative of specified {@code order} with respect to specified variable (new instance created)
     *
     * @param variable the variable
     * @param order    derivative order
     * @return partial derivative of specified {@code order} with respect to specified variable
     */
    public abstract Poly derivative(int variable, int order);

    /**
     * Gives (unevaluated) coefficient of Taylor series expansion for specified variable that is {@code derivative(poly,
     * variable, order) / order! }, where the derivative is formal derivative and calculated with arithmetic performed
     * in Z ring (to overcome possible zeros in Zp).
     *
     * @param variable the variable
     * @param order    derivative order
     * @return {@code derivative(poly, variable, order) / order! }, where the derivative is formal derivative and
     * calculated with arithmetic performed in Z ring (to overcome possible zeros in Zp)
     */
    public abstract Poly seriesCoefficient(int variable, int order);

    /**
     * Substitutes {@code 0} for {@code variable} (new instance created).
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code 0} substituted for {@code variable}
     */
    public final Poly evaluateAtZero(int variable) {
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        for (Term el : terms)
            if (el.exponents[variable] == 0)
                newData.add(el);
        return create(newData);
    }

    /**
     * Substitutes {@code 0} for all specified {@code variables} (new instance created).
     *
     * @param variables the variables
     * @return a new multivariate polynomial with {@code 0} substituted for all specified {@code variables}
     */
    public final Poly evaluateAtZero(int[] variables) {
        if (variables.length == 0)
            return clone();
        MonomialSet<Term> newData = new MonomialSet<>(ordering);
        out:
        for (Term el : terms) {
            for (int variable : variables)
                if (el.exponents[variable] != 0)
                    continue out;
            newData.add(el);
        }
        return create(newData);
    }

    /**
     * Gives the derivative vector
     *
     * @return derivative vector
     */
    public final Poly[] derivative() {
        Poly[] result = createArray(nVariables);
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
    public abstract Poly evaluateAtRandom(int variable, RandomGenerator rnd);

    /**
     * Evaluates {@code poly} at random point chosen in such way that the skeleton of evaluated version is the same as
     * of the original {@code poly} with respect to all except {@code variable} variables
     */
    public abstract Poly evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd);

    /**
     * Collector which collects stream of element to a UnivariatePolynomial
     */
    public static final class PolynomialCollector
            <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
            implements Collector<Term, Poly, Poly> {
        final Supplier<Poly> supplier;
        final BiConsumer<Poly, Term> accumulator = Poly::add;
        final BinaryOperator<Poly> combiner = (l, r) -> {
            l.add(r);
            return l;
        };
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

    @Override
    public final String toString() {
        return toString(WithVariables.defaultVars(nVariables));
    }
}
