package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Function;
import java.util.function.ToLongFunction;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import static cc.r2.core.poly.Integers.Integers;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariatePolynomial<E> extends AMultivariatePolynomial<MonomialTerm<E>, MultivariatePolynomial<E>> {
    /** the domain */
    final Domain<E> domain;

    MultivariatePolynomial(int nVariables, Domain<E> domain, Comparator<DegreeVector> ordering, MonomialsSet<MonomialTerm<E>> terms) {
        super(nVariables, ordering, terms);
        this.domain = domain;
    }

    /* ============================================ Factory methods ============================================ */

    static <E> void add(MonomialsSet<MonomialTerm<E>> polynomial, MonomialTerm<E> term, Domain<E> domain) {
        if (domain.isZero(term.coefficient))
            return;
        MonomialTerm<E> pTerm = polynomial.get(term);
        if (pTerm == null)
            polynomial.add(term);
        else {
            E r = domain.add(pTerm.coefficient, term.coefficient);
            if (domain.isZero(r))
                polynomial.remove(pTerm);
            else
                polynomial.add(pTerm.setCoefficient(r));
        }
    }

    static <E> void subtract(MonomialsSet<MonomialTerm<E>> polynomial, MonomialTerm<E> term, Domain<E> domain) {
        if (domain.isZero(term.coefficient))
            return;
        MonomialTerm<E> pTerm = polynomial.get(term);
        if (pTerm == null)
            polynomial.add(term.negate(domain));
        else {
            E r = domain.subtract(pTerm.coefficient, term.coefficient);
            if (domain.isZero(r))
                polynomial.remove(pTerm);
            else
                polynomial.add(pTerm.setCoefficient(r));
        }
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param domain   the domain
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> create(int nVariables, Domain<E> domain, Comparator<DegreeVector> ordering, MonomialTerm<E>... terms) {
        if (terms.length == 0)
            throw new IllegalArgumentException("empty");
        MonomialsSet<MonomialTerm<E>> map = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : terms)
            add(map, term.setDomain(domain), domain);

        return new MultivariatePolynomial<>(nVariables, domain, ordering, map);
    }

    /**
     * Creates zero.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return zero
     */
    public static <E> MultivariatePolynomial<E> zero(int nVariables, Domain<E> domain, Comparator<DegreeVector> ordering) {
        return new MultivariatePolynomial<>(nVariables, domain, ordering, new MonomialsSet<>(ordering));
    }

    /**
     * Creates unit.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return unit
     */
    public static <E> MultivariatePolynomial<E> one(int nVariables, Domain<E> domain, Comparator<DegreeVector> ordering) {
        return create(nVariables, domain, ordering, MonomialTerm.withZeroExponents(nVariables, domain.getOne()));
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
        return Parser.parse(string, DegreeVector.LEX, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    string polynomials
     * @param domain    the domain
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Domain<E> domain, String... variables) {
        return Parser.parse(string, domain, DegreeVector.LEX, variables);
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
     * Parse multivariate polynomial from string.
     *
     * @param string string polynomials
     * @param domain domain
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Domain<E> domain) {
        return Parser.parse(string, domain);
    }

    /**
     * Converts multivariate polynomial over BigIntegers to multivariate polynomial over machine sized modular integers
     *
     * @param poly the polynomial
     * @return multivariate polynomial over machine sized modular integers
     * @throws ArithmeticException if some of coefficients will not exactly fit in a {@code long}.
     */
    public static lMultivariatePolynomialZp asLongPolyZp(MultivariatePolynomial<BigInteger> poly) {
        if (!(poly.domain instanceof IntegersModulo))
            throw new IllegalArgumentException("Poly is not over modular domain: " + poly.domain);
        IntegersModulo domain = (IntegersModulo) poly.domain;
        MonomialsSet<lMonomialTerm> terms = new MonomialsSet<>(poly.ordering);
        for (MonomialTerm<BigInteger> term : poly.terms)
            terms.add(new lMonomialTerm(term.exponents, term.totalDegree, term.coefficient.longValueExact()));
        return lMultivariatePolynomialZp.create(poly.nVariables, domain.asLong(), poly.ordering, terms);
    }

    /**
     * Converts univariate polynomial to multivariate.
     *
     * @param poly       univariate polynomial
     * @param nVariables number of variables in the result
     * @param variable   variable that will be used as a primary variable
     * @param ordering   ordering
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> asMultivariate(UnivariatePolynomial<E> poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        MonomialsSet<MonomialTerm<E>> map = new MonomialsSet<>(ordering);
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;
            map.add(new MonomialTerm<>(degreeVector, i, poly.get(i)));
        }
        return new MultivariatePolynomial<>(nVariables, poly.domain, ordering, map);
    }

    /**
     * Converts this polynomial to a univariate polynomial.
     *
     * @return univariate polynomial
     * @throws IllegalArgumentException if this is not effectively a univariate polynomial
     */
    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E> asUnivariate() {
        if (isConstant())
            return UnivariatePolynomial.create(domain, lc());
        int[] degrees = degrees();
        int theVar = -1;
        for (int i = 0; i < degrees.length; i++) {
            if (degrees[i] != 0) {
                if (theVar != -1)
                    throw new IllegalArgumentException("not a univariate polynomial: " + this);
                theVar = i;
            }
        }
        if (theVar == -1)
            throw new IllegalStateException("Not a univariate polynomial: " + this);
        E[] univarData = domain.createZeroesArray(degrees[theVar] + 1);
        for (MonomialTerm<E> e : terms)
            univarData[e.exponents[theVar]] = e.coefficient;
        return UnivariatePolynomial.createUnsafe(domain, univarData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomial<E>> asOverUnivariate(int variable) {
        UnivariatePolynomial<E> factory = UnivariatePolynomial.zero(domain);
        UnivariatePolynomials<UnivariatePolynomial<E>> pDomain = new UnivariatePolynomials<>(factory);
        MonomialsSet<MonomialTerm<UnivariatePolynomial<E>>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> e : terms) {
            add(newData, new MonomialTerm<>(
                            e.set(variable, 0).exponents,
                            factory.createMonomial(e.coefficient, e.exponents[variable])),
                    pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomial<E>> asOverUnivariateEliminate(int variable) {
        UnivariatePolynomial<E> factory = UnivariatePolynomial.zero(domain);
        UnivariatePolynomials<UnivariatePolynomial<E>> pDomain = new UnivariatePolynomials<>(factory);
        MonomialsSet<MonomialTerm<UnivariatePolynomial<E>>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> e : terms) {
            add(newData, new MonomialTerm<>(
                            e.without(variable).exponents,
                            factory.createMonomial(e.coefficient, e.exponents[variable])),
                    pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomial<E>> asOverMultivariate(int... variables) {
        Domain<MultivariatePolynomial<E>> domain = new MultivariatePolynomials<>(this);
        MonomialsSet<MonomialTerm<MultivariatePolynomial<E>>> terms = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : this) {
            int[] coeffExponents = new int[nVariables];
            for (int var : variables)
                coeffExponents[var] = term.exponents[var];
            MonomialTerm<E> restTerm = term.setZero(variables, this.domain.getOne());
            MonomialTerm<MultivariatePolynomial<E>> newTerm
                    = new MonomialTerm<>(
                    restTerm.exponents,
                    restTerm.totalDegree,
                    create(new MonomialTerm<>(coeffExponents, term.coefficient)));

            add(terms, newTerm, domain);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomial<E>> asOverMultivariateEliminate(int... variables) {
        variables = variables.clone();
        Arrays.sort(variables);
        int[] restVariables = ArraysUtil.intSetDifference(ArraysUtil.sequence(nVariables), variables);
        Domain<MultivariatePolynomial<E>> domain = new MultivariatePolynomials<>(create(variables.length, new MonomialsSet<>(ordering)));
        MonomialsSet<MonomialTerm<MultivariatePolynomial<E>>> terms = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : this) {
            int i = 0;
            int[] coeffExponents = new int[variables.length];
            for (int var : variables)
                coeffExponents[i++] = term.exponents[var];

            i = 0;
            int[] termExponents = new int[restVariables.length];
            for (int var : restVariables)
                termExponents[i++] = term.exponents[var];

            MonomialTerm<MultivariatePolynomial<E>> newTerm
                    = new MonomialTerm<>(
                    termExponents,
                    create(variables.length, this.domain, this.ordering, new MonomialTerm<>(coeffExponents, term.coefficient)));

            add(terms, newTerm, domain);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms);
    }

    /**
     * Converts multivariate polynomial over univariate polynomial domain to a multivariate polynomial over
     * coefficient domain
     *
     * @param poly     the polynomial
     * @param variable the variable to insert
     * @return multivariate polynomial over normal coefficient domain
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(MultivariatePolynomial<UnivariatePolynomial<E>> poly, int variable) {
        Domain<E> domain = poly.domain.getZero().domain;
        int nVariables = poly.nVariables + 1;
        MultivariatePolynomial<E> result = zero(nVariables, domain, poly.ordering);
        for (MonomialTerm<UnivariatePolynomial<E>> entry : poly.terms) {
            UnivariatePolynomial<E> uPoly = entry.coefficient;
            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new MonomialTerm<>(cdv, uPoly.get(i)));
            }
        }
        return result;
    }

    /**
     * Converts multivariate polynomial over multivariate polynomial domain to a multivariate polynomial over
     * coefficient domain
     *
     * @param poly the polynomial
     * @return multivariate polynomial over normal coefficient domain
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(MultivariatePolynomial<MultivariatePolynomial<E>> poly) {
        Domain<E> domain = poly.domain.getZero().domain;
        int nVariables = poly.nVariables;
        MultivariatePolynomial<E> result = zero(nVariables, domain, poly.ordering);
        for (MonomialTerm<MultivariatePolynomial<E>> term : poly.terms) {
            MultivariatePolynomial<E> uPoly = term.coefficient;
            result.add(uPoly.clone().multiply(new MonomialTerm<>(term.exponents, term.totalDegree, domain.getOne())));
        }
        return result;
    }

    /**
     * Converts multivariate polynomial over multivariate polynomial domain to a multivariate polynomial over
     * coefficient domain
     *
     * @param poly the polynomial
     * @return multivariate polynomial over normal coefficient domain
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(
            MultivariatePolynomial<MultivariatePolynomial<E>> poly,
            int[] coefficientVariables, int[] mainVariables) {
        Domain<E> domain = poly.domain.getZero().domain;
        int nVariables = coefficientVariables.length + mainVariables.length;
        MultivariatePolynomial<E> result = zero(nVariables, domain, poly.ordering);
        for (MonomialTerm<MultivariatePolynomial<E>> term : poly.terms) {
            MultivariatePolynomial<E> coefficient =
                    term.coefficient.joinNewVariables(nVariables, coefficientVariables);
            MonomialTerm<MultivariatePolynomial<E>> t = term.joinNewVariables(nVariables, mainVariables);
            result.add(coefficient.multiply(new MonomialTerm<>(t.exponents, t.totalDegree, domain.getOne())));
        }
        return result;
    }

    /**
     * Returns Z[X] polynomial formed from the coefficients of the poly.
     *
     * @param poly the polynomial
     * @param copy whether to copy the internal data
     * @return Z[X] version of the poly
     */
    public static MultivariatePolynomial<BigInteger> asPolyZ(MultivariatePolynomial<BigInteger> poly, boolean copy) {
        return new MultivariatePolynomial<>(poly.nVariables, Integers, poly.ordering, copy ? poly.terms.clone() : poly.terms);
    }

    /**
     * Converts Zp[x] polynomial to Z[x] polynomial formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @param poly Zp polynomial
     * @return Z[x] version of the poly with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     * @throws IllegalArgumentException is {@code poly.domain} is not a {@link IntegersModulo}
     */
    public static MultivariatePolynomial<BigInteger> asPolyZSymmetric(MultivariatePolynomial<BigInteger> poly) {
        if (!(poly.domain instanceof IntegersModulo))
            throw new IllegalArgumentException("Not a modular domain: " + poly.domain);
        IntegersModulo domain = (IntegersModulo) poly.domain;
        MonomialsSet<MonomialTerm<BigInteger>> newTerms = new MonomialsSet<>(poly.ordering);
        for (MonomialTerm<BigInteger> term : poly)
            newTerms.add(term.setCoefficient(domain.symmetricForm(term.coefficient)));

        return new MultivariatePolynomial<>(poly.nVariables, Integers, poly.ordering, newTerms);
    }


    /* ============================================ Main methods ============================================ */

    @Override
    public MultivariatePolynomial<E> contentAsPoly() {
        return createConstant(content());
    }

    @Override
    public MultivariatePolynomial<E> lcAsPoly() {
        return createConstant(lc());
    }

    @Override
    public MultivariatePolynomial<E> ccAsPoly() {
        return createConstant(cc());
    }

    @Override
    MultivariatePolynomial<E> create(int nVariables, MonomialsSet<MonomialTerm<E>> monomialTerms) {
        return new MultivariatePolynomial<>(nVariables, domain, ordering, monomialTerms);
    }

    @Override
    MonomialTerm<E> createTermWithUnitCoefficient(int[] exponents) {
        return new MonomialTerm<E>(exponents, domain.getOne());
    }

    @Override
    MonomialTerm<E> getUnitTerm() {
        return MonomialTerm.withZeroExponents(nVariables, domain.getOne());
    }

    @Override
    MonomialTerm<E> getZeroTerm() {
        return MonomialTerm.withZeroExponents(nVariables, domain.getZero());
    }

    @Override
    public boolean isOverField() {return domain.isField();}

    @Override
    public boolean isOverFiniteField() {return domain.isFiniteField();}

    /** {@inheritDoc} */
    @Override
    public boolean isOverZ() {return domain.equals(Integers);}

    @Override
    public BigInteger coefficientDomainCardinality() {return domain.cardinality();}

    @Override
    public BigInteger coefficientDomainCharacteristics() {return domain.characteristics();}

    @Override
    public boolean isOverPerfectPower() {
        return domain.isPerfectPower();
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerBase() {
        return domain.perfectPowerBase();
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerExponent() {
        return domain.perfectPowerExponent();
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[] arrayNewInstance(int length) {return new MultivariatePolynomial[length];}

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[][] arrayNewInstance2D(int length) {
        return new MultivariatePolynomial[length][];
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[][] arrayNewInstance2D(int length1, int length2) {
        return new MultivariatePolynomial[length1][length2];
    }

    @Override
    public boolean sameDomainWith(MultivariatePolynomial<E> oth) {
        return nVariables == oth.nVariables && domain.equals(oth.domain);
    }

    @Override
    boolean isZeroMonomial(MonomialTerm<E> a) {
        return domain.isZero(a.coefficient);
    }

    @Override
    MonomialTerm<E> multiply(MonomialTerm<E> a, MonomialTerm<E> b) {
        return a.multiply(b, domain.multiply(a.coefficient, b.coefficient));
    }

    @Override
    MonomialTerm<E> divideOrNull(MonomialTerm<E> a, MonomialTerm<E> b) {
        E e = domain.divideOrNull(a.coefficient, b.coefficient);
        if (e == null)
            return null;
        return a.divide(b, e);
    }

    /** release caches */
    @Override
    protected void release() {
        super.release();
        /* add cache in the future */
    }


    /**
     * Switches to another domain specified by {@code newDomain}
     *
     * @param newDomain the new domain
     * @return a copy of this reduced to the domain specified by {@code newDomain}
     */
    public MultivariatePolynomial<E> setDomain(Domain<E> newDomain) {
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> e : terms)
            add(newData, e.setDomain(newDomain));
        return new MultivariatePolynomial<>(nVariables, newDomain, ordering, newData);
    }

    public MultivariatePolynomial<E> setDomainUnsafe(Domain<E> newDomain) {
        return new MultivariatePolynomial<>(nVariables, newDomain, ordering, terms);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public MultivariatePolynomial<E> createConstant(E val) {
        MonomialsSet<MonomialTerm<E>> data = new MonomialsSet<>(ordering);
        if (!domain.isZero(val))
            data.add(MonomialTerm.withZeroExponents(nVariables, val));
        return new MultivariatePolynomial<>(nVariables, domain, ordering, data);
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
     * @param lc       the leading coefficient
     * @return linear polynomial {@code cc + lc * variable}
     */
    public MultivariatePolynomial<E> createLinear(int variable, E cc, E lc) {
        MonomialsSet<MonomialTerm<E>> data = new MonomialsSet<>(ordering);
        if (!domain.isZero(cc))
            data.add(MonomialTerm.withZeroExponents(nVariables, cc));
        if (!domain.isZero(lc)) {
            int[] lcDegreeVector = new int[nVariables];
            lcDegreeVector[variable] = 1;
            data.add(new MonomialTerm<>(lcDegreeVector, 1, lc));
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, data);
    }


    @Override
    public boolean isMonic() {
        return domain.isOne(lc());
    }

    @Override
    public int signum() {
        return domain.signum(lc());
    }

    @Override
    public boolean isZero() {
        return terms.size() == 0;
    }

    @Override
    public boolean isOne() {
        if (size() != 1)
            return false;
        MonomialTerm<E> lt = terms.first();
        return lt.isZeroVector() && domain.isOne(lt.coefficient);
    }

    @Override
    public boolean isUnitCC() {
        return domain.isOne(cc());
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && terms.first().isZeroVector());
    }


    /** {@inheritDoc} */
    @Override
    public MonomialTerm<E> lt() {
        return size() == 0 ? MonomialTerm.withZeroExponents(nVariables, domain.getZero()) : terms.last();
    }

    /**
     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the ordering.
     *
     * @return leading coefficient of this polynomial
     */
    public E lc() {
        return lt().coefficient;
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
        terms.add(lt().setCoefficient(domain.valueOf(val)));
        return this;
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public E cc() {
        MonomialTerm<E> zero = MonomialTerm.withZeroExponents(nVariables, domain.getZero());
        return terms.getOrDefault(zero, zero).coefficient;
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public E content() {
        return domain.gcd(coefficients());
    }

    /**
     * Returns iterable over polynomial coefficients
     *
     * @return iterable over polynomial coefficients
     */
    public Iterable<E> coefficients() {
        return () -> new It<>(terms.iterator());
    }

    /**
     * Returns polynomial coefficients
     *
     * @return polynomial coefficients
     */
    public E[] coefficientsArray() {
        if (isZero())
            return domain.createZeroesArray(1);

        E[] array = domain.createArray(size());
        int i = 0;
        for (MonomialTerm<E> term : this)
            array[i++] = term.coefficient;
        return array;
    }

    private static class It<V> implements Iterator<V> {
        final Iterator<MonomialTerm<V>> inner;

        public It(Iterator<MonomialTerm<V>> inner) {
            this.inner = inner;
        }

        @Override
        public boolean hasNext() {
            return inner.hasNext();
        }

        @Override
        public V next() {
            return inner.next().coefficient;
        }
    }

    @Override
    public MultivariatePolynomial<E> primitivePart(int variable) {
        return asNormalMultivariate(asOverUnivariateEliminate(variable).primitivePart(), variable);
    }

    @Override
    public UnivariatePolynomial<E> contentUnivariate(int variable) {
        return asOverUnivariate(variable).content();
    }

    @Override
    public MultivariatePolynomial<E> primitivePart() {
        MultivariatePolynomial<E> r = divideOrNull(content());
        assert r != null;
        return r;
    }

    @Override
    public MultivariatePolynomial<E> primitivePartSameSign() {
        E c = content();
        if (domain.signum(c) < 0)
            c = domain.negate(c);
        MultivariatePolynomial<E> r = divideOrNull(c);
        assert r != null;
        return r;
    }

    @Override
    public MultivariatePolynomial<E> divideByLC(MultivariatePolynomial<E> other) {
        return divideOrNull(other.lc());
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
        if (domain.isField())
            return multiply(domain.reciprocal(factor)); // <- this is typically faster than the division
        for (Entry<DegreeVector, MonomialTerm<E>> entry : terms.entrySet()) {
            MonomialTerm<E> term = entry.getValue();
            E quot = domain.divideOrNull(term.coefficient, factor);
            if (quot == null)
                return null;
            entry.setValue(term.setCoefficient(quot));
        }
        release();
        return this;
    }

    /** {@inheritDoc} */
    @Override
    public MultivariatePolynomial<E> divideOrNull(MonomialTerm<E> monomial) {
        if (monomial.isZeroVector())
            return divideOrNull(monomial.coefficient);
        MonomialsSet<MonomialTerm<E>> map = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : terms) {
            MonomialTerm<E> dv = divideOrNull(term, monomial);
            if (dv == null)
                return null;
            map.add(dv);
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
    @Override
    public MultivariatePolynomial<E> monic() {
        if (isZero())
            return this;
        return divideOrNull(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public MultivariatePolynomial<E> monic(E factor) {
        E lc = lc();
        return multiply(factor).divideOrNull(lc);
    }

    @Override
    public MultivariatePolynomial<E> monicWithLC(MultivariatePolynomial<E> other) {
        return monic(other.lc());
    }

    /**
     * Substitutes {@code 0} for {@code variable}.
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code 0} substituted for {@code variable}
     */
    public MultivariatePolynomial<E> evaluateAtZero(int variable) {
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> el : terms)
            if (el.exponents[variable] == 0)
                newData.add(el);
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newData);
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
        if (domain.isZero(value))
            return evaluateAtZero(variable);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, domain);
        return evaluate(variable, powers);
    }

    /**
     * Substitutes {@code value} for {@code variable}.
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    MultivariatePolynomial<E> evaluate(int variable, PrecomputedPowers<E> powers) {
        if (domain.isZero(powers.value))
            return evaluateAtZero(variable);
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> el : terms) {
            E val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable, val));
        }
        return new MultivariatePolynomial<E>(nVariables, domain, ordering, newData);
    }

    /**
     * Substitutes {@code 0} for all specified {@code variables}.
     *
     * @param variables the variables
     * @return a new multivariate polynomial with {@code 0} substituted for all specified {@code variables}
     */
    public MultivariatePolynomial<E> evaluateAtZero(int[] variables) {
        if (variables.length == 0)
            return this.clone();
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        out:
        for (MonomialTerm<E> el : terms) {
            for (int variable : variables)
                if (el.exponents[variable] != 0)
                    continue out;
            newData.add(el);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newData);
    }

    UnivariatePolynomial<E> evaluateAtZeroAllExcept(int variable) {
        E[] uData = domain.createArray(degree(variable) + 1);
        out:
        for (MonomialTerm<E> el : terms) {
            if (el.totalDegree != 0 && el.exponents[variable] == 0)
                continue;
            for (int i = 0; i < nVariables; ++i)
                if (i != variable && el.exponents[i] != 0)
                    continue out;
            int uExp = el.exponents[variable];
            uData[uExp] = domain.add(uData[uExp], el.coefficient);
        }
        return UnivariatePolynomial.createUnsafe(domain, uData);
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
        for (E value : values)
            if (!domain.isZero(value))
                return evaluate(new PrecomputedPowersHolder<>(nVariables, variables, values, domain), variables, ones(nVariables));

        // <- all values are zero
        return evaluateAtZero(variables);
    }

    /* cached array of units */
    private static int[][] ones = new int[16][];

    private static int[] ones(int len) {
        if (len < ones.length)
            return ones[len];
        else
            return ArraysUtil.arrayOf(1, len);
    }

    static {
        for (int i = 0; i < ones.length; i++)
            ones[i] = ArraysUtil.arrayOf(1, i);
    }

    /** substitutes {@code values} for {@code variables} */
    @SuppressWarnings("unchecked")
    MultivariatePolynomial<E> evaluate(PrecomputedPowersHolder<E> powers, int[] variables) {
        return evaluate(powers, variables, ones(variables.length));
    }

    /** substitutes {@code values} for {@code variables} */
    @SuppressWarnings("unchecked")
    MultivariatePolynomial<E> evaluate(PrecomputedPowersHolder<E> powers, int[] variables, int[] raiseFactors) {
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> el : terms) {
            MonomialTerm<E> r = el;
            E value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = domain.multiply(value, powers.pow(variables[i], raiseFactors[i] * el.exponents[variables[i]]));
                r = r.setZero(variables[i], value);
            }

            add(newData, r);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newData);
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
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, domain);
        for (MonomialTerm<E> el : terms) {
            E val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.without(variable, val));
        }
        return new MultivariatePolynomial<>(nVariables - 1, domain, ordering, newData);
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
        // TODO: store map here
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
            if (exponent >= precomputedPowers.length)
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
        final int cacheSize;
        final Domain<E> domain;
        final PrecomputedPowers<E>[] powers;

        PrecomputedPowersHolder(int nVariables, int[] variables, E[] points, Domain<E> domain) {
            this(SIZE_OF_POWERS_CACHE, nVariables, variables, points, domain);
        }

        @SuppressWarnings("unchecked")
        PrecomputedPowersHolder(int cacheSize, int nVariables, int[] variables, E[] points, Domain<E> domain) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = new PrecomputedPowers[nVariables];
            for (int i = 0; i < points.length; i++)
                powers[variables[i]] = points[i] == null ? null : new PrecomputedPowers<>(cacheSize, points[i], domain);
        }

        void set(int i, E point) {
            if (powers[i] == null || !powers[i].value.equals(point))
                powers[i] = new PrecomputedPowers<>(cacheSize, point, domain);
        }

        E pow(int variable, int exponent) {
            return powers[variable].pow(exponent);
        }
    }


    /**
     * Substitutes {@code variable -> poly} and returns the result (copied)
     *
     * @param variable the variable
     * @param poly     the replacement for the variable
     * @return a copy of this with  {@code variable -> poly}
     */
    public MultivariatePolynomial<E> substitute(int variable, MultivariatePolynomial<E> poly) {
        if (poly.isConstant())
            return evaluate(variable, poly.cc());
        PrecomputedSubstitution<E> subsPowers;
        if (poly.isEffectiveUnivariate())
            subsPowers = new USubstitution<>(poly.asUnivariate(), poly.univariateVariable(), nVariables, ordering);
        else
            subsPowers = new MSubstitution<>(poly);

        MultivariatePolynomial<E> result = createZero();
        for (MonomialTerm<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent == 0) {
                result.add(term);
                continue;
            }

            result.add(subsPowers.pow(exponent).multiply(term.setZero(variable)));
        }
        return result;
    }

    /**
     * Substitutes {@code variable -> variable + shift} and returns the result (copied)
     *
     * @param variable the variable
     * @param shift    shift amount
     * @return a copy of this with  {@code variable -> variable + shift}
     */
    public MultivariatePolynomial<E> shift(int variable, long shift) {
        return shift(variable, domain.valueOf(shift));
    }

    /**
     * Substitutes {@code variable -> variable + shift} and returns the result (copied)
     *
     * @param variable the variable
     * @param shift    shift amount
     * @return a copy of this with  {@code variable -> variable + shift}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> shift(int variable, E shift) {
        if (domain.isZero(shift))
            return clone();

        USubstitution<E> shifts = new USubstitution<>(UnivariatePolynomial.create(domain, shift, domain.getOne()), variable, nVariables, ordering);
        MultivariatePolynomial<E> result = createZero();
        for (MonomialTerm<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent == 0) {
                result.add(term);
                continue;
            }

            result.add(shifts.pow(exponent).multiply(term.setZero(variable)));
        }
        return result;
    }

    /**
     * Substitutes {@code variable -> variable + shift} for each variable from {@code variables} array
     *
     * @param variables the variables
     * @param shifts    the corresponding shifts
     * @return a copy of this with  {@code variable -> variable + shift}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> shift(int[] variables, E[] shifts) {
        PrecomputedSubstitution<E>[] powers = new PrecomputedSubstitution[nVariables];
        boolean allShiftsAreZero = true;
        for (int i = 0; i < variables.length; ++i) {
            if (!domain.isZero(shifts[i]))
                allShiftsAreZero = false;
            powers[variables[i]] = new USubstitution<>(UnivariatePolynomial.create(domain, shifts[i], domain.getOne()), variables[i], nVariables, ordering);
        }

        if (allShiftsAreZero)
            return clone();

        PrecomputedSubstitutions<E> calculatedShifts = new PrecomputedSubstitutions<>(powers);

        MultivariatePolynomial<E> result = createZero();
        for (MonomialTerm<E> term : this) {
            MultivariatePolynomial<E> temp = createOne();
            for (int variable : variables) {
                if (term.exponents[variable] != 0) {
                    temp = temp.multiply(calculatedShifts.getSubstitutionPower(variable, term.exponents[variable]));
                    term = term.setZero(variable);
                }
            }
            if (temp.isOne()) {
                result.add(term);
                continue;
            }
            result.add(temp.multiply(term));
        }
        return result;
    }

    static final class PrecomputedSubstitutions<E> {
        final PrecomputedSubstitution<E>[] subs;

        public PrecomputedSubstitutions(PrecomputedSubstitution<E>[] subs) {
            this.subs = subs;
        }

        MultivariatePolynomial<E> getSubstitutionPower(int var, int exponent) {
            if (subs[var] == null)
                throw new IllegalArgumentException();

            return subs[var].pow(exponent);
        }
    }

    interface PrecomputedSubstitution<E> {
        MultivariatePolynomial<E> pow(int exponent);
    }

    static final class USubstitution<E> implements PrecomputedSubstitution<E> {
        final int variable;
        final int nVariables;
        final Comparator<DegreeVector> ordering;
        final UnivariatePolynomial<E> base;
        final TIntObjectHashMap<UnivariatePolynomial<E>> uCache = new TIntObjectHashMap<>();
        final TIntObjectHashMap<MultivariatePolynomial<E>> mCache = new TIntObjectHashMap<>();

        USubstitution(UnivariatePolynomial<E> base, int variable, int nVariables, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.variable = variable;
            this.ordering = ordering;
            this.base = base;
        }

        @Override
        public MultivariatePolynomial<E> pow(int exponent) {
            MultivariatePolynomial<E> cached = mCache.get(exponent);
            if (cached != null)
                return cached.clone();
            UnivariatePolynomial<E> r = CommonPolynomialsArithmetics.polyPow(base, exponent, true, uCache);
            mCache.put(exponent, cached = asMultivariate(r, nVariables, variable, ordering));
            return cached.clone();
        }
    }

    static final class MSubstitution<E> implements PrecomputedSubstitution<E> {
        final MultivariatePolynomial<E> base;
        final TIntObjectHashMap<MultivariatePolynomial<E>> cache = new TIntObjectHashMap<>();

        MSubstitution(MultivariatePolynomial<E> base) {
            this.base = base;
        }

        @Override
        public MultivariatePolynomial<E> pow(int exponent) {
            return CommonPolynomialsArithmetics.polyPow(base, exponent, true, cache);
        }
    }

    @Override
    public MultivariatePolynomial<E> negate() {
        for (Entry<DegreeVector, MonomialTerm<E>> entry : terms.entrySet()) {
            MonomialTerm<E> term = entry.getValue();
            entry.setValue(term.negate(domain));
        }
        release();
        return this;
    }

    @Override
    void add(MonomialsSet<MonomialTerm<E>> terms, MonomialTerm<E> term) {
        add(terms, term, domain);
    }

    @Override
    void subtract(MonomialsSet<MonomialTerm<E>> terms, MonomialTerm<E> term) {
        subtract(terms, term, domain);
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
        add(terms, MonomialTerm.withZeroExponents(nVariables, oth));
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
        return add(domain.negate(domain.valueOf(oth)));
    }

    @Override
    public MultivariatePolynomial<E> increment() {
        return add(domain.getOne());
    }

    @Override
    public MultivariatePolynomial<E> decrement() {
        return subtract(domain.getOne());
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
        for (Entry<DegreeVector, MonomialTerm<E>> entry : terms.entrySet()) {
            MonomialTerm<E> term = entry.getValue();
            E val = domain.multiply(term.coefficient, factor);
            if (!domain.isZero(val))
                entry.setValue(term.setCoefficient(val));
        }
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> multiplyByLC(MultivariatePolynomial<E> other) {
        return multiply(other.lc());
    }

    /** {@inheritDoc} */
    @Override
    public MultivariatePolynomial<E> multiply(MonomialTerm<E> term) {
        checkSameDomainWith(term);
        if (term.isZeroVector())
            return multiply(term.coefficient);
        if (domain.isZero(term.coefficient))
            return toZero();

        MonomialsSet<MonomialTerm<E>> newMap = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> thisElement : terms) {
            E val = domain.multiply(thisElement.coefficient, term.coefficient);
            if (!domain.isZero(val))
                newMap.add(thisElement.multiply(term, val));
        }

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> multiply(long factor) {
        return multiply(domain.valueOf(factor));
    }

    @Override
    public MultivariatePolynomial<E> multiplyByBigInteger(BigInteger factor) {
        return multiply(domain.valueOfBigInteger(factor));
    }

    @Override
    public MultivariatePolynomial<E> multiply(MultivariatePolynomial<E> oth) {
        checkSameDomainWith(oth);
        MonomialsSet<MonomialTerm<E>> newMap = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> othElement : oth.terms)
            for (MonomialTerm<E> thisElement : terms)
                add(newMap, thisElement.multiply(othElement, domain.multiply(thisElement.coefficient, othElement.coefficient)));

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> square() {
        return multiply(this);
    }

    @Override
    MultivariatePolynomial<E> evaluateAtRandom(int variable, RandomGenerator rnd) {
        return evaluate(variable, domain.randomElement(rnd));
    }

    @Override
    MultivariatePolynomial<E> evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd) {
        //desired skeleton
        Set<DegreeVector> skeleton = getSkeletonExcept(variable);
        MultivariatePolynomial<E> tmp;
        do {
            E randomPoint = domain.randomElement(rnd);
            tmp = evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    @Override
    public MultivariatePolynomial<E> derivative(int variable, int order) {
        MonomialsSet<MonomialTerm<E>> newTerms = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;
            E newCoefficient = term.coefficient;
            for (int i = 0; i < order; ++i)
                newCoefficient = domain.multiply(newCoefficient, domain.valueOf(exponent - i));
            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            add(newTerms, new MonomialTerm<>(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newTerms);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static BigInteger seriesCoefficientFactor0(int exponent, int order, IntegersModulo domain) {
        if (!domain.modulus.isInt() || order < domain.modulus.intValueExact())
            return seriesCoefficientFactor1(exponent, order, domain);
        return BigInteger.valueOf(lMultivariatePolynomialZp.seriesCoefficientFactor(exponent, order, domain.asMachineSizedDomain()));
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static <E> E seriesCoefficientFactor1(int exponent, int order, Domain<E> domain) {
        E factor = domain.getOne();
        for (int i = 0; i < order; ++i)
            factor = domain.multiply(factor, domain.valueOf(exponent - i));
        factor = domain.divideExact(factor, domain.factorial(order));
        return factor;
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static <E> E seriesCoefficientFactor2(int exponent, int order, Domain<E> domain) {
        BigInteger factor = BigInteger.ONE;
        for (int i = 0; i < order; ++i)
            factor = factor.multiply(BigInteger.valueOf(exponent - i));
        factor = factor.divideExact(Integers.factorial(order));
        return domain.valueOfBigInteger(factor);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    @SuppressWarnings("unchecked")
    static <E> E seriesCoefficientFactor(int exponent, int order, Domain<E> domain) {
        if (domain instanceof IntegersModulo)
            return (E) seriesCoefficientFactor0(exponent, order, (IntegersModulo) domain);
        BigInteger characteristics = domain.characteristics();
        if (characteristics == null || !characteristics.isInt() || characteristics.intValueExact() > order)
            return seriesCoefficientFactor1(exponent, order, domain);

        return seriesCoefficientFactor2(exponent, order, domain);
    }

    @Override
    public MultivariatePolynomial<E> seriesCoefficient(int variable, int order) {
        if (order == 0)
            return clone();
        if (isConstant())
            return createZero();

        MonomialsSet<MonomialTerm<E>> newTerms = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;

            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            E newCoefficient = domain.multiply(term.coefficient, seriesCoefficientFactor(exponent, order, domain));
            add(newTerms, new MonomialTerm<>(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newTerms);
    }

    /**
     * Returns a stream of coefficients of this
     *
     * @return stream of coefficients
     */
    public Stream<E> stream() {
        return terms.values().stream().map(e -> e.coefficient);
    }

    /**
     * Maps terms of this using specified mapping function
     *
     * @param newDomain the new domain
     * @param mapper    mapping
     * @param <T>       new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms
     */
    public <T> MultivariatePolynomial<T> mapTerms(Domain<T> newDomain, Function<MonomialTerm<E>, MonomialTerm<T>> mapper) {
        return terms.values()
                .stream()
                .map(mapper)
                .collect(new PolynomialCollector<>(() -> zero(nVariables, newDomain, ordering)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newDomain the new domain
     * @param mapper    mapping
     * @param <T>       new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public <T> MultivariatePolynomial<T> mapCoefficients(Domain<T> newDomain, Function<E, T> mapper) {
        return mapTerms(newDomain, t -> new MonomialTerm<>(t.exponents, t.totalDegree, mapper.apply(t.coefficient)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newDomain the new domain
     * @param mapper    mapping
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public lMultivariatePolynomialZp mapCoefficients(lIntegersModulo newDomain, ToLongFunction<E> mapper) {
        return terms.values()
                .stream()
                .map(t -> new lMonomialTerm(t.exponents, t.totalDegree, mapper.applyAsLong(t.coefficient)))
                .collect(new PolynomialCollector<>(() -> lMultivariatePolynomialZp.zero(nVariables, newDomain, ordering)));
    }

    @Override
    public int compareTo(MultivariatePolynomial<E> oth) {
        int c = Integer.compare(size(), oth.size());
        if (c != 0)
            return c;
        Iterator<MonomialTerm<E>>
                thisIt = iterator(),
                othIt = oth.iterator();

        while (thisIt.hasNext() && othIt.hasNext()) {
            MonomialTerm<E>
                    a = thisIt.next(),
                    b = othIt.next();

            if ((c = ordering.compare(a, b)) != 0)
                return c;

            if ((c = domain.compare(a.coefficient, b.coefficient)) != 0)
                return c;
        }
        return 0;
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> clone() {
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms.clone());
    }

    @Override
    public MultivariatePolynomial<E> parsePoly(String string) {
        return parse(string, domain, ordering);
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
        for (MonomialTerm<E> term : terms) {
            E coeff = term.coefficient;
            if (domain.isZero(coeff))
                continue;
            String monomialString = term.toString(vars);
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

    @Override
    public String coefficientDomainToString() {
        return domain.toString();
    }

    static String[] defaultVars(int nVars) {
        char v = 'a';
        String[] vars = new String[nVars];
        for (int i = 0; i < nVars; i++)
            vars[i] = Character.toString(v++);
        return vars;
    }

}
