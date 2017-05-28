package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.IntegersModulo;
import cc.r2.core.poly.UnivariatePolynomials;
import cc.r2.core.poly.univar.UnivariatePolynomial;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.regex.Pattern;

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
     * @param string    string polynomials
     * @param domain    domain
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

    /**
     * Converts this to a multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     *
     * @param variable variable
     * @return multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     */
    public MultivariatePolynomial<UnivariatePolynomial<E>> asOverUnivariate(int variable) {
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
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[] arrayNewInstance(int length) {return new MultivariatePolynomial[length];}

    @Override
    public boolean sameDomainWith(MultivariatePolynomial<E> oth) {
        return nVariables == oth.nVariables && domain.equals(oth.domain);
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
        for (MonomialTerm<E> e : terms) {
            MonomialTerm<E> term = e.setDomain(newDomain);
            if (!newDomain.isZero(term.coefficient))
                newData.add(term);
        }
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
        return size() == 1 && domain.isOne(terms.first().coefficient);
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
            E quot = domain.divideOrNull(term.coefficient, monomial.coefficient);
            if (quot == null)
                return null;
            MonomialTerm<E> dv = term.divide(monomial, quot);
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
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, domain);
        for (MonomialTerm<E> el : terms) {
            E val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable, val));
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, newData);
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
        return evaluate(new PrecomputedPowersHolder<>(values, domain), variables, ones(nVariables));
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
    MultivariatePolynomial<E> evaluate(PrecomputedPowersHolder<E> powers, int[] variables, int[] raiseFactors) {
        MonomialsSet<MonomialTerm<E>> newData = new MonomialsSet<>(ordering);
        for (MonomialTerm<E> el : terms) {
            MonomialTerm<E> r = el;
            E value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = domain.multiply(value, powers.pow(i, raiseFactors[i] * el.exponents[variables[i]]));
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

    private static final int SIZE_OF_POWERS_CACHE = 332;

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
                powers[i] = new PrecomputedPowers<>(cacheSize, point, domain);
        }

        E pow(int i, int exponent) {
            return powers[i].pow(exponent);
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
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> clone() {
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms.clone());
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

    static String[] defaultVars(int nVars) {
        char v = 'a';
        String[] vars = new String[nVars];
        for (int i = 0; i < nVars; i++)
            vars[i] = Character.toString(v++);
        return vars;
    }

}
