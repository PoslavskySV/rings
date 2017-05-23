package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Integers;
import cc.r2.core.poly.LongArithmetics;
import cc.r2.core.poly.UnivariatePolynomials;
import cc.r2.core.poly.lIntegersModulo;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Map.Entry;
import java.util.Set;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class lMultivariatePolynomial extends AMultivariatePolynomial<lMonomialTerm, lMultivariatePolynomial> {
    final lIntegersModulo domain;

    private lMultivariatePolynomial(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, MonomialsSet<lMonomialTerm> lMonomialTerms) {
        super(nVariables, ordering, lMonomialTerms);
        this.domain = domain;
    }

    private static void add(MonomialsSet<lMonomialTerm> polynomial, lMonomialTerm term, lIntegersModulo domain) {
        if (term.coefficient == 0)
            return;
        lMonomialTerm pTerm = polynomial.get(term);
        if (pTerm == null)
            polynomial.add(term);
        else {
            long r = domain.add(pTerm.coefficient, term.coefficient);
            if (r == 0)
                polynomial.remove(pTerm);
            else
                polynomial.add(pTerm.setCoefficient(r));
        }
    }

    private static void subtract(MonomialsSet<lMonomialTerm> polynomial, lMonomialTerm term, lIntegersModulo domain) {
        if (term.coefficient == 0)
            return;
        lMonomialTerm pTerm = polynomial.get(term);
        if (pTerm == null)
            polynomial.add(term.negate(domain));
        else {
            long r = domain.subtract(pTerm.coefficient, term.coefficient);
            if (r == 0)
                polynomial.remove(pTerm);
            else
                polynomial.add(pTerm.setCoefficient(r));
        }
    }

    /**
     * Creates multivariate polynomial from a set of monomial terms
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   term ordering
     * @param terms      the monomial terms
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomial create(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, MonomialsSet<lMonomialTerm> terms) {
        return new lMultivariatePolynomial(nVariables, domain, ordering, terms);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param domain   the domain
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomial create(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, lMonomialTerm... terms) {
        if (terms.length == 0)
            throw new IllegalArgumentException("empty");
        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
        for (lMonomialTerm term : terms)
            add(map, term.setDomain(domain), domain);

        return new lMultivariatePolynomial(nVariables, domain, ordering, map);
    }

    /**
     * Creates zero.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return zero
     */
    public static lMultivariatePolynomial zero(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
        return new lMultivariatePolynomial(nVariables, domain, ordering, new MonomialsSet<>(ordering));
    }

    /**
     * Creates unit.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return unit
     */
    public static lMultivariatePolynomial one(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
        return create(nVariables, domain, ordering, lMonomialTerm.withZeroExponents(nVariables, 1L));
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    string polynomials
     * @param domain    domain
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomial parse(String string, lIntegersModulo domain, String... variables) {
        return MultivariatePolynomial.asLongPolyZp(Parser.parse(string, domain.asDomain(), DegreeVector.LEX, variables));
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
    public static lMultivariatePolynomial parse(String string, lIntegersModulo domain, Comparator<DegreeVector> ordering, String... variables) {
        return MultivariatePolynomial.asLongPolyZp(Parser.parse(string, domain.asDomain(), ordering, variables));
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
    public static lMultivariatePolynomial asMultivariate(lUnivariatePolynomialZp poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;
            map.add(new lMonomialTerm(degreeVector, i, poly.get(i)));
        }
        return new lMultivariatePolynomial(nVariables, poly.domain, ordering, map);
    }

    /**
     * Converts this polynomial to a univariate polynomial.
     *
     * @return univariate polynomial
     * @throws IllegalArgumentException if this is not effectively a univariate polynomial
     */
    @Override
    public lUnivariatePolynomialZp asUnivariate() {
        if (isConstant())
            return lUnivariatePolynomialZp.createUnsafe(domain, new long[]{lc()});
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
        long[] univarData = new long[degrees[theVar] + 1];
        for (lMonomialTerm e : terms)
            univarData[e.exponents[theVar]] = e.coefficient;
        return lUnivariatePolynomialZp.createUnsafe(domain, univarData);
    }

    /**
     * Converts this to a multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     *
     * @param variable variable
     * @return multivariate polynomial with coefficients being univariate polynomials over {@code variable}
     */
    public MultivariatePolynomial<lUnivariatePolynomialZp> asOverUnivariate(int variable) {
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(domain);
        UnivariatePolynomials<lUnivariatePolynomialZp> pDomain = new UnivariatePolynomials<>(factory);
        MonomialsSet<MonomialTerm<lUnivariatePolynomialZp>> newData = new MonomialsSet<>(ordering);
        for (lMonomialTerm e : terms) {
            MonomialTerm<lUnivariatePolynomialZp> eMonomialTerm = new MonomialTerm<>(
                    e.without(variable).exponents,
                    factory.createMonomial(e.coefficient, e.exponents[variable]));
            MultivariatePolynomial.add(newData, eMonomialTerm, pDomain);
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
    public static lMultivariatePolynomial asNormalMultivariate(MultivariatePolynomial<lUnivariatePolynomialZp> poly, int variable) {
        lIntegersModulo domain = poly.domain.getZero().domain;
        int nVariables = poly.nVariables + 1;
        lMultivariatePolynomial result = zero(nVariables, domain, poly.ordering);
        for (MonomialTerm<lUnivariatePolynomialZp> entry : poly.terms) {
            lUnivariatePolynomialZp uPoly = entry.coefficient;
            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new lMonomialTerm(cdv, uPoly.get(i)));
            }
        }
        return result;
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     * represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[X] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     */
    public MultivariatePolynomial<BigInteger> asPolyZSymmetric() {
        MonomialsSet<MonomialTerm<BigInteger>> bTerms = new MonomialsSet<>(ordering);
        for (lMonomialTerm t : this)
            bTerms.add(new MonomialTerm<>(t.exponents, t.totalDegree, BigInteger.valueOf(domain.symmetricForm(t.coefficient))));
        return new MultivariatePolynomial<>(nVariables, Integers.Integers, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> asPolyZ() {
        MonomialsSet<MonomialTerm<BigInteger>> bTerms = new MonomialsSet<>(ordering);
        for (lMonomialTerm t : this)
            bTerms.add(t.asBigInteger());
        return new MultivariatePolynomial<>(nVariables, Integers.Integers, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> toBigPoly() {
        MonomialsSet<MonomialTerm<BigInteger>> bTerms = new MonomialsSet<>(ordering);
        for (lMonomialTerm t : this)
            bTerms.add(t.asBigInteger());
        return new MultivariatePolynomial<>(nVariables, domain.asDomain(), ordering, bTerms);
    }

    /* ============================================ Main methods ============================================ */

    @Override
    lMultivariatePolynomial create(int nVariables, MonomialsSet<lMonomialTerm> lMonomialTerms) {
        return new lMultivariatePolynomial(nVariables, domain, ordering, lMonomialTerms);
    }

    @Override
    lMonomialTerm createTermWithUnitCoefficient(int[] exponents) {
        return new lMonomialTerm(exponents, 1L);
    }

    @Override
    lMonomialTerm getUnitTerm() {
        return lMonomialTerm.withZeroExponents(nVariables, 1L);
    }

    @Override
    lMonomialTerm getZeroTerm() {
        return lMonomialTerm.withZeroExponents(nVariables, 0L);
    }

    @Override
    public boolean isOverField() {return true;}

    @Override
    public boolean isOverFiniteField() {return true;}

    @Override
    public boolean isOverZ() {return false;}

    @Override
    public BigInteger coefficientDomainCardinality() {return BigInteger.valueOf(domain.modulus);}

    @Override
    public BigInteger coefficientDomainCharacteristics() {return BigInteger.valueOf(domain.modulus);}

    @Override
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial[] arrayNewInstance(int length) {return new lMultivariatePolynomial[length];}

    @Override
    public boolean sameDomainWith(lMultivariatePolynomial oth) {
        return nVariables == oth.nVariables && domain.equals(oth.domain);
    }

    /** release caches */
    @Override
    protected void release() {
        super.release();
        /* add cache in the future */
    }

    /**
     * Switches to another domain specified by {@code newModulus}
     *
     * @param newModulus the new modulus
     * @return a copy of this reduced to the domain specified by {@code newModulus}
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial setDomain(long newModulus) {
        return setDomain(new lIntegersModulo(newModulus));
    }

    /**
     * Switches to another domain specified by {@code newDomain}
     *
     * @param newDomain the new domain
     * @return a copy of this reduced to the domain specified by {@code newDomain}
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial setDomain(lIntegersModulo newDomain) {
        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
        for (lMonomialTerm e : terms)
            newData.add(e.setDomain(newDomain));
        return new lMultivariatePolynomial(nVariables, newDomain, ordering, newData);
    }

    public lMultivariatePolynomial setDomainUnsafe(lIntegersModulo newDomain) {
        return new lMultivariatePolynomial(nVariables, newDomain, ordering, terms);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public lMultivariatePolynomial createConstant(long val) {
        MonomialsSet<lMonomialTerm> data = new MonomialsSet<>(ordering);
        if (val != 0)
            data.add(lMonomialTerm.withZeroExponents(nVariables, val));
        return new lMultivariatePolynomial(nVariables, domain, ordering, data);
    }

    @Override
    public lMultivariatePolynomial createZero() {
        return createConstant(0L);
    }

    @Override
    public lMultivariatePolynomial createOne() {
        return createConstant(1L);
    }

    /**
     * Creates linear polynomial of the form {@code cc + lc * variable}
     *
     * @param variable the variable
     * @param cc       the constant coefficient
     * @param lc       the leading coefficient
     * @return linear polynomial {@code cc + lc * variable}
     */
    public lMultivariatePolynomial createLinear(int variable, long cc, long lc) {
        MonomialsSet<lMonomialTerm> data = new MonomialsSet<>(ordering);

        lc = domain.modulus(lc);
        cc = domain.modulus(cc);

        if (cc != 0L)
            data.add(lMonomialTerm.withZeroExponents(nVariables, cc));
        if (lc != 0L) {
            int[] lcDegreeVector = new int[nVariables];
            lcDegreeVector[variable] = 1;
            data.add(new lMonomialTerm(lcDegreeVector, 1, lc));
        }
        return new lMultivariatePolynomial(nVariables, domain, ordering, data);
    }


    @Override
    public boolean isMonic() {
        return lc() == 1L;
    }

    @Override
    public int signum() {
        return Long.signum(lc());
    }

    @Override
    public boolean isZero() {
        return terms.size() == 0;
    }

    @Override
    public boolean isOne() {
        return size() == 1 && terms.first().coefficient == 1;
    }

    @Override
    public boolean isUnitCC() {
        return cc() == 1;
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && terms.first().isZeroVector());
    }


    /** {@inheritDoc} */
    @Override
    public lMonomialTerm lt() {
        return size() == 0 ? lMonomialTerm.withZeroExponents(nVariables, 0L) : terms.last();
    }

    /**
     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the ordering.
     *
     * @return leading coefficient of this polynomial
     */
    public long lc() {
        return lt().coefficient;
    }

    /**
     * Sets the leading coefficient to the specified value
     *
     * @param val new value for the lc
     * @return the leading coefficient to the specified value
     */
    public lMultivariatePolynomial setLC(long val) {
        if (isZero())
            return add(val);
        terms.add(lt().setCoefficient(domain.modulus(val)));
        return this;
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public long cc() {
        lMonomialTerm zero = lMonomialTerm.withZeroExponents(nVariables, 0L);
        return terms.getOrDefault(zero, zero).coefficient;
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public long content() {
        long gcd = -1;
        for (lMonomialTerm term : terms) {
            if (gcd == -1)
                gcd = term.coefficient;
            else
                gcd = LongArithmetics.gcd(gcd, term.coefficient);
        }
        return gcd;
    }

    /**
     * Returns array of polynomial coefficients
     *
     * @return array of polynomial coefficients
     */
    public long[] coefficients() {
        return terms.values().stream().mapToLong(x -> x.coefficient).toArray();
    }

    @Override
    public lMultivariatePolynomial primitivePart() {
        return divide(content());
    }

    @Override
    public lMultivariatePolynomial primitivePartSameSign() {
        return primitivePart();
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public lMultivariatePolynomial divide(long factor) {
        if (factor == 1)
            return this;
        return multiply(domain.reciprocal(factor)); // <- this is typically faster than the division
    }

    /** {@inheritDoc} */
    @Override
    public lMultivariatePolynomial divideOrNull(lMonomialTerm monomial) {
        if (monomial.isZeroVector())
            return divide(monomial.coefficient);
        long reciprocal = domain.reciprocal(monomial.coefficient);
        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
        for (lMonomialTerm term : terms) {
            lMonomialTerm dv = term.divide(monomial, domain.multiply(term.coefficient, reciprocal));
            if (dv == null)
                return null;
            map.add(dv);
        }
        loadFrom(map);
        release();
        return this;
    }

    /**
     * Makes this polynomial monic
     *
     * @return monic this
     */
    @Override
    public lMultivariatePolynomial monic() {
        return divide(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public lMultivariatePolynomial monic(long factor) {
        return multiply(domain.multiply(domain.modulus(factor), domain.reciprocal(lc())));
    }

    /**
     * Substitutes {@code value} for {@code variable}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    public lMultivariatePolynomial evaluate(int variable, long value) {
        value = domain.modulus(value);
        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, domain);
        for (lMonomialTerm el : terms) {
            long val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable, val));
        }
        return new lMultivariatePolynomial(nVariables, domain, ordering, newData);
    }

    /**
     * Substitutes {@code values} for {@code variables}.
     *
     * @param variables the variables
     * @param values    the values
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial evaluate(int[] variables, long[] values) {
        return evaluate(new lPrecomputedPowersHolder(values, domain), variables, ones(nVariables));
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
    lMultivariatePolynomial evaluate(lPrecomputedPowersHolder powers, int[] variables, int[] raiseFactors) {
        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
        for (lMonomialTerm el : terms) {
            lMonomialTerm r = el;
            long value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = domain.multiply(value, powers.pow(i, raiseFactors[i] * el.exponents[variables[i]]));
                r = r.setZero(variables[i], value);
            }

            add(newData, r);
        }
        return new lMultivariatePolynomial(nVariables, domain, ordering, newData);
    }

    /**
     * Evaluates this polynomial at specified points
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial[] evaluate(int variable, long... values) {
        return Arrays.stream(values).mapToObj(p -> evaluate(variable, p)).toArray(lMultivariatePolynomial[]::new);
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
     */
    public lMultivariatePolynomial eliminate(int variable, long value) {
        value = domain.modulus(value);
        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, domain);
        for (lMonomialTerm el : terms) {
            long val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.without(variable, val));
        }
        return new lMultivariatePolynomial(nVariables - 1, domain, ordering, newData);
    }

    private static final int SIZE_OF_POWERS_CACHE = 32;

    /** cached powers used to save some time */
    static final class lPrecomputedPowers {
        private final long value;
        private final lIntegersModulo domain;
        // TODO: store map here
        private final long[] precomputedPowers;

        lPrecomputedPowers(long value, lIntegersModulo domain) {
            this(SIZE_OF_POWERS_CACHE, value, domain);
        }

        lPrecomputedPowers(int cacheSize, long value, lIntegersModulo domain) {
            this.value = domain.modulus(value);
            this.domain = domain;
            this.precomputedPowers = new long[cacheSize];
            Arrays.fill(precomputedPowers, -1);
        }

        long pow(int exponent) {
            if (exponent >= SIZE_OF_POWERS_CACHE)
                return domain.powMod(value, exponent);

            if (precomputedPowers[exponent] != -1)
                return precomputedPowers[exponent];

            long result = 1L;
            long k2p = value;
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
    static final class lPrecomputedPowersHolder {
        private final int cacheSize;
        private final lIntegersModulo domain;
        private final lPrecomputedPowers[] powers;

        lPrecomputedPowersHolder(long[] points, lIntegersModulo domain) {
            this(SIZE_OF_POWERS_CACHE, points, domain);
        }

        @SuppressWarnings("unchecked")
        lPrecomputedPowersHolder(int cacheSize, long[] points, lIntegersModulo domain) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = new lPrecomputedPowers[points.length];
            for (int i = 0; i < points.length; i++)
                powers[i] = new lPrecomputedPowers(cacheSize, points[i], domain);
        }

        lPrecomputedPowersHolder(int cacheSize, lIntegersModulo domain, lPrecomputedPowers[] powers) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = powers;
        }

        void set(int i, long point) {
            if (powers[i] == null || powers[i].value != point)
                powers[i] = new lPrecomputedPowers(cacheSize, point, domain);
        }

        long pow(int i, int exponent) {
            return powers[i].pow(exponent);
        }
    }

    @Override
    public lMultivariatePolynomial negate() {
        for (Entry<DegreeVector, lMonomialTerm> entry : terms.entrySet()) {
            lMonomialTerm term = entry.getValue();
            entry.setValue(term.negate(domain));
        }
        release();
        return this;
    }

    @Override
    void add(MonomialsSet<lMonomialTerm> terms, lMonomialTerm term) {
        add(terms, term, domain);
    }

    @Override
    void subtract(MonomialsSet<lMonomialTerm> terms, lMonomialTerm term) {
        subtract(terms, term, domain);
    }

    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public lMultivariatePolynomial add(long oth) {
        oth = domain.modulus(oth);
        if (oth == 0)
            return this;
        add(terms, lMonomialTerm.withZeroExponents(nVariables, oth));
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public lMultivariatePolynomial subtract(long oth) {
        return add(domain.negate(domain.modulus(oth)));
    }

    @Override
    public lMultivariatePolynomial increment() {
        return add(1L);
    }

    @Override
    public lMultivariatePolynomial decrement() {
        return subtract(1L);
    }

    /**
     * Raises {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    @Override
    public lMultivariatePolynomial multiply(long factor) {
        factor = domain.modulus(factor);
        if (factor == 1)
            return this;
        if (factor == 0)
            return toZero();
        for (Entry<DegreeVector, lMonomialTerm> entry : terms.entrySet()) {
            lMonomialTerm term = entry.getValue();
            long val = domain.multiply(term.coefficient, factor);
            if (val != 0)
                entry.setValue(term.setCoefficient(val));
        }
        release();
        return this;
    }


    /** {@inheritDoc} */
    @Override
    public lMultivariatePolynomial multiply(lMonomialTerm term) {
        checkSameDomainWith(term);
        if (term.isZeroVector())
            return multiply(term.coefficient);
        if (term.coefficient == 0)
            return toZero();

        MonomialsSet<lMonomialTerm> newMap = new MonomialsSet<>(ordering);
        for (lMonomialTerm thisElement : terms) {
            long val = domain.multiply(thisElement.coefficient, term.coefficient);
            if (val != 0)
                newMap.add(thisElement.multiply(term, val));
        }

        return loadFrom(newMap);
    }

    @Override
    public lMultivariatePolynomial multiply(lMultivariatePolynomial oth) {
        checkSameDomainWith(oth);
        MonomialsSet<lMonomialTerm> newMap = new MonomialsSet<>(ordering);
        for (lMonomialTerm othElement : oth.terms)
            for (lMonomialTerm thisElement : terms)
                add(newMap, thisElement.multiply(othElement, domain.multiply(thisElement.coefficient, othElement.coefficient)));

        return loadFrom(newMap);
    }

    @Override
    public lMultivariatePolynomial square() {
        return multiply(this);
    }

    @Override
    lMultivariatePolynomial evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd) {
        //desired skeleton
        Set<DegreeVector> skeleton = getSkeletonExcept(variable);
        lMultivariatePolynomial tmp;
        do {
            long randomPoint = domain.randomElement(rnd);
            tmp = evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    @Override
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomial clone() {
        return new lMultivariatePolynomial(nVariables, domain, ordering, terms.clone());
    }

    public String toString(String... vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (lMonomialTerm term : terms) {
            long coeff = term.coefficient;
            if (coeff == 0)
                continue;
            String monomialString = term.toString(vars);
            if (first) {
                if (coeff != 1 || monomialString.isEmpty()) {
                    sb.append(coeff);
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
                first = false;
            } else {
                if (coeff > 0)
                    sb.append("+");
                else {
                    sb.append("-");
                    coeff = domain.negate(coeff);
                }

                if (coeff != 1 || monomialString.isEmpty()) {
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

    static String[] defaultVars(int nVars) {
        char v = 'a';
        String[] vars = new String[nVars];
        for (int i = 0; i < nVars; i++)
            vars[i] = Character.toString(v++);
        return vars;
    }
}
