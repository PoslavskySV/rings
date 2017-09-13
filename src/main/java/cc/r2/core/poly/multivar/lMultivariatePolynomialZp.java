package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.*;
import cc.r2.core.poly.univar.lUnivariatePolynomialZ;
import cc.r2.core.poly.univar.lUnivariatePolynomialZp;
import cc.r2.core.util.ArraysUtil;
import cc.redberry.libdivide4j.FastDivision;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Function;
import java.util.function.LongFunction;

/**
 * Multivariate polynomial over Zp domain with the modulus in the range (0, 2^62) (see {@link MachineArithmetic#MAX_SUPPORTED_MODULUS}).
 * For details on polynomial data structure and properties see {@link AMultivariatePolynomial}.
 *
 * @author Stanislav Poslavsky
 * @see AMultivariatePolynomial
 * @since 1.0
 */
public final class lMultivariatePolynomialZp extends AMultivariatePolynomial<lMonomialZp, lMultivariatePolynomialZp> {
    private static final long serialVersionUID = 1L;
    /**
     * The domain.
     */
    public final lIntegersModulo domain;

    private lMultivariatePolynomialZp(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, MonomialSet<lMonomialZp> lMonomialTerms) {
        super(nVariables, ordering, lMonomialTerms);
        this.domain = domain;
    }

    private static void add(MonomialSet<lMonomialZp> polynomial, lMonomialZp term, lIntegersModulo domain) {
        if (term.coefficient == 0)
            return;
        lMonomialZp pTerm = polynomial.get(term);
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

    private static void subtract(MonomialSet<lMonomialZp> polynomial, lMonomialZp term, lIntegersModulo domain) {
        if (term.coefficient == 0)
            return;
        lMonomialZp pTerm = polynomial.get(term);
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
     * Creates multivariate polynomial from a set of monomials
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   term ordering
     * @param terms      the monomials
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp create(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, MonomialSet<lMonomialZp> terms) {
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, terms);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param domain   the domain
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp create(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, Iterable<lMonomialZp> terms) {
        MonomialSet<lMonomialZp> map = new MonomialSet<>(ordering);
        for (lMonomialZp term : terms)
            add(map, term.setDomain(domain), domain);

        return new lMultivariatePolynomialZp(nVariables, domain, ordering, map);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param domain   the domain
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp create(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering, lMonomialZp... terms) {
        return create(nVariables, domain, ordering, Arrays.asList(terms));
    }

    /**
     * Creates zero polynomial.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return zero
     */
    public static lMultivariatePolynomialZp zero(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, new MonomialSet<>(ordering));
    }

    /**
     * Creates unit polynomial.
     *
     * @param nVariables number of variables
     * @param domain     the domain
     * @param ordering   the ordering
     * @return unit
     */
    public static lMultivariatePolynomialZp one(int nVariables, lIntegersModulo domain, Comparator<DegreeVector> ordering) {
        return create(nVariables, domain, ordering, lMonomialZp.withZeroExponents(nVariables, 1L));
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    the string
     * @param domain    the domain
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp parse(String string, lIntegersModulo domain, String... variables) {
        return MultivariatePolynomial.asLongPolyZp(Parser.parse(string, domain.asDomain(), MonomialOrder.LEX, variables));
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    the string
     * @param domain    the domain
     * @param ordering  monomial order
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp parse(String string, lIntegersModulo domain, Comparator<DegreeVector> ordering, String... variables) {
        return MultivariatePolynomial.asLongPolyZp(Parser.parse(string, domain.asDomain(), ordering, variables));
    }

    /**
     * Converts univariate polynomial to multivariate.
     *
     * @param poly       univariate polynomial
     * @param nVariables number of variables in the result
     * @param variable   variable that will be used as a primary variable
     * @param ordering   the ordering
     * @return multivariate polynomial
     */
    public static lMultivariatePolynomialZp asMultivariate(lUnivariatePolynomialZp poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        MonomialSet<lMonomialZp> map = new MonomialSet<>(ordering);
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;
            map.add(new lMonomialZp(degreeVector, i, poly.get(i)));
        }
        return new lMultivariatePolynomialZp(nVariables, poly.domain, ordering, map);
    }

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
        for (lMonomialZp e : terms)
            univarData[e.exponents[theVar]] = e.coefficient;
        return lUnivariatePolynomialZp.createUnsafe(domain, univarData);
    }

    @Override
    public MultivariatePolynomial<lUnivariatePolynomialZp> asOverUnivariate(int variable) {
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(domain);
        UnivariatePolynomials<lUnivariatePolynomialZp> pDomain = new UnivariatePolynomials<>(factory);
        MonomialSet<Monomial<lUnivariatePolynomialZp>> newData = new MonomialSet<>(ordering);
        for (lMonomialZp e : terms) {
            Monomial<lUnivariatePolynomialZp> eMonomial = new Monomial<>(
                    e.set(variable, 0).exponents,
                    factory.createMonomial(e.coefficient, e.exponents[variable]));
            MultivariatePolynomial.add(newData, eMonomial, pDomain);
        }
        return new MultivariatePolynomial<>(nVariables, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<lUnivariatePolynomialZp> asOverUnivariateEliminate(int variable) {
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(domain);
        UnivariatePolynomials<lUnivariatePolynomialZp> pDomain = new UnivariatePolynomials<>(factory);
        MonomialSet<Monomial<lUnivariatePolynomialZp>> newData = new MonomialSet<>(ordering);
        for (lMonomialZp e : terms) {
            Monomial<lUnivariatePolynomialZp> eMonomial = new Monomial<>(
                    e.without(variable).exponents,
                    factory.createMonomial(e.coefficient, e.exponents[variable]));
            MultivariatePolynomial.add(newData, eMonomial, pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<lMultivariatePolynomialZp> asOverMultivariate(int... variables) {
        Domain<lMultivariatePolynomialZp> domain = new MultivariatePolynomials<>(this);
        MonomialSet<Monomial<lMultivariatePolynomialZp>> terms = new MonomialSet<>(ordering);
        for (lMonomialZp term : this) {
            int[] coeffExponents = new int[nVariables];
            for (int var : variables)
                coeffExponents[var] = term.exponents[var];
            lMonomialZp restTerm = term.setZero(variables, 1);
            Monomial<lMultivariatePolynomialZp> newTerm
                    = new Monomial<>(
                    restTerm.exponents,
                    restTerm.totalDegree,
                    create(new lMonomialZp(coeffExponents, term.coefficient)));

            MultivariatePolynomial.add(terms, newTerm, domain);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms);
    }

    @Override
    public MultivariatePolynomial<lMultivariatePolynomialZp> asOverMultivariateEliminate(int... variables) {
        variables = variables.clone();
        Arrays.sort(variables);
        int[] restVariables = ArraysUtil.intSetDifference(ArraysUtil.sequence(nVariables), variables);
        Domain<lMultivariatePolynomialZp> domain = new MultivariatePolynomials<>(create(variables.length, new MonomialSet<>(ordering)));
        MonomialSet<Monomial<lMultivariatePolynomialZp>> terms = new MonomialSet<>(ordering);
        for (lMonomialZp term : this) {
            int i = 0;
            int[] coeffExponents = new int[variables.length];
            for (int var : variables)
                coeffExponents[i++] = term.exponents[var];

            i = 0;
            int[] termExponents = new int[restVariables.length];
            for (int var : restVariables)
                termExponents[i++] = term.exponents[var];

            Monomial<lMultivariatePolynomialZp> newTerm
                    = new Monomial<>(
                    termExponents,
                    create(variables.length, this.domain, this.ordering, new lMonomialZp(coeffExponents, term.coefficient)));

            MultivariatePolynomial.add(terms, newTerm, domain);
        }
        return new MultivariatePolynomial<>(nVariables, domain, ordering, terms);
    }

    /**
     * Converts multivariate polynomial over univariate polynomial domain (Zp[variable][other_variables]) to a
     * multivariate polynomial over coefficient domain (Zp[all_variables])
     *
     * @param poly     the polynomial
     * @param variable the variable to insert
     * @return multivariate polynomial over normal coefficient domain
     */
    public static lMultivariatePolynomialZp asNormalMultivariate(MultivariatePolynomial<lUnivariatePolynomialZp> poly, int variable) {
        lIntegersModulo domain = poly.domain.getZero().domain;
        int nVariables = poly.nVariables + 1;
        lMultivariatePolynomialZp result = zero(nVariables, domain, poly.ordering);
        for (Monomial<lUnivariatePolynomialZp> entry : poly.terms) {
            lUnivariatePolynomialZp uPoly = entry.coefficient;
            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new lMonomialZp(cdv, uPoly.get(i)));
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
    public static lMultivariatePolynomialZp asNormalMultivariate(MultivariatePolynomial<lMultivariatePolynomialZp> poly) {
        lIntegersModulo domain = poly.domain.getZero().domain;
        int nVariables = poly.nVariables;
        lMultivariatePolynomialZp result = zero(nVariables, domain, poly.ordering);
        for (Monomial<lMultivariatePolynomialZp> term : poly.terms) {
            lMultivariatePolynomialZp uPoly = term.coefficient;
            result.add(uPoly.clone().multiply(new lMonomialZp(term.exponents, term.totalDegree, 1)));
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
    public static lMultivariatePolynomialZp asNormalMultivariate(
            MultivariatePolynomial<lMultivariatePolynomialZp> poly,
            int[] coefficientVariables, int[] mainVariables) {
        lIntegersModulo domain = poly.domain.getZero().domain;
        int nVariables = coefficientVariables.length + mainVariables.length;
        lMultivariatePolynomialZp result = zero(nVariables, domain, poly.ordering);
        for (Monomial<lMultivariatePolynomialZp> term : poly.terms) {
            lMultivariatePolynomialZp coefficient =
                    term.coefficient.joinNewVariables(nVariables, coefficientVariables);
            Monomial<lMultivariatePolynomialZp> t = term.joinNewVariables(nVariables, mainVariables);
            result.add(coefficient.multiply(new lMonomialZp(t.exponents, t.totalDegree, 1)));
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
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (lMonomialZp t : this)
            bTerms.add(new Monomial<>(t.exponents, t.totalDegree, BigInteger.valueOf(domain.symmetricForm(t.coefficient))));
        return new MultivariatePolynomial<>(nVariables, Domains.Z, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> asPolyZ() {
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (lMonomialZp t : this)
            bTerms.add(t.toBigMonomial());
        return new MultivariatePolynomial<>(nVariables, Domains.Z, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> toBigPoly() {
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (lMonomialZp t : this)
            bTerms.add(t.toBigMonomial());
        return new MultivariatePolynomial<>(nVariables, domain.asDomain(), ordering, bTerms);
    }

    /* ============================================ Main methods ============================================ */

    @Override
    public lMultivariatePolynomialZp contentAsPoly() {
        return createConstant(content());
    }

    @Override
    public lMultivariatePolynomialZp lcAsPoly() {
        return createConstant(lc());
    }

    @Override
    public lMultivariatePolynomialZp ccAsPoly() {
        return createConstant(cc());
    }

    @Override
    lMultivariatePolynomialZp create(int nVariables, MonomialSet<lMonomialZp> lMonomialTerms) {
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, lMonomialTerms);
    }

    @Override
    public lMonomialZp createTermWithUnitCoefficient(int[] exponents) {
        return new lMonomialZp(exponents, 1L);
    }

    @Override
    public lMonomialZp createUnitTerm() {
        return lMonomialZp.withZeroExponents(nVariables, 1L);
    }

    @Override
    public lMonomialZp createZeroTerm() {
        return lMonomialZp.withZeroExponents(nVariables, 0L);
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
    public boolean isOverPerfectPower() {
        return domain.isPerfectPower();
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerBase() {
        return BigInteger.valueOf(domain.perfectPowerBase());
    }

    @Override
    public BigInteger coefficientDomainPerfectPowerExponent() {
        return BigInteger.valueOf(domain.perfectPowerExponent());
    }

    @Override
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomialZp[] createArray(int length) {return new lMultivariatePolynomialZp[length];}

    @Override
    public lMultivariatePolynomialZp[][] createArray2d(int length) {
        return new lMultivariatePolynomialZp[length][];
    }

    @Override
    public lMultivariatePolynomialZp[][] createArray2d(int length1, int length2) {
        return new lMultivariatePolynomialZp[length1][length2];
    }

    @Override
    public boolean sameDomainWith(lMultivariatePolynomialZp oth) {
        return nVariables == oth.nVariables && domain.equals(oth.domain);
    }

    @Override
    public lMultivariatePolynomialZp setDomainFrom(lMultivariatePolynomialZp lMonomialTerms) {
        return clone();
    }

    @Override
    boolean isZeroMonomial(lMonomialZp a) {
        return a.coefficient == 0L;
    }

    @Override
    lMonomialZp multiply(lMonomialZp a, lMonomialZp b) {
        return a.multiply(b, domain.multiply(a.coefficient, b.coefficient));
    }

    @Override
    lMonomialZp divideOrNull(lMonomialZp a, lMonomialZp b) {
        return a.divide(b, domain.divide(a.coefficient, b.coefficient));
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
    public lMultivariatePolynomialZp setDomain(long newModulus) {
        return setDomain(new lIntegersModulo(newModulus));
    }

    /**
     * Switches to another domain specified by {@code newDomain}
     *
     * @param newDomain the new domain
     * @return a copy of this reduced to the domain specified by {@code newDomain}
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomialZp setDomain(lIntegersModulo newDomain) {
        MonomialSet<lMonomialZp> newData = new MonomialSet<>(ordering);
        for (lMonomialZp e : terms)
            add(newData, e.setDomain(newDomain));
        return new lMultivariatePolynomialZp(nVariables, newDomain, ordering, newData);
    }

    /**
     * Switches to another domain specified by {@code newDomain}
     *
     * @param newDomain the new domain
     * @return a copy of this reduced to the domain specified by {@code newDomain}
     */
    @SuppressWarnings("unchecked")
    public <E> MultivariatePolynomial<E> setDomain(Domain<E> newDomain) {
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        for (lMonomialZp e : terms)
            MultivariatePolynomial.add(newData, e.setDomain(newDomain), newDomain);
        return new MultivariatePolynomial(nVariables, newDomain, ordering, newData);
    }

    /** internal API */
    public lMultivariatePolynomialZp setDomainUnsafe(lIntegersModulo newDomain) {
        return new lMultivariatePolynomialZp(nVariables, newDomain, ordering, terms);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public lMultivariatePolynomialZp createConstant(long val) {
        MonomialSet<lMonomialZp> data = new MonomialSet<>(ordering);
        if (val != 0)
            data.add(lMonomialZp.withZeroExponents(nVariables, val));
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, data);
    }

    @Override
    public lMultivariatePolynomialZp createZero() {
        return createConstant(0L);
    }

    @Override
    public lMultivariatePolynomialZp createOne() {
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
    public lMultivariatePolynomialZp createLinear(int variable, long cc, long lc) {
        MonomialSet<lMonomialZp> data = new MonomialSet<>(ordering);

        lc = domain.modulus(lc);
        cc = domain.modulus(cc);

        if (cc != 0L)
            data.add(lMonomialZp.withZeroExponents(nVariables, cc));
        if (lc != 0L) {
            int[] lcDegreeVector = new int[nVariables];
            lcDegreeVector[variable] = 1;
            data.add(new lMonomialZp(lcDegreeVector, 1, lc));
        }
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, data);
    }


    @Override
    public boolean isMonic() {
        return lc() == 1L;
    }

    @Override
    public int signumOfLC() {
        return Long.signum(lc());
    }

    @Override
    public boolean isZero() {
        return terms.size() == 0;
    }

    @Override
    public boolean isOne() {
        if (size() != 1)
            return false;
        lMonomialZp lt = terms.first();
        return lt.isZeroVector() && lt.coefficient == 1L;
    }

    @Override
    public boolean isUnitCC() {
        return cc() == 1;
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && terms.first().isZeroVector());
    }


    @Override
    public lMonomialZp lt() {
        return size() == 0 ? lMonomialZp.withZeroExponents(nVariables, 0L) : terms.last();
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
    public lMultivariatePolynomialZp setLC(long val) {
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
        lMonomialZp zero = lMonomialZp.withZeroExponents(nVariables, 0L);
        return terms.getOrDefault(zero, zero).coefficient;
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public long content() {
        long gcd = -1;
        for (lMonomialZp term : terms) {
            if (gcd == -1)
                gcd = term.coefficient;
            else
                gcd = MachineArithmetic.gcd(gcd, term.coefficient);
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
    public lMultivariatePolynomialZp primitivePart(int variable) {
        return asNormalMultivariate(asOverUnivariateEliminate(variable).primitivePart(), variable);
    }

    @Override
    public lUnivariatePolynomialZp contentUnivariate(int variable) {
        return asOverUnivariate(variable).content();
    }

    @Override
    public lMultivariatePolynomialZp primitivePart() {
        return divide(content());
    }

    @Override
    public lMultivariatePolynomialZp primitivePartSameSign() {
        return primitivePart();
    }

    @Override
    public lMultivariatePolynomialZp divideByLC(lMultivariatePolynomialZp other) {
        return divide(other.lc());
    }

    /**
     * Divides this polynomial by a {@code factor}
     *
     * @param factor the factor
     * @return {@code this / factor}
     */
    public lMultivariatePolynomialZp divide(long factor) {
        if (factor == 1)
            return this;
        return multiply(domain.reciprocal(factor)); // <- this is typically faster than the division
    }

    @Override
    public lMultivariatePolynomialZp divideOrNull(lMonomialZp monomial) {
        if (monomial.isZeroVector())
            return divide(monomial.coefficient);
        long reciprocal = domain.reciprocal(monomial.coefficient);
        MonomialSet<lMonomialZp> map = new MonomialSet<>(ordering);
        for (lMonomialZp term : terms) {
            lMonomialZp dv = term.divide(monomial, domain.multiply(term.coefficient, reciprocal));
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
    public lMultivariatePolynomialZp monic() {
        if (isMonic())
            return this;
        if (isZero())
            return this;
        return divide(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is
     * {@code monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public lMultivariatePolynomialZp monic(long factor) {
        return multiply(domain.multiply(domain.modulus(factor), domain.reciprocal(lc())));
    }

    @Override
    public lMultivariatePolynomialZp monicWithLC(lMultivariatePolynomialZp other) {
        return monic(other.lc());
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    public lMultivariatePolynomialZp evaluate(int variable, long value) {
        value = domain.modulus(value);
        if (value == 0)
            return evaluateAtZero(variable);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, domain);
        return evaluate(variable, powers);
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    lMultivariatePolynomialZp evaluate(int variable, lPrecomputedPowers powers) {
        if (powers.value == 0)
            return evaluateAtZero(variable);
        MonomialSet<lMonomialZp> newData = new MonomialSet<>(ordering);
        for (lMonomialZp el : terms) {
            long val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable, val));
        }
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, newData);
    }

    lUnivariatePolynomialZp evaluateAtZeroAllExcept(int variable) {
        long[] uData = new long[degree(variable) + 1];
        for (lMonomialZp el : terms) {
            if (el.totalDegree != el.exponents[variable])
                continue;
            int uExp = el.exponents[variable];
            uData[uExp] = domain.add(uData[uExp], el.coefficient);
        }
        return lUnivariatePolynomialZp.createUnsafe(domain, uData);
    }

    /**
     * Returns a copy of this with {@code values} substituted for {@code variables}
     *
     * @param variables the variables
     * @param values    the values
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomialZp evaluate(int[] variables, long[] values) {
        for (long value : values)
            if (value != 0)
                return evaluate(new lPrecomputedPowersHolder(nVariables, variables, values, domain), variables, ones(nVariables));

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
    lMultivariatePolynomialZp evaluate(lPrecomputedPowersHolder powers, int[] variables) {
        return evaluate(powers, variables, ones(variables.length));
    }

    /** substitutes {@code values} for {@code variables} */
    @SuppressWarnings("unchecked")
    lMultivariatePolynomialZp evaluate(lPrecomputedPowersHolder powers, int[] variables, int[] raiseFactors) {
        MonomialSet<lMonomialZp> newData = new MonomialSet<>(ordering);
        for (lMonomialZp el : terms) {
            lMonomialZp r = el;
            long value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = domain.multiply(value, powers.pow(variables[i], raiseFactors[i] * el.exponents[variables[i]]));
                r = r.setZero(variables[i], value);
            }

            add(newData, r);
        }
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, newData);
    }

    /**
     * Evaluates this polynomial at specified points
     */
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomialZp[] evaluate(int variable, long... values) {
        return Arrays.stream(values).mapToObj(p -> evaluate(variable, p)).toArray(lMultivariatePolynomialZp[]::new);
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
     */
    public lMultivariatePolynomialZp eliminate(int variable, long value) {
        value = domain.modulus(value);
        MonomialSet<lMonomialZp> newData = new MonomialSet<>(ordering);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, domain);
        for (lMonomialZp el : terms) {
            long val = domain.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.without(variable, val));
        }
        return new lMultivariatePolynomialZp(nVariables - 1, domain, ordering, newData);
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
            if (exponent >= precomputedPowers.length)
                return domain.powMod(value, exponent);

            if (precomputedPowers[exponent] != -1)
                return precomputedPowers[exponent];

            long result = 1L;
            long k2p = value;
            int rExp = 0, kExp = 1;
            for (; ; ) {
                if ((exponent & 1) != 0)
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
        final int cacheSize;
        final lIntegersModulo domain;
        final lPrecomputedPowers[] powers;

        lPrecomputedPowersHolder(int nVariables, int[] variables, long[] points, lIntegersModulo domain) {
            this(SIZE_OF_POWERS_CACHE, nVariables, variables, points, domain);
        }

        @SuppressWarnings("unchecked")
        lPrecomputedPowersHolder(int cacheSize, int nVariables, int[] variables, long[] points, lIntegersModulo domain) {
            this.cacheSize = cacheSize;
            this.domain = domain;
            this.powers = new lPrecomputedPowers[nVariables];
            for (int i = 0; i < variables.length; i++)
                powers[variables[i]] = new lPrecomputedPowers(cacheSize, points[i], domain);
        }

        void set(int i, long point) {
            if (powers[i] == null || powers[i].value != point)
                powers[i] = new lPrecomputedPowers(cacheSize, point, domain);
        }

        long pow(int variable, int exponent) {
            return powers[variable].pow(exponent);
        }
    }

    /**
     * Returns a copy of this with {@code poly} substituted for {@code variable}
     *
     * @param variable the variable
     * @param poly     the replacement for the variable
     * @return a copy of this with  {@code variable -> poly}
     */
    public lMultivariatePolynomialZp substitute(int variable, lMultivariatePolynomialZp poly) {
        if (poly.isConstant())
            return evaluate(variable, poly.cc());
        lPrecomputedSubstitution subsPowers;
        if (poly.isEffectiveUnivariate())
            subsPowers = new lUSubstitution(poly.asUnivariate(), poly.univariateVariable(), nVariables, ordering);
        else
            subsPowers = new lMSubstitution(poly);

        lMultivariatePolynomialZp result = createZero();
        for (lMonomialZp term : this) {
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
     * Returns a copy of this with {@code variable -> variable + shift}
     *
     * @param variable the variable
     * @param shift    shift amount
     * @return a copy of this with  {@code variable -> variable + shift}
     */
    public lMultivariatePolynomialZp shift(int variable, long shift) {
        if (shift == 0)
            return clone();

        lUSubstitution shifts = new lUSubstitution(lUnivariatePolynomialZ.create(shift, 1).modulus(domain), variable, nVariables, ordering);
        lMultivariatePolynomialZp result = createZero();
        for (lMonomialZp term : this) {
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
    public lMultivariatePolynomialZp shift(int[] variables, long[] shifts) {
        lPrecomputedSubstitution[] powers = new lPrecomputedSubstitution[nVariables];
        boolean allShiftsAreZero = true;
        for (int i = 0; i < variables.length; ++i) {
            if (shifts[i] != 0)
                allShiftsAreZero = false;
            powers[variables[i]] = new lUSubstitution(lUnivariatePolynomialZ.create(shifts[i], 1).modulus(domain, false), variables[i], nVariables, ordering);
        }

        if (allShiftsAreZero)
            return clone();

        lPrecomputedSubstitutions calculatedShifts = new lPrecomputedSubstitutions(powers);

        lMultivariatePolynomialZp result = createZero();
        for (lMonomialZp term : this) {
            lMultivariatePolynomialZp temp = createOne();
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

    static final class lPrecomputedSubstitutions {
        final lPrecomputedSubstitution[] subs;

        public lPrecomputedSubstitutions(lPrecomputedSubstitution[] subs) {
            this.subs = subs;
        }

        lMultivariatePolynomialZp getSubstitutionPower(int var, int exponent) {
            if (subs[var] == null)
                throw new IllegalArgumentException();

            return subs[var].pow(exponent);
        }
    }

    interface lPrecomputedSubstitution {
        lMultivariatePolynomialZp pow(int exponent);
    }

    static final class lUSubstitution implements lPrecomputedSubstitution {
        final int variable;
        final int nVariables;
        final Comparator<DegreeVector> ordering;
        final lUnivariatePolynomialZp base;
        final TIntObjectHashMap<lUnivariatePolynomialZp> uCache = new TIntObjectHashMap<>();
        final TIntObjectHashMap<lMultivariatePolynomialZp> mCache = new TIntObjectHashMap<>();

        lUSubstitution(lUnivariatePolynomialZp base, int variable, int nVariables, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.variable = variable;
            this.ordering = ordering;
            this.base = base;
        }

        @Override
        public lMultivariatePolynomialZp pow(int exponent) {
            lMultivariatePolynomialZp cached = mCache.get(exponent);
            if (cached != null)
                return cached.clone();
            lUnivariatePolynomialZp r = PolynomialMethods.polyPow(base, exponent, true, uCache);
            mCache.put(exponent, cached = asMultivariate(r, nVariables, variable, ordering));
            return cached.clone();
        }
    }

    static final class lMSubstitution implements lPrecomputedSubstitution {
        final lMultivariatePolynomialZp base;
        final TIntObjectHashMap<lMultivariatePolynomialZp> cache = new TIntObjectHashMap<>();

        lMSubstitution(lMultivariatePolynomialZp base) {
            this.base = base;
        }

        @Override
        public lMultivariatePolynomialZp pow(int exponent) {
            return PolynomialMethods.polyPow(base, exponent, true, cache);
        }
    }

    @Override
    public lMultivariatePolynomialZp negate() {
        for (Entry<DegreeVector, lMonomialZp> entry : terms.entrySet()) {
            lMonomialZp term = entry.getValue();
            entry.setValue(term.negate(domain));
        }
        release();
        return this;
    }

    @Override
    void add(MonomialSet<lMonomialZp> terms, lMonomialZp term) {
        add(terms, term, domain);
    }

    @Override
    void subtract(MonomialSet<lMonomialZp> terms, lMonomialZp term) {
        subtract(terms, term, domain);
    }

    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public lMultivariatePolynomialZp add(long oth) {
        oth = domain.modulus(oth);
        if (oth == 0)
            return this;
        add(terms, lMonomialZp.withZeroExponents(nVariables, oth));
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public lMultivariatePolynomialZp subtract(long oth) {
        return add(domain.negate(domain.modulus(oth)));
    }

    @Override
    public lMultivariatePolynomialZp increment() {
        return add(1L);
    }

    @Override
    public lMultivariatePolynomialZp decrement() {
        return subtract(1L);
    }

    @Override
    public lMultivariatePolynomialZp multiply(long factor) {
        factor = domain.modulus(factor);
        if (factor == 1)
            return this;
        if (factor == 0)
            return toZero();
        for (Entry<DegreeVector, lMonomialZp> entry : terms.entrySet()) {
            lMonomialZp term = entry.getValue();
            long val = domain.multiply(term.coefficient, factor);
            if (val != 0)
                entry.setValue(term.setCoefficient(val));
        }
        release();
        return this;
    }

    @Override
    public lMultivariatePolynomialZp multiplyByLC(lMultivariatePolynomialZp other) {
        return multiply(other.lc());
    }

    @Override
    public lMultivariatePolynomialZp multiplyByBigInteger(BigInteger factor) {
        return multiply(factor.mod(BigInteger.valueOf(domain.modulus)).longValueExact());
    }

    @Override
    public lMultivariatePolynomialZp multiply(lMonomialZp monomial) {
        checkSameDomainWith(monomial);
        if (monomial.isZeroVector())
            return multiply(monomial.coefficient);
        if (monomial.coefficient == 0)
            return toZero();

        MonomialSet<lMonomialZp> newMap = new MonomialSet<>(ordering);
        for (lMonomialZp thisElement : terms) {
            long val = domain.multiply(thisElement.coefficient, monomial.coefficient);
            if (val != 0)
                newMap.add(thisElement.multiply(monomial, val));
        }

        return loadFrom(newMap);
    }

    @Override
    public lMultivariatePolynomialZp multiply(lMultivariatePolynomialZp oth) {
        assertSameDomainWith(oth);
        if (oth.isZero())
            return toZero();
        if (isZero())
            return this;
        if (oth.isConstant())
            return multiply(oth.cc());
        MonomialSet<lMonomialZp> newMap = new MonomialSet<>(ordering);
        for (lMonomialZp othElement : oth.terms)
            for (lMonomialZp thisElement : terms)
                add(newMap, thisElement.multiply(othElement, domain.multiply(thisElement.coefficient, othElement.coefficient)));

        return loadFrom(newMap);
    }

    @Override
    public lMultivariatePolynomialZp square() {
        return multiply(this);
    }

    @Override
    public lMultivariatePolynomialZp evaluateAtRandom(int variable, RandomGenerator rnd) {
        return evaluate(variable, domain.randomElement(rnd));
    }

    @Override
    public lMultivariatePolynomialZp evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd) {
        //desired skeleton
        Set<DegreeVector> skeleton = getSkeletonExcept(variable);
        lMultivariatePolynomialZp tmp;
        do {
            long randomPoint = domain.randomElement(rnd);
            tmp = evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    @Override
    public lMultivariatePolynomialZp derivative(int variable, int order) {
        if (order == 0)
            return this.clone();
        if (isConstant())
            return createZero();
        MonomialSet<lMonomialZp> newTerms = new MonomialSet<>(ordering);
        for (lMonomialZp term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;
            long newCoefficient = term.coefficient;
            for (int i = 0; i < order; ++i)
                newCoefficient = domain.multiply(newCoefficient, exponent - i);
            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            add(newTerms, new lMonomialZp(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, newTerms);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static long seriesCoefficientFactor(int exponent, int order, lIntegersModulo domain) {
        lIntegersModulo baseDomain = domain.perfectPowerBaseDomain();
        if (order < baseDomain.modulus) {
            long factor = 1;
            for (int i = 0; i < order; ++i)
                factor = domain.multiply(factor, exponent - i);
            factor = domain.divide(factor, domain.factorial(order));
            return factor;
        }

        long numerator = 1, denominator = 1;
        int numZeros = 0, denZeros = 0;
        for (int i = 1; i <= order; ++i) {
            long
                    num = exponent - i + 1,
                    numMod = baseDomain.modulus(num);

            while (num > 1 && numMod == 0) {
                num = FastDivision.divideSignedFast(num, baseDomain.magic);
                numMod = baseDomain.modulus(num);
                ++numZeros;
            }

            if (numMod != 0)
                numerator = domain.multiply(numerator, domain == baseDomain ? numMod : domain.modulus(num));

            long
                    den = i,
                    denMod = baseDomain.modulus(i);

            while (den > 1 && denMod == 0) {
                den = FastDivision.divideSignedFast(den, baseDomain.magic);
                denMod = baseDomain.modulus(den);
                ++denZeros;
            }

            if (denMod != 0)
                denominator = domain.multiply(denominator, domain == baseDomain ? denMod : domain.modulus(den));
        }

        if (numZeros > denZeros)
            numerator = domain.multiply(numerator, domain.powMod(baseDomain.modulus, numZeros - denZeros));
        else if (denZeros < numZeros)
            denominator = domain.multiply(denominator, domain.powMod(baseDomain.modulus, denZeros - numZeros));

        if (numerator == 0)
            return numerator;
        return domain.divide(numerator, denominator);
    }

    @Override
    public lMultivariatePolynomialZp seriesCoefficient(int variable, int order) {
        if (order == 0)
            return this.clone();
        if (isConstant())
            return createZero();

        MonomialSet<lMonomialZp> newTerms = new MonomialSet<>(ordering);
        for (lMonomialZp term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;

            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            long newCoefficient = domain.multiply(term.coefficient, seriesCoefficientFactor(exponent, order, domain));
            add(newTerms, new lMonomialZp(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, newTerms);
    }

    /**
     * Maps terms of this using specified mapping function
     *
     * @param newDomain the new domain
     * @param mapper    mapping
     * @param <T>       new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms
     */
    public <T> MultivariatePolynomial<T> mapTerms(Domain<T> newDomain, Function<lMonomialZp, Monomial<T>> mapper) {
        return terms.values()
                .stream()
                .map(mapper)
                .collect(new PolynomialCollector<>(() -> MultivariatePolynomial.zero(nVariables, newDomain, ordering)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newDomain the new domain
     * @param mapper    mapping
     * @param <T>       new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public <T> MultivariatePolynomial<T> mapCoefficients(Domain<T> newDomain, LongFunction<T> mapper) {
        return mapTerms(newDomain, t -> new Monomial<>(t.exponents, t.totalDegree, mapper.apply(t.coefficient)));
    }

    @Override
    public int compareTo(lMultivariatePolynomialZp oth) {
        int c = Integer.compare(size(), oth.size());
        if (c != 0)
            return c;
        Iterator<lMonomialZp>
                thisIt = iterator(),
                othIt = oth.iterator();

        while (thisIt.hasNext() && othIt.hasNext()) {
            lMonomialZp
                    a = thisIt.next(),
                    b = othIt.next();

            if ((c = ordering.compare(a, b)) != 0)
                return c;

            if ((c = Long.compare(a.coefficient, b.coefficient)) != 0)
                return c;
        }
        return 0;
    }

    @Override
    @SuppressWarnings("unchecked")
    public lMultivariatePolynomialZp clone() {
        return new lMultivariatePolynomialZp(nVariables, domain, ordering, terms.clone());
    }

    @Override
    public lMultivariatePolynomialZp parsePoly(String string) {
        lMultivariatePolynomialZp r = parse(string, domain, ordering, defaultVars(nVariables));
        if (r.nVariables != nVariables)
            throw new IllegalArgumentException("not from this field");
        return r;
    }

    @Override
    public lMultivariatePolynomialZp parsePoly(String string, String[] variables) {
        lMultivariatePolynomialZp r = parse(string, domain, ordering, variables);
        if (r.nVariables != nVariables)
            throw new IllegalArgumentException("not from this field");
        return r;
    }

    @Override
    public String toString(String... vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (lMonomialZp term : terms) {
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
    public String coefficientDomainToString() {
        return domain.toString();
    }
}
