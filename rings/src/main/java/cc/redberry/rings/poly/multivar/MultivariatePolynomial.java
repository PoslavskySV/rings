package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Function;
import java.util.function.ToLongFunction;
import java.util.stream.Stream;

/**
 * @since 1.0
 */
public final class MultivariatePolynomial<E> extends AMultivariatePolynomial<Monomial<E>, MultivariatePolynomial<E>> {
    private static final long serialVersionUID = 1L;
    /** The coefficient ring */
    public final Ring<E> ring;

    MultivariatePolynomial(int nVariables, Ring<E> ring, Comparator<DegreeVector> ordering, MonomialSet<Monomial<E>> terms) {
        super(nVariables, ordering, terms);
        this.ring = ring;
    }

    /* ============================================ Factory methods ============================================ */

    static <E> void add(MonomialSet<Monomial<E>> polynomial, Monomial<E> term, Ring<E> ring) {
        if (ring.isZero(term.coefficient))
            return;
        polynomial.merge(term, term, (o, n) -> {
            E r = ring.add(o.coefficient, n.coefficient);
            if (ring.isZero(r))
                return null;
            else {
                return o.setCoefficient(r);
            }
        });
    }

    static <E> void subtract(MonomialSet<Monomial<E>> polynomial, Monomial<E> term, Ring<E> ring) {
        if (ring.isZero(term.coefficient))
            return;
        Monomial<E> pTerm = polynomial.get(term);
        if (pTerm == null)
            polynomial.add(term.negate(ring));
        else {
            E r = ring.subtract(pTerm.coefficient, term.coefficient);
            if (ring.isZero(r))
                polynomial.remove(pTerm);
            else
                polynomial.add(pTerm.setCoefficient(r));
        }
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param ring     the ring
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> create(int nVariables, Ring<E> ring, Comparator<DegreeVector> ordering, Iterable<Monomial<E>> terms) {
        MonomialSet<Monomial<E>> map = new MonomialSet<>(ordering);
        for (Monomial<E> term : terms)
            add(map, term.setRing(ring), ring);

        return new MultivariatePolynomial<>(nVariables, ring, ordering, map);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param ring     the ring
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> create(int nVariables, Ring<E> ring, Comparator<DegreeVector> ordering, Monomial<E>... terms) {
        return create(nVariables, ring, ordering, Arrays.asList(terms));
    }

    /**
     * Creates zero polynomial.
     *
     * @param nVariables number of variables
     * @param ring       the ring
     * @param ordering   the ordering
     * @return zero polynomial
     */
    public static <E> MultivariatePolynomial<E> zero(int nVariables, Ring<E> ring, Comparator<DegreeVector> ordering) {
        return new MultivariatePolynomial<>(nVariables, ring, ordering, new MonomialSet<>(ordering));
    }

    /**
     * Creates unit polynomial.
     *
     * @param nVariables number of variables
     * @param ring       the ring
     * @param ordering   the ordering
     * @return unit polynomial
     */
    public static <E> MultivariatePolynomial<E> one(int nVariables, Ring<E> ring, Comparator<DegreeVector> ordering) {
        return create(nVariables, ring, ordering, Monomial.withZeroExponents(nVariables, ring.getOne()));
    }

    /**
     * Parse multivariate Z[X] polynomial from string.
     *
     * @param string    the string
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate Z[X] polynomial
     */
    public static MultivariatePolynomial<BigInteger> parse(String string, String... variables) {
        return Parser.parse(string, MonomialOrder.LEX, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    the string
     * @param ring      the ring
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Ring<E> ring, String... variables) {
        return Parser.parse(string, ring, ring, MonomialOrder.LEX, variables);
    }

    /**
     * Parse multivariate Z[X] polynomial from string.
     *
     * @param string    the string
     * @param ordering  monomial order
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return Z[X] multivariate polynomial
     */
    public static MultivariatePolynomial<BigInteger> parse(String string, Comparator<DegreeVector> ordering, String... variables) {
        return Parser.parse(string, ordering, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    the string
     * @param ring      the ring
     * @param ordering  monomial order
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are different,
     *                  since the first one is considered as Z[a], while the second as Z[a,b]
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Ring<E> ring, Comparator<DegreeVector> ordering, String... variables) {
        return Parser.parse(string, ring, ring, ordering, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string the string
     * @param ring   the ring
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> parse(String string, Ring<E> ring) {
        return Parser.parse(string, ring, ring);
    }

    /**
     * Converts multivariate polynomial over BigIntegers to multivariate polynomial over machine modular integers
     *
     * @param poly the polynomial
     * @return multivariate polynomial over machine sized modular integers
     * @throws IllegalArgumentException if poly.ring is not Zp
     * @throws ArithmeticException      if some of coefficients will not exactly fit in a {@code long}.
     */
    public static MultivariatePolynomialZp64 asOverZp64(MultivariatePolynomial<BigInteger> poly) {
        if (!(poly.ring instanceof IntegersZp))
            throw new IllegalArgumentException("Poly is not over modular ring: " + poly.ring);
        IntegersZp ring = (IntegersZp) poly.ring;
        MonomialSet<MonomialZp64> terms = new MonomialSet<>(poly.ordering);
        for (Monomial<BigInteger> term : poly.terms)
            terms.add(new MonomialZp64(term.exponents, term.totalDegree, term.coefficient.longValueExact()));
        return MultivariatePolynomialZp64.create(poly.nVariables, ring.asMachineRing(), poly.ordering, terms);
    }

    /**
     * Converts univariate polynomial to multivariate.
     *
     * @param poly       univariate polynomial
     * @param nVariables number of variables in the result
     * @param variable   variable that will be used as a primary variable
     * @param ordering   monomial order
     * @return multivariate polynomial
     */
    public static <E> MultivariatePolynomial<E> asMultivariate(UnivariatePolynomial<E> poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        MonomialSet<Monomial<E>> map = new MonomialSet<>(ordering);
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;
            map.add(new Monomial<>(degreeVector, i, poly.get(i)));
        }
        return new MultivariatePolynomial<>(nVariables, poly.ring, ordering, map);
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<E> asUnivariate() {
        if (isConstant())
            return UnivariatePolynomial.constant(ring, lc());
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
        E[] univarData = ring.createZeroesArray(degrees[theVar] + 1);
        for (Monomial<E> e : terms)
            univarData[e.exponents[theVar]] = e.coefficient;
        return UnivariatePolynomial.createUnsafe(ring, univarData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomial<E>> asOverUnivariate(int variable) {
        UnivariatePolynomial<E> factory = UnivariatePolynomial.zero(ring);
        UnivariateRing<UnivariatePolynomial<E>> pDomain = new UnivariateRing<>(factory);
        MonomialSet<Monomial<UnivariatePolynomial<E>>> newData = new MonomialSet<>(ordering);
        for (Monomial<E> e : terms) {
            add(newData, new Monomial<>(
                            e.set(variable, 0).exponents,
                            factory.createMonomial(e.coefficient, e.exponents[variable])),
                    pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomial<E>> asOverUnivariateEliminate(int variable) {
        UnivariatePolynomial<E> factory = UnivariatePolynomial.zero(ring);
        UnivariateRing<UnivariatePolynomial<E>> pDomain = new UnivariateRing<>(factory);
        MonomialSet<Monomial<UnivariatePolynomial<E>>> newData = new MonomialSet<>(ordering);
        for (Monomial<E> e : terms) {
            add(newData, new Monomial<>(
                            e.without(variable).exponents,
                            factory.createMonomial(e.coefficient, e.exponents[variable])),
                    pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomial<E>> asOverMultivariate(int... variables) {
        Ring<MultivariatePolynomial<E>> ring = new MultivariateRing<>(this);
        MonomialSet<Monomial<MultivariatePolynomial<E>>> terms = new MonomialSet<>(ordering);
        for (Monomial<E> term : this) {
            int[] coeffExponents = new int[nVariables];
            for (int var : variables)
                coeffExponents[var] = term.exponents[var];
            Monomial<E> restTerm = term.setZero(variables, this.ring.getOne());
            Monomial<MultivariatePolynomial<E>> newTerm
                    = new Monomial<>(
                    restTerm.exponents,
                    restTerm.totalDegree,
                    create(new Monomial<>(coeffExponents, term.coefficient)));

            add(terms, newTerm, ring);
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, terms);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomial<E>> asOverMultivariateEliminate(int... variables) {
        variables = variables.clone();
        Arrays.sort(variables);
        int[] restVariables = ArraysUtil.intSetDifference(ArraysUtil.sequence(nVariables), variables);
        Ring<MultivariatePolynomial<E>> ring = new MultivariateRing<>(create(variables.length, new MonomialSet<>(ordering)));
        MonomialSet<Monomial<MultivariatePolynomial<E>>> terms = new MonomialSet<>(ordering);
        for (Monomial<E> term : this) {
            int i = 0;
            int[] coeffExponents = new int[variables.length];
            for (int var : variables)
                coeffExponents[i++] = term.exponents[var];

            i = 0;
            int[] termExponents = new int[restVariables.length];
            for (int var : restVariables)
                termExponents[i++] = term.exponents[var];

            Monomial<MultivariatePolynomial<E>> newTerm
                    = new Monomial<>(
                    termExponents,
                    create(variables.length, this.ring, this.ordering, new Monomial<>(coeffExponents, term.coefficient)));

            add(terms, newTerm, ring);
        }
        return new MultivariatePolynomial<>(restVariables.length, ring, ordering, terms);
    }

    /**
     * Converts multivariate polynomial over univariate polynomial ring (R[variable][other_variables]) to a multivariate
     * polynomial over coefficient ring (R[variables])
     *
     * @param poly     the polynomial
     * @param variable the variable to insert
     * @return multivariate polynomial over normal coefficient ring
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(MultivariatePolynomial<UnivariatePolynomial<E>> poly, int variable) {
        Ring<E> ring = poly.ring.getZero().ring;
        int nVariables = poly.nVariables + 1;
        MultivariatePolynomial<E> result = zero(nVariables, ring, poly.ordering);
        for (Monomial<UnivariatePolynomial<E>> entry : poly.terms) {
            UnivariatePolynomial<E> uPoly = entry.coefficient;
            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new Monomial<>(cdv, uPoly.get(i)));
            }
        }
        return result;
    }

    /**
     * Converts multivariate polynomial over multivariate polynomial ring to a multivariate polynomial over coefficient
     * ring
     *
     * @param poly the polynomial
     * @return multivariate polynomial over normal coefficient ring
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(MultivariatePolynomial<MultivariatePolynomial<E>> poly) {
        Ring<E> ring = poly.ring.getZero().ring;
        int nVariables = poly.nVariables;
        MultivariatePolynomial<E> result = zero(nVariables, ring, poly.ordering);
        for (Monomial<MultivariatePolynomial<E>> term : poly.terms) {
            MultivariatePolynomial<E> uPoly = term.coefficient;
            result.add(uPoly.clone().multiply(new Monomial<>(term.exponents, term.totalDegree, ring.getOne())));
        }
        return result;
    }

    /**
     * Converts multivariate polynomial over multivariate polynomial ring to a multivariate polynomial over coefficient
     * ring
     *
     * @param poly the polynomial
     * @return multivariate polynomial over normal coefficient ring
     */
    public static <E> MultivariatePolynomial<E> asNormalMultivariate(
            MultivariatePolynomial<MultivariatePolynomial<E>> poly,
            int[] coefficientVariables, int[] mainVariables) {
        Ring<E> ring = poly.ring.getZero().ring;
        int nVariables = coefficientVariables.length + mainVariables.length;
        MultivariatePolynomial<E> result = zero(nVariables, ring, poly.ordering);
        for (Monomial<MultivariatePolynomial<E>> term : poly.terms) {
            MultivariatePolynomial<E> coefficient =
                    term.coefficient.joinNewVariables(nVariables, coefficientVariables);
            Monomial<MultivariatePolynomial<E>> t = term.joinNewVariables(nVariables, mainVariables);
            result.add(coefficient.multiply(new Monomial<>(t.exponents, t.totalDegree, ring.getOne())));
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
        return new MultivariatePolynomial<>(poly.nVariables, Rings.Z, poly.ordering, copy ? poly.terms.clone() : poly.terms);
    }

    /**
     * Converts Zp[x] polynomial to Z[x] polynomial formed from the coefficients of this represented in symmetric
     * modular form ({@code -modulus/2 <= cfx <= modulus/2}).
     *
     * @param poly Zp polynomial
     * @return Z[x] version of the poly with coefficients represented in symmetric modular form ({@code -modulus/2 <=
     * cfx <= modulus/2}).
     * @throws IllegalArgumentException is {@code poly.ring} is not a {@link IntegersZp}
     */
    public static MultivariatePolynomial<BigInteger> asPolyZSymmetric(MultivariatePolynomial<BigInteger> poly) {
        if (!(poly.ring instanceof IntegersZp))
            throw new IllegalArgumentException("Not a modular ring: " + poly.ring);
        IntegersZp ring = (IntegersZp) poly.ring;
        MonomialSet<Monomial<BigInteger>> newTerms = new MonomialSet<>(poly.ordering);
        for (Monomial<BigInteger> term : poly)
            newTerms.add(term.setCoefficient(ring.symmetricForm(term.coefficient)));

        return new MultivariatePolynomial<>(poly.nVariables, Rings.Z, poly.ordering, newTerms);
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
    MultivariatePolynomial<E> create(int nVariables, Comparator<DegreeVector> ordering, MonomialSet<Monomial<E>> monomialTerms) {
        return new MultivariatePolynomial<>(nVariables, ring, ordering, monomialTerms);
    }

    @Override
    public Monomial<E> createTermWithUnitCoefficient(int[] exponents) {
        return new Monomial<E>(exponents, ring.getOne());
    }

    @Override
    public Monomial<E> createUnitTerm() {
        return Monomial.withZeroExponents(nVariables, ring.getOne());
    }

    @Override
    public Monomial<E> createZeroTerm() {
        return Monomial.withZeroExponents(nVariables, ring.getZero());
    }

    @Override
    public boolean isOverField() {return ring.isField();}

    @Override
    public boolean isOverFiniteField() {return ring.isFiniteField();}

    @Override
    public boolean isOverZ() {return ring.equals(Rings.Z);}

    @Override
    public BigInteger coefficientRingCardinality() {return ring.cardinality();}

    @Override
    public BigInteger coefficientRingCharacteristic() {return ring.characteristic();}

    @Override
    public boolean isOverPerfectPower() {
        return ring.isPerfectPower();
    }

    @Override
    public BigInteger coefficientRingPerfectPowerBase() {
        return ring.perfectPowerBase();
    }

    @Override
    public BigInteger coefficientRingPerfectPowerExponent() {
        return ring.perfectPowerExponent();
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[] createArray(int length) {return new MultivariatePolynomial[length];}

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[][] createArray2d(int length) {
        return new MultivariatePolynomial[length][];
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[][] createArray2d(int length1, int length2) {
        return new MultivariatePolynomial[length1][length2];
    }

    @Override
    public boolean sameCoefficientRingWith(MultivariatePolynomial<E> oth) {
        return nVariables == oth.nVariables && ring.equals(oth.ring);
    }

    @Override
    public MultivariatePolynomial<E> setCoefficientRingFrom(MultivariatePolynomial<E> poly) {
        return setRing(poly.ring);
    }

    @Override
    boolean isZeroMonomial(Monomial<E> a) {
        return ring.isZero(a.coefficient);
    }

    @Override
    Monomial<E> multiply(Monomial<E> a, Monomial<E> b) {
        return a.multiply(b, ring.multiply(a.coefficient, b.coefficient));
    }

    @Override
    Monomial<E> divideOrNull(Monomial<E> a, Monomial<E> b) {
        E e = ring.divideOrNull(a.coefficient, b.coefficient);
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
     * Returns a copy of this with coefficient reduced to a {@code newRing}
     *
     * @param newRing the new ring
     * @return a copy of this reduced to the ring specified by {@code newRing}
     */
    public MultivariatePolynomial<E> setRing(Ring<E> newRing) {
        if (ring == newRing)
            return clone();
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        for (Monomial<E> e : terms)
            add(newData, e.setRing(newRing));
        return new MultivariatePolynomial<>(nVariables, newRing, ordering, newData);
    }

    /** internal API */
    public MultivariatePolynomial<E> setRingUnsafe(Ring<E> newRing) {
        return new MultivariatePolynomial<>(nVariables, newRing, ordering, terms);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public MultivariatePolynomial<E> createConstant(E val) {
        MonomialSet<Monomial<E>> data = new MonomialSet<>(ordering);
        if (!ring.isZero(val))
            data.add(Monomial.withZeroExponents(nVariables, val));
        return new MultivariatePolynomial<>(nVariables, ring, ordering, data);
    }

    @Override
    public MultivariatePolynomial<E> createZero() {
        return createConstant(ring.getZero());
    }

    @Override
    public MultivariatePolynomial<E> createOne() {
        return createConstant(ring.getOne());
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
        MonomialSet<Monomial<E>> data = new MonomialSet<>(ordering);
        if (!ring.isZero(cc))
            data.add(Monomial.withZeroExponents(nVariables, cc));
        if (!ring.isZero(lc)) {
            int[] lcDegreeVector = new int[nVariables];
            lcDegreeVector[variable] = 1;
            data.add(new Monomial<>(lcDegreeVector, 1, lc));
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, data);
    }

    @Override
    public boolean isMonic() {
        return ring.isOne(lc());
    }

    @Override
    public int signumOfLC() {
        return ring.signum(lc());
    }

    @Override
    public boolean isZero() {
        return terms.size() == 0;
    }

    @Override
    public boolean isOne() {
        if (size() != 1)
            return false;
        Monomial<E> lt = terms.first();
        return lt.isZeroVector() && ring.isOne(lt.coefficient);
    }

    @Override
    public boolean isUnitCC() {
        return ring.isOne(cc());
    }

    @Override
    public boolean isConstant() {
        return size() == 0 || (size() == 1 && terms.first().isZeroVector());
    }

    @Override
    public Monomial<E> lt() {
        return size() == 0 ? Monomial.withZeroExponents(nVariables, ring.getZero()) : terms.last();
    }

    /**
     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the
     * ordering.
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
        terms.add(lt().setCoefficient(ring.valueOf(val)));
        return this;
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public E cc() {
        Monomial<E> zero = Monomial.withZeroExponents(nVariables, ring.getZero());
        return terms.getOrDefault(zero, zero).coefficient;
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public E content() {
        return ring.gcd(coefficients());
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
     * Returns array of polynomial coefficients
     *
     * @return array of polynomial coefficients
     */
    public E[] coefficientsArray() {
        if (isZero())
            return ring.createZeroesArray(1);

        E[] array = ring.createArray(size());
        int i = 0;
        for (Monomial<E> term : this)
            array[i++] = term.coefficient;
        return array;
    }

    private static class It<V> implements Iterator<V> {
        final Iterator<Monomial<V>> inner;

        public It(Iterator<Monomial<V>> inner) {
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
        if (signumOfLC() < 0)
            c = ring.negate(c);
        MultivariatePolynomial<E> r = divideOrNull(c);
        assert r != null;
        return r;
    }

    @Override
    public MultivariatePolynomial<E> divideByLC(MultivariatePolynomial<E> other) {
        return divideOrNull(other.lc());
    }

    /**
     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of
     * the elements can't be exactly divided by the {@code factor}. NOTE: if {@code null} is returned, the content of
     * {@code this} is destroyed.
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor} or {@code null}
     */
    public MultivariatePolynomial<E> divideOrNull(E factor) {
        if (ring.isOne(factor))
            return this;
        if (ring.isField())
            return multiply(ring.reciprocal(factor)); // <- this is typically faster than the division
        for (Entry<DegreeVector, Monomial<E>> entry : terms.entrySet()) {
            Monomial<E> term = entry.getValue();
            E quot = ring.divideOrNull(term.coefficient, factor);
            if (quot == null)
                return null;
            entry.setValue(term.setCoefficient(quot));
        }
        release();
        return this;
    }

    /**
     * Divides this polynomial by a {@code factor} or throws exception if exact division is not possible
     *
     * @param factor the factor
     * @return {@code this} divided by the {@code factor}
     * @throws ArithmeticException if exact division is not possible
     */
    public MultivariatePolynomial<E> divideExact(E factor) {
        MultivariatePolynomial<E> r = divideOrNull(factor);
        if (r == null)
            throw new ArithmeticException("not divisible " + this + " / " + factor);
        return r;
    }

    @Override
    public MultivariatePolynomial<E> divideOrNull(Monomial<E> monomial) {
        if (monomial.isZeroVector())
            return divideOrNull(monomial.coefficient);
        MonomialSet<Monomial<E>> map = new MonomialSet<>(ordering);
        for (Monomial<E> term : terms) {
            Monomial<E> dv = divideOrNull(term, monomial);
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
     * @return monic this or null if the ring does not support exact division by lc
     */
    @Override
    public MultivariatePolynomial<E> monic() {
        if (isZero())
            return this;
        return divideOrNull(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is {@code
     * monic(modulus).multiply(factor)} ).
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
     * Returns a copy of this with {@code value} substituted for {@code variable}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, Object)})
     * @see #eliminate(int, Object)
     */
    public MultivariatePolynomial<E> evaluate(int variable, E value) {
        value = ring.valueOf(value);
        if (ring.isZero(value))
            return evaluateAtZero(variable);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, ring);
        return evaluate(variable, powers);
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}.
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, long)})
     */
    MultivariatePolynomial<E> evaluate(int variable, PrecomputedPowers<E> powers) {
        if (degree(variable) == 0)
            return clone();
        if (ring.isZero(powers.value))
            return evaluateAtZero(variable);
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        for (Monomial<E> el : terms) {
            E val = ring.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable, val));
        }
        return new MultivariatePolynomial<E>(nVariables, ring, ordering, newData);
    }

    UnivariatePolynomial<E> evaluateAtZeroAllExcept(int variable) {
        E[] uData = ring.createArray(degree(variable) + 1);
        out:
        for (Monomial<E> el : terms) {
            if (el.totalDegree != 0 && el.exponents[variable] == 0)
                continue;
            for (int i = 0; i < nVariables; ++i)
                if (i != variable && el.exponents[i] != 0)
                    continue out;
            int uExp = el.exponents[variable];
            uData[uExp] = ring.add(uData[uExp], el.coefficient);
        }
        return UnivariatePolynomial.createUnsafe(ring, uData);
    }

    /**
     * Returns a copy of this with {@code values} substituted for {@code variables}.
     *
     * @param variables the variables
     * @param values    the values
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, Object)})
     * @see #eliminate(int, Object)
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> evaluate(int[] variables, E[] values) {
        for (E value : values)
            if (!ring.isZero(value))
                return evaluate(new PrecomputedPowersHolder<>(nVariables, variables, values, ring), variables, ones(nVariables));

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
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        for (Monomial<E> el : terms) {
            Monomial<E> r = el;
            E value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = ring.multiply(value, powers.pow(variables[i], raiseFactors[i] * el.exponents[variables[i]]));
                r = r.setZero(variables[i], value);
            }

            add(newData, r);
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, newData);
    }

    /**
     * Evaluates this polynomial at specified points
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E>[] evaluate(int variable, E... values) {
        return Arrays.stream(values).map(p -> evaluate(variable, p)).toArray(MultivariatePolynomial[]::new);
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, long)})
     * @see #eliminate(int, long)
     */
    public MultivariatePolynomial<E> evaluate(int variable, long value) {
        return evaluate(variable, ring.valueOf(value));
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so that
     * the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables
     * = nVariables - 1})
     * @see #evaluate(int, Object)
     */
    public MultivariatePolynomial<E> eliminate(int variable, E value) {
        value = ring.valueOf(value);
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        PrecomputedPowers<E> powers = new PrecomputedPowers<>(value, ring);
        for (Monomial<E> el : terms) {
            E val = ring.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.without(variable, val));
        }
        return new MultivariatePolynomial<>(nVariables - 1, ring, ordering, newData);
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so that
     * the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables
     * = nVariables - 1})
     * @see #evaluate(int, long)
     */
    public MultivariatePolynomial<E> eliminate(int variable, long value) {
        return eliminate(variable, ring.valueOf(value));
    }

    private static final int SIZE_OF_POWERS_CACHE = 32;

    /** cached powers used to save some time */
    static final class PrecomputedPowers<E> {
        private final E value;
        private final Ring<E> ring;
        // TODO: store map here
        private final E[] precomputedPowers;

        PrecomputedPowers(E value, Ring<E> ring) {
            this(SIZE_OF_POWERS_CACHE, value, ring);
        }

        PrecomputedPowers(int cacheSize, E value, Ring<E> ring) {
            this.value = ring.valueOf(value);
            this.ring = ring;
            this.precomputedPowers = ring.createArray(cacheSize);
        }

        E pow(int exponent) {
            if (exponent >= precomputedPowers.length)
                return ring.pow(value, exponent);

            if (precomputedPowers[exponent] != null)
                return precomputedPowers[exponent];

            E result = ring.getOne();
            E k2p = value;
            int rExp = 0, kExp = 1;
            for (; ; ) {
                if ((exponent & 1) != 0)
                    precomputedPowers[rExp += kExp] = result = ring.multiply(result, k2p);
                exponent = exponent >> 1;
                if (exponent == 0)
                    return precomputedPowers[rExp] = result;
                precomputedPowers[kExp *= 2] = k2p = ring.multiply(k2p, k2p);
            }
        }
    }

    /** holds an array of precomputed powers */
    static final class PrecomputedPowersHolder<E> {
        final int cacheSize;
        final Ring<E> ring;
        final PrecomputedPowers<E>[] powers;

        PrecomputedPowersHolder(int nVariables, int[] variables, E[] points, Ring<E> ring) {
            this(SIZE_OF_POWERS_CACHE, nVariables, variables, points, ring);
        }

        @SuppressWarnings("unchecked")
        PrecomputedPowersHolder(int cacheSize, int nVariables, int[] variables, E[] points, Ring<E> ring) {
            this.cacheSize = cacheSize;
            this.ring = ring;
            this.powers = new PrecomputedPowers[nVariables];
            for (int i = 0; i < points.length; i++)
                powers[variables[i]] = points[i] == null ? null : new PrecomputedPowers<>(cacheSize, points[i], ring);
        }

        void set(int i, E point) {
            if (powers[i] == null || !powers[i].value.equals(point))
                powers[i] = new PrecomputedPowers<>(cacheSize, point, ring);
        }

        E pow(int variable, int exponent) {
            return powers[variable].pow(exponent);
        }
    }


    /**
     * Returns a copy of this with {@code poly} substituted for {@code variable}.
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
        for (Monomial<E> term : this) {
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
    public MultivariatePolynomial<E> shift(int variable, long shift) {
        return shift(variable, ring.valueOf(shift));
    }

    /**
     * Returns a copy of this with {@code variable -> variable + shift}
     *
     * @param variable the variable
     * @param shift    shift amount
     * @return a copy of this with  {@code variable -> variable + shift}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> shift(int variable, E shift) {
        if (ring.isZero(shift))
            return clone();
        shift = ring.valueOf(shift);

        USubstitution<E> shifts = new USubstitution<>(UnivariatePolynomial.createUnsafe(ring, ring.createArray(shift, ring.getOne())), variable, nVariables, ordering);
        MultivariatePolynomial<E> result = createZero();
        for (Monomial<E> term : this) {
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
     * Returns a copy of this with {@code variables -> variables + shifts}
     *
     * @param variables the variables
     * @param shifts    the corresponding shifts
     * @return a copy of this with {@code variables -> variables + shifts}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> shift(int[] variables, E[] shifts) {
        PrecomputedSubstitution<E>[] powers = new PrecomputedSubstitution[nVariables];
        boolean allShiftsAreZero = true;
        for (int i = 0; i < variables.length; ++i) {
            if (!ring.isZero(shifts[i]))
                allShiftsAreZero = false;
            powers[variables[i]] = new USubstitution<>(UnivariatePolynomial.create(ring, ring.createArray(shifts[i], ring.getOne())), variables[i], nVariables, ordering);
        }

        if (allShiftsAreZero)
            return clone();

        PrecomputedSubstitutions<E> calculatedShifts = new PrecomputedSubstitutions<>(powers);

        MultivariatePolynomial<E> result = createZero();
        for (Monomial<E> term : this) {
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
            UnivariatePolynomial<E> r = PolynomialMethods.polyPow(base, exponent, true, uCache);
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
            return PolynomialMethods.polyPow(base, exponent, true, cache);
        }
    }

    @Override
    public MultivariatePolynomial<E> negate() {
        for (Entry<DegreeVector, Monomial<E>> entry : terms.entrySet()) {
            Monomial<E> term = entry.getValue();
            entry.setValue(term.negate(ring));
        }
        release();
        return this;
    }

    @Override
    void add(MonomialSet<Monomial<E>> terms, Monomial<E> term) {
        add(terms, term, ring);
    }

    @Override
    void subtract(MonomialSet<Monomial<E>> terms, Monomial<E> term) {
        subtract(terms, term, ring);
    }

    /**
     * Adds {@code oth} to this polynomial
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomial<E> add(E oth) {
        oth = ring.valueOf(oth);
        if (ring.isZero(oth))
            return this;
        add(terms, Monomial.withZeroExponents(nVariables, oth));
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public MultivariatePolynomial<E> subtract(E oth) {
        return add(ring.negate(ring.valueOf(oth)));
    }

    @Override
    public MultivariatePolynomial<E> increment() {
        return add(ring.getOne());
    }

    @Override
    public MultivariatePolynomial<E> decrement() {
        return subtract(ring.getOne());
    }

    /**
     * Multiplies {@code this} by the {@code factor}
     *
     * @param factor the factor
     * @return {@code} this multiplied by the {@code factor}
     */
    public MultivariatePolynomial<E> multiply(E factor) {
        factor = ring.valueOf(factor);
        if (ring.isOne(factor))
            return this;
        if (ring.isZero(factor))
            return toZero();
        for (Entry<DegreeVector, Monomial<E>> entry : terms.entrySet()) {
            Monomial<E> term = entry.getValue();
            E val = ring.multiply(term.coefficient, factor);
            if (!ring.isZero(val))
                entry.setValue(term.setCoefficient(val));
        }
        release();
        return this;
    }

    @Override
    public MultivariatePolynomial<E> multiplyByLC(MultivariatePolynomial<E> other) {
        return multiply(other.lc());
    }

    @Override
    public MultivariatePolynomial<E> multiply(Monomial<E> monomial) {
        checkSameDomainWith(monomial);
        if (monomial.isZeroVector())
            return multiply(monomial.coefficient);
        if (ring.isZero(monomial.coefficient))
            return toZero();

        MonomialSet<Monomial<E>> newMap = new MonomialSet<>(ordering);
        for (Monomial<E> thisElement : terms) {
            E val = ring.multiply(thisElement.coefficient, monomial.coefficient);
            if (!ring.isZero(val))
                newMap.add(thisElement.multiply(monomial, val));
        }

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> multiply(long factor) {
        return multiply(ring.valueOf(factor));
    }

    @Override
    public MultivariatePolynomial<E> multiplyByBigInteger(BigInteger factor) {
        return multiply(ring.valueOfBigInteger(factor));
    }

    @Override
    public MultivariatePolynomial<E> multiply(MultivariatePolynomial<E> oth) {
        assertSameCoefficientRingWith(oth);
        MonomialSet<Monomial<E>> newMap = new MonomialSet<>(ordering);
        for (Monomial<E> othElement : oth.terms)
            for (Monomial<E> thisElement : terms)
                add(newMap, thisElement.multiply(othElement, ring.multiply(thisElement.coefficient, othElement.coefficient)));

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomial<E> square() {
        return multiply(this);
    }

    @Override
    public MultivariatePolynomial<E> evaluateAtRandom(int variable, RandomGenerator rnd) {
        return evaluate(variable, ring.randomElement(rnd));
    }

    @Override
    public MultivariatePolynomial<E> evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd) {
        if (degree(variable) == 0)
            return clone();
        //desired skeleton
        Set<DegreeVector> skeleton = getSkeletonExcept(variable);
        MultivariatePolynomial<E> tmp;
        do {
            E randomPoint = ring.randomElement(rnd);
            tmp = evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    @Override
    public MultivariatePolynomial<E> derivative(int variable, int order) {
        MonomialSet<Monomial<E>> newTerms = new MonomialSet<>(ordering);
        for (Monomial<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;
            E newCoefficient = term.coefficient;
            for (int i = 0; i < order; ++i)
                newCoefficient = ring.multiply(newCoefficient, ring.valueOf(exponent - i));
            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            add(newTerms, new Monomial<>(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, newTerms);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static BigInteger seriesCoefficientFactor0(int exponent, int order, IntegersZp ring) {
        if (!ring.modulus.isInt() || order < ring.modulus.intValueExact())
            return seriesCoefficientFactor1(exponent, order, ring);
        return BigInteger.valueOf(MultivariatePolynomialZp64.seriesCoefficientFactor(exponent, order, ring.asZp64()));
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static <E> E seriesCoefficientFactor1(int exponent, int order, Ring<E> ring) {
        E factor = ring.getOne();
        for (int i = 0; i < order; ++i)
            factor = ring.multiply(factor, ring.valueOf(exponent - i));
        factor = ring.divideExact(factor, ring.factorial(order));
        return factor;
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static <E> E seriesCoefficientFactor2(int exponent, int order, Ring<E> ring) {
        BigInteger factor = BigInteger.ONE;
        for (int i = 0; i < order; ++i)
            factor = factor.multiply(BigInteger.valueOf(exponent - i));
        factor = factor.divideExact(Rings.Z.factorial(order));
        return ring.valueOfBigInteger(factor);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    @SuppressWarnings("unchecked")
    static <E> E seriesCoefficientFactor(int exponent, int order, Ring<E> ring) {
        if (ring instanceof IntegersZp)
            return (E) seriesCoefficientFactor0(exponent, order, (IntegersZp) ring);
        BigInteger characteristics = ring.characteristic();
        if (characteristics == null || !characteristics.isInt() || characteristics.intValueExact() > order)
            return seriesCoefficientFactor1(exponent, order, ring);

        return seriesCoefficientFactor2(exponent, order, ring);
    }

    @Override
    public MultivariatePolynomial<E> seriesCoefficient(int variable, int order) {
        if (order == 0)
            return clone();
        if (isConstant())
            return createZero();

        MonomialSet<Monomial<E>> newTerms = new MonomialSet<>(ordering);
        for (Monomial<E> term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;

            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            E newCoefficient = ring.multiply(term.coefficient, seriesCoefficientFactor(exponent, order, ring));
            add(newTerms, new Monomial<>(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, newTerms);
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
     * @param newRing the new ring
     * @param mapper  mapping
     * @param <T>     new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms
     */
    public <T> MultivariatePolynomial<T> mapTerms(Ring<T> newRing, Function<Monomial<E>, Monomial<T>> mapper) {
        return terms.values()
                .stream()
                .map(mapper)
                .collect(new PolynomialCollector<>(() -> zero(nVariables, newRing, ordering)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newRing the new ring
     * @param mapper  mapping
     * @param <T>     new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public <T> MultivariatePolynomial<T> mapCoefficients(Ring<T> newRing, Function<E, T> mapper) {
        return mapTerms(newRing, t -> new Monomial<>(t.exponents, t.totalDegree, mapper.apply(t.coefficient)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newDomain the new ring
     * @param mapper    mapping
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public MultivariatePolynomialZp64 mapCoefficients(IntegersZp64 newDomain, ToLongFunction<E> mapper) {
        return terms.values()
                .stream()
                .map(t -> new MonomialZp64(t.exponents, t.totalDegree, mapper.applyAsLong(t.coefficient)))
                .collect(new PolynomialCollector<>(() -> MultivariatePolynomialZp64.zero(nVariables, newDomain, ordering)));
    }

    @Override
    public int compareTo(MultivariatePolynomial<E> oth) {
        int c = Integer.compare(size(), oth.size());
        if (c != 0)
            return c;
        Iterator<Monomial<E>>
                thisIt = iterator(),
                othIt = oth.iterator();

        while (thisIt.hasNext() && othIt.hasNext()) {
            Monomial<E>
                    a = thisIt.next(),
                    b = othIt.next();

            if ((c = ordering.compare(a, b)) != 0)
                return c;

            if ((c = ring.compare(a.coefficient, b.coefficient)) != 0)
                return c;
        }
        return 0;
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomial<E> clone() {
        return new MultivariatePolynomial<>(nVariables, ring, ordering, terms.clone());
    }

    @Override
    public MultivariatePolynomial<E> parsePoly(String string) {
        MultivariatePolynomial<E> r = parse(string, ring, ordering);
        if (r.nVariables != nVariables)
            return parsePoly(string, WithVariables.defaultVars(nVariables));
        return r;
    }

    @Override
    public MultivariatePolynomial<E> parsePoly(String string, String[] variables) {
        MultivariatePolynomial<E> r = parse(string, ring, ordering, variables);
        if (r.nVariables != nVariables)
            throw new IllegalArgumentException("not from this field");
        return r;
    }

    public MultivariatePolynomial<E> parsePoly(String string, ElementParser<E> eParser, String[] variables) {
        MultivariatePolynomial<E> r = Parser.parse(string, ring, eParser, ordering, variables);
        if (r.nVariables != nVariables)
            throw new IllegalArgumentException("not from this field");
        return r;
    }

    private static <E> String coeffToString(ToStringSupport<E> cfxToString, E coeff) {
        String cfs = cfxToString.toString(coeff);
        if (cfs.matches("\\d+"))
            return cfs;
        if (cfs.startsWith("(") && cfs.endsWith(")"))
            return cfs;
        else
            return "(" + cfs + ")";
    }

    @Override
    public String toString(String[] vars) {
        return toString(ring, vars);
    }

    public String toString(ToStringSupport<E> cfxToString, String[] vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (Monomial<E> term : terms) {
            E coeff = term.coefficient;
            if (ring.isZero(coeff))
                continue;
            String monomialString = term.toString(vars);
            if (first) {
                if (!ring.isOne(coeff) || monomialString.isEmpty()) {
                    sb.append(coeffToString(cfxToString, coeff));
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
                first = false;
            } else {
                if (ring.signum(coeff) > 0)
                    sb.append("+");
                else {
                    sb.append("-");
                    coeff = ring.negate(coeff);
                }

                if (!ring.isOne(coeff) || monomialString.isEmpty()) {
                    sb.append(coeffToString(cfxToString, coeff));
                    if (!monomialString.isEmpty())
                        sb.append("*");
                }
                sb.append(monomialString);
            }
        }
        return sb.length() == 0 ? "0" : sb.toString();
    }

    @Override
    public String coefficientRingToString() {
        return ring.toString();
    }
}
