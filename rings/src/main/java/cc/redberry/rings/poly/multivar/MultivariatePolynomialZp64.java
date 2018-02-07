package cc.redberry.rings.poly.multivar;

import cc.redberry.libdivide4j.FastDivision;
import cc.redberry.rings.*;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZ64;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.iterator.TLongObjectIterator;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.function.Function;
import java.util.function.LongFunction;

/**
 * Multivariate polynomial over Zp ring with the modulus in the range (0, 2^62) (see {@link
 * MachineArithmetic#MAX_SUPPORTED_MODULUS}). For details on polynomial data structure and properties see {@link
 * AMultivariatePolynomial}.
 *
 * @see AMultivariatePolynomial
 * @since 1.0
 */
public final class MultivariatePolynomialZp64 extends AMultivariatePolynomial<MonomialZp64, MultivariatePolynomialZp64> {
    private static final long serialVersionUID = 1L;
    /**
     * The ring.
     */
    public final IntegersZp64 ring;

    private MultivariatePolynomialZp64(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering, MonomialSet<MonomialZp64> lMonomialTerms) {
        super(nVariables, ordering, new IMonomialAlgebra.MonomialAlgebraZp64(ring), lMonomialTerms);
        this.ring = ring;
    }

    private static void add(MonomialSet<MonomialZp64> polynomial, MonomialZp64 term, IntegersZp64 ring) {
        if (term.coefficient == 0)
            return;
        polynomial.merge(term, term, (o, n) -> {
            long r = ring.add(o.coefficient, n.coefficient);
            if (r == 0)
                return null;
            else
                return o.setCoefficient(r);
        });
    }

    private static void subtract(MonomialSet<MonomialZp64> polynomial, MonomialZp64 term, IntegersZp64 ring) {
        add(polynomial, term.setCoefficient(ring.negate(term.coefficient)), ring);
    }

    /**
     * Creates multivariate polynomial from a set of monomials
     *
     * @param nVariables number of variables
     * @param ring       the ring
     * @param ordering   term ordering
     * @param terms      the monomials
     * @return multivariate polynomial
     */
    public static MultivariatePolynomialZp64 create(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering, MonomialSet<MonomialZp64> terms) {
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, terms);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param ring     the ring
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static MultivariatePolynomialZp64 create(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering, Iterable<MonomialZp64> terms) {
        MonomialSet<MonomialZp64> map = new MonomialSet<>(ordering);
        for (MonomialZp64 term : terms)
            add(map, term.setCoefficient(ring.modulus(term.coefficient)), ring);

        return new MultivariatePolynomialZp64(nVariables, ring, ordering, map);
    }

    /**
     * Creates multivariate polynomial from a list of monomial terms
     *
     * @param ring     the ring
     * @param ordering term ordering
     * @param terms    the monomial terms
     * @return multivariate polynomial
     */
    public static MultivariatePolynomialZp64 create(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering, MonomialZp64... terms) {
        return create(nVariables, ring, ordering, Arrays.asList(terms));
    }

    /**
     * Creates zero polynomial.
     *
     * @param nVariables number of variables
     * @param ring       the ring
     * @param ordering   the ordering
     * @return zero
     */
    public static MultivariatePolynomialZp64 zero(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering) {
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, new MonomialSet<>(ordering));
    }

    /**
     * Creates unit polynomial.
     *
     * @param nVariables number of variables
     * @param ring       the ring
     * @param ordering   the ordering
     * @return unit
     */
    public static MultivariatePolynomialZp64 one(int nVariables, IntegersZp64 ring, Comparator<DegreeVector> ordering) {
        return create(nVariables, ring, ordering, new MonomialZp64(nVariables, 1L));
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
    public static MultivariatePolynomialZp64 parse(String string, IntegersZp64 ring, String... variables) {
        return parse(string, ring, MonomialOrder.LEX, variables);
    }

    /**
     * Parse multivariate polynomial from string.
     *
     * @param string    the string
     * @param ring      the ring
     * @param ordering  monomial order
     * @param variables string variables that should be taken into account. For examples: {@code parse("a", LEX)} and
     *                  {@code parse("a", LEX, "a", "b")} although give the same mathematical expressions are differ,
     *                  since the first one is considered as Z[x], while the second as Z[x1,x2]
     * @return multivariate polynomial
     */
    public static MultivariatePolynomialZp64 parse(String string, IntegersZp64 ring, Comparator<DegreeVector> ordering, String... variables) {
        IntegersZp lDomain = ring.asGenericRing();
        return MultivariatePolynomial.asOverZp64(Parser.parse(string, lDomain, lDomain, ordering, variables));
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
    public static MultivariatePolynomialZp64 asMultivariate(UnivariatePolynomialZp64 poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
        MonomialSet<MonomialZp64> map = new MonomialSet<>(ordering);
        for (int i = poly.degree(); i >= 0; --i) {
            if (poly.isZeroAt(i))
                continue;
            int[] degreeVector = new int[nVariables];
            degreeVector[variable] = i;
            map.add(new MonomialZp64(degreeVector, i, poly.get(i)));
        }
        return new MultivariatePolynomialZp64(nVariables, poly.ring, ordering, map);
    }

    @Override
    public UnivariatePolynomialZp64 asUnivariate() {
        if (isConstant())
            return UnivariatePolynomialZp64.createUnsafe(ring, new long[]{lc()});
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
        for (MonomialZp64 e : terms)
            univarData[e.exponents[theVar]] = e.coefficient;
        return UnivariatePolynomialZp64.createUnsafe(ring, univarData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomialZp64> asOverUnivariate(int variable) {
        UnivariatePolynomialZp64 factory = UnivariatePolynomialZp64.zero(ring);
        UnivariateRing<UnivariatePolynomialZp64> pDomain = new UnivariateRing<>(factory);
        MonomialSet<Monomial<UnivariatePolynomialZp64>> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 e : terms) {
            Monomial<UnivariatePolynomialZp64> eMonomial = new Monomial<>(
                    e.dvSetZero(variable),
                    factory.createMonomial(e.coefficient, e.exponents[variable]));
            MultivariatePolynomial.add(newData, eMonomial, pDomain);
        }
        return new MultivariatePolynomial<>(nVariables, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<UnivariatePolynomialZp64> asOverUnivariateEliminate(int variable) {
        UnivariatePolynomialZp64 factory = UnivariatePolynomialZp64.zero(ring);
        UnivariateRing<UnivariatePolynomialZp64> pDomain = new UnivariateRing<>(factory);
        MonomialSet<Monomial<UnivariatePolynomialZp64>> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 e : terms) {
            Monomial<UnivariatePolynomialZp64> eMonomial = new Monomial<>(
                    e.dvWithout(variable),
                    factory.createMonomial(e.coefficient, e.exponents[variable]));
            MultivariatePolynomial.add(newData, eMonomial, pDomain);
        }
        return new MultivariatePolynomial<>(nVariables - 1, pDomain, ordering, newData);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomialZp64> asOverMultivariate(int... variables) {
        Ring<MultivariatePolynomialZp64> ring = new MultivariateRing<>(this);
        MonomialSet<Monomial<MultivariatePolynomialZp64>> terms = new MonomialSet<>(ordering);
        for (MonomialZp64 term : this) {
            int[] coeffExponents = new int[nVariables];
            for (int var : variables)
                coeffExponents[var] = term.exponents[var];

            Monomial<MultivariatePolynomialZp64> newTerm =
                    new Monomial<>(
                            term.dvSetZero(variables),
                            create(new MonomialZp64(coeffExponents, term.coefficient)));

            MultivariatePolynomial.add(terms, newTerm, ring);
        }
        return new MultivariatePolynomial<>(nVariables, ring, ordering, terms);
    }

    @Override
    public MultivariatePolynomial<MultivariatePolynomialZp64> asOverMultivariateEliminate(int... variables) {
        variables = variables.clone();
        Arrays.sort(variables);
        int[] restVariables = ArraysUtil.intSetDifference(ArraysUtil.sequence(nVariables), variables);
        Ring<MultivariatePolynomialZp64> ring = new MultivariateRing<>(create(variables.length, new MonomialSet<>(ordering)));
        MonomialSet<Monomial<MultivariatePolynomialZp64>> terms = new MonomialSet<>(ordering);
        for (MonomialZp64 term : this) {
            int i = 0;
            int[] coeffExponents = new int[variables.length];
            for (int var : variables)
                coeffExponents[i++] = term.exponents[var];

            i = 0;
            int[] termExponents = new int[restVariables.length];
            for (int var : restVariables)
                termExponents[i++] = term.exponents[var];

            Monomial<MultivariatePolynomialZp64> newTerm =
                    new Monomial<>(
                            termExponents,
                            create(variables.length, this.ring, this.ordering, new MonomialZp64(coeffExponents, term.coefficient)));

            MultivariatePolynomial.add(terms, newTerm, ring);
        }
        return new MultivariatePolynomial<>(restVariables.length, ring, ordering, terms);
    }

    /**
     * Converts multivariate polynomial over univariate polynomial ring (Zp[variable][other_variables]) to a
     * multivariate polynomial over coefficient ring (Zp[all_variables])
     *
     * @param poly     the polynomial
     * @param variable the variable to insert
     * @return multivariate polynomial over normal coefficient ring
     */
    public static MultivariatePolynomialZp64 asNormalMultivariate(MultivariatePolynomial<UnivariatePolynomialZp64> poly, int variable) {
        IntegersZp64 ring = poly.ring.getZero().ring;
        int nVariables = poly.nVariables + 1;
        MultivariatePolynomialZp64 result = zero(nVariables, ring, poly.ordering);
        for (Monomial<UnivariatePolynomialZp64> entry : poly.terms) {
            UnivariatePolynomialZp64 uPoly = entry.coefficient;
            DegreeVector dv = entry.dvInsert(variable);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                result.add(new MonomialZp64(dv.dvSet(variable, i), uPoly.get(i)));
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
    public static MultivariatePolynomialZp64 asNormalMultivariate(MultivariatePolynomial<MultivariatePolynomialZp64> poly) {
        IntegersZp64 ring = poly.ring.getZero().ring;
        int nVariables = poly.nVariables;
        MultivariatePolynomialZp64 result = zero(nVariables, ring, poly.ordering);
        for (Monomial<MultivariatePolynomialZp64> term : poly.terms) {
            MultivariatePolynomialZp64 uPoly = term.coefficient;
            result.add(uPoly.clone().multiply(new MonomialZp64(term.exponents, term.totalDegree, 1)));
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
    public static MultivariatePolynomialZp64 asNormalMultivariate(
            MultivariatePolynomial<MultivariatePolynomialZp64> poly,
            int[] coefficientVariables, int[] mainVariables) {
        IntegersZp64 ring = poly.ring.getZero().ring;
        int nVariables = coefficientVariables.length + mainVariables.length;
        MultivariatePolynomialZp64 result = zero(nVariables, ring, poly.ordering);
        for (Monomial<MultivariatePolynomialZp64> term : poly.terms) {
            MultivariatePolynomialZp64 coefficient =
                    term.coefficient.joinNewVariables(nVariables, coefficientVariables);
            Monomial<MultivariatePolynomialZp64> t = term.joinNewVariables(nVariables, mainVariables);
            result.add(coefficient.multiply(new MonomialZp64(t.exponents, t.totalDegree, 1)));
        }
        return result;
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this represented in symmetric modular form ({@code
     * -modulus/2 <= cfx <= modulus/2}).
     *
     * @return Z[X] version of this with coefficients represented in symmetric modular form ({@code -modulus/2 <= cfx <=
     * modulus/2}).
     */
    public MultivariatePolynomial<BigInteger> asPolyZSymmetric() {
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (MonomialZp64 t : this)
            bTerms.add(new Monomial<>(t.exponents, t.totalDegree, BigInteger.valueOf(ring.symmetricForm(t.coefficient))));
        return new MultivariatePolynomial<>(nVariables, Rings.Z, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> asPolyZ() {
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (MonomialZp64 t : this)
            bTerms.add(t.toBigMonomial());
        return new MultivariatePolynomial<>(nVariables, Rings.Z, ordering, bTerms);
    }

    /**
     * Returns polynomial over Z formed from the coefficients of this
     *
     * @return Z[X] version of this
     */
    public MultivariatePolynomial<BigInteger> toBigPoly() {
        MonomialSet<Monomial<BigInteger>> bTerms = new MonomialSet<>(ordering);
        for (MonomialZp64 t : this)
            bTerms.add(t.toBigMonomial());
        return new MultivariatePolynomial<>(nVariables, ring.asGenericRing(), ordering, bTerms);
    }

    /* ============================================ Main methods ============================================ */

    @Override
    public MultivariatePolynomialZp64 contentAsPoly() {
        return createConstant(content());
    }

    @Override
    public MultivariatePolynomialZp64 lcAsPoly() {
        return createConstant(lc());
    }

    @Override
    public MultivariatePolynomialZp64 ccAsPoly() {
        return createConstant(cc());
    }

    @Override
    MultivariatePolynomialZp64 create(int nVariables, Comparator<DegreeVector> ordering, MonomialSet<MonomialZp64> lMonomialTerms) {
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, lMonomialTerms);
    }

    @Override
    public boolean isOverField() {return true;}

    @Override
    public boolean isOverFiniteField() {return true;}

    @Override
    public boolean isOverZ() {return false;}

    @Override
    public BigInteger coefficientRingCardinality() {return BigInteger.valueOf(ring.modulus);}

    @Override
    public BigInteger coefficientRingCharacteristic() {return BigInteger.valueOf(ring.modulus);}

    @Override
    public boolean isOverPerfectPower() {
        return ring.isPerfectPower();
    }

    @Override
    public BigInteger coefficientRingPerfectPowerBase() {
        return BigInteger.valueOf(ring.perfectPowerBase());
    }

    @Override
    public BigInteger coefficientRingPerfectPowerExponent() {
        return BigInteger.valueOf(ring.perfectPowerExponent());
    }

    @Override
    @SuppressWarnings("unchecked")
    public MultivariatePolynomialZp64[] createArray(int length) {return new MultivariatePolynomialZp64[length];}

    @Override
    public MultivariatePolynomialZp64[][] createArray2d(int length) {
        return new MultivariatePolynomialZp64[length][];
    }

    @Override
    public MultivariatePolynomialZp64[][] createArray2d(int length1, int length2) {
        return new MultivariatePolynomialZp64[length1][length2];
    }

    @Override
    public boolean sameCoefficientRingWith(MultivariatePolynomialZp64 oth) {
        return nVariables == oth.nVariables && ring.equals(oth.ring);
    }

    @Override
    public MultivariatePolynomialZp64 setCoefficientRingFrom(MultivariatePolynomialZp64 lMonomialTerms) {
        return clone();
    }

    /** release caches */
    @Override
    protected void release() {
        super.release();
        /* add cache in the future */
    }

    /**
     * Switches to another ring specified by {@code newModulus}
     *
     * @param newModulus the new modulus
     * @return a copy of this reduced to the ring specified by {@code newModulus}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomialZp64 setRing(long newModulus) {
        return setRing(new IntegersZp64(newModulus));
    }

    /**
     * Switches to another ring specified by {@code newDomain}
     *
     * @param newDomain the new ring
     * @return a copy of this reduced to the ring specified by {@code newDomain}
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomialZp64 setRing(IntegersZp64 newDomain) {
        MonomialSet<MonomialZp64> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 e : terms)
            add(newData, e.setCoefficient(newDomain.modulus(e.coefficient)));
        return new MultivariatePolynomialZp64(nVariables, newDomain, ordering, newData);
    }

    /**
     * Switches to another ring specified by {@code newRing}
     *
     * @param newRing the new ring
     * @return a copy of this reduced to the ring specified by {@code newRing}
     */
    @SuppressWarnings("unchecked")
    public <E> MultivariatePolynomial<E> setRing(Ring<E> newRing) {
        MonomialSet<Monomial<E>> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 e : terms)
            MultivariatePolynomial.add(newData, new Monomial<>(e, newRing.valueOf(e.coefficient)), newRing);
        return new MultivariatePolynomial(nVariables, newRing, ordering, newData);
    }

    /** internal API */
    public MultivariatePolynomialZp64 setRingUnsafe(IntegersZp64 newDomain) {
        return new MultivariatePolynomialZp64(nVariables, newDomain, ordering, terms);
    }

    /**
     * Creates constant polynomial with specified value
     *
     * @param val value
     * @return constant polynomial with specified value
     */
    public MultivariatePolynomialZp64 createConstant(long val) {
        MonomialSet<MonomialZp64> data = new MonomialSet<>(ordering);
        if (val != 0)
            data.add(new MonomialZp64(nVariables, val));
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, data);
    }

    @Override
    public MultivariatePolynomialZp64 createZero() {
        return createConstant(0L);
    }

    @Override
    public MultivariatePolynomialZp64 createOne() {
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
    public MultivariatePolynomialZp64 createLinear(int variable, long cc, long lc) {
        MonomialSet<MonomialZp64> data = new MonomialSet<>(ordering);

        lc = ring.modulus(lc);
        cc = ring.modulus(cc);

        if (cc != 0L)
            data.add(new MonomialZp64(nVariables, cc));
        if (lc != 0L) {
            int[] lcDegreeVector = new int[nVariables];
            lcDegreeVector[variable] = 1;
            data.add(new MonomialZp64(lcDegreeVector, 1, lc));
        }
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, data);
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
        MonomialZp64 lt = terms.first();
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

    /**
     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the
     * ordering.
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
    public MultivariatePolynomialZp64 setLC(long val) {
        if (isZero())
            return add(val);
        terms.add(lt().setCoefficient(ring.modulus(val)));
        return this;
    }

    /**
     * Returns the constant coefficient of this polynomial.
     *
     * @return constant coefficient of this polynomial
     */
    public long cc() {
        MonomialZp64 zero = new MonomialZp64(nVariables, 0L);
        return terms.getOrDefault(zero, zero).coefficient;
    }

    /**
     * Returns the content of this polynomial.
     *
     * @return content of this polynomial
     */
    public long content() {
        long gcd = -1;
        for (MonomialZp64 term : terms) {
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
    public MultivariatePolynomialZp64 primitivePart(int variable) {
        return asNormalMultivariate(asOverUnivariateEliminate(variable).primitivePart(), variable);
    }

    @Override
    public UnivariatePolynomialZp64 contentUnivariate(int variable) {
        return asOverUnivariate(variable).content();
    }

    @Override
    public MultivariatePolynomialZp64 primitivePart() {
        return divide(content());
    }

    @Override
    public MultivariatePolynomialZp64 primitivePartSameSign() {
        return primitivePart();
    }

    @Override
    public MultivariatePolynomialZp64 divideByLC(MultivariatePolynomialZp64 other) {
        return divide(other.lc());
    }

    /**
     * Divides this polynomial by a {@code factor}
     *
     * @param factor the factor
     * @return {@code this / factor}
     */
    public MultivariatePolynomialZp64 divide(long factor) {
        if (factor == 1)
            return this;
        return multiply(ring.reciprocal(factor)); // <- this is typically faster than the division
    }

    @Override
    public MultivariatePolynomialZp64 divideOrNull(MonomialZp64 monomial) {
        if (monomial.isZeroVector())
            return divide(monomial.coefficient);
        long reciprocal = ring.reciprocal(monomial.coefficient);
        MonomialSet<MonomialZp64> map = new MonomialSet<>(ordering);
        for (MonomialZp64 term : terms) {
            DegreeVector dv = term.dvDivideOrNull(monomial);
            if (dv == null)
                return null;
            map.add(new MonomialZp64(dv, ring.multiply(term.coefficient, reciprocal)));
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
    public MultivariatePolynomialZp64 monic() {
        if (isMonic())
            return this;
        if (isZero())
            return this;
        return divide(lc());
    }

    /**
     * Sets {@code this} to its monic part multiplied by the {@code factor} modulo {@code modulus} (that is {@code
     * monic(modulus).multiply(factor)} ).
     *
     * @param factor the factor
     * @return {@code this}
     */
    public MultivariatePolynomialZp64 monic(long factor) {
        return multiply(ring.multiply(ring.modulus(factor), ring.reciprocal(lc())));
    }

    @Override
    public MultivariatePolynomialZp64 monicWithLC(MultivariatePolynomialZp64 other) {
        if (lc() == other.lc())
            return this;
        return monic(other.lc());
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, long)})
     */
    public MultivariatePolynomialZp64 evaluate(int variable, long value) {
        value = ring.modulus(value);
        if (value == 0)
            return evaluateAtZero(variable);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, ring);
        return evaluate(variable, powers);
    }

    /**
     * Returns a copy of this with {@code value} substituted for {@code variable}
     *
     * @param variable the variable
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, long)})
     */
    MultivariatePolynomialZp64 evaluate(int variable, lPrecomputedPowers powers) {
        if (degree(variable) == 0)
            return clone();
        if (powers.value == 0)
            return evaluateAtZero(variable);
        MonomialSet<MonomialZp64> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 el : terms) {
            long val = ring.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.setZero(variable).setCoefficient(val));
        }
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, newData);
    }

    UnivariatePolynomialZp64 evaluateAtZeroAllExcept(int variable) {
        long[] uData = new long[degree(variable) + 1];
        for (MonomialZp64 el : terms) {
            if (el.totalDegree != el.exponents[variable])
                continue;
            int uExp = el.exponents[variable];
            uData[uExp] = ring.add(uData[uExp], el.coefficient);
        }
        return UnivariatePolynomialZp64.createUnsafe(ring, uData);
    }

    /**
     * Returns a copy of this with {@code values} substituted for {@code variables}
     *
     * @param variables the variables
     * @param values    the values
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the same
     * {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link
     * #eliminate(int, long)})
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomialZp64 evaluate(int[] variables, long[] values) {
        for (long value : values)
            if (value != 0)
                return evaluate(new lPrecomputedPowersHolder(nVariables, variables, values, ring), variables, ones(nVariables));

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
    MultivariatePolynomialZp64 evaluate(lPrecomputedPowersHolder powers, int[] variables) {
        return evaluate(powers, variables, ones(variables.length));
    }

    /** substitutes {@code values} for {@code variables} */
    @SuppressWarnings("unchecked")
    MultivariatePolynomialZp64 evaluate(lPrecomputedPowersHolder powers, int[] variables, int[] raiseFactors) {
        MonomialSet<MonomialZp64> newData = new MonomialSet<>(ordering);
        for (MonomialZp64 el : terms) {
            MonomialZp64 r = el;
            long value = el.coefficient;
            for (int i = 0; i < variables.length; ++i) {
                value = ring.multiply(value, powers.pow(variables[i], raiseFactors[i] * el.exponents[variables[i]]));
                r = r.setZero(variables[i]).setCoefficient(value);
            }

            add(newData, r);
        }
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, newData);
    }

    /**
     * Evaluates this polynomial at specified points
     */
    @SuppressWarnings("unchecked")
    public MultivariatePolynomialZp64[] evaluate(int variable, long... values) {
        return Arrays.stream(values).mapToObj(p -> evaluate(variable, p)).toArray(MultivariatePolynomialZp64[]::new);
    }

    /**
     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so that
     * the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
     *
     * @param variable the variable
     * @param value    the value
     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables
     * = nVariables - 1})
     */
    public MultivariatePolynomialZp64 eliminate(int variable, long value) {
        value = ring.modulus(value);
        MonomialSet<MonomialZp64> newData = new MonomialSet<>(ordering);
        lPrecomputedPowers powers = new lPrecomputedPowers(value, ring);
        for (MonomialZp64 el : terms) {
            long val = ring.multiply(el.coefficient, powers.pow(el.exponents[variable]));
            add(newData, el.without(variable).setCoefficient(val));
        }
        return new MultivariatePolynomialZp64(nVariables - 1, ring, ordering, newData);
    }

    private static final int SIZE_OF_POWERS_CACHE = 32;

    /** cached powers used to save some time */
    static final class lPrecomputedPowers {
        private final long value;
        private final IntegersZp64 ring;
        // TODO: store map here
        private final long[] precomputedPowers;

        lPrecomputedPowers(long value, IntegersZp64 ring) {
            this(SIZE_OF_POWERS_CACHE, value, ring);
        }

        lPrecomputedPowers(int cacheSize, long value, IntegersZp64 ring) {
            this.value = ring.modulus(value);
            this.ring = ring;
            this.precomputedPowers = new long[cacheSize];
            Arrays.fill(precomputedPowers, -1);
        }

        long pow(int exponent) {
            if (exponent >= precomputedPowers.length)
                return ring.powMod(value, exponent);

            if (precomputedPowers[exponent] != -1)
                return precomputedPowers[exponent];

            long result = 1L;
            long k2p = value;
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
    static final class lPrecomputedPowersHolder {
        final int cacheSize;
        final IntegersZp64 ring;
        final lPrecomputedPowers[] powers;

        lPrecomputedPowersHolder(int nVariables, int[] variables, long[] points, IntegersZp64 ring) {
            this(SIZE_OF_POWERS_CACHE, nVariables, variables, points, ring);
        }

        @SuppressWarnings("unchecked")
        lPrecomputedPowersHolder(int cacheSize, int nVariables, int[] variables, long[] points, IntegersZp64 ring) {
            this.cacheSize = cacheSize;
            this.ring = ring;
            this.powers = new lPrecomputedPowers[nVariables];
            for (int i = 0; i < variables.length; i++)
                powers[variables[i]] = new lPrecomputedPowers(cacheSize, points[i], ring);
        }

        void set(int i, long point) {
            if (powers[i] == null || powers[i].value != point)
                powers[i] = new lPrecomputedPowers(cacheSize, point, ring);
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
    public MultivariatePolynomialZp64 substitute(int variable, MultivariatePolynomialZp64 poly) {
        if (poly.isConstant())
            return evaluate(variable, poly.cc());
        lPrecomputedSubstitution subsPowers;
        if (poly.isEffectiveUnivariate())
            subsPowers = new lUSubstitution(poly.asUnivariate(), poly.univariateVariable(), nVariables, ordering);
        else
            subsPowers = new lMSubstitution(poly);

        MultivariatePolynomialZp64 result = createZero();
        for (MonomialZp64 term : this) {
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
    public MultivariatePolynomialZp64 shift(int variable, long shift) {
        if (shift == 0)
            return clone();

        lUSubstitution shifts = new lUSubstitution(UnivariatePolynomialZ64.create(shift, 1).modulus(ring), variable, nVariables, ordering);
        MultivariatePolynomialZp64 result = createZero();
        for (MonomialZp64 term : this) {
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
    public MultivariatePolynomialZp64 shift(int[] variables, long[] shifts) {
        lPrecomputedSubstitution[] powers = new lPrecomputedSubstitution[nVariables];
        boolean allShiftsAreZero = true;
        for (int i = 0; i < variables.length; ++i) {
            if (shifts[i] != 0)
                allShiftsAreZero = false;
            powers[variables[i]] = new lUSubstitution(UnivariatePolynomialZ64.create(shifts[i], 1).modulus(ring, false), variables[i], nVariables, ordering);
        }

        if (allShiftsAreZero)
            return clone();

        lPrecomputedSubstitutions calculatedShifts = new lPrecomputedSubstitutions(powers);

        MultivariatePolynomialZp64 result = createZero();
        for (MonomialZp64 term : this) {
            MultivariatePolynomialZp64 temp = createOne();
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

        MultivariatePolynomialZp64 getSubstitutionPower(int var, int exponent) {
            if (subs[var] == null)
                throw new IllegalArgumentException();

            return subs[var].pow(exponent);
        }
    }

    interface lPrecomputedSubstitution {
        MultivariatePolynomialZp64 pow(int exponent);
    }

    static final class lUSubstitution implements lPrecomputedSubstitution {
        final int variable;
        final int nVariables;
        final Comparator<DegreeVector> ordering;
        final UnivariatePolynomialZp64 base;
        final TIntObjectHashMap<UnivariatePolynomialZp64> uCache = new TIntObjectHashMap<>();
        final TIntObjectHashMap<MultivariatePolynomialZp64> mCache = new TIntObjectHashMap<>();

        lUSubstitution(UnivariatePolynomialZp64 base, int variable, int nVariables, Comparator<DegreeVector> ordering) {
            this.nVariables = nVariables;
            this.variable = variable;
            this.ordering = ordering;
            this.base = base;
        }

        @Override
        public MultivariatePolynomialZp64 pow(int exponent) {
            MultivariatePolynomialZp64 cached = mCache.get(exponent);
            if (cached != null)
                return cached.clone();
            UnivariatePolynomialZp64 r = PolynomialMethods.polyPow(base, exponent, true, uCache);
            mCache.put(exponent, cached = asMultivariate(r, nVariables, variable, ordering));
            return cached.clone();
        }
    }

    static final class lMSubstitution implements lPrecomputedSubstitution {
        final MultivariatePolynomialZp64 base;
        final TIntObjectHashMap<MultivariatePolynomialZp64> cache = new TIntObjectHashMap<>();

        lMSubstitution(MultivariatePolynomialZp64 base) {
            this.base = base;
        }

        @Override
        public MultivariatePolynomialZp64 pow(int exponent) {
            return PolynomialMethods.polyPow(base, exponent, true, cache);
        }
    }

    @Override
    void add(MonomialSet<MonomialZp64> terms, MonomialZp64 term) {
        add(terms, term, ring);
    }

    @Override
    void subtract(MonomialSet<MonomialZp64> terms, MonomialZp64 term) {
        subtract(terms, term, ring);
    }

    /**
     * Adds {@code oth} to this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this + oth}
     */
    public MultivariatePolynomialZp64 add(long oth) {
        oth = ring.modulus(oth);
        if (oth == 0)
            return this;
        add(terms, new MonomialZp64(nVariables, oth));
        release();
        return this;
    }

    /**
     * Subtracts {@code oth} from this polynomial and returns it
     *
     * @param oth other polynomial
     * @return {@code this - oth}
     */
    public MultivariatePolynomialZp64 subtract(long oth) {
        return add(ring.negate(ring.modulus(oth)));
    }

    @Override
    public MultivariatePolynomialZp64 increment() {
        return add(1L);
    }

    @Override
    public MultivariatePolynomialZp64 decrement() {
        return subtract(1L);
    }

    @Override
    public MultivariatePolynomialZp64 multiply(long factor) {
        factor = ring.modulus(factor);
        if (factor == 1)
            return this;
        if (factor == 0)
            return toZero();
        for (Entry<DegreeVector, MonomialZp64> entry : terms.entrySet()) {
            MonomialZp64 term = entry.getValue();
            long val = ring.multiply(term.coefficient, factor);
            if (val != 0)
                entry.setValue(term.setCoefficient(val));
        }
        release();
        return this;
    }

    @Override
    public MultivariatePolynomialZp64 multiplyByLC(MultivariatePolynomialZp64 other) {
        return multiply(other.lc());
    }

    @Override
    public MultivariatePolynomialZp64 multiplyByBigInteger(BigInteger factor) {
        return multiply(factor.mod(BigInteger.valueOf(ring.modulus)).longValueExact());
    }

    @Override
    public MultivariatePolynomialZp64 multiply(MonomialZp64 monomial) {
        checkSameDomainWith(monomial);
        if (monomial.isZeroVector())
            return multiply(monomial.coefficient);
        if (monomial.coefficient == 0)
            return toZero();

        MonomialSet<MonomialZp64> newMap = new MonomialSet<>(ordering);
        for (MonomialZp64 thisElement : terms) {
            MonomialZp64 mul = monomialAlgebra.multiply(thisElement, monomial);
            if (mul.coefficient != 0)
                newMap.add(mul);
        }

        return loadFrom(newMap);
    }

    @Override
    public MultivariatePolynomialZp64 multiply(MultivariatePolynomialZp64 oth) {
        assertSameCoefficientRingWith(oth);
        if (oth.isZero())
            return toZero();
        if (isZero())
            return this;
        if (oth.isConstant())
            return multiply(oth.cc());

        if (size() > KRONECKER_THRESHOLD && oth.size() > KRONECKER_THRESHOLD)
            return multiplyKronecker(oth);
        else
            return multiplyClassic(oth);
    }

    private MultivariatePolynomialZp64 multiplyClassic(MultivariatePolynomialZp64 oth) {
        MonomialSet<MonomialZp64> newMap = new MonomialSet<>(ordering);
        for (MonomialZp64 othElement : oth.terms)
            for (MonomialZp64 thisElement : terms)
                add(newMap, monomialAlgebra.multiply(thisElement, othElement), ring);

        return loadFrom(newMap);
    }

    private MultivariatePolynomialZp64 multiplyKronecker(MultivariatePolynomialZp64 oth) {
        int[] resultDegrees = new int[nVariables];
        int[] thisDegrees = degrees();
        int[] othDegrees = oth.degrees();
        for (int i = 0; i < resultDegrees.length; i++)
            resultDegrees[i] = thisDegrees[i] + othDegrees[i];

        long[] map = KroneckerMap(resultDegrees);
        if (map == null)
            return multiplyClassic(oth);

        // check that degrees fit long
        double threshold = 0.;
        for (int i = 0; i < nVariables; i++)
            threshold += 1.0 * resultDegrees[i] * map[i];
        threshold *= 2;

        if (threshold > Long.MAX_VALUE)
            return multiplyClassic(oth);

        return fromKronecker(multiplySparseUnivariate(ring, toKronecker(map), oth.toKronecker(map)), map);
    }

    /**
     * Convert to Kronecker's representation
     */
    private TLongObjectHashMap<CfHolder> toKronecker(long[] kroneckerMap) {
        TLongObjectHashMap<CfHolder> result = new TLongObjectHashMap<>(size());
        for (MonomialZp64 term : this) {
            long exponent = term.exponents[0];
            for (int i = 1; i < term.exponents.length; i++)
                exponent += term.exponents[i] * kroneckerMap[i];
            assert !result.contains(exponent);
            result.put(exponent, new CfHolder(term.coefficient));
        }
        return result;
    }

    private static TLongObjectHashMap<CfHolder> multiplySparseUnivariate(IntegersZp64 ring,
                                                                         TLongObjectHashMap<CfHolder> a,
                                                                         TLongObjectHashMap<CfHolder> b) {
        TLongObjectHashMap<CfHolder> result = new TLongObjectHashMap<>(a.size() + b.size());
        TLongObjectIterator<CfHolder> ait = a.iterator();
        while (ait.hasNext()) {
            ait.advance();
            TLongObjectIterator<CfHolder> bit = b.iterator();
            while (bit.hasNext()) {
                bit.advance();

                long deg = ait.key() + bit.key();
                long val = ring.multiply(ait.value().coefficient, bit.value().coefficient);

                CfHolder r = result.putIfAbsent(deg, new CfHolder(val));
                if (r != null)
                    r.coefficient = ring.add(r.coefficient, val);
            }
        }
        return result;
    }

    private MultivariatePolynomialZp64 fromKronecker(TLongObjectHashMap<CfHolder> p,
                                                     long[] kroneckerMap) {
        terms.clear();
        TLongObjectIterator<CfHolder> it = p.iterator();
        while (it.hasNext()) {
            it.advance();
            if (it.value().coefficient == 0)
                continue;
            long exponent = it.key();
            int[] exponents = new int[nVariables];
            for (int i = 0; i < nVariables; i++) {
                long div = exponent / kroneckerMap[nVariables - i - 1];
                exponent = exponent - (div * kroneckerMap[nVariables - i - 1]);
                exponents[nVariables - i - 1] = MachineArithmetic.safeToInt(div);
            }
            terms.add(new MonomialZp64(exponents, it.value().coefficient));
        }
        return this;
    }

    static final class CfHolder {
        long coefficient = 0;

        CfHolder(long coefficient) {
            this.coefficient = coefficient;
        }
    }

    @Override
    public MultivariatePolynomialZp64 square() {
        return multiply(this);
    }

    @Override
    public MultivariatePolynomialZp64 evaluateAtRandom(int variable, RandomGenerator rnd) {
        return evaluate(variable, ring.randomElement(rnd));
    }

    @Override
    public MultivariatePolynomialZp64 evaluateAtRandomPreservingSkeleton(int variable, RandomGenerator rnd) {
        //desired skeleton
        Set<DegreeVector> skeleton = getSkeletonExcept(variable);
        MultivariatePolynomialZp64 tmp;
        do {
            long randomPoint = ring.randomElement(rnd);
            tmp = evaluate(variable, randomPoint);
        } while (!skeleton.equals(tmp.getSkeleton()));
        return tmp;
    }

    @Override
    public MultivariatePolynomialZp64 derivative(int variable, int order) {
        if (order == 0)
            return this.clone();
        if (isConstant())
            return createZero();
        MonomialSet<MonomialZp64> newTerms = new MonomialSet<>(ordering);
        for (MonomialZp64 term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;
            long newCoefficient = term.coefficient;
            for (int i = 0; i < order; ++i)
                newCoefficient = ring.multiply(newCoefficient, exponent - i);
            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            add(newTerms, new MonomialZp64(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, newTerms);
    }

    /** exp * (exp - 1) * ... * (exp - order + 1) / (1 * 2 * ... * order) mod modulus */
    static long seriesCoefficientFactor(int exponent, int order, IntegersZp64 ring) {
        IntegersZp64 baseDomain = ring.perfectPowerBaseDomain();
        if (order < baseDomain.modulus) {
            long factor = 1;
            for (int i = 0; i < order; ++i)
                factor = ring.multiply(factor, exponent - i);
            factor = ring.divide(factor, ring.factorial(order));
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
                numerator = ring.multiply(numerator, ring == baseDomain ? numMod : ring.modulus(num));

            long
                    den = i,
                    denMod = baseDomain.modulus(i);

            while (den > 1 && denMod == 0) {
                den = FastDivision.divideSignedFast(den, baseDomain.magic);
                denMod = baseDomain.modulus(den);
                ++denZeros;
            }

            if (denMod != 0)
                denominator = ring.multiply(denominator, ring == baseDomain ? denMod : ring.modulus(den));
        }

        if (numZeros > denZeros)
            numerator = ring.multiply(numerator, ring.powMod(baseDomain.modulus, numZeros - denZeros));
        else if (denZeros < numZeros)
            denominator = ring.multiply(denominator, ring.powMod(baseDomain.modulus, denZeros - numZeros));

        if (numerator == 0)
            return numerator;
        return ring.divide(numerator, denominator);
    }

    @Override
    public MultivariatePolynomialZp64 seriesCoefficient(int variable, int order) {
        if (order == 0)
            return this.clone();
        if (isConstant())
            return createZero();

        MonomialSet<MonomialZp64> newTerms = new MonomialSet<>(ordering);
        for (MonomialZp64 term : this) {
            int exponent = term.exponents[variable];
            if (exponent < order)
                continue;

            int[] newExponents = term.exponents.clone();
            newExponents[variable] -= order;

            long newCoefficient = ring.multiply(term.coefficient, seriesCoefficientFactor(exponent, order, ring));
            add(newTerms, new MonomialZp64(newExponents, term.totalDegree - order, newCoefficient));
        }
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, newTerms);
    }

    /**
     * Maps terms of this using specified mapping function
     *
     * @param newRing the new ring
     * @param mapper  mapping
     * @param <T>     new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms
     */
    public <T> MultivariatePolynomial<T> mapTerms(Ring<T> newRing, Function<MonomialZp64, Monomial<T>> mapper) {
        return terms.values()
                .stream()
                .map(mapper)
                .collect(new PolynomialCollector<>(() -> MultivariatePolynomial.zero(nVariables, newRing, ordering)));
    }

    /**
     * Maps coefficients of this using specified mapping function
     *
     * @param newRing the new ring
     * @param mapper  mapping
     * @param <T>     new element type
     * @return a new polynomial with terms obtained by applying mapper to this terms (only coefficients are changed)
     */
    public <T> MultivariatePolynomial<T> mapCoefficients(Ring<T> newRing, LongFunction<T> mapper) {
        return mapTerms(newRing, t -> new Monomial<>(t.exponents, t.totalDegree, mapper.apply(t.coefficient)));
    }

    @Override
    public int compareTo(MultivariatePolynomialZp64 oth) {
        int c = Integer.compare(size(), oth.size());
        if (c != 0)
            return c;
        Iterator<MonomialZp64>
                thisIt = iterator(),
                othIt = oth.iterator();

        while (thisIt.hasNext() && othIt.hasNext()) {
            MonomialZp64
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
    public MultivariatePolynomialZp64 clone() {
        return new MultivariatePolynomialZp64(nVariables, ring, ordering, terms.clone());
    }

    @Override
    public MultivariatePolynomialZp64 parsePoly(String string) {
        MultivariatePolynomialZp64 r = parse(string, ring, ordering);
        if (r.nVariables != nVariables)
            return parsePoly(string, WithVariables.defaultVars(nVariables));
        return r;
    }

    @Override
    public MultivariatePolynomialZp64 parsePoly(String string, String[] variables) {
        MultivariatePolynomialZp64 r = parse(string, ring, ordering, variables);
        if (r.nVariables != nVariables)
            throw new IllegalArgumentException("not from this field");
        return r;
    }

    @Override
    public String toString(String... vars) {
        StringBuilder sb = new StringBuilder();
        boolean first = true;
        for (MonomialZp64 term : terms) {
            long coeff = term.coefficient;
            if (coeff == 0)
                continue;
            String monomialString = term.dvToString(vars);
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
                    coeff = ring.negate(coeff);
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
    public String coefficientRingToString() {
        return ring.toString();
    }
}
