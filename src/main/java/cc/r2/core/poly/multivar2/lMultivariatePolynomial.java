//package cc.r2.core.poly.multivar2;
//
//import cc.r2.core.poly.Domain;
//import cc.r2.core.poly.LongArithmetics;
//import cc.r2.core.util.ArraysUtil;
//import cc.redberry.libdivide4j.FastDivision;
//import cc.redberry.libdivide4j.FastDivision.Magic;
//
//import java.util.Arrays;
//import java.util.Comparator;
//import java.util.Map.Entry;
//import java.util.Set;
//import java.util.TreeSet;
//import java.util.regex.Pattern;
//import java.util.stream.Collectors;
//
//import static cc.redberry.libdivide4j.FastDivision.modUnsignedFast;
//import static cc.redberry.libdivide4j.FastDivision.multiplyMod128Unsigned;
//
//long;
//
///**
// * @author Stanislav Poslavsky
// * @since 1.0
// */
//public class lMultivariatePolynomial extends AbstractMultivariatepolynomial<lMonomialTerm, lMultivariatePolynomial> {
//    /** the modulus */
//    final long modulus;
//    /** magic **/
//    private final Magic magic, magic32MulMod;
//    /** whether modulus less then 2^32 (if so, faster mulmod available) **/
//    private final boolean modulusFits32;
//
//    private lMultivariatePolynomial(int nVariables, long modulus, Magic magic, Magic magic32MulMod, boolean modulusFits32, Comparator<DegreeVector> ordering, MonomialsSet<lMonomialTerm> data) {
//        super(nVariables, ordering, data);
//        this.modulus = modulus;
//        this.magic = magic;
//        this.magic32MulMod = magic32MulMod;
//        this.modulusFits32 = modulusFits32;
//    }
//
//    private lMultivariatePolynomial(int nVariables, long modulus, Comparator<DegreeVector> ordering, MonomialsSet<lMonomialTerm> data) {
//        super(nVariables, ordering, data);
//        this.modulus = modulus;
//        this.magic = FastDivision.magicUnsigned(modulus);
//        this.magic32MulMod = FastDivision.magic32ForMultiplyMod(modulus);
//        this.modulusFits32 = LongArithmetics.fits32bitWord(modulus);
//    }
//
//    /* =========================== Factory methods =========================== */
//
//    @Override
//    lMultivariatePolynomial create(MonomialsSet<lMonomialTerm> set) {
//        return new lMultivariatePolynomial(nVariables, modulus, magic, magic32MulMod, modulusFits32, ordering, set);
//    }
//
//    /**
//     * Creates multivariate polynomial from a list of monomial terms
//     *
//     * @param modulus  the modulus
//     * @param ordering term ordering
//     * @param terms    the monomial terms
//     * @return multivariate polynomial
//     */
//    public static lMultivariatePolynomial create(long modulus, Comparator<DegreeVector> ordering, lMonomialTerm... terms) {
//        if (terms.length == 0)
//            throw new IllegalArgumentException("empty");
//        lMultivariatePolynomial result = new lMultivariatePolynomial(terms[0].exponents.length, modulus, ordering, new MonomialsSet<>(ordering));
//        for (lMonomialTerm term : terms)
//            result.add(result.data, term);
//        return result;
//    }
//
//    /**
//     * Creates zero
//     *
//     * @param modulus    the modulus
//     * @param ordering   the ordering
//     * @param nVariables number of variables
//     * @return zero
//     */
//    public static lMultivariatePolynomial zero(long modulus, Comparator<DegreeVector> ordering, int nVariables) {
//        return new lMultivariatePolynomial(nVariables, modulus, ordering, new MonomialsSet<>(ordering));
//    }
//
//
//    /*=========================== Main methods ===========================*/
//
//
//    /** modulus operation */
//    long mod(long val) {
//        return modUnsignedFast(val, magic);
//    }
//
//    /** multiplyMod operation */
//    long multiplyMod(long a, long b) {
//        return modulusFits32 ? mod(a * b) : multiplyMod128Unsigned(a, b, modulus, magic32MulMod);
//    }
//
//    /** addMod operation */
//    long addMod(long a, long b) {
//        long r = a + b;
//        return r - modulus >= 0 ? r - modulus : r;
//    }
//
//    /** subtractMod operation */
//    long subtractMod(long a, long b) {
//        long r = a - b;
//        return r + ((r >> 63)&modulus);
//    }
//
//    /** negateMod operation */
//    long negateMod(long val) {
//        return val == 0 ? val : modulus - val;
//    }
//
//    /** to symmetric modulus */
//    long symmetricForm(long value) {
//        return value <= modulus / 2 ? value : value - modulus;
//    }
//
//
//    /**
//     * Switches to another modulus specified by {@code newModulus}
//     *
//     * @param newModulus the new modulus
//     * @return a copy of this with new modulus
//     */
//    @SuppressWarnings("unchecked")
//    public lMultivariatePolynomial setModulus(long newModulus) {
//        Magic magic = FastDivision.magicUnsigned(newModulus);
//        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
//        for (lMonomialTerm e : data)
//            newData.add(e.setCoefficient(FastDivision.modUnsignedFast(e.coefficient, magic)));
//        return new lMultivariatePolynomial(nVariables, newModulus, magic, FastDivision.magic32ForMultiplyMod(newModulus), LongArithmetics.fits32bitWord(newModulus), ordering, newData);
//    }
//
//    @Override
//    public lMultivariatePolynomial createOne() {
//        return createConstant(1L);
//    }
//
//    /**
//     * Creates constant polynomial with specified value
//     *
//     * @param val value
//     * @return constant polynomial with specified value
//     */
//    public lMultivariatePolynomial createConstant(long val) {
//        MonomialsSet<lMonomialTerm> data = new MonomialsSet<>(ordering);
//        if (val != 0)
//            data.add(lMonomialTerm.withZeroExponents(nVariables, val));
//        return new lMultivariatePolynomial(nVariables, modulus, magic, magic32MulMod, modulusFits32, ordering, data);
//    }
//
//    /**
//     * Creates multivariate polynomial from a list of monomial terms
//     *
//     * @param terms the monomial terms
//     * @return multivariate polynomial
//     */
//    public lMultivariatePolynomial create(lMonomialTerm... terms) {
//        if (terms.length == 0)
//            throw new IllegalArgumentException("empty");
//        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
//        for (lMonomialTerm term : terms)
//            add(map, term);
//        return new lMultivariatePolynomial(terms[0].exponents.length, modulus, magic, magic32MulMod, modulusFits32, ordering, map);
//    }
//
//    /**
//     * Creates linear polynomial of the form {@code cc + lc * variable}
//     *
//     * @param variable the variable
//     * @param cc       the constant coefficient
//     * @param lc       the leading coefficient
//     * @return linear polynomial {@code cc + lc * variable}
//     */
//    public lMultivariatePolynomial createLinear(int variable, long cc, long lc) {
//        MonomialsSet<lMonomialTerm> data = new MonomialsSet<>(ordering);
//        if (cc != 0)
//            data.add(lMonomialTerm.withZeroExponents(nVariables, cc));
//        if (lc != 0) {
//            int[] lcDegreeVector = new int[nVariables];
//            lcDegreeVector[variable] = 1;
//            lMonomialTerm lcTerm = new lMonomialTerm(lcDegreeVector, 1, lc);
//
//            data.add(lcTerm);
//        }
//        return new lMultivariatePolynomial(nVariables, modulus, magic, magic32MulMod, modulusFits32, ordering, data);
//    }
//
//
//    @Override
//    public boolean isOne() {
//        return size() == 1 && data.first().coefficient == 1L;
//    }
//
//    @Override
//    public boolean isUnitCC() {
//        return cc() == 1L;
//    }
//
//    /**
//     * Returns the leading term in this polynomial according to ordering
//     *
//     * @return the leading term in this polynomial according to ordering
//     */
//    public lMonomialTerm lt() {
//        return size() == 0 ? lMonomialTerm.withZeroExponents(nVariables, 0) : data.last();
//    }
//
//    /**
//     * Returns the leading coefficient of this polynomial that is coefficient of the largest term according to the ordering.
//     *
//     * @return leading coefficient of this polynomial
//     */
//    public long lc() {
//        return lt().coefficient;
//    }
//
//    /**
//     * Sets the leading coefficient to the specified value
//     *
//     * @param val new value for the lc
//     * @return the leading coefficient to the specified value
//     */
//    public lMultivariatePolynomial setLC(long val) {
//        val = mod(val);
//        if (isZero())
//            return add(val);
//        data.add(lt().setCoefficient(val));
//        return this;
//    }
//
//    /**
//     * Returns the constant coefficient of this polynomial.
//     *
//     * @return constant coefficient of this polynomial
//     */
//    public long cc() {
//        lMonomialTerm zero = lMonomialTerm.withZeroExponents(nVariables, 0);
//        return data.getOrDefault(zero, zero).coefficient;
//    }
//
//    /**
//     * Returns the content of this polynomial.
//     *
//     * @return content of this polynomial
//     */
//    public long content() {
//        long gcd = 1;
//
//        for (lMonomialTerm t : data) {
//            long cf = t.coefficient;
//            if (gcd == 1) {
//                gcd = cf;
//                continue;
//            } else
//                gcd = LongArithmetics.gcd(gcd, cf);
//            if (gcd == 1)
//                break;
//        }
//        return gcd;
//    }
//
//    /**
//     * Returns the monomial content of this polynomial
//     *
//     * @return the monomial content of this polynomial
//     */
//    //todo rename!!
//    public lMonomialTerm monomialContent() {
//        return commonContent(null);
//    }
//
//    /**
//     * Returns common content of {@code this} and {@code monomial}
//     *
//     * @param monomial the monomial
//     * @return common monomial factor of {@code this} and {@code monomial}
//     */
//    lMonomialTerm commonContent(lMonomialTerm monomial) {
//        int[] exponents = monomial == null ? null : monomial.exponents.clone();
//        for (lMonomialTerm degreeVector : data)
//            if (exponents == null)
//                exponents = degreeVector.exponents.clone();
//            else
//                setMin(degreeVector, exponents);
//        if (exponents == null)
//            return lMonomialTerm.withZeroExponents(nVariables, 1);
//        return new lMonomialTerm(exponents, 1);
//    }
//
//    static void setMin(lMonomialTerm degreeVector, int[] exponents) {
//        int[] dv = degreeVector.exponents;
//        for (int i = 0; i < exponents.length; ++i)
//            if (dv[i] < exponents[i])
//                exponents[i] = dv[i];
//    }
//
//    @Override
//    public lMultivariatePolynomial primitivePart() {
//        lMultivariatePolynomial r = divide(content());
//        assert r != null;
//        release();
//        return r;
//    }
//
//    @Override
//    public lMultivariatePolynomial primitivePartSameSign() {
//        return primitivePart();
//    }
//
//    /**
//     * Divides this polynomial by a {@code factor} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
//     * divided by the {@code factor}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
//     *
//     * @param factor the factor
//     * @return {@code this} divided by the {@code factor} or {@code null}
//     */
//    public lMultivariatePolynomial divide(long factor) {
//        return multiply(LongArithmetics.modInverse(factor, modulus));
//    }
//
//    /**
//     * Divides this polynomial by a {@code monomial} or returns {@code null} (causing loss of internal data) if some of the elements can't be exactly
//     * divided by the {@code monomial}. NOTE: is {@code null} is returned, the content of {@code this} is destroyed.
//     *
//     * @param monomial monomial degrees
//     * @return {@code this} divided by the {@code factor * monomial} or {@code null}
//     */
//    public lMultivariatePolynomial divideOrNull(lMonomialTerm monomial) {
//        if (monomial.isZeroVector())
//            return divide(monomial.coefficient);
//        long cfInverse = LongArithmetics.modInverse(monomial.coefficient, modulus);
//        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
//        for (lMonomialTerm term : data) {
//            lMonomialTerm dv = term.divide(monomial, multiplyMod(monomial.coefficient, cfInverse));
//            if (dv == null)
//                return null;
//            map.add(dv);
//        }
//        loadFrom(map);
//        release();
//        return this;
//    }
//
//    /**
//     * Makes this polynomial monic if possible, if not -- destroys this and returns null
//     *
//     * @return monic this or null if the domain does not support exact division by lc
//     */
//    public lMultivariatePolynomial monic() {
//        return divide(lc());
//    }
//
//    /**
//     * Substitutes {@code value} for {@code variable}.
//     *
//     * @param variable the variable
//     * @param value    the value
//     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
//     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, Object)})
//     * @see #eliminate(int, Object)
//     */
//    public lMultivariatePolynomial evaluate(int variable, long value) {
//        value = mod(value);
//        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
//        PrecomputedPowers powers = new PrecomputedPowers<>(value, domain);
//        for (lMonomialTerm el : data) {
//            add(newData, el.setZero(variable, domain.multiply(el.coefficient, powers.pow(el.exponents[variable]))));
//        }
//        return new lMultivariatePolynomial(nVariables, domain, ordering, newData);
//    }
//
//    /**
//     * Substitutes {@code values} for {@code variables}.
//     *
//     * @param variables the variables
//     * @param values    the values
//     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
//     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, Object)})
//     * @see #eliminate(int, Object)
//     */
//    @SuppressWarnings("unchecked")
//    public lMultivariatePolynomial evaluate(int[] variables, E[] values) {
//        return evaluate(new PrecomputedPowersHolder<>(values, domain), ones == null ? ones = ArraysUtil.arrayOf(1, nVariables) : ones, variables);
//    }
//
//    private int[] ones = null;
//
//    /** substitutes {@code values} for {@code variables} */
//    @SuppressWarnings("unchecked")
//    lMultivariatePolynomial evaluate(PrecomputedPowersHolder powers, int[] variables, int[] raiseFactors) {
//        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
//        for (lMonomialTerm el : data) {
//            lMonomialTerm r = el;
//            E value = el.coefficient;
//            for (int i = 0; i < variables.length; ++i) {
//                value = domain.multiply(value, powers.pow(i, raiseFactors[i] * el.exponents[variables[i]]));
//                r = r.setZero(variables[i], value);
//            }
//
//            add(newData, r);
//        }
//        return new lMultivariatePolynomial(nVariables, domain, ordering, newData);
//    }
//
//    /**
//     * Evaluates this polynomial at specified points
//     */
//    @SuppressWarnings("unchecked")
//    public lMultivariatePolynomial[] evaluate(int variable, E... values) {
//        return Arrays.stream(values).map(p -> evaluate(variable, p)).toArray(MultivariatePolynomial[]::new);
//    }
//
//    /**
//     * Substitutes {@code value} for {@code variable}. NOTE: the resulting polynomial will
//     *
//     * @param variable the variable
//     * @param value    the value
//     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} but still with the
//     * same {@link #nVariables} (though the effective number of variables is {@code nVariables - 1}, compare to {@link #eliminate(int, long)})
//     * @see #eliminate(int, long)
//     */
//    public lMultivariatePolynomial evaluate(int variable, long value) {
//        return evaluate(variable, domain.valueOf(value));
//    }
//
//    /**
//     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
//     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
//     *
//     * @param variable the variable
//     * @param value    the value
//     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
//     * @see #evaluate(int, Object)
//     */
//    public lMultivariatePolynomial eliminate(int variable, E value) {
//        value = domain.valueOf(value);
//        MonomialsSet<lMonomialTerm> newData = new MonomialsSet<>(ordering);
//        PrecomputedPowers powers = new PrecomputedPowers<>(value, domain);
//        for (lMonomialTerm el : data) {
//            add(newData, el.without(variable, domain.multiply(el.coefficient, powers.pow(el.exponents[variable]))));
//        }
//        return new lMultivariatePolynomial(nVariables - 1, domain, ordering, newData);
//    }
//
//    /**
//     * Substitutes {@code value} for {@code variable} and eliminates {@code variable} from the list of variables so
//     * that the resulting polynomial has {@code result.nVariables = this.nVariables - 1}.
//     *
//     * @param variable the variable
//     * @param value    the value
//     * @return a new multivariate polynomial with {@code value} substituted for {@code variable} and  {@code nVariables = nVariables - 1})
//     * @see #evaluate(int, long)
//     */
//    public lMultivariatePolynomial eliminate(int variable, long value) {
//        return eliminate(variable, domain.valueOf(value));
//    }
//
//    private static final int SIZE_OF_POWERS_CACHE = 32;
//
//    /** cached powers used to save some time */
//    static final class PrecomputedPowers {
//        private final long value;
//        private final Magic magic, magic32MulMod;
//        private final boolean modulusFits32;
//        private final long[] precomputedPowers;
//
//        PrecomputedPowers(long value, Magic magic) {
//            this(SIZE_OF_POWERS_CACHE, value, magic);
//        }
//
//        public PrecomputedPowers(int cacheSize, long value, Magic magic, Magic magic32MulMod, boolean modulusFits32) {
//            this.value = value;
//            this.magic = magic;
//            this.magic32MulMod = magic32MulMod;
//            this.modulusFits32 = modulusFits32;
//            this.precomputedPowers = new long[cacheSize];
//        }
//
//        long pow(int exponent) {
//            if (exponent >= SIZE_OF_POWERS_CACHE)
//                return LongArithmetics.powModUnsigned(value, exponent, magic);
//
//            if (precomputedPowers[exponent] != 0)
//                return precomputedPowers[exponent];
//
//            long result = 1;
//            long k2p = value;
//            int rExp = 0, kExp = 1;
//            for (; ; ) {
//                if ((exponent&1) != 0)
//                    precomputedPowers[rExp += kExp] = result = domain.multiply(result, k2p);
//                exponent = exponent >> 1;
//                if (exponent == 0)
//                    return precomputedPowers[rExp] = result;
//                precomputedPowers[kExp *= 2] = k2p = domain.multiply(k2p, k2p);
//            }
//        }
//    }
//
//    /** holds an array of precomputed powers */
//    static final class PrecomputedPowersHolder {
//        private final int cacheSize;
//        private final Domain domain;
//        private final PrecomputedPowers[] powers;
//
//        PrecomputedPowersHolder(E[] points, Domain domain) {
//            this(SIZE_OF_POWERS_CACHE, points, domain);
//        }
//
//        @SuppressWarnings("unchecked")
//        PrecomputedPowersHolder(int cacheSize, E[] points, Domain domain) {
//            this.cacheSize = cacheSize;
//            this.domain = domain;
//            this.powers = new PrecomputedPowers[points.length];
//            for (int i = 0; i < points.length; i++)
//                powers[i] = points[i] == null ? null : new PrecomputedPowers(cacheSize, points[i], domain);
//        }
//
//        PrecomputedPowersHolder(int cacheSize, Domain domain, PrecomputedPowers[] powers) {
//            this.cacheSize = cacheSize;
//            this.domain = domain;
//            this.powers = powers;
//        }
//
//        void set(int i, E point) {
//            if (powers[i] == null || !powers[i].value.equals(point))
//                powers[i] = new PrecomputedPowers(cacheSize, point, domain);
//        }
//
//        E pow(int i, int exponent) {
//            return powers[i].pow(exponent);
//        }
//    }
//
//    @Override
//    public lMultivariatePolynomial negate() {
//        for (Entry<DegreeVector, lMonomialTerm> entry : data.entrySet()) {
//            lMonomialTerm term = entry.getValue();
//            entry.setValue(term.setCoefficient(domain.negate(term.coefficient)));
//        }
//        release();
//        return this;
//    }
//
//    @Override
//    public lMultivariatePolynomial add(lMultivariatePolynomial oth) {
//        if (data == oth.data)
//            return multiply(2);
//        ensureCompatible(oth);
//        if (oth.isZero())
//            return this;
//        oth.data.forEach((key, value) -> add(data, value));
//        release();
//        return this;
//    }
//
//    @Override
//    public lMultivariatePolynomial subtract(lMultivariatePolynomial oth) {
//        if (data == oth.data)
//            return toZero();
//        ensureCompatible(oth);
//        oth.data.forEach((key, value) -> subtract(data, value));
//        release();
//        return this;
//    }
//
//    /**
//     * Adds {@code term} to this polynomial and returns it
//     *
//     * @param term some term
//     * @return {@code this + oth}
//     */
//    lMultivariatePolynomial add(lMonomialTerm term) {
//        ensureCompatible(term);
//        add(data, term.setDomain(domain));
//        release();
//        return this;
//    }
//
//    /**
//     * Adds terms to this polynomial and returns it
//     *
//     * @param terms terms
//     * @return {@code this + terms}
//     */
//    lMultivariatePolynomial add(lMonomialTerm... terms) {
//        if (terms.length == 0)
//            throw new IllegalArgumentException("empty");
//        for (lMonomialTerm term : terms)
//            add(term);
//        return this;
//    }
//
//    /**
//     * Adds {@code oth} to this polynomial and returns it
//     *
//     * @param oth other polynomial
//     * @return {@code this + oth}
//     */
//    public lMultivariatePolynomial add(E oth) {
//        oth = domain.valueOf(oth);
//        if (domain.isZero(oth))
//            return this;
//        add(data, MonomialTerm.withZeroExponents(nVariables, oth));
//        release();
//        return this;
//    }
//
//    /**
//     * Subtracts {@code term} to this polynomial and returns it
//     *
//     * @param term some term
//     * @return {@code this - oth}
//     */
//    lMultivariatePolynomial subtract(lMonomialTerm term) {
//        ensureCompatible(term);
//        subtract(data, term);
//        release();
//        return this;
//    }
//
//    /**
//     * Subtracts {@code oth} from this polynomial and returns it
//     *
//     * @param oth other polynomial
//     * @return {@code this - oth}
//     */
//    public lMultivariatePolynomial subtract(E oth) {
//        return add(domain.negate(domain.valueOf(oth)));
//    }
//
//    /**
//     * Removes the leading term from this polynomial
//     *
//     * @return this - this.lt()
//     */
//    public lMultivariatePolynomial subtractLt() {
//        data.pollLastEntry();
//        release();
//        return this;
//    }
//
//    @Override
//    public lMultivariatePolynomial increment() {
//        return add(domain.getOne());
//    }
//
//    @Override
//    public lMultivariatePolynomial decrement() {
//        return subtract(domain.getOne());
//    }
//
//    private void add(MonomialsSet<lMonomialTerm> map, lMonomialTerm term) {
//        if (term.coefficient == 0)
//            return;
//        map.compute(term, (thisVector, thisValue) -> {
//            if (thisValue == null)
//                return term;
//            long r = addMod(thisValue.coefficient, term.coefficient);
//            return r == 0 ? null : thisValue.setCoefficient(r);
//        });
//    }
//
//    private void subtract(MonomialsSet<lMonomialTerm> map, lMonomialTerm term) {
//        map.compute(term, (thisVector, thisValue) -> {
//            if (thisValue == null)
//                return term;
//            long r = subtractMod(thisValue.coefficient, term.coefficient);
//            return r == 0 ? null : thisValue.setCoefficient(r);
//        });
//    }
//
//    /**
//     * Raises {@code this} by the {@code factor}
//     *
//     * @param factor the factor
//     * @return {@code} this multiplied by the {@code factor}
//     */
//    public lMultivariatePolynomial multiply(E factor) {
//        factor = domain.valueOf(factor);
//        if (domain.isOne(factor))
//            return this;
//        if (domain.isZero(factor))
//            return toZero();
//        for (Entry<DegreeVector, lMonomialTerm> entry : data.entrySet()) {
//            lMonomialTerm term = entry.getValue();
//            entry.setValue(term.setCoefficient(domain.multiply(term.coefficient, factor)));
//        }
//        release();
//        return this;
//    }
//
//    private lMultivariatePolynomial loadFrom(MonomialsSet<lMonomialTerm> map) {
//        data.clear();
//        data.putAll(map);
//        release();
//        return this;
//    }
//
//    /**
//     * Multiplies {@code this} by the {@code term}
//     *
//     * @param term the factor
//     * @return {@code} this multiplied by the {@code term}
//     */
//    public lMultivariatePolynomial multiply(lMonomialTerm term) {
//        ensureCompatible(term);
//        if (term.isZeroVector())
//            return multiply(term.coefficient);
//        if (domain.isZero(term.coefficient))
//            return toZero();
//
//        MonomialsSet<lMonomialTerm> newMap = new MonomialsSet<>(ordering);
//        for (lMonomialTerm thisElement : data) {
//            E m = domain.multiply(thisElement.coefficient, term.coefficient);
//            if (!domain.isZero(m))
//                newMap.add(thisElement.multiply(term, m));
//        }
//
//        return loadFrom(newMap);
//    }
//
//    @Override
//    public lMultivariatePolynomial multiply(long factor) {
//        return multiply(domain.valueOf(factor));
//    }
//
//    @Override
//    public lMultivariatePolynomial multiply(lMultivariatePolynomial oth) {
//        ensureCompatible(oth);
//        MonomialsSet<lMonomialTerm> newMap = new MonomialsSet<>(ordering);
//        for (lMonomialTerm othElement : oth.data)
//            for (lMonomialTerm thisElement : data)
//                add(newMap, thisElement.multiply(othElement, domain.multiply(thisElement.coefficient, othElement.coefficient)));
//
//        return loadFrom(newMap);
//    }
//
//    @Override
//    public lMultivariatePolynomial square() {
//        return multiply(this);
//    }
//
//    @Override
//    @SuppressWarnings("unchecked")
//    public lMultivariatePolynomial clone() {
//        return new lMultivariatePolynomial(nVariables, domain, ordering, data.clone());
//    }
//
//    /**
//     * Returns skeleton of this poly
//     *
//     * @return skeleton of this poly
//     */
//    public Set<DegreeVector> getSkeleton() {
//        return data.keySet();
//    }
//
//    /**
//     * Returns skeleton of this poly with respect to specified {@code variables}
//     *
//     * @param variables the variables
//     * @return skeleton of this poly with respect to specified {@code variables}
//     */
//    public Set<DegreeVector> getSkeleton(int... variables) {
//        return data.keySet().stream().map(dv -> dv.of(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
//    }
//
//    /**
//     * Returns skeleton of this poly with respect to all except specified {@code variables}
//     *
//     * @param variables the variables to exclude
//     * @return skeleton of this poly with respect to all except specified {@code variables}
//     */
//    public Set<DegreeVector> getSkeletonExcept(int... variables) {
//        return data.keySet().stream().map(dv -> dv.except(variables)).collect(Collectors.toCollection(() -> new TreeSet<>(ordering)));
//    }
//
//    /**
//     * Tests whether {@code this} and {@code oth} have the same skeleton
//     *
//     * @param oth other multivariate polynomial
//     * @return {@code true} if {@code this} and {@code oth} have the same skeleton and {@code false} otherwise
//     */
//    public boolean sameSkeleton(lMultivariatePolynomial oth) {
//        return getSkeleton().equals(oth.getSkeleton());
//    }
//
//    /**
//     * Tests whether {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables}
//     *
//     * @param oth       other multivariate polynomial
//     * @param variables variables to test
//     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to specified {@code variables} and {@code false} otherwise
//     */
//    public boolean sameSkeleton(lMultivariatePolynomial oth, int... variables) {
//        return getSkeleton(variables).equals(oth.getSkeleton(variables));
//    }
//
//    /**
//     * Tests whether {@code this} and {@code oth} have the same skeleton with respect all except specified {@code variables}
//     *
//     * @param oth       other multivariate polynomial
//     * @param variables variables to exclude
//     * @return {@code true} if {@code this} and {@code oth} have the same skeleton with respect to all except specified  {@code variables} and {@code false} otherwise
//     */
//    public boolean sameSkeletonExcept(lMultivariatePolynomial oth, int... variables) {
//        return getSkeletonExcept(variables).equals(oth.getSkeletonExcept(variables));
//    }
//
//    @Override
//    public boolean equals(Object o) {
//        if (this == o) return true;
//        if (o == null || getClass() != o.getClass()) return false;
//
//        lMultivariatePolynomial that = (lMultivariatePolynomial) o;
//
//        if (nVariables != that.nVariables)
//            return false;
//        return data.equals(that.data);
//    }
//
//    @Override
//    public int hashCode() {
//        int result = nVariables;
//        result = 31 * result + data.hashCode();
//        return result;
//    }
//
//    private static final Pattern nonTrivialCoefficientString = Pattern.compile("[\\+\\-\\*]");
//
//    private static String coeffToString(E coeff) {
//        String cfs = coeff.toString();
//        if (coeff instanceof long)
//            return cfs;
//        if (nonTrivialCoefficientString.matcher(cfs).find())
//            return "(" + cfs + ")";
//        else return cfs;
//    }
//
//    public String toString(String... vars) {
//        StringBuilder sb = new StringBuilder();
//        boolean first = true;
//        for (lMonomialTerm term : data) {
//            E coeff = term.coefficient;
//            if (domain.isZero(coeff))
//                continue;
//            String monomialString = term.toString(vars);
//            if (first) {
//                if (!domain.isOne(coeff) || monomialString.isEmpty()) {
//                    sb.append(coeffToString(coeff));
//                    if (!monomialString.isEmpty())
//                        sb.append("*");
//                }
//                sb.append(monomialString);
//                first = false;
//            } else {
//                if (domain.signum(coeff) > 0)
//                    sb.append("+");
//                else {
//                    sb.append("-");
//                    coeff = domain.negate(coeff);
//                }
//
//                if (!domain.isOne(coeff) || monomialString.isEmpty()) {
//                    sb.append(coeffToString(coeff));
//                    if (!monomialString.isEmpty())
//                        sb.append("*");
//                }
//                sb.append(monomialString);
//            }
//        }
//        return sb.length() == 0 ? "0" : sb.toString();
//    }
//
//    @Override
//    public String toString() {
//        return toString(defaultVars(nVariables));
//    }
//
//    static String[] defaultVars(int nVars) {
//        char v = 'a';
//        String[] vars = new String[nVars];
//        for (int i = 0; i < nVars; i++)
//            vars[i] = Character.toString(v++);
//        return vars;
//    }
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////
////    /**
////     * Converts multivariate polynomial to univariate if it is effectively univariate
////     *
////     * @param poly the multivariate polynomial
////     * @return the univariate polynomial
////     * @throws IllegalArgumentException if {@code poly} is not actually a univariate polynomial
////     */
////    public static lMutablePolynomialZp asUnivariateZp(lMultivariatePolynomial poly) {
////        long[] data = asUnivariate(poly);
////        return lMutablePolynomialZp.createUnsafe(poly.modulus, data);
////    }
////
////    private static long[] asUnivariate(lMultivariatePolynomial poly) {
////        int[] degrees = poly.degrees();
////        int theVar = -1;
////        for (int i = 0; i < degrees.length; i++) {
////            if (degrees[i] != 0) {
////                if (theVar != -1)
////                    throw new IllegalArgumentException("not a univariate polynomial: " + poly);
////                theVar = i;
////            }
////        }
////
////        long[] data = new long[degrees[theVar] + 1];
////        for (lMonomialTerm e : poly.data)
////            data[e.exponents[theVar]] = e.coefficient;
////        return data;
////    }
////
////    /**
////     * Converts univariate polynomial to multivariate over Zp domain.
////     *
////     * @param poly       univariate polynomial
////     * @param nVariables number of variables in the result
////     * @param variable   variable that will be used as a primary variable
////     * @param ordering   ordering
////     * @return multivariate polynomial
////     */
////    public static lMultivariatePolynomial asMultivariate(lMutablePolynomialZp poly, int nVariables, int variable, Comparator<DegreeVector> ordering) {
////        return asMultivariate(poly.getDataReferenceUnsafe(), poly.degree(), nVariables, variable, poly, ordering);
////    }
////
////    private static lMultivariatePolynomial asMultivariate(long[] data, int degree, int nVariables, int variable, lMutablePolynomialZp upoly, Comparator<DegreeVector> ordering) {
////        MonomialsSet<lMonomialTerm> map = new MonomialsSet<>(ordering);
////        for (int i = 0; i <= degree; i++) {
////            if (data[i] == 0)
////                continue;
////            int[] degreeVector = new int[nVariables];
////            degreeVector[variable] = i;
////
////            map.add(new lMonomialTerm(degreeVector, i, data[i]));
////        }
////        return new lMultivariatePolynomial(nVariables, upoly.modulus, upoly.magic, upoly.magic32MulMod, upoly.modulusFits32, ordering, map);
////    }
////
////    public static MultivariatePolynomial<lMutablePolynomialZp> convertZp(lMultivariatePolynomial poly, int variable) {
////        long modulus = poly.modulus;
////        MonomialsSet<MonomialTerm<lMutablePolynomialZp>> map = new MonomialsSet<>(poly.ordering);
////        UnivariatePolynomialDomain<lMutablePolynomialZp> domain = new UnivariatePolynomialDomain<>(lMutablePolynomialZp.zero(modulus));
////        for (MonomialTerm<long> e : poly.data) {
////            MultivariatePolynomial.add(map, new MonomialTerm<>(
////                            e.without(variable).exponents,
////                            lMutablePolynomialZp.createMonomial(modulus, e.coefficient, e.exponents[variable])),
////                    domain);
////        }
////        return new MultivariatePolynomial<>(poly.nVariables - 1, domain, poly.ordering, map);
////    }
////
////    public static lMultivariatePolynomial fromZp(MultivariatePolynomial<lMutablePolynomialZp> poly, long modulus, int variable) {
////        int nVariables = poly.nVariables + 1;
////        lMultivariatePolynomial result = zero(modulus, poly.ordering, nVariables);
////        for (MonomialTerm<lMutablePolynomialZp> entry : poly.data) {
////            lMutablePolynomialZp upoly = entry.coefficient;
////            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
////            for (int i = 0; i <= upoly.degree(); ++i) {
////                if (upoly.get(i) == 0)
////                    continue;
////                int[] cdv = dv.clone();
////                cdv[variable] = i;
////                result.add(new lMonomialTerm(cdv, upoly.get(i)));
////            }
////        }
////        return result;
////    }
//}
