package cc.r2.core.poly.multivar;


import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.number.ChineseRemainders;
import cc.r2.core.number.primes.PrimesIterator;
import cc.r2.core.poly.*;
import cc.r2.core.poly.multivar.HenselLifting.IEvaluation;
import cc.r2.core.poly.multivar.LinearAlgebra.SystemInfo;
import cc.r2.core.poly.multivar.MultivariateInterpolation.Interpolation;
import cc.r2.core.poly.multivar.MultivariateInterpolation.lInterpolation;
import cc.r2.core.poly.multivar.lMultivariatePolynomialZp.lPrecomputedPowersHolder;
import cc.r2.core.poly.univar.*;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.set.hash.TLongHashSet;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;

import static cc.r2.core.poly.multivar.MultivariatePolynomial.*;
import static cc.r2.core.poly.multivar.MultivariateReduction.divideExact;
import static cc.r2.core.poly.multivar.MultivariateReduction.dividesQ;
import static cc.r2.core.util.ArraysUtil.negate;


/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class MultivariateGCD {
    private MultivariateGCD() {}

    /**
     * Calculates greatest common divisor for array of polynomials
     *
     * @param arr set of polynomials
     * @return the gcd
     */
    public static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Poly... arr) {
        return PolynomialGCD(arr, MultivariateGCD::PolynomialGCD);
    }

    /**
     * Calculates greatest common divisor for array of polynomials
     *
     * @param arr       set of polynomials
     * @param algorithm gcd algorithm
     * @return the gcd
     */
    @SuppressWarnings("unchecked")
    static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Poly[] arr, BiFunction<Poly, Poly, Poly> algorithm) {
        //choose poly of minimal total degree
        int iMin = 0, degSum = arr[0].degreeSum();
        for (int i = 1; i < arr.length; ++i) {
            int t = arr[i].degreeSum();
            if (t < degSum) {
                iMin = i;
                degSum = t;
            }
        }

        ArraysUtil.swap(arr, 0, iMin);
        AMultivariatePolynomial
                a = arr[0],
                b = arr[1].clone();

        RandomGenerator rnd = PrivateRandom.getRandom();
        for (int i = 2; i < arr.length; i++) {
            long factor;
            do { factor = rnd.nextInt(2048); } while (factor == 0);
            b = b.add((AMultivariatePolynomial) arr[i].clone().multiply(factor));
        }

        AMultivariatePolynomial gcd = algorithm.apply((Poly) a, (Poly) b);
        if (gcd.isConstant()) {
            ArraysUtil.swap(arr, 0, iMin); // <-restore
            return (Poly) gcd;
        }

        ArrayList<AMultivariatePolynomial> remainders = new ArrayList<>();
        for (int i = 1; i < arr.length; i++) {
            if (!dividesQ(arr[i], gcd))
                remainders.add(arr[i]);
        }

        if (!remainders.isEmpty()) {
            remainders.add(gcd);
            gcd = PolynomialGCD((Poly[]) remainders.toArray((Poly[]) a.arrayNewInstance(remainders.size())), algorithm);
        }

        ArraysUtil.swap(arr, 0, iMin); // <-restore
        return (Poly) gcd;
    }

    /**
     * Calculates greatest common divisor for array of polynomials
     *
     * @param arr set of polynomials
     * @return the gcd
     */
    public static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Iterable<Poly> arr) {
        return PolynomialGCD(arr, MultivariateGCD::PolynomialGCD);
    }

    /**
     * Calculates greatest common divisor for array of polynomials
     *
     * @param arr       set of polynomials
     * @param algorithm gcd algorithm
     * @return the gcd
     */
    @SuppressWarnings("unchecked")
    static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Iterable<Poly> arr, BiFunction<Poly, Poly, Poly> algorithm) {
        Iterator<Poly> iterator = arr.iterator();
        ArrayList<Poly> list = new ArrayList<>();
        while (iterator.hasNext())
            list.add(iterator.next());
        if (list.isEmpty())
            throw new IllegalArgumentException("Empty iterable");
        Poly[] array = (Poly[]) list.get(0).arrayNewInstance(list.size());
        return PolynomialGCD(list.toArray(array), algorithm);
    }

    /**
     * Calculates greatest common divisor for two multivariate polynomials
     *
     * @param a the first poly
     * @param b the second poly
     * @return the gcd
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial> Poly PolynomialGCD(Poly a, Poly b) {
        a.checkSameDomainWith(b);
        if (a instanceof lMultivariatePolynomialZp)
            return (Poly) ZippelGCD((lMultivariatePolynomialZp) a, (lMultivariatePolynomialZp) b);
        else if (a instanceof MultivariatePolynomial) {
            Domain domain = ((MultivariatePolynomial) a).domain;
            if (Integers.Integers.equals(domain))
                return (Poly) ModularGCD((MultivariatePolynomial<BigInteger>) a, (MultivariatePolynomial<BigInteger>) b);
            else if (domain.isField())
                return (Poly) ZippelGCD((MultivariatePolynomial) a, (MultivariatePolynomial) b);
            else
                throw new RuntimeException("GCD over " + domain + " domain not supported.");
        } else
            throw new RuntimeException();
    }

    /* ============================================== Auxiliary methods ============================================= */

    /** calculates the inverse permutation */
    private static int[] inversePermutation(int[] permutation) {
        final int[] inv = new int[permutation.length];
        for (int i = permutation.length - 1; i >= 0; --i)
            inv[permutation[i]] = i;
        return inv;
    }

    /** structure with required input for GCD algorithms */
    private static final class GCDInput<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** input polynomials (with variables renamed) and earlyGCD if possible (in trivial cases) */
        final Poly aReduced, bReduced, earlyGCD;
        /** gcd degree bounds, mapping used to rename variables so that degreeBounds are in descending order */
        final int[] degreeBounds, mapping;
        /** last present variable */
        final int lastPresentVariable;
        /** domain cardinality (or -1, if cardinality is greater than Integer.MAX_VALUE) */
        final int evaluationStackLimit;
        /** GCD of monomial content of a and b **/
        final Term monomialGCD;
        /**
         * degree of irreducible univariate polynomial used to construct field extension q^n
         * (if the coefficient domain has so small cardinality so that modular algorithm will fail)
         */
        final int finiteExtensionDegree;

        GCDInput(Poly earlyGCD) {
            this.earlyGCD = earlyGCD;
            aReduced = bReduced = null;
            degreeBounds = mapping = null;
            lastPresentVariable = evaluationStackLimit = -1;
            monomialGCD = null;
            finiteExtensionDegree = -1;
        }

        GCDInput(Poly aReduced, Poly bReduced, Term monomialGCD,
                 int evaluationStackLimit, int[] degreeBounds, int[] mapping, int lastPresentVariable,
                 int finiteExtensionDegree) {
            //assert monomialGCD == null || aReduced.domain.isOne(monomialGCD.coefficient);
            this.aReduced = aReduced;
            this.bReduced = bReduced;
            this.monomialGCD = monomialGCD;
            this.earlyGCD = null;
            this.evaluationStackLimit = evaluationStackLimit;
            this.degreeBounds = degreeBounds;
            this.mapping = inversePermutation(mapping);
            this.lastPresentVariable = lastPresentVariable;
            this.finiteExtensionDegree = finiteExtensionDegree;
        }

        /** recover initial order of variables in the result */
        Poly restoreGCD(Poly result) {
            return renameVariables(result, mapping).multiply(monomialGCD);
        }
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly trivialGCD(Poly a, Poly b) {
        if (a.isZero())
            return b.clone();
        if (b.isZero())
            return a.clone();
        if (a.isConstant() || b.isConstant())
            return a.createOne();
        if (a.size() == 1)
            return gcdWithMonomial(a.lt(), b);
        if (b.size() == 1)
            return gcdWithMonomial(b.lt(), a);
        return null;
    }

    /** prepare input for modular GCD algorithms (Brown, Zippel, LinZip) */
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    GCDInput<Term, Poly> preparedGCDInput(Poly a, Poly b, BiFunction<Poly, Poly, Poly> gcdAlgorithm) {
        Poly trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return new GCDInput<>(trivialGCD);

        BigInteger domainSize = a.coefficientDomainCardinality();
        // domain cardinality, i.e. number of possible random choices
        int evaluationStackLimit = domainSize == null ? -1 : (domainSize.isInt() ? domainSize.intValue() : -1);

        // find monomial GCD
        // and remove monomial content from a and b
        a = a.clone(); b = b.clone(); // prevent rewriting original data
        Term monomialGCD = reduceMonomialContent(a, b);

        int
                nVariables = a.nVariables,
                nUsedVariables = 0, // number of really present variables in gcd
                nRedundantVariables = 0, // number of variables that present either in a or b but do not appear in the gcd
                lastGCDVariable = -1, // last variable that present in both input polynomials
                aDegrees[] = a.degrees(),
                bDegrees[] = b.degrees(),
                degreeBounds[] = new int[nVariables]; // degree bounds for gcd

        // determine initial degree bounds for gcd
        for (int i = 0; i < nVariables; i++) {
            if ((aDegrees[i] == 0 && bDegrees[i] != 0) || (aDegrees[i] != 0 && bDegrees[i] == 0))
                ++nRedundantVariables;

            degreeBounds[i] = Math.min(aDegrees[i], bDegrees[i]);
            if (degreeBounds[i] != 0) {
                ++nUsedVariables;
                lastGCDVariable = i;
            }
        }

        if (nUsedVariables == 0)
            // gcd is constant => 1
            return new GCDInput<>(a.create(monomialGCD));

        if (nUsedVariables == 1)
        // switch to univariate gcd
        {
            @SuppressWarnings("unchecked")
            IUnivariatePolynomial
                    uaContent = a.asOverUnivariate(lastGCDVariable).content(),
                    ubContent = b.asOverUnivariate(lastGCDVariable).content(),
                    iUnivar = UnivariateGCD.PolynomialGCD(uaContent, ubContent);

            Poly poly = AMultivariatePolynomial.asMultivariate(iUnivar, nVariables, lastGCDVariable, a.ordering);
            return new GCDInput<>(poly.multiply(monomialGCD));
        }

        if (nRedundantVariables > 0) {
            // some of the variables are present only in either a or b but not in both simultaneously
            // call this vars {dummies} and the rest vars as {used}
            // => then we can consider polynomials as over R[{used}][{dummies}] and calculate
            // gcd via recursive call as GCD[content(a) in R[{used}], content(b) in R[{used}]]

            int[] usedVariables = new int[nUsedVariables];
            int counter = 0;
            for (int i = 0; i < nVariables; ++i)
                if (degreeBounds[i] != 0)
                    usedVariables[counter++] = i;

            // coefficients in R[{used}][{dummies}]
            Poly[]
                    aEffective = a.asOverMultivariateEliminate(usedVariables).coefficientsArray(),
                    bEffective = b.asOverMultivariateEliminate(usedVariables).coefficientsArray();

            Poly[] all = ArraysUtil.addAll(aEffective, bEffective);
            // recursive call in PolynomialGCD
            Poly gcd = PolynomialGCD(all, gcdAlgorithm);

            gcd = gcd.joinNewVariables(nVariables, usedVariables);
            return new GCDInput<>(gcd.multiply(monomialGCD));
        }

        // now swap variables so that the first variable will have the maximal degree (univariate gcd is fast),
        // and all non-used variables are at the end of poly's

        int[] variables = ArraysUtil.sequence(nVariables);
        //sort in descending order
        ArraysUtil.quickSort(negate(degreeBounds), variables);
        negate(degreeBounds);//recover degreeBounds

        lastGCDVariable = 0; //recalculate lastPresentVariable
        for (; lastGCDVariable < degreeBounds.length; ++lastGCDVariable)
            if (degreeBounds[lastGCDVariable] == 0)
                break;
        --lastGCDVariable;

        a = renameVariables(a, variables);
        b = renameVariables(b, variables);

        // check whether coefficient domain cardinality is large enough
        int finiteExtensionDegree = 1;
        int cardinalityBound = 5 * ArraysUtil.max(degreeBounds);
        if (domainSize != null && domainSize.isInt() && domainSize.intValueExact() < cardinalityBound) {
            long ds = domainSize.intValueExact();
            finiteExtensionDegree = 2;
            long tmp = ds;
            for (; tmp < cardinalityBound; ++finiteExtensionDegree)
                tmp = tmp * ds;
        }
        return new GCDInput<>(a, b, monomialGCD, evaluationStackLimit, degreeBounds, variables, lastGCDVariable, finiteExtensionDegree);
    }

    /** gcd with monomial */
    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly gcdWithMonomial(Term monomial, Poly poly) {
        return poly.create(poly.commonContent(monomial));
    }

    /**
     * Removes monomial content from a and b and returns monomial gcd
     */
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Term reduceMonomialContent(Poly a, Poly b) {

        Term aMonomialContent = a.monomialContent();
        int[] exponentsGCD = b.monomialContent().exponents;
        AMultivariatePolynomial.setMin(aMonomialContent.exponents, exponentsGCD);
        Term monomialGCD = a.createTermWithUnitCoefficient(exponentsGCD);

        a = a.divideDegreeVectorOrNull(monomialGCD);
        b = b.divideDegreeVectorOrNull(monomialGCD);
        assert a != null && b != null;

        // remove the rest of monomial content
        a = a.divideDegreeVectorOrNull(a.monomialContent());
        b = b.divideDegreeVectorOrNull(b.monomialContent());
        assert a != null && b != null;

        return monomialGCD;
    }

    /**
     * Primitive part and content of multivariate polynomial considered as polynomial over Zp[x_i][x_1, ..., x_N]
     */
    private static final class UnivariateContent<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        final MultivariatePolynomial<uPoly> poly;
        final Poly primitivePart;
        final uPoly content;

        public UnivariateContent(MultivariatePolynomial<uPoly> poly, Poly primitivePart, uPoly content) {
            this.poly = poly;
            this.primitivePart = primitivePart;
            this.content = content;
        }
    }

    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    MultivariatePolynomial<uPoly> asOverUnivariate(Poly poly, int variable) {
        if (poly instanceof MultivariatePolynomial)
            return (MultivariatePolynomial<uPoly>) ((MultivariatePolynomial) poly).asOverUnivariateEliminate(variable);
        else if (poly instanceof lMultivariatePolynomialZp)
            return (MultivariatePolynomial<uPoly>) ((lMultivariatePolynomialZp) poly).asOverUnivariateEliminate(variable);
        else
            throw new RuntimeException();
    }

    /**
     * Returns primitive part and content of {@code poly} considered as polynomial over Zp[variable][x_1, ..., x_N]
     *
     * @param poly     the polynomial
     * @param variable the variable
     * @return primitive part and content
     */
    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    UnivariateContent<Term, Poly, uPoly> univariateContent(
            Poly poly, int variable) {
        //convert poly to Zp[var][x...]

        MultivariatePolynomial<uPoly> conv = asOverUnivariate(poly, variable);
        //univariate content
        uPoly uContent = UnivariateGCD.PolynomialGCD(conv.coefficients());
        Poly mContent = AMultivariatePolynomial.asMultivariate(uContent, poly.nVariables, variable, poly.ordering);
        Poly primitivePart = divideExact(poly, mContent);
        return new UnivariateContent(conv, primitivePart, uContent);
    }

    /** holds primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
    private static final class PrimitiveInput<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>> {
        /** primitive parts of a and b viewed as polynomials in Zp[x_k][x_1 ... x_{k-1}] */
        final Poly aPrimitive, bPrimitive;
        /** gcd of content and leading coefficient of a and b given as Zp[x_k][x_1 ... x_{k-1}] */
        final uPoly contentGCD, lcGCD;

        public PrimitiveInput(Poly aPrimitive, Poly bPrimitive,
                              uPoly contentGCD, uPoly lcGCD) {
            this.aPrimitive = aPrimitive;
            this.bPrimitive = bPrimitive;
            this.contentGCD = contentGCD;
            this.lcGCD = lcGCD;
        }
    }

    /**
     * Primitive parts, contents and leading coefficients of {@code a} and {@code b} viewed as polynomials in
     * Zp[x_k][x_1 ... x_{k-1}] with {@code x_k = variable}
     */
    @SuppressWarnings("unchecked")
    private static <
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>,
            uPoly extends IUnivariatePolynomial<uPoly>>
    PrimitiveInput<Term, Poly, uPoly> makePrimitive(Poly a, Poly b, int variable) {
        // a and b as Zp[x_k][x_1 ... x_{k-1}]
        UnivariateContent<Term, Poly, uPoly>
                aContent = univariateContent(a, variable),
                bContent = univariateContent(b, variable);

        a = aContent.primitivePart;
        b = bContent.primitivePart;

        // gcd of Zp[x_k] content and lc
        uPoly
                contentGCD = UnivariateGCD.PolynomialGCD(aContent.content, bContent.content),
                lcGCD = UnivariateGCD.PolynomialGCD(aContent.poly.lc(), bContent.poly.lc());
        return new PrimitiveInput(a, b, contentGCD, lcGCD);
    }

    /** coefficients specified variable as multivariate polynomials */
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly[] multivariateCoefficients(Poly poly, int variable) {
        return Arrays.stream(poly.degrees(variable)).mapToObj(d -> poly.coefficientOf(variable, d)).toArray(poly::arrayNewInstance);
    }

    /** gcd of content of a and b with respect to specified variable */
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly contentGCD(Poly a, Poly b, int variable, BiFunction<Poly, Poly, Poly> algorithm) {
        Poly[] aCfs = multivariateCoefficients(a, variable);
        Poly[] bCfs = multivariateCoefficients(b, variable);

        Poly[] all = a.arrayNewInstance(aCfs.length + bCfs.length);
        System.arraycopy(aCfs, 0, all, 0, aCfs.length);
        System.arraycopy(bCfs, 0, all, aCfs.length, bCfs.length);

        Poly gcd = PolynomialGCD(all, algorithm);
        assert dividesQ(a, gcd);
        assert dividesQ(b, gcd);
        return gcd;
    }

    /* =========================================== Multivariate GCD over Z ========================================== */

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("ConstantConditions")
    public static MultivariatePolynomial<BigInteger> ModularGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        CommonUtils.ensureZDomain(a, b);
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree() < b.degree())
            return ModularGCD(b, a);
        BigInteger aContent = a.content(), bContent = b.content();
        BigInteger contentGCD = BigIntegerArithmetics.gcd(aContent, bContent);
        if (a.isConstant() || b.isConstant())
            return a.createConstant(contentGCD);

        return ModularGCD0(a.clone().divideOrNull(aContent), b.clone().divideOrNull(bContent)).multiply(contentGCD);
    }

    static MultivariatePolynomial<BigInteger> ModularGCD0(
            MultivariatePolynomial<BigInteger> a,
            MultivariatePolynomial<BigInteger> b) {

        GCDInput<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>>
                gcdInput = preparedGCDInput(a, b, MultivariateGCD::ModularGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

        BigInteger lcGCD = BigIntegerArithmetics.gcd(a.lc(), b.lc());

        RandomGenerator random = PrivateRandom.getRandom();
        PrimesIterator primesLoop = new PrimesIterator(1031);
        main_loop:
        while (true) {
            // prepare the skeleton
            long basePrime = primesLoop.take();
            BigInteger bPrime = BigInteger.valueOf(basePrime);
            assert basePrime != -1 : "long overflow";

            IntegersModulo domain = new IntegersModulo(basePrime);
            // reduce Z -> Zp
            MultivariatePolynomial<BigInteger>
                    abMod = a.setDomain(domain),
                    bbMod = b.setDomain(domain);
            if (!abMod.sameSkeleton(a) || !bbMod.sameSkeleton(b))
                continue;

            lMultivariatePolynomialZp
                    aMod = asLongPolyZp(abMod),
                    bMod = asLongPolyZp(bbMod);

            // the base image
            // accumulator to update coefficients via Chineese remainding
            lMultivariatePolynomialZp base = PolynomialGCD(aMod, bMod);
            long lLcGCD = lcGCD.mod(bPrime).longValueExact();
            // scale to correct l.c.
            base = base.monic(lLcGCD);

            if (base.isConstant())
                return gcdInput.restoreGCD(a.createOne());

            // cache the previous base
            MultivariatePolynomial<BigInteger> previousBase = null;

            // over all primes
            while (true) {
                long prime = primesLoop.take();
                bPrime = BigInteger.valueOf(prime);
                domain = new IntegersModulo(bPrime);

                // reduce Z -> Zp
                abMod = a.setDomain(domain);
                bbMod = b.setDomain(domain);
                if (!abMod.sameSkeleton(a) || !bbMod.sameSkeleton(b))
                    continue;

                aMod = asLongPolyZp(abMod);
                bMod = asLongPolyZp(bbMod);

                lIntegersModulo lDomain = new lIntegersModulo(prime);

                // calculate new GCD using previously calculated skeleton via sparse interpolation
                lMultivariatePolynomialZp modularGCD = interpolateGCD(aMod, bMod, base.setDomainUnsafe(lDomain), random);
                if (modularGCD == null) {
                    // interpolation failed => assumed form is wrong => start over
                    continue main_loop;
                }

                assert dividesQ(aMod, modularGCD);
                assert dividesQ(bMod, modularGCD);

                if (modularGCD.isConstant())
                    return gcdInput.restoreGCD(a.createOne());

                // better degree bound found -> start over
                if (modularGCD.degree(0) < base.degree(0)) {
                    lLcGCD = lcGCD.mod(bPrime).longValueExact();
                    // scale to correct l.c.
                    base = modularGCD.monic(lLcGCD);
                    basePrime = prime;
                    continue;
                }

                //skip unlucky prime
                if (modularGCD.degree(0) > base.degree(0))
                    continue;

                if (LongArithmetics.isOverflowMultiply(basePrime, prime) || basePrime * prime > LongArithmetics.MAX_SUPPORTED_MODULUS)
                    break;

                //lifting
                long newBasePrime = basePrime * prime;
                long monicFactor = modularGCD.domain.multiply(
                        LongArithmetics.modInverse(modularGCD.lc(), prime),
                        lcGCD.mod(bPrime).longValueExact());

                ChineseRemainders.ChineseRemaindersMagic magic = ChineseRemainders.createMagic(basePrime, prime);
                PairIterator<lMonomialTerm, lMultivariatePolynomialZp> iterator = new PairIterator<>(base, modularGCD);
                while (iterator.hasNext()) {
                    iterator.advance();

                    lMonomialTerm
                            baseTerm = iterator.aTerm,
                            imageTerm = iterator.bTerm;

                    if (baseTerm.coefficient == 0)
                        // term is absent in the base
                        continue;

                    if (imageTerm.coefficient == 0) {
                        // term is absent in the modularGCD => remove it from the base
                        base.subtract(baseTerm);
                        continue;
                    }

                    long oth = lDomain.multiply(imageTerm.coefficient, monicFactor);

                    // update base term
                    long newCoeff = ChineseRemainders.ChineseRemainders(magic, baseTerm.coefficient, oth);
                    base.terms.add(baseTerm.setCoefficient(newCoeff));
                }

                base = base.setDomainUnsafe(new lIntegersModulo(newBasePrime));
                basePrime = newBasePrime;

                // two trials didn't change the result, probably we are done
                MultivariatePolynomial<BigInteger> candidate = base.asPolyZSymmetric().primitivePart();
                if (previousBase != null && candidate.equals(previousBase)) {
                    previousBase = candidate;
                    //first check b since b is less degree
                    if (!dividesQ(b, candidate))
                        continue;

                    if (!dividesQ(a, candidate))
                        continue;

                    return gcdInput.restoreGCD(candidate);
                }
                previousBase = candidate;
            }

            //continue lifting with multi-precision integers
            MultivariatePolynomial<BigInteger> bBase = base.toBigPoly();
            BigInteger bBasePrime = BigInteger.valueOf(basePrime);
            // over all primes
            while (true) {
                long prime = primesLoop.take();
                bPrime = BigInteger.valueOf(prime);
                domain = new IntegersModulo(bPrime);

                // reduce Z -> Zp
                abMod = a.setDomain(domain);
                bbMod = b.setDomain(domain);
                if (!abMod.sameSkeleton(a) || !bbMod.sameSkeleton(b))
                    continue;

                aMod = asLongPolyZp(abMod);
                bMod = asLongPolyZp(bbMod);

                lIntegersModulo lDomain = new lIntegersModulo(prime);

                // calculate new GCD using previously calculated skeleton via sparse interpolation
                lMultivariatePolynomialZp modularGCD = interpolateGCD(aMod, bMod, base.setDomainUnsafe(lDomain), random);
                if (modularGCD == null) {
                    // interpolation failed => assumed form is wrong => start over
                    continue main_loop;
                }

                assert dividesQ(aMod, modularGCD);
                assert dividesQ(bMod, modularGCD);

                if (modularGCD.isConstant())
                    return gcdInput.restoreGCD(a.createOne());

                // better degree bound found -> start over
                if (modularGCD.degree(0) < bBase.degree(0)) {
                    lLcGCD = lcGCD.mod(bPrime).longValueExact();
                    // scale to correct l.c.
                    bBase = modularGCD.monic(lLcGCD).toBigPoly();
                    bBasePrime = bPrime;
                    continue;
                }

                //skip unlucky prime
                if (modularGCD.degree(0) > bBase.degree(0))
                    continue;

                //lifting
                BigInteger newBasePrime = bBasePrime.multiply(bPrime);
                long monicFactor = lDomain.multiply(
                        lDomain.reciprocal(modularGCD.lc()),
                        lcGCD.mod(bPrime).longValueExact());

                PairIterator<MonomialTerm<BigInteger>, MultivariatePolynomial<BigInteger>> iterator = new PairIterator<>(bBase, modularGCD.toBigPoly());
                while (iterator.hasNext()) {
                    iterator.advance();

                    MonomialTerm<BigInteger>
                            baseTerm = iterator.aTerm,
                            imageTerm = iterator.bTerm;

                    if (baseTerm.coefficient.isZero())
                        // term is absent in the base
                        continue;

                    if (imageTerm.coefficient.isZero()) {
                        // term is absent in the modularGCD => remove it from the base
                        bBase.subtract(baseTerm);
                        continue;
                    }

                    long oth = lDomain.multiply(imageTerm.coefficient.longValueExact(), monicFactor);

                    // update base term
                    BigInteger newCoeff = ChineseRemainders.ChineseRemainders(bBasePrime, bPrime, baseTerm.coefficient, BigInteger.valueOf(oth));
                    bBase.terms.add(baseTerm.setCoefficient(newCoeff));
                }

                bBase = bBase.setDomainUnsafe(new IntegersModulo(newBasePrime));
                bBasePrime = newBasePrime;

                // two trials didn't change the result, probably we are done
                MultivariatePolynomial<BigInteger> candidate = MultivariatePolynomial.asPolyZSymmetric(bBase).primitivePart();
                if (previousBase != null && candidate.equals(previousBase)) {
                    previousBase = candidate;
                    //first check b since b is less degree
                    if (!dividesQ(b, candidate))
                        continue;

                    if (!dividesQ(a, candidate))
                        continue;

                    return gcdInput.restoreGCD(candidate);
                }
                previousBase = candidate;
            }
        }
    }

    static lMultivariatePolynomialZp interpolateGCD(lMultivariatePolynomialZp a, lMultivariatePolynomialZp b, lMultivariatePolynomialZp skeleton, RandomGenerator rnd) {
        a.checkSameDomainWith(b);
        a.checkSameDomainWith(skeleton);

        lMultivariatePolynomialZp content = contentGCD(a, b, 0, MultivariateGCD::PolynomialGCD);
        a = divideExact(a, content);
        b = divideExact(b, content);
        skeleton = divideSkeletonExact(skeleton, content);

        lSparseInterpolation interpolation = createInterpolation(-1, a, b, skeleton, rnd);
        if (interpolation == null)
            return null;
        lMultivariatePolynomialZp gcd = interpolation.evaluate();
        if (gcd == null)
            return null;

        return gcd.multiply(content);
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> Poly divideSkeletonExact(Poly dividend, Poly divider) {
        if (divider.isConstant())
            return dividend;
        if (divider.isMonomial())
            return dividend.clone().divideDegreeVectorOrNull(divider.lt());

        dividend = dividend.clone().setAllCoefficientsToUnit();
        divider = divider.clone().setAllCoefficientsToUnit();

        Poly quotient = dividend.createZero();
        dividend = dividend.clone();
        while (!dividend.isZero()) {
            Term dlt = dividend.lt();
            Term ltDiv = dlt.divide(divider.lt());
            if (ltDiv == null)
                throw new RuntimeException();
            quotient = quotient.add(ltDiv);
            dividend = dividend.subtract(divider.clone().multiply(ltDiv));
        }
        return quotient;
    }

    static final class PairIterator<Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>> {
        final Poly factory;
        final Term zeroTerm;
        final Comparator<DegreeVector> ordering;
        final Iterator<Term> aIterator, bIterator;

        PairIterator(Poly a, Poly b) {
            this.factory = a;
            this.zeroTerm = factory.getZeroTerm();
            this.ordering = a.ordering;
            this.aIterator = a.iterator();
            this.bIterator = b.iterator();
        }

        boolean hasNext() {
            return aIterator.hasNext() || bIterator.hasNext();
        }

        Term aTerm = null, bTerm = null;
        Term aTermCached = null, bTermCached = null;

        void advance() {
            if (aTermCached != null) {
                aTerm = aTermCached;
                aTermCached = null;
            } else
                aTerm = aIterator.hasNext() ? aIterator.next() : zeroTerm;
            if (bTermCached != null) {
                bTerm = bTermCached;
                bTermCached = null;
            } else
                bTerm = bIterator.hasNext() ? bIterator.next() : zeroTerm;

            int c = ordering.compare(aTerm, bTerm);
            if (c < 0) {
                aTermCached = aTerm;
                aTerm = zeroTerm;
            } else if (c > 0) {
                bTermCached = bTerm;
                bTerm = zeroTerm;
            } else {
                aTermCached = null;
                bTermCached = null;
            }
        }
    }


    /* ======================== Multivariate GCD over finite fields with small cardinality ========================== */

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("ConstantConditions")
    public static <E> MultivariatePolynomial<E> ModularGCDFiniteField(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b) {
        CommonUtils.ensureFiniteFieldDomain(a, b);
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree() < b.degree())
            return ModularGCDFiniteField(b, a);

        GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>>
                gcdInput = preparedGCDInput(a, b, MultivariateGCD::ModularGCDFiniteField);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;
        return ModularGCDFiniteField(gcdInput);
    }

    private static <E> MultivariatePolynomial<E> ModularGCDFiniteField(GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>> gcdInput) {
        MultivariatePolynomial<E> a = gcdInput.aReduced;
        MultivariatePolynomial<E> b = gcdInput.bReduced;

        int uVariable = a.nVariables - 1;
        if (a.domain instanceof IntegersModulo && a.coefficientDomainCardinality().isLong()) {
            // use machine integers
            @SuppressWarnings("unchecked")
            MultivariatePolynomial<lUnivariatePolynomialZp>
                    ua = asOverUnivariate0((MultivariatePolynomial<BigInteger>) a, uVariable),
                    ub = asOverUnivariate0((MultivariatePolynomial<BigInteger>) b, uVariable);

            lUnivariatePolynomialZp aContent = ua.content(), bContent = ub.content();
            lUnivariatePolynomialZp contentGCD = ua.domain.gcd(aContent, bContent);

            ua = ua.divideOrNull(aContent);
            ub = ub.divideOrNull(bContent);

            MultivariatePolynomial<lUnivariatePolynomialZp> ugcd =
                    ModularGCDFiniteField0(ua, ub, gcdInput.degreeBounds[uVariable], gcdInput.finiteExtensionDegree);
            ugcd = ugcd.multiply(contentGCD);

            @SuppressWarnings("unchecked")
            MultivariatePolynomial<E> r = gcdInput.restoreGCD((MultivariatePolynomial<E>) asNormalMultivariate0(ugcd, uVariable));
            return r;
        } else {
            MultivariatePolynomial<UnivariatePolynomial<E>>
                    ua = a.asOverUnivariateEliminate(uVariable),
                    ub = b.asOverUnivariateEliminate(uVariable);

            UnivariatePolynomial<E> aContent = ua.content(), bContent = ub.content();
            UnivariatePolynomial<E> contentGCD = ua.domain.gcd(aContent, bContent);

            ua = ua.divideOrNull(aContent);
            ub = ub.divideOrNull(bContent);

            MultivariatePolynomial<UnivariatePolynomial<E>> ugcd =
                    ModularGCDFiniteField0(ua, ub, gcdInput.degreeBounds[uVariable], gcdInput.finiteExtensionDegree);
            ugcd = ugcd.multiply(contentGCD);

            return gcdInput.restoreGCD(MultivariatePolynomial.asNormalMultivariate(ugcd, uVariable));
        }
    }

    private static MultivariatePolynomial<lUnivariatePolynomialZp> asOverUnivariate0(MultivariatePolynomial<BigInteger> poly, int variable) {
        lIntegersModulo domain = new lIntegersModulo(((IntegersModulo) poly.domain).modulus.longValueExact());
        lUnivariatePolynomialZp factory = lUnivariatePolynomialZp.zero(domain);
        UnivariatePolynomials<lUnivariatePolynomialZp> pDomain = new UnivariatePolynomials<>(factory);
        MonomialsSet<MonomialTerm<lUnivariatePolynomialZp>> newData = new MonomialsSet<>(poly.ordering);
        for (MonomialTerm<BigInteger> e : poly) {
            add(newData, new MonomialTerm<>(
                            e.without(variable).exponents,
                            factory.createMonomial(e.coefficient.longValueExact(), e.exponents[variable])),
                    pDomain);
        }
        return new MultivariatePolynomial<>(poly.nVariables - 1, pDomain, poly.ordering, newData);
    }

    private static MultivariatePolynomial<BigInteger> asNormalMultivariate0(MultivariatePolynomial<lUnivariatePolynomialZp> poly, int variable) {
        Domain<BigInteger> domain = poly.domain.getZero().domain.asDomain();
        int nVariables = poly.nVariables + 1;
        MultivariatePolynomial<BigInteger> result = zero(nVariables, domain, poly.ordering);
        for (MonomialTerm<lUnivariatePolynomialZp> entry : poly.terms) {
            lUnivariatePolynomialZp uPoly = entry.coefficient;
            int[] dv = ArraysUtil.insert(entry.exponents, variable, 0);
            for (int i = 0; i <= uPoly.degree(); ++i) {
                if (uPoly.isZeroAt(i))
                    continue;
                int[] cdv = dv.clone();
                cdv[variable] = i;
                result.add(new MonomialTerm<>(cdv, BigInteger.valueOf(uPoly.get(i))));
            }
        }
        return result;
    }

    /**
     * Modular GCD algorithm for polynomials over Z.
     *
     * @param a the first polynomial
     * @param b the second polynomial
     * @return GCD of two polynomials
     */
    @SuppressWarnings("ConstantConditions")
    public static lMultivariatePolynomialZp ModularGCDFiniteField(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b) {
        CommonUtils.ensureFiniteFieldDomain(a, b);
        if (a == b)
            return a.clone();
        if (a.isZero()) return b.clone();
        if (b.isZero()) return a.clone();

        if (a.degree() < b.degree())
            return ModularGCDFiniteField(b, a);

        GCDInput<lMonomialTerm, lMultivariatePolynomialZp>
                gcdInput = preparedGCDInput(a, b, MultivariateGCD::ModularGCDFiniteField);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;
        return lModularGCDFiniteField(gcdInput);
    }

    private static lMultivariatePolynomialZp lModularGCDFiniteField(GCDInput<lMonomialTerm, lMultivariatePolynomialZp> gcdInput) {
        lMultivariatePolynomialZp a = gcdInput.aReduced;
        lMultivariatePolynomialZp b = gcdInput.bReduced;

        int uVariable = a.nVariables - 1;
        MultivariatePolynomial<lUnivariatePolynomialZp>
                ua = a.asOverUnivariateEliminate(uVariable),
                ub = b.asOverUnivariateEliminate(uVariable);

        lUnivariatePolynomialZp aContent = ua.content(), bContent = ub.content();
        lUnivariatePolynomialZp contentGCD = ua.domain.gcd(aContent, bContent);

        ua = ua.divideOrNull(aContent);
        ub = ub.divideOrNull(bContent);

        MultivariatePolynomial<lUnivariatePolynomialZp> ugcd =
                ModularGCDFiniteField0(ua, ub, gcdInput.degreeBounds[uVariable], gcdInput.finiteExtensionDegree);
        ugcd = ugcd.multiply(contentGCD);

        return gcdInput.restoreGCD(lMultivariatePolynomialZp.asNormalMultivariate(ugcd, uVariable));
    }

    private static final int MAX_OVER_ITERATIONS = 16;

    static <uPoly extends IUnivariatePolynomial<uPoly>>
    MultivariatePolynomial<uPoly> ModularGCDFiniteField0(
            MultivariatePolynomial<uPoly> a,
            MultivariatePolynomial<uPoly> b,
            int uDegreeBound,
            int finiteExtensionDegree) {

        uPoly lcGCD = UnivariateGCD.PolynomialGCD(a.lc(), b.lc());

        Domain<uPoly> univariateDomain = a.domain;
        RandomGenerator random = PrivateRandom.getRandom();
        IrreduciblePolynomialsIterator<uPoly> primesLoop = new IrreduciblePolynomialsIterator<>(univariateDomain.getOne(), finiteExtensionDegree);
        main_loop:
        while (true) {
            // prepare the skeleton
            uPoly basePrime = primesLoop.next();

            FiniteField<uPoly> fField = new FiniteField<>(basePrime);

            // reduce Zp[x] -> GF
            MultivariatePolynomial<uPoly>
                    aMod = a.setDomain(fField),
                    bMod = b.setDomain(fField);
            if (!aMod.sameSkeleton(a) || !bMod.sameSkeleton(b))
                continue;

            //accumulator to update coefficients via Chineese remainding
            MultivariatePolynomial<uPoly> base = PolynomialGCD(aMod, bMod);
            uPoly lLcGCD = fField.valueOf(lcGCD);

            if (base.isConstant())
                return a.createOne();

            // scale to correct l.c.
            base = base.monic(lLcGCD);

            // cache the previous base
            MultivariatePolynomial<uPoly> previousBase = null;

            // over all primes
            while (true) {
                if (basePrime.degree() >= uDegreeBound + MAX_OVER_ITERATIONS)
                    throw new RuntimeException();

                uPoly prime = primesLoop.next();
                fField = new FiniteField<>(prime);

                // reduce Zp[x] -> GF
                aMod = a.setDomain(fField);
                bMod = b.setDomain(fField);
                if (!aMod.sameSkeleton(a) || !bMod.sameSkeleton(b))
                    continue;

                // calculate new GCD using previously calculated skeleton via sparse interpolation
                MultivariatePolynomial<uPoly> modularGCD = interpolateGCD(aMod, bMod, base.setDomainUnsafe(fField), random);
                if (modularGCD == null) {
                    // interpolation failed => assumed form is wrong => start over
                    continue main_loop;
                }

                assert dividesQ(aMod, modularGCD);
                assert dividesQ(bMod, modularGCD);

                if (modularGCD.isConstant())
                    return a.createOne();

                // better degree bound found -> start over
                if (modularGCD.degree(0) < base.degree(0)) {
                    lLcGCD = fField.valueOf(lcGCD);
                    // scale to correct l.c.
                    base = modularGCD.monic(lLcGCD);
                    basePrime = prime;
                    continue;
                }

                //skip unlucky prime
                if (modularGCD.degree(0) > base.degree(0))
                    continue;

                //lifting
                uPoly newBasePrime = basePrime.clone().multiply(prime);
                uPoly monicFactor = fField.divideExact(fField.valueOf(lcGCD), modularGCD.lc());

                PairIterator<MonomialTerm<uPoly>, MultivariatePolynomial<uPoly>> iterator = new PairIterator<>(base, modularGCD);
                while (iterator.hasNext()) {
                    iterator.advance();

                    MonomialTerm<uPoly>
                            baseTerm = iterator.aTerm,
                            imageTerm = iterator.bTerm;

                    if (baseTerm.coefficient.isZero())
                        // term is absent in the base
                        continue;

                    if (imageTerm.coefficient.isZero()) {
                        // term is absent in the modularGCD => remove it from the base
                        base.subtract(baseTerm);
                        continue;
                    }

                    uPoly oth = fField.multiply(imageTerm.coefficient, monicFactor);

                    // update base term
                    uPoly newCoeff = ChineseRemainders.ChineseRemainders(univariateDomain, basePrime, prime, baseTerm.coefficient, oth);
                    base.terms.add(baseTerm.setCoefficient(newCoeff));
                }

                basePrime = newBasePrime;

                // set domain back to the normal univariate domain
                base = base.setDomainUnsafe(univariateDomain);
                // two trials didn't change the result, probably we are done
                MultivariatePolynomial<uPoly> candidate = base.clone().primitivePart();
                if (basePrime.degree() >= uDegreeBound || (previousBase != null && candidate.equals(previousBase))) {
                    previousBase = candidate;
                    //first check b since b is less degree
                    if (!dividesQ(b, candidate))
                        continue;

                    if (!dividesQ(a, candidate))
                        continue;


                    return candidate;
                }
                previousBase = candidate;
            }
        }
    }

    private static final class IrreduciblePolynomialsIterator<uPoly extends IUnivariatePolynomial<uPoly>> {
        final uPoly factory;
        int degree;

        IrreduciblePolynomialsIterator(uPoly factory, int degree) {
            this.factory = factory;
            this.degree = degree;
        }

        uPoly next() {
            return IrreduciblePolynomials.randomIrreduciblePolynomial(factory, degree++, PrivateRandom.getRandom());
        }
    }

    static <E> MultivariatePolynomial<E> interpolateGCD(MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, MultivariatePolynomial<E> skeleton, RandomGenerator rnd) {
        a.checkSameDomainWith(b);
        a.checkSameDomainWith(skeleton);

        MultivariatePolynomial<E> content = contentGCD(a, b, 0, MultivariateGCD::ZippelGCD);
        a = divideExact(a, content);
        b = divideExact(b, content);
        skeleton = divideSkeletonExact(skeleton, content);

        SparseInterpolation<E> interpolation = createInterpolation(-1, a, b, skeleton, rnd);
        if (interpolation == null)
            return null;
        MultivariatePolynomial<E> gcd = interpolation.evaluate();
        if (gcd == null)
            return null;

        return gcd.multiply(content);
    }

    /* ===================================== Multivariate GCD over finite fields ==================================== */

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Brown's algorithm with dense interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static <E> MultivariatePolynomial<E> BrownGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b) {
        CommonUtils.ensureFieldDomain(a, b);

        // prepare input and test for early termination
        GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>> gcdInput = preparedGCDInput(a, b, MultivariateGCD::BrownGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        if (gcdInput.finiteExtensionDegree > 1)
            return ModularGCDFiniteField(gcdInput);

        MultivariatePolynomial<E> result = BrownGCD(
                gcdInput.aReduced, gcdInput.bReduced, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            // ground fill is too small for modular algorithm
            return ModularGCDFiniteField(gcdInput);

        return gcdInput.restoreGCD(result);
    }

    /**
     * Actual implementation of dense interpolation
     *
     * @param variable             current variable (all variables {@code v > variable} are fixed so far)
     * @param degreeBounds         degree bounds for gcd
     * @param evaluationStackLimit domain cardinality
     */
    @SuppressWarnings("unchecked")
    private static <E> MultivariatePolynomial<E> BrownGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        MultivariatePolynomial<E> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        MultivariatePolynomial<E> factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0)
        // switch to univariate gcd
        {
            UnivariatePolynomial<E> gcd = UnivariateGCD.PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return AMultivariatePolynomial.asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        PrimitiveInput<MonomialTerm<E>, MultivariatePolynomial<E>, UnivariatePolynomial<E>> primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        UnivariatePolynomial<E>
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            MultivariatePolynomial<E> poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        Domain<E> domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //dense interpolation
        Interpolation<E> interpolation = null;
        //previous interpolation (used to detect whether update doesn't change the result)
        MultivariatePolynomial<E> previousInterpolation;
        //store points that were already used in interpolation
        Set<E> evaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // all elements of the domain are tried
                // do division check (last chance) and return
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            //pickup the next random element for variable
            E randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;
            evaluationStack.add(randomPoint);

            E lcVal = lcGCD.evaluate(randomPoint);
            if (domain.isZero(lcVal))
                continue;

            // evaluate a and b at variable = randomPoint
            MultivariatePolynomial<E>
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<E> cVal = BrownGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            // normalize gcd
            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            if (currExponent < prevVarExponent) {
                //better degree bound detected => start over
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
                continue;
            }

            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new Interpolation<>(variable, randomPoint, cVal);
                continue;
            }

            // Cache previous interpolation. NOTE: clone() is important, since the poly will
            // be modified inplace by the update() method
            previousInterpolation = interpolation.getInterpolatingPolynomial().clone();
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] <= interpolation.numberOfPoints()
                    || previousInterpolation.equals(interpolation.getInterpolatingPolynomial())) {
                MultivariatePolynomial<E> result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    /** division test **/
    private static <E> MultivariatePolynomial<E> doDivisionCheck(
            MultivariatePolynomial<E> a, MultivariatePolynomial<E> b,
            UnivariatePolynomial<E> contentGCD, Interpolation<E> interpolation, int variable) {
        if (interpolation == null)
            return null;
        MultivariatePolynomial<E> interpolated =
                MultivariatePolynomial.asNormalMultivariate(interpolation.getInterpolatingPolynomial().asOverUnivariateEliminate(variable).primitivePart(), variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        if (contentGCD == null)
            return interpolated;
        MultivariatePolynomial<E> poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
        return interpolated.multiply(poly);
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Zippel's algorithm with sparse interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static <E> MultivariatePolynomial<E> ZippelGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b) {
        CommonUtils.ensureFieldDomain(a, b);

        // prepare input and test for early termination
        GCDInput<MonomialTerm<E>, MultivariatePolynomial<E>> gcdInput = preparedGCDInput(a, b, MultivariateGCD::ZippelGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        if (gcdInput.finiteExtensionDegree > 1)
            return ModularGCDFiniteField(gcdInput);

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

        // content in the main variable => avoid raise condition in LINZIP!
        // see Example 4 in "Algorithms for the non-monic case of the sparse modular GCD algorithm"
        MultivariatePolynomial<E> content = contentGCD(a, b, 0, MultivariateGCD::ZippelGCD);
        a = divideExact(a, content);
        b = divideExact(b, content);

        MultivariatePolynomial<E> result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            // ground fill is too small for modular algorithm
            return ModularGCDFiniteField(gcdInput);

        result = result.multiply(content);
        return gcdInput.restoreGCD(result);
    }

    /** Maximal number of fails before switch to a new homomorphism */
    private static final int MAX_SPARSE_INTERPOLATION_FAILS = 1000;
    /** Maximal number of sparse interpolations after interpolation.numberOfPoints() > degreeBounds[variable] */
    private static final int ALLOWED_OVER_INTERPOLATED_ATTEMPTS = 32;

    @SuppressWarnings("unchecked")
    private static <E> MultivariatePolynomial<E> ZippelGCD(
            MultivariatePolynomial<E> a,
            MultivariatePolynomial<E> b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        MultivariatePolynomial<E> trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        MultivariatePolynomial<E> factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0)
        // switch to univariate gcd
        {
            UnivariatePolynomial<E> gcd = UnivariateGCD.PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return AMultivariatePolynomial.asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

//        MultivariatePolynomial<E> content = ZippelContentGCD(a, b, variable);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        PrimitiveInput<
                MonomialTerm<E>,
                MultivariatePolynomial<E>,
                UnivariatePolynomial<E>>
                primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        UnivariatePolynomial<E>
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            MultivariatePolynomial<E> poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        Domain<E> domain = factory.domain;
        //store points that were already used in interpolation
        Set<E> globalEvaluationStack = new HashSet<>();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        int failedSparseInterpolations = 0;

        int[] tmpDegreeBounds = degreeBounds.clone();
        main:
        while (true) {
            if (evaluationStackLimit == globalEvaluationStack.size())
                return null;

            E seedPoint = domain.randomElement(rnd);
            if (globalEvaluationStack.contains(seedPoint))
                continue;

            globalEvaluationStack.add(seedPoint);

            E lcVal = lcGCD.evaluate(seedPoint);
            if (domain.isZero(lcVal))
                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            MultivariatePolynomial<E>
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            MultivariatePolynomial<E> cVal = ZippelGCD(aVal, bVal, rnd, variable - 1, tmpDegreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > tmpDegreeBounds[variable - 1])
                //unlucky homomorphism
                continue;

            if (currExponent < tmpDegreeBounds[variable - 1]) {
                //better degree bound detected
                tmpDegreeBounds[variable - 1] = currExponent;
            }

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc().equals(lcVal);

            SparseInterpolation sparseInterpolator = createInterpolation(variable, a, b, cVal, rnd);
            if (sparseInterpolator == null)
                //unlucky homomorphism
                continue;

            // we are applying dense interpolation for univariate skeleton coefficients
            Interpolation<E> denseInterpolation = new Interpolation<>(variable, seedPoint, cVal);
            //previous interpolation (used to detect whether update doesn't change the result)
            MultivariatePolynomial<E> previousInterpolation;
            //local evaluation stack for points that are calculated via sparse interpolation (but not gcd evaluation) -> always same skeleton
            HashSet<E> localEvaluationStack = new HashSet<>(globalEvaluationStack);
            while (true) {
                if (evaluationStackLimit == localEvaluationStack.size())
                    return null;

                if (denseInterpolation.numberOfPoints() > tmpDegreeBounds[variable] + ALLOWED_OVER_INTERPOLATED_ATTEMPTS) {
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                E randomPoint = domain.randomElement(rnd);
                if (localEvaluationStack.contains(randomPoint))
                    continue;
                localEvaluationStack.add(randomPoint);

                lcVal = lcGCD.evaluate(randomPoint);
                if (domain.isZero(lcVal))
                    continue;

                cVal = sparseInterpolator.evaluate(randomPoint);
                if (cVal == null) {
                    ++failedSparseInterpolations;
                    if (failedSparseInterpolations == MAX_SPARSE_INTERPOLATION_FAILS)
                        throw new RuntimeException("Sparse interpolation failed");
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc().equals(lcVal);

                // Cache previous interpolation. NOTE: clone() is important, since the poly will
                // be modified inplace by the update() method
                previousInterpolation = denseInterpolation.getInterpolatingPolynomial().clone();
                denseInterpolation.update(randomPoint, cVal);

                // do division test
                if (tmpDegreeBounds[variable] <= denseInterpolation.numberOfPoints()
                        || previousInterpolation.equals(denseInterpolation.getInterpolatingPolynomial())) {
                    MultivariatePolynomial<E> result = doDivisionCheck(a, b, contentGCD, denseInterpolation, variable);
                    if (result != null)
                        return result;
                }
            }
        }
    }

    static boolean ALWAYS_LINZIP = false;
    /** maximal number of attempts to choose a good evaluation point for sparse interpolation */
    private static final int MAX_FAILED_SUBSTITUTIONS = 32;

    static <E> SparseInterpolation<E> createInterpolation(int variable,
                                                          MultivariatePolynomial<E> a,
                                                          MultivariatePolynomial<E> b,
                                                          MultivariatePolynomial<E> skeleton,
                                                          RandomGenerator rnd) {
        skeleton = skeleton.clone().setAllCoefficientsToUnit();
        if (skeleton.size() == 1)
            return new TrivialSparseInterpolation<>(skeleton);

        boolean monic = a.coefficientOf(0, a.degree(0)).isConstant() && b.coefficientOf(0, a.degree(0)).isConstant();
        Set<DegreeVector> globalSkeleton = skeleton.terms.keySet();
        TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        Domain<E> domain = a.domain;

        int lastVariable = variable == -1 ? a.nVariables - 1 : variable;
        int[] evaluationVariables = ArraysUtil.sequence(1, lastVariable + 1);//variable inclusive
        E[] evaluationPoint = domain.createArray(evaluationVariables.length);

        PrecomputedPowersHolder<E> powers;
        int fails = 0;
        search_for_good_evaluation_point:
        while (true) {
            if (fails >= MAX_FAILED_SUBSTITUTIONS)
                return null;
            //avoid zero evaluation points
            for (int i = lastVariable - 1; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (domain.isZero(evaluationPoint[i]));

            powers = new PrecomputedPowersHolder<>(a.nVariables, evaluationVariables, evaluationPoint, domain);
            int[] raiseFactors = ArraysUtil.arrayOf(1, evaluationVariables.length);

            for (MultivariatePolynomial<E> p : Arrays.asList(a, b, skeleton))
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables, raiseFactors).getSkeleton())) {
                    ++fails;
                    continue search_for_good_evaluation_point;
                }
            break;
        }

        int requiredNumberOfEvaluations = -1, monicScalingExponent = -1;
        for (TIntObjectIterator<MultivariatePolynomial<E>> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            MultivariatePolynomial<E> v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
            if (v.size() == 1)
                monicScalingExponent = it.key();
        }

        if (!ALWAYS_LINZIP) {
            if (monic)
                monicScalingExponent = -1;

            if (monic || monicScalingExponent != -1)
                return new MonicInterpolation<>(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                        evaluationVariables, evaluationPoint, powers, rnd, requiredNumberOfEvaluations, monicScalingExponent);
        }

        return new LinZipInterpolation<>(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, rnd);
    }

    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
    private static <E> TIntObjectHashMap<MultivariatePolynomial<E>> getSkeleton(MultivariatePolynomial<E> poly) {
        TIntObjectHashMap<MultivariatePolynomial<E>> skeleton = new TIntObjectHashMap<>();
        for (MonomialTerm<E> term : poly.terms) {
            MonomialTerm<E> newDV = term.setZero(0);
            MultivariatePolynomial<E> coeff = skeleton.get(term.exponents[0]);
            if (coeff != null)
                coeff.add(newDV);
            else
                skeleton.put(term.exponents[0], MultivariatePolynomial.create(poly.nVariables,
                        poly.domain, poly.ordering, newDV));
        }
        return skeleton;
    }

    interface SparseInterpolation<E> {
        /**
         * Returns interpolating gcd
         *
         * @return gcd
         */
        MultivariatePolynomial<E> evaluate();

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        MultivariatePolynomial<E> evaluate(E newPoint);
    }

    static final class TrivialSparseInterpolation<E> implements SparseInterpolation<E> {
        final MultivariatePolynomial<E> val;

        TrivialSparseInterpolation(MultivariatePolynomial<E> val) {
            this.val = val;
        }

        @Override
        public MultivariatePolynomial<E> evaluate() {
            return val;
        }

        @Override
        public MultivariatePolynomial<E> evaluate(E newPoint) {
            return val;
        }
    }

    static abstract class ASparseInterpolation<E> implements SparseInterpolation<E> {
        /** the domain */
        final Domain<E> domain;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final MultivariatePolynomial<E> a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final E[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final PrecomputedPowersHolder<E> powers;
        /** random */
        final RandomGenerator rnd;

        ASparseInterpolation(Domain<E> domain, int variable,
                             MultivariatePolynomial<E> a, MultivariatePolynomial<E> b,
                             Set<DegreeVector> globalSkeleton,
                             TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton,
                             int[] sparseUnivarDegrees, int[] evaluationVariables,
                             E[] evaluationPoint,
                             PrecomputedPowersHolder<E> powers, RandomGenerator rnd) {
            this.domain = domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.globalSkeleton = globalSkeleton;
            this.univarSkeleton = univarSkeleton;
            this.sparseUnivarDegrees = sparseUnivarDegrees;
            this.evaluationVariables = evaluationVariables;
            this.evaluationPoint = evaluationPoint;
            this.powers = powers;
            this.rnd = rnd;
        }

        @SuppressWarnings("unchecked")
        @Override
        public final MultivariatePolynomial<E> evaluate(E newPoint) {
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationVariables[evaluationVariables.length - 1], newPoint);
            return evaluate();
        }
    }

    /** Number of retries when raise condition occurs; then drop up with new homomorphism */
    private static final int NUMBER_OF_UNDER_DETERMINED_RETRIES = 8;

    static final class LinZipInterpolation<E> extends ASparseInterpolation<E> {
        LinZipInterpolation(Domain<E> domain, int variable, MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, E[] evaluationPoint, PrecomputedPowersHolder<E> powers, RandomGenerator rnd) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
        }

        @Override
        public MultivariatePolynomial<E> evaluate() {
            @SuppressWarnings("unchecked")
            LinZipSystem<E>[] systems = new LinZipSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new LinZipSystem<>(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? a.nVariables - 1 : variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            if (variable != -1)
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;

            int nUnknowns = globalSkeleton.size(), nUnknownScalings = -1;
            int raiseFactor = 0;

            for (int nTries = 0; nTries < NUMBER_OF_UNDER_DETERMINED_RETRIES; ++nTries) {
                int previousFreeVars = -1, underDeterminedTries = 0;
                for (; ; ) {
                    // increment at each loop!
                    ++nUnknownScalings;
                    ++raiseFactor;
                    // sequential powers of evaluation point
                    Arrays.fill(raiseFactors, 0, variable == -1 ? raiseFactors.length : raiseFactors.length - 1, raiseFactor);
                    // evaluate a and b to univariate and calculate gcd
                    UnivariatePolynomial<E>
                            aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            gcdUnivar = UnivariateGCD.PolynomialGCD(aUnivar, bUnivar);

                    if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                        // unlucky main homomorphism or bad evaluation point
                        return null;

                    assert gcdUnivar.isMonic();
                    if (!univarSkeleton.keySet().containsAll(gcdUnivar.exponents()))
                        // univariate gcd contains terms that are not present in the skeleton
                        // again unlucky main homomorphism
                        return null;

                    int totalEquations = 0;
                    for (LinZipSystem<E> system : systems) {
                        E rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs, nUnknownScalings != 0);
                        totalEquations += system.nEquations();
                    }
                    if (nUnknowns + nUnknownScalings <= totalEquations)
                        break;


                    if (underDeterminedTries > NUMBER_OF_UNDER_DETERMINED_RETRIES) {
                        // raise condition: new equations does not fix enough variables
                        return null;
                    }

                    int freeVars = nUnknowns + nUnknownScalings - totalEquations;
                    if (freeVars >= previousFreeVars)
                        ++underDeterminedTries;
                    else
                        underDeterminedTries = 0;

                    previousFreeVars = freeVars;
                }

                MultivariatePolynomial<E> result = a.createZero();
                SystemInfo info = solveLinZip(a, systems, nUnknownScalings, result);
                if (info == SystemInfo.UnderDetermined) {
                    //try to generate more equations
                    continue;
                }

                if (info == SystemInfo.Consistent)
                    //well done
                    return result;
                if (info == SystemInfo.Inconsistent)
                    //inconsistent system => unlucky homomorphism
                    return null;
            }

            // the system is still under-determined => bad evaluation homomorphism
            return null;
        }
    }

    private static <E> SystemInfo solveLinZip(MultivariatePolynomial<E> factory, LinZipSystem<E>[] subSystems, int nUnknownScalings, MultivariatePolynomial<E> destination) {
        ArrayList<MonomialTerm<E>> unknowns = new ArrayList<>();
        for (LinZipSystem<E> system : subSystems)
            for (MonomialTerm<E> degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = nUnknownsMonomials + nUnknownScalings;
        ArrayList<E[]> lhsGlobal = new ArrayList<>();
        ArrayList<E> rhsGlobal = new ArrayList<>();
        int offset = 0;
        Domain<E> domain = factory.domain;
        for (LinZipSystem<E> system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                E[] row = domain.createZeroesArray(nUnknownsTotal);
                E[] subRow = system.matrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                if (j > 0)
                    row[nUnknownsMonomials + j - 1] = system.scalingMatrix.get(j);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        E[] solution = domain.createArray(nUnknownsTotal);
        SystemInfo info = LinearAlgebra.solve(domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent) {
            @SuppressWarnings("unchecked")
            MonomialTerm<E>[] terms = new MonomialTerm[unknowns.size()];
            for (int i = 0; i < terms.length; i++)
                terms[i] = unknowns.get(i).setCoefficient(solution[i]);
            destination.add(terms);
        }
        return info;
    }


    static final class MonicInterpolation<E> extends ASparseInterpolation<E> {
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** univar exponent with monomial factor that can be used for scaling */
        final int monicScalingExponent;

        MonicInterpolation(Domain<E> domain, int variable, MultivariatePolynomial<E> a, MultivariatePolynomial<E> b, Set<DegreeVector> globalSkeleton, TIntObjectHashMap<MultivariatePolynomial<E>> univarSkeleton, int[] sparseUnivarDegrees, int[] evaluationVariables, E[] evaluationPoint, PrecomputedPowersHolder<E> powers, RandomGenerator rnd, int requiredNumberOfEvaluations, int monicScalingExponent) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
            this.requiredNumberOfEvaluations = requiredNumberOfEvaluations;
            this.monicScalingExponent = monicScalingExponent;
        }

        @Override
        public MultivariatePolynomial<E> evaluate() {
            @SuppressWarnings("unchecked")
            VandermondeSystem<E>[] systems = new VandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new VandermondeSystem<>(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? a.nVariables - 1 : variable - 1);

            int[] raiseFactors = new int[evaluationVariables.length];
            if (variable != -1)
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;
            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                Arrays.fill(raiseFactors, 0, variable == -1 ? raiseFactors.length : raiseFactors.length - 1, i + 1);
                // evaluate a and b to univariate and calculate gcd
                UnivariatePolynomial<E>
                        aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        gcdUnivar = UnivariateGCD.PolynomialGCD(aUnivar, bUnivar);

                if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                    // unlucky main homomorphism or bad evaluation point
                    return null;

                assert gcdUnivar.isMonic();
                if (!univarSkeleton.keySet().containsAll(gcdUnivar.exponents()))
                    // univariate gcd contains terms that are not present in the skeleton
                    // again unlucky main homomorphism
                    return null;

                if (monicScalingExponent != -1) {
                    // single scaling factor
                    // scale the system according to it

                    if (gcdUnivar.degree() < monicScalingExponent || domain.isZero(gcdUnivar.get(monicScalingExponent)))
                        // unlucky homomorphism
                        return null;

                    E normalization = evaluateExceptFirst(domain, powers, domain.getOne(), univarSkeleton.get(monicScalingExponent).lt(), i + 1, variable == -1 ? a.nVariables - 1 : variable - 1);
                    //normalize univariate gcd in order to reconstruct leading coefficient polynomial
                    normalization = domain.multiply(domain.reciprocal(gcdUnivar.get(monicScalingExponent)), normalization);
                    gcdUnivar = gcdUnivar.multiply(normalization);
                }

                boolean allDone = true;
                for (VandermondeSystem<E> system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        E rhs = gcdUnivar.degree() < system.univarDegree ? domain.getZero() : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (VandermondeSystem<E> system : systems) {
                //solve each system
                SystemInfo info = system.solve();
                if (info != SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            MultivariatePolynomial<E> gcdVal = a.createZero();
            for (VandermondeSystem<E> system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || domain.isOne(system.solution[0]);
                for (int i = 0; i < system.skeleton.length; i++) {
                    MonomialTerm<E> degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    E value = system.solution[i];
                    gcdVal.add(degreeVector.setCoefficient(value));
                }
            }

            return gcdVal;
        }
    }

    private static abstract class LinearSystem<E> {
        final int univarDegree;
        /** the domain */
        final Domain<E> domain;
        /** the skeleton */
        final MonomialTerm<E>[] skeleton;
        /** the lhs matrix */
        final ArrayList<E[]> matrix;
        /** the rhs values */
        final ArrayList<E> rhs = new ArrayList<>();
        /** precomputed powers */
        final PrecomputedPowersHolder<E> powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        final int nVars;

        @SuppressWarnings("unchecked")
        LinearSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            //todo refactor generics
            this.skeleton = skeleton.getSkeleton().toArray(new MonomialTerm[skeleton.size()]);
            this.powers = powers;
            this.nVars = nVars;
            this.matrix = new ArrayList<>();
        }

        final int nUnknownVariables() {
            return skeleton.length;
        }

        final int nEquations() {
            return matrix.size();
        }

        @Override
        public String toString() {
            return "{" + matrix.stream().map(Arrays::toString).collect(Collectors.joining(",")) + "} = " + rhs;
        }
    }

    /** Vandermonde system builder */
    private static final class LinZipSystem<E> extends LinearSystem<E> {
        public LinZipSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        private final ArrayList<E> scalingMatrix = new ArrayList<>();

        public void oneMoreEquation(E rhsVal, boolean newScalingIntroduced) {
            E[] row = domain.createArray(skeleton.length);
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);

            if (newScalingIntroduced) {
                scalingMatrix.add(domain.negate(rhsVal));
                rhsVal = domain.getZero();
            } else
                scalingMatrix.add(domain.getZero());
            rhs.add(rhsVal);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < matrix.size(); i++) {
                E[] row = matrix.get(i);
                for (int j = 0; j < row.length; j++) {
                    if (j != 0)
                        sb.append("+");
                    sb.append(row[j]).append("*c" + j);
                }
                if (i != 0)
                    sb.append("+").append(scalingMatrix.get(i)).append("*m" + (i - 1));

                sb.append("=").append(rhs.get(i)).append("\n");
            }
            return sb.toString();
        }
    }

    /** Vandermonde system builder */
    private static final class VandermondeSystem<E> extends LinearSystem<E> {
        public VandermondeSystem(int univarDegree, MultivariatePolynomial<E> skeleton, PrecomputedPowersHolder<E> powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        E[] solution = null;

        SystemInfo solve() {
            if (solution == null)
                solution = domain.createArray(nUnknownVariables());

            if (nUnknownVariables() <= 8)
                // for small systems Gaussian elimination is indeed faster
                return LinearAlgebra.solve(domain, matrix.toArray(domain.createArray2d(matrix.size())), rhs.toArray(domain.createArray(rhs.size())), solution);

            // solve vandermonde system
            E[] vandermondeRow = matrix.get(0);
            SystemInfo info = LinearAlgebra.solveVandermondeT(domain, vandermondeRow, rhs.toArray(domain.createArray(rhs.size())), solution);
            if (info == SystemInfo.Consistent)
                for (int i = 0; i < solution.length; ++i)
                    solution[i] = domain.divideExact(solution[i], vandermondeRow[i]);

            return info;
        }

        public VandermondeSystem<E> oneMoreEquation(E rhsVal) {
            E[] row = domain.createArray(skeleton.length);
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, domain.getOne(), skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);
            rhs.add(rhsVal);
            return this;
        }
    }

    private static <E> E evaluateExceptFirst(Domain<E> domain,
                                             PrecomputedPowersHolder<E> powers,
                                             E coefficient,
                                             MonomialTerm<E> skeleton,
                                             int raiseFactor,
                                             int nVars) {
        E tmp = coefficient;
        for (int k = 1; k <= nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[k]));
        return tmp;
    }

    private static <E> boolean isVandermonde(E[][] lhs, Domain<E> domain) {
        for (int i = 1; i < lhs.length; i++) {
            for (int j = 0; j < lhs[0].length; j++) {
                if (!lhs[i][j].equals(domain.pow(lhs[0][j], i + 1)))
                    return false;
            }
        }
        return true;
    }


    /* ================================================ Machine numbers ============================================= */

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Brown's algorithm with dense interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static lMultivariatePolynomialZp BrownGCD(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b) {

        // prepare input and test for early termination
        GCDInput<lMonomialTerm, lMultivariatePolynomialZp> gcdInput = preparedGCDInput(a, b, MultivariateGCD::BrownGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        if (gcdInput.finiteExtensionDegree > 1)
            return lModularGCDFiniteField(gcdInput);

        lMultivariatePolynomialZp result = BrownGCD(
                gcdInput.aReduced, gcdInput.bReduced, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            // ground fill is too small for modular algorithm
            return lModularGCDFiniteField(gcdInput);

        return gcdInput.restoreGCD(result);
    }

    /**
     * Actual implementation of dense interpolation
     *
     * @param variable             current variable (all variables {@code v > variable} are fixed so far)
     * @param degreeBounds         degree bounds for gcd
     * @param evaluationStackLimit domain cardinality
     */
    @SuppressWarnings("unchecked")
    private static lMultivariatePolynomialZp BrownGCD(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        lMultivariatePolynomialZp trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        lMultivariatePolynomialZp factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0)
        // switch to univariate gcd
        {
            lUnivariatePolynomialZp gcd = UnivariateGCD.PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return AMultivariatePolynomial.asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

        PrimitiveInput<lMonomialTerm, lMultivariatePolynomialZp, lUnivariatePolynomialZp> primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        lUnivariatePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            lMultivariatePolynomialZp poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        lIntegersModulo domain = factory.domain;
        //degree bound for the previous variable
        int prevVarExponent = degreeBounds[variable - 1];
        //dense interpolation
        lInterpolation interpolation = null;
        //previous interpolation (used to detect whether update doesn't change the result)
        lMultivariatePolynomialZp previousInterpolation;
        //store points that were already used in interpolation
        TLongHashSet evaluationStack = new TLongHashSet();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        main:
        while (true) {
            if (evaluationStackLimit == evaluationStack.size())
                // all elements of the domain are tried
                // do division check (last chance) and return
                return doDivisionCheck(a, b, contentGCD, interpolation, variable);

            //pickup the next random element for variable
            long randomPoint = domain.randomElement(rnd);
            if (evaluationStack.contains(randomPoint))
                continue;
            evaluationStack.add(randomPoint);

            long lcVal = lcGCD.evaluate(randomPoint);
            if (lcVal == 0)
                continue;

            // evaluate a and b at variable = randomPoint
            lMultivariatePolynomialZp
                    aVal = a.evaluate(variable, randomPoint),
                    bVal = b.evaluate(variable, randomPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            // calculate gcd of the result by the recursive call
            lMultivariatePolynomialZp cVal = BrownGCD(aVal, bVal, rnd, variable - 1, degreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > prevVarExponent)
                //unlucky homomorphism
                continue;

            // normalize gcd
            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc() == lcVal;

            if (currExponent < prevVarExponent) {
                //better degree bound detected => start over
                interpolation = new lInterpolation(variable, randomPoint, cVal);
                degreeBounds[variable - 1] = prevVarExponent = currExponent;
                continue;
            }

            if (interpolation == null) {
                //first successful homomorphism
                interpolation = new lInterpolation(variable, randomPoint, cVal);
                continue;
            }

            // Cache previous interpolation. NOTE: clone() is important, since the poly will
            // be modified inplace by the update() method
            previousInterpolation = interpolation.getInterpolatingPolynomial().clone();
            // interpolate
            interpolation.update(randomPoint, cVal);

            // do division test
            if (degreeBounds[variable] <= interpolation.numberOfPoints()
                    || previousInterpolation.equals(interpolation.getInterpolatingPolynomial())) {
                lMultivariatePolynomialZp result = doDivisionCheck(a, b, contentGCD, interpolation, variable);
                if (result != null)
                    return result;
            }
        }
    }

    /** division test **/
    private static lMultivariatePolynomialZp doDivisionCheck(
            lMultivariatePolynomialZp a, lMultivariatePolynomialZp b,
            lUnivariatePolynomialZp contentGCD, lInterpolation interpolation, int variable) {
        if (interpolation == null)
            return null;
        lMultivariatePolynomialZp interpolated =
                lMultivariatePolynomialZp.asNormalMultivariate(interpolation.getInterpolatingPolynomial().asOverUnivariateEliminate(variable).primitivePart(), variable);
        if (!dividesQ(a, interpolated) || !dividesQ(b, interpolated))
            return null;

        if (contentGCD == null)
            return interpolated;
        lMultivariatePolynomialZp poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
        return interpolated.multiply(poly);
    }

    /**
     * Calculates GCD of two multivariate polynomials over Zp using Zippel's algorithm with sparse interpolation.
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static lMultivariatePolynomialZp ZippelGCD(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b) {

        // prepare input and test for early termination
        GCDInput<lMonomialTerm, lMultivariatePolynomialZp> gcdInput = preparedGCDInput(a, b, MultivariateGCD::ZippelGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        if (gcdInput.finiteExtensionDegree > 1)
            return lModularGCDFiniteField(gcdInput);

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

        // content in the main variable => avoid raise condition in LINZIP!
        // see Example 4 in "Algorithms for the non-monic case of the sparse modular GCD algorithm"
        lMultivariatePolynomialZp content = contentGCD(a, b, 0, MultivariateGCD::ZippelGCD);
        a = divideExact(a, content);
        b = divideExact(b, content);

        lMultivariatePolynomialZp result = ZippelGCD(a, b, PrivateRandom.getRandom(),
                gcdInput.lastPresentVariable, gcdInput.degreeBounds, gcdInput.evaluationStackLimit);
        if (result == null)
            // ground fill is too small for modular algorithm
            return lModularGCDFiniteField(gcdInput);

        result = result.multiply(content);
        return gcdInput.restoreGCD(result);
    }

    @SuppressWarnings("unchecked")
    private static lMultivariatePolynomialZp ZippelGCD(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b,
            RandomGenerator rnd,
            int variable,
            int[] degreeBounds,
            int evaluationStackLimit) {

        //check for trivial gcd
        lMultivariatePolynomialZp trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null)
            return trivialGCD;

        lMultivariatePolynomialZp factory = a;
        int nVariables = factory.nVariables;
        if (variable == 0)
        // switch to univariate gcd
        {
            lUnivariatePolynomialZp gcd = UnivariateGCD.PolynomialGCD(a.asUnivariate(), b.asUnivariate());
            if (gcd.degree() == 0)
                return factory.createOne();
            return AMultivariatePolynomial.asMultivariate(gcd, nVariables, variable, factory.ordering);
        }

//        lMultivariatePolynomialZp content = ZippelContentGCD(a, b, variable);
//        a = divideExact(a, content);
//        b = divideExact(b, content);

        PrimitiveInput<
                lMonomialTerm,
                lMultivariatePolynomialZp,
                lUnivariatePolynomialZp>
                primitiveInput = makePrimitive(a, b, variable);
        // primitive parts of a and b as Zp[x_k][x_1 ... x_{k-1}]
        a = primitiveInput.aPrimitive;
        b = primitiveInput.bPrimitive;
        // gcd of Zp[x_k] content and lc
        lUnivariatePolynomialZp
                contentGCD = primitiveInput.contentGCD,
                lcGCD = primitiveInput.lcGCD;

        //check again for trivial gcd
        trivialGCD = trivialGCD(a, b);
        if (trivialGCD != null) {
            lMultivariatePolynomialZp poly = AMultivariatePolynomial.asMultivariate(contentGCD, a.nVariables, variable, a.ordering);
            return trivialGCD.multiply(poly);
        }

        lIntegersModulo domain = factory.domain;
        //store points that were already used in interpolation
        TLongHashSet globalEvaluationStack = new TLongHashSet();

        int[] aDegrees = a.degrees(), bDegrees = b.degrees();
        int failedSparseInterpolations = 0;

        int[] tmpDegreeBounds = degreeBounds.clone();
        main:
        while (true) {
            if (evaluationStackLimit == globalEvaluationStack.size())
                return null;

            long seedPoint = domain.randomElement(rnd);
            if (globalEvaluationStack.contains(seedPoint))
                continue;

            globalEvaluationStack.add(seedPoint);

            long lcVal = lcGCD.evaluate(seedPoint);
            if (lcVal == 0)
                continue;

            // evaluate a and b at variable = randomPoint
            // calculate gcd of the result by the recursive call
            lMultivariatePolynomialZp
                    aVal = a.evaluate(variable, seedPoint),
                    bVal = b.evaluate(variable, seedPoint);

            // check for unlucky substitution
            int[] aValDegrees = aVal.degrees(), bValDegrees = bVal.degrees();
            for (int i = variable - 1; i >= 0; --i)
                if (aDegrees[i] != aValDegrees[i] || bDegrees[i] != bValDegrees[i])
                    continue main;

            lMultivariatePolynomialZp cVal = ZippelGCD(aVal, bVal, rnd, variable - 1, tmpDegreeBounds, evaluationStackLimit);
            if (cVal == null)
                //unlucky homomorphism
                continue;

            int currExponent = cVal.degree(variable - 1);
            if (currExponent > tmpDegreeBounds[variable - 1])
                //unlucky homomorphism
                continue;

            if (currExponent < tmpDegreeBounds[variable - 1]) {
                //better degree bound detected
                tmpDegreeBounds[variable - 1] = currExponent;
            }

            cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
            assert cVal.lc() == lcVal;

            lSparseInterpolation sparseInterpolator = createInterpolation(variable, a, b, cVal, rnd);
            if (sparseInterpolator == null)
                //unlucky homomorphism
                continue;

            // we are applying dense interpolation for univariate skeleton coefficients
            lInterpolation denseInterpolation = new lInterpolation(variable, seedPoint, cVal);
            //previous interpolation (used to detect whether update doesn't change the result)
            lMultivariatePolynomialZp previousInterpolation;
            //local evaluation stack for points that are calculated via sparse interpolation (but not gcd evaluation) -> always same skeleton
            TLongHashSet localEvaluationStack = new TLongHashSet(globalEvaluationStack);
            while (true) {
                if (evaluationStackLimit == localEvaluationStack.size())
                    return null;

                if (denseInterpolation.numberOfPoints() > tmpDegreeBounds[variable] + ALLOWED_OVER_INTERPOLATED_ATTEMPTS) {
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                long randomPoint = domain.randomElement(rnd);
                if (localEvaluationStack.contains(randomPoint))
                    continue;
                localEvaluationStack.add(randomPoint);

                lcVal = lcGCD.evaluate(randomPoint);
                if (lcVal == 0)
                    continue;

                cVal = sparseInterpolator.evaluate(randomPoint);
                if (cVal == null) {
                    ++failedSparseInterpolations;
                    if (failedSparseInterpolations == MAX_SPARSE_INTERPOLATION_FAILS)
                        throw new RuntimeException("Sparse interpolation failed");
                    // restore original degree bounds, since unlucky homomorphism may destruct correct bounds
                    tmpDegreeBounds = degreeBounds.clone();
                    continue main;
                }
                cVal = cVal.multiply(domain.multiply(domain.reciprocal(cVal.lc()), lcVal));
                assert cVal.lc() == lcVal;

                // Cache previous interpolation. NOTE: clone() is important, since the poly will
                // be modified inplace by the update() method
                previousInterpolation = denseInterpolation.getInterpolatingPolynomial().clone();
                denseInterpolation.update(randomPoint, cVal);

                // do division test
                if (tmpDegreeBounds[variable] <= denseInterpolation.numberOfPoints()
                        || previousInterpolation.equals(denseInterpolation.getInterpolatingPolynomial())) {
                    lMultivariatePolynomialZp result = doDivisionCheck(a, b, contentGCD, denseInterpolation, variable);
                    if (result != null)
                        return result;
                }
            }
        }
    }

    static lSparseInterpolation createInterpolation(int variable,
                                                    lMultivariatePolynomialZp a,
                                                    lMultivariatePolynomialZp b,
                                                    lMultivariatePolynomialZp skeleton,
                                                    RandomGenerator rnd) {
        skeleton = skeleton.clone().setAllCoefficientsToUnit();
        if (skeleton.size() == 1)
            return new lTrivialSparseInterpolation(skeleton);

        boolean monic = a.coefficientOf(0, a.degree(0)).isConstant() && b.coefficientOf(0, a.degree(0)).isConstant();

        Set<DegreeVector> globalSkeleton = skeleton.getSkeleton();
        TIntObjectHashMap<lMultivariatePolynomialZp> univarSkeleton = getSkeleton(skeleton);
        int[] sparseUnivarDegrees = univarSkeleton.keys();

        lIntegersModulo domain = a.domain;

        int lastVariable = variable == -1 ? a.nVariables - 1 : variable;
        int[] evaluationVariables = ArraysUtil.sequence(1, lastVariable + 1);//variable inclusive
        long[] evaluationPoint = new long[evaluationVariables.length];

        lPrecomputedPowersHolder powers;
        int fails = 0;
        search_for_good_evaluation_point:
        while (true) {
            if (fails >= MAX_FAILED_SUBSTITUTIONS)
                return null;
            //avoid zero evaluation points
            for (int i = lastVariable - 1; i >= 0; --i)
                do {
                    evaluationPoint[i] = domain.randomElement(rnd);
                } while (evaluationPoint[i] == 0);

            powers = new lPrecomputedPowersHolder(a.nVariables, evaluationVariables, evaluationPoint, domain);
            int[] raiseFactors = ArraysUtil.arrayOf(1, evaluationVariables.length);

            for (lMultivariatePolynomialZp p : Arrays.asList(a, b, skeleton))
                if (!p.getSkeleton(0).equals(p.evaluate(powers, evaluationVariables, raiseFactors).getSkeleton())) {
                    ++fails;
                    continue search_for_good_evaluation_point;
                }
            break;
        }

        int requiredNumberOfEvaluations = -1, monicScalingExponent = -1;
        for (TIntObjectIterator<lMultivariatePolynomialZp> it = univarSkeleton.iterator(); it.hasNext(); ) {
            it.advance();
            lMultivariatePolynomialZp v = it.value();
            if (v.size() > requiredNumberOfEvaluations)
                requiredNumberOfEvaluations = v.size();
            if (v.size() == 1)
                monicScalingExponent = it.key();
        }

        if (!ALWAYS_LINZIP) {
            if (monic)
                monicScalingExponent = -1;

            if (monic || monicScalingExponent != -1)
                return new lMonicInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                        evaluationVariables, evaluationPoint, powers, rnd, requiredNumberOfEvaluations, monicScalingExponent);
        }

        return new lLinZipInterpolation(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees,
                evaluationVariables, evaluationPoint, powers, rnd);
    }

    /**
     * view multivariate polynomial as a univariate in Zp[x_1, ... x_N][x_0] and
     * return the map (x_0)^exponent -> coefficient in Zp[x_1, ... x_N]
     */
    private static TIntObjectHashMap<lMultivariatePolynomialZp> getSkeleton(lMultivariatePolynomialZp poly) {
        TIntObjectHashMap<lMultivariatePolynomialZp> skeleton = new TIntObjectHashMap<>();
        for (lMonomialTerm term : poly.terms) {
            lMonomialTerm newDV = term.setZero(0);
            lMultivariatePolynomialZp coeff = skeleton.get(term.exponents[0]);
            if (coeff != null)
                coeff.add(newDV);
            else
                skeleton.put(term.exponents[0], lMultivariatePolynomialZp.create(poly.nVariables, poly.domain, poly.ordering, newDV));
        }
        return skeleton;
    }

    interface lSparseInterpolation {
        /**
         * Returns interpolating gcd
         *
         * @return gcd
         */
        lMultivariatePolynomialZp evaluate();

        /**
         * Returns interpolating gcd at specified point
         *
         * @param newPoint evaluation point
         * @return gcd at {@code variable = newPoint}
         */
        lMultivariatePolynomialZp evaluate(long newPoint);
    }

    static final class lTrivialSparseInterpolation implements lSparseInterpolation {
        final lMultivariatePolynomialZp val;

        lTrivialSparseInterpolation(lMultivariatePolynomialZp val) {
            this.val = val;
        }

        @Override
        public lMultivariatePolynomialZp evaluate() {
            return val;
        }

        @Override
        public lMultivariatePolynomialZp evaluate(long newPoint) {
            return val;
        }
    }

    static abstract class lASparseInterpolation implements lSparseInterpolation {
        /** the domain */
        final lIntegersModulo domain;
        /** variable which we are evaluating */
        final int variable;
        /** initial polynomials */
        final lMultivariatePolynomialZp a, b;
        /** global skeleton of the result */
        final Set<DegreeVector> globalSkeleton;
        /** skeleton of each of the coefficients of polynomial viewed as Zp[x_1,...,x_N][x_0] */
        final TIntObjectHashMap<lMultivariatePolynomialZp> univarSkeleton;
        /** univariate degrees of {@code univarSkeleton} with respect to x_0 */
        final int[] sparseUnivarDegrees;
        /** variables that will be substituted with random values for sparse interpolation, i.e. {@code [1, 2 ... variable] } */
        final int[] evaluationVariables;
        /**
         * values that will be subsituted for {@code evaluationVariables};
         * values {@code V} for variables {@code [1, 2 ... variable-1]} are fixed and successive
         * powers {@code V, V^2, V^3 ...} will be used to form Vandermonde matrix
         */
        final long[] evaluationPoint;
        /** cached powers for {@code evaluationPoint} */
        final lPrecomputedPowersHolder powers;
        /** random */
        final RandomGenerator rnd;

        lASparseInterpolation(lIntegersModulo domain, int variable,
                              lMultivariatePolynomialZp a, lMultivariatePolynomialZp b,
                              Set<DegreeVector> globalSkeleton,
                              TIntObjectHashMap<lMultivariatePolynomialZp> univarSkeleton,
                              int[] sparseUnivarDegrees, int[] evaluationVariables,
                              long[] evaluationPoint,
                              lPrecomputedPowersHolder powers, RandomGenerator rnd) {
            this.domain = domain;
            this.variable = variable;
            this.a = a;
            this.b = b;
            this.globalSkeleton = globalSkeleton;
            this.univarSkeleton = univarSkeleton;
            this.sparseUnivarDegrees = sparseUnivarDegrees;
            this.evaluationVariables = evaluationVariables;
            this.evaluationPoint = evaluationPoint;
            this.powers = powers;
            this.rnd = rnd;
        }

        @Override
        public final lMultivariatePolynomialZp evaluate(long newPoint) {
            // constant is constant
            if (globalSkeleton.size() == 1)
                return a.create(((lMonomialTerm) globalSkeleton.iterator().next()).setCoefficient(1));
            // variable = newPoint
            evaluationPoint[evaluationPoint.length - 1] = newPoint;
            powers.set(evaluationVariables[evaluationVariables.length - 1], newPoint);
            return evaluate();
        }
    }

    static final class lLinZipInterpolation extends lASparseInterpolation {
        lLinZipInterpolation(lIntegersModulo domain, int variable, lMultivariatePolynomialZp a,
                             lMultivariatePolynomialZp b, Set<DegreeVector> globalSkeleton,
                             TIntObjectHashMap<lMultivariatePolynomialZp> univarSkeleton, int[] sparseUnivarDegrees,
                             int[] evaluationVariables, long[] evaluationPoint, lPrecomputedPowersHolder powers,
                             RandomGenerator rnd) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
        }

        @Override
        public lMultivariatePolynomialZp evaluate() {
            @SuppressWarnings("unchecked")
            lLinZipSystem[] systems = new lLinZipSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new lLinZipSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? a.nVariables - 1 : variable - 1);


            int[] raiseFactors = new int[evaluationVariables.length];
            if (variable != -1)
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;

            int nUnknowns = globalSkeleton.size(), nUnknownScalings = -1;
            int raiseFactor = 0;

            for (int nTries = 0; nTries < NUMBER_OF_UNDER_DETERMINED_RETRIES; ++nTries) {
                int previousFreeVars = -1, underDeterminedTries = 0;
                for (; ; ) {
                    // increment at each loop!
                    ++nUnknownScalings;
                    ++raiseFactor;
                    // sequential powers of evaluation point
                    Arrays.fill(raiseFactors, 0, variable == -1 ? raiseFactors.length : raiseFactors.length - 1, raiseFactor);
                    // evaluate a and b to univariate and calculate gcd
                    lUnivariatePolynomialZp
                            aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                            gcdUnivar = UnivariateGCD.PolynomialGCD(aUnivar, bUnivar);

                    if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                        // unlucky main homomorphism or bad evaluation point
                        return null;

                    assert gcdUnivar.isMonic();
                    if (!univarSkeleton.keySet().containsAll(gcdUnivar.exponents()))
                        // univariate gcd contains terms that are not present in the skeleton
                        // again unlucky main homomorphism
                        return null;

                    int totalEquations = 0;
                    for (lLinZipSystem system : systems) {
                        long rhs = gcdUnivar.degree() < system.univarDegree ? 0 : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs, nUnknownScalings != 0);
                        totalEquations += system.nEquations();
                    }
                    if (nUnknowns + nUnknownScalings <= totalEquations)
                        break;


                    if (underDeterminedTries > NUMBER_OF_UNDER_DETERMINED_RETRIES) {
                        // raise condition: new equations does not fix enough variables
                        return null;
                    }

                    int freeVars = nUnknowns + nUnknownScalings - totalEquations;
                    if (freeVars >= previousFreeVars)
                        ++underDeterminedTries;
                    else
                        underDeterminedTries = 0;

                    previousFreeVars = freeVars;
                }

                lMultivariatePolynomialZp result = a.createZero();
                SystemInfo info = solveLinZip(a, systems, nUnknownScalings, result);
                if (info == SystemInfo.UnderDetermined) {
                    //try to generate more equations
                    continue;
                }

                if (info == SystemInfo.Consistent)
                    //well done
                    return result;
                if (info == SystemInfo.Inconsistent)
                    //inconsistent system => unlucky homomorphism
                    return null;
            }

            // the system is still under-determined => bad evaluation homomorphism
            return null;
        }
    }

    private static SystemInfo solveLinZip(lMultivariatePolynomialZp factory, lLinZipSystem[] subSystems, int nUnknownScalings, lMultivariatePolynomialZp destination) {
        ArrayList<lMonomialTerm> unknowns = new ArrayList<>();
        for (lLinZipSystem system : subSystems)
            for (lMonomialTerm degreeVector : system.skeleton)
                unknowns.add(degreeVector.set(0, system.univarDegree));

        int nUnknownsMonomials = unknowns.size();
        int nUnknownsTotal = nUnknownsMonomials + nUnknownScalings;
        ArrayList<long[]> lhsGlobal = new ArrayList<>();
        TLongArrayList rhsGlobal = new TLongArrayList();
        int offset = 0;
        lIntegersModulo domain = factory.domain;
        for (lLinZipSystem system : subSystems) {
            for (int j = 0; j < system.matrix.size(); j++) {
                long[] row = new long[nUnknownsTotal];
                long[] subRow = system.matrix.get(j);

                System.arraycopy(subRow, 0, row, offset, subRow.length);
                if (j > 0)
                    row[nUnknownsMonomials + j - 1] = system.scalingMatrix.get(j);
                lhsGlobal.add(row);
                rhsGlobal.add(system.rhs.get(j));
            }

            offset += system.skeleton.length;
        }

        long[] solution = new long[nUnknownsTotal];
        SystemInfo info = LinearAlgebra.solve(domain, lhsGlobal, rhsGlobal, solution);
        if (info == SystemInfo.Consistent) {
            @SuppressWarnings("unchecked")
            lMonomialTerm[] terms = new lMonomialTerm[unknowns.size()];
            for (int i = 0; i < terms.length; i++)
                terms[i] = unknowns.get(i).setCoefficient(solution[i]);
            destination.add(terms);
        }
        return info;
    }


    static final class lMonicInterpolation extends lASparseInterpolation {
        /** required number of sparse evaluations to reconstruct all coefficients in skeleton */
        final int requiredNumberOfEvaluations;
        /** univar exponent with monomial factor that can be used for scaling */
        final int monicScalingExponent;

        lMonicInterpolation(lIntegersModulo domain, int variable, lMultivariatePolynomialZp a, lMultivariatePolynomialZp b,
                            Set<DegreeVector> globalSkeleton, TIntObjectHashMap<lMultivariatePolynomialZp> univarSkeleton,
                            int[] sparseUnivarDegrees, int[] evaluationVariables, long[] evaluationPoint,
                            lPrecomputedPowersHolder powers, RandomGenerator rnd, int requiredNumberOfEvaluations,
                            int monicScalingExponent) {
            super(domain, variable, a, b, globalSkeleton, univarSkeleton, sparseUnivarDegrees, evaluationVariables, evaluationPoint, powers, rnd);
            this.requiredNumberOfEvaluations = requiredNumberOfEvaluations;
            this.monicScalingExponent = monicScalingExponent;
        }

        @Override
        public lMultivariatePolynomialZp evaluate() {
            @SuppressWarnings("unchecked")
            lVandermondeSystem[] systems = new lVandermondeSystem[sparseUnivarDegrees.length];
            for (int i = 0; i < sparseUnivarDegrees.length; i++)
                systems[i] = new lVandermondeSystem(sparseUnivarDegrees[i], univarSkeleton.get(sparseUnivarDegrees[i]), powers, variable == -1 ? a.nVariables - 1 : variable - 1);

            int[] raiseFactors = new int[evaluationVariables.length];
            if (variable != -1)
                // the last variable (the variable) is the same for all evaluations = newPoint
                raiseFactors[raiseFactors.length - 1] = 1;
            for (int i = 0; i < requiredNumberOfEvaluations; ++i) {
                // sequential powers of evaluation point
                Arrays.fill(raiseFactors, 0, variable == -1 ? raiseFactors.length : raiseFactors.length - 1, i + 1);
                // evaluate a and b to univariate and calculate gcd
                lUnivariatePolynomialZp
                        aUnivar = a.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        bUnivar = b.evaluate(powers, evaluationVariables, raiseFactors).asUnivariate(),
                        gcdUnivar = UnivariateGCD.PolynomialGCD(aUnivar, bUnivar);

                if (a.degree(0) != aUnivar.degree() || b.degree(0) != bUnivar.degree())
                    // unlucky main homomorphism or bad evaluation point
                    return null;

                assert gcdUnivar.isMonic();
                if (!univarSkeleton.keySet().containsAll(gcdUnivar.exponents()))
                    // univariate gcd contains terms that are not present in the skeleton
                    // again unlucky main homomorphism
                    return null;

                if (monicScalingExponent != -1) {
                    // single scaling factor
                    // scale the system according to it

                    if (gcdUnivar.degree() < monicScalingExponent || gcdUnivar.get(monicScalingExponent) == 0)
                        // unlucky homomorphism
                        return null;

                    long normalization = evaluateExceptFirst(domain, powers, 1, univarSkeleton.get(monicScalingExponent).lt(), i + 1, variable == -1 ? a.nVariables - 1 : variable - 1);
                    //normalize univariate gcd in order to reconstruct leading coefficient polynomial
                    normalization = domain.multiply(domain.reciprocal(gcdUnivar.get(monicScalingExponent)), normalization);
                    gcdUnivar = gcdUnivar.multiply(normalization);
                }

                boolean allDone = true;
                for (lVandermondeSystem system : systems)
                    if (system.nEquations() < system.nUnknownVariables()) {
                        long rhs = gcdUnivar.degree() < system.univarDegree ? 0 : gcdUnivar.get(system.univarDegree);
                        system.oneMoreEquation(rhs);
                        allDone = false;
                    }

                if (allDone)
                    break;
            }

            for (lVandermondeSystem system : systems) {
                //solve each system
                SystemInfo info = system.solve();
                if (info != SystemInfo.Consistent)
                    // system is inconsistent or under determined
                    // unlucky homomorphism
                    return null;
            }

            lMultivariatePolynomialZp gcdVal = a.createZero();
            for (lVandermondeSystem system : systems) {
                assert monicScalingExponent == -1 || system.univarDegree != monicScalingExponent || system.solution[0] == 1;
                for (int i = 0; i < system.skeleton.length; i++) {
                    lMonomialTerm degreeVector = system.skeleton[i].set(0, system.univarDegree);
                    long value = system.solution[i];
                    gcdVal.add(degreeVector.setCoefficient(value));
                }
            }

            return gcdVal;
        }
    }

    private static abstract class lLinearSystem {
        final int univarDegree;
        /** the domain */
        final lIntegersModulo domain;
        /** the skeleton */
        final lMonomialTerm[] skeleton;
        /** the lhs matrix */
        final ArrayList<long[]> matrix;
        /** the rhs values */
        final TLongArrayList rhs = new TLongArrayList();
        /** precomputed powers */
        final lPrecomputedPowersHolder powers;
        /** number of non-fixed variables, i.e. variables that will be substituted */
        final int nVars;

        @SuppressWarnings("unchecked")
        lLinearSystem(int univarDegree, lMultivariatePolynomialZp skeleton, lPrecomputedPowersHolder powers, int nVars) {
            this.univarDegree = univarDegree;
            this.domain = skeleton.domain;
            //todo refactor generics
            this.skeleton = skeleton.getSkeleton().toArray(new lMonomialTerm[skeleton.size()]);
            this.powers = powers;
            this.nVars = nVars;
            this.matrix = new ArrayList<>();
        }

        final int nUnknownVariables() {
            return skeleton.length;
        }

        final int nEquations() {
            return matrix.size();
        }

        @Override
        public String toString() {
            return "{" + matrix.stream().map(Arrays::toString).collect(Collectors.joining(",")) + "} = " + rhs;
        }
    }

    /** Vandermonde system builder */
    private static final class lLinZipSystem extends lLinearSystem {
        public lLinZipSystem(int univarDegree, lMultivariatePolynomialZp skeleton, lPrecomputedPowersHolder powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        private final TLongArrayList scalingMatrix = new TLongArrayList();

        public void oneMoreEquation(long rhsVal, boolean newScalingIntroduced) {
            long[] row = new long[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, 1, skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);

            if (newScalingIntroduced) {
                scalingMatrix.add(domain.negate(rhsVal));
                rhsVal = 0;
            } else
                scalingMatrix.add(0);
            rhs.add(rhsVal);
        }
    }

    /** Vandermonde system builder */
    private static final class lVandermondeSystem extends lLinearSystem {
        public lVandermondeSystem(int univarDegree, lMultivariatePolynomialZp skeleton, lPrecomputedPowersHolder powers, int nVars) {
            super(univarDegree, skeleton, powers, nVars);
        }

        long[] solution = null;

        SystemInfo solve() {
            if (solution == null)
                solution = new long[nUnknownVariables()];

            if (nUnknownVariables() <= 8)
                // for small systems Gaussian elimination is indeed faster
                return LinearAlgebra.solve(domain, matrix.toArray(new long[matrix.size()][]), rhs.toArray(), solution);

            // solve vandermonde system
            long[] vandermondeRow = matrix.get(0);
            SystemInfo info = LinearAlgebra.solveVandermondeT(domain, vandermondeRow, rhs.toArray(), solution);
            if (info == SystemInfo.Consistent)
                for (int i = 0; i < solution.length; ++i)
                    solution[i] = domain.divide(solution[i], vandermondeRow[i]);

            return info;
        }

        public lVandermondeSystem oneMoreEquation(long rhsVal) {
            long[] row = new long[skeleton.length];
            for (int i = 0; i < skeleton.length; i++)
                row[i] = evaluateExceptFirst(domain, powers, 1, skeleton[i], matrix.size() + 1, nVars);
            matrix.add(row);
            rhs.add(rhsVal);
            return this;
        }
    }

    private static long evaluateExceptFirst(lIntegersModulo domain,
                                            lPrecomputedPowersHolder powers,
                                            long coefficient,
                                            lMonomialTerm skeleton,
                                            int raiseFactor,
                                            int nVars) {
        long tmp = coefficient;
        for (int k = 1; k <= nVars; k++)
            tmp = domain.multiply(tmp, powers.pow(k, raiseFactor * skeleton.exponents[k]));
        return tmp;
    }

    /* =============================================== EZ-GCD algorithm ============================================ */


    /**
     * Calculates GCD of two multivariate polynomials over Zp using EZ algorithm
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static lMultivariatePolynomialZp EZGCD(
            lMultivariatePolynomialZp a,
            lMultivariatePolynomialZp b) {

        // prepare input and test for early termination
        GCDInput<lMonomialTerm, lMultivariatePolynomialZp> gcdInput = preparedGCDInput(a, b, MultivariateGCD::EZGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

        // remove content in each variable
        lMultivariatePolynomialZp content = a.createOne();
        for (int i = 0; i < a.nVariables; ++i) {
            if (a.degree(i) == 0 || b.degree(i) == 0)
                continue;
            lMultivariatePolynomialZp tmpContent = contentGCD(a, b, i, MultivariateGCD::EZGCD);
            a = divideExact(a, tmpContent);
            b = divideExact(b, tmpContent);
            content = content.multiply(tmpContent);
        }

        // one more reduction; removing of content may shuffle required variables order
        GCDInput<lMonomialTerm, lMultivariatePolynomialZp> gcdInput2 = preparedGCDInput(a, b, MultivariateGCD::EZGCD);
        a = gcdInput2.aReduced; b = gcdInput2.bReduced;

        lMultivariatePolynomialZp result = gcdInput2.earlyGCD != null
                ? gcdInput2.earlyGCD
                : gcdInput2.restoreGCD(EZGCD0(a, b, gcdInput2.lastPresentVariable + 1, PrivateRandom.getRandom()));

        result = result.multiply(content);
        return gcdInput.restoreGCD(result);
    }

    /** actual EZ-GCD implementation */
    private static lMultivariatePolynomialZp EZGCD0(
            lMultivariatePolynomialZp a, lMultivariatePolynomialZp b, int nUsedVariables, RandomGenerator rnd) {

        // degree of univariate gcd
        int ugcdDegree = Integer.MAX_VALUE;

        EZGCDEvaluations evaluations = new EZGCDEvaluations(a, b, nUsedVariables, rnd);
        choose_evaluation:
        while (true) {
            // set new evaluation point (true returned if base variable has changed)
            boolean swapped = evaluations.nextEvaluation();
            if (swapped)
                ugcdDegree = Integer.MAX_VALUE;
            a = evaluations.aReduced;
            b = evaluations.bReduced;

            // degrees of a and b as Z[y1, ... ,yN][x]
            int
                    uaDegree = a.degree(0),
                    ubDegree = b.degree(0);

            lUnivariatePolynomialZp
                    ua = a.evaluateAtZeroAllExcept(0),
                    ub = b.evaluateAtZeroAllExcept(0);

            assert ua.degree() == uaDegree;
            assert ub.degree() == ubDegree;

            // gcd of a mod I and b mod I (univariate)
            lUnivariatePolynomialZp ugcd = UnivariateGCD.PolynomialGCD(ua, ub);

            if (ugcd.degree() == 0) {
                // coprime polynomials
                return evaluations.reconstruct(a.createOne());
            }

            if (ugcd.degree() > ugcdDegree)
                // unlucky evaluation
                continue choose_evaluation;

            ugcdDegree = ugcd.degree();

            if (ugcdDegree == uaDegree) {
                // a is a divisor of b
                if (dividesQ(b, a))
                    return evaluations.reconstruct(a);
                //continue choose_evaluation;
            }

            if (ugcdDegree == ubDegree) {
                // b is a divisor of a
                if (dividesQ(a, b))
                    return evaluations.reconstruct(b);
                //continue choose_evaluation;
            }

            // base polynomial to lift (either a or b)
            lMultivariatePolynomialZp base = null;
            // a cofactor with gcd to lift, i.e. base = ugcd * uCoFactor mod I
            lUnivariatePolynomialZp uCoFactor = null;
            lMultivariatePolynomialZp coFactorLC = null;

            // deg(b) < deg(a), so it is better to try to lift b
            lUnivariatePolynomialZp ubCoFactor = DivisionWithRemainder.quotient(ub, ugcd, true);
            if (UnivariateGCD.PolynomialGCD(ugcd, ubCoFactor).isConstant()) {
                // b splits into coprime factors
                base = b;
                uCoFactor = ubCoFactor;
                coFactorLC = b.lc(0);
            }

            if (base == null) {
                // b does not split into coprime factors => try a
                lUnivariatePolynomialZp uaCoFactor = DivisionWithRemainder.quotient(ua, ugcd, true);
                if (UnivariateGCD.PolynomialGCD(ugcd, uaCoFactor).isConstant()) {
                    base = a;
                    uCoFactor = uaCoFactor;
                    coFactorLC = a.lc(0);
                }
            }

            if (base == null) {
                // neither a nor b does not split into coprime factors => square free decomposition required

                lMultivariatePolynomialZp
                        bRepeatedFactor = EZGCD(b, b.derivative(0, 1)),
                        squareFreeGCD = EZGCD(divideExact(b, bRepeatedFactor), a),
                        aRepeatedFactor = divideExact(a, squareFreeGCD),
                        gcd = squareFreeGCD.clone();

                lUnivariatePolynomialZp
                        uSquareFreeGCD = squareFreeGCD.evaluateAtZeroAllExcept(0),
                        uaRepeatedFactor = aRepeatedFactor.evaluateAtZeroAllExcept(0),
                        ubRepeatedFactor = bRepeatedFactor.evaluateAtZeroAllExcept(0);
                while (true) {
                    ugcd = UnivariateGCD.PolynomialGCD(
                            uSquareFreeGCD, uaRepeatedFactor, ubRepeatedFactor);

                    if (ugcd.degree() == 0)
                        return evaluations.reconstruct(gcd);

                    if (!uSquareFreeGCD.clone().monic().equals(ugcd.clone().monic())) {
                        lMultivariatePolynomialZp
                                mgcd = lMultivariatePolynomialZp.asMultivariate(ugcd, a.nVariables, 0, a.ordering),
                                gcdCoFactor = lMultivariatePolynomialZp.asMultivariate(
                                        DivisionWithRemainder.divideExact(uSquareFreeGCD, ugcd, false), a.nVariables, 0, a.ordering);

                        HenselLifting.liftPair(squareFreeGCD, mgcd, gcdCoFactor);

                        squareFreeGCD = mgcd.clone();
                        uSquareFreeGCD = squareFreeGCD.evaluateAtZeroAllExcept(0);
                    }

                    gcd = gcd.multiply(squareFreeGCD);
                    uaRepeatedFactor = DivisionWithRemainder.divideExact(uaRepeatedFactor, ugcd, false);
                    ubRepeatedFactor = DivisionWithRemainder.divideExact(ubRepeatedFactor, ugcd, false);
                }
            }

            lMultivariatePolynomialZp gcd = lMultivariatePolynomialZp.asMultivariate(ugcd, a.nVariables, 0, a.ordering);
            lMultivariatePolynomialZp coFactor = lMultivariatePolynomialZp.asMultivariate(uCoFactor, a.nVariables, 0, a.ordering);

            // impose the leading coefficient
            lMultivariatePolynomialZp lcCorrection = EZGCD(a.lc(0), b.lc(0));
            assert ZippelGCD(a.lc(0), b.lc(0)).monic(lcCorrection.lc()).equals(lcCorrection) : "\n" + a.lc(0) + "  \n " + b.lc(0);

            if (lcCorrection.isOne()) {
                HenselLifting.liftPair0(base, gcd, coFactor, null, null);
            } else {
                long lcCorrectionMod = lcCorrection.cc(); // substitute all zeros
                assert lcCorrectionMod != 0;

                coFactor = coFactor.multiply(gcd.lc());
                gcd = gcd.monic(lcCorrectionMod);

                HenselLifting.liftPair0(base.clone().multiply(lcCorrection), gcd, coFactor, lcCorrection, coFactorLC);
                assert gcd.lc(0).equals(lcCorrection);
                gcd = HenselLifting.primitivePart(gcd);
            }

            if (dividesQ(b, gcd) && dividesQ(a, gcd))
                return evaluations.reconstruct(gcd);
        }
    }

    static ZeroVariables commonPossibleZeroes(lMultivariatePolynomialZp a, lMultivariatePolynomialZp b, int nVariables) {
        return commonPossibleZeroes(possibleZeros(a, nVariables), possibleZeros(b, nVariables));
    }

    static ZeroVariables commonPossibleZeroes(ZeroVariables a, ZeroVariables b) {
        if (a.allZeroes)
            return b;
        if (b.allZeroes)
            return a;

        ZeroVariables result = new ZeroVariables(a.nVariables);
        for (BitSet az : a.pZeros)
            for (BitSet bz : b.pZeros) {
                BitSet tmp = (BitSet) az.clone();
                tmp.and(bz);
                result.add(tmp);
            }

        return result;
    }

    /** list of available zero substitutions for poly */
    static ZeroVariables possibleZeros(lMultivariatePolynomialZp poly, int nVariables) {
        ZeroVariables result = new ZeroVariables(nVariables);
        if (poly.cc() != 0) {
            BitSet s = new BitSet(nVariables);
            s.set(0, nVariables);
            result.add(s);
            return result;
        }
        for (lMonomialTerm term : poly.terms) {
            BitSet zeroes = new BitSet(nVariables);
            for (int i = 0; i < nVariables; i++)
                if (term.exponents[i] == 0)
                    zeroes.set(i);
            result.add(zeroes);
        }
        return result;
    }

    /** a list of variables that can be replaced with zeroes */
    static final class ZeroVariables {
        final List<BitSet> pZeros = new ArrayList<>();
        final int nVariables;

        public ZeroVariables(int nVariables) {
            this.nVariables = nVariables;
        }

        int iMaxPossibleZeros = -1, maxPossibleZeros = -1;
        boolean allZeroes = false;

        void add(BitSet zeros) {
            if (allZeroes)
                return;
            int nZeros = zeros.cardinality();
            if (nVariables == nZeros) {
                maxPossibleZeros = zeros.length();
                iMaxPossibleZeros = 0;
                allZeroes = true;
                pZeros.clear();
                pZeros.add(zeros);
                return;
            }
            Iterator<BitSet> it = pZeros.iterator();
            while (it.hasNext()) {
                BitSet e = it.next();
                if (contains(e, zeros))
                    return;
                if (contains(zeros, e))
                    it.remove();
            }
            pZeros.add(zeros);
            if (nZeros > maxPossibleZeros) {
                maxPossibleZeros = nZeros;
                iMaxPossibleZeros = pZeros.size() - 1;
            }
        }

        BitSet maxZeros() {
            return pZeros.isEmpty() ? new BitSet(nVariables) : pZeros.get(iMaxPossibleZeros);
        }

        @Override
        public String toString() {
            return pZeros.toString();
        }
    }

    /** test whether major contains all elements of minor */
    static boolean contains(BitSet major, BitSet minor) {
        if (major.equals(minor))
            return true;
        BitSet tmp = (BitSet) major.clone();
        tmp.and(minor);
        return tmp.equals(minor);
    }

    private static final int
            SAME_BASE_VARIABLE_ATTEMPTS = 8,
            ATTEMPTS_WITH_MAX_ZEROS = 2;

    static final class EZGCDEvaluations {
        final lMultivariatePolynomialZp a, b;
        final RandomGenerator rnd;
        // actual number of variables in use [0, ... nVariables] (may be less than a.nVariables)
        final int nVariables;
        // variables in which both cc and lc in a and b are non zeroes
        final TIntArrayList perfectVariables = new TIntArrayList();

        EZGCDEvaluations(lMultivariatePolynomialZp a,
                         lMultivariatePolynomialZp b,
                         int nVariables,
                         RandomGenerator rnd) {
            this.a = a;
            this.b = b;
            this.rnd = rnd;
            this.nVariables = nVariables;

            // search for main variables
            for (int var = 0; var < nVariables; ++var) {
                lUnivariatePolynomialZp ub = b.evaluateAtZeroAllExcept(var);
                if (ub.degree() == 0 || ub.cc() == 0 || ub.degree() != b.degree(var))
                    continue;

                lUnivariatePolynomialZp ua = a.evaluateAtZeroAllExcept(var);
                if (ua.degree() == 0 || ua.cc() == 0 || ua.degree() != a.degree(var))
                    continue;

                perfectVariables.add(var);
            }
        }

        // modified a and b so that all except first variables
        // may be evaluated to zero
        lMultivariatePolynomialZp aReduced, bReduced;

        // pointers in perfectVariables list
        int
                currentPerfectVariable = -1,
                nextPerfectVariable = 0;

        int baseVariable = -1;
        int[]
                zeroVariables, // variables that may be replaced with zeros ( 0 <> baseVariable are swaped )
                shiftingVariables; // variables that must be shifted before substituting zeros ( 0 <> baseVariable are swaped )
        // shift values for shiftingVariables
        long[] shifts;

        // number of tries keeping as many zero variables as possible
        int nAttemptsWithZeros = 0;
        // variable that were used as main variables
        final TIntHashSet usedBaseVariables = new TIntHashSet();
        int nAttemptsWithSameBaseVariable = 0;

        // returns whether the main variable was changed
        boolean nextEvaluation() {
            // first try to evaluate all but one variables to zeroes

            if (nextPerfectVariable < perfectVariables.size()) {
                // we have the main variable
                currentPerfectVariable = perfectVariables.get(nextPerfectVariable);
                ++nextPerfectVariable;

                aReduced = AMultivariatePolynomial.swapVariables(a, 0, currentPerfectVariable);
                bReduced = AMultivariatePolynomial.swapVariables(b, 0, currentPerfectVariable);
                return true;
            }

            // no possible zero substitutions more
            // => shift some variables and try to
            // minimize intermediate expression swell
            // due to "bad zero problem"

            boolean changedMainVariable = false;
            currentPerfectVariable = -1;

            int previousBaseVariable = baseVariable;
            // try to substitute as many zeroes as possible
            if (baseVariable == -1 || nAttemptsWithSameBaseVariable > SAME_BASE_VARIABLE_ATTEMPTS) {
                if (usedBaseVariables.size() == nVariables)
                    usedBaseVariables.clear();
                nAttemptsWithSameBaseVariable = nAttemptsWithZeros = 0;
                ZeroVariables maxZeros = null;
                for (int var = 0; var < nVariables; ++var) {
                    if (usedBaseVariables.contains(var))
                        continue;
                    ZeroVariables pZeros;
                    if (a.degree(var) == 0 || b.degree(var) == 0)
                        continue;
                    pZeros = possibleZeros(b.coefficientOf(var, 0), nVariables);
                    if (maxZeros != null && pZeros.maxPossibleZeros < maxZeros.maxPossibleZeros)
                        continue;

                    pZeros = commonPossibleZeroes(pZeros, possibleZeros(b.coefficientOf(var, b.degree(var)), nVariables));
                    if (maxZeros != null && pZeros.maxPossibleZeros < maxZeros.maxPossibleZeros)
                        continue;

                    pZeros = commonPossibleZeroes(pZeros, possibleZeros(a.coefficientOf(var, 0), nVariables));
                    if (maxZeros != null && pZeros.maxPossibleZeros < maxZeros.maxPossibleZeros)
                        continue;

                    pZeros = commonPossibleZeroes(pZeros, possibleZeros(a.coefficientOf(var, a.degree(var)), nVariables));
                    if (maxZeros != null && pZeros.maxPossibleZeros < maxZeros.maxPossibleZeros)
                        continue;

                    maxZeros = pZeros;
                    baseVariable = var;
                }
                if (maxZeros == null) {
                    // all free variables are bad
                    usedBaseVariables.clear();
                    return nextEvaluation();
                }
                usedBaseVariables.add(baseVariable);

                // <- we chose the main variable and those that can be safely substituted with zeros

                BitSet zeroSubstitutions = maxZeros.maxZeros();
                // swap 0 and baseVariable
                zeroSubstitutions.set(baseVariable, zeroSubstitutions.get(0));
                zeroSubstitutions.clear(0);

                zeroVariables = new int[zeroSubstitutions.cardinality()];
                shiftingVariables = new int[nVariables - zeroSubstitutions.cardinality() - 1];
                int iCounter = 0, jCounter = 0;
                for (int j = 0; j < nVariables; j++)
                    if (zeroSubstitutions.get(j))
                        zeroVariables[iCounter++] = j;
                    else if (j != 0) // don't shift the main variable
                        shiftingVariables[jCounter++] = j;

                if (shiftingVariables.length == 0) { // <- perfect variable exists; just pickup some var for random substitutions
                    shiftingVariables = new int[]{zeroVariables[zeroVariables.length - 1]};
                    zeroVariables = Arrays.copyOf(zeroVariables, zeroVariables.length - 1);
                }
                shifts = new long[shiftingVariables.length];

                changedMainVariable = previousBaseVariable != baseVariable;
            }

            aReduced = AMultivariatePolynomial.swapVariables(a, 0, baseVariable);
            bReduced = AMultivariatePolynomial.swapVariables(b, 0, baseVariable);

            if (nAttemptsWithZeros > ATTEMPTS_WITH_MAX_ZEROS && zeroVariables.length != 0) {
                // fill with some value some zero variable
                shiftingVariables = Arrays.copyOf(shiftingVariables, shiftingVariables.length + 1);
                shiftingVariables[shiftingVariables.length - 1] = zeroVariables[zeroVariables.length - 1];
                zeroVariables = Arrays.copyOf(zeroVariables, zeroVariables.length - 1);
                shifts = new long[shiftingVariables.length];
                nAttemptsWithZeros = 0;
            }

            ++nAttemptsWithSameBaseVariable;
            ++nAttemptsWithZeros;
            int
                    ubDegree = bReduced.degree(0),
                    uaDegree = aReduced.degree(0);

            // choosing random shifts
            for (int nFails = 0; ; ++nFails) {
                if (nFails > 8) {
                    // base variable is unlucky -> start over
                    nAttemptsWithSameBaseVariable = SAME_BASE_VARIABLE_ATTEMPTS + 1;
                    return nextEvaluation();
                }
                for (int i = 0; i < shiftingVariables.length; ++i)
                    shifts[i] = a.domain.randomElement(rnd);

                lUnivariatePolynomialZp ub = bReduced
                        .evaluateAtZero(zeroVariables)
                        .evaluate(shiftingVariables, shifts)
                        .asUnivariate();

                if (ub.isZeroAt(0) || ub.isZeroAt(ubDegree))
                    continue;

                lUnivariatePolynomialZp ua = aReduced
                        .evaluateAtZero(zeroVariables)
                        .evaluate(shiftingVariables, shifts)
                        .asUnivariate();

                if (ua.isZeroAt(0) || ua.isZeroAt(uaDegree))
                    continue;

                break;
            }
            aReduced = aReduced.shift(shiftingVariables, shifts);
            bReduced = bReduced.shift(shiftingVariables, shifts);

            return changedMainVariable;
        }

        lMultivariatePolynomialZp convert(lMultivariatePolynomialZp poly) {
            if (currentPerfectVariable != -1)
                // just swap the vars
                return currentPerfectVariable == 0 ?
                        poly.clone() : lMultivariatePolynomialZp.swapVariables(poly, 0, currentPerfectVariable);
            return AMultivariatePolynomial.swapVariables(poly, 0, baseVariable).shift(shiftingVariables, shifts);
        }

        lMultivariatePolynomialZp reconstruct(lMultivariatePolynomialZp poly) {
            if (currentPerfectVariable != -1)
                // just swap the vars
                return currentPerfectVariable == 0 ?
                        poly : lMultivariatePolynomialZp.swapVariables(poly, 0, currentPerfectVariable);

            poly = poly.shift(shiftingVariables, ArraysUtil.negate(shifts.clone()));
            poly = lMultivariatePolynomialZp.swapVariables(poly, 0, baseVariable);
            return poly;
        }
    }


    /* =============================================== EEZ-GCD algorithm ============================================ */

    /**
     * Calculates GCD of two multivariate polynomials over Zp using enhanced EZ algorithm
     *
     * @param a the first multivariate polynomial
     * @param b the second multivariate polynomial
     * @return greatest common divisor of {@code a} and {@code b}
     */
    @SuppressWarnings("unchecked")
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly EEZGCD(Poly a, Poly b) {

        // prepare input and test for early termination
        GCDInput<Term, Poly> gcdInput = preparedGCDInput(a, b, MultivariateGCD::EEZGCD);
        if (gcdInput.earlyGCD != null)
            return gcdInput.earlyGCD;

        a = gcdInput.aReduced;
        b = gcdInput.bReduced;

        // remove content in each variable
        Poly content = a.createOne();
        for (int i = 0; i < a.nVariables; ++i) {
            if (a.degree(i) == 0 || b.degree(i) == 0)
                continue;
            Poly tmpContent = contentGCD(a, b, i, MultivariateGCD::EEZGCD);
            a = divideExact(a, tmpContent);
            b = divideExact(b, tmpContent);
            content = content.multiply(tmpContent);
        }

        // one more reduction; removing of content may shuffle required variables order
        GCDInput<Term, Poly> gcdInput2 = preparedGCDInput(a, b, MultivariateGCD::EEZGCD);
        a = gcdInput2.aReduced; b = gcdInput2.bReduced;

        Poly result = gcdInput2.earlyGCD != null
                ? gcdInput2.earlyGCD
                : gcdInput2.restoreGCD(EEZGCD0(a, b, PrivateRandom.getRandom()));

        result = result.multiply(content);
        return gcdInput.restoreGCD(result);
    }

    @SuppressWarnings("unchecked")
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    IEvaluation<Term, Poly> createEvaluation(Poly factory, RandomGenerator rnd) {
        if (factory instanceof lMultivariatePolynomialZp) {
            lMultivariatePolynomialZp lFactory = (lMultivariatePolynomialZp) factory;
            // set new evaluation point
            long[] substitutions = new long[factory.nVariables - 1];
            for (int j = 0; j < substitutions.length; j++)
                substitutions[j] = lFactory.domain.randomNonZeroElement(rnd);

            return (IEvaluation<Term, Poly>) new HenselLifting.lEvaluation(lFactory.nVariables, substitutions, lFactory.domain, lFactory.ordering);
        } else if (factory instanceof MultivariatePolynomial) {
            MultivariatePolynomial lFactory = (MultivariatePolynomial) factory;
            // set new evaluation point
            Object[] substitutions = lFactory.domain.createArray(lFactory.nVariables - 1);
            for (int j = 0; j < substitutions.length; j++)
                substitutions[j] = lFactory.domain.randomNonZeroElement(rnd);

            return (IEvaluation<Term, Poly>) new HenselLifting.Evaluation(lFactory.nVariables, substitutions, lFactory.domain, lFactory.ordering);
        } else
            throw new RuntimeException();
    }

    @SuppressWarnings("unchecked")
    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly EEZGCD0(Poly a, Poly b, RandomGenerator rnd) {

        // degree of univariate gcd
        int ugcdDegree = Integer.MAX_VALUE;
        // degrees of a and b as Z[y1, ... ,yN][x]
        int
                uaDegree = a.degree(0),
                ubDegree = b.degree(0);

        Set<DegreeVector> aSkeleton = a.getSkeleton(0), bSkeleton = b.getSkeleton(0);
        choose_evaluation:
        while (true) {
            IEvaluation<Term, Poly> evaluation = createEvaluation(a, rnd);
            Poly
                    mua = evaluation.evaluateFrom(a, 1),
                    mub = evaluation.evaluateFrom(b, 1);

            if (!aSkeleton.equals(mua.getSkeleton()) || !bSkeleton.equals(mub.getSkeleton()))
                continue;

            IUnivariatePolynomial
                    ua = mua.asUnivariate(),
                    ub = mub.asUnivariate();

            assert ua.degree() == uaDegree;
            assert ub.degree() == ubDegree;

            // gcd of a mod I and b mod I (univariate)
            IUnivariatePolynomial ugcd = UnivariateGCD.PolynomialGCD(ua, ub);

            if (ugcd.degree() == 0)
                // coprime polynomials
                return a.createOne();

            if (ugcd.degree() > ugcdDegree)
                // unlucky evaluation
                continue choose_evaluation;

            ugcdDegree = ugcd.degree();

            if (ugcdDegree == uaDegree) {
                // a is a divisor of b
                if (dividesQ(b, a))
                    return a;
                //continue choose_evaluation;
            }

            if (ugcdDegree == ubDegree) {
                // b is a divisor of a
                if (dividesQ(a, b))
                    return b;
                //continue choose_evaluation;
            }

            // base polynomial to lift (either a or b)
            Poly base = null;
            // a cofactor with gcd to lift, i.e. base = ugcd * uCoFactor mod I
            IUnivariatePolynomial uCoFactor = null;
            Poly coFactorLC = null;

            // deg(b) < deg(a), so it is better to try to lift b
            IUnivariatePolynomial ubCoFactor = DivisionWithRemainder.quotient(ub, ugcd, true);
            if (UnivariateGCD.PolynomialGCD(ugcd, ubCoFactor).isConstant()) {
                // b splits into coprime factors
                base = b;
                uCoFactor = ubCoFactor;
                coFactorLC = b.lc(0);
            }

            if (base == null) {
                // b does not split into coprime factors => try a
                IUnivariatePolynomial uaCoFactor = DivisionWithRemainder.quotient(ua, ugcd, true);
                if (UnivariateGCD.PolynomialGCD(ugcd, uaCoFactor).isConstant()) {
                    base = a;
                    uCoFactor = uaCoFactor;
                    coFactorLC = a.lc(0);
                }
            }

            if (base == null) {
                // neither a nor b does not split into coprime factors => square free decomposition required
                Poly
                        bRepeatedFactor = EEZGCD(b, b.derivative(0, 1)),
                        squareFreeGCD = EEZGCD(divideExact(b, bRepeatedFactor), a),
                        aRepeatedFactor = divideExact(a, squareFreeGCD),
                        gcd = squareFreeGCD.clone();

                IUnivariatePolynomial
                        uSquareFreeGCD = evaluation.evaluateFrom(squareFreeGCD, 1).asUnivariate(),
                        uaRepeatedFactor = evaluation.evaluateFrom(aRepeatedFactor, 1).asUnivariate(),
                        ubRepeatedFactor = evaluation.evaluateFrom(bRepeatedFactor, 1).asUnivariate();
                while (true) {
                    ugcd = UnivariateGCD.PolynomialGCD(
                            uSquareFreeGCD, uaRepeatedFactor, ubRepeatedFactor);

                    if (ugcd.degree() == 0)
                        return gcd;

                    if (!uSquareFreeGCD.clone().monic().equals(ugcd.clone().monic())) {
                        Poly
                                mgcd = (Poly) AMultivariatePolynomial.asMultivariate(ugcd, a.nVariables, 0, a.ordering),
                                gcdCoFactor = (Poly) AMultivariatePolynomial.asMultivariate(
                                        DivisionWithRemainder.divideExact(uSquareFreeGCD, ugcd, false), a.nVariables, 0, a.ordering);

                        HenselLifting.liftWang(squareFreeGCD, mgcd, gcdCoFactor, evaluation);

                        squareFreeGCD = mgcd.clone();
                        uSquareFreeGCD = evaluation.evaluateFrom(squareFreeGCD, 1).asUnivariate();
                    }

                    gcd = gcd.multiply(squareFreeGCD);
                    uaRepeatedFactor = DivisionWithRemainder.divideExact(uaRepeatedFactor, ugcd, false);
                    ubRepeatedFactor = DivisionWithRemainder.divideExact(ubRepeatedFactor, ugcd, false);
                }
            }

            Poly gcd = AMultivariatePolynomial.asMultivariate(ugcd, a.nVariables, 0, a.ordering);
            Poly coFactor = AMultivariatePolynomial.asMultivariate(uCoFactor, a.nVariables, 0, a.ordering);

            // impose the leading coefficient
            Poly lcCorrection = EEZGCD(a.lc(0), b.lc(0));
            //assert ZippelGCD(a.lc(0), b.lc(0)).monic(lcCorrection.lc()).equals(lcCorrection) : "\n" + a.lc(0) + "  \n " + b.lc(0);

            if (lcCorrection.isOne()) {
                HenselLifting.liftWang(base, gcd, coFactor, null, base.lc(0), evaluation);
            } else {
                Poly lcCorrectionMod = evaluation.evaluateFrom(lcCorrection, 1);
                assert lcCorrectionMod.isConstant();

                coFactor = coFactor.multiplyByLC(gcd);
                gcd = gcd.monicWithLC(lcCorrectionMod);

                HenselLifting.liftWang(base.clone().multiply(lcCorrection), gcd, coFactor, lcCorrection, coFactorLC, evaluation);
                assert gcd.lc(0).equals(lcCorrection);
                gcd = HenselLifting.primitivePart(gcd);
            }

            if (dividesQ(b, gcd) && dividesQ(a, gcd))
                return gcd;
        }
    }
}
