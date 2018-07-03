package cc.redberry.rings.poly.univar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.poly.Util;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

import static cc.redberry.rings.poly.univar.Conversions64bit.asOverZp64;
import static cc.redberry.rings.poly.univar.Conversions64bit.canConvertToZp64;

/**
 * Irreducibility tests and generators for random irreducible polynomials.
 *
 * @since 1.0
 */
public final class IrreduciblePolynomials {
    private IrreduciblePolynomials() {}

    /**
     * Tests whether {@code poly} is irreducible
     *
     * @param poly the polynomial
     * @return whether {@code poly} is an irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> boolean irreducibleQ(Poly poly) {
        if (poly.isOverFiniteField())
            return finiteFieldIrreducibleQ(poly);
        else
            return UnivariateFactorization.Factor(poly).isTrivial();
    }

    /**
     * Tests whether {@code poly} is irreducible over the finite field
     *
     * @param poly the polynomial over finite field
     * @return whether {@code poly} is an irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> boolean finiteFieldIrreducibleQ(Poly poly) {
        return finiteFieldIrreducibleBenOr(poly);
    }

    /**
     * Tests whether {@code poly} is irreducible over the finite field
     *
     * @param poly the polynomial over finite field
     * @return whether {@code poly} is an irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> boolean
    finiteFieldIrreducibleViaModularComposition(Poly poly) {
        Util.ensureOverFiniteField(poly);

        if (canConvertToZp64(poly))
            return finiteFieldIrreducibleViaModularComposition(asOverZp64(poly));

        if (poly.degree() <= 1)
            return true;

        poly = poly.clone().monic();
        UnivariateDivision.InverseModMonomial<Poly> invMod = UnivariateDivision.fastDivisionPreConditioning(poly);
        // x^q
        Poly xq = UnivariatePolynomialArithmetic.createMonomialMod(poly.coefficientRingCardinality(), poly, invMod);

        // cached powers x^(q^i) for different i
        TIntObjectMap<Poly> cache = new TIntObjectHashMap<>();

        int degree = poly.degree();

        // x^(q^n)
        Poly xqn = composition(xq, degree, poly, invMod, cache);
        assert
                xqn.equals(UnivariatePolynomialArithmetic.createMonomialMod(BigIntegerUtil.pow(poly.coefficientRingCardinality(), degree), poly, invMod))
                : "\n" + xqn + "\n" + UnivariatePolynomialArithmetic.createMonomialMod(BigIntegerUtil.pow(poly.coefficientRingCardinality(), degree), poly, invMod);

        Poly xMonomial = poly.createMonomial(1);
        if (!xqn.equals(xMonomial))
            // x^(q^n) != x
            return false;

        // primes factors of poly.degree
        int[] primes = ArraysUtil.getSortedDistinct(SmallPrimes.primeFactors(degree));
        for (int prime : primes) {
            // x^(q^(n/p))
            Poly xqp = composition(xq, degree / prime, poly, invMod, cache);
            if (!UnivariateGCD.PolynomialGCD(xqp.subtract(xMonomial), poly).isOne())
                return false;
        }
        return true;
    }

    /* fives xq^exponent using repeated compositions */
    static <Poly extends IUnivariatePolynomial<Poly>> Poly composition(
            Poly xq, int exponent,
            Poly poly, UnivariateDivision.InverseModMonomial<Poly> invMod,
            TIntObjectMap<Poly> cache) {
        assert exponent > 0;

        Poly cached = cache.get(exponent);
        if (cached != null)
            return cached;

        Poly result = xq.createMonomial(1);
        Poly k2p = xq;
        int rExp = 0, kExp = 1;
        for (; ; ) {
            if ((exponent & 1) != 0)
                cache.put(rExp += kExp, result = ModularComposition.composition(k2p, result, poly, invMod));
            exponent = exponent >> 1;
            if (exponent == 0) {
                cache.put(rExp, result);
                return result;
            }
            cache.put(kExp *= 2, k2p = ModularComposition.composition(k2p, k2p, poly, invMod));
        }
    }

    /**
     * Tests whether {@code poly} is irreducible over the finite field
     *
     * @param poly the polynomial over finite field
     * @return whether {@code poly} is an irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> boolean
    finiteFieldIrreducibleBenOr(Poly poly) {
        Util.ensureOverFiniteField(poly);

        if (canConvertToZp64(poly))
            return finiteFieldIrreducibleBenOr(asOverZp64(poly));

        if (!poly.isMonic())
            poly = poly.clone().monic();

        UnivariateDivision.InverseModMonomial<Poly> invMod = UnivariateDivision.fastDivisionPreConditioning(poly);
        // x
        Poly xMonomial = poly.createMonomial(1);
        // x^q
        Poly xq = xMonomial.clone();
        // primes factors of poly.degree
        for (int i = 1; i <= (poly.degree() / 2); ++i) {
            // x^(q^i)
            xq = UnivariatePolynomialArithmetic.polyPowMod(xq, poly.coefficientRingCardinality(), poly, invMod, false);
            if (!UnivariateGCD.PolynomialGCD(xq.clone().subtract(xMonomial), poly).isOne())
                return false;
        }
        return true;
    }

    /**
     * Generated random irreducible Zp polynomial of degree {@code degree}
     *
     * @param modulus the modulus
     * @param degree  the degree
     * @param rnd     random source
     * @return irreducible polynomial
     */
    public static UnivariatePolynomialZp64 randomIrreduciblePolynomial(long modulus, int degree, RandomGenerator rnd) {
        return randomIrreduciblePolynomial(UnivariatePolynomialZp64.zero(modulus), degree, rnd);
    }


    /**
     * Generated random irreducible polynomial over specified ring of degree {@code degree}
     *
     * @param ring   coefficient ring
     * @param degree the degree
     * @param rnd    random source
     * @return irreducible polynomial
     */
    public static <E> UnivariatePolynomial<E> randomIrreduciblePolynomial(Ring<E> ring, int degree, RandomGenerator rnd) {
        UnivariatePolynomial<E> poly;
        do {
            poly = RandomUnivariatePolynomials.randomPoly(degree, ring, rnd);
            if (ring.isField())
                poly.monic();
        } while (!irreducibleQ(poly));
        assert poly.degree == degree;
        return poly;
    }

    /**
     * Generated random irreducible polynomial over Z
     *
     * @param degree the degree
     * @param rnd    random source
     * @return irreducible polynomial over Z
     */
    public static UnivariatePolynomial<BigInteger> randomIrreduciblePolynomialOverZ(int degree, RandomGenerator rnd) {
        // some prime number
        long modulus = SmallPrimes.nextPrime(1 << 25);
        return randomIrreduciblePolynomial(modulus, degree, rnd).toBigPoly().setRing(Rings.Z);
    }

    /**
     * Generated random irreducible polynomial of degree {@code degree}
     *
     * @param factory type marker
     * @param degree  the degree
     * @param rnd     random source
     * @return irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly randomIrreduciblePolynomial(Poly factory, int degree, RandomGenerator rnd) {
        Poly poly;
        do {
            poly = RandomUnivariatePolynomials.randomPoly(factory, degree, rnd);
        } while (!irreducibleQ(poly));
        assert poly.degree() == degree;
        return poly;
    }
}
