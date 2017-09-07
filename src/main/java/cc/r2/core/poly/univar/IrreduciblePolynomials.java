package cc.r2.core.poly.univar;

import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.number.primes.SmallPrimes;
import cc.r2.core.poly.Util;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.univar.DivisionWithRemainder.InverseModMonomial;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class IrreduciblePolynomials {
    private IrreduciblePolynomials() {}

    /**
     * Rests whether {@code poly} is an irreducible polynomial over the finite field
     *
     * @param poly the polynomials over finite field
     * @return whether {@code poly} is an irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> boolean irreducibleQ(Poly poly) {
        Util.ensureFiniteFieldDomain(poly);
        if (poly.degree() <= 1)
            return true;

        poly = poly.clone().monic();
        InverseModMonomial<Poly> invMod = DivisionWithRemainder.fastDivisionPreConditioning(poly);
        // x^q
        Poly xq = PolynomialArithmetics.createMonomialMod(poly.coefficientDomainCardinality(), poly, invMod);

        // cached powers x^(q^i) for different i
        TIntObjectMap<Poly> cache = new TIntObjectHashMap<>();

        int degree = poly.degree();

        // x^(q^n)
        Poly xqn = composition(xq, degree, poly, invMod, cache);
        assert
                xqn.equals(PolynomialArithmetics.createMonomialMod(BigIntegerArithmetics.safePow(poly.coefficientDomainCardinality(), degree), poly, invMod))
                : "\n" + xqn + "\n" + PolynomialArithmetics.createMonomialMod(BigIntegerArithmetics.safePow(poly.coefficientDomainCardinality(), degree), poly, invMod);

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

    /**
     * Generated random irreducible Zp polynomial of degree {@code degree}
     *
     * @param modulus the modulus
     * @param degree  the degree
     * @param rnd     random source
     * @return irreducible polynomial
     */
    public static lUnivariatePolynomialZp randomIrreduciblePolynomial(long modulus, int degree, RandomGenerator rnd) {
        lUnivariatePolynomialZp poly;
        do {
            poly = RandomPolynomials.randomMonicPoly(degree, modulus, rnd);
        } while (!irreducibleQ(poly));
        assert poly.degree == degree;
        return poly;
    }

    /**
     * Generated random irreducible polynomial over finite field of degree {@code degree}
     *
     * @param domain finite field domain
     * @param degree the degree
     * @param rnd    random source
     * @return irreducible polynomial
     */
    public static <E> UnivariatePolynomial<E> randomIrreduciblePolynomial(Domain<E> domain, int degree, RandomGenerator rnd) {
        if (!domain.isFiniteField())
            throw new IllegalArgumentException("Not a finite domain.");
        UnivariatePolynomial<E> poly;
        do {
            poly = RandomPolynomials.randomMonicPoly(degree, domain, rnd);
        } while (!irreducibleQ(poly));
        assert poly.degree == degree;
        return poly;
    }

    /**
     * Generated random irreducible Zp polynomial of degree {@code degree}
     *
     * @param factory tyep marker
     * @param degree  the degree
     * @param rnd     random source
     * @return irreducible polynomial
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly randomIrreduciblePolynomial(Poly factory, int degree, RandomGenerator rnd) {
        if (!factory.isOverFiniteField())
            throw new IllegalArgumentException("Not a finite domain.");
        Poly poly;
        do {
            poly = RandomPolynomials.randomPoly(factory, degree, rnd);
        } while (!irreducibleQ(poly));
        assert poly.degree() == degree;
        return poly;
    }

    /* fives xq^exponent using repeated compositions */
    static <Poly extends IUnivariatePolynomial<Poly>> Poly composition(
            Poly xq, int exponent,
            Poly poly, InverseModMonomial<Poly> invMod,
            TIntObjectMap<Poly> cache) {
        assert exponent > 0;

        Poly cached = cache.get(exponent);
        if (cached != null)
            return cached;

        Poly result = xq.createMonomial(1);
        Poly k2p = xq;
        int rExp = 0, kExp = 1;
        for (; ; ) {
            if ((exponent&1) != 0)
                cache.put(rExp += kExp, result = ModularComposition.composition(k2p, result, poly, invMod));
            exponent = exponent >> 1;
            if (exponent == 0) {
                cache.put(rExp, result);
                return result;
            }
            cache.put(kExp *= 2, k2p = ModularComposition.composition(k2p, k2p, poly, invMod));
        }
    }
}
