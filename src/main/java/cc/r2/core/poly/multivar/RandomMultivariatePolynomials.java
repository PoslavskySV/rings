package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Domains;
import cc.r2.core.poly.IntegersZp64;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Comparator;

/**
 * Methods to generate random multivariate polynomials.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class RandomMultivariatePolynomials {
    private RandomMultivariatePolynomials() { }

    /**
     * Generates random Z[X] polynomial with coefficients bounded by {@code bound}
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param bound    coefficient bound
     * @param ordering monomial order
     * @param rnd      random source
     * @return random Z[X] polynomial
     */
    public static MultivariatePolynomial<BigInteger> randomPolynomial(int nVars, int degree, int size, BigInteger bound, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        @SuppressWarnings("unchecked")
        Monomial<BigInteger>[] terms = new Monomial[size];
        for (int i = 0; i < size; i++) {
            BigInteger cfx = RandomUtil.randomInt(bound, rnd);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                cfx = cfx.negate();
            terms[i] = new Monomial<>(RandomUtil.randomIntArray(nVars, 0, nd, rnd), cfx);
        }
        return MultivariatePolynomial.create(nVars, Domains.Z, ordering, terms);
    }

    /**
     * Generates random Z[X] polynomial
     *
     * @param nVars  number of variables
     * @param degree maximal degree of the result
     * @param size   number of elements in the result
     * @param rnd    random source
     * @return random polynomial
     */
    public static MultivariatePolynomial<BigInteger> randomPolynomial(int nVars, int degree, int size, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, BigInteger.TEN, MonomialOrder.LEX, rnd);
    }

    /**
     * Generates random polynomial
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param domain   the coefficient domain
     * @param ordering monomial order
     * @param rnd      random source
     * @return random polynomial
     */
    public static <E> MultivariatePolynomial<E> randomPolynomial(int nVars, int degree, int size, Domain<E> domain, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        @SuppressWarnings("unchecked")
        Monomial<E>[] terms = new Monomial[size];
        for (int i = 0; i < size; i++) {
            E cfx = domain.randomElement(rnd);
            terms[i] = new Monomial<>(RandomUtil.randomIntArray(nVars, 0, nd, rnd), cfx);
        }
        return MultivariatePolynomial.create(nVars, domain, ordering, terms);
    }

    /**
     * Generates random Zp[X] polynomial over machine integers
     *
     * @param nVars  number of variables
     * @param degree maximal degree of the result
     * @param size   number of elements in the result
     * @param domain the coefficient domain
     * @param rnd    random source
     * @return random polynomial
     */
    public static MultivariatePolynomialZp64 randomPolynomial(int nVars, int degree, int size, IntegersZp64 domain, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, domain, MonomialOrder.LEX, rnd);
    }

    /**
     * Generates random Zp[X] polynomial over machine integers
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param domain   the coefficient domain
     * @param ordering monomial order
     * @param rnd      random source
     * @return random polynomial
     */
    public static MultivariatePolynomialZp64 randomPolynomial(int nVars, int degree, int size, IntegersZp64 domain, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        @SuppressWarnings("unchecked")
        MonomialZp64[] terms = new MonomialZp64[size];
        for (int i = 0; i < size; i++) {
            long cfx = domain.randomElement(rnd);
            terms[i] = new MonomialZp64(RandomUtil.randomIntArray(nVars, 0, nd, rnd), cfx);
        }
        return MultivariatePolynomialZp64.create(nVars, domain, ordering, terms);
    }
}
