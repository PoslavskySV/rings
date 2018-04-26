package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Comparator;
import java.util.function.Function;

/**
 * Methods to generate random multivariate polynomials.
 *
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
            int[] exponents = RandomUtil.randomIntArray(nVars, 0, nd, rnd);
            terms[i] = new Monomial<>(exponents, cfx);
        }
        return MultivariatePolynomial.create(nVars, Rings.Z, ordering, terms);
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
        return randomPolynomial(nVars, degree, size, BigInteger.TEN, MonomialOrder.DEFAULT, rnd);
    }

    /**
     * Generates random polynomial
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param ring     the coefficient ring
     * @param ordering monomial order
     * @param method   method for generating random coefficients
     * @param rnd      random source
     * @return random polynomial
     */
    public static <E> MultivariatePolynomial<E> randomPolynomial(int nVars, int degree, int size, Ring<E> ring,
                                                                 Comparator<DegreeVector> ordering,
                                                                 Function<RandomGenerator, E> method, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        @SuppressWarnings("unchecked")
        Monomial<E>[] terms = new Monomial[size];
        for (int i = 0; i < size; i++) {
            E cfx = method.apply(rnd);
            int[] exponents = RandomUtil.randomIntArray(nVars, 0, nd, rnd);
            terms[i] = new Monomial<>(exponents, cfx);
        }
        return MultivariatePolynomial.create(nVars, ring, ordering, terms);
    }

    /**
     * Generates random polynomial
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param ring     the coefficient ring
     * @param ordering monomial order
     * @param rnd      random source
     * @return random polynomial
     */
    public static <E> MultivariatePolynomial<E> randomPolynomial(int nVars, int degree, int size, Ring<E> ring, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, ring, ordering, ring::randomElement, rnd);
    }

    /**
     * Generates random Zp[X] polynomial over machine integers
     *
     * @param nVars  number of variables
     * @param degree maximal degree of the result
     * @param size   number of elements in the result
     * @param ring   the coefficient ring
     * @param rnd    random source
     * @return random polynomial
     */
    public static MultivariatePolynomialZp64 randomPolynomial(int nVars, int degree, int size, IntegersZp64 ring, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, ring, MonomialOrder.DEFAULT, rnd);
    }

    /**
     * Generates random Zp[X] polynomial over machine integers
     *
     * @param nVars    number of variables
     * @param degree   maximal degree of the result
     * @param size     number of elements in the result
     * @param ring     the coefficient ring
     * @param ordering monomial order
     * @param rnd      random source
     * @return random polynomial
     */
    public static MultivariatePolynomialZp64 randomPolynomial(int nVars, int degree, int size, IntegersZp64 ring, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        @SuppressWarnings("unchecked")
        MonomialZp64[] terms = new MonomialZp64[size];
        for (int i = 0; i < size; i++) {
            long cfx = ring.randomElement(rnd);
            int[] exponents = RandomUtil.randomIntArray(nVars, 0, nd, rnd);
            terms[i] = new MonomialZp64(exponents, cfx);
        }
        return MultivariatePolynomialZp64.create(nVars, ring, ordering, terms);
    }

    /**
     * Generates random multivariate polynomial
     *
     * @param factory factory polynomial
     * @param degree  maximal degree of the result
     * @param size    number of elements in the result
     * @param rnd     random source
     * @return random polynomial
     */
    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly randomPolynomial(Poly factory, int degree, int size, RandomGenerator rnd) {
        if (factory instanceof MultivariatePolynomialZp64)
            return (Poly) randomPolynomial(((MultivariatePolynomialZp64) factory).nVariables, degree, size, ((MultivariatePolynomialZp64) factory).ring, ((MultivariatePolynomialZp64) factory).ordering, rnd);
        else
            return (Poly) randomPolynomial(((MultivariatePolynomial) factory).nVariables, degree, size, ((MultivariatePolynomial) factory).ring, ((MultivariatePolynomial) factory).ordering, rnd);
    }
}
