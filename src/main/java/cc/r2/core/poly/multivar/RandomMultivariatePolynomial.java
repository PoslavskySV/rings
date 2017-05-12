package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Integers;
import cc.r2.core.poly.multivar.MultivariatePolynomial.DegreeVector;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Comparator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class RandomMultivariatePolynomial {
    public static MultivariatePolynomial<BigInteger> randomPolynomial(int nVars, int degree, int size, BigInteger bound, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, bound, Integers.Integers, ordering, rnd);
    }

    public static MultivariatePolynomial<BigInteger> randomPolynomial(int nVars, int degree, int size, BigInteger bound, Domain<BigInteger> domain, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        BigInteger[] cfx = new BigInteger[size];
        DegreeVector[] dvs = new DegreeVector[size];
        for (int i = 0; i < size; i++) {
            dvs[i] = new DegreeVector(RandomUtil.randomIntArray(nVars, 0, nd, rnd));
            cfx[i] = domain.valueOf(RandomUtil.randomInt(bound, rnd));
            if (rnd.nextBoolean() && rnd.nextBoolean())
                cfx[i] = domain.negate(cfx[i]);
        }
        return MultivariatePolynomial.create(domain, ordering, dvs, cfx);
    }

    public static MultivariatePolynomial<BigInteger> randomPolynomial(int nVars, int degree, int size, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, BigInteger.TEN, MultivariatePolynomial.LEX, rnd);
    }
}
