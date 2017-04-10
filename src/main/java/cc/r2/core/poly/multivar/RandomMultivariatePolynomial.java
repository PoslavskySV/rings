package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.multivar.MultivariatePolynomial.DegreeVector;
import cc.r2.core.util.RandomUtil;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Comparator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class RandomMultivariatePolynomial {
    public static MultivariatePolynomial randomPolynomial(int nVars, int degree, int size, BigInteger bound, Comparator<DegreeVector> ordering, RandomGenerator rnd) {
        int nd = 3 * degree / 2;
        BigInteger[] cfx = new BigInteger[size];
        DegreeVector[] dvs = new DegreeVector[size];
        for (int i = 0; i < size; i++) {
            dvs[i] = new DegreeVector(RandomUtil.randomIntArray(nVars, 0, nd, rnd));
            cfx[i] = RandomUtil.randomInt(bound, rnd);
            if (rnd.nextBoolean() && rnd.nextBoolean())
                cfx[i] = cfx[i].negate();
        }
        return MultivariatePolynomial.create(cfx, dvs, ordering);
    }

    public static MultivariatePolynomial randomPolynomial(int nVars, int degree, int size, RandomGenerator rnd) {
        return randomPolynomial(nVars, degree, size, BigInteger.TEN, MultivariatePolynomial.LEX, rnd);
    }
}
