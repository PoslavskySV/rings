package cc.redberry.rings.poly.univar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.Rings;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.test.APolynomialTest;
import cc.redberry.rings.util.RandomDataGenerator;
import org.apache.commons.math3.random.RandomGenerator;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class DiophantineEquationsTest extends APolynomialTest {
    @Test
    public void test1() throws Exception {
        RandomGenerator rnd = getRandom();
        RandomDataGenerator rndd = getRandomData();
        for (long modulus : Arrays.asList(2, 3, 17, Integer.MAX_VALUE)) {
            IntegersZp64 ring = Rings.Zp64(modulus);
            UnivariateRing<UnivariatePolynomialZp64> uring = Rings.UnivariateRingZp64(ring);
            for (int i = 0; i < 1000; i++) {
                int nPolynomials = rndd.nextInt(2, 10);
                UnivariatePolynomialZp64[] polys = new UnivariatePolynomialZp64[nPolynomials];
                for (int j = 0; j < polys.length; j++)
                    polys[j] = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(15, 30), ring.modulus, rnd);

                UnivariatePolynomialZp64 gcd = UnivariateGCD.PolynomialGCD(polys);
                UnivariatePolynomialZp64 rhs = RandomUnivariatePolynomials.randomMonicPoly(rndd.nextInt(0, 10), ring.modulus, rnd);
                rhs.multiply(gcd);

                DiophantineEquations.DiophantineSolver<UnivariatePolynomialZp64> solver = new DiophantineEquations.DiophantineSolver<>(polys);
                UnivariatePolynomialZp64 g = uring.getZero();
                for (int l = 0; l < solver.solution.length; l++)
                    g.add(solver.solution[l].clone().multiply(polys[l]));

                UnivariatePolynomialZp64[] solve = solver.solve(rhs);
                for (int j = 0; j < solve.length; j++)
                    rhs.subtract(solve[j].multiply(polys[j]));
                Assert.assertTrue(rhs.isZero());
            }
        }
    }
}