package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariateFactorization;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class Rings_vs_Singular_vs_Mathematica_Factor2 {

    public static void main(String[] args) throws Exception {
        for (Ring<BigInteger> ring : new Ring[]{Rings.Zp(2), Rings.Z, Rings.Zp(524287)}) {
            MultivariatePolynomial<BigInteger> poly = MultivariatePolynomial.parse("1 + 3*a + 5*b + 7*c + 9*d + 11*e + 13*f + 15*g", ring);
            poly = PolynomialMethods.polyPow(poly, 15);
            poly.decrement();

            System.out.println("Ring: " + ring);

            for (int i = 0; i < 5; ++i) {
                long start = System.nanoTime();
                int size = MultivariateFactorization.Factor(poly).size();
                long ringsTime = System.nanoTime() - start;
                if (i == 0)
                    System.out.println("#factors = " + size);

                System.out.println(TimeUnits.nanosecondsToString(ringsTime));

                // MMA and Singular should be called directly
                //
                // long singularTime = Bench.singularFactor(poly).nanoTime;
                // System.out.println(singularTime);
                //
                // long mmaTime = Bench.singularFactor(poly).nanoTime;
                // System.out.println(mmaTime);
                //
                // System.out.println(i + "\t" + TimeUnits.nanosecondsToString(ringsTime) + "\t" + TimeUnits.nanosecondsToString(singularTime) + "\t" + TimeUnits.nanosecondsToString(mmaTime));
            }
        }
    }
}
