package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MultivariateFactorization;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

import static cc.redberry.rings.poly.PolynomialMethods.polyPow;
import static cc.redberry.rings.poly.multivar.MultivariatePolynomial.parse;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class Rings_vs_Singular_vs_Mathematica_Factor3 {
    public static void main(String[] args) throws Exception {
        for (Ring<BigInteger> ring : new Ring[]{Rings.Z, Rings.Zp(2), Rings.Zp(524287)}) {
            MultivariatePolynomial<BigInteger>
                    p1 = polyPow(parse("1 + 3*a*b + 5*b*c + 7*c*d + 9*d*e + 11*e*f + 13*f*g + 15*g*a", ring), 3),
                    p2 = polyPow(parse("1 + 3*a*c + 5*b*d + 7*c*e + 9*f*e + 11*g*f + 13*f*a + 15*g*b", ring), 3),
                    p3 = polyPow(parse("1 + 3*a*d + 5*b*e + 7*c*f + 9*f*g + 11*g*a + 13*f*b + 15*g*c", ring), 3),
                    poly = p1.multiply(p2, p3);
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
