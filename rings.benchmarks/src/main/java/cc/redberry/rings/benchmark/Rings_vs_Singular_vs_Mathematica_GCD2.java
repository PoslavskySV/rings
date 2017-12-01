package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateGCD;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 2.2
 */
public class Rings_vs_Singular_vs_Mathematica_GCD2 {
    public static void main(String[] args) throws Exception {

        for (Ring<BigInteger> ring : new Ring[]{Rings.Z, Rings.Zp((1 << 19) - 1)}) {
            System.out.println("Ring: " + ring);

            for (int exp = 4; exp < 7; exp++) {
                System.out.println("\n\n\n\n");
                System.out.println("<><><><><><><><><><><><><><><><><><><><><><><><><><><>");
                System.out.println("exp: " + exp);

                MultivariatePolynomial<BigInteger>
                        a = MultivariatePolynomial.parse("1 + 3*a + 5*b + 7*c + 9*d + 11*e + 13*f + 15*g", ring),
                        b = MultivariatePolynomial.parse("1 - 3*a + 5*b - 7*c + 9*d - 11*e + 13*f - 15*g", ring),
                        g = MultivariatePolynomial.parse("1 + 3*a - 5*b + 7*c - 9*d + 11*e - 13*f + 15*g", ring);

                a = PolynomialMethods.polyPow(a, exp);
                a.decrement();
                b = PolynomialMethods.polyPow(b, exp);
                g = PolynomialMethods.polyPow(g, exp);
                g.increment();

                MultivariatePolynomial<BigInteger> ag = a.clone().multiply(g);
                MultivariatePolynomial<BigInteger> bg = b.clone().multiply(g);

//                info(a, "a");
//                info(b, "b");
//                info(g, "g");
//
//                info(ag, "ag");
//                info(bg, "bg");

                System.out.println("\n=================\n");
                for (int i = 0; i < 5; ++i) {
                    long start = System.nanoTime();
                    int size = MultivariateGCD.PolynomialGCD(ag, bg).size();
                    long ringsTime = System.nanoTime() - start;

                    System.out.println("Rings: " + TimeUnits.nanosecondsToString(ringsTime));

                    long singularTime = Bench.singularGCD(ag, bg).nanoTime;
                    System.out.println("Singular: " + TimeUnits.nanosecondsToString(singularTime));

                    long mmaTime = Bench.mathematicaGCD(ag, bg).nanoTime;
                    System.out.println("MMA: " + TimeUnits.nanosecondsToString(mmaTime));

                    System.out.println();
                }
            }
        }
    }

    private static void info(AMultivariatePolynomial p) {
        info(p, "");
    }

    private static void info(AMultivariatePolynomial p, String poly) {
        System.out.println();
        if (!poly.isEmpty())
            System.out.println("Info for " + poly);
        System.out.println("Ring:      " + p.coefficientRingToString());
        System.out.println("Degrees:   " + Arrays.toString(p.degrees()));
        System.out.println("Size:      " + p.size());
        System.out.println("Sparsity:  " + p.sparsity());
        System.out.println("Sparsity2: " + p.sparsity2());
    }
}
