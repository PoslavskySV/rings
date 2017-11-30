package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariateGCD;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class Rings_vs_Singular_vs_Mathematica_GCD2 {


//    public static void main(String[] args) throws Exception {
//        // warm up
////        run(3, 5, 3, 100, Rings.Z, false);
//        System.out.println("warmed");
//        silent = false;
//        TDoubleHashSet db = new TDoubleHashSet();
//        for (int nVariables = 5; nVariables <= 7; ++nVariables) {
//
//            int nIterations = 100;
//
//            db.addAll(run(nVariables, 20, 40, nIterations, Rings.Z, false));
//            db.addAll(run(nVariables, 20, 40, nIterations, Rings.Z, true));
//            db.addAll(run(nVariables, 20, 40, nIterations, Rings.Zp(2), false));
//            db.addAll(run(nVariables, 20, 40, nIterations, Rings.Zp(2), true));
//        }
//
//        DescriptiveStatistics ds = new DescriptiveStatistics();
//        for (double v : db.toArray()) {
//            ds.addValue(v);
//        }
//        System.out.println(ds);
//    }
//
//    static double[] run(
//            int nVariables,
//            int degree,
//            int size,
//            int nIterations,
//            Ring<BigInteger> cfRing,
//            boolean coprime) throws Exception {
//        TDoubleHashSet db = new TDoubleHashSet();
//        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(nVariables, cfRing);
//        for (int i = 0; i < nIterations; i++) {
//            MultivariatePolynomial<BigInteger>
//                    a = ring.randomElement(degree, size),
//                    b = ring.randomElement(degree, size),
//                    gcd = ring.randomElement(degree, size);
//            a = a.multiply(gcd);
//            b = b.multiply(gcd);
//            db.add(a.sparsity2());
//            db.add(b.sparsity2());
//        }
//        return db.toArray();
//    }


    static BigInteger formula(int integer, int nPartitions) {
        return BigIntegerUtil.binomial(integer + nPartitions - 1, integer);
    }

    public static void main(String[] args) throws Exception {

        for (Ring<BigInteger> ring : new Ring[]{Rings.Z, Rings.Zp(1031)}) {
            MultivariatePolynomial<BigInteger>
                    a = MultivariatePolynomial.parse("1 + 3*a + 5*b + 7*c + 9*d + 11*e + 13*f + 15*g", ring),
                    b = MultivariatePolynomial.parse("1 - 3*a + 5*b - 7*c + 9*d - 11*e + 13*f - 15*g", ring),
                    g = MultivariatePolynomial.parse("1 + 3*a - 5*b + 7*c - 9*d + 11*e - 13*f + 15*g", ring);
            int exp = 5;
            a = PolynomialMethods.polyPow(a, exp);
            a.decrement();
            b = PolynomialMethods.polyPow(b, exp);
            g = PolynomialMethods.polyPow(g, exp);
            g.increment();

            MultivariatePolynomial<BigInteger> ag = a.clone().multiply(g);
            MultivariatePolynomial<BigInteger> bg = b.clone().multiply(g);


            info(a, "a");
            info(b, "b");
            info(g, "g");

            info(ag, "ag");
            info(bg, "bg");

            for (int i = 0; i < 5; ++i) {
                long start = System.nanoTime();
                int size = MultivariateGCD.PolynomialGCD(ag, bg).size();
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
