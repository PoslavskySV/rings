package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

import java.nio.file.Paths;

import static cc.redberry.rings.benchmark.Bench.*;

public class Rings_vs_Singular_vs_Mathematica_GCD {
    static long SINGULAR_TIMEOUT = 300_000;

    public static void main(String[] args) throws Exception {
        // warm up
        run(3, 5, 3, 100, Rings.Z, false);
        System.out.println("warmed");
        silent = false;
        long[][] timings;
        for (int nVariables = 6; nVariables <= 7; ++nVariables) {
            MATHEMATICA_TIMEOUT_SECONDS = 300;

            int nIterations = 100;

            doFactorMathematica = true;
            doRunSingular = false;
            timings = run(nVariables, 20, 40, nIterations, Rings.Z, false);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings.benchmarks", "results", String.format("gcd_in_z_%s.tsv", nVariables)), timings);
            timings = run(nVariables, 20, 40, nIterations, Rings.Z, true);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings.benchmarks", "results", String.format("gcd_in_z_coprime_%s.tsv", nVariables)), timings);

//            doFactorMathematica = false;
//            timings = run(nVariables, 20, 40, nIterations, Rings.Zp(2), false);
//            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings.benchmarks", "results", String.format("gcd_in_z2_%s.tsv", nVariables)), timings);
//            timings = run(nVariables, 20, 40, nIterations, Rings.Zp(2), true);
//            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings.benchmarks", "results", String.format("gcd_in_z2_coprime_%s.tsv", nVariables)), timings);
        }
    }

    static boolean silent = true;

    /**
     * @param degree      degree of polynomials to test
     * @param size        number of terms in factors
     * @param nIterations number of iterations
     * @param cfRing      coefficient ring
     * @param coprime     whether to make input coprime
     * @return array of timings {rings, mathematica, singular}
     */
    static long[][] run(
            int nVariables,
            int degree,
            int size,
            int nIterations,
            Ring<BigInteger> cfRing,
            boolean coprime) throws Exception {
        long[][] timings = new long[nIterations][];
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(nVariables, cfRing);
        System.out.println("Ring: " + ring);
        for (int i = 0; i < nIterations; i++) {
            if (!silent)
                System.out.println("\n-> " + i);
            MultivariatePolynomial<BigInteger>
                    a = ring.randomElement(degree, size),
                    b = ring.randomElement(degree, size),
                    gcd = ring.randomElement(degree, size);

            a = a.multiply(gcd);
            b = b.multiply(gcd);

            if (coprime)
                b.increment();

            long start = System.nanoTime();
            MultivariatePolynomial<BigInteger> ringsResult = PolynomialMethods.PolynomialGCD(a, b);
            long ringsTime = System.nanoTime() - start;

            if (!silent)
                System.out.println("Rings:" + TimeUnits.nanosecondsToString(ringsTime));

            ExternalResult singularResult = new SingularGCD(a, b).run(SINGULAR_TIMEOUT);
            long singularTime = singularResult.nanoTime;
            if (!silent)
                System.out.println("Singular:" + TimeUnits.nanosecondsToString(singularTime));

            ExternalResult mmaResult = mathematicaGCD(a, b);
            long mmaTime = mmaResult.nanoTime;
            if (!silent)
                System.out.println("MMA:" + TimeUnits.nanosecondsToString(mmaTime));

            timings[i] = new long[]{ringsTime, singularTime, mmaTime};
        }
        return timings;
    }
}
