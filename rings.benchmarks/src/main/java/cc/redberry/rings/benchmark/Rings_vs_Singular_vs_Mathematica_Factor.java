package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.primes.SmallPrimes;
import cc.redberry.rings.util.TimeUnits;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static cc.redberry.rings.benchmark.Bench.*;

public class Rings_vs_Singular_vs_Mathematica_Factor {
    public static void main(String[] args) throws Exception {
        // warm up
        run(3, 5, 3, 2, 100, Rings.Z, false);
        System.out.println("warmed");
        long[][] timings;
        silent = false;


        for (int nVariables = 3; nVariables <= 5; nVariables++) {
            MATHEMATICA_TIMEOUT_SECONDS = 100;
            doFactorMathematica = nVariables <= 4;

            System.out.println("#variables = " + nVariables);

            int nIterations = 900 / nVariables / nVariables;
            int size = 20;
            int degree = 10;

            timings = run(nVariables, degree, size, 2, nIterations, Rings.Z, false);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z_%s.tsv", nVariables)), timings);
            timings = run(nVariables, degree, size, 2, nIterations, Rings.Z, true);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z_coprime_%s.tsv", nVariables)), timings);

            for (long prime : new long[]{(1 << 19) - 1, 2}) {
                System.out.println("Modulus: " + prime);
                timings = run(nVariables, degree, size, 2, nIterations, Rings.Zp(prime), false);
                writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z%s_%s.tsv", prime, nVariables)), timings);
                timings = run(nVariables, degree, size, 2, nIterations, Rings.Zp(prime), true);
                writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z%s_coprime_%s.tsv", prime, nVariables)), timings);
            }
        }
    }

    static void writeTimingsTSV(Path path, long[][] timings) throws IOException {
        Files.write(path, (Iterable<String>)
                () -> Arrays.stream(timings)
                        .map(row -> Arrays.stream(row)
                                .mapToObj(String::valueOf)
                                .collect(Collectors.joining("\t")))
                        .iterator(), StandardOpenOption.CREATE);
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
    @SuppressWarnings("unchecked")
    static long[][] run(int nVariables,
                        int degree,
                        int size,
                        int nFactors,
                        int nIterations,
                        Ring<BigInteger> cfRing,
                        boolean coprime) throws Exception {
        System.out.println("\tRings\tSingular\tMathematica");
        long[][] timings = new long[nIterations][];
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(nVariables, cfRing);
        for (int i = 0; i < nIterations; i++) {
            MultivariatePolynomial<BigInteger>[] factors
                    = IntStream.range(0, nFactors)
                    .mapToObj(__ -> ring.randomElement(degree, size))
                    .filter(__ -> !__.isZero())
                    .map(__ -> __.setRing(Rings.Zp(1000)).setRing(Rings.Z))// reduce coefficients
                    .toArray(MultivariatePolynomial[]::new);

            MultivariatePolynomial<BigInteger> poly = ring.multiply(factors);

            if (coprime)
                poly.increment();

            long start = System.nanoTime();
            FactorDecomposition<MultivariatePolynomial<BigInteger>> ringsResult = PolynomialMethods.Factor(poly);
            long ringsTime = System.nanoTime() - start;

            ExternalResult singularResult = singularFactor(poly);
            long singularTime = singularResult.nanoTime;

            ExternalResult mmaResult = mathematicaFactor(poly);
            long mmaTime = mmaResult.nanoTime;

            timings[i] = new long[]{ringsTime, singularTime, mmaTime};
            if (!silent)
                System.out.println(i + "\t" + TimeUnits.nanosecondsToString(ringsTime) + "\t" + TimeUnits.nanosecondsToString(singularTime) + "\t" + TimeUnits.nanosecondsToString(mmaTime));
        }
        return timings;
    }
}
