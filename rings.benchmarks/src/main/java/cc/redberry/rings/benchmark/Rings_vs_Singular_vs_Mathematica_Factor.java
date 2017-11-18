package cc.redberry.rings.benchmark;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.util.TimeUnits;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static cc.redberry.rings.benchmark.Bench.*;

public class Rings_vs_Singular_vs_Mathematica_Factor {
    public static void main(String[] args) throws Exception {
        // warm up
        run(3, 5, 3, 2, 100, Rings.Z, false);
        System.out.println("warmed");
        long[][] timings;
        silent = false;


        for (int nVariables = 3; nVariables <= 7; nVariables++) {
            MATHEMATICA_TIMEOUT_SECONDS = 300;
            doFactorMathematica = nVariables <= 3;

            System.out.println("#variables = " + nVariables);

            int nIterations = 900 / nVariables / nVariables;
            int size = 20;
            int degree = 10;

            timings = run(nVariables, degree, size, 3, nIterations, Rings.Z, false);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z_%s.tsv", nVariables)), timings);
            timings = run(nVariables, degree, size, 3, nIterations, Rings.Z, true);
            writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z_coprime_%s.tsv", nVariables)), timings);

            for (long prime : new long[]{(1 << 19) - 1, 2}) {
                System.out.println("Modulus: " + prime);
                timings = run(nVariables, degree, size, 3, nIterations, Rings.Zp(prime), false);
                writeTimingsTSV(Paths.get(System.getProperty("user.dir"), "rings", "target", String.format("factor_z%s_%s.tsv", prime, nVariables)), timings);
                timings = run(nVariables, degree, size, 3, nIterations, Rings.Zp(prime), true);
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

    static final String log = Paths.get(System.getProperty("user.dir"), "rings", "target", "polynomials.log").toString();

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
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(log))) {
            System.out.println("\tRings\tSingular\tMathematica");
            long[][] timings = new long[nIterations][];
            MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(nVariables, cfRing);
            for (int i = 0; i < nIterations; i++) {
                Stream<MultivariatePolynomial<BigInteger>> ss = IntStream.range(0, nFactors)
                        .mapToObj(__ -> ring.randomElement(degree, size))
                        .filter(__ -> !__.isZero());
                if (cfRing.equals(Rings.Z))
                    ss = ss.map(__ -> __.setRing(Rings.Zp(1000)).setRing(cfRing));// reduce coefficients

                MultivariatePolynomial<BigInteger>[] factors
                        = ss.toArray(MultivariatePolynomial[]::new);

                MultivariatePolynomial<BigInteger> poly = ring.multiply(factors);

                if (coprime)
                    poly.increment();

                //System.out.println(poly.sparsity());

                writer.write(poly.toString());
                writer.newLine();

                long start = System.nanoTime();
                FactorDecomposition<MultivariatePolynomial<BigInteger>> ringsResult = PolynomialMethods.Factor(poly);
                long ringsTime = System.nanoTime() - start;
                //System.out.println(ringsTime);

                ExternalResult singularResult = singularFactor(poly);
                long singularTime = singularResult.nanoTime;
                //System.out.println(singularTime);

                ExternalResult mmaResult = mathematicaFactor(poly);
                long mmaTime = mmaResult.nanoTime;
                //System.out.println(mmaTime );

                timings[i] = new long[]{ringsTime, singularTime, mmaTime};
                if (!silent)
                    System.out.println(i + "\t" + TimeUnits.nanosecondsToString(ringsTime) + "\t" + TimeUnits.nanosecondsToString(singularTime) + "\t" + TimeUnits.nanosecondsToString(mmaTime));
            }
            return timings;
        }
    }
}
