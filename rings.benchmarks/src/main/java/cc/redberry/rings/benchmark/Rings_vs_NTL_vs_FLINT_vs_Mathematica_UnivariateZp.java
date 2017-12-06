package cc.redberry.rings.benchmark;

import cc.redberry.rings.benchmark.Bench.ExternalResult;
import cc.redberry.rings.poly.PolynomialFactorDecomposition;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.util.TimeUnits;

import java.io.*;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import static cc.redberry.rings.benchmark.Bench.mathematicaFactor;

public class Rings_vs_NTL_vs_FLINT_vs_Mathematica_UnivariateZp {
    //NTL and FLINT data are generated from C++ binaries!

    public static void main(String[] args) throws IOException {
        // warm up
        run(6, true);
        System.out.println("warmed");
        run(12, false);
    }

    static int dummy = 0;
    static int nIterations = 3;

    public static long[][] run(int maxL, boolean silent) throws IOException {
        Path path = Paths.get(System.getProperty("user.dir"), "rings.benchmarks", "results", "factor_uni_rings_mathematica.tsv");
        try(FileOutputStream writer = new FileOutputStream(path.toFile())) {
            System.out.println("Degree\tRings\tMathematica");
            ArrayList<long[]> timings = new ArrayList<>();
            for (long l = 4; l <= maxL; l += 1)
                for (long idx = 0; idx < 10; idx++) {
                    long degree = (1 << l) + idx * (((1 << (l + 1)) - (1 << l)) / 10);

                    UnivariatePolynomialZp64 poly = UnivariatePolynomialZp64.zero(17);
                    for (int j = 0; j <= degree; j++)
                        poly.set(j, j);

                    poly.set(0, 1);
                    poly.monic();


                    long start = System.nanoTime();
                    for (int i = 0; i < nIterations; i++) {
                        PolynomialFactorDecomposition<UnivariatePolynomialZp64> rFactors = PolynomialMethods.Factor(poly);
                        dummy += rFactors.signum();
                    }
                    long ringsTime = (System.nanoTime() - start);

                    long mmaTime = 0;
                    for (int i = 0; i < nIterations; i++) {
                        ExternalResult mathematicaResult = mathematicaFactor(poly);
                        dummy += mathematicaResult.nanoTime;
                        mmaTime += mathematicaResult.nanoTime;
                    }

                    writer.write(String.format("%s\t%s\t%s\n", degree, ringsTime, mmaTime).getBytes());
                    timings.add(new long[]{degree, ringsTime, mmaTime});
                    if (!silent)
                        System.out.println(String.format("%s\t%s\t%s", degree, TimeUnits.nanosecondsToString(ringsTime), TimeUnits.nanosecondsToString(mmaTime)));
                }
            return timings.toArray(new long[timings.size()][]);
        }
    }
}
