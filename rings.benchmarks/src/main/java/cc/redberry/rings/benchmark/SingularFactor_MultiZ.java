package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64;
import cc.redberry.rings.util.TimeUnits;
import org.apache.commons.math3.random.Well44497b;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class SingularFactor_MultiZ {

    public static void main(String[] args) throws Exception {
        // warm up
        run(3, 3, 2, 100, true);
        System.out.println("warmed");
        long[][] timings = run(10, 10, 2, 100, false);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    public static long[][] run(int degree, int size, int nFactors, int nIterations, boolean silent) throws Exception {
        long[][] timings = new long[nIterations][];
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(5, Rings.Zp(17));
        Well44497b rnd = new Well44497b();
        for (int i = 0; i < nIterations; i++) {
            MultivariatePolynomial<BigInteger>[] factors
                    = IntStream.range(0, nFactors)
                    .mapToObj(__ -> ring.randomElement(degree, size, rnd))
                    .filter(__ -> !__.isZero())
                    .toArray(MultivariatePolynomial[]::new);

            MultivariatePolynomial<BigInteger> poly = ring.multiply(factors);

            long start = System.nanoTime();
            FactorDecomposition<MultivariatePolynomial<BigInteger>> rFactors = PolynomialMethods.Factor(poly);
            long ringsTiming = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(ringsTiming));
            }

            String[] mString = singularFactor(poly);
            long mmaTime = Long.valueOf(mString[1]);
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(mmaTime));
            }

            if (!silent)
                System.out.println(i + "\t" + TimeUnits.nanosecondsToString(ringsTiming) + "\t" + TimeUnits.nanosecondsToString(1000 * 1000 * mmaTime));

            timings[i] = new long[]{ringsTiming, mmaTime};
        }
        return timings;
    }


    public static String[] singularFactor(MultivariatePolynomial<BigInteger> poly) throws Exception {
        String[] singularCmds = {
                "system(\"--ticks-per-sec\",1000);",
                "ring r = 17,(x1,x2,x3,x4,x5),dp;",
                String.format("poly p = %s;", poly),
                "int t = timer;",
                "list factors = factorize(p);",
                "int elapsed = timer-t;",
                "print(g);",
                "print(\"SEPARATOR\");",
                "print(elapsed);",
                "exit;"
        };
        String cmd = Arrays.stream(singularCmds).reduce((l, r) -> l + r).get();
        Process process = new ProcessBuilder("/Applications/Singular.app/Contents/bin/Singular",
                "-q", "-c", cmd)
                .redirectErrorStream(true)
                .start();

        process.getOutputStream().write(cmd.getBytes());
        process.getOutputStream().flush();
        process.getOutputStream().close();

        process.waitFor();
        String singularOut = new BufferedReader(new InputStreamReader(process.getInputStream())).lines().reduce((l, r) -> l + r).get();
        return singularOut.split("SEPARATOR");
    }
}
