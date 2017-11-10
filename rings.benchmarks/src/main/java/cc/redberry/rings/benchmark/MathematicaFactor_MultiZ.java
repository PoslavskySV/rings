package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64;
import cc.redberry.rings.util.TimeUnits;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;
import org.apache.commons.math3.random.Well44497b;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MathematicaFactor_MultiZ {

    static final KernelLink ml;

    static {
        try {
            ml = MathLinkFactory.createKernelLink("-linkmode launch -linkname '\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\" -mathlink'");
            ml.discardAnswer();
        } catch (MathLinkException e) {
            System.out.println("Fatal error opening link: " + e.getMessage());
            throw new RuntimeException(e);
        }

    }

    public static String mFactor(MultivariatePolynomial<BigInteger> poly) {
        return ml.evaluateToInputForm(String.format("Factor[%s]", poly), Integer.MAX_VALUE);
    }

    public static void main(String[] args) {
        // warm up
        run(3, 3, 2, 100, true);
        System.out.println("warmed");
        long[][] timings = run(5, 10, 2, 100, false);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    @SuppressWarnings("unchecked")
    public static long[][] run(int degree, int size, int nFactors, int nIterations, boolean silent) {
        long[][] timings = new long[nIterations][];
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(4, Rings.Z);
        Well44497b rnd = new Well44497b();
        for (int i = 0; i < nIterations; i++) {
            MultivariatePolynomial<BigInteger>[] factors
                    = IntStream.range(0, nFactors)
                    .mapToObj(__ -> ring.randomElement(degree, size))
                    .filter(__ -> !__.isZero())
                    .toArray(MultivariatePolynomial[]::new);

            MultivariatePolynomial<BigInteger> poly = ring.multiply(factors);

            long start = System.nanoTime();
            FactorDecomposition<MultivariatePolynomial<BigInteger>> rFactors = PolynomialMethods.Factor(poly);
            long ringsTiming = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(ringsTiming));
            }

            start = System.nanoTime();
            String mString = mFactor(poly);
            long mmaTimings = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(mmaTimings));
            }

            int mFactorsNum = mString.split("\\*\\(").length;

            if (rFactors.size() != mFactorsNum) {
//                System.out.println();
//                System.out.println(rFactors.size());
//                System.out.println(mFactorsNum);
//                System.out.println();
//                throw new AssertionError();
            }
            if (!silent)
                System.out.println();

            timings[i] = new long[]{ringsTiming, mmaTimings};
        }
        return timings;
    }


    public static String[] singularFactor(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b) throws Exception {
        String[] singularCmds = {
                "system(\"--ticks-per-sec\",1000);",
                "ring r = 2,(x,y,z),dp;",
                String.format("poly poly1 = %s;", a),
                String.format("poly poly2 = %s;", b),
                "int t = timer;",
                "poly g = gcd(poly1, poly2);",
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

        process.waitFor();
        String singularOut = new BufferedReader(new InputStreamReader(process.getInputStream())).lines().reduce((l, r) -> l + r).get();
        return singularOut.split("SEPARATOR");
    }
}
