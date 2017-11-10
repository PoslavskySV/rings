package cc.redberry.rings.benchmark;

import cc.redberry.rings.WithVariables;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.Arrays;
import java.util.stream.Collectors;

class Bench {
    static final class ExternalResult {
        final String result;
        final long nanoTime;

        ExternalResult(String result, long nanoTime) {
            this.result = result;
            this.nanoTime = nanoTime;
        }
    }

    /**
     * How to run Mathematica
     */
    static String MATHEMATICA_CMD_LINE = "-linkmode launch -linkname '\"/Applications/Mathematica.app/Contents/MacOS/MathKernel\" -mathlink'";
    static int MATHEMATICA_TIMEOUT_SECONDS = 3 * 60;

    // mathematica kernel
    static final KernelLink mathKernel;

    static {
        try {
            mathKernel = MathLinkFactory.createKernelLink(MATHEMATICA_CMD_LINE);
            mathKernel.discardAnswer();
        } catch (MathLinkException e) {
            System.out.println("Fatal error opening link: " + e.getMessage());
            throw new RuntimeException(e);
        }
        Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
            @Override
            public void run() {
                mathKernel.abortEvaluation();
                mathKernel.terminateKernel();
                mathKernel.close();
            }
        }));
    }

    /**
     * Path to Singular executable
     */
    public static String SINGULAR = "/Applications/Singular.app/Contents/bin/Singular";

    /**
     * Calculate gcd of two polynomials with Mathematica
     *
     * @return {gcd, timing (in seconds)}
     */
    static ExternalResult mathematicaGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        if (!doFactorMathematica)
            return new ExternalResult("", -1);
        String mResult;
        if (a.isOverFiniteField())
            mResult = mathKernel.evaluateToInputForm(String.format("Timing[TimeConstrained[PolynomialGCD[%s, %s, Modulus->%s], %s]]", a, b, a.ring.cardinality(), MATHEMATICA_TIMEOUT_SECONDS), Integer.MAX_VALUE);
        else
            mResult = mathKernel.evaluateToInputForm(String.format("Timing[TimeConstrained[PolynomialGCD[%s, %s], %s]]", a, b, MATHEMATICA_TIMEOUT_SECONDS), Integer.MAX_VALUE);
        String[] split = mResult.replace("{", "").replace("}", "").split(",");
        return new ExternalResult(split[1], (long) (Double.valueOf(split[0]) * 1000_000_000));
    }

    static boolean doRunSingular = true;

    /**
     * Calculate gcd of two polynomials with Singular
     */
    static ExternalResult singularGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) throws Exception {
        if (!doRunSingular)
            return new ExternalResult("", -1);
        String[] singularCmds = {
                "system(\"--ticks-per-sec\",1000);",
                String.format("ring r = %s,(%s),dp;", a.isOverFiniteField() ? a.coefficientRingCardinality() : "0", String.join(",", WithVariables.defaultVars(a.nVariables))),
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
        Process process = new ProcessBuilder(SINGULAR,
                "-q")
                .redirectErrorStream(true)
                .start();

        process.getOutputStream().write(cmd.getBytes());
        process.getOutputStream().flush();
        process.getOutputStream().close();

        process.waitFor();
        String singularOut = new BufferedReader(new InputStreamReader(process.getInputStream())).lines().reduce((l, r) -> l + r).get();
        String[] split = singularOut.split("SEPARATOR");
        return new ExternalResult(split[0], (long) (Double.valueOf(split[1]) * 1000_000));
    }

    static boolean doFactorMathematica = true;

    /**
     * Calculate gcd of two polynomials with Mathematica
     *
     * @return {gcd, timing (in seconds)}
     */
    static ExternalResult mathematicaFactor(IPolynomial poly) {
        if (!doFactorMathematica)
            return new ExternalResult("", -1);

        if (poly.isOverFiniteField() && (poly instanceof AMultivariatePolynomial && !((AMultivariatePolynomial) poly).isEffectiveUnivariate()))
            return new ExternalResult("", -1);

        String mResult;
        if (poly.isOverField())
            mResult = mathKernel.evaluateToInputForm(String.format("Timing[TimeConstrained[Factor[%s, Modulus->%s], %s]]", poly, poly.coefficientRingCardinality(), MATHEMATICA_TIMEOUT_SECONDS), Integer.MAX_VALUE);
        else
            mResult = mathKernel.evaluateToInputForm(String.format("Timing[TimeConstrained[Factor[%s], %s]]", poly, MATHEMATICA_TIMEOUT_SECONDS), Integer.MAX_VALUE);
        String[] split = mResult.replace("{", "").replace("}", "").split(",");
        return new ExternalResult(split[1], (long) (Double.valueOf(split[0]) * 1000_000_000));
    }

    /**
     * Calculate gcd of two polynomials with Singular
     */
    static ExternalResult singularFactor(AMultivariatePolynomial poly) throws Exception {

        String[] singularCmds = {
                "system(\"--ticks-per-sec\",1000);",
                String.format("ring r = %s,(%s),dp;", poly.isOverFiniteField() ? poly.coefficientRingCardinality() : "0", String.join(",", WithVariables.defaultVars(poly.nVariables))),
                String.format("poly p = %s;", poly),
                "int t = timer;",
                "list g = factorize(p);",
                "int elapsed = timer-t;",
                "print(g);",
                "print(\"SEPARATOR\");",
                "print(elapsed);",
                "exit;"
        };
        String cmd = Arrays.stream(singularCmds).reduce((l, r) -> l + r).get();
        Process process = new ProcessBuilder(SINGULAR,
                "-q")
                .redirectErrorStream(true)
                .start();

        process.getOutputStream().write(cmd.getBytes());
        process.getOutputStream().flush();
        process.getOutputStream().close();

        process.waitFor();
        String singularOut = new BufferedReader(new InputStreamReader(process.getInputStream())).lines().reduce((l, r) -> l + r).get();
        String[] split = singularOut.split("SEPARATOR");
        return new ExternalResult(split[0], (long) (Double.valueOf(split[1]) * 1000_000));
    }


    static void writeTimingsTSV(Path path, long[][] timings) throws IOException {
        Files.write(path, (Iterable<String>)
                () -> Arrays.stream(timings)
                        .map(row -> Arrays.stream(row)
                                .mapToObj(String::valueOf)
                                .collect(Collectors.joining("\t")))
                        .iterator(), StandardOpenOption.CREATE);
    }
}
