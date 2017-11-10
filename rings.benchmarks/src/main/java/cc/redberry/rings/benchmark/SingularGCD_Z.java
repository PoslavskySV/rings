package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64;
import cc.redberry.rings.util.TimeUnits;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class SingularGCD_Z {

    public static void main(String[] args) throws Exception {
        // warm up
        run(5, 3, 1000);
        System.out.println("warmed");
        silent = false;
        long[][] timings = run(30, 50, 100);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    static boolean silent = true;

    public static long[][] run(int degree, int size, int nIterations) throws Exception {
        long[][] timings = new long[nIterations][];
        String[] vars = {"a", "b", "c"};
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(3, Rings.Z);
        for (int i = 0; i < 100; i++) {
            MultivariatePolynomial<BigInteger>
                    a = ring.randomElement(degree, size),
                    b = ring.randomElement(degree, size),
                    gcd = ring.randomElement(degree, size);

            a = a.multiply(gcd);
            b = b.multiply(gcd);

            String[] singular = singularGCD(a, b);
            long mmaTime = Long.valueOf(singular[1]);
            System.out.println(mmaTime);

            long start = System.nanoTime();
            MultivariatePolynomial<BigInteger> rGCD = PolynomialMethods.PolynomialGCD(a, b);
            long ringsTime = System.nanoTime() - start;
            System.out.println(ringsTime);


            assert rGCD != null;

            timings[i] = new long[]{ringsTime, mmaTime};
            if (!silent)
                System.out.println(i + "\t" + TimeUnits.nanosecondsToString(ringsTime) + "\t" + TimeUnits.nanosecondsToString(1000 * 1000 * mmaTime));
        }
        return timings;
    }


    public static String[] singularGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) throws Exception {
        String[] singularCmds = {
                "system(\"--ticks-per-sec\",1000);",
                "ring r = 17,(x,y,z),dp;",
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
                "-q")
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
