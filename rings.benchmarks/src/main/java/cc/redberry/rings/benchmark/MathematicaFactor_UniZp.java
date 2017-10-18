package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import cc.redberry.rings.primes.BigPrimes;
import cc.redberry.rings.util.TimeUnits;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;
import org.apache.commons.math3.random.Well44497b;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MathematicaFactor_UniZp {

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

    public static String mFactor(UnivariatePolynomialZp64 poly) {
        return ml.evaluateToInputForm(String.format("Factor[%s, Modulus->%s]", poly, poly.ring.modulus), Integer.MAX_VALUE);
    }

    public static void main(String[] args) {
        // warm up
        run(5, 5, 10000, true);
        System.out.println("warmed");
        long[][] timings = run(50, 5, 100, false);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    public static long[][] run(int fDegree, int nFactors, int nIterations, boolean silent) {
        long[][] timings = new long[nIterations][];
        UnivariateRing<UnivariatePolynomialZp64> ring = Rings.UnivariateRingZp64(BigPrimes.nextPrime(1L << 15));
        Well44497b rnd = new Well44497b();
        for (int i = 0; i < nIterations; i++) {
            UnivariatePolynomialZp64[] factors
                    = IntStream.range(0, nFactors)
                    .mapToObj(__ -> ring.randomElement(fDegree, fDegree, rnd))
                    .filter(__ -> !__.isZero())
                    .toArray(UnivariatePolynomialZp64[]::new);

            UnivariatePolynomialZp64 poly = ring.multiply(factors).monic();

            long start = System.nanoTime();
            FactorDecomposition<UnivariatePolynomialZp64> rFactors = PolynomialMethods.Factor(poly);
            long ringsTime = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(ringsTime));
            }

            start = System.nanoTime();
            String mString = mFactor(poly);
            long mmaTime = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(mmaTime));
            }

            int mFactorsNum = mString.split("\\*\\(").length;

            if (rFactors.size() != mFactorsNum) {
                System.out.println(rFactors);
                System.out.println(mString);
                throw new AssertionError();
            }
            if (!silent)
                System.out.println();

            timings[i] = new long[]{ringsTime, mmaTime};
        }
        return timings;
    }
}
