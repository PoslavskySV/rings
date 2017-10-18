package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
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
public class MathematicaFactor_UniZpLarge {

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

    public static String mFactor(UnivariatePolynomial<BigInteger> poly) {
        return ml.evaluateToInputForm(String.format("Factor[%s, Modulus->%s]", poly, poly.ring.characteristic()), Integer.MAX_VALUE);
    }

    public static void main(String[] args) {
        // warm up
        run(3, 5, 100, false);
        System.out.println("warmed");
        long[][] timings = run(50, 5, 100, false);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    public static long[][] run(int fDegree, int nFactors, int nIterations, boolean silent) {
        long[][] timings = new long[nIterations][];
        UnivariateRing<UnivariatePolynomial<BigInteger>> ring = Rings.UnivariateRingZp(BigPrimes.nextPrime(Rings.Z.valueOf(1).shiftLeft(100)));
        Well44497b rnd = new Well44497b();
        for (int i = 0; i < nIterations; i++) {
            UnivariatePolynomial<BigInteger>[] factors
                    = IntStream.range(0, nFactors)
                    .mapToObj(__ -> ring.randomElement(fDegree, fDegree, rnd))
                    .filter(__ -> !__.isZero())
                    .toArray(UnivariatePolynomial[]::new);

            UnivariatePolynomial<BigInteger> poly = ring.multiply(factors).monic();

            long start = System.nanoTime();
            FactorDecomposition<UnivariatePolynomial<BigInteger>> rFactors = PolynomialMethods.Factor(poly);
            long ringsTiming = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(ringsTiming));
            }

            start = System.nanoTime();
            String mString = mFactor(poly);
            long mmaTiming = System.nanoTime() - start;
            if (!silent) {
                System.out.println(TimeUnits.nanosecondsToString(mmaTiming));
            }

            int mFactorsNum = mString.split("\\*\\(").length;

            if (rFactors.size() != mFactorsNum) {
                System.out.println(rFactors);
                System.out.println(mString);
                throw new AssertionError();
            }
            if (!silent)
                System.out.println();

            timings[i] = new long[]{ringsTiming, mmaTiming};
        }
        return timings;
    }
}
