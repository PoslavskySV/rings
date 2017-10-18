package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomialZp64;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MathematicaGCD_Z2 {


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

    public static String mGCD(MultivariatePolynomialZp64 a, MultivariatePolynomialZp64 b) {
        return ml.evaluateToInputForm(String.format("TimeConstrained[PolynomialGCD[%s, %s, Modulus->%s], 3]", a, b, a.ring.modulus), Integer.MAX_VALUE);
    }

    public static void main(String[] args) {
        // warm up
        run(5, 3, 1000);
        System.out.println("warmed");
        silent = false;
        long[][] timings = run(10, 10, 100);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    static boolean silent = true;

    public static long[][] run(int degree, int size, int nIterations) {
        long[][] timings = new long[nIterations][];
        String[] vars = {"a", "b", "c"};
        MultivariateRing<MultivariatePolynomialZp64> ring = Rings.MultivariateRingZp64(3, 2);
        for (int i = 0; i < 100; i++) {
            MultivariatePolynomialZp64
                    a = ring.randomElement(degree, size),
                    b = ring.randomElement(degree, size),
                    gcd = ring.randomElement(degree, size);

            a = a.multiply(gcd);
            b = b.multiply(gcd);

            long start = System.nanoTime();
            MultivariatePolynomialZp64 rGCD = PolynomialMethods.PolynomialGCD(a, b);
            long ringsTime = System.nanoTime() - start;

            start = System.nanoTime();
            String mString = mGCD(a, b);
            long mmaTime = System.nanoTime() - start;
            MultivariatePolynomialZp64 mGCD = null;
            try {
                mGCD = ring.parse(mString, vars);
            } catch (Exception e) {}

            if (mGCD != null && !rGCD.monic().equals(mGCD.monic())) {
                System.out.println(rGCD);
                System.out.println(mGCD);
                System.out.println(gcd);
                throw new AssertionError();
            }

            timings[i] = new long[]{ringsTime, mmaTime};
            if (!silent)
                System.out.println(i);
        }
        return timings;
    }
}
