package cc.redberry.rings.benchmark;

import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MathematicaGCD_Z {


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

    public static String mGCD(MultivariatePolynomial<BigInteger> a, MultivariatePolynomial<BigInteger> b) {
        return ml.evaluateToInputForm(String.format("PolynomialGCD[%s, %s]", a, b), Integer.MAX_VALUE);
    }

    public static void main(String[] args) {
        // warm up
        run(5, 5, 100);
        long[][] timings = run(10, 10, 100);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    public static long[][] run(int degree, int size, int nIterations) {
        long[][] timings = new long[nIterations][];
        MultivariateRing<MultivariatePolynomial<BigInteger>> ring = Rings.MultivariateRing(3, Rings.Z);
        for (int i = 0; i < nIterations; i++) {
            MultivariatePolynomial<BigInteger>
                    a = ring.randomElement(degree, size),
                    b = ring.randomElement(degree, size),
                    gcd = ring.randomElement(degree, size);

            a = a.multiply(gcd);
            b = b.multiply(gcd);

            long start = System.nanoTime();
            MultivariatePolynomial<BigInteger> rGCD = PolynomialMethods.PolynomialGCD(a, b);
            long ringsTime = System.nanoTime() - start;


            start = System.nanoTime();
            MultivariatePolynomial<BigInteger> mGCD = ring.parse(mGCD(a, b));
            long mmaTime = System.nanoTime() - start;

            if (rGCD.signumOfLC() != mGCD.signumOfLC())
                rGCD.negate();

            if (!rGCD.primitivePart().equals(mGCD.primitivePart())) {
                System.out.println(rGCD);
                System.out.println(mGCD);
                System.out.println(gcd);
                throw new AssertionError();
            }

            timings[i] = new long[]{ringsTime, mmaTime};
        }
        return timings;
    }
}
