package cc.redberry.rings.benchmark;

import cc.redberry.rings.poly.FactorDecomposition;
import cc.redberry.rings.poly.PolynomialMethods;
import cc.redberry.rings.poly.univar.UnivariatePolynomialZp64;
import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkException;
import com.wolfram.jlink.MathLinkFactory;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public class MathematicaFactor_UniZp_2 {

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
        run(6, true);
        System.out.println("warmed");
        System.out.println(String.format("degree\tRings\tMathematica"));
        double[][] timings = run(12, false);
        System.out.println(Arrays.deepToString(timings).replace("[", "{").replace("]", "}"));
    }

    static int dummy = 0;
    static int nIterations = 3;

    public static double[][] run(int maxL, boolean silent) {
        ArrayList<double[]> timings = new ArrayList<>();
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
                    FactorDecomposition<UnivariatePolynomialZp64> rFactors = PolynomialMethods.Factor(poly);
                    dummy += rFactors.signum();
                }
                double ringsTime = (System.nanoTime() - start) / 1000_000_000.0;

                start = System.nanoTime();
                for (int i = 0; i < nIterations; i++) {
                    String mString = mFactor(poly);
                    dummy += mString.length();
                }
                double mmaTime = (System.nanoTime() - start) / 1000_000_000.0;

                timings.add(new double[]{degree, ringsTime, mmaTime});
                if (!silent)
                    System.out.println(String.format("%s\t%5.4f\t%5.4f", degree, ringsTime, mmaTime));
            }
        return timings.toArray(new double[timings.size()][]);
    }
}
