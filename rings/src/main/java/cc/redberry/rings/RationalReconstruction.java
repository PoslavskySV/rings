package cc.redberry.rings;

import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.bigint.BigIntegerUtil;
import cc.redberry.rings.poly.MachineArithmetic;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariateDivision;
import cc.redberry.rings.poly.univar.UnivariateGCD;

/**
 * @since 2.3
 */
public final class RationalReconstruction {
    private RationalReconstruction() {}

    /**
     * Performs a rational number reconstruction. If the answer is not unique, {@code null} is returned.
     *
     * @param n                num * den^(-1) mod modulus
     * @param modulus          the modulus
     * @param numeratorBound   numerator bound
     * @param denominatorBound denominator bound
     */
    public static long[] reconstruct(long n, long modulus,
                                     long numeratorBound,
                                     long denominatorBound) {
        long[] v = {modulus, 0};
        long[] w = {n, 1};

        while (w[0] > numeratorBound) {
            long q = v[0] / w[0];
            long[] z = {v[0] - q * w[0], v[1] - q * w[1]}; // this is safe (no long overflow)
            v = w;
            w = z;
        }
        if (w[1] < 0) {
            w[0] = -w[0];
            w[1] = -w[1];
        }
        if (w[1] <= denominatorBound && MachineArithmetic.gcd(w[0], w[1]) == 1)
            return w;
        return null;
    }

    /**
     * Performs a rational number reconstruction. If the answer is not unique, {@code null} is returned.
     */
    public static BigInteger[] reconstruct(BigInteger n, BigInteger modulus,
                                           BigInteger numeratorBound,
                                           BigInteger denominatorBound) {
        BigInteger[] v = {modulus, BigInteger.ZERO};
        BigInteger[] w = {n, BigInteger.ONE};

        while (w[0].compareTo(numeratorBound) > 0) {
            BigInteger q = v[0].divide(w[0]);
            BigInteger[] z = {v[0].subtract(q.multiply(w[0])), v[1].subtract(q.multiply(w[1]))};
            v = w;
            w = z;
        }
        if (w[1].signum() < 0) {
            w[0] = w[0].negate();
            w[1] = w[1].negate();
        }
        if (w[1].compareTo(denominatorBound) <= 0 && BigIntegerUtil.gcd(w[0], w[1]).isOne())
            return w;
        return null;
    }

    /**
     * Performs a rational number reconstruction. If the answer is not unique, {@code null} is returned.
     *
     * @param n                num * den^(-1) mod modulus
     * @param modulus          the modulus
     * @param numeratorBound   numerator bound
     * @param denominatorBound denominator bound
     */
    public static <Poly extends IUnivariatePolynomial<Poly>> Poly[] reconstruct(Poly n, Poly modulus,
                                                                                int numeratorBound,
                                                                                int denominatorBound) {
        Poly[] v = n.createArray(modulus, n.createZero());
        Poly[] w = n.createArray(n, n.createOne());

        while (w[0].degree() > numeratorBound) {
            Poly q = UnivariateDivision.quotient(v[0], w[0], true);
            Poly[] z = n.createArray(
                    v[0].clone().subtract(q.clone().multiply(w[0])),
                    v[1].clone().subtract(q.clone().multiply(w[1])));
            v = w;
            w = z;
        }
        if (w[1].signumOfLC() < 0) {
            w[0].negate();
            w[1].negate();
        }
        if (w[1].degree() <= denominatorBound && UnivariateGCD.PolynomialGCD(w[0], w[1]).isConstant())
            return w;
        return null;
    }
}
