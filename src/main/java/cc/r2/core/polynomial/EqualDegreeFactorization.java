package cc.r2.core.polynomial;

import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well1024a;

import static cc.r2.core.polynomial.DivideAndRemainder.fastDivisionPreConditioning;
import static cc.r2.core.polynomial.PolynomialGCD.PolynomialGCD;

/**
 * Created by poslavsky on 21/01/2017.
 */
public final class EqualDegreeFactorization {
    private EqualDegreeFactorization() {}

    private static final RandomGenerator CZRANDOM = new Well1024a();

    public MutablePolynomial CantorZassenhaus(MutablePolynomial poly, long modulus, int d) {
        MutablePolynomial a = RandomPolynomials.randomMonicPoly(poly.degree - 1, modulus, CZRANDOM);
        MutablePolynomial gcd = PolynomialGCD(a, poly, modulus);
        if (!gcd.isConstant())
            return gcd;

        InverseModMonomial invMod = fastDivisionPreConditioning(poly, modulus);

        // we have to power a^(p^d-1)/2
        // (p^d-1)/2 = (1 + p + p^2 + ... + p^(d-1)) * (p-1)/2
        return null;
    }

}
