package cc.r2.core.polynomial;

import cc.r2.core.polynomial.DivideAndRemainder.InverseModMonomial;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.MutablePolynomial.createMonomial;
import static cc.r2.core.polynomial.PolynomialArithmetics.polyMod;

public final class ModularComposition {
    private ModularComposition() {
    }

    /**
     * Returns {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree.
     *
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the prime modulus
     * @return {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static ArrayList<MutablePolynomial> xPowers(MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus) {
        return polyPowers(createMonomial(1, (int) modulus), polyModulus, invMod, modulus, polyModulus.degree);
    }

    /**
     * Returns {@code poly^{i} mod polyModulus} for i in {@code [0...nIterations]}
     *
     * @param poly        the polynomial
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the prime modulus
     * @return {@code poly^{i} mod polyModulus} for i in {@code [0...nIterations]}
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static ArrayList<MutablePolynomial> polyPowers(MutablePolynomial poly, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, int nIterations) {
        ArrayList<MutablePolynomial> exponents = new ArrayList<>();
        polyPowers(poly, polyModulus, invMod, modulus, nIterations, exponents);
        return exponents;
    }

    /** writes poly^{i} mod polyModulus for i in [0...nIterations] to exponents */
    static void polyPowers(MutablePolynomial poly, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus, int nIterations, ArrayList<MutablePolynomial> exponents) {
        exponents.add(MutablePolynomial.one());
        MutablePolynomial base = polyMod(poly, polyModulus, invMod, modulus, true);
        exponents.add(base);
        MutablePolynomial prev = base;
        for (int i = 0; i < nIterations; i++)
            exponents.add(prev = polyMod(prev.clone().multiply(base, modulus), polyModulus, invMod, modulus, false));
    }

    /**
     * Returns {@code poly^modulus mod polyModulus} using precomputed monomial powers {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree(poly)]}
     *
     * @param poly        the polynomial
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @param xPowers     precomputed monomial powers {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree(poly)]}
     * @return {@code poly^modulus mod polyModulus}
     * @see #xPowers(MutablePolynomial, InverseModMonomial, long)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     **/
    public static MutablePolynomial powModulusMod(MutablePolynomial poly,
                                                  MutablePolynomial polyModulus,
                                                  InverseModMonomial invMod,
                                                  long modulus,
                                                  ArrayList<MutablePolynomial> xPowers) {
        poly = polyMod(poly, polyModulus, invMod, modulus, true);
        return powModulusMod0(poly, polyModulus, invMod, modulus, xPowers);
    }

    /** doesn't do poly mod polyModulus first */
    private static MutablePolynomial powModulusMod0(MutablePolynomial poly,
                                                    MutablePolynomial polyModulus,
                                                    InverseModMonomial invMod,
                                                    long modulus,
                                                    ArrayList<MutablePolynomial> xPowers) {
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] == 0)
                continue;
            res.addMul(xPowers.get(i), poly.data[i], modulus);
        }
        return polyMod(res, polyModulus, invMod, modulus, false);
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated using Brent & Kung algorithm for modular composition.
     *
     * @param main        the polynomial
     * @param pointPowers precomputed powers of evaluation point {@code point^{i} mod polyModulus}
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @param tBrentKung  Brent-Kung splitting parameter (optimal choice is ~sqrt(main.degree))
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see #polyPowers(MutablePolynomial, MutablePolynomial, InverseModMonomial, long, int)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static MutablePolynomial compositionBrentKung(
            MutablePolynomial main,
            ArrayList<MutablePolynomial> pointPowers,
            MutablePolynomial polyModulus,
            InverseModMonomial invMod,
            long modulus,
            int tBrentKung) {
        if (main.isConstant())
            return main;
        ArrayList<MutablePolynomial> gj = new ArrayList<>();
        for (int i = 0; i <= main.degree; ) {
            int to = i + tBrentKung;
            if (to > (main.degree + 1))
                to = main.degree + 1;
            MutablePolynomial g = MutablePolynomial.create(Arrays.copyOfRange(main.data, i, to));
            gj.add(powModulusMod0(g, polyModulus, invMod, modulus, pointPowers));
            i = to;
        }
        MutablePolynomial pt = pointPowers.get(tBrentKung);
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = gj.size() - 1; i >= 0; --i)
            res = polyMod(res.multiply(pt, modulus).add(gj.get(i)), polyModulus, invMod, modulus, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated using Brent & Kung algorithm for modular composition.
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static MutablePolynomial compositionBrentKung(MutablePolynomial poly, MutablePolynomial point, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus) {
        if (poly.isConstant())
            return poly;
        int t = (int) Math.sqrt(poly.degree);
        ArrayList<MutablePolynomial> hPowers = polyPowers(point, polyModulus, invMod, modulus, t);
        return compositionBrentKung(poly, hPowers, polyModulus, invMod, modulus, t);
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated with plain Horner scheme.
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static MutablePolynomial compositionHorner(MutablePolynomial poly, MutablePolynomial point, MutablePolynomial polyModulus, InverseModMonomial invMod, long modulus) {
        if (poly.isConstant())
            return poly;
        MutablePolynomial res = MutablePolynomial.zero();
        for (int i = poly.degree; i >= 0; --i)
            res = polyMod(res.multiply(point, modulus).addMonomial(poly.data[i], 0, modulus), polyModulus, invMod, modulus, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus}. Brent & Kung algorithm used
     * ({@link #compositionBrentKung(MutablePolynomial, MutablePolynomial, MutablePolynomial, InverseModMonomial, long)}).
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param modulus     the modulus
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see #polyPowers(MutablePolynomial, MutablePolynomial, InverseModMonomial, long, int)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomial, long)
     */
    public static MutablePolynomial composition(MutablePolynomial poly, MutablePolynomial point, MutablePolynomial polyModulus, long modulus) {
        return compositionBrentKung(poly, point, polyModulus, DivideAndRemainder.fastDivisionPreConditioning(point, modulus), modulus);
    }
}
