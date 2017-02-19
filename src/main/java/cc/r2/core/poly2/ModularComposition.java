package cc.r2.core.poly2;


import cc.r2.core.poly2.DivideAndRemainder.InverseModMonomial;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.poly2.PolynomialArithmetics.polyMod;


public final class ModularComposition {
    private ModularComposition() {
    }

    /**
     * Returns {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree.
     *
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @return {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static ArrayList<MutablePolynomialMod> xPowers(MutablePolynomialMod polyModulus, InverseModMonomial invMod) {
        return polyPowers(MutablePolynomialMod.createMonomial(polyModulus.modulus, 1, (int) polyModulus.modulus), polyModulus, invMod, polyModulus.degree);
    }

    /**
     * Returns {@code poly^{i} mod polyModulus} for i in {@code [0...nIterations]}
     *
     * @param poly        the polynomial
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @return {@code poly^{i} mod polyModulus} for i in {@code [0...nIterations]}
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static ArrayList<MutablePolynomialMod> polyPowers(MutablePolynomialMod poly, MutablePolynomialMod polyModulus, InverseModMonomial invMod, int nIterations) {
        ArrayList<MutablePolynomialMod> exponents = new ArrayList<>();
        polyPowers(poly, polyModulus, invMod, nIterations, exponents);
        return exponents;
    }

    /** writes poly^{i} mod polyModulus for i in [0...nIterations] to exponents */
    static void polyPowers(MutablePolynomialMod poly, MutablePolynomialMod polyModulus, InverseModMonomial invMod, int nIterations, ArrayList<MutablePolynomialMod> exponents) {
        exponents.add(poly.createOne());
        MutablePolynomialMod base = polyMod(poly, polyModulus, invMod, true);
        exponents.add(base);
        MutablePolynomialMod prev = base;
        for (int i = 0; i < nIterations; i++)
            exponents.add(prev = polyMod(prev.clone().multiply(base), polyModulus, invMod, false));
    }

    /**
     * Returns {@code poly^modulus mod polyModulus} using precomputed monomial powers {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree(poly)]}
     *
     * @param poly        the polynomial
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param xPowers     precomputed monomial powers {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree(poly)]}
     * @return {@code poly^modulus mod polyModulus}
     * @see #xPowers(MutablePolynomialMod, InverseModMonomial)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     **/
    public static MutablePolynomialMod powModulusMod(MutablePolynomialMod poly,
                                                     MutablePolynomialMod polyModulus,
                                                     InverseModMonomial invMod,
                                                     ArrayList<MutablePolynomialMod> xPowers) {
        poly = polyMod(poly, polyModulus, invMod, true);
        return powModulusMod0(poly, polyModulus, invMod, xPowers);
    }

    /** doesn't do poly mod polyModulus first */
    private static MutablePolynomialMod powModulusMod0(MutablePolynomialMod poly,
                                                       MutablePolynomialMod polyModulus,
                                                       InverseModMonomial invMod,
                                                       ArrayList<MutablePolynomialMod> xPowers) {
        MutablePolynomialMod res = poly.createZero();
        for (int i = poly.degree; i >= 0; --i) {
            if (poly.data[i] == 0)
                continue;
            res.addMul(xPowers.get(i), poly.data[i]);
        }
        return polyMod(res, polyModulus, invMod, false);
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated using Brent & Kung algorithm for modular composition.
     *
     * @param poly        the polynomial
     * @param pointPowers precomputed powers of evaluation point {@code point^{i} mod polyModulus}
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param tBrentKung  Brent-Kung splitting parameter (optimal choice is ~sqrt(main.degree))
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see #polyPowers(MutablePolynomialMod, MutablePolynomialMod, InverseModMonomial, int)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod compositionBrentKung(
            MutablePolynomialMod poly,
            ArrayList<MutablePolynomialMod> pointPowers,
            MutablePolynomialMod polyModulus,
            InverseModMonomial invMod,
            int tBrentKung) {
        if (poly.isConstant())
            return poly;
        ArrayList<MutablePolynomialMod> gj = new ArrayList<>();
        for (int i = 0; i <= poly.degree; ) {
            int to = i + tBrentKung;
            if (to > (poly.degree + 1))
                to = poly.degree + 1;
            MutablePolynomialMod g = poly.createFromArray(Arrays.copyOfRange(poly.data, i, to));
            gj.add(powModulusMod0(g, polyModulus, invMod, pointPowers));
            i = to;
        }
        MutablePolynomialMod pt = pointPowers.get(tBrentKung);
        MutablePolynomialMod res = poly.createZero();
        for (int i = gj.size() - 1; i >= 0; --i)
            res = polyMod(res.multiply(pt).add(gj.get(i)), polyModulus, invMod, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated using Brent & Kung algorithm for modular composition.
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod compositionBrentKung(MutablePolynomialMod poly, MutablePolynomialMod point, MutablePolynomialMod polyModulus, InverseModMonomial invMod) {
        if (poly.isConstant())
            return poly;
        int t = (int) Math.sqrt(poly.degree);
        ArrayList<MutablePolynomialMod> hPowers = polyPowers(point, polyModulus, invMod, t);
        return compositionBrentKung(poly, hPowers, polyModulus, invMod, t);
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated with plain Horner scheme.
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod compositionHorner(MutablePolynomialMod poly, MutablePolynomialMod point, MutablePolynomialMod polyModulus, InverseModMonomial invMod) {
        if (poly.isConstant())
            return poly;
        MutablePolynomialMod res = poly.createZero();
        for (int i = poly.degree; i >= 0; --i)
            res = polyMod(res.multiply(point).addMonomial(poly.data[i], 0), polyModulus, invMod, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus}. Brent & Kung algorithm used
     * ({@link #compositionBrentKung(MutablePolynomialMod, MutablePolynomialMod, MutablePolynomialMod, InverseModMonomial)}).
     *
     * @param poly        the polynomial
     * @param point       the evaluation point
     * @param polyModulus the monic polynomial modulus
     * @return modular composition {@code poly(point) mod polyModulus }
     * @see #polyPowers(MutablePolynomialMod, MutablePolynomialMod, InverseModMonomial, int)
     * @see DivideAndRemainder#fastDivisionPreConditioning(MutablePolynomialMod)
     */
    public static MutablePolynomialMod composition(MutablePolynomialMod poly, MutablePolynomialMod point, MutablePolynomialMod polyModulus) {
        return compositionBrentKung(poly, point, polyModulus, DivideAndRemainder.fastDivisionPreConditioning(point));
    }
}
