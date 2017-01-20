package cc.r2.core.polynomial;

import cc.r2.core.polynomial.SmallPolynomialsDivideAndRemainder.InverseModMonomial;

import java.util.ArrayList;
import java.util.Arrays;

import static cc.r2.core.polynomial.MutableLongPoly.createMonomial;
import static cc.r2.core.polynomial.SmallPolynomialArithmetics.polyMod;

public final class SmallPolynomialsModularComposition {
    private SmallPolynomialsModularComposition() {
    }

    /**
     * Returns {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree.
     *
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the prime modulus
     * @return {@code x^{i*modulus} mod polyModulus} for i in {@code [0...degree]}, where {@code degree} is {@code polyModulus} degree
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static ArrayList<MutableLongPoly> xPowers(MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus) {
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
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static ArrayList<MutableLongPoly> polyPowers(MutableLongPoly poly, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, int nIterations) {
        ArrayList<MutableLongPoly> exponents = new ArrayList<>();
        polyPowers(poly, polyModulus, invMod, modulus, nIterations, exponents);
        return exponents;
    }

    /** writes poly^{i} mod polyModulus for i in [0...nIterations] to exponents */
    private static void polyPowers(MutableLongPoly poly, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus, int nIterations, ArrayList<MutableLongPoly> exponents) {
        exponents.add(MutableLongPoly.one());
        MutableLongPoly base = polyMod(poly, polyModulus, invMod, modulus, true);
        exponents.add(base);
        MutableLongPoly prev = base;
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
     * @see #xPowers(MutableLongPoly, InverseModMonomial, long)
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     **/
    public static MutableLongPoly powModulusMod(MutableLongPoly poly,
                                                MutableLongPoly polyModulus,
                                                InverseModMonomial invMod,
                                                long modulus,
                                                ArrayList<MutableLongPoly> xPowers) {
        poly = polyMod(poly, polyModulus, invMod, modulus, true);
        return powModulusMod0(poly, polyModulus, invMod, modulus, xPowers);
    }

    /** doesn't do poly mod polyModulus first */
    private static MutableLongPoly powModulusMod0(MutableLongPoly poly,
                                                  MutableLongPoly polyModulus,
                                                  InverseModMonomial invMod,
                                                  long modulus,
                                                  ArrayList<MutableLongPoly> xPowers) {
        MutableLongPoly res = MutableLongPoly.zero();
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
     * @param main        the main polynomial
     * @param pointPowers precomputed powers of evaluation point {@code point^{i} mod polyModulus}
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @param tBrentKung  Brent-Kung splitting parameter (optimal choice is ~sqrt(main.degree))
     * @return modular composition {@code main(point) mod polyModulus }
     * @see #polyPowers(MutableLongPoly, MutableLongPoly, InverseModMonomial, long, int)
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static MutableLongPoly compositionBrentKung(
            MutableLongPoly main,
            ArrayList<MutableLongPoly> pointPowers,
            MutableLongPoly polyModulus,
            InverseModMonomial invMod,
            long modulus,
            int tBrentKung) {
        if (main.isConstant())
            return main;
        ArrayList<MutableLongPoly> gj = new ArrayList<>();
        for (int i = 0; i <= main.degree; ) {
            int to = i + tBrentKung;
            if (to > (main.degree + 1))
                to = main.degree + 1;
            MutableLongPoly g = MutableLongPoly.create(Arrays.copyOfRange(main.data, i, to));
            gj.add(powModulusMod0(g, polyModulus, invMod, modulus, pointPowers));
            i = to;
        }
        MutableLongPoly pt = pointPowers.get(tBrentKung);
        MutableLongPoly res = MutableLongPoly.zero();
        for (int i = gj.size() - 1; i >= 0; --i)
            res = polyMod(res.multiply(pt, modulus).add(gj.get(i)), polyModulus, invMod, modulus, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated using Brent & Kung algorithm for modular composition.
     *
     * @param poly        the main polynomial
     * @param point       evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @return modular composition {@code main(point) mod polyModulus }
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static MutableLongPoly compositionBrentKung(MutableLongPoly poly, MutableLongPoly point, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus) {
        if (poly.isConstant())
            return poly;
        int t = (int) Math.sqrt(poly.degree);
        ArrayList<MutableLongPoly> hPowers = polyPowers(point, polyModulus, invMod, modulus, t);
        return compositionBrentKung(poly, hPowers, polyModulus, invMod, modulus, t);
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus } calculated with plain Horner scheme.
     *
     * @param poly        the main polynomial
     * @param point       evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param invMod      precomputed inversed {@code rev[polyModulus]}
     * @param modulus     the modulus
     * @return modular composition {@code main(point) mod polyModulus }
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static MutableLongPoly compositionHorner(MutableLongPoly poly, MutableLongPoly point, MutableLongPoly polyModulus, InverseModMonomial invMod, long modulus) {
        if (poly.isConstant())
            return poly;
        MutableLongPoly res = MutableLongPoly.zero();
        for (int i = poly.degree; i >= 0; --i)
            res = polyMod(res.multiply(point, modulus).addMonomial(poly.data[i], 0, modulus), polyModulus, invMod, modulus, false);
        return res;
    }

    /**
     * Returns modular composition {@code poly(point) mod polyModulus }
     *
     * @param poly        the main polynomial
     * @param point       evaluation point
     * @param polyModulus the monic polynomial modulus
     * @param modulus     the modulus
     * @return modular composition {@code main(point) mod polyModulus }
     * @see #polyPowers(MutableLongPoly, MutableLongPoly, InverseModMonomial, long, int)
     * @see SmallPolynomialsDivideAndRemainder#fastDivisionPreConditioning(MutableLongPoly, long)
     */
    public static MutableLongPoly composition(MutableLongPoly poly, MutableLongPoly point, MutableLongPoly polyModulus, long modulus) {
        return compositionBrentKung(poly, point, polyModulus, SmallPolynomialsDivideAndRemainder.fastDivisionPreConditioning(point, modulus), modulus);
    }
}
