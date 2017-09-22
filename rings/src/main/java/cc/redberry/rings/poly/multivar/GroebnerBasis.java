package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;

import java.util.ArrayList;
import java.util.List;

/**
 * Basic Groebner basis.
 *
 * @since 1.0
 */
public final class GroebnerBasis {
    private GroebnerBasis() {}

    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> ideal) {
        Poly factory = ideal.get(0);
        List<Poly> groebner = new ArrayList<>(ideal);
        List<Poly> temporary = new ArrayList<>();
        while (true) {
            temporary.clear();
            Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));
            for (int i = 0; i < groebner.size() - 1; ++i) {
                for (int j = i + 1; j < groebner.size(); ++j) {
                    Poly syzygy = MultivariateDivision.remainder(syzygy(groebner.get(i), groebner.get(j)), groebnerArray);
                    if (!syzygy.isZero())
                        temporary.add(syzygy);
                }
            }
            if (temporary.isEmpty()) {
                minimizeGroebnerBases(groebner);
                reduceGroebnerBases(groebner);
                return groebner;
            }
            groebner.addAll(temporary);
        }
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void minimizeGroebnerBases(List<Poly> basis) {
        outer:
        for (int i = basis.size() - 1; i >= 1; --i) {
            for (int j = i - 1; j >= 0; --j) {
                Poly pi = basis.get(i), pj = basis.get(j);
                if (pi.lt().dividesQ(pj.lt())) {
                    basis.remove(i);
                    continue outer;
                }
                if (pj.lt().dividesQ(pi.lt())) {
                    basis.remove(j);
                    --i;
                    continue;
                }
            }
        }
        for (Poly el : basis)
            el.monic();
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void reduceGroebnerBases(List<Poly> basis) {
        for (int i = 0, size = basis.size(); i < size; ++i) {
            Poly el = basis.remove(i);
            Poly r = MultivariateDivision.remainder(el, basis);
            basis.add(i, r);
        }
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Poly a, Poly b) {
        int[] gamma = ArraysUtil.max(a.multidegree(), b.multidegree());
        Term xGamma = a.createTermWithUnitCoefficient(gamma);
        Poly
                aReduced = a.clone().multiply(a.divideOrNull(xGamma, a.lt())),
                bReduced = b.clone().multiply(b.divideOrNull(xGamma, b.lt())),
                syzygy = aReduced.subtract(bReduced);
        return syzygy;
    }
}
