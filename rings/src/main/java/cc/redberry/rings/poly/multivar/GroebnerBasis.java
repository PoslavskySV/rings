package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.set.hash.TLongHashSet;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;

/**
 * Basic Groebner basis.
 *
 * @since 1.0
 */
public final class GroebnerBasis {
    private GroebnerBasis() {}

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> ideal) {
        Poly factory = ideal.get(0);
        Comparator<DegreeVector> ordering = factory.ordering;
        Comparator<Poly> ltOrder = (a, b) -> ordering.compare(a.lt(), b.lt()); // <- this gives 2x speed up

        List<Poly> groebner = new ArrayList<>(ideal);
        removeRedundand(groebner);
        groebner.sort(ltOrder);
        List<Poly> temporary = new ArrayList<>();
        while (true) {
            temporary.clear();
            Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));
            for (int i = 0; i < groebner.size() - 1; ++i) {
                for (int j = i + 1; j < groebner.size(); ++j) {
                    Poly
                            fi = groebner.get(i),
                            fj = groebner.get(j);
                    if (!shareVariables(fi.lt(), fj.lt()))
                        continue;
                    Poly syzygy = MultivariateDivision.remainder(syzygy(fi, fj), groebnerArray);
                    if (syzygy.isZero())
                        continue;
                    temporary.add(syzygy);
                }
            }
            if (temporary.isEmpty()) {
                minimizeGroebnerBases(groebner);
                removeRedundand(groebner);
                return groebner;
            }
            groebner.addAll(temporary);
            groebner.sort(ltOrder);
            removeRedundand(groebner);
            // groebner.sort(ltOrder); <- this make things slower...
        }
    }

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB2(List<Poly> ideal) {
        Poly factory = ideal.get(0);
        Comparator<DegreeVector> ordering = factory.ordering;

        List<Poly> groebner = new ArrayList<>(ideal);
        removeRedundand(groebner);

        groebner.sort((a, b) -> ordering.compare(a.lt(), b.lt())); // <- this gives 2x speed up

        // fill all initial pairs
        TLongHashSet pairs = new TLongHashSet(ideal.size() * ideal.size());
        for (int i = 0; i < groebner.size() - 1; i++)
            for (int j = i + 1; j < groebner.size(); ++j)
                pairs.add(pack(i, j));

        Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));
        int its = 0;
        int nontriv = 0;
        main:
        while (!pairs.isEmpty()) {
            ++its;
            long el = pairs.iterator().next();
            pairs.remove(el);

            int i = (int) (0xFFFFFFFF & el);
            int j = (int) (el >> 32);

            Poly
                    fi = groebner.get(i),
                    fj = groebner.get(j);

            if (!shareVariables(fi.lt(), fj.lt()))
                continue;

            // criterion
            int[] ltlcm = ArraysUtil.max(fi.lt().exponents, fj.lt().exponents);
            for (int k = 0; k < groebner.size(); ++k) {
                if (groebner.get(k) == null)
                    continue;
                if (k == i || k == j)
                    continue;
                if (pairs.contains(packOrdered(i, k)) || pairs.contains(packOrdered(j, k)))
                    continue;
                if (dividesQ(ltlcm, groebner.get(k).lt().exponents))
                    continue main;
            }

            Poly syzygy = MultivariateDivision.remainder(syzygy(fi, fj), groebnerArray);
            if (syzygy.isZero())
                continue;

            ++nontriv;
            groebner.add(syzygy);
            // recompute array
            groebnerArray = groebner.stream().filter(Objects::nonNull).toArray(factory::createArray);

            pairs.ensureCapacity(groebner.size());
            for (int k = 0; k < groebner.size() - 1; k++)
                if (groebner.get(k) != null)
                    pairs.add(pack(k, groebner.size() - 1));

            // remove redundant elements from GB
            for (int k = 0; k < groebner.size() - 1; ++k) {
                Poly fk = groebner.get(k);
                if (fk == null)
                    continue;

                // proceed only if new syzygy can reduce fk
                if (!MultivariateDivision.nontrivialQuotientQ(fk, syzygy))
                    continue;

                groebner.remove(k);
                Poly rem = MultivariateDivision.remainder(fk, groebner.stream().filter(Objects::nonNull).toArray(factory::createArray));
                if (rem.isZero())
                    rem = null;
                groebner.add(k, rem);

                if (!fk.equals(rem)) {
                    if (rem == null) {
                        // remove all pairs with k
                        for (int l = 0; l < groebner.size(); l++)
                            if (l != k && groebner.get(l) != null)
                                pairs.remove(packOrdered(l, k));
                    } else
                        // update all pairs with k
                        for (int l = 0; l < groebner.size(); l++)
                            if (l != k && groebner.get(l) != null)
                                pairs.add(packOrdered(l, k));
                }
            }
        }

        System.out.println("nitss = " + its);
        System.out.println("nontr = " + nontriv);

        // remove null entries
        while (groebner.remove(null)) ;

        minimizeGroebnerBases(groebner);
        removeRedundand(groebner);
        return groebner;
    }

    private static boolean dividesQ(int[] dividend, int[] divider) {
        for (int i = 0; i < dividend.length; i++)
            if (dividend[i] < divider[i])
                return false;
        return true;
    }

    private static long pack(int i, int j) {
        assert i < j;
        return ((long) j) << 32 | (long) i;
    }

    private static long packOrdered(int i, int j) { return i < j ? pack(i, j) : pack(j, i); }

    private static boolean shareVariables(DegreeVector a, DegreeVector b) {
        for (int i = 0; i < a.exponents.length; i++)
            if (a.exponents[i] != 0 && b.exponents[i] != 0)
                return true;
        return false;
    }

    /**
     * Minimizes Groebner basis
     */
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void minimizeGroebnerBases(List<Poly> basis) {
        outer:
        for (int i = basis.size() - 1; i >= 1; --i) {
            for (int j = i - 1; j >= 0; --j) {
                Poly pi = basis.get(i), pj = basis.get(j);
                if (pi.lt().divisibleBy(pj.lt())) {
                    basis.remove(i);
                    continue outer;
                }
                if (pj.lt().divisibleBy(pi.lt())) {
                    basis.remove(j);
                    --i;
                    continue;
                }
            }
        }
        for (Poly el : basis)
            el.monic();
    }

    /**
     * Computes reduced Groebner basis
     */
    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void removeRedundand(List<Poly> basis) {
        for (int i = 0, size = basis.size(); i < size; ++i) {
            Poly el = basis.remove(i);
            Poly r = MultivariateDivision.remainder(el, basis);
            if (r.isZero()) {
                --i;
                --size;
            } else
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
