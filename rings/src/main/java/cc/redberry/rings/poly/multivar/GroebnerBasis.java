package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.map.hash.TLongObjectHashMap;

import java.util.*;

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
        removeRedundant(groebner);
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
                    if (!shareVariablesQ(fi.lt(), fj.lt()))
                        continue;
                    Poly syzygy = MultivariateDivision.remainder(syzygy(fi, fj), groebnerArray);
                    if (syzygy.isZero())
                        continue;
                    temporary.add(syzygy);
                }
            }
            if (temporary.isEmpty()) {
                minimizeGroebnerBases(groebner);
                removeRedundant(groebner);
                return groebner;
            }
            groebner.addAll(temporary);
            groebner.sort(ltOrder);
            removeRedundant(groebner);
            // groebner.sort(ltOrder); <- this make things slower...
        }
    }

    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    void setMonomialOrder(List<Poly> list, Comparator<DegreeVector> order) {
        for (int i = 0; i < list.size(); i++) {
            Poly p = list.get(i);
            if (!order.equals(p.ordering))
                list.set(i, p.setOrdering(order));
        }
    }

    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    void canonicalize(List<Poly> list) {
        list.forEach(Poly::monic);
        list.sort(Comparable::compareTo);
    }

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm.
     *
     * NOTE: this algorithm is optimized to perform fast with GREVLEX order
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB3(List<Poly> ideal, Comparator<DegreeVector> monomialOrder) {
        return BuchbergerGB3(ideal, monomialOrder, normalSelectionStrategy(ideal.get(0).ordering), (i, j) -> false);
    }

    public interface MinimizationStrategy {
        boolean doMinimize(int previousSize, int currentSize);
    }

    public static MinimizationStrategy NO_MINIMIZATION = (prev, curr) -> false;

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm.
     *
     * NOTE: this algorithm is optimized to perform fast with GREVLEX order
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB3(List<Poly> ideal,
                             Comparator<DegreeVector> monomialOrder,
                             Comparator<SyzygyPair> selectionStrategy,
                             MinimizationStrategy minimizationStrategy) {
        setMonomialOrder(ideal, monomialOrder);

        Poly factory = ideal.get(0);

        // stack of GB candidates
        List<Poly> groebner = new ArrayList<>(ideal);
        // remove redundant elements from the basis
        removeRedundant(groebner);

        // sort polynomials in basis in the ascending order, to achieve faster divisions
        // this gives 2x performance boost for GREVLEX
        Comparator<Poly> polyOrder = (a, b) -> monomialOrder.compare(a.lt(), b.lt());
        groebner.sort(polyOrder);

        // pairs are ordered like (0, 1), (0, 2), ... (0, N), (1, 2), (1, 3), ...
        TreeSet<SyzygyPair<Term, Poly>> sPairs = new TreeSet<>(selectionStrategy);
        TLongObjectHashMap<SyzygyPair<Term, Poly>> ijPairs = new TLongObjectHashMap<>();
        for (int i = 0; i < groebner.size() - 1; i++)
            for (int j = i + 1; j < groebner.size(); ++j) {
                SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(i, j, groebner);
                sPairs.add(sPair);
                ijPairs.put(pack(i, j), sPair);
            }

        // cache array used in divisions (little performance improvement actually)
        Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));

        int sizeAfterReduction = groebner.size();
        main:
        while (!sPairs.isEmpty()) {
            // pick up the pair with "smallest" syzygy
            SyzygyPair<Term, Poly> pair = sPairs.pollFirst();
            ijPairs.remove(pair.index());
            int
                    i = pair.i,
                    j = pair.j;

            Poly
                    fi = pair.fi,
                    fj = pair.fj;

            if (!shareVariablesQ(fi.lt(), fj.lt()))
                continue;

            // test criterion
            // l.c.m. of lts
            int[] lcm = ArraysUtil.max(fi.lt().exponents, fj.lt().exponents);
            for (int k = 0; k < groebner.size(); ++k) {
                if (groebner.get(k) == null)
                    continue;
                if (k == i || k == j)
                    continue;
                if (ijPairs.contains(pack(i, k)) || ijPairs.contains(pack(j, k)))
                    continue;
                if (dividesQ(lcm, groebner.get(k).lt().exponents))
                    continue main;
            }

            Poly syzygy = MultivariateDivision.remainder(pair.computeSyzygy(), groebnerArray);
            if (syzygy.isZero())
                continue;

            groebner.add(syzygy);
            // recompute array
            groebnerArray = groebner.stream().filter(Objects::nonNull).toArray(factory::createArray);
            //Arrays.sort(groebnerArray, (a, b) -> a.ordering.compare(a.lt(), b.lt()));

            for (int k = 0; k < groebner.size() - 1; k++)
                if (groebner.get(k) != null) {
                    SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(k, groebner.size() - 1, groebner);
                    sPairs.add(sPair);
                    ijPairs.put(sPair.index(), sPair);
                }

            if (minimizationStrategy.doMinimize(sizeAfterReduction, groebner.size())) {
                reduceAndMinimizeGroebnerBases(groebner, sPairs, ijPairs, sizeAfterReduction);
                sizeAfterReduction = groebner.size();
            }
        }

        // remove null entries
        groebner.removeAll(Collections.<Poly>singleton(null)); // batch remove all nulls

        // minimize Groebner basis & canonicalize it
        minimizeGroebnerBases(groebner);
        removeRedundant(groebner);
        return groebner;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void reduceAndMinimizeGroebnerBases(List<Poly> groebner,
                                        TreeSet<SyzygyPair<Term, Poly>> sPairs,
                                        TLongObjectHashMap<SyzygyPair<Term, Poly>> ijPairs,
                                        int from) {
        for (int i = from; i < groebner.size(); i++) {
            // this are newly added elements of Groebner basis
            Poly fi = groebner.get(i);
            if (fi == null)
                continue;

            for (int j = 0; j < i; ++j) { // search for fj which is divisible by fi
                Poly fj = groebner.get(j);
                if (fj == null)
                    continue;

                // proceed only if new syzygy can reduce fk
                if (!MultivariateDivision.nontrivialQuotientQ(fj, fi))
                    continue;

                groebner.remove(j);
                Poly reduced = remainder0(fj, groebner);
                if (reduced.equals(fj))
                    continue;

                if (reduced.isZero()) // remove element
                    reduced = null;

                groebner.add(j, reduced);

                if (!fj.equals(reduced)) {
                    if (reduced == null) {
                        // remove all pairs with k
                        for (int l = 0; l < groebner.size(); l++)
                            if (l != j && groebner.get(l) != null) {
                                SyzygyPair<Term, Poly> sPair = ijPairs.remove(pack(l, j));
                                if (sPair != null)
                                    sPairs.remove(sPair);
                            }
                    } else
                        // update all pairs with k
                        for (int l = 0; l < groebner.size(); l++)
                            if (l != j && groebner.get(l) != null) {
                                SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(l, j, groebner);
                                sPairs.add(sPair);
                                ijPairs.put(sPair.index(), sPair);
                            }
                }
            }
        }
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder0(Poly dividend, List<Poly> dividers) {
        Poly[] dividersArr = dividers.stream().filter(Objects::nonNull).toArray(dividend::createArray);
        //Arrays.sort(dividersArr, (a, b) -> a.ordering.compare(a.lt(), b.lt()));
        return MultivariateDivision.remainder(dividend, dividersArr);
    }

    static Comparator<SyzygyPair> normalSelectionStrategy(Comparator<DegreeVector> monomialOrder) {
        return (a, b) -> {
            int c = monomialOrder.compare(a.syzygyGamma, b.syzygyGamma);
            if (c != 0)
                return c;
            c = Integer.compare(a.j, b.j);
            if (c != 0)
                return c;
            return Integer.compare(a.i, b.i);
        };
    }

    static Comparator<SyzygyPair> withSugar(Comparator<SyzygyPair> initial) {
        return (sa, sb) -> {
            int c = Integer.compare(sa.sugar, sb.sugar);
            if (c != 0)
                return c;
            return initial.compare(sa, sb);
        };
    }

    public static final class SyzygyPair<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        final int i, j;
        final Poly fi, fj;
        final Term syzygyGamma;
        // the sugar
        final int sugar;

        SyzygyPair(int i, int j, List<Poly> basis) {
            this(i, j, basis.get(i), basis.get(j));
        }

        SyzygyPair(int i, int j, Poly fi, Poly fj) {
            if (i > j) {
                int s = i; i = j; j = s;
                Poly fs = fi; fi = fj; fj = fs;
            }
            assert i != j;
            assert fi != fj;
            this.i = i; this.fi = fi;
            this.j = j; this.fj = fj;
            this.syzygyGamma = fi.createTermWithUnitCoefficient(ArraysUtil.max(fi.multidegree(), fj.multidegree()));
            this.sugar = Math.max(
                    fi.degreeSum() - fi.lt().totalDegree,
                    fj.degreeSum() - fj.lt().totalDegree) + syzygyGamma.totalDegree;
        }

        long index() { return pack(i, j);}

        Poly computeSyzygy() {
            return syzygy(syzygyGamma, fi, fj);
        }
    }

    private static boolean dividesQ(int[] dividend, int[] divider) {
        for (int i = 0; i < dividend.length; i++)
            if (dividend[i] < divider[i])
                return false;
        return true;
    }

    private static long pack(int i, int j) {
        if (i > j)
            return pack(j, i);
        return ((long) j) << 32 | (long) i;
    }

    private static boolean shareVariablesQ(DegreeVector a, DegreeVector b) {
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
    void removeRedundant(List<Poly> basis) {
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
        return syzygy(a.createTermWithUnitCoefficient(ArraysUtil.max(a.multidegree(), b.multidegree())), a, b);
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Term xGamma, Poly a, Poly b) {
        Poly
                aReduced = a.clone().multiply(a.divideOrNull(xGamma, a.lt())),
                bReduced = b.clone().multiply(b.divideOrNull(xGamma, b.lt())),
                syzygy = aReduced.subtract(bReduced);
        return syzygy;
    }
}
