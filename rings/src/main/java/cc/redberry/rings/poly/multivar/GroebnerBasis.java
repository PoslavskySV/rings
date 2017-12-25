package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.map.hash.TLongObjectHashMap;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Groebner basis computation.
 *
 * @since 1.0
 */
public final class GroebnerBasis {
    private GroebnerBasis() {}


    /* **************************************** Util methods ************************************************ */

    /** Set monomial order inplace */
    @SuppressWarnings("unchecked")
    static <Poly extends AMultivariatePolynomial<?, Poly>>
    void setMonomialOrder(List<Poly> list, Comparator<DegreeVector> order) {
        for (int i = 0; i < list.size(); i++) {
            Poly p = list.get(i);
            if (!order.equals(p.ordering))
                list.set(i, p.setOrdering(order));
        }
    }

    /** Make each poly monic and sort list */
    static <Poly extends AMultivariatePolynomial<?, Poly>>
    List<Poly> canonicalize(List<Poly> list) {
        // todo replace monic with `to common denominator` for Q[X]
        list.forEach(Poly::monic);
        // fixme sort with respect to minimal total degree first
        list.sort(Comparable::compareTo);
        return list;
    }

    /** whether monomial order is graded */
    static boolean isGradedOrder(Comparator<DegreeVector> monomialOrder) {
        return monomialOrder == MonomialOrder.GREVLEX || monomialOrder == MonomialOrder.GRLEX;
    }


    /**
     * Computes Groebner basis (minimized and reduced) of a given ideal represented as a list of generators.
     *
     * @param generators    generators of the ideal
     * @param monomialOrder monomial order
     * @return Groebner basis
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> GroebnerBasis(List<Poly> generators,
                             Comparator<DegreeVector> monomialOrder) {
        return BuchbergerGB(generators, monomialOrder);
    }

    private static <Poly extends AMultivariatePolynomial<?, Poly>>
    boolean containsNonZeroConstant(List<Poly> generators) {
        return generators.stream().anyMatch(p -> p.isConstant() && !p.isZero());
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> prepareGenerators(List<Poly> generators, Comparator<DegreeVector> monomialOrder) {
        generators = new ArrayList<>(generators);
        setMonomialOrder(generators, monomialOrder);
        // remove redundant elements from the basis
        removeRedundant(generators);
        generators = generators.stream().filter(p -> !p.isZero()).collect(Collectors.toList());
        if (containsNonZeroConstant(generators))
            return Collections.singletonList(generators.get(0).createOne());
        if (generators.get(0).nVariables == 1)
            // univariate case
            return canonicalize(Collections.singletonList(MultivariateGCD.PolynomialGCD(generators)));
        return generators;
    }

    /* ************************************** Buchberger algorithm ********************************************** */

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm. It uses normal strategy
     * for graded orders and sugar strategy for lexicographic.
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> ideal, Comparator<DegreeVector> monomialOrder) {
        Comparator<SyzygyPair> selectionStrategy = normalSelectionStrategy(ideal.get(0).ordering);
        // fixme use sugar always?
        if (!isGradedOrder(monomialOrder))
            // add sugar for non-graded orders
            selectionStrategy = withSugar(selectionStrategy);
        return BuchbergerGB(ideal, monomialOrder, selectionStrategy, NO_MINIMIZATION);
    }

    /** Strategy used to reduce and minimize basis in the intermediate steps of Buchberger algorithm */
    public interface MinimizationStrategy {
        /**
         * true means "yes, do minimization and reduction", false means "just keep all generators as is"
         *
         * @param previousSize size of generators list on the previous invocation of minimization
         * @param currentSize  current size of generators list
         */
        boolean doMinimize(int previousSize, int currentSize);
    }

    /** no any minimization at intermediate steps, just keep all track of generators as is */
    public static MinimizationStrategy NO_MINIMIZATION = (prev, curr) -> false;

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm.
     *
     * @param generators           generators of the ideal
     * @param monomialOrder        monomial order to use
     * @param selectionStrategy    critical pair selection strategy
     * @param minimizationStrategy how to minimize Groebner basis at intermediate steps
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> generators,
                            Comparator<DegreeVector> monomialOrder,
                            Comparator<SyzygyPair> selectionStrategy,
                            MinimizationStrategy minimizationStrategy) {
        Poly factory = generators.get(0);
        // stack of GB candidates
        List<Poly> groebner = prepareGenerators(generators, monomialOrder);
        if (groebner.size() == 1)
            return groebner;

        // sort polynomials in basis to achieve faster divisions
        Comparator<Poly> polyOrder = isGradedOrder(monomialOrder)
                ? (a, b) -> monomialOrder.compare(a.lt(), b.lt())
                : new EcartComparator<>();
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

        // cache size of basis after each minimization
        int sizeAfterMinimization = groebner.size();
        main:
        while (!sPairs.isEmpty()) {
            // pick up the pair with "smallest" syzygy
            SyzygyPair<Term, Poly> pair = sPairs.pollFirst();
            ijPairs.remove(pair.hash());
            int
                    i = pair.i,
                    j = pair.j;

            Poly
                    fi = pair.fi,
                    fj = pair.fj;

            if (!shareVariablesQ(fi.lt(), fj.lt()))
                // don't test for relatively prime lts
                continue;

            // test criterion (Gebauer-Moller)
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
            // don't tail reduce
            // syzygy = syzygy.ltAsPoly().add(MultivariateDivision.remainder(syzygy.subtractLt(), groebnerArray));
            if (syzygy.isZero())
                continue;

            if (syzygy.isConstant())
                // ideal = ring
                return Collections.singletonList(factory.createOne());

            groebner.add(syzygy);
            // recompute array
            groebnerArray = groebner.stream().filter(Objects::nonNull).toArray(factory::createArray);
            // don't sort here, not practical actually
            // Arrays.sort(groebnerArray, polyOrder);

            for (int k = 0; k < groebner.size() - 1; k++)
                if (groebner.get(k) != null) {
                    SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(k, groebner.size() - 1, groebner);
                    sPairs.add(sPair);
                    ijPairs.put(sPair.hash(), sPair);
                }

            if (minimizationStrategy.doMinimize(sizeAfterMinimization, groebner.size())) {
                reduceAndMinimizeGroebnerBases(groebner, sPairs, ijPairs, sizeAfterMinimization);
                sizeAfterMinimization = groebner.size();
            }
        }

        // batch remove all nulls
        groebner.removeAll(Collections.<Poly>singleton(null));

        // minimize Groebner basis
        minimizeGroebnerBases(groebner);
        // speed up final reduction
        groebner.sort(polyOrder);
        // reduce Groebner basis
        removeRedundant(groebner);
        // canonicalize Groebner basis
        canonicalize(groebner);

        return groebner;
    }

    private static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void reduceAndMinimizeGroebnerBases(List<Poly> generators,
                                        TreeSet<SyzygyPair<Term, Poly>> sPairs,
                                        TLongObjectHashMap<SyzygyPair<Term, Poly>> ijPairs,
                                        int from) {
        for (int i = from; i < generators.size(); i++) {
            // this are newly added elements of Groebner basis
            Poly fi = generators.get(i);
            if (fi == null)
                continue;

            for (int j = 0; j < i; ++j) { // search for fj which is divisible by fi
                Poly fj = generators.get(j);
                if (fj == null)
                    continue;

                // proceed only if new syzygy can reduce fk
                if (!MultivariateDivision.nontrivialQuotientQ(fj, fi))
                    continue;

                generators.remove(j);
                Poly reduced = remainder0(fj, generators);
                if (reduced.equals(fj))
                    continue;

                if (reduced.isZero()) // remove element
                    reduced = null;

                generators.add(j, reduced);

                if (!fj.equals(reduced)) {
                    if (reduced == null) {
                        // remove all pairs with k
                        for (int l = 0; l < generators.size(); l++)
                            if (l != j && generators.get(l) != null) {
                                SyzygyPair<Term, Poly> sPair = ijPairs.remove(pack(l, j));
                                if (sPair != null)
                                    sPairs.remove(sPair);
                            }
                    } else
                        // update all pairs with k
                        for (int l = 0; l < generators.size(); l++)
                            if (l != j && generators.get(l) != null) {
                                SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(l, j, generators);
                                sPairs.add(sPair);
                                ijPairs.put(sPair.hash(), sPair);
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

    /**
     * Normal selection strategy: chose syzygy with the less lcm(fi.lt(), fj.lt()) with respect to monomialOrder
     */
    public static Comparator<SyzygyPair> normalSelectionStrategy(Comparator<DegreeVector> monomialOrder) {
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

    /**
     * Add sugar to selection strategy: pick syzygy with less sugar first, break tie with initial strategy
     */
    public static Comparator<SyzygyPair> withSugar(Comparator<SyzygyPair> initial) {
        return (sa, sb) -> {
            int c = Integer.compare(sa.sugar, sb.sugar);
            if (c != 0)
                return c;
            return initial.compare(sa, sb);
        };
    }

    private static final class EcartComparator<Poly extends AMultivariatePolynomial<?, Poly>>
            implements Comparator<Poly> {
        @Override
        public int compare(Poly a, Poly b) {
            int c = Integer.compare(a.ecart(), b.ecart());
            if (c != 0)
                return c;
            return a.ordering.compare(a.lt(), b.lt());
        }
    }

    /** Syzygy pair */
    public static final class SyzygyPair<
            Term extends DegreeVector<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> {
        /** Positions of polynomials {@code fi} and {@code fj} in the list of generators */
        final int i, j;
        /** Polynomials to form a syzygy */
        final Poly fi, fj;
        /** {@code lcm(fi.lt(), fj.lt())} */
        final Term syzygyGamma;
        /** The sugar */
        final int sugar;

        public SyzygyPair(int i, int j, List<Poly> generators) {
            this(i, j, generators.get(i), generators.get(j));
        }

        public SyzygyPair(int i, int j, Poly fi, Poly fj) {
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

        long hash() { return pack(i, j);}

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

    static boolean shareVariablesQ(DegreeVector a, DegreeVector b) {
        for (int i = 0; i < a.exponents.length; i++)
            if (a.exponents[i] != 0 && b.exponents[i] != 0)
                return true;
        return false;
    }

    /**
     * Minimizes Groebner basis
     */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

    /** Computes syzygy of given polynomials */
    public static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Poly a, Poly b) {
        return syzygy(a.createTermWithUnitCoefficient(ArraysUtil.max(a.multidegree(), b.multidegree())), a, b);
    }

    static <Term extends DegreeVector<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Term xGamma, Poly a, Poly b) {
        Poly
                aReduced = a.clone().multiply(a.divideOrNull(xGamma, a.lt())),
                bReduced = b.clone().multiply(b.divideOrNull(xGamma, b.lt())),
                syzygy = aReduced.subtract(bReduced);
        return syzygy;
    }
}
