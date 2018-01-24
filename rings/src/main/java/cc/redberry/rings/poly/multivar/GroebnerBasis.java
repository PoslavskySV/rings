package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.IntegersZp64;
import cc.redberry.rings.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;

import java.util.*;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Groebner basis computation.
 *
 * @since 1.0
 */
public final class GroebnerBasis {
    private GroebnerBasis() {}

    /**
     * Computes Groebner basis (minimized and reduced) of a given ideal represented as a list of generators.
     *
     * @param generators    generators of the ideal
     * @param monomialOrder monomial order
     * @return Groebner basis
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> GroebnerBasis(List<Poly> generators,
                             Comparator<DegreeVector> monomialOrder) {
        return BuchbergerGB(generators, monomialOrder);
    }

    /* **************************************** Common methods ************************************************ */

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

    /** set monomial order, monicize, sort etc. */
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> prepareGenerators(List<Poly> generators, Comparator<DegreeVector> monomialOrder) {
        // todo monicize

        // clone generators
        generators = generators.stream().map(Poly::clone).collect(Collectors.toList());
        // set monomial order
        setMonomialOrder(generators, monomialOrder);

        Poly factory = generators.get(0);
        if (factory.nVariables == 1)
            // univariate case
            return canonicalize(Collections.singletonList(MultivariateGCD.PolynomialGCD(generators)));

        // remove zeroes
        generators.removeIf(Poly::isZero);
        if (generators.isEmpty()) {
            // empty ideal
            generators.add(factory.createZero());
            return generators;
        }

        // remove redundant elements from the basis
        removeRedundant(generators);

        // remove zeroes (may occur after reduction)
        generators.removeIf(Poly::isZero);
        if (generators.isEmpty()) {
            // empty ideal
            generators.add(factory.createZero());
            return generators;
        }

        if (generators.stream().anyMatch(Poly::isConstant))
            // contains non zero constant => ideal == whole ring
            return Collections.singletonList(factory.createOne());

        return generators;
    }


    /** Minimizes Groebner basis. The input should be a Groebner basis */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void minimizeGroebnerBases(List<Poly> basis) {
        outer:
        for (int i = basis.size() - 1; i >= 1; --i) {
            for (int j = i - 1; j >= 0; --j) {
                Poly pi = basis.get(i), pj = basis.get(j);
                if (pi.lt().dvDivisibleBy(pj.lt())) {
                    basis.remove(i);
                    continue outer;
                }
                if (pj.lt().dvDivisibleBy(pi.lt())) {
                    basis.remove(j);
                    --i;
                    continue;
                }
            }
        }
        for (Poly el : basis)
            el.monic();
    }

    /** Computes reduced Groebner basis */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

    /** Abstract critical pair: used with different Poly type for Buchberger and F4 algorithms */
    public static class SyzygyPair<Term extends AMonomial<Term>, Poly extends MonomialSetView<Term>> {
        /** Positions of polynomials {@code fi} and {@code fj} in the list of generators */
        final int i, j;
        /** Polynomials to form a syzygy */
        final Poly fi, fj;
        /** {@code lcm(fi.lt(), fj.lt())} */
        final DegreeVector syzygyGamma;
        /** The sugar */
        final int sugar;

        SyzygyPair(int i, int j, Poly fi, Poly fj) {
            if (i > j) {
                int s = i; i = j; j = s;
                Poly fs = fi; fi = fj; fj = fs;
            }
            assert i != j && fi != fj;
            this.i = i; this.j = j;
            this.fi = fi; this.fj = fj;
            this.syzygyGamma = lcm(fi.lt(), fj.lt());
            this.sugar = Math.max(
                    fi.degreeSum() - fi.lt().totalDegree,
                    fj.degreeSum() - fj.lt().totalDegree) + syzygyGamma.totalDegree;
        }

        SyzygyPair(int i, int j, List<Poly> generators) {
            this(i, j, generators.get(i), generators.get(j));
        }

        /** syzygy degree */
        final int degree() { return syzygyGamma.totalDegree;}
    }

    /** Computes syzygy of given polynomials */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Poly a, Poly b) {
        return syzygy(a.monomialAlgebra.create(lcm(a.multidegree(), b.multidegree())), a, b);
    }

    /** Computes syzygy of given polynomials */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(SyzygyPair<Term, Poly> sPair) {
        return syzygy(sPair.fi.monomialAlgebra.create(sPair.syzygyGamma), sPair.fi, sPair.fj);
    }

    /** Computes syzygy of given polynomials */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Term lcm, Poly a, Poly b) {
        IMonomialAlgebra<Term> mAlgebra = a.monomialAlgebra;
        Poly
                aReduced = a.clone().multiply(mAlgebra.divideExact(lcm, a.lt())),
                bReduced = b.clone().multiply(mAlgebra.divideExact(lcm, b.lt())),
                syzygy = aReduced.subtract(bReduced);
        return syzygy;
    }

    /**
     * Data structure which holds current set of critical pairs. This set may be just a Set&lt;SyzygyPair&gt; (TreeSet
     * actually), or some more complicated structure (map of sets, e.g. TreeMap&lt;Integer,
     * TreeSet&lt;SyzygyPair&lt;Term, Poly&gt;&gt;&gt; of degree -&gt; set of pairs with this degree) for graded
     * algorithms (F4, homogeneous Buchberger etc).
     */
    private interface SyzygySet<Term extends AMonomial<Term>, Poly extends MonomialSetView<Term>> {
        /** add new syzygy to set */
        void add(SyzygyPair<Term, Poly> sPair);

        /** remove syzygy from set */
        void remove(SyzygyPair<Term, Poly> sPair);

        /** is empty */
        boolean isEmpty();

        /** remove a bunch of syzygies to process */
        Collection<SyzygyPair<Term, Poly>> getAndRemoveNextBunch();

        /** a view of all sets of syzygies */
        Collection<TreeSet<SyzygyPair<Term, Poly>>> allSets();
    }

    /**
     * A plain set of critical pairs --- wrapper around TreeSet&lt;SyzygyPair&lt;Term, Poly&gt;&gt;. Elements are
     * ordered according to the comparator of TreeSet instance.
     */
    private static final class SyzygyTreeSet<Term extends AMonomial<Term>, Poly extends MonomialSetView<Term>>
            implements SyzygySet<Term, Poly> {
        final TreeSet<SyzygyPair<Term, Poly>> sPairs;

        SyzygyTreeSet(TreeSet<SyzygyPair<Term, Poly>> sPairs) { this.sPairs = sPairs; }

        @Override
        public void add(SyzygyPair<Term, Poly> sPair) { sPairs.add(sPair); }

        @Override
        public void remove(SyzygyPair<Term, Poly> sPair) { sPairs.remove(sPair); }

        @Override
        public boolean isEmpty() { return sPairs.isEmpty(); }

        /**
         * Returns the single element from the top of TreeSet
         */
        @Override
        public Collection<SyzygyPair<Term, Poly>> getAndRemoveNextBunch() {return Collections.singleton(sPairs.pollFirst());}

        /**
         * Returns the single TreeSet (this)
         */
        @Override
        public Collection<TreeSet<SyzygyPair<Term, Poly>>> allSets() { return Collections.singleton(sPairs);}
    }

    /**
     * Graded set of syzygies: critical pairs are collected in the sets each one containing critical pairs with the same
     * degree. Wrapper around TreeMap&lt;Integer, TreeSet&lt;SyzygyPair&lt;Term, Poly&gt;&gt;&gt;
     */
    private static final class GradedSyzygyTreeSet<Term extends AMonomial<Term>, Poly extends MonomialSetView<Term>>
            implements SyzygySet<Term, Poly> {
        final TreeMap<Integer, TreeSet<SyzygyPair<Term, Poly>>> sPairs;
        final Comparator<SyzygyPair> selectionStrategy;

        GradedSyzygyTreeSet(TreeMap<Integer, TreeSet<SyzygyPair<Term, Poly>>> sPairs,
                            Comparator<SyzygyPair> selectionStrategy) {
            this.sPairs = sPairs;
            this.selectionStrategy = selectionStrategy;
        }

        @Override
        public void add(SyzygyPair<Term, Poly> sPair) {
            sPairs.computeIfAbsent(sPair.degree(), __ -> new TreeSet<>(selectionStrategy)).add(sPair);
        }

        @Override
        public void remove(SyzygyPair<Term, Poly> sPair) {
            TreeSet<SyzygyPair<Term, Poly>> set = sPairs.get(sPair.degree());
            if (set != null)
                set.remove(sPair);
        }

        @Override
        public boolean isEmpty() { return sPairs.isEmpty(); }

        /**
         * Returns all syzygies with the minimal degree
         */
        @Override
        public Collection<SyzygyPair<Term, Poly>> getAndRemoveNextBunch() {return sPairs.pollFirstEntry().getValue();}

        /**
         * Returns sets of syzygies each one with fixed degree
         */
        @Override
        public Collection<TreeSet<SyzygyPair<Term, Poly>>> allSets() { return sPairs.values();}
    }

    /**
     * Generic Buchberger criteria via Gebauer-Moller installation. Algorithm from page 230 of Becker & Weispfenning,
     * 1993. Note: redundant generators in the {@code basis} are marked as null.
     */
    private static <Term extends AMonomial<Term>, Poly extends MonomialSetView<Term>>
    void updateBasis(List<Poly> basis,
                     SyzygySet<Term, Poly> sPairs,
                     Poly newElement) {
        // array of lcm( lt(fi), lt(newElement) ) <- cache once for performance
        int[][] lcm = new int[basis.size()][];

        Term newLeadTerm = newElement.lt();
        for (int i = 0; i < basis.size(); i++) {
            Poly fi = basis.get(i);
            if (fi == null)
                continue;
            lcm[i] = lcm(fi.lt().exponents, newLeadTerm.exponents);
        }

        // first indices of new critical pairs to add
        TIntArrayList pairsToAdd = new TIntArrayList();
        // find new critical pairs that should be definitely added
        filter:
        for (int iIndex = 0; iIndex < basis.size(); iIndex++) {
            Poly fi = basis.get(iIndex);
            if (fi == null)
                continue;
            if (!shareVariablesQ(fi.lt(), newLeadTerm)) {
                pairsToAdd.add(iIndex); // add disjoint elements (will be removed in the next step)
                continue;
            }

            // ruling out redundant Buchberger triplets: 1st pass
            for (int i = 0; i < pairsToAdd.size(); i++) {
                int jIndex = pairsToAdd.get(i);
                Poly fj = basis.get(jIndex);
                assert fj != null;
                if (dividesQ(lcm[jIndex], lcm[iIndex]))
                    continue filter;
            }

            // ruling out redundant Buchberger triplets: 2st pass
            for (int jIndex = iIndex + 1; jIndex < basis.size(); ++jIndex) {
                Poly fj = basis.get(jIndex);
                if (fj == null)
                    continue;
                if (dividesQ(lcm[jIndex], lcm[iIndex]))
                    continue filter;
            }

            // no any redundant Buchberger triplets found -> add new critical pair
            pairsToAdd.add(iIndex);
        }

        // now rule out disjoint elements
        for (int i = pairsToAdd.size() - 1; i >= 0; --i)
            if (!shareVariablesQ(basis.get(pairsToAdd.get(i)).lt(), newLeadTerm))
                pairsToAdd.removeAt(i);

        // ruling out redundant Buchberger triplets from the old set of critical pairs
        for (TreeSet<SyzygyPair<Term, Poly>> c : sPairs.allSets())
            // remove redundant critical pairs
            c.removeIf(sPair -> dividesQ(newLeadTerm, sPair.syzygyGamma)
                    && !Arrays.equals(lcm[sPair.i] != null ? lcm[sPair.i] : lcm(sPair.fi.lt().exponents, newLeadTerm.exponents), sPair.syzygyGamma.exponents)
                    && !Arrays.equals(lcm[sPair.j] != null ? lcm[sPair.j] : lcm(sPair.fj.lt().exponents, newLeadTerm.exponents), sPair.syzygyGamma.exponents));

        // now add new element to the basis
        basis.add(newElement);

        // update set of critical pairs with pairsToAdd
        int oldSize = basis.size() - 1;
        for (int i = 0; i < pairsToAdd.size(); i++) {
            int iIndex = pairsToAdd.get(i);
            assert basis.get(iIndex) != null;
            sPairs.add(new SyzygyPair<>(iIndex, oldSize, basis));
        }

        // remove old basis elements that are now redundant
        // note: do that only after update of critical pair set
        for (int iIndex = 0; iIndex < oldSize; ++iIndex) {
            Poly fi = basis.get(iIndex);
            if (fi == null)
                continue;
            if (dividesQ(newLeadTerm, fi.lt()))
                basis.set(iIndex, null);
        }
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
     * Add sugar to selection strategy: pick syzygy with less sugar first, break tie with the initial strategy
     */
    public static Comparator<SyzygyPair> withSugar(Comparator<SyzygyPair> initial) {
        return (sa, sb) -> {
            int c = Integer.compare(sa.sugar, sb.sugar);
            if (c != 0)
                return c;
            return initial.compare(sa, sb);
        };
    }

    /** Ecart comparator */
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

    /** whether monomial order is graded */
    static boolean isGradedOrder(Comparator<DegreeVector> monomialOrder) {
        return monomialOrder == MonomialOrder.GREVLEX || monomialOrder == MonomialOrder.GRLEX;
    }

    /** l.c.m. of degree vectors */
    static DegreeVector lcm(DegreeVector a, DegreeVector b) {
        return new DegreeVector(lcm(a.exponents, b.exponents));
    }

    /** l.c.m. of degree vectors */
    static int[] lcm(int[] a, int[] b) {
        return ArraysUtil.max(a, b);
    }

    /** whether dividend degree vector can be divided by divider degree vector */
    private static boolean dividesQ(int[] divider, int[] dividend) {
        for (int i = 0; i < dividend.length; i++)
            if (dividend[i] < divider[i])
                return false;
        return true;
    }

    /** whether dividend monomial can be divided by divider monomial */
    private static boolean dividesQ(DegreeVector divider, DegreeVector dividend) {
        return dividend.dvDivisibleBy(divider);
    }

    /** have variables in common */
    static boolean shareVariablesQ(DegreeVector a, DegreeVector b) {
        for (int i = 0; i < a.exponents.length; i++)
            if (a.exponents[i] != 0 && b.exponents[i] != 0)
                return true;
        return false;
    }

    /* ************************************** Buchberger algorithm ********************************************** */

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm. It uses normal strategy
     * for graded orders and sugar strategy for lexicographic.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> ideal, Comparator<DegreeVector> monomialOrder) {
        Comparator<SyzygyPair> selectionStrategy = normalSelectionStrategy(ideal.get(0).ordering);
        // fixme use sugar always?
        if (!isGradedOrder(monomialOrder))
            // add sugar for non-graded orders
            selectionStrategy = withSugar(selectionStrategy);
        return BuchbergerGB(ideal, monomialOrder, selectionStrategy, NO_MINIMIZATION);
    }

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm.
     *
     * @param generators           generators of the ideal
     * @param monomialOrder        monomial order to use
     * @param selectionStrategy    critical pair selection strategy
     * @param minimizationStrategy how to minimize Groebner basis at intermediate steps
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> generators,
                            Comparator<DegreeVector> monomialOrder,
                            Comparator<SyzygyPair> selectionStrategy,
                            MinimizationStrategy minimizationStrategy) {
        return BuchbergerGB(generators, monomialOrder, minimizationStrategy,
                () -> new SyzygyTreeSet<>(new TreeSet<>(selectionStrategy)));
    }

    /**
     * Computes minimized and reduced Groebner basis of a given homogeneous ideal via Buchberger algorithm.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerHomogeneousGB(List<Poly> generators,
                                       Comparator<DegreeVector> monomialOrder) {
        Comparator<SyzygyPair> selectionStrategy = normalSelectionStrategy(generators.get(0).ordering);
        // fixme use sugar always?
        if (!isGradedOrder(monomialOrder))
            // add sugar for non-graded orders
            selectionStrategy = withSugar(selectionStrategy);
        return BuchbergerHomogeneousGB(generators, monomialOrder, selectionStrategy, NO_MINIMIZATION);
    }

    /**
     * Computes minimized and reduced Groebner basis of a given homogeneous ideal via Buchberger algorithm.
     *
     * @param generators        generators of the ideal
     * @param monomialOrder     monomial order to use
     * @param selectionStrategy critical pair selection strategy
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerHomogeneousGB(List<Poly> generators,
                                       Comparator<DegreeVector> monomialOrder,
                                       Comparator<SyzygyPair> selectionStrategy,
                                       MinimizationStrategy minimizationStrategy) {
        return BuchbergerGB(generators, monomialOrder, minimizationStrategy,
                () -> new GradedSyzygyTreeSet<>(new TreeMap<>(), selectionStrategy));
    }

    /**
     * Generic implementation of Buchberger algorithm.
     *
     * @param generators           generators of the ideal
     * @param monomialOrder        monomial order to use
     * @param minimizationStrategy minimization strategy
     * @param syzygySetSupplier    supplies data structure for managing syzygies
     */
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerGB(List<Poly> generators,
                            Comparator<DegreeVector> monomialOrder,
                            MinimizationStrategy minimizationStrategy,
                            Supplier<SyzygySet<Term, Poly>> syzygySetSupplier) {
        Poly factory = generators.get(0);
        // simplify generators as much as possible
        generators = prepareGenerators(generators, monomialOrder);
        if (generators.size() == 1)
            return canonicalize(generators);

        // sort polynomials in the basis to achieve faster divisions
        Comparator<Poly> polyOrder = isGradedOrder(monomialOrder)
                ? (a, b) -> monomialOrder.compare(a.lt(), b.lt())
                : new EcartComparator<>();
        generators.sort(polyOrder);

        // set of syzygies
        SyzygySet<Term, Poly> sPairs = syzygySetSupplier.get();

        // Groebner basis that will be computed
        List<Poly> groebner = new ArrayList<>();
        // update Groebner basis with initial generators
        generators.forEach(g -> updateBasis(groebner, sPairs, g));

        // cache array used in divisions (little performance improvement actually)
        Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));

        // cache size of basis after each minimization
        int sizeAfterMinimization = groebner.size();
        while (!sPairs.isEmpty())
            // pick up (and remove) a bunch of critical pairs
            for (SyzygyPair<Term, Poly> pair : sPairs.getAndRemoveNextBunch()) {
                // compute actual syzygy
                Poly syzygy = MultivariateDivision.remainder(syzygy(pair), groebnerArray);
                if (syzygy.isZero())
                    continue;

                if (syzygy.isConstant())
                    // ideal = ring
                    return Collections.singletonList(factory.createOne());

                // add syzygy to basis
                updateBasis(groebner, sPairs, syzygy);
                // recompute array
                groebnerArray = groebner.stream().filter(Objects::nonNull).toArray(factory::createArray);
                // don't sort here, not practical actually
                // Arrays.sort(groebnerArray, polyOrder);

                if (minimizationStrategy.doMinimize(sizeAfterMinimization, groebner.size())) {
                    reduceAndMinimizeGroebnerBases(groebner, sPairs, sizeAfterMinimization);
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

    /** reduce & minimize current basis */
    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    void reduceAndMinimizeGroebnerBases(List<Poly> generators,
                                        SyzygySet<Term, Poly> sPairs,
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
                            if (l != j && generators.get(l) != null)
                                sPairs.remove(new SyzygyPair<>(l, j, generators));
                    } else
                        // update all pairs with k
                        for (int l = 0; l < generators.size(); l++)
                            if (l != j && generators.get(l) != null)
                                sPairs.add(new SyzygyPair<>(l, j, generators));
                }
            }
        }
    }

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly remainder0(Poly dividend, List<Poly> dividers) {
        Poly[] dividersArr = dividers.stream().filter(Objects::nonNull).toArray(dividend::createArray);
        //Arrays.sort(dividersArr, (a, b) -> a.ordering.compare(a.lt(), b.lt()));
        return MultivariateDivision.remainder(dividend, dividersArr);
    }

    /* ************************************************** F4 ******************************************************* */

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via F4 algorithm.
     *
     * @param generators    generators of the ideal
     * @param monomialOrder monomial order to use
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> F4GB(List<Poly> generators,
                    Comparator<DegreeVector> monomialOrder) {
        return F4GB(generators, monomialOrder, () -> new GradedSyzygyTreeSet<>(new TreeMap<>(), normalSelectionStrategy(monomialOrder)));
    }

    /*
     * In F4 algorithm we don't need tree-based structure for multivariate polynomials, since no any additions/divisions
     * performed (all operations are reduced to matrix operations). The only arithmetic operation used under the hood
     * of F4 is multiplication by monomial, which still preserves the order of terms. So, we switch here to array-based
     * (with sorted array) representation of polynomials which is much faster than tree-based.
     */

    /**
     * Array-based polynomial representation. The data array is sorted in descending order, so the first element is the
     * leading term.
     */
    private static final class ArrayBasedPoly<Term extends AMonomial<Term>> implements MonomialSetView<Term> {
        /** sorted (in reverse order, so lt is the first) array of terms */
        final Term[] data;
        final int[] degrees;

        ArrayBasedPoly(AMultivariatePolynomial<Term, ?> poly) {
            // retrieve data array from poly in descending order
            this.data = poly.terms.descendingMap().values().toArray(poly.monomialAlgebra.createArray(poly.size()));
            this.degrees = poly.degrees();
        }

        ArrayBasedPoly(Term[] data, int[] degrees) {
            this.data = data;
            this.degrees = degrees;
        }

        ArrayBasedPoly(Term[] data, int nVariables) {
            this.data = data;
            this.degrees = new int[nVariables];
            for (Term db : data)
                for (int i = 0; i < nVariables; i++)
                    if (db.exponents[i] > degrees[i])
                        degrees[i] = db.exponents[i];
        }

        @Override
        public Iterator<Term> ascendingIterator() { throw new UnsupportedOperationException(); }

        @Override
        public Iterator<Term> descendingIterator() { return Arrays.asList(data).iterator(); }

        @Override
        public Iterator<Term> iterator() { return descendingIterator(); }

        @Override
        public int[] degrees() { return degrees; }

        @Override
        public Term first() { return data[data.length - 1]; }

        @Override
        public Term last() { return data[0]; }

        int size() { return data.length; }

        Term get(int i) { return data[i]; }

        Term find(Term dv, Comparator<DegreeVector> ordering) {
            int i = Arrays.binarySearch(data, dv, (a, b) -> ordering.compare(b, a));
            if (i < 0)
                return null;
            return data[i];
        }

        @Override
        public Collection<Term> collection() { return Arrays.asList(data); }

        @Override
        public ArrayBasedPoly<Term> clone() { return new ArrayBasedPoly<>(data.clone(), degrees.clone()); }

        @Override
        public String toString() {
            return Arrays.toString(data);
        }
    }

    /** Multiply array-based poly by a monomial */
    static <Term extends AMonomial<Term>>
    void multiply(ArrayBasedPoly<Term> poly, Term term, IMonomialAlgebra<Term> mAlgebra) {
        Term[] data = poly.data;
        for (int i = data.length - 1; i >= 0; --i)
            data[i] = mAlgebra.multiply(data[i], term);
        int[] degrees = poly.degrees;
        for (int i = 0; i < degrees.length; ++i)
            degrees[i] += term.exponents[i];
    }

    /** Multiply array-based poly by a monomial */
    static <Term extends AMonomial<Term>>
    void monic(ArrayBasedPoly<Term> poly, IMonomialAlgebra<Term> mAlgebra) {
        Term unit = mAlgebra.getUnitTerm(poly.lt().exponents.length);
        multiply(poly, mAlgebra.divideExact(unit, poly.lt().setDegreeVector(unit)), mAlgebra);
    }

    /** Minimal size of selected set of critical pairs */
    private static final int F4_MIN_SELECTION_SIZE = 0;

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Faugère's F4 algorithm.
     *
     * Simplification routine used is due to Joux & Vitse, "A variant of F4 algorithm",
     * https://eprint.iacr.org/2010/158.pdf, 2011
     *
     * @param generators    generators of the ideal
     * @param monomialOrder monomial order to use
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> F4GB(List<Poly> generators,
                    Comparator<DegreeVector> monomialOrder,
                    Supplier<SyzygySet<Term, ArrayBasedPoly<Term>>> syzygySetSupplier) {
        Poly factory = generators.get(0);
        IMonomialAlgebra<Term> mAlgebra = factory.monomialAlgebra;

        // simplify generators as much as possible
        generators = prepareGenerators(generators, monomialOrder);
        if (generators.size() == 1)
            return generators;

        // sort polynomials in the basis to achieve in general faster inital reductions
        Comparator<Poly> polyOrder = isGradedOrder(monomialOrder)
                ? (a, b) -> monomialOrder.compare(a.lt(), b.lt())
                : new EcartComparator<>();
        generators.sort(polyOrder);

        // set of all syzygies
        SyzygySet<Term, ArrayBasedPoly<Term>> sPairs = syzygySetSupplier.get();
        // Groebner basis that will be computed
        List<ArrayBasedPoly<Term>> groebner = new ArrayList<>();

        // a history of performed reductions: index of generator -> list of available reductions
        // this is extensively used in simplify routine; the data is updated each time when a new row-echelon
        // form is calculated; for details see Joux & Vitse, "A variant of F4 algorithm" (TabSimplify structure)
        List<List<ArrayBasedPoly<Term>>> f4reductions = new ArrayList<>();

        // update Groebner basis with initial generators
        for (Poly generator : generators) {
            ArrayBasedPoly<Term> g = new ArrayBasedPoly<>(generator);
            updateBasis(groebner, sPairs, g);
            f4reductions.add(new ArrayList<>(Collections.singletonList(g)));
        }

        while (!sPairs.isEmpty()) {
            // pick up and remove a bunch of pairs (select at least F4_MIN_SELECTION_SIZE pairs)
            Collection<SyzygyPair<Term, ArrayBasedPoly<Term>>> subset = sPairs.getAndRemoveNextBunch();
            while (!sPairs.isEmpty() && subset.size() < F4_MIN_SELECTION_SIZE)
                subset.addAll(sPairs.getAndRemoveNextBunch());


            // the list of polynomials to reduce (H-polynomials)
            List<HPolynomial<Term>> hPolynomials = new ArrayList<>();
            // all monomials that occur in H-polynomials
            TreeSet<DegreeVector> hMonomials = new TreeSet<>(monomialOrder);
            // monomials that were annihilated
            TreeSet<DegreeVector> hAnnihilated = new TreeSet<>(monomialOrder);

            // form initial set of H-polynomials from the current set of critical pairs
            for (SyzygyPair<Term, ArrayBasedPoly<Term>> syzygy : subset) {
                // correction factors
                Term
                        fiQuot = mAlgebra.divideExact(syzygy.syzygyGamma, syzygy.fi.lt()),
                        fjQuot = mAlgebra.divideExact(syzygy.syzygyGamma, syzygy.fj.lt());

                // a pair of H-polynomials
                ArrayBasedPoly<Term>
                        fiHPoly = simplify(syzygy.fi, fiQuot, f4reductions.get(syzygy.i), mAlgebra),
                        fjHPoly = simplify(syzygy.fj, fjQuot, f4reductions.get(syzygy.j), mAlgebra);

                assert fiHPoly.lt().equals(fjHPoly.lt()) : fiHPoly.lt() + " " + fjHPoly.lt();

                // add to set of H-polynomials
                hPolynomials.add(new HPolynomial<>(fiHPoly, syzygy.i));
                hPolynomials.add(new HPolynomial<>(fjHPoly, syzygy.j));

                // store all monomials that we have in H
                hMonomials.addAll(fiHPoly.collection());
                hMonomials.addAll(fjHPoly.collection());

                // lts will be annihilated
                hAnnihilated.add(fiHPoly.lt());
            }

            // the diff = Mon(H-polynomials) / annihilated
            TreeSet<DegreeVector> diff = new TreeSet<>(monomialOrder);
            diff.addAll(hMonomials);
            diff.removeAll(hAnnihilated);

            // compute a whole set of required H-polynomials for reductions
            // ("Symbolic Preprocessing" routine from Faugère's original paper)
            while (!diff.isEmpty()) {
                // pick the "highest" term from diff
                DegreeVector dv = diff.pollLast();
                hAnnihilated.add(dv);

                // select some polynomial from basis which can reduce the "highest" term (lt-divisor)
                // todo selection possible
                OptionalInt divisorOpt = IntStream.range(0, groebner.size())
                        .filter(i -> groebner.get(i) != null)
                        .filter(i -> dv.dvDivisibleBy(groebner.get(i).lt()))
                        .findAny();
                if (divisorOpt.isPresent()) {
                    // index in the basis
                    int iIndex = divisorOpt.getAsInt();
                    // correction factor
                    DegreeVector quot = dv.dvDivideExact(groebner.get(iIndex).lt());
                    // new H-polynomial to add
                    ArrayBasedPoly<Term> newH = simplify(groebner.get(iIndex), mAlgebra.create(quot), f4reductions.get(iIndex), mAlgebra);

                    // append newH to H-polynomials
                    hPolynomials.add(new HPolynomial<>(newH, iIndex));
                    // update monomials set
                    TreeSet<DegreeVector> newMonomials = new TreeSet<>(monomialOrder);
                    newMonomials.addAll(newH.collection());
                    hMonomials.addAll(newMonomials);

                    // update the diff
                    newMonomials.removeAll(hAnnihilated);
                    diff.addAll(newMonomials);
                }
            }

            // all monomials occurring in H sorted in descending order, i.e. lead term first
            DegreeVector[] hMonomialsArray = hMonomials.descendingSet().toArray(new DegreeVector[hMonomials.size()]);
            // sort rows to make the initial matrix upper-right-triangular as max as possible
            hPolynomials.sort((a, b) -> {
                int c = monomialOrder.compare(b.hPoly.lt(), a.hPoly.lt());
                if (c != 0)
                    return c;
                return Integer.compare(a.hPoly.size(), b.hPoly.size());
            });
            // reduce all H-polynomials with linear algebra and compute new basis elements (N+)
            List<ArrayBasedPoly<Term>> nPlus = nPlus(factory, groebner, hPolynomials, hMonomialsArray, f4reductions);
            // enlarge the basis
            nPlus.forEach(g -> updateBasis(groebner, sPairs, g));
        }
        // batch remove all nulls
        groebner.removeAll(Collections.<ArrayBasedPoly<Term>>singleton(null));

        // convert from array-based polynomials to normal
        List<Poly> result = groebner.stream().map(factory::create).collect(Collectors.toList());
        // minimize Groebner basis
        minimizeGroebnerBases(result);
        // speed up final reduction
        result.sort(polyOrder);
        // reduce Groebner basis
        removeRedundant(result);
        // canonicalize Groebner basis
        canonicalize(result);

        return result;
    }

    /**
     * A simple container for polynomial and the index of this polynomial (or initial basis polynomial which was reduced
     * in F4 to this polynomial) in the list of generators
     */
    static final class HPolynomial<Term extends AMonomial<Term>> {
        /** polynomial */
        final ArrayBasedPoly<Term> hPoly;
        /** index in generators list of the original polynomial */
        final int indexOfOrigin;

        HPolynomial(ArrayBasedPoly<Term> hPoly, int indexOfOrigin) {
            this.hPoly = hPoly;
            this.indexOfOrigin = indexOfOrigin;
        }
    }

    /**
     * Simplifies h-polynomial by replacing it with its previously computed reduction. For details see Joux & Vitse, "A
     * variant of F4 algorithm", https://eprint.iacr.org/2010/158.pdf, 2011
     *
     * @param generator  the generator (H-polynomial to add is a product of factor * generator)
     * @param factor     correction factor for generator (H-polynomial to add is a product of factor * generator)
     * @param reductions a list of previously obtained reductions of the generator
     * @param mAlgebra   monomial algebra
     */
    static <Term extends AMonomial<Term>>
    ArrayBasedPoly<Term> simplify(ArrayBasedPoly<Term> generator,
                                  Term factor,
                                  List<ArrayBasedPoly<Term>> reductions,
                                  IMonomialAlgebra<Term> mAlgebra) {
        // the desired leading term of H-polynomial
        Term desiredLeadTerm = mAlgebra.multiply(generator.lt(), factor);

        // iterate from the last (most recently added, thus most simplified) to the first (the initial generator) reductions
        for (int i = reductions.size() - 1; i >= 0; --i) {
            ArrayBasedPoly<Term> reduction = reductions.get(i);
            // leading term of previously computed reduction
            Term rLeadTerm = reduction.lt();

            if (rLeadTerm.dvEquals(desiredLeadTerm))
                // we just have the same calculation in the map (no any appropriate reductions were calculated)
                return setLeadCoeffFrom(reduction, desiredLeadTerm, mAlgebra);

            if (desiredLeadTerm.dvDivisibleBy(rLeadTerm)) {
                // <- nontrivial appropriate reduction

                DegreeVector quot = desiredLeadTerm.dvDivideExact(rLeadTerm);
                ArrayBasedPoly<Term> g = reduction.clone();
                multiply(g, mAlgebra.create(quot), mAlgebra);
                g = setLeadCoeffFrom(g, desiredLeadTerm, mAlgebra);

                // cache this reduction too
                reductions.add(g);
                return g;
            }
        }

        // <- this point is unreachable since the initial generator is already added to the list (in the begining)
        throw new RuntimeException();
    }

//    /**
//     * Simplifies h-polynomial by replacing it with its previously computed reduction. For details see Joux & Vitse, "A
//     * variant of F4 algorithm", https://eprint.iacr.org/2010/158.pdf, 2011
//     *
//     * @param generator  the generator (H-polynomial to add is a product of factor * generator)
//     * @param factor     correction factor for generator (H-polynomial to add is a product of factor * generator)
//     * @param reductions a list of previously obtained reductions of the generator
//     * @param mAlgebra   monomial algebra
//     */
//    static <Term extends AMonomial<Term>>
//    ArrayBasedPoly<Term> simplify(ArrayBasedPoly<Term> generator,
//                                  Term factor,
//                                  List<ArrayBasedPoly<Term>> reductions,
//                                  IMonomialAlgebra<Term> mAlgebra) {
//        // the desired leading term of H-polynomial
//        Term desiredLeadTerm = mAlgebra.multiply(generator.lt(), factor);
//
//        ArrayBasedPoly<Term> bestReduction = null;
//        // iterate from the last (most recently added, thus most simplified) to the first (the initial generator) reductions
//        for (int i = reductions.size() - 1; i >= 0; --i) {
//            ArrayBasedPoly<Term> reduction = reductions.get(i);
//            // leading term of previously computed reduction
//            Term rLeadTerm = reduction.lt();
//
//            if (rLeadTerm.dvEquals(desiredLeadTerm) || desiredLeadTerm.dvDivisibleBy(rLeadTerm)) {
//                if (bestReduction == null
//                        || reduction.size() < bestReduction.size()
////                        || MonomialOrder.GREVLEX.compare(reduction.lt(), bestReduction.lt()) > 0
//                        )
//                    bestReduction = reduction;
//            }
//        }
//
//        ArrayBasedPoly<Term> reduction = bestReduction;
//        // leading term of previously computed reduction
//        Term rLeadTerm = reduction.lt();
//        if (rLeadTerm.dvEquals(desiredLeadTerm)) {
//            // we just have the same calculation in the map (no any appropriate reductions were calculated)
//            return setLeadCoeffFrom(reduction, desiredLeadTerm, mAlgebra);
//        }
//
//        if (desiredLeadTerm.dvDivisibleBy(rLeadTerm)) {
//            // <- nontrivial appropriate reduction
//            DegreeVector quot = desiredLeadTerm.dvDivideExact(rLeadTerm);
//            ArrayBasedPoly<Term> g = reduction.clone();
//            multiply(g, mAlgebra.create(quot), mAlgebra);
//            g = setLeadCoeffFrom(g, desiredLeadTerm, mAlgebra);
//
//            // cache this reduction too
//            reductions.add(g);
//            return g;
//        }
//
//        // <- this point is unreachable since the initial generator is already added to the list (in the begining)
//        throw new RuntimeException();
//    }


    /** Sets the leading coefficient of poly from the coefficient of term (poly is copied only if necessary) */
    static <Term extends AMonomial<Term>>
    ArrayBasedPoly<Term> setLeadCoeffFrom(ArrayBasedPoly<Term> poly, Term term, IMonomialAlgebra<Term> mAlgebra) {
        if (mAlgebra.haveSameCoefficients(poly.lt(), term))
            return poly;
        poly = poly.clone();
        multiply(poly, mAlgebra.divideExact(term.toZero(), poly.lt().toZero()), mAlgebra);
        return poly;
    }

    /**
     * Computes row echelon form of the matrix and returns new elements to add to the basis
     *
     * @param factory      factory polynomials
     * @param basis        the basis computed so far
     * @param hPolynomials list of prepared H-polynomials
     * @param hMonomials   all the monomials occurring in the H-polynomials list (sorted in descending order)
     * @param f4reductions reductions computed previously for all basis elements
     */
    @SuppressWarnings("unchecked")
    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<ArrayBasedPoly<Term>> nPlus(Poly factory,
                                     List<ArrayBasedPoly<Term>> basis,
                                     List<HPolynomial<Term>> hPolynomials,
                                     DegreeVector[] hMonomials,
                                     List<List<ArrayBasedPoly<Term>>> f4reductions) {
        if (hPolynomials.isEmpty())
            return new ArrayList<>();

        if (factory instanceof MultivariatePolynomialZp64)
            return (List<ArrayBasedPoly<Term>>) nPlusZp64_1(
                    (MultivariatePolynomialZp64) factory,
                    (List) basis,
                    (List) hPolynomials,
                    hMonomials,
                    (List) f4reductions);
        else
            throw new RuntimeException();
    }

    private static final double DENSE_FILLING_THRESHOLD = 0.1;

    static void axsy(IntegersZp64 ring, long[] denseArray, long factor, SparseArray arr) {
        for (int i = 0; i < arr.densePositions.length; i++)
            denseArray[arr.densePositions[i]] = ring.subtract(denseArray[arr.densePositions[i]], ring.multiply(factor, arr.denseValues[i]));
        for (int i = 0; i < arr.sparsePositions.length; i++)
            denseArray[arr.sparsePositions[i]] = ring.subtract(denseArray[arr.sparsePositions[i]], ring.multiply(factor, arr.sparseValues[i]));
    }

    static long STEP0= 0, STEP1= 0, STEP2 = 0, STEP3 = 0, STEP4= 0, STEP5 = 0, STEP6 = 0;
    // for Zp64
    @SuppressWarnings("unchecked")
    private static List<ArrayBasedPoly<MonomialZp64>> nPlusZp64_1(
            MultivariatePolynomialZp64 factory,
            List<ArrayBasedPoly<MonomialZp64>> basis,
            List<HPolynomial<MonomialZp64>> hPolynomials,
            DegreeVector[] hMonomials,
            List<List<ArrayBasedPoly<MonomialZp64>>> f4reductions) {
        IntegersZp64 ring = factory.ring;

        // We use the structured Gaussian elimination strategy as described in
        // J.-C. Faugere & S. Lachartre, PASCO'10 https://doi.org/10.1145/1837210.1837225
        // "Parallel Gaussian Elimination for Gröbner bases computations in finite fields"
        //
        // By swapping the rows and non-pivoting columns, each matrix in F4
        // can be rewritten in the following form (F4-form):
        //
        //    \x0x00xx000x | x0x00xx0000x0x0x0xx000xx000xx0
        //    0\x0x0x00xx0 | 0xx0000x0xxx0x00xx0x0000xx0000
        //    00\x00x000x0 | 0x0x0000x0x0xx00xx0x00x0x0xx00
        //    000\xx0x0x00 | xx0xxx00x000x0x0xx00x0x0xx000x
        //    0000\xx0x0x0 | x0000xx0x00x0xxx0xx0000x000xx0
        //    00000\x0000x | 00x0000x0x0x0xx0xx0xx000xx0000
        //    000000\xx00x | 0x0x000x00x0xxx0xx00xxx0x0xx00
        //    0000000\x0x0 | xx00xx00xx00x000xx0xx00x0x000x
        //    ............ | ..............................
        //    -------------+-------------------------------
        //    0xx000x0x0xx | xxxxxx0xxxxxxx0xxxxxxxxxxxxxxx
        //    x0xx000x0x00 | xxxx0xxxxxxxxxxxxxx0xxxxxxxxxx
        //    00x00x0000xx | xxxxxxx0xxxxxxxxxxxxxxx0xxxxxx
        //    x0000x00xx0x | xxxxxxxxxxxxxxxxx0xxxxxxx0xxxx
        //    ............ | ..............................
        //
        // We denote:
        //
        // A - upper left  block (very sparse, triangular)         -- pivoting rows
        // B - upper right block (partially sparse, rectangular)   -- pivoting rows
        // C -  down left  block (partially  sparse, rectangular)  -- non-pivoting rows
        // D -  down right block (dense, rectangular)              -- non-pivoting rows
        //
        // The algorithm to reduce the matrix is then very simple:
        //
        // 1) row reduce A (B is still partially sparse)
        // 2) annihilate C (D is now almost certainly dense)
        // 3) row echelon & row reduce D
        // 4) row reduce B

        // reverse order for binary searching
        Comparator<DegreeVector> reverseOrder = (a, b) -> factory.ordering.compare(b, a);
        int
                nRows = hPolynomials.size(),
                nColumns = hMonomials.length,
                iRow, iColumn;

        long start;
        /* ======  STEP 0: bring matrix to F4-form  ====== */

        start = System.nanoTime();
        // number of lead-terms in each column (columns with no any lead
        // terms are non-pivoting and can be rearranged)
        int[] columnsLeadTermsFilling = new int[nColumns];
        // detect non-pivoting columns
        int iOldColumnPrev = 0;
        for (HPolynomial<MonomialZp64> hPoly : hPolynomials) {
            iColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, hPoly.hPoly.lt(), reverseOrder);
            iOldColumnPrev = iColumn;
            assert iColumn >= 0;
            ++columnsLeadTermsFilling[iColumn];
        }

        // find non pivoting columns
        TIntArrayList nonPivotColumns = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows + nonPivotColumns.size() && iColumn < nColumns; ++iColumn)
            if (columnsLeadTermsFilling[iColumn] == 0)
                nonPivotColumns.add(iColumn);
        // now we move non-pivoting columns to the right
        // mapping between old and new columns numeration
        int[] nonPivotColumnsArr = nonPivotColumns.toArray();
        int[] columnsRearrangement = ArraysUtil.addAll(
                ArraysUtil.intSetDifference(ArraysUtil.sequence(0, nColumns), nonPivotColumnsArr),
                nonPivotColumnsArr);

        // back mapping between new and old columns numeration
        int[] columnsBackRearrangement = new int[nColumns];
        for (int i = 0; i < nColumns; ++i)
            columnsBackRearrangement[columnsRearrangement[i]] = i;

        // number of non-zero entries in each column
        int[] columnsFilling = new int[nColumns];
        // index of each term in hPolynomials in the hMonomials array
        int[][] mapping = new int[nRows][];
        // estimated row filling of B matrix (number of nonzero elements in each row)
        int[] bRowsFilling = new int[nRows];
        // first iteration: gather info about matrix pattern
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;
            mapping[iRow] = new int[hPoly.size()];
            iOldColumnPrev = 0;
            for (int i = 0; i < hPoly.size(); ++i) {
                MonomialZp64 term = hPoly.get(i);
                // column in old numeration
                int iOldColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, term, reverseOrder);
                iOldColumnPrev = iOldColumn;
                assert iOldColumn >= 0;
                // column in new numeration
                iColumn = columnsBackRearrangement[iOldColumn];
                mapping[iRow][i] = iColumn;

                ++columnsFilling[iColumn];
                if (iColumn >= nRows)
                    ++bRowsFilling[iRow];
            }
        }

        // choose pivoting rows, so that B matrix is maximally sparse (rows with minimal bFillIns)
        TIntArrayList pivots = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows; ++iColumn) {
            int minFillIn = Integer.MAX_VALUE;
            int pivot = -1;
            for (iRow = iColumn; iRow < nRows; ++iRow)
                if (mapping[iRow][0] == iColumn && bRowsFilling[iRow] < minFillIn) {
                    minFillIn = bRowsFilling[iRow];
                    pivot = iRow;
                } else if (pivot != -1 && mapping[iRow][0] != iColumn)
                    break;
            if (pivot == -1)
                break;
            pivots.add(pivot);
        }
        bRowsFilling = null; // prevent further use

        // rearrange rows: move pivots up and non-pivots down
        int nPivotRows = pivots.size();
        for (int i = 0; i < nPivotRows; ++i) {
            int pivot = pivots.get(i);
            Collections.swap(hPolynomials, i, pivot);
            ArraysUtil.swap(mapping, i, pivot);
        }

        for (int i = 0; i < nPivotRows; ++i)
            monic(hPolynomials.get(i).hPoly, factory.monomialAlgebra);

        // <- the matrix is now is in the desired F4-form

        STEP0 += System.nanoTime() - start;
        /* ======  STEP 1: prepare data structures  ====== */

        start = System.nanoTime();
        // dense columns in matrices B & D
        int[] bDenseColumns = IntStream.range(nPivotRows, nColumns)
                .filter(i -> 1.0 * columnsFilling[i] / nRows > DENSE_FILLING_THRESHOLD)
                .map(i -> i - nPivotRows)
                .toArray(); // <- it is sorted (for below binary searching)

        // A & C are very sparse
        SparseRowMatrix
                aMatrix = new SparseRowMatrix(ring, nPivotRows, nPivotRows, new int[0]),
                cMatrix = new SparseRowMatrix(ring, nRows - nPivotRows, nPivotRows, new int[0]);

        // sparse matrices B & D
        SparseRowMatrix
                bMatrix = new SparseRowMatrix(ring, nPivotRows, nColumns - nPivotRows, bDenseColumns),
                dMatrix = new SparseRowMatrix(ring, nRows - nPivotRows, nColumns - nPivotRows, bDenseColumns);

        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;

            TIntArrayList
                    acSparseCols = new TIntArrayList(),
                    bdSparseCols = new TIntArrayList();
            TLongArrayList
                    acSparseVals = new TLongArrayList(),
                    bdSparseVals = new TLongArrayList();
            long[] bdDenseVals = new long[bDenseColumns.length];
            for (int i = 0; i < hPoly.size(); i++) {
                iColumn = mapping[iRow][i];
                long coefficient = hPoly.get(i).coefficient;
                if (iColumn < nPivotRows) {
                    // element of matrix A or C
                    acSparseCols.add(iColumn);
                    acSparseVals.add(coefficient);
                } else {
                    // element of matrix B or D
                    iColumn -= nPivotRows;
                    int iDense;
                    if ((iDense = Arrays.binarySearch(bDenseColumns, iColumn)) >= 0)
                        // this is dense column (in B or D)
                        bdDenseVals[iDense] = coefficient;
                    else {
                        bdSparseCols.add(iColumn);
                        bdSparseVals.add(coefficient);
                    }
                }
            }

            int[] bdSparseColumns = bdSparseCols.toArray();
            long[] bdSparseValues = bdSparseVals.toArray();
            ArraysUtil.quickSort(bdSparseColumns, bdSparseValues);
            int[] acSparseColumns = acSparseCols.toArray();
            long[] acSparseValues = acSparseVals.toArray();
            ArraysUtil.quickSort(acSparseColumns, acSparseValues);

            if (iRow < nPivotRows) {
                aMatrix.rows[iRow] = new SparseArray(ring, nPivotRows, new int[0], acSparseColumns, new long[0], acSparseValues);
                bMatrix.rows[iRow] = new SparseArray(ring, nColumns - nPivotRows, bDenseColumns, bdSparseColumns, bdDenseVals, bdSparseValues);
            } else {
                cMatrix.rows[iRow - nPivotRows] = new SparseArray(ring, nPivotRows, new int[0], acSparseColumns, new long[0], acSparseValues);
                dMatrix.rows[iRow - nPivotRows] = new SparseArray(ring, nColumns - nPivotRows, bDenseColumns, bdSparseColumns, bdDenseVals, bdSparseValues);
            }
        }

        STEP1 += System.nanoTime() - start;
        /* ======  STEP 2: row reduce matrix A  ====== */

        start = System.nanoTime();
        // we start from the last column in matrix A
        for (iRow = nPivotRows - 2; iRow >= 0; --iRow) {
            SparseArray aRow = aMatrix.rows[iRow];
            long[] bRow = bMatrix.rows[iRow].toDense();

            for (int i = 0; i < aRow.sparsePositions.length; ++i) {
                iColumn = aRow.sparsePositions[i];
                if (iColumn == iRow) // diagonal
                    continue;
                axsy(ring, bRow, aRow.sparseValues[i], bMatrix.rows[iColumn]);
            }
            bMatrix.rows[iRow] = new SparseArray(ring, bMatrix.densePositions, bRow);
        }

        STEP2 += System.nanoTime() - start;

        /* ======  STEP 3: annihilate matrix C  ====== */

        start = System.nanoTime();
        for (iRow = 0; iRow < nRows - nPivotRows; ++iRow) {
            SparseArray cRow = cMatrix.rows[iRow];
            long[] dRow = dMatrix.rows[iRow].toDense();
            for (int i = 0; i < cRow.sparsePositions.length; ++i) {
                iColumn = cRow.sparsePositions[i];
                axsy(ring, dRow, cRow.sparseValues[i], bMatrix.rows[iColumn]);
            }
            dMatrix.rows[iRow] = new SparseArray(ring, dMatrix.densePositions, dRow);
        }
        STEP3 += System.nanoTime() - start;
//
//        for (iColumn = 0; iColumn < nPivotRows; ++iColumn) {
//            List<ColumnElement> column = cColumns[iColumn];
//            // annihilate each element in the column
//            for (int i = column.size() - 1; i >= 0; --i) {
//                ColumnElement cElement = column.get(i);
//                dMatrix.rows[cElement.iRow].subtract(bMatrix.rows[iColumn], cElement.value);
//            }
//        }

        /* ======  STEP 4: compute row reduced echelon form of matrix D  ====== */
        start = System.nanoTime();

        // this can be optimized (use dense structures, use same structured elimination), but actually
        // it doesn't take too much time

        // make D maximally triangular
        Arrays.sort(dMatrix.rows, Comparator.comparingInt(SparseArray::firstNonZeroPosition));
        dMatrix.rowReduce();

        STEP4 += System.nanoTime() - start;

        /* ======  STEP 5: row reduce B  ====== */

        start = System.nanoTime();
        int dShift = 0;
        for (iRow = 0; iRow < dMatrix.nRows; iRow++) {
            SparseArray dRow = dMatrix.rows[iRow];
            iColumn = iRow + dShift;
            if (iColumn >= nColumns)
                break;
            if (dRow.coefficient(iColumn) == 0) {
                --iRow;
                ++dShift;
                continue;
            }
            for (int i = 0; i < bMatrix.nRows; ++i)
                if (bMatrix.rows[i].coefficient(iColumn) != 0)
                    bMatrix.rows[i].subtract(dRow, bMatrix.rows[i].coefficient(iColumn));
        }


        STEP5 += System.nanoTime() - start;
        /* ======  STEP 6: finally form N+ polynomials ====== */

        start = System.nanoTime();

        // leading monomials of H-polynomials
        TreeSet<DegreeVector> hLeadMonomials = hPolynomials.stream().map(p -> p.hPoly.lt())
                .collect(Collectors.toCollection(() -> new TreeSet<>(factory.ordering)));

        List<ArrayBasedPoly<MonomialZp64>> nPolynomials = new ArrayList<>();
        // collect from A and B
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayList<MonomialZp64> candidateList = new ArrayList<>();

            if (iRow < nPivotRows)
                candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[iRow]], 1L));

            SparseArray row = iRow < nPivotRows ? bMatrix.rows[iRow] : dMatrix.rows[iRow - nPivotRows];
            for (int i = 0; i < row.sparsePositions.length; i++) {
                long val = row.sparseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + row.sparsePositions[i]]], val));
            }

            for (int i = 0; i < bDenseColumns.length; i++) {
                long val = row.denseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + bDenseColumns[i]]], val));
            }

            candidateList.sort(reverseOrder);
            if (!candidateList.isEmpty()) {
                ArrayBasedPoly<MonomialZp64> poly = new ArrayBasedPoly<>(
                        candidateList.toArray(factory.monomialAlgebra.createArray(candidateList.size())),
                        factory.nVariables);
                if (poly.lt().coefficient != 1) {
                    long factor = ring.reciprocal(poly.lt().coefficient);
                    for (int i = 0; i < poly.data.length; i++)
                        poly.data[i] = poly.data[i].setCoefficient(ring.multiply(factor, poly.data[i].coefficient));
                }
                assert poly.lt().coefficient == 1;
                nPolynomials.add(poly);
            }
        }
        nPolynomials.sort((a, b) -> reverseOrder.compare(a.lt(), b.lt()));
        // resulting N+ set
        List<ArrayBasedPoly<MonomialZp64>> nPlusPolynomials = new ArrayList<>();
        for (ArrayBasedPoly<MonomialZp64> candidate : nPolynomials) {
            if (!hLeadMonomials.contains(candidate.lt())) {
                // lt is new -> just add
                nPlusPolynomials.add(candidate);
                f4reductions.add(new ArrayList<>(Collections.singletonList(candidate)));
            } else
                // update f4reductions
                for (int iIndex = 0; iIndex < basis.size(); ++iIndex) {
                    ArrayBasedPoly<MonomialZp64> g = basis.get(iIndex);
                    if (g == null)
                        continue;
                    if (candidate.lt().dvDivisibleBy(g.lt())) {
                        List<ArrayBasedPoly<MonomialZp64>> reductions = f4reductions.get(iIndex);
                        boolean reduced = false;
                        for (int i = 0; i < reductions.size(); ++i) {
                            ArrayBasedPoly<MonomialZp64> red = reductions.get(i);
                            if (red.lt().dvEquals(candidate.lt())) {
                                reductions.set(i, candidate);
                                reduced = true;
                                break;
                            }
                        }
                        if (!reduced)
                            reductions.add(candidate);
                    }
                }
        }
        STEP6 += System.nanoTime() - start;

        return nPlusPolynomials;
    }

    // for Zp64
    @SuppressWarnings("unchecked")
    private static List<ArrayBasedPoly<MonomialZp64>> nPlusZp64_0(
            MultivariatePolynomialZp64 factory,
            List<ArrayBasedPoly<MonomialZp64>> basis,
            List<HPolynomial<MonomialZp64>> hPolynomials,
            DegreeVector[] hMonomials,
            List<List<ArrayBasedPoly<MonomialZp64>>> f4reductions) {
        IntegersZp64 ring = factory.ring;

        // We use the structured Gaussian elimination strategy as described in
        // J.-C. Faugere & S. Lachartre, PASCO'10 https://doi.org/10.1145/1837210.1837225
        // "Parallel Gaussian Elimination for Gröbner bases computations in finite fields"
        //
        // By swapping the rows and non-pivoting columns, each matrix in F4
        // can be rewritten in the following form (F4-form):
        //
        //    \x0x00xx000x | x0x00xx0000x0x0x0xx000xx000xx0
        //    0\x0x0x00xx0 | 0xx0000x0xxx0x00xx0x0000xx0000
        //    00\x00x000x0 | 0x0x0000x0x0xx00xx0x00x0x0xx00
        //    000\xx0x0x00 | xx0xxx00x000x0x0xx00x0x0xx000x
        //    0000\xx0x0x0 | x0000xx0x00x0xxx0xx0000x000xx0
        //    00000\x0000x | 00x0000x0x0x0xx0xx0xx000xx0000
        //    000000\xx00x | 0x0x000x00x0xxx0xx00xxx0x0xx00
        //    0000000\x0x0 | xx00xx00xx00x000xx0xx00x0x000x
        //    ............ | ..............................
        //    -------------+-------------------------------
        //    0xx000x0x0xx | xxxxxx0xxxxxxx0xxxxxxxxxxxxxxx
        //    x0xx000x0x00 | xxxx0xxxxxxxxxxxxxx0xxxxxxxxxx
        //    00x00x0000xx | xxxxxxx0xxxxxxxxxxxxxxx0xxxxxx
        //    x0000x00xx0x | xxxxxxxxxxxxxxxxx0xxxxxxx0xxxx
        //    ............ | ..............................
        //
        // We denote:
        //
        // A - upper left  block (very sparse, triangular)         -- pivoting rows
        // B - upper right block (partially sparse, rectangular)   -- pivoting rows
        // C -  down left  block (partially  sparse, rectangular)  -- non-pivoting rows
        // D -  down right block (dense, rectangular)              -- non-pivoting rows
        //
        // The algorithm to reduce the matrix is then very simple:
        //
        // 1) row reduce A (B is still partially sparse)
        // 2) annihilate C (D is now almost certainly dense)
        // 3) row echelon & row reduce D
        // 4) row reduce B

        // reverse order for binary searching
        Comparator<DegreeVector> reverseOrder = (a, b) -> factory.ordering.compare(b, a);
        int
                nRows = hPolynomials.size(),
                nColumns = hMonomials.length,
                iRow, iColumn;

        /* ======  STEP 0: bring matrix to F4-form  ====== */

        // number of lead-terms in each column (columns with no any lead
        // terms are non-pivoting and can be rearranged)
        int[] columnsLeadTermsFilling = new int[nColumns];
        // detect non-pivoting columns
        int iOldColumnPrev = 0;
        for (HPolynomial<MonomialZp64> hPoly : hPolynomials) {
            iColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, hPoly.hPoly.lt(), reverseOrder);
            iOldColumnPrev = iColumn;
            assert iColumn >= 0;
            ++columnsLeadTermsFilling[iColumn];
        }

        // find non pivoting columns
        TIntArrayList nonPivotColumns = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows + nonPivotColumns.size() && iColumn < nColumns; ++iColumn)
            if (columnsLeadTermsFilling[iColumn] == 0)
                nonPivotColumns.add(iColumn);
        // now we move non-pivoting columns to the right
        // mapping between old and new columns numeration
        int[] nonPivotColumnsArr = nonPivotColumns.toArray();
        int[] columnsRearrangement = ArraysUtil.addAll(
                ArraysUtil.intSetDifference(ArraysUtil.sequence(0, nColumns), nonPivotColumnsArr),
                nonPivotColumnsArr);

        // back mapping between new and old columns numeration
        int[] columnsBackRearrangement = new int[nColumns];
        for (int i = 0; i < nColumns; ++i)
            columnsBackRearrangement[columnsRearrangement[i]] = i;

        // number of non-zero entries in each column
        int[] columnsFilling = new int[nColumns];
        // index of each term in hPolynomials in the hMonomials array
        int[][] mapping = new int[nRows][];
        // estimated row filling of B matrix (number of nonzero elements in each row)
        int[] bRowsFilling = new int[nRows];
        // first iteration: gather info about matrix pattern
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;
            mapping[iRow] = new int[hPoly.size()];
            iOldColumnPrev = 0;
            for (int i = 0; i < hPoly.size(); ++i) {
                MonomialZp64 term = hPoly.get(i);
                // column in old numeration
                int iOldColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, term, reverseOrder);
                iOldColumnPrev = iOldColumn;
                assert iOldColumn >= 0;
                // column in new numeration
                iColumn = columnsBackRearrangement[iOldColumn];
                mapping[iRow][i] = iColumn;

                ++columnsFilling[iColumn];
                if (iColumn >= nRows)
                    ++bRowsFilling[iRow];
            }
        }

        // choose pivoting rows, so that B matrix is maximally sparse (rows with minimal bFillIns)
        TIntArrayList pivots = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows; ++iColumn) {
            int minFillIn = Integer.MAX_VALUE;
            int pivot = -1;
            for (iRow = iColumn; iRow < nRows; ++iRow)
                if (mapping[iRow][0] == iColumn && bRowsFilling[iRow] < minFillIn) {
                    minFillIn = bRowsFilling[iRow];
                    pivot = iRow;
                } else if (pivot != -1 && mapping[iRow][0] != iColumn)
                    break;
            if (pivot == -1)
                break;
            pivots.add(pivot);
        }
        bRowsFilling = null; // prevent further use

        // rearrange rows: move pivots up and non-pivots down
        int nPivotRows = pivots.size();
        for (int i = 0; i < nPivotRows; ++i) {
            int pivot = pivots.get(i);
            Collections.swap(hPolynomials, i, pivot);
            ArraysUtil.swap(mapping, i, pivot);
        }
        // <- the matrix is now is in the desired F4-form

        /* ======  STEP 1: prepare data structures  ====== */

        // dense columns in matrices B & D
        int[] bDenseColumns = IntStream.range(nPivotRows, nColumns)
                .filter(i -> 1.0 * columnsFilling[i] / nRows > DENSE_FILLING_THRESHOLD)
                .map(i -> i - nPivotRows)
                .toArray(); // <- it is sorted (for below binary searching)

        // we represent A & C as a list of columns
        // column lists for A & C
        List<ColumnElement>[]
                aColumns = new List[nPivotRows],
                cColumns = new List[nPivotRows];
        for (int i = 0; i < nPivotRows; ++i) {
            aColumns[i] = new ArrayList<>();
            cColumns[i] = new ArrayList<>();
        }

        // sparse matrices B & D
        SparseRowMatrix
                bMatrix = new SparseRowMatrix(ring, nPivotRows, nColumns - nPivotRows, bDenseColumns),
                dMatrix = new SparseRowMatrix(ring, nRows - nPivotRows, nColumns - nPivotRows, bDenseColumns);

        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;

            TIntArrayList sparseCols = new TIntArrayList();
            TLongArrayList sparseVals = new TLongArrayList();
            long[] denseVals = new long[bDenseColumns.length];
            for (int i = 0; i < hPoly.size(); i++) {
                iColumn = mapping[iRow][i];
                long coefficient = hPoly.get(i).coefficient;
                if (iColumn < nPivotRows) {
                    // element of matrix A or C
                    if (iRow < nPivotRows)
                        aColumns[iColumn].add(new ColumnElement(iRow, coefficient));
                    else
                        cColumns[iColumn].add(new ColumnElement(iRow - nPivotRows, coefficient));
                } else {
                    // element of matrix B or D
                    iColumn -= nPivotRows;
                    int iDense;
                    if ((iDense = Arrays.binarySearch(bDenseColumns, iColumn)) >= 0)
                        // this is dense column (in B or D)
                        denseVals[iDense] = coefficient;
                    else {
                        sparseCols.add(iColumn);
                        sparseVals.add(coefficient);
                    }
                }
            }

            int[] sparseColumns = sparseCols.toArray();
            long[] sparseValues = sparseVals.toArray();
            ArraysUtil.quickSort(sparseColumns, sparseValues);
            if (iRow < nPivotRows)
                bMatrix.rows[iRow] = new SparseArray(ring, nColumns, bDenseColumns, sparseColumns, denseVals, sparseValues);
            else
                dMatrix.rows[iRow - nPivotRows] = new SparseArray(ring, nColumns, bDenseColumns, sparseColumns, denseVals, sparseValues);
        }

        /* ======  STEP 2: row reduce matrix A  ====== */

        // we start from the last column in matrix A
        for (iColumn = nPivotRows - 1; iColumn >= 0; --iColumn) {
            List<ColumnElement> column = aColumns[iColumn];
            // monicize diagonal
            ColumnElement diagonal = column.get(column.size() - 1);
            assert diagonal.iRow == iColumn;
            bMatrix.rows[diagonal.iRow].multiply(ring.reciprocal(diagonal.value));

            // annihilate each non-diagonal element
            for (int i = column.size() - 2; i >= 0; --i) {
                ColumnElement nonDiagonal = column.get(i);
                bMatrix.rows[nonDiagonal.iRow].subtract(bMatrix.rows[diagonal.iRow], nonDiagonal.value);
            }
        }

        /* ======  STEP 3: annihilate matrix C  ====== */

        for (iColumn = 0; iColumn < nPivotRows; ++iColumn) {
            List<ColumnElement> column = cColumns[iColumn];
            // annihilate each element in the column
            for (int i = column.size() - 1; i >= 0; --i) {
                ColumnElement cElement = column.get(i);
                dMatrix.rows[cElement.iRow].subtract(bMatrix.rows[iColumn], cElement.value);
            }
        }

        /* ======  STEP 4: compute row reduced echelon form of matrix D  ====== */

        // this can be optimized (use dense structures, use same structured elimination), but actually
        // it doesn't take too much time

        // make D maximally triangular
        Arrays.sort(dMatrix.rows, Comparator.comparingInt(SparseArray::firstNonZeroPosition));
        dMatrix.rowReduce();

        /* ======  STEP 5: row reduce B  ====== */

        int dShift = 0;
        for (iRow = 0; iRow < dMatrix.nRows; iRow++) {
            SparseArray dRow = dMatrix.rows[iRow];
            iColumn = iRow + dShift;
            if (iColumn >= nColumns)
                break;
            if (dRow.coefficient(iColumn) == 0) {
                --iRow;
                ++dShift;
                continue;
            }
            for (int i = 0; i < bMatrix.nRows; ++i)
                if (bMatrix.rows[i].coefficient(iColumn) != 0)
                    bMatrix.rows[i].subtract(dRow, bMatrix.rows[i].coefficient(iColumn));
        }

        /* ======  STEP 6: finally form N+ polynomials ====== */

        // leading monomials of H-polynomials
        TreeSet<DegreeVector> hLeadMonomials = hPolynomials.stream().map(p -> p.hPoly.lt())
                .collect(Collectors.toCollection(() -> new TreeSet<>(factory.ordering)));

        List<ArrayBasedPoly<MonomialZp64>> nPolynomials = new ArrayList<>();
        // collect from A and B
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayList<MonomialZp64> candidateList = new ArrayList<>();

            if (iRow < nPivotRows)
                candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[iRow]], 1L));

            SparseArray row = iRow < nPivotRows ? bMatrix.rows[iRow] : dMatrix.rows[iRow - nPivotRows];
            for (int i = 0; i < row.sparsePositions.length; i++) {
                long val = row.sparseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + row.sparsePositions[i]]], val));
            }

            for (int i = 0; i < bDenseColumns.length; i++) {
                long val = row.denseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + bDenseColumns[i]]], val));
            }

            candidateList.sort(reverseOrder);
            if (!candidateList.isEmpty()) {
                ArrayBasedPoly<MonomialZp64> poly = new ArrayBasedPoly<>(
                        candidateList.toArray(factory.monomialAlgebra.createArray(candidateList.size())),
                        factory.nVariables);
                if (poly.lt().coefficient != 1) {
                    long factor = ring.reciprocal(poly.lt().coefficient);
                    for (int i = 0; i < poly.data.length; i++)
                        poly.data[i] = poly.data[i].setCoefficient(ring.multiply(factor, poly.data[i].coefficient));
                }
                assert poly.lt().coefficient == 1;
                nPolynomials.add(poly);
            }
        }
        nPolynomials.sort((a, b) -> reverseOrder.compare(a.lt(), b.lt()));
        // resulting N+ set
        List<ArrayBasedPoly<MonomialZp64>> nPlusPolynomials = new ArrayList<>();
        for (ArrayBasedPoly<MonomialZp64> candidate : nPolynomials) {
            if (!hLeadMonomials.contains(candidate.lt())) {
                // lt is new -> just add
                nPlusPolynomials.add(candidate);
                f4reductions.add(new ArrayList<>(Collections.singletonList(candidate)));
            } else
                // update f4reductions
                for (int iIndex = 0; iIndex < basis.size(); ++iIndex) {
                    ArrayBasedPoly<MonomialZp64> g = basis.get(iIndex);
                    if (g == null)
                        continue;
                    if (candidate.lt().dvDivisibleBy(g.lt())) {
                        List<ArrayBasedPoly<MonomialZp64>> reductions = f4reductions.get(iIndex);
                        boolean reduced = false;
                        for (int i = 0; i < reductions.size(); ++i) {
                            ArrayBasedPoly<MonomialZp64> red = reductions.get(i);
                            if (red.lt().dvEquals(candidate.lt())) {
                                reductions.set(i, candidate);
                                reduced = true;
                                break;
                            }
                        }
                        if (!reduced)
                            reductions.add(candidate);
                    }
                }
        }

        return nPlusPolynomials;
    }

    // for Zp64
    @SuppressWarnings("unchecked")
    private static List<ArrayBasedPoly<MonomialZp64>> nPlusZp64(
            MultivariatePolynomialZp64 factory,
            List<ArrayBasedPoly<MonomialZp64>> basis,
            List<HPolynomial<MonomialZp64>> hPolynomials,
            DegreeVector[] hMonomials,
            List<List<ArrayBasedPoly<MonomialZp64>>> f4reductions) {
        IntegersZp64 ring = factory.ring;

        // We use the structured Gaussian elimination strategy as described in
        // J.-C. Faugere & S. Lachartre, PASCO'10 https://doi.org/10.1145/1837210.1837225
        // "Parallel Gaussian Elimination for Gröbner bases computations in finite fields"
        //
        // By swapping the rows and non-pivoting columns, each matrix in F4
        // can be rewritten in the following form (F4-form):
        //
        //    \x0x00xx000x | x0x00xx0000x0x0x0xx000xx000xx0
        //    0\x0x0x00xx0 | 0xx0000x0xxx0x00xx0x0000xx0000
        //    00\x00x000x0 | 0x0x0000x0x0xx00xx0x00x0x0xx00
        //    000\xx0x0x00 | xx0xxx00x000x0x0xx00x0x0xx000x
        //    0000\xx0x0x0 | x0000xx0x00x0xxx0xx0000x000xx0
        //    00000\x0000x | 00x0000x0x0x0xx0xx0xx000xx0000
        //    000000\xx00x | 0x0x000x00x0xxx0xx00xxx0x0xx00
        //    0000000\x0x0 | xx00xx00xx00x000xx0xx00x0x000x
        //    ............ | ..............................
        //    -------------+-------------------------------
        //    0xx000x0x0xx | xxxxxx0xxxxxxx0xxxxxxxxxxxxxxx
        //    x0xx000x0x00 | xxxx0xxxxxxxxxxxxxx0xxxxxxxxxx
        //    00x00x0000xx | xxxxxxx0xxxxxxxxxxxxxxx0xxxxxx
        //    x0000x00xx0x | xxxxxxxxxxxxxxxxx0xxxxxxx0xxxx
        //    ............ | ..............................
        //
        // We denote:
        //
        // A - upper left  block (very sparse, triangular)         -- pivoting rows
        // B - upper right block (partially sparse, rectangular)   -- pivoting rows
        // C -  down left  block (partially  sparse, rectangular)  -- non-pivoting rows
        // D -  down right block (dense, rectangular)              -- non-pivoting rows
        //
        // The algorithm to reduce the matrix is then very simple:
        //
        // 1) row reduce A (B is still partially sparse)
        // 2) annihilate C (D is now almost certainly dense)
        // 3) row echelon & row reduce D
        // 4) row reduce B

        // reverse order for binary searching
        Comparator<DegreeVector> reverseOrder = (a, b) -> factory.ordering.compare(b, a);
        int
                nRows = hPolynomials.size(),
                nColumns = hMonomials.length,
                iRow, iColumn;

        /* ======  STEP 0: bring matrix to F4-form  ====== */

        // number of lead-terms in each column (columns with no any lead
        // terms are non-pivoting and can be rearranged)
        int[] columnsLeadTermsFilling = new int[nColumns];
        // detect non-pivoting columns
        int iOldColumnPrev = 0;
        for (HPolynomial<MonomialZp64> hPoly : hPolynomials) {
            iColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, hPoly.hPoly.lt(), reverseOrder);
            iOldColumnPrev = iColumn;
            assert iColumn >= 0;
            ++columnsLeadTermsFilling[iColumn];
        }

        // find non pivoting columns
        TIntArrayList nonPivotColumns = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows + nonPivotColumns.size() && iColumn < nColumns; ++iColumn)
            if (columnsLeadTermsFilling[iColumn] == 0)
                nonPivotColumns.add(iColumn);
        // now we move non-pivoting columns to the right
        // mapping between old and new columns numeration
        int[] nonPivotColumnsArr = nonPivotColumns.toArray();
        int[] columnsRearrangement = ArraysUtil.addAll(
                ArraysUtil.intSetDifference(ArraysUtil.sequence(0, nColumns), nonPivotColumnsArr),
                nonPivotColumnsArr);

        // back mapping between new and old columns numeration
        int[] columnsBackRearrangement = new int[nColumns];
        for (int i = 0; i < nColumns; ++i)
            columnsBackRearrangement[columnsRearrangement[i]] = i;

        // number of non-zero entries in each column
        int[] columnsFilling = new int[nColumns];
        // index of each term in hPolynomials in the hMonomials array
        int[][] mapping = new int[nRows][];
        // estimated row filling of B matrix (number of nonzero elements in each row)
        int[] bRowsFilling = new int[nRows];
        // first iteration: gather info about matrix pattern
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;
            mapping[iRow] = new int[hPoly.size()];
            iOldColumnPrev = 0;
            for (int i = 0; i < hPoly.size(); ++i) {
                MonomialZp64 term = hPoly.get(i);
                // column in old numeration
                int iOldColumn = Arrays.binarySearch(hMonomials, iOldColumnPrev, hMonomials.length, term, reverseOrder);
                iOldColumnPrev = iOldColumn;
                assert iOldColumn >= 0;
                // column in new numeration
                iColumn = columnsBackRearrangement[iOldColumn];
                int rrr = ArraysUtil.firstIndexOf(iColumn, mapping[iRow]);
                if (iColumn != 0 && rrr >= 0) {
                    System.out.println(hPoly.get(rrr));
                    System.out.println(hPoly.get(i));
                    System.out.println(Arrays.toString(hPoly.data));
                }
                assert iColumn == 0 || ArraysUtil.firstIndexOf(iColumn, mapping[iRow]) < 0;
                mapping[iRow][i] = iColumn;


                ++columnsFilling[iColumn];
                if (iColumn >= nRows)
                    ++bRowsFilling[iRow];
            }
        }

        // choose pivoting rows, so that B matrix is maximally sparse (rows with minimal bFillIns)
        TIntArrayList pivots = new TIntArrayList();
        for (iColumn = 0; iColumn < nRows; ++iColumn) {
            int minFillIn = Integer.MAX_VALUE;
            int pivot = -1;
            for (iRow = iColumn; iRow < nRows; ++iRow)
                if (mapping[iRow][0] == iColumn && bRowsFilling[iRow] < minFillIn) {
                    minFillIn = bRowsFilling[iRow];
                    pivot = iRow;
                } else if (pivot != -1 && mapping[iRow][0] != iColumn)
                    break;
            if (pivot == -1)
                break;
            pivots.add(pivot);
        }

        // rearrange rows: move pivots up and non-pivots down
        int nPivotRows = pivots.size();
        for (int i = 0; i < nPivotRows; ++i) {
            int pivot = pivots.get(i);
            Collections.swap(hPolynomials, i, pivot);
            ArraysUtil.swap(mapping, i, pivot);
            ArraysUtil.swap(bRowsFilling, i, pivot);
        }
        // <- the matrix is now is in the desired F4-form

        /* ======  STEP 1: prepare data structures  ====== */

        // dense columns in matrices B & D
        int[] bDenseColumns = IntStream.range(nPivotRows, nColumns)
                .filter(i -> 1.0 * columnsFilling[i] / nRows > DENSE_FILLING_THRESHOLD)
                .map(i -> i - nPivotRows)
                .toArray(); // <- it is sorted (for below binary searching)

        int nbdColumns = nColumns - nPivotRows;

        // dense rows in matrices B & D
        int[] bdDenseRows = IntStream.range(0, nRows)
                .filter(i -> 1.0 * bRowsFilling[i] / nbdColumns > DENSE_FILLING_THRESHOLD)
                .toArray(); // <- it is sorted (for below binary searching)

        // we represent A & C as a list of columns
        // column lists for A & C
        List<ColumnElement>[]
                aColumns = new List[nPivotRows],
                cColumns = new List[nPivotRows];
        for (int i = 0; i < nPivotRows; ++i) {
            aColumns[i] = new ArrayList<>();
            cColumns[i] = new ArrayList<>();
        }

        TIntArrayList[] bdSparseRows = new TIntArrayList[nbdColumns];
        TLongArrayList[] bdSparseRowsValues = new TLongArrayList[nbdColumns];
        for (int i = 0; i < nbdColumns; ++i) {
            bdSparseRows[i] = new TIntArrayList();
            bdSparseRowsValues[i] = new TLongArrayList();
        }
        long[][] bdDenseValues = new long[nbdColumns][bdDenseRows.length];

        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayBasedPoly<MonomialZp64> hPoly = hPolynomials.get(iRow).hPoly;

            int iDenseRow = Arrays.binarySearch(bdDenseRows, iRow);
            for (int i = 0; i < hPoly.size(); i++) {
                iColumn = mapping[iRow][i];
                long coefficient = hPoly.get(i).coefficient;
                if (iColumn < nPivotRows) {
                    // element of matrix A or C
                    if (iRow < nPivotRows)
                        aColumns[iColumn].add(new ColumnElement(iRow, coefficient));
                    else
                        cColumns[iColumn].add(new ColumnElement(iRow - nPivotRows, coefficient));
                } else {
                    // element of matrix B or D
                    iColumn -= nPivotRows;
                    if (iDenseRow >= 0)
                        bdDenseValues[iColumn][iDenseRow] = coefficient;
                    else {
                        bdSparseRows[iColumn].add(iRow);
                        bdSparseRowsValues[iColumn].add(coefficient);
                    }
                }
            }
        }

        // sparse B-D matrix
        SparseColumnMatrix bdMatrix = new SparseColumnMatrix(ring, nRows, nbdColumns, bdDenseRows);
        for (iColumn = 0; iColumn < bdMatrix.nColumns; iColumn++) {
            int[] sparsePositions = bdSparseRows[iColumn].toArray();
            long[] sparseValues = bdSparseRowsValues[iColumn].toArray();
            ArraysUtil.quickSort(sparsePositions, sparseValues);
            bdMatrix.columns[iColumn] = new SparseArray(ring, nRows, bdDenseRows, sparsePositions, bdDenseValues[iColumn], sparseValues);
        }

        /* ======  STEP 2: row reduce matrix A  ====== */

        // we start from the last column in matrix A
        for (iColumn = nPivotRows - 1; iColumn >= 0; --iColumn) {
            List<ColumnElement> column = aColumns[iColumn];
            // monicize diagonal
            ColumnElement diagonal = column.get(column.size() - 1);
            assert diagonal.iRow == iColumn;
            bdMatrix.multiplyRow(diagonal.iRow, ring.reciprocal(diagonal.value));

            long[] pivot = new long[bdMatrix.nColumns];
            for (int jColumn = 0; jColumn < pivot.length; ++jColumn)
                pivot[jColumn] = bdMatrix.columns[jColumn].coefficient(diagonal.iRow);

            // annihilate each non-diagonal element
            for (int i = column.size() - 2; i >= 0; --i) {
                ColumnElement nonDiagonal = column.get(i);
                assert nonDiagonal.value != 0;
                for (int jColumn = 0; jColumn < bdMatrix.nColumns; ++jColumn)
                    bdMatrix.columns[jColumn].subtract(nonDiagonal.iRow, ring.multiply(pivot[jColumn], nonDiagonal.value));
            }
        }


        /* ======  STEP 3: annihilate matrix C  ====== */

        for (iColumn = nPivotRows - 1; iColumn >= 0; --iColumn) {
            List<ColumnElement> column = cColumns[iColumn];
            long[] pivot = new long[bdMatrix.nColumns];
            for (int jColumn = 0; jColumn < pivot.length; ++jColumn)
                pivot[jColumn] = bdMatrix.columns[jColumn].coefficient(iColumn);
            // annihilate each element in the column
            for (int i = column.size() - 1; i >= 0; --i) {
                ColumnElement nonDiagonal = column.get(i);
                assert nonDiagonal.value != 0;
                for (int jColumn = 0; jColumn < bdMatrix.nColumns; ++jColumn)
                    bdMatrix.columns[jColumn].subtract(nPivotRows + nonDiagonal.iRow, ring.multiply(pivot[jColumn], nonDiagonal.value));
            }
        }

        SparseRowMatrix bdRowMatrix = bdMatrix.toRowMatrix(bDenseColumns);
        SparseRowMatrix bMatrix = bdRowMatrix.range(0, nPivotRows);
        SparseRowMatrix dMatrix = bdRowMatrix.range(nPivotRows, nRows);

//        System.out.println(prettyMatrix(bMatrix.denseMatrix()));
//        System.out.println(prettyMatrix(dMatrix.denseMatrix()));
        /* ======  STEP 4: compute row reduced echelon form of matrix D  ====== */

        // this can be optimized (use dense structures, use same structured elimination), but actually
        // it doesn't take too much time

        // make D maximally triangular
        Arrays.sort(dMatrix.rows, Comparator.comparingInt(SparseArray::firstNonZeroPosition));
        dMatrix.rowReduce();

        /* ======  STEP 5: row reduce B  ====== */

        int dShift = 0;
        for (iRow = 0; iRow < dMatrix.nRows; iRow++) {
            SparseArray dRow = dMatrix.rows[iRow];
            iColumn = iRow + dShift;
            if (iColumn >= nColumns)
                break;
            if (dRow.coefficient(iColumn) == 0) {
                --iRow;
                ++dShift;
                continue;
            }
            for (int i = 0; i < bMatrix.nRows; ++i)
                if (bMatrix.rows[i].coefficient(iColumn) != 0)
                    bMatrix.rows[i].subtract(dRow, bMatrix.rows[i].coefficient(iColumn));
        }

        /* ======  STEP 6: finally form N+ polynomials ====== */

        // leading monomials of H-polynomials
        TreeSet<DegreeVector> hLeadMonomials = hPolynomials.stream().map(p -> p.hPoly.lt())
                .collect(Collectors.toCollection(() -> new TreeSet<>(factory.ordering)));

        List<ArrayBasedPoly<MonomialZp64>> nPolynomials = new ArrayList<>();
        // collect from A and B
        for (iRow = 0; iRow < nRows; ++iRow) {
            ArrayList<MonomialZp64> candidateList = new ArrayList<>();

            if (iRow < nPivotRows)
                candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[iRow]], 1L));

            SparseArray row = iRow < nPivotRows ? bMatrix.rows[iRow] : dMatrix.rows[iRow - nPivotRows];
            for (int i = 0; i < row.sparsePositions.length; i++) {
                long val = row.sparseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + row.sparsePositions[i]]], val));
            }

            for (int i = 0; i < bDenseColumns.length; i++) {
                long val = row.denseValues[i];
                if (val != 0)
                    candidateList.add(new MonomialZp64(hMonomials[columnsRearrangement[nPivotRows + bDenseColumns[i]]], val));
            }

            candidateList.sort(reverseOrder);
            if (!candidateList.isEmpty()) {
                ArrayBasedPoly<MonomialZp64> poly = new ArrayBasedPoly<>(
                        candidateList.toArray(factory.monomialAlgebra.createArray(candidateList.size())),
                        factory.nVariables);
                if (poly.lt().coefficient != 1) {
                    long factor = ring.reciprocal(poly.lt().coefficient);
                    for (int i = 0; i < poly.data.length; i++)
                        poly.data[i] = poly.data[i].setCoefficient(ring.multiply(factor, poly.data[i].coefficient));
                }
                assert poly.lt().coefficient == 1;
                nPolynomials.add(poly);
            }
        }
        nPolynomials.sort((a, b) -> reverseOrder.compare(a.lt(), b.lt()));
        // resulting N+ set
        List<ArrayBasedPoly<MonomialZp64>> nPlusPolynomials = new ArrayList<>();
        for (ArrayBasedPoly<MonomialZp64> candidate : nPolynomials) {
            if (!hLeadMonomials.contains(candidate.lt())) {
                // lt is new -> just add
                nPlusPolynomials.add(candidate);
                f4reductions.add(new ArrayList<>(Collections.singletonList(candidate)));
            } else
                // update f4reductions
                for (int iIndex = 0; iIndex < basis.size(); ++iIndex) {
                    ArrayBasedPoly<MonomialZp64> g = basis.get(iIndex);
                    if (g == null)
                        continue;
                    if (candidate.lt().dvDivisibleBy(g.lt())) {
                        List<ArrayBasedPoly<MonomialZp64>> reductions = f4reductions.get(iIndex);
                        boolean reduced = false;
                        for (int i = 0; i < reductions.size(); ++i) {
                            ArrayBasedPoly<MonomialZp64> red = reductions.get(i);
                            if (red.lt().dvEquals(candidate.lt())) {
                                reductions.set(i, candidate);
                                reduced = true;
                                break;
                            }
                        }
                        if (!reduced)
                            reductions.add(candidate);
                    }
                }
        }

        return nPlusPolynomials;
    }

    static final class SparseArray {
        final IntegersZp64 ring;
        final int length;
        final int[] densePositions;
        final long[] denseValues;
        int[] sparsePositions;
        long[] sparseValues;

        SparseArray(IntegersZp64 ring,
                    int length,
                    int[] densePositions,
                    int[] sparsePositions,
                    long[] denseValues,
                    long[] sparseValues) {
            this.ring = ring;
            this.length = length;
            this.densePositions = densePositions;
            this.sparsePositions = sparsePositions;
            this.denseValues = denseValues;
            this.sparseValues = sparseValues;
            assert isSorted(sparsePositions);
        }

        SparseArray(IntegersZp64 ring,
                    int[] densePositions,
                    long[] denseArray) {
            this.ring = ring;
            this.densePositions = densePositions;
            this.length = denseArray.length;
            this.denseValues = new long[densePositions.length];
            for (int i = 0; i < densePositions.length; i++)
                denseValues[i] = denseArray[densePositions[i]];

            TIntArrayList sparsePositions = new TIntArrayList();
            TLongArrayList sparseValues = new TLongArrayList();
            for (int i = 0; i < denseArray.length; ++i) {
                if (denseArray[i] == 0)
                    continue;
                if (Arrays.binarySearch(densePositions, i) >= 0)
                    continue;
                sparsePositions.add(i);
                sparseValues.add(denseArray[i]);
            }
            this.sparsePositions = sparsePositions.toArray();
            this.sparseValues = sparseValues.toArray();
        }

        long[] toDense() {
            long[] result = new long[length];
            for (int i = 0; i < densePositions.length; ++i)
                result[densePositions[i]] = denseValues[i];
            for (int i = 0; i < sparsePositions.length; ++i)
                result[sparsePositions[i]] = sparseValues[i];
            return result;
        }

        void subtract(SparseArray pivot, long factor, int iColumn, int firstDense) {
            if (factor == 0)
                return;
            assert pivot.densePositions == densePositions;
            long negFactor = ring.negate(factor);

            // adding dense parts
            for (int i = firstDense; i < denseValues.length; i++)
                denseValues[i] = ring.subtract(denseValues[i], ring.multiply(factor, pivot.denseValues[i]));

            // subtracting sparse parts
            int[] pivCols = pivot.sparsePositions;
            long[] pivVals = pivot.sparseValues;

            int firstSparse = ArraysUtil.binarySearch1(sparsePositions, iColumn);

            // resulting non-zero columns
            TIntArrayList resCols = new TIntArrayList(sparsePositions.length + pivCols.length);
            TLongArrayList resVals = new TLongArrayList(sparsePositions.length + pivCols.length);

            resCols.add(sparsePositions, 0, firstSparse);
            resVals.add(sparseValues, 0, firstSparse);

            int iSel = firstSparse, iPiv = 0;
            while (iSel < sparsePositions.length && iPiv < pivCols.length) {
                int
                        selCol = sparsePositions[iSel],
                        othCol = pivCols[iPiv];

                if (selCol == othCol) {
                    long subtract = ring.subtract(sparseValues[iSel], ring.multiply(factor, pivVals[iPiv]));
                    if (subtract != 0) {
                        resCols.add(selCol);
                        resVals.add(subtract);
                    }

                    ++iSel;
                    ++iPiv;
                } else if (selCol < othCol) {
                    resCols.add(selCol);
                    resVals.add(sparseValues[iSel]);

                    ++iSel;
                } else if (selCol > othCol) {
                    resCols.add(othCol);
                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));

                    ++iPiv;
                }
            }

            if (iSel < sparsePositions.length)
                for (; iSel < sparsePositions.length; ++iSel) {
                    resCols.add(sparsePositions[iSel]);
                    resVals.add(sparseValues[iSel]);
                }

            if (iPiv < pivCols.length)
                for (; iPiv < pivCols.length; ++iPiv) {
                    resCols.add(pivCols[iPiv]);
                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));
                }

            sparsePositions = resCols.toArray();
            sparseValues = resVals.toArray();
        }

        void subtract(SparseArray pivot, long factor) {
            subtract(pivot, factor, 0, 0);
        }

        void multiply(long factor) {
            if (factor == 1)
                return;
            for (int i = 0; i < sparseValues.length; ++i)
                sparseValues[i] = ring.multiply(sparseValues[i], factor);
            for (int i = 0; i < denseValues.length; ++i)
                denseValues[i] = ring.multiply(denseValues[i], factor);
        }

        long coefficient(int iRow) {
            int index = Arrays.binarySearch(sparsePositions, iRow);
            if (index >= 0)
                return sparseValues[index];
            index = Arrays.binarySearch(densePositions, iRow);
            if (index >= 0)
                return denseValues[index];
            return 0;
        }

        void multiply(int iRow, long value) {
            int index = Arrays.binarySearch(sparsePositions, iRow);
            if (index >= 0)
                sparseValues[index] = ring.multiply(sparseValues[index], value);
            else {
                index = Arrays.binarySearch(densePositions, iRow);
                if (index >= 0)
                    denseValues[index] = ring.multiply(denseValues[index], value);
            }
        }

        void subtract(int iRow, long value) {
            if (value == 0)
                return;

            int index = Arrays.binarySearch(densePositions, iRow);
            if (index >= 0)
                denseValues[index] = ring.subtract(denseValues[index], value);
            else {
                index = Arrays.binarySearch(sparsePositions, iRow);
                if (index >= 0) {
                    sparseValues[index] = ring.subtract(sparseValues[index], value);
                    if (sparseValues[index] == 0) {
                        sparsePositions = ArraysUtil.remove(sparsePositions, index);
                        sparseValues = ArraysUtil.remove(sparseValues, index);
                    }
                } else {
                    index = ~index;
//                    TIntHashSet rr = new TIntHashSet(sparsePositions);
//                    rr.add(iRow);
                    sparsePositions = ArraysUtil.insert(sparsePositions, index, iRow);
//                    assert rr.equals(new TIntHashSet(sparsePositions));

//                    TLongHashSet ll = new TLongHashSet(sparseValues);
//                    ll.add(ring.negate( value));
                    sparseValues = ArraysUtil.insert(sparseValues, index, ring.negate(value));
//                    assert ll.equals(new TLongHashSet(sparseValues));

//                    assert sparseValues[Arrays.binarySearch(sparsePositions, iRow)] == ring.negate(value);
//                    assert isSorted(sparsePositions);
                }
            }
        }

        int firstNonZeroPosition() {
            int firstSparse = sparsePositions.length != 0 ? sparsePositions[0] : Integer.MAX_VALUE;
            for (int i = 0; i < densePositions.length; ++i)
                if (denseValues[i] != 0)
                    return Math.min(densePositions[i], firstSparse);
            return firstSparse;
        }
    }

    static final class SparseColumnMatrix {
        final IntegersZp64 ring;
        final int nRows, nColumns;
        final int[] densePositions;
        final SparseArray[] columns;

        SparseColumnMatrix(IntegersZp64 ring, int nRows, int nColumns, int[] densePositions) {
            this.ring = ring;
            this.nRows = nRows;
            this.nColumns = nColumns;
            this.densePositions = densePositions;
            this.columns = new SparseArray[nColumns];
        }

        void multiplyRow(int iRow, long value) {
            if (value != 1)
                for (int i = 0; i < nColumns; ++i)
                    columns[i].multiply(iRow, value);
        }

        SparseRowMatrix toRowMatrix(int[] densePositions) {
            SparseRowMatrix rowMatrix = new SparseRowMatrix(ring, nRows, nColumns, densePositions);
            TIntArrayList[] sparseColumns = new TIntArrayList[nRows];
            TLongArrayList[] sparseValues = new TLongArrayList[nRows];
            long[][] denseValues = new long[nRows][densePositions.length];

            for (int iRow = 0; iRow < nRows; ++iRow) {
                sparseColumns[iRow] = new TIntArrayList();
                sparseValues[iRow] = new TLongArrayList();
            }

            for (int iColumn = 0; iColumn < nColumns; ++iColumn) {
                SparseArray column = columns[iColumn];
                int iDenseColumn = Arrays.binarySearch(densePositions, iColumn);
                if (iDenseColumn >= 0) {
                    for (int i = 0; i < column.densePositions.length; i++)
                        denseValues[column.densePositions[i]][iDenseColumn] = column.denseValues[i];
                    for (int i = 0; i < column.sparsePositions.length; ++i)
                        denseValues[column.sparsePositions[i]][iDenseColumn] = column.sparseValues[i];
                } else {
                    for (int i = 0; i < column.densePositions.length; i++) {
                        sparseColumns[column.densePositions[i]].add(iColumn);
                        sparseValues[column.densePositions[i]].add(column.denseValues[i]);
                    } for (int i = 0; i < column.sparsePositions.length; ++i) {
                        sparseColumns[column.sparsePositions[i]].add(iColumn);
                        sparseValues[column.sparsePositions[i]].add(column.sparseValues[i]);
                    }
                }
            }

            for (int iRow = 0; iRow < nRows; ++iRow)
                rowMatrix.rows[iRow] = new SparseArray(ring, nColumns, densePositions, sparseColumns[iRow].toArray(), denseValues[iRow], sparseValues[iRow].toArray());

            return rowMatrix;
        }
    }

    static final class SparseRowMatrix {
        final IntegersZp64 ring;
        final int nRows, nColumns;
        final int[] densePositions;
        final SparseArray[] rows;

        SparseRowMatrix(IntegersZp64 ring, int nRows, int nColumns, int[] densePositions) {
            this.ring = ring;
            this.nRows = nRows;
            this.nColumns = nColumns;
            this.densePositions = densePositions;
            this.rows = new SparseArray[nRows];
        }

        SparseRowMatrix(IntegersZp64 ring, int nRows, int nColumns, int[] densePositions, SparseArray[] rows) {
            this.ring = ring;
            this.nRows = nRows;
            this.nColumns = nColumns;
            this.densePositions = densePositions;
            this.rows = rows;
        }

        SparseRowMatrix range(int from, int to) {
            return new SparseRowMatrix(ring, to - from, nColumns, densePositions, Arrays.copyOfRange(rows, from, to));
        }

        long[][] denseMatrix() {
            long[][] result = new long[rows.length][nColumns];
            for (int iRow = 0; iRow < rows.length; ++iRow) {
                SparseArray row = rows[iRow];
                for (int i = 0; i < row.sparsePositions.length; i++)
                    result[iRow][row.sparsePositions[i]] = row.sparseValues[i];
                for (int i = 0; i < densePositions.length; i++)
                    result[iRow][densePositions[i]] = row.denseValues[i];
            }
            return result;
        }


        // return rank
        int rowReduce() {
            // Gaussian elimination
            int nZeroRows = 0;
            for (int iCol = 0, to = rows.length; iCol < to; ++iCol) {
                int iRow = iCol - nZeroRows;
                int pivotIndex = -1;

                for (int i = iRow; i < rows.length; ++i) {
                    if (rows[i].coefficient(iCol) == 0)
                        continue;
                    if (pivotIndex == -1) {
                        pivotIndex = i;
                        continue;
                    }
                    if (rows[i].sparsePositions.length < rows[pivotIndex].sparsePositions.length)
                        pivotIndex = i;
                }

                if (pivotIndex == -1) {
                    ++nZeroRows;
                    to = Math.min(nColumns, rows.length + nZeroRows);
                    continue;
                }

                ArraysUtil.swap(rows, pivotIndex, iRow);
                SparseArray pivot = rows[iRow];
                long diagonalValue = pivot.coefficient(iCol);

                // row-reduction
                pivot.multiply(ring.reciprocal(diagonalValue));

                int firstDense = ArraysUtil.binarySearch1(densePositions, iCol);
                for (int row = 0; row < rows.length; ++row) {
                    if (row == iRow)
                        continue;

                    long value = rows[row].coefficient(iCol);
                    if (value == 0)
                        continue;

                    rows[row].subtract(pivot, value, iCol, firstDense);
                }
            }
            return Math.min(nRows, nColumns) - nZeroRows;
        }
    }

    static final class SparseColumn {
        final int[] rows;
        final long[] values;

        SparseColumn(int[] rows, long[] values) {
            this.rows = rows;
            this.values = values;
        }
    }

    static final class ColumnElement {
        final int iRow;
        final long value;

        ColumnElement(int iRow, long value) { this.iRow = iRow; this.value = value; }

        @Override
        public String toString() { return "(" + iRow + " -> " + value + ")"; }
    }

    private static boolean isSorted(int[] arr) {
        for (int i = 1; i < arr.length; ++i)
            if (arr[i - 1] >= arr[i])
                return false;
        return true;
    }
//
//    static final class SparseMatrixZp64 {
//        final IntegersZp64 ring;
//        // indexes of dense columns
//        final int[] denseColumns;
//        // the rows (initially nulls)
//        final SparseRow[] matrix;
//        final int nRows, nColumns;
//
//        SparseMatrixZp64(IntegersZp64 ring, int[] denseColumns, int nRows, int nColumns) {
//            this.ring = ring;
//            this.denseColumns = denseColumns;
//            this.matrix = new SparseRow[nRows];
//            this.nRows = nRows;
//            this.nColumns = nColumns;
//        }
//
//        /* copy constructor */
//        SparseMatrixZp64(IntegersZp64 ring, int[] denseColumns, SparseRow[] matrix, int nRows, int nColumns) {
//            this.ring = ring;
//            this.denseColumns = denseColumns;
//            this.matrix = matrix;
//            this.nRows = nRows;
//            this.nColumns = nColumns;
//        }
//
//        long[][] denseMatrix() {
//            long[][] result = new long[matrix.length][nColumns];
//            for (int iRow = 0; iRow < matrix.length; ++iRow) {
//                SparseRow row = this.matrix[iRow];
//                for (int i = 0; i < row.sparseColumns.length; i++)
//                    result[iRow][row.sparseColumns[i]] = row.sparse[i];
//                for (int i = 0; i < denseColumns.length; i++)
//                    result[iRow][denseColumns[i]] = row.dense[i];
//            }
//            return result;
//        }
//
//        void clear() {
//            for (SparseRow row : matrix)
//                row.clear();
//        }
//
//        void normalize() {
//            for (SparseRow row : matrix)
//                row.normalize();
//        }
//
//        final class SparseRow {
//            // non-zero columns
//            int[] sparseColumns;
//            // values in non-zero columns
//            long[] sparse;
//            // values in dense columns
//            long[] dense;
//
//            SparseRow(int[] sparseColumns, long[] sparse, long[] dense) {
//                assert isSorted(sparseColumns);
//                assert sparse.length == 0 || LongStream.of(sparse).noneMatch(s -> s == 0);
//                this.sparseColumns = sparseColumns;
//                this.sparse = sparse;
//                this.dense = dense;
//            }
//
//            void clear() {
//                Arrays.fill(dense, 0);
//                Arrays.fill(sparse, 0);
//            }
//
//            void normalize() {
//                TIntArrayList zeros = new TIntArrayList();
//                for (int i = 0; i < sparse.length; ++i)
//                    if (sparse[i] == 0)
//                        zeros.add(i);
//                int[] arr = zeros.toArray();
//                sparseColumns = ArraysUtil.remove(sparseColumns, arr);
//                sparse = ArraysUtil.remove(sparse, arr);
//            }
//
//            void setDenseValue(int iDenseColumn, long value) {
//                dense[iDenseColumn] = value;
//            }
//
//            void setSparseValue(int iColumn, long value) {
//                int ind = Arrays.binarySearch(sparseColumns, iColumn);
//                if (ind >= 0)
//                    sparse[ind] = value;
//                else {
//                    ind = ~ind;
//                    sparseColumns = ArraysUtil.insert(sparseColumns, ind, iColumn);
//                    sparse = ArraysUtil.insert(sparse, ind, value);
//                    assert isSorted(sparseColumns);
//                }
//            }
//
//            void multiply(long factor) {
//                if (factor == 1)
//                    return;
//                for (int i = 0; i < sparse.length; ++i)
//                    sparse[i] = ring.multiply(sparse[i], factor);
//                for (int i = 0; i < dense.length; ++i)
//                    dense[i] = ring.multiply(dense[i], factor);
//            }
//
//            long coefficient(int iColumn) {
//                int index = Arrays.binarySearch(sparseColumns, iColumn);
//                if (index >= 0)
//                    return sparse[index];
//                index = Arrays.binarySearch(denseColumns, iColumn);
//                if (index >= 0)
//                    return dense[index];
//                return 0;
//            }
//
//            int firstColumn() {
//                if (sparseColumns.length != 0)
//                    return sparseColumns[0];
//                for (int i = 0; i < dense.length; ++i)
//                    if (dense[i] != 0)
//                        return denseColumns[i];
//                return Integer.MAX_VALUE;
//            }
//
//            @Override
//            public SparseRow clone() {
//                return new SparseRow(sparseColumns.clone(), sparse.clone(), dense.clone());
//            }
//
//            @Override
//            public boolean equals(Object o) {
//                if (this == o) return true;
//                if (o == null || getClass() != o.getClass()) return false;
//
//                SparseRow sparseRow = (SparseRow) o;
//
//                if (!Arrays.equals(sparseColumns, sparseRow.sparseColumns)) return false;
//                if (!Arrays.equals(sparse, sparseRow.sparse)) return false;
//                return Arrays.equals(dense, sparseRow.dense);
//            }
//
//            @Override
//            public int hashCode() {
//                int result = Arrays.hashCode(sparseColumns);
//                result = 31 * result + Arrays.hashCode(sparse);
//                result = 31 * result + Arrays.hashCode(dense);
//                return result;
//            }
//
//        }
//
//        // return rank
//        int rowReduce() {
//            // Gaussian elimination
//            int nZeroRows = 0;
//            for (int iCol = 0, to = matrix.length; iCol < to; ++iCol) {
//                int iRow = iCol - nZeroRows;
//                int pivotIndex = -1;
//
//                for (int i = iRow; i < matrix.length; ++i) {
//                    if (matrix[i].coefficient(iCol) == 0)
//                        continue;
//                    if (pivotIndex == -1) {
//                        pivotIndex = i;
//                        continue;
//                    }
//                    if (matrix[i].sparse.length < matrix[pivotIndex].sparse.length)
//                        pivotIndex = i;
//                }
//
//                if (pivotIndex == -1) {
//                    ++nZeroRows;
//                    to = Math.min(nColumns, matrix.length + nZeroRows);
//                    continue;
//                }
//
//                ArraysUtil.swap(matrix, pivotIndex, iRow);
//                SparseRow pivot = matrix[iRow];
//                long diagonalValue = pivot.coefficient(iCol);
//
//                // row-reduction
//                pivot.multiply(ring.reciprocal(diagonalValue));
//
//                int firstDense = ArraysUtil.binarySearch1(denseColumns, iCol);
//                for (int row = 0; row < matrix.length; ++row) {
//                    if (row == iRow)
//                        continue;
//
//                    long value = matrix[row].coefficient(iCol);
//                    if (value == 0)
//                        continue;
//
//                    matrix[row] = subtract(matrix[row], pivot, value, iCol, firstDense);
//                }
//            }
//            return Math.min(nRows, nColumns) - nZeroRows;
//        }
//
//        /** self - factor * amount */
//        SparseRow subtract(SparseRow self, SparseRow pivot, long factor, int iColumn, int firstDense) {
//            long negFactor = ring.negate(factor);
//
//            // adding dense parts
//            long[]
//                    selDense = self.dense,
//                    pivDense = pivot.dense;
//
//            for (int i = firstDense; i < selDense.length; i++)
//                selDense[i] = ring.subtract(selDense[i], ring.multiply(factor, pivDense[i]));
//
//            // subtracting sparse parts
//            int[]
//                    selCols = self.sparseColumns,
//                    pivCols = pivot.sparseColumns;
//            long[]
//                    selVals = self.sparse,
//                    pivVals = pivot.sparse;
//
//            int firstSparse = ArraysUtil.binarySearch1(selCols, iColumn);
//
//            // resulting non-zero columns
//            TIntArrayList resCols = new TIntArrayList(selCols.length + pivCols.length);
//            TLongArrayList resVals = new TLongArrayList(selCols.length + pivCols.length);
//
//            resCols.add(selCols, 0, firstSparse);
//            resVals.add(selVals, 0, firstSparse);
//
//            int iSel = firstSparse, iPiv = 0;
//            while (iSel < selCols.length && iPiv < pivCols.length) {
//                int
//                        selCol = selCols[iSel],
//                        othCol = pivCols[iPiv];
//
//                if (selCol == othCol) {
//                    long subtract = ring.subtract(selVals[iSel], ring.multiply(factor, pivVals[iPiv]));
//                    if (subtract != 0) {
//                        resCols.add(selCol);
//                        resVals.add(subtract);
//                    }
//
//                    ++iSel;
//                    ++iPiv;
//                } else if (selCol < othCol) {
//                    resCols.add(selCol);
//                    resVals.add(selVals[iSel]);
//
//                    ++iSel;
//                } else if (selCol > othCol) {
//                    resCols.add(othCol);
//                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));
//
//                    ++iPiv;
//                }
//            }
//
//            if (iSel < selCols.length)
//                for (; iSel < selCols.length; ++iSel) {
//                    resCols.add(selCols[iSel]);
//                    resVals.add(selVals[iSel]);
//                }
//
//            if (iPiv < pivCols.length)
//                for (; iPiv < pivCols.length; ++iPiv) {
//                    resCols.add(pivCols[iPiv]);
//                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));
//                }
//
//            self.sparseColumns = resCols.toArray();
//            self.sparse = resVals.toArray();
//            return self;
//        }
//
//        /** self - factor * amount */
//        SparseRow subtract(SparseRow self, SparseRow pivot, long factor) {
//            long negFactor = ring.negate(factor);
//
//            // adding dense parts
//            long[]
//                    selDense = self.dense,
//                    pivDense = pivot.dense;
//
//            for (int i = 0; i < selDense.length; i++)
//                selDense[i] = ring.subtract(selDense[i], ring.multiply(factor, pivDense[i]));
//
//            // subtracting sparse parts
//            int[]
//                    selCols = self.sparseColumns,
//                    pivCols = pivot.sparseColumns;
//            long[]
//                    selVals = self.sparse,
//                    pivVals = pivot.sparse;
//
//            // resulting non-zero columns
//            TIntArrayList resCols = new TIntArrayList(selCols.length + pivCols.length);
//            TLongArrayList resVals = new TLongArrayList(selCols.length + pivCols.length);
//
//            int iSel = 0, iPiv = 0;
//            while (iSel < selCols.length && iPiv < pivCols.length) {
//                int
//                        selCol = selCols[iSel],
//                        othCol = pivCols[iPiv];
//
//                if (selCol == othCol) {
//                    long subtract = ring.subtract(selVals[iSel], ring.multiply(factor, pivVals[iPiv]));
//                    if (subtract != 0) {
//                        resCols.add(selCol);
//                        resVals.add(subtract);
//                    }
//
//                    ++iSel;
//                    ++iPiv;
//                } else if (selCol < othCol) {
//                    resCols.add(selCol);
//                    resVals.add(selVals[iSel]);
//
//                    ++iSel;
//                } else if (selCol > othCol) {
//                    resCols.add(othCol);
//                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));
//
//                    ++iPiv;
//                }
//            }
//
//            if (iSel < selCols.length)
//                for (; iSel < selCols.length; ++iSel) {
//                    resCols.add(selCols[iSel]);
//                    resVals.add(selVals[iSel]);
//                }
//
//            if (iPiv < pivCols.length)
//                for (; iPiv < pivCols.length; ++iPiv) {
//                    resCols.add(pivCols[iPiv]);
//                    resVals.add(ring.multiply(negFactor, pivVals[iPiv]));
//                }
//
//            return new SparseRow(resCols.toArray(), resVals.toArray(), selDense);
//        }
//
//        @Override
//        public SparseMatrixZp64 clone() {
//            SparseRow[] newMatrix = new SparseRow[nRows];
//            for (int i = 0; i < nRows; ++i)
//                newMatrix[i] = matrix[i].clone();
//            return new SparseMatrixZp64(ring, denseColumns.clone(), newMatrix, nRows, nColumns);
//        }
//
//        @Override
//        public boolean equals(Object o) {
//            if (this == o) return true;
//            if (o == null || getClass() != o.getClass()) return false;
//
//            SparseMatrixZp64 that = (SparseMatrixZp64) o;
//
//            if (nRows != that.nRows) return false;
//            if (nColumns != that.nColumns) return false;
//            if (!Arrays.equals(denseColumns, that.denseColumns)) return false;
//            // Probably incorrect - comparing Object[] arrays with Arrays.equals
//            return Arrays.equals(matrix, that.matrix);
//        }
//
//        @Override
//        public int hashCode() {
//            int result = Arrays.hashCode(denseColumns);
//            result = 31 * result + Arrays.hashCode(matrix);
//            result = 31 * result + nRows;
//            result = 31 * result + nColumns;
//            return result;
//        }
//    }

    public static String prettyMatrix(long[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = Long.toString(matrix[i][j]);
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.deleteCharAt(sb.length() - 1);
        return sb.append("\n}").toString();
    }


    private static String padding(char c, int len) {
        return new String(ArraysUtil.arrayOf(c, len));
    }

    private static String padd(String str, int newLen) {
        return padding(' ', newLen - str.length()) + str;
    }

    public static String prettyMatrix(Object[][] matrix) {
        int maxLength = 0;

        String[][] strings = new String[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            strings[i] = new String[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                strings[i][j] = matrix[i][j].toString();
                maxLength = Math.max(maxLength, strings[i][j].length());
            }
        }
        ++maxLength;
        for (int i = 0; i < matrix.length; i++)
            for (int j = 0; j < matrix[i].length; j++)
                strings[i][j] = padd(strings[i][j], maxLength);


        StringBuilder sb = new StringBuilder().append("{").append("\n");
        String sep = "    ";
        for (int i = 0; i < strings.length; i++) {
            sb.append(sep)
                    .append("{").append(Arrays.stream(strings[i]).collect(Collectors.joining(","))).append("}")
                    .append(",\n");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.deleteCharAt(sb.length() - 1);
        return sb.append("\n}").toString();
    }
}
