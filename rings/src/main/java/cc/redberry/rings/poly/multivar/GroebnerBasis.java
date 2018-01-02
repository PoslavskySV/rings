package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.Ring;
import cc.redberry.rings.linear.LinearSolver;
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> GroebnerBasis(List<Poly> generators,
                             Comparator<DegreeVector> monomialOrder) {
        return BuchbergerGB(generators, monomialOrder);
    }

    private static <Poly extends AMultivariatePolynomial<?, Poly>>
    boolean containsNonZeroConstant(List<Poly> generators) {
        return generators.stream().anyMatch(p -> p.isConstant() && !p.isZero());
    }

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

    /**
     * Computes minimized and reduced Groebner basis of a given homogenious ideal via Buchberger algorithm.
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerHomogeneousGB(List<Poly> generators,
                                       Comparator<DegreeVector> monomialOrder) {
        Comparator<SyzygyPair> selectionStrategy = normalSelectionStrategy(generators.get(0).ordering);
        // fixme use sugar always?
        if (!isGradedOrder(monomialOrder))
            // add sugar for non-graded orders
            selectionStrategy = withSugar(selectionStrategy);
        return BuchbergerHomogeneousGB(generators, monomialOrder, selectionStrategy);
    }

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via Buchberger algorithm.
     *
     * @param generators        generators of the ideal
     * @param monomialOrder     monomial order to use
     * @param selectionStrategy critical pair selection strategy
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> BuchbergerHomogeneousGB(List<Poly> generators,
                                       Comparator<DegreeVector> monomialOrder,
                                       Comparator<SyzygyPair> selectionStrategy) {
        assert generators.stream().allMatch(Poly::isHomogeneous);

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

        // degree -> pairs
        TreeMap<Integer, TreeSet<SyzygyPair<Term, Poly>>> sPairs = new TreeMap<>();
        TLongObjectHashMap<SyzygyPair<Term, Poly>> ijPairs = new TLongObjectHashMap<>();

        for (int i = 0; i < groebner.size() - 1; i++)
            for (int j = i + 1; j < groebner.size(); ++j) {
                SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(i, j, groebner);
                sPairs.computeIfAbsent(sPair.degree(), __ -> new TreeSet<>(selectionStrategy)).add(sPair);
                ijPairs.put(pack(i, j), sPair);
            }

        // cache array used in divisions (little performance improvement actually)
        Poly[] groebnerArray = groebner.toArray(factory.createArray(groebner.size()));

        while (!sPairs.isEmpty()) {
            // pick up (and remove) all pairs with the smallest degree
            main:
            for (SyzygyPair<Term, Poly> pair : sPairs.pollFirstEntry().getValue()) {
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
                        sPairs.computeIfAbsent(sPair.degree(), __ -> new TreeSet<>(selectionStrategy)).add(sPair);
                        ijPairs.put(sPair.hash(), sPair);
                    }
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

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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

    private static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
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
            Term extends AMonomial<Term>,
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
            this.syzygyGamma = fi.monomialAlgebra.createTermWithUnitCoefficient(ArraysUtil.max(fi.multidegree(), fj.multidegree()));
            this.sugar = Math.max(
                    fi.degreeSum() - fi.lt().totalDegree,
                    fj.degreeSum() - fj.lt().totalDegree) + syzygyGamma.totalDegree;
        }

        int degree() { return syzygyGamma.totalDegree;}

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

    /**
     * Computes reduced Groebner basis
     */
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

    /** Computes syzygy of given polynomials */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Poly a, Poly b) {
        return syzygy(a.monomialAlgebra.createTermWithUnitCoefficient(ArraysUtil.max(a.multidegree(), b.multidegree())), a, b);
    }

    static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Poly syzygy(Term xGamma, Poly a, Poly b) {
        IMonomialAlgebra<Term> mAlgebra = a.monomialAlgebra;
        Poly
                aReduced = a.clone().multiply(mAlgebra.divideExact(xGamma, a.lt())),
                bReduced = b.clone().multiply(mAlgebra.divideExact(xGamma, b.lt())),
                syzygy = aReduced.subtract(bReduced);
        return syzygy;
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
        return F4GB(generators, monomialOrder, normalSelectionStrategy(monomialOrder));
    }

    private static final int F4_MIN_SELECTION_SIZE = 16;

    /**
     * Computes minimized and reduced Groebner basis of a given ideal via F4 algorithm.
     *
     * @param generators        generators of the ideal
     * @param monomialOrder     monomial order to use
     * @param selectionStrategy critical pair selection strategy
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> F4GB(List<Poly> generators,
                    Comparator<DegreeVector> monomialOrder,
                    Comparator<SyzygyPair> selectionStrategy) {
        Poly factory = generators.get(0);
        IMonomialAlgebra<Term> mAlgebra = factory.monomialAlgebra;

        // stack of GB candidates
        List<Poly> groebner = prepareGenerators(generators, monomialOrder);
        if (groebner.size() == 1)
            return groebner;

        // sort polynomials in basis to achieve faster divisions
        Comparator<Poly> polyOrder = isGradedOrder(monomialOrder)
                ? (a, b) -> monomialOrder.compare(a.lt(), b.lt())
                : new EcartComparator<>();
        groebner.sort(polyOrder);

        // degree -> pairs (not actually computed)
        TreeMap<Integer, TreeSet<SyzygyPair<Term, Poly>>> sPairs = new TreeMap<>();
        TLongObjectHashMap<SyzygyPair<Term, Poly>> ijPairs = new TLongObjectHashMap<>();

        for (int i = 0; i < groebner.size() - 1; i++)
            for (int j = i + 1; j < groebner.size(); ++j) {
                SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(i, j, groebner);
                sPairs.computeIfAbsent(sPair.degree(), __ -> new TreeSet<>(selectionStrategy)).add(sPair);
                ijPairs.put(pack(i, j), sPair);
            }

        while (!sPairs.isEmpty()) {
            // pick up (and remove) all pairs with the smallest degree (normal selection strategy)
            TreeSet<SyzygyPair<Term, Poly>> subset = sPairs.pollFirstEntry().getValue();

            // ! todo if subset.size == 1 ! ? mb just Buchberger ?

            while (!sPairs.isEmpty() && subset.size() < F4_MIN_SELECTION_SIZE)
                subset.addAll(sPairs.pollFirstEntry().getValue());


            subset.forEach(s -> ijPairs.remove(s.hash()));

            System.out.println(subset.size());


            // H: the list of polynomials to reduce
            List<Poly> hPolynomials = new ArrayList<>();
            // all possible monomials in H
            TreeSet<DegreeVector> hMonomials = new TreeSet<>(monomialOrder);
            // monomials that were annihilated
            TreeSet<DegreeVector> hAnnihilated = new TreeSet<>(monomialOrder);

            // populate L with initial values
            for (SyzygyPair<Term, Poly> syzygy : subset) {
                Poly
                        fij = syzygy.fi.clone().multiply(mAlgebra.divideExact(syzygy.syzygyGamma, syzygy.fi.lt())),
                        fji = syzygy.fj.clone().multiply(mAlgebra.divideExact(syzygy.syzygyGamma, syzygy.fj.lt()));

                assert fij.lt().equals(fji.lt());

                // add to set H
                hPolynomials.add(fij);
                hPolynomials.add(fji);

                // store all monomials that we have in H
                hMonomials.addAll(fij.terms.keySet());
                hMonomials.addAll(fji.terms.keySet());

                // lts will be annihilated
                hAnnihilated.add(fij.lt());
            }

            // the diff = Mon(H) / annihilated
            TreeSet<DegreeVector> diff = new TreeSet<>(monomialOrder);
            diff.addAll(hMonomials);
            diff.removeAll(hAnnihilated);

            while (!diff.isEmpty()) {
                // pick the "highest" term from diff
                DegreeVector dv = diff.pollLast();
                hAnnihilated.add(dv);

                // todo <- selection possible
                Optional<Poly> divisorOpt = groebner.stream().filter(g -> dv.dvDivisibleBy(g.lt())).findAny();
                if (divisorOpt.isPresent()) {
                    Poly divisor = divisorOpt.get();
                    Poly newH = divisor.clone().multiplyByDegreeVector(dv.dvDivideExact(divisor.lt()));

                    // append newH to H-polynomials
                    hPolynomials.add(newH);

                    // update monomials set
                    TreeSet<DegreeVector> newMonomials = new TreeSet<>(monomialOrder);
                    newMonomials.addAll(newH.terms.keySet());
                    hMonomials.addAll(newMonomials);

                    // update diff
                    newMonomials.removeAll(hAnnihilated);
                    diff.addAll(newMonomials);
                }
            }

            // all monomials occurring in H
            // columns in decreasing order
            DegreeVector[] hMonomialsArray = hMonomials.descendingSet().toArray(new DegreeVector[hMonomials.size()]);

            // reduce all polys with linear algebra
            List<Poly> nPlus = nPlus(hPolynomials, hMonomialsArray);
            for (Poly g : nPlus) {
                int gSize = groebner.size();
                groebner.add(g);
                for (int i = 0; i < gSize; i++) {
                    SyzygyPair<Term, Poly> sPair = new SyzygyPair<>(i, gSize, groebner);
                    sPairs.computeIfAbsent(sPair.degree(), __ -> new TreeSet<>(selectionStrategy)).add(sPair);
                    ijPairs.put(pack(i, gSize), sPair);
                }
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

    @SuppressWarnings("unchecked")
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    List<Poly> nPlus(List<Poly> hPolynomials,
                     DegreeVector[] hMonomials) {
        if (hPolynomials.get(0) instanceof MultivariatePolynomialZp64)
            return (List<Poly>) nPlusZp64((List) hPolynomials, hMonomials);
        else
            return (List<Poly>) nPlusE((List) hPolynomials, hMonomials);
    }

    /** compute N+ polynomials */
    private static List<MultivariatePolynomialZp64> nPlusZp64(List<MultivariatePolynomialZp64> hPolynomials,
                                                              DegreeVector[] hMonomials) {
        MultivariatePolynomialZp64 factory = hPolynomials.get(0);

        // leading monomials of H-polynomials
        TreeSet<DegreeVector> hLeadMonomials = new TreeSet<>(factory.ordering);
        hPolynomials.forEach(p -> hLeadMonomials.add(p.lt()));

        int nRows = hPolynomials.size();
        int nColumns = hMonomials.length;
        long[][] matrix = new long[nRows][nColumns];
        for (int iCol = 0; iCol < nColumns; ++iCol) {
            DegreeVector dv = hMonomials[iCol];
            for (int iRow = 0; iRow < nRows; ++iRow) {
                MultivariatePolynomialZp64 poly = hPolynomials.get(iRow);
                MonomialZp64 term = poly.terms.get(dv);
                long value = 0;
                if (term != null)
                    value = term.coefficient;
                matrix[iRow][iCol] = value;
            }
        }

        // row reduced echelon form
        LinearSolver.rowEchelonForm(factory.ring, matrix, true);

        // form N+
        List<MultivariatePolynomialZp64> nPlus = new ArrayList<>();
        for (int iRow = 0; iRow < nRows; ++iRow) {
            // todo don't calculate candidate first
            MultivariatePolynomialZp64 candidate = factory.createZero();
            // find lt
            for (int iCol = iRow; iCol < nColumns; ++iCol) {
                long v = matrix[iRow][iCol];
                if (v != 0)
                    candidate.add(new MonomialZp64(hMonomials[iCol], v));
            }

            if (candidate.isZero())
                continue;

            // check whether lt is not reducible by initial lts
            MonomialZp64 lt = candidate.lt();
            if (hLeadMonomials.stream().anyMatch(lt::dvDivisibleBy))
                continue;

            nPlus.add(candidate);
        }

        return nPlus;
    }

    private static void assertRREF(long[][] m) {
        int nRows = m.length;
        int nCols = m[0].length;
        for (int col = 0; col < Math.min(nRows, nCols); col++) {
            if (m[col][col] == 0)
                continue;
            for (int row = 0; row < nRows; row++) {
                if (col != row)
                    assert m[row][col] == 0 : row + " " + col + "\n" + prettyMatrix(m);
            }
            assert m[col][col] == 0 || m[col][col] == 1 : prettyMatrix(m);
        }
    }

    /** compute N+ polynomials */
    private static <E> List<MultivariatePolynomial<E>> nPlusE(List<MultivariatePolynomial<E>> hPolynomials,
                                                              DegreeVector[] hMonomials) {
        DegreeVector[] hLeadMonomials = null;

        MultivariatePolynomial<E> factory = hPolynomials.get(0);
        Ring<E> ring = factory.ring;

        int nRows = hPolynomials.size();
        int nColumns = hMonomials.length;
        E[][] matrix = ring.createArray2d(nRows, nColumns);
        for (int iCol = 0; iCol < nColumns; ++iCol) {
            DegreeVector dv = hMonomials[iCol];
            for (int iRow = 0; iRow < nRows; ++iRow) {
                MultivariatePolynomial<E> poly = hPolynomials.get(iRow);
                Monomial<E> term = poly.terms.get(dv);
                E value;
                if (term != null)
                    value = term.coefficient;
                else
                    value = ring.getZero();
                matrix[iRow][iCol] = value;
            }
        }

        //System.out.println(prettyMatrix( matrix));
        // row reduced echelon form
        LinearSolver.rowEchelonForm(ring, matrix, true);

        //System.out.println(prettyMatrix( matrix));

        // form N+
        List<MultivariatePolynomial<E>> nPlus = new ArrayList<>();
        for (int iRow = 0; iRow < nRows; ++iRow) {
            // todo don't calculate candidate first
            MultivariatePolynomial<E> candidate = factory.createZero();
            // find lt
            for (int iCol = iRow; iCol < nColumns; ++iCol) {
                E v = matrix[iRow][iCol];
                if (!ring.isZero(v))
                    candidate.add(new Monomial<>(hMonomials[iCol], v));
            }

            if (candidate.isZero())
                continue;

            // check whether lt is not reducible by initial lts
            Monomial<E> lt = candidate.lt();
            if (Arrays.stream(hLeadMonomials).anyMatch(lt::dvDivisibleBy))
                continue;

            nPlus.add(candidate);
        }

        return nPlus;
    }

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
