package cc.redberry.rings.poly.multivar;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Common monomial orderings.
 *
 * @since 1.0
 */
public final class MonomialOrder {
    private MonomialOrder() { }

    /**
     * Lexicographic monomial order.
     */
    public static final Comparator<DegreeVector> LEX = Lex.instance;

    private static final class Lex implements Comparator<DegreeVector>, Serializable {
        private static final Lex instance = new Lex();

        private Lex() {}

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            for (int i = 0; i < a.exponents.length; ++i) {
                int c = Integer.compare(a.exponents[i], b.exponents[i]);
                if (c != 0)
                    return c;
            }
            return 0;
        }

        protected Object readResolve() { return instance; }
    }

    /**
     * Graded lexicographic monomial order.
     */
    public static final Comparator<DegreeVector> GRLEX = Grlex.instance;

    private static final class Grlex implements Comparator<DegreeVector>, Serializable {
        private static final Grlex instance = new Grlex();

        private Grlex() {}

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            int c = Integer.compare(a.totalDegree, b.totalDegree);
            return c != 0 ? c : LEX.compare(a, b);
        }

        protected Object readResolve() { return instance; }
    }


    /**
     * Antilexicographic monomial order.
     */
    public static final Comparator<DegreeVector> ALEX = Alex.instance;

    private static final class Alex implements Comparator<DegreeVector>, Serializable {
        private static final Alex instance = new Alex();

        private Alex() {}

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            return LEX.compare(b, a);
        }

        protected Object readResolve() { return instance; }
    }

    /**
     * Graded reverse lexicographic monomial order
     */
    public static final Comparator<DegreeVector> GREVLEX = Grevlex.instance;

    private static final class Grevlex implements Comparator<DegreeVector>, Serializable {
        private static final Grevlex instance = new Grevlex();

        private Grevlex() {}

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            int c = Integer.compare(a.totalDegree, b.totalDegree);
            if (c != 0)
                return c;
            for (int i = a.exponents.length - 1; i >= 0; --i) {
                c = Integer.compare(b.exponents[i], a.exponents[i]);
                if (c != 0)
                    return c;
            }
            return 0;
        }

        protected Object readResolve() { return instance; }
    }

    /** Default monomial order (GREVLEX) */
    public static final Comparator<DegreeVector> DEFAULT = parse(System.getProperty("defaultMonomialOrder", "grevlex").toLowerCase());

    static Comparator<DegreeVector> parse(String string) {
        switch (string.toLowerCase()) {
            case "lex": return LEX;
            case "grlex": return GRLEX;
            case "grevlex": return GREVLEX;
            case "alex": return ALEX;
            default: throw new RuntimeException("unknown: " + string);
        }
    }

    /**
     * Block product of orderings
     */
    public static Comparator<DegreeVector> product(Comparator<DegreeVector> orderings[], int[] nVariables) {
        return new ProductOrder(orderings, nVariables);
    }

    /**
     * Block product of orderings
     */
    @SuppressWarnings("unchecked")
    public static Comparator<DegreeVector> product(Comparator<DegreeVector> a, int anVariables, Comparator<DegreeVector> b, int bnVariable) {
        return new ProductOrder(new Comparator[]{a, b}, new int[]{anVariables, bnVariable});
    }

    /** whether monomial order is graded */
    public static boolean isGradedOrder(Comparator<DegreeVector> monomialOrder) {
        return monomialOrder == GREVLEX
                || monomialOrder == GRLEX
                || monomialOrder instanceof GrevLexWithPermutation
                || (monomialOrder instanceof EliminationOrder && isGradedOrder(((EliminationOrder) monomialOrder).baseOrder));
    }

    static final class ProductOrder implements Comparator<DegreeVector>, Serializable {
        final Comparator<DegreeVector> orderings[];
        final int[] nVariables;

        ProductOrder(Comparator<DegreeVector>[] orderings, int[] nVariables) {
            this.orderings = orderings;
            this.nVariables = nVariables;
        }

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            int prev = 0;
            for (int i = 0; i < nVariables.length; i++) {
                // for each block
                DegreeVector
                        aBlock = a.dvRange(prev, prev + nVariables[i]),
                        bBlock = b.dvRange(prev, prev + nVariables[i]);

                int c = orderings[i].compare(aBlock, bBlock);
                if (c != 0)
                    return c;

                prev += nVariables[i];
            }
            return 0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            ProductOrder that = (ProductOrder) o;

            // Probably incorrect - comparing Object[] arrays with Arrays.equals
            if (!Arrays.equals(orderings, that.orderings)) return false;
            return Arrays.equals(nVariables, that.nVariables);
        }

        @Override
        public int hashCode() {
            int result = Arrays.hashCode(orderings);
            result = 31 * result + Arrays.hashCode(nVariables);
            return result;
        }
    }

    public static final class GrevLexWithPermutation implements Comparator<DegreeVector>, Serializable {
        final int[] permutation;

        GrevLexWithPermutation(int[] permutation) {
            this.permutation = permutation;
        }

        @Override
        public int compare(DegreeVector a, DegreeVector b) {
            int c = Integer.compare(a.totalDegree, b.totalDegree);
            if (c != 0)
                return c;
            for (int i = a.exponents.length - 1; i >= 0; --i) {
                c = Integer.compare(b.exponents[permutation[i]], a.exponents[permutation[i]]);
                if (c != 0)
                    return c;
            }
            return 0;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            GrevLexWithPermutation that = (GrevLexWithPermutation) o;
            return Arrays.equals(permutation, that.permutation);
        }

        @Override
        public int hashCode() {
            return Arrays.hashCode(permutation);
        }
    }

    public static final class EliminationOrder implements Comparator<DegreeVector>, Serializable {
        final Comparator<DegreeVector> baseOrder;
        final int variable;

        public EliminationOrder(Comparator<DegreeVector> baseOrder, int variable) {
            this.baseOrder = baseOrder;
            this.variable = variable;
        }

        @Override
        public int compare(DegreeVector o1, DegreeVector o2) {
            int c = Integer.compare(o1.exponents[variable], o2.exponents[variable]);
            if (c != 0)
                return c;
            return baseOrder.compare(o1, o2);
        }
    }
}