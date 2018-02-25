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
    public static final Comparator<DegreeVector> LEX = (Comparator<DegreeVector> & Serializable)
            (DegreeVector a, DegreeVector b) -> {
                for (int i = 0; i < a.exponents.length; ++i) {
                    int c = Integer.compare(a.exponents[i], b.exponents[i]);
                    if (c != 0)
                        return c;
                }
                return 0;
            };
    /**
     * Graded lexicographic monomial order.
     */
    public static final Comparator<DegreeVector> GRLEX = (Comparator<DegreeVector> & Serializable)
            (DegreeVector a, DegreeVector b) -> {
                int c = Integer.compare(a.totalDegree, b.totalDegree);
                return c != 0 ? c : LEX.compare(a, b);
            };
    /**
     * Antilexicographic monomial order.
     */
    public static final Comparator<DegreeVector> ALEX = (Comparator<DegreeVector> & Serializable)
            (DegreeVector a, DegreeVector b) -> LEX.compare(b, a);
    /**
     * Graded reverse lexicographic monomial order
     */
    public static final Comparator<DegreeVector> GREVLEX = (Comparator<DegreeVector> & Serializable)
            (Comparator<DegreeVector> & Serializable) (DegreeVector a, DegreeVector b) -> {
                int c = Integer.compare(a.totalDegree, b.totalDegree);
                if (c != 0)
                    return c;
                for (int i = a.exponents.length - 1; i >= 0; --i) {
                    c = Integer.compare(b.exponents[i], a.exponents[i]);
                    if (c != 0)
                        return c;
                }
                return 0;
            };

    /** Default monomial order (GREVLEX) */
    public static final Comparator<DegreeVector> DEFAULT = GREVLEX;

    /**
     * Block product of orderings
     */
    public static Comparator<DegreeVector> product(Comparator<DegreeVector> orderings[], int[] nVariables) {
        return new ProductOrdering(orderings, nVariables);
    }

    /**
     * Block product of orderings
     */
    @SuppressWarnings("unchecked")
    public static Comparator<DegreeVector> product(Comparator<DegreeVector> a, int anVariables, Comparator<DegreeVector> b, int bnVariable) {
        return new ProductOrdering(new Comparator[]{a, b}, new int[]{anVariables, bnVariable});
    }

    static final class ProductOrdering implements Comparator<DegreeVector>, Serializable {
        final Comparator<DegreeVector> orderings[];
        final int[] nVariables;

        ProductOrdering(Comparator<DegreeVector>[] orderings, int[] nVariables) {
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

            ProductOrdering that = (ProductOrdering) o;

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
    }
}
