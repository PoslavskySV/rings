package cc.r2.core.poly.multivar;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Different monomial orderings.
 *
 * @author Stanislav Poslavsky
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
}
