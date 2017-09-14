package cc.r2.core.poly;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface PolynomialDomain<Poly extends IPolynomial<Poly>> extends Domain<Poly>, WithVariables {
    /**
     * Parse string into polynomial
     *
     * @param string    string
     * @param variables names of variables
     * @return domain element
     */
    Poly parse(String string, String[] variables);

    /**
     * String representation of this domain
     *
     * @param coefficientDomain coefficient domain
     * @param vars              variable names
     */
    String toString(String coefficientDomain, String[] vars);
}
