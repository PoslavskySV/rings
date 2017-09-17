package cc.r2.core.poly;

/**
 * Polynomial domain.
 *
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public interface PolynomialDomain<Poly extends IPolynomial<Poly>> extends Domain<Poly>, WithVariables {
    /**
     * Number of polynomial variables
     */
    int nVariables();

    /**
     * Factory polynomial
     */
    Poly factory();

    /**
     * Parse string into polynomial
     *
     * @param string    string
     * @param variables names of variables
     */
    Poly parse(String string, String[] variables);

    /**
     * Creates monomial corresponding to specified variable
     */
    Poly variable(int variable);

    /**
     * String representation of this domain.
     *
     * @param coefficientDomain coefficient domain
     * @param vars              variable names
     */
    String toString(String coefficientDomain, String[] vars);
}
