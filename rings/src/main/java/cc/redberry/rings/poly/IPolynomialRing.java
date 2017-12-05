package cc.redberry.rings.poly;

import cc.redberry.rings.Ring;
import cc.redberry.rings.WithVariables;

/**
 * Polynomial ring.
 *
 * @since 1.0
 */
public interface IPolynomialRing<Poly extends IPolynomial<Poly>> extends Ring<Poly>, WithVariables {
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
     * String representation of this ring.
     *
     * @param coefficientDomain coefficient ring
     * @param vars              variable names
     */
    String toString(String coefficientDomain, String[] vars);
}
