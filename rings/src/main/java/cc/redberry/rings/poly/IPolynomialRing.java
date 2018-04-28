package cc.redberry.rings.poly;

import cc.redberry.rings.Ring;
import cc.redberry.rings.io.Coder;

/**
 * Polynomial ring.
 *
 * @since 1.0
 */
public interface IPolynomialRing<Poly extends IPolynomial<Poly>> extends Ring<Poly> {
    /**
     * Number of polynomial variables
     */
    int nVariables();

    /**
     * Factory polynomial
     */
    Poly factory();

    /**
     * Creates poly representing a single specified variable
     */
    Poly variable(int variable);

    /**
     * Parse poly from string using specified variables representation
     */
    default Poly parse(String string, String... variables) {
        return mkCoder(variables).parse(string);
    }

    /**
     * Simple coder for this ring
     */
    default Coder<Poly, ?, ?> mkCoder(String... variables) {
        return Coder.mkPolynomialCoder(this, variables);
    }
}
