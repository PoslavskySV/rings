package cc.redberry.rings.poly.multivar;


import cc.redberry.rings.Rings;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;
import cc.redberry.rings.util.ArraysUtil;

import java.util.Arrays;

/**
 * @author Stanislav Poslavsky
 * @since 2.2
 */
public final class MultivariateConversions {
    private MultivariateConversions() {}

    /**
     * Given poly in R[x1,x2,...,xN] converts to poly in R[variables][other_variables]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    MultivariatePolynomial<Poly> split(Poly poly, int... variables) {
        return poly.asOverMultivariateEliminate(variables);
    }

    /**
     * Given poly in R[variables][other_variables] converts it to poly in R[x1,x2,...,xN]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    Poly mergeSplit(MultivariatePolynomial<Poly> poly, int... variables) {
        variables = variables.clone();
        Arrays.sort(variables);
        int[] mainVariables = ArraysUtil.intSetDifference(
                ArraysUtil.sequence(0, poly.nVariables + variables.length),
                variables);
        if (poly.cc() instanceof MultivariatePolynomial)
            return (Poly) MultivariatePolynomial.asNormalMultivariate((MultivariatePolynomial) poly, variables, mainVariables);
        else
            return (Poly) MultivariatePolynomialZp64.asNormalMultivariate((MultivariatePolynomial) poly, variables, mainVariables);
    }

    /**
     * Given poly in R[x1,x2,...,xN] converts to poly in R[variables][other_variables]
     */
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    IPolynomialRing<MultivariatePolynomial<Poly>> split(IPolynomialRing<Poly> ring, int... variables) {
        return Rings.MultivariateRing(split(ring.factory(), variables));
    }

    /**
     * Given poly in R[x1,x2,...,xN] converts to poly in R[variables][other_variables]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    IPolynomialRing<Poly> mergeSplit(IPolynomialRing<MultivariatePolynomial<Poly>> ring, int... variables) {
        return (IPolynomialRing<Poly>) Rings.MultivariateRing((AMultivariatePolynomial) mergeSplit(ring.factory(), variables));
    }

    /**
     * Given poly in R[x1,x2,...,xN] converts to poly in R[other_variables][variable]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    UnivariatePolynomial<Poly> splitUnivariate(Poly poly, int variable) {
        return poly.asUnivariateEliminate(variable);
    }

    /**
     * Given poly in R[variables][other_variables] converts it to poly in R[x1,x2,...,xN]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    Poly mergeUnivariateSplit(UnivariatePolynomial<Poly> poly, int variable) {
        if (poly.cc() instanceof MultivariatePolynomial)
            return (Poly) MultivariatePolynomial.asMultivariate((UnivariatePolynomial) poly, variable, true);
        else
            return (Poly) MultivariatePolynomialZp64.asMultivariate((UnivariatePolynomial) poly, variable, true);
    }

    /**
     * Given poly in R[x1,x2,...,xN] converts to poly in R[other_variables][variable]
     */
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    IPolynomialRing<UnivariatePolynomial<Poly>> splitUnivariate(IPolynomialRing<Poly> ring, int variable) {
        return Rings.UnivariateRing(splitUnivariate(ring.factory(), variable));
    }

    /**
     * Given poly in R[variables][other_variables] converts it to poly in R[x1,x2,...,xN]
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends AMultivariatePolynomial<?, Poly>>
    IPolynomialRing<Poly> mergeUnivariateSplit(IPolynomialRing<UnivariatePolynomial<Poly>> ring, int variable) {
        return (IPolynomialRing<Poly>) Rings.MultivariateRing((AMultivariatePolynomial) mergeUnivariateSplit(ring.factory(), variable));
    }
}
