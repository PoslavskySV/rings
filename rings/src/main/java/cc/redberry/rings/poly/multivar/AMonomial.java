package cc.redberry.rings.poly.multivar;

import cc.redberry.rings.util.ArraysUtil;

/**
 * Abstract monomial (degree vector + coefficient). The parent class for {@link MonomialZp64} and {@link Monomial}.
 * Instances are immutable. Algebraic operations on monomials (multiplication and division) are specified in {@link
 * IMonomialAlgebra}.
 *
 * @see DegreeVector
 * @see IMonomialAlgebra
 * @since 2.3
 */
public abstract class AMonomial<Term extends AMonomial<Term>> extends DegreeVector {
    /**
     * @param exponents   exponents
     * @param totalDegree total degree (sum of exponents)
     */
    protected AMonomial(int[] exponents, int totalDegree) {
        super(exponents, totalDegree);
    }

    /**
     * @param exponents exponents
     */
    protected AMonomial(int[] exponents) {
        this(exponents, ArraysUtil.sum(exponents));
    }

    /**
     * @param degreeVector DegreeVector
     */
    protected AMonomial(DegreeVector degreeVector) {
        this(degreeVector.exponents, degreeVector.totalDegree);
    }

    /** Sets coefficient of this with coefficient of oth */
    public abstract Term setCoefficientFrom(Term oth);

    /** Sets the degree vector */
    public abstract Term setDegreeVector(DegreeVector oth);

    /** Sets the degree vector */
    public abstract Term setDegreeVector(int[] exponents, int totalDegree);

    /** Sets the degree vector */
    public final Term setDegreeVector(int[] exponents) { return setDegreeVector(exponents, ArraysUtil.sum(exponents));}

    /** Multiplies this by oth */
    public final Term multiply(DegreeVector oth) {return setDegreeVector(dvMultiply(oth));}

    /** Multiplies this by oth */
    public final Term multiply(int[] oth) {return setDegreeVector(dvMultiply(oth));}

    /** Gives quotient {@code this / oth } or null if exact division is not possible (e.g. a^2*b^3 / a^3*b^5) */
    public final Term divideOrNull(DegreeVector divider) {return setDegreeVector(dvDivideOrNull(divider));}

    /** Gives quotient {@code this / oth } or null if exact division is not possible (e.g. a^2*b^3 / a^3*b^5) */
    public final Term divideOrNull(int[] divider) {return setDegreeVector(dvDivideOrNull(divider));}

    /** Joins new variable (with zero exponent) to degree vector */
    public final Term joinNewVariable() {return setDegreeVector(dvJoinNewVariable());}

    /** Joins new variables (with zero exponents) to degree vector */
    public final Term joinNewVariables(int n) {return setDegreeVector(dvJoinNewVariables(n));}

    /** internal API */
    public final Term joinNewVariables(int newNVariables, int[] mapping) {return setDegreeVector(dvJoinNewVariables(newNVariables, mapping));}

    /** Sets the number of variables */
    public final Term setNVariables(int n) {return setDegreeVector(dvSetNVariables(n));}

    /** Sets exponents of all variables except the specified variable to zero */
    public final Term select(int var) {return setDegreeVector(dvSelect(var));}

    /** Set's exponents of all variables except specified variables to zero */
    public final Term select(int[] variables) {return setDegreeVector(dvSelect(variables));}

    /**
     * Selects range from this
     *
     * @param from from inclusive
     * @param to   to exclusive
     */
    public final Term range(int from, int to) {return setDegreeVector(dvRange(from, to));}

    /** Set exponent of specified {@code var} to zero */
    public final Term setZero(int var) {return setDegreeVector(dvSetZero(var));}

    /** Set exponents of specified variables to zero */
    public final Term setZero(int[] variables) {return setDegreeVector(dvSetZero(variables));}

    /** Drops specified variable (number of variables will be reduced) */
    public final Term without(int variable) {return setDegreeVector(dvWithout(variable));}

    /** Drops specified variables (number of variables will be reduced) */
    public final Term without(int[] variables) {return setDegreeVector(dvWithout(variables));}

    /** Inserts new variable (with zero exponent) */
    public final Term insert(int variable) {return setDegreeVector(dvInsert(variable));}

    /**
     * Set's exponent of specified variable to specified value
     *
     * @param variable the variable
     * @param exponent new exponent
     */
    public final Term set(int variable, int exponent) {return setDegreeVector(dvSet(variable, exponent));}
}
