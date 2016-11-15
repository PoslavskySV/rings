package cc.r2.core.polynomial;

import cc.r2.core.number.Ring;
import cc.r2.core.number.RingElement;

import java.util.Arrays;

import static cc.r2.core.polynomial.Polynomials.newArray;

final class MutableUnivariatePolynomial<R extends RingElement<R>> {
    final Ring<R> ring;
    R[] coefficients;
    int internalDegree;

    public MutableUnivariatePolynomial(Ring<R> ring, int degree, R... coefficients) {
        this.ring = ring;
        this.coefficients = coefficients;
        this.internalDegree = degree + 1;
    }

    public R lc() {
        return coefficients[internalDegree - 1];
    }

    public int degree() {
        return internalDegree - 1;
    }

    @SuppressWarnings("unchecked")
    private void ensureCapacity(int newInternalDegree) {
        if (internalDegree < newInternalDegree) {

            R[] newData;
            if (coefficients.length < newInternalDegree)
                newData = newArray(ring, newInternalDegree);
            else
                newData = coefficients;

            Arrays.fill(newData, internalDegree, newInternalDegree, ring.getZero());
            System.arraycopy(coefficients, 0, newData, 0, internalDegree);
            this.coefficients = newData;
            this.internalDegree = newInternalDegree;
        }
    }

    private void updateDegree() {
        int i = internalDegree - 1;
        while (i >= 1 && coefficients[i].isZero()) --i;
        internalDegree = i + 1;
    }


    MutableUnivariatePolynomial<R> add(R[] oth) {
        ensureCapacity(oth.length);
        for (int i = 0; i < oth.length; i++)
            coefficients[i] = coefficients[i].add(oth[i]);
        updateDegree();
        return this;
    }

    MutableUnivariatePolynomial<R> subtract(R[] oth) {
        ensureCapacity(oth.length);
        for (int i = 0; i < oth.length; i++)
            coefficients[i] = coefficients[i].subtract(oth[i]);
        updateDegree();
        return this;
    }

    MutableUnivariatePolynomial<R> subtract(R[] oth, R monomialCoeff, int monomialDegree) {
        ensureCapacity(oth.length + monomialDegree);
        for (int i = monomialDegree; i < oth.length + monomialDegree; i++)
            coefficients[i] = coefficients[i].subtract(monomialCoeff.multiply(oth[i - monomialDegree]));
        updateDegree();
        return this;
    }

    MutableUnivariatePolynomial<R> negate(R[] oth) {
        for (int i = 0; i < internalDegree; i++)
            coefficients[i] = coefficients[i].negate();
        return this;
    }

    UnivariatePolynomial<R> toPoly() {
        return new UnivariatePolynomial<>(ring, Arrays.copyOf(coefficients, internalDegree));
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < coefficients.length; i++) {
            if (coefficients[i].isZero())
                continue;
            String c = coefficients[i].toString();
            if (!c.startsWith("-") && sb.length() != 0)
                sb.append("+");
            sb.append(c);
            if (i != 0)
                sb.append("x^").append(i);
        }

        if (sb.length() == 0)
            return "0";
        return sb.toString();
    }
}
