package cc.r2.core.polynomial;

import cc.r2.core.number.Ring;
import cc.r2.core.number.RingElement;

import java.util.Arrays;

import static cc.r2.core.polynomial.Polynomials.newArray;

public class UnivariatePolynomial<R extends RingElement<R>>
        implements RingElement<UnivariatePolynomial<R>> {
    final Ring<R> ring;
    final R[] coefficients;
    final int internalDegree;

    public UnivariatePolynomial(Ring<R> ring, R c) {
        this.ring = ring;

        this.coefficients = newArray(ring, 1);
        this.coefficients[0] = c;
        this.internalDegree = coefficients.length;
    }

    public UnivariatePolynomial(Ring<R> ring, R[] coefficients) {
        this.ring = ring;

        //truncate zero coeffs
        int p = coefficients.length - 1;
        while (p >= 1 && coefficients[p].isZero()) --p;
        if (p != coefficients.length - 1)
            coefficients = Arrays.copyOf(coefficients, p + 1);

        this.coefficients = coefficients;
        this.internalDegree = coefficients.length;
    }

    public R lc() {
        return coefficients[coefficients.length - 1];
    }

    public int degree() {
        return internalDegree - 1;
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> add(UnivariatePolynomial<R> oth) {
        if (oth.isZero())
            return this;
        if (this.internalDegree < oth.internalDegree)
            return oth.add(this);

        final R[] newData = newArray(ring, internalDegree);
        int i = 0;
        for (; i < oth.internalDegree; i++)
            newData[i] = coefficients[i].add(oth.coefficients[i]);
        for (; i < internalDegree; i++)
            newData[i] = coefficients[i];

        return new UnivariatePolynomial<>(ring, newData);
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> subtract(UnivariatePolynomial<R> oth) {
        if (oth.isZero())
            return this;
        if (internalDegree < oth.internalDegree)
            return oth.subtract(this).negate();

        final R[] newData = newArray(ring, internalDegree);
        int i = 0;
        for (; i < oth.internalDegree; i++)
            newData[i] = coefficients[i].subtract(oth.coefficients[i]);
        for (; i < internalDegree; i++)
            newData[i] = coefficients[i];

        return new UnivariatePolynomial<>(ring, newData);
    }

    @Override
    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> negate() {
        final R[] newData = newArray(ring, internalDegree);
        for (int i = 0; i < coefficients.length; i++)
            newData[i] = coefficients[i].negate();
        return new UnivariatePolynomial<R>(ring, newData);
    }

    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> multiplyByFactor(R oth) {
        if (oth.isZero())
            return getZero();
        if (isZero())
            return getZero();
        if (oth.isOne())
            return this;
        final R[] newData = newArray(ring, internalDegree);
        for (int i = 0; i < internalDegree; i++)
            newData[i] = coefficients[i].multiply(oth);
        return new UnivariatePolynomial<>(ring, newData);
    }

    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> multiplyByMonomial(R monomialCoefficient, int monomialDegree) {
        if (monomialDegree == 0)
            return multiplyByFactor(monomialCoefficient);
        if (monomialCoefficient.isZero())
            return getZero();

        final R[] newData = newArray(ring, internalDegree + monomialDegree);
        Arrays.fill(newData, 0, monomialDegree, ring.getZero());
        for (int i = monomialDegree; i < monomialDegree + internalDegree; i++)
            newData[i] = coefficients[i - monomialDegree].multiply(monomialCoefficient);
        return new UnivariatePolynomial<R>(ring, newData);
    }

    @SuppressWarnings("unchecked")
    public UnivariatePolynomial<R> multiply(UnivariatePolynomial<R> oth) {
        if (oth.isZero() || isZero())
            return getZero();
        if (oth.isOne())
            return this;
        if (isOne())
            return oth;

        final R[] newData = newArray(ring, this.internalDegree + oth.internalDegree);
        Arrays.fill(newData, ring.getZero());

        for (int i = 0; i < this.internalDegree; i++)
            for (int j = 0; j < oth.internalDegree; j++)
                newData[i + j] = newData[i + j].add(coefficients[i].multiply(oth.coefficients[j]));

        return new UnivariatePolynomial<>(ring, newData);
    }

    @Override
    public boolean isZero() {
        return internalDegree == 1 && coefficients[0].isZero();
    }

    @Override
    public boolean isOne() {
        return internalDegree == 1 && coefficients[0].isOne();
    }

    @Override
    public UnivariatePolynomial<R> getZero() {
        return new UnivariatePolynomial<>(ring, ring.getZero());
    }

    @Override
    public UnivariatePolynomial<R> getOne() {
        return new UnivariatePolynomial<>(ring, ring.getOne());
    }

    @Override
    public Ring<UnivariatePolynomial<R>> getRing() {
        return null;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        UnivariatePolynomial<?> that = (UnivariatePolynomial<?>) o;
        if (internalDegree != that.internalDegree) return false;
        for (int i = 0; i < internalDegree; i++)
            if (!coefficients[i].equals(that.coefficients[i]))
                return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = Arrays.hashCode(coefficients);
        result = 31 * result + internalDegree;
        return result;
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
