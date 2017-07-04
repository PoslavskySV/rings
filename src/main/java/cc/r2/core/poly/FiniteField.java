package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.poly.univar.*;
import org.apache.commons.math3.random.RandomGenerator;

import java.lang.reflect.Array;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class FiniteField<Poly extends IUnivariatePolynomial<Poly>> extends ADomain<Poly> {
    /** GF(3^3) */
    public static final FiniteField<lUnivariatePolynomialZp> GF27 = new FiniteField<>(lUnivariatePolynomialZ.create(-1, -1, 0, 1).modulus(3));
    /** GF(17^5) */
    public static final FiniteField<lUnivariatePolynomialZp> GF17p5 = new FiniteField<>(lUnivariatePolynomialZ.create(11, 11, 0, 3, 9, 9).modulus(17).monic());

    /** Irreducible polynomial that generates this finite field */
    public final Poly irreducible;
    private final DivisionWithRemainder.InverseModMonomial<Poly> inverseMod;
    private final BigInteger cardinality;

    /**
     * Constructs finite field from the specified irreducible polynomial. NOTE: irreducibility test for the input
     * polynomial is not performed here, use {@link cc.r2.core.poly.univar.IrreduciblePolynomials#irreducibleQ(IUnivariatePolynomial)}
     * to test irreducibility
     *
     * @param irreducible irreducible polynomial
     */
    public FiniteField(Poly irreducible) {
        if (!irreducible.isOverField())
            throw new IllegalArgumentException();
        this.irreducible = irreducible;
        this.inverseMod = DivisionWithRemainder.fastDivisionPreConditioning(irreducible);
        this.cardinality = BigIntegerArithmetics.pow(irreducible.coefficientDomainCardinality(), irreducible.degree());
    }

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public BigInteger cardinality() {
        return cardinality;
    }

    @Override
    public BigInteger characteristics() {
        return irreducible.coefficientDomainCardinality();
    }

    @Override
    public Poly add(Poly a, Poly b) {
        return PolynomialArithmetics.polyAddMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly subtract(Poly a, Poly b) {
        return PolynomialArithmetics.polySubtractMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly multiply(Poly a, Poly b) {
        return PolynomialArithmetics.polyMultiplyMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public Poly negate(Poly val) {
        return PolynomialArithmetics.polyNegateMod(val, irreducible, inverseMod, true);
    }

    @Override
    public int signum(Poly a) {
        return a.signum();
    }

    @Override
    public Poly[] divideAndRemainder(Poly a, Poly b) {
        return a.arrayNewInstance(multiply(a, reciprocal(b)), getZero());
    }

    @Override
    public Poly remainder(Poly dividend, Poly divider) {
        return getZero();
    }

    @Override
    public Poly reciprocal(Poly a) {
        Poly
                t = getZero(),
                newt = getOne(),
                r = irreducible.clone(),
                newr = a.clone();
        Poly tmp;
        while (!newr.isZero()) {
            Poly quotient = DivisionWithRemainder.quotient(r, newr, true);

            tmp = r;
            r = newr;
            newr = tmp.clone().subtract(quotient.clone().multiply(newr));

            tmp = t;
            t = newt;
            newt = tmp.clone().subtract(quotient.clone().multiply(newt));
        }
        if (r.degree() > 0)
            throw new IllegalArgumentException("Either p is not irreducible or a is a multiple of p");
        return t.divideByLC(r);
    }

    @Override
    public Poly gcd(Poly a, Poly b) {
        return UnivariateGCD.PolynomialGCD(a, b);
    }

    @Override
    public Poly getZero() {
        return irreducible.createZero();
    }

    @Override
    public Poly getOne() {
        return irreducible.createOne();
    }

    @Override
    public boolean isZero(Poly poly) {
        return poly.isZero();
    }

    @Override
    public boolean isOne(Poly poly) {
        return poly.isOne();
    }

    @Override
    public boolean isUnit(Poly poly) {
        return !isZero(poly);
    }

    @Override
    public Poly valueOf(long val) {
        return getOne().multiply(val);
    }

    @Override
    public Poly valueOf(Poly val) {
        return PolynomialArithmetics.polyMod(val, irreducible, inverseMod, true);
    }

    @Override
    public Poly[] createArray(int length) {
        return irreducible.arrayNewInstance(length);
    }

    @Override
    @SuppressWarnings("unchecked")
    public Poly[][] createArray2d(int length) {
        Poly[] array = createArray(0);
        return (Poly[][]) Array.newInstance(array.getClass(), length);
    }

    @Override
    public Poly[][] createArray2d(int m, int n) {
        Poly[][] arr = createArray2d(m);
        for (int i = 0; i < arr.length; i++)
            arr[i] = createArray(n);
        return arr;
    }

    @Override
    public int compare(Poly o1, Poly o2) {
        return o1.compareTo(o2);
    }

    @Override
    public Poly randomElement(RandomGenerator rnd) {
        return valueOf(RandomPolynomials.randomPoly(irreducible, 2 * irreducible.degree(), rnd));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FiniteField<?> that = (FiniteField<?>) o;

        if (!irreducible.equals(that.irreducible)) return false; // todo: need to check this?
        return cardinality.equals(that.cardinality);
    }

    @Override
    public int hashCode() {
        int result = irreducible.hashCode();
        result = 31 * result + cardinality.hashCode();
        return result;
    }
}
