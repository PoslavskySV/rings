package cc.r2.core.poly;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import cc.r2.core.poly.univar.*;
import org.apache.commons.math3.random.RandomGenerator;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class FiniteField implements Domain<lUnivariatePolynomialZp> {
    public static final FiniteField GF27 = new FiniteField(lUnivariatePolynomialZ.create(-1, -1, 0, 1).modulus(3));
    public static final FiniteField GF17p5 = new FiniteField(lUnivariatePolynomialZ.create(11, 11, 0, 3, 9, 9).modulus(17).monic());

    public final lUnivariatePolynomialZp irreducible;
    private final DivisionWithRemainder.InverseModMonomial<lUnivariatePolynomialZp> inverseMod;
    private final BigInteger cardinality;

    public FiniteField(lUnivariatePolynomialZp irreducible) {
        this.irreducible = irreducible;
        this.inverseMod = DivisionWithRemainder.fastDivisionPreConditioning(irreducible);
        this.cardinality = BigIntegerArithmetics.pow(irreducible.domain.modulus, irreducible.degree());
    }

    public long modulus() {return irreducible.domain.modulus;}

    @Override
    public boolean isField() {
        return true;
    }

    @Override
    public BigInteger cardinality() {
        return cardinality;
    }

    @Override
    public BigInteger characteristics() {return BigInteger.valueOf(modulus());}

    @Override
    public lUnivariatePolynomialZp add(lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
        return PolynomialArithmetics.polyAddMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public lUnivariatePolynomialZp subtract(lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
        return PolynomialArithmetics.polySubtractMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public lUnivariatePolynomialZp multiply(lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
        return PolynomialArithmetics.polyMultiplyMod(a, b, irreducible, inverseMod, true);
    }

    @Override
    public lUnivariatePolynomialZp negate(lUnivariatePolynomialZp val) {
        return PolynomialArithmetics.polyNegateMod(val, irreducible, inverseMod, true);
    }

    @Override
    public int signum(lUnivariatePolynomialZp a) {
        return Long.signum(a.lc());
    }

    @Override
    public lUnivariatePolynomialZp[] divideAndRemainder(lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
        return new lUnivariatePolynomialZp[]{multiply(a, reciprocal(b)), getZero()};
    }

    @Override
    public lUnivariatePolynomialZp reciprocal(lUnivariatePolynomialZp a) {
        lUnivariatePolynomialZp p = irreducible;
        lUnivariatePolynomialZp
                t = irreducible.createZero(),
                newt = irreducible.createOne(),
                r = p.clone(),
                newr = a.clone();
        lUnivariatePolynomialZp tmp;
        while (!newr.isZero()) {
            lUnivariatePolynomialZp quotient = DivisionWithRemainder.quotient(r, newr, true);

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
    public lUnivariatePolynomialZp gcd(lUnivariatePolynomialZp a, lUnivariatePolynomialZp b) {
        return UnivariateGCD.PolynomialGCD(a, b);
    }

    @Override
    public lUnivariatePolynomialZp getZero() {
        return irreducible.createZero();
    }

    @Override
    public lUnivariatePolynomialZp getOne() {
        return irreducible.createOne();
    }

    @Override
    public boolean isZero(lUnivariatePolynomialZp poly) {
        return poly.isZero();
    }

    @Override
    public boolean isOne(lUnivariatePolynomialZp poly) {
        return poly.isOne();
    }

    @Override
    public lUnivariatePolynomialZp valueOf(long val) {
        return irreducible.createConstant(val);
    }

    @Override
    public lUnivariatePolynomialZp valueOf(lUnivariatePolynomialZp val) {
        return PolynomialArithmetics.polyMod(val, irreducible, inverseMod, true);
    }

    public lUnivariatePolynomialZp valueOf(lUnivariatePolynomialZ val) {
        return PolynomialArithmetics.polyMod(val.modulus(irreducible.domain), irreducible, inverseMod, true);
    }

    @Override
    public lUnivariatePolynomialZp[] createArray(int length) {
        return new lUnivariatePolynomialZp[length];
    }

    @Override
    public lUnivariatePolynomialZp[][] createArray2d(int length) {
        return new lUnivariatePolynomialZp[length][];
    }

    @Override
    public lUnivariatePolynomialZp[][] createArray2d(int m, int n) {
        return new lUnivariatePolynomialZp[m][n];
    }

    @Override
    public int compare(lUnivariatePolynomialZp o1, lUnivariatePolynomialZp o2) {
        return o1.compareTo(o2);
    }

    @Override
    public lUnivariatePolynomialZp randomElement(RandomGenerator rnd) {
        return valueOf(RandomPolynomials.randomPoly(2 * irreducible.degree(), rnd));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        FiniteField that = (FiniteField) o;

        if (!irreducible.equals(that.irreducible)) return false;
        if (!inverseMod.equals(that.inverseMod)) return false;
        return cardinality.equals(that.cardinality);
    }

    @Override
    public int hashCode() {
        int result = irreducible.hashCode();
        result = 31 * result + inverseMod.hashCode();
        result = 31 * result + cardinality.hashCode();
        return result;
    }
}
