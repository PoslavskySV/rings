package cc.redberry.rings.poly;

import cc.redberry.rings.IntegersZp;
import cc.redberry.rings.Rational;
import cc.redberry.rings.Ring;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.poly.multivar.MultivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;

/**
 * @since 1.0
 */
public final class Util {
    private Util() {}

    public static void ensureOverFiniteField(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverFiniteField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureOverField(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureOverZ(IPolynomial... polys) {
        for (IPolynomial poly : polys)
            if (!poly.isOverZ())
                throw new IllegalArgumentException("Polynomial over Z is expected, but got " + poly.getClass());
    }

    /**
     * Test whether poly is over Zp with modulus less then 2^63
     */
    public static boolean canConvertToZp64(IPolynomial poly) {
        Ring ring;
        if (poly instanceof UnivariatePolynomial)
            ring = ((UnivariatePolynomial) poly).ring;
        else if (poly instanceof MultivariatePolynomial)
            ring = ((MultivariatePolynomial) poly).ring;
        else
            return false;

        return ring instanceof IntegersZp && ((IntegersZp) ring).modulus.bitLength() < MachineArithmetic.MAX_SUPPORTED_MODULUS_BITS;
    }

    /** Whether coefficient domain is rationals */
    public static <T extends IPolynomial<T>> boolean isOverRationals(T poly) {
        if (poly instanceof UnivariatePolynomial && ((UnivariatePolynomial) poly).ring.getOne() instanceof Rational)
            return true;
        else if (poly instanceof MultivariatePolynomial && ((MultivariatePolynomial) poly).ring.getOne() instanceof Rational)
            return true;
        else
            return false;
    }

    /** Whether coefficient domain is Q */
    public static <T extends IPolynomial<T>> boolean isOverQ(T poly) {
        Object rep;

        if (poly instanceof UnivariatePolynomial)
            rep = ((UnivariatePolynomial) poly).ring.getOne();
        else if (poly instanceof MultivariatePolynomial)
            rep = ((MultivariatePolynomial) poly).ring.getOne();
        else
            return false;

        if (!(rep instanceof Rational))
            return false;

        return ((Rational) rep).numerator() instanceof BigInteger;
    }

    public static final class Tuple2<A, B> {
        public final A _1;
        public final B _2;

        public Tuple2(A _1, B _2) {
            this._1 = _1;
            this._2 = _2;
        }
    }

    /**
     * Brings polynomial with rational coefficients to common denominator
     *
     * @param poly the polynomial
     * @return (reduced poly, common denominator)
     */
    public static <E> Tuple2<UnivariatePolynomial<E>, E> toCommonDenominator(UnivariatePolynomial<Rational<E>> poly) {
        Ring<Rational<E>> field = poly.ring;
        Ring<E> integralRing = field.getOne().ring;
        E denominator = integralRing.getOne();
        for (int i = 0; i <= poly.degree(); i++)
            if (!poly.isZeroAt(i))
                denominator = integralRing.lcm(denominator, poly.get(i).denominator());

        E[] data = integralRing.createArray(poly.degree() + 1);
        for (int i = 0; i <= poly.degree(); i++) {
            Rational<E> cf = poly.get(i).multiply(denominator);
            assert cf.isIntegral();
            data[i] = cf.numerator();
        }
        return new Tuple2<>(UnivariatePolynomial.createUnsafe(integralRing, data), denominator);
    }

    /**
     * Brings polynomial with rational coefficients to common denominator
     *
     * @param poly the polynomial
     * @return (reduced poly, common denominator)
     */
    public static <E> Tuple2<MultivariatePolynomial<E>, E> toCommonDenominator(MultivariatePolynomial<Rational<E>> poly) {
        Ring<Rational<E>> field = poly.ring;
        Ring<E> integralRing = field.getOne().ring;
        E denominator = integralRing.getOne();
        for (Rational<E> cf : poly.coefficients())
            denominator = integralRing.lcm(denominator, cf.denominator());

        final E d = denominator;
        MultivariatePolynomial<E> integral = poly.mapCoefficients(integralRing, cf -> {
            Rational<E> r = cf.multiply(d);
            assert integralRing.isOne(r.denominator());
            return r.numerator();
        });
        return new Tuple2<>(integral, denominator);
    }

    public static <E> UnivariatePolynomial<Rational<E>> asOverRationals(Ring<Rational<E>> field, UnivariatePolynomial<E> poly) {
        return poly.mapCoefficients(field, cf -> new Rational<>(poly.ring, cf));
    }

    public static <E> MultivariatePolynomial<Rational<E>> asOverRationals(Ring<Rational<E>> field, MultivariatePolynomial<E> poly) {
        return poly.mapCoefficients(field, cf -> new Rational<>(poly.ring, cf));
    }

    public static <E> UnivariatePolynomial<Rational<E>> divideOverRationals(Ring<Rational<E>> field, UnivariatePolynomial<E> poly, E denominator) {
        return poly.mapCoefficients(field, cf -> new Rational<>(poly.ring, cf, denominator));
    }

    public static <E> MultivariatePolynomial<Rational<E>> divideOverRationals(Ring<Rational<E>> field, MultivariatePolynomial<E> poly, E denominator) {
        return poly.mapCoefficients(field, cf -> new Rational<>(poly.ring, cf, denominator));
    }
}
