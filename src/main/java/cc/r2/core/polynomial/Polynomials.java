package cc.r2.core.polynomial;

import cc.r2.core.number.EuclideanRingElement;
import cc.r2.core.number.Ring;
import cc.r2.core.number.RingElement;
import gnu.trove.list.array.TIntArrayList;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by poslavsky on 04/11/2016.
 */
public final class Polynomials {
    private Polynomials() {
    }

    @SuppressWarnings("unchecked")
    static <R extends RingElement<R>> R[] newArray(Ring<R> ring, int length) {
        return (R[]) Array.newInstance(ring.getElementType(), length);
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    R content(UnivariatePolynomial<R> poly) {
        if (poly.degree() == 0)
            return poly.coefficients[0];
        return cc.r2.core.number.Util.gcd(poly.coefficients);
    }

    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    primitivePart(UnivariatePolynomial<R> poly) {
        R content = content(poly);
        if (poly.lc().compareTo(poly.ring.getZero()) < 0)
            content = content.negate();
        return divide(poly, content);
    }

    private static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    divide(UnivariatePolynomial<R> poly, R content) {
        if (content.isOne())
            return poly;
        final R[] newData = newArray(poly.ring, poly.internalDegree);
        for (int i = 0; i < newData.length; i++) {
            R[] rs = poly.coefficients[i].divideAndRemainder(content);
            if (!rs[1].isZero())
                throw new IllegalArgumentException("" + rs[1]);
            newData[i] = rs[0];
        }
        return new UnivariatePolynomial<R>(poly.ring, newData);
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>[] pseudoDivideAndRemainder(
            final Ring<R> ring,
            UnivariatePolynomial<R> a,
            final UnivariatePolynomial<R> b) {
        return divideAndRemainder(ring, a.multiplyByFactor(b.lc().pow(a.degree() - b.degree() + 1)), b);
    }


    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>[] divideAndRemainder(
            final Ring<R> ring,
            final UnivariatePolynomial<R> a,
            final UnivariatePolynomial<R> b) {
        if (a.internalDegree < b.internalDegree)
            return null;

        MutableUnivariatePolynomial<R> r = new MutableUnivariatePolynomial<>(ring, a.degree(), a.coefficients.clone());

        final R[] q = newArray(ring, a.degree() - b.degree() + 1);
        for (int i = a.degree() - b.degree(); i >= 0; --i) {
            if (r.degree() == b.degree() + i) {
                final R[] qc = r.lc().divideAndRemainder(b.lc());
                if (!qc[1].isZero())
                    return null;
                q[i] = qc[0];
                r.subtract(b.coefficients, q[i], i);
            } else q[i] = ring.getZero();
        }
        return new UnivariatePolynomial[]{new UnivariatePolynomial<>(ring, q), r.toPoly()};
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> EuclideanPRS(Ring<R> ring,
                                               final UnivariatePolynomial<R> a,
                                               final UnivariatePolynomial<R> b) {
        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        UnivariatePolynomial<R> x = a, y = b, r;
        while (!y.isZero()) {
            UnivariatePolynomial<R>[] tmp = divideAndRemainder(ring, x, y);
            r = tmp[1];
            prs.add(r);
            x = y;
            y = r;
        }
        prs.add(x);
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> naivePRS(Ring<R> ring,
                                           final UnivariatePolynomial<R> a,
                                           final UnivariatePolynomial<R> b) {
        R aContent = content(a), bContent = content(b), contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        UnivariatePolynomial<R> x = aPP, y = bPP, r;
        while (!y.isZero()) {
            UnivariatePolynomial<R>[] tmp = pseudoDivideAndRemainder(ring, x, y);
            r = tmp[1];
            if (!r.isZero())
                r = primitivePart(r);
            prs.add(r);
            x = y;
            y = r;
        }
        if (x.degree() == 0)
            x = new UnivariatePolynomial<R>(ring, contentGCD);
        else
            x = primitivePart(x).multiplyByFactor(contentGCD);
        prs.add(x);
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> primitivePRS(final Ring<R> ring,
                                               final UnivariatePolynomial<R> a,
                                               final UnivariatePolynomial<R> b) {
        if (b.degree() > a.degree())
            return primitivePRS(ring, b, a);

        R aContent = content(a),
                bContent = content(b),
                contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        for (int i = 0; ; i++) {
            UnivariatePolynomial<R> ith = prs.get(i);
            UnivariatePolynomial<R> ithp1 = prs.get(i + 1);

            UnivariatePolynomial<R> q = divideAndRemainder(ring, ith.multiplyByFactor(ithp1.lc().pow(ith.degree() - ithp1.degree() + 1)), ithp1)[1];
            if (q.isZero())
                break;
            q = primitivePart(q);
            prs.add(q);
        }
        prs.set(prs.size() - 1, prs.get(prs.size() - 1).multiplyByFactor(contentGCD));
        return prs;
    }

    @SuppressWarnings({"unchecked", "ConstantConditions"})
    public static <R extends EuclideanRingElement<R>>
    List<UnivariatePolynomial<R>> subResultantPRS(final Ring<R> ring,
                                                  final UnivariatePolynomial<R> a,
                                                  final UnivariatePolynomial<R> b) {
        if (b.degree() > a.degree())
            return primitivePRS(ring, b, a);

        R aContent = content(a),
                bContent = content(b),
                contentGCD = aContent.gcd(bContent);
        UnivariatePolynomial<R> aPP = divide(a, aContent), bPP = divide(b, bContent);

        ArrayList<UnivariatePolynomial<R>> prs = new ArrayList<>();
        prs.add(aPP);
        prs.add(bPP);

        ArrayList<R> beta = new ArrayList<>(), psi = new ArrayList<>();
        TIntArrayList deltas = new TIntArrayList();

        R cBeta, cPsi;
        for (int i = 0; ; i++) {
            UnivariatePolynomial<R> ith = prs.get(i);
            UnivariatePolynomial<R> ithp1 = prs.get(i + 1);
            int delta = ith.degree() - ithp1.degree();
            if (i == 0) {
                cBeta = (delta + 1) % 2 == 0 ? ring.getOne() : ring.getOne().negate();
                cPsi = ring.getOne().negate();
            } else {
                cPsi = ith.lc().negate().pow(deltas.get(i - 1))
                ;//.multiply((R) BigInteger.valueOf(i + 1));
                if (deltas.get(i - 1) < 1)
                    cPsi = cPsi.multiply(psi.get(i - 1).pow(-deltas.get(i - 1) + 1));
                else {
                    R[] rs = cPsi.divideAndRemainder(psi.get(i - 1).pow(deltas.get(i - 1) - 1));
                    assert rs[1].isZero();
                    cPsi = rs[0];
                }
                cBeta = ith.lc().negate().multiply(cPsi.pow(delta));
            }

            UnivariatePolynomial<R> q = divideAndRemainder(ring, ith.multiplyByFactor(ithp1.lc().pow(ith.degree() - ithp1.degree() + 1)), ithp1)[1];
            if (q.isZero())
                break;

            q = divide(q, cBeta);
            prs.add(q);

            deltas.add(delta);
            beta.add(cBeta);
            psi.add(cPsi);
        }
        prs.set(prs.size() - 1, prs.get(prs.size() - 1).multiplyByFactor(contentGCD));
        return prs;
    }

    @SuppressWarnings("unchecked")
    public static <R extends EuclideanRingElement<R>> UnivariatePolynomial<R>
    EuclideanGCD(final Ring<R> ring,
                 final UnivariatePolynomial<R> a,
                 final UnivariatePolynomial<R> b) {
        UnivariatePolynomial<R> x = a, y = b, r;
        while (!y.isZero()) {
            r = divideAndRemainder(ring, x, y)[1];
            x = y;
            y = r;
        }
        return x;
    }
}
