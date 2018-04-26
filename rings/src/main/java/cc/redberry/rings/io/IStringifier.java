package cc.redberry.rings.io;

import cc.redberry.rings.Ring;
import cc.redberry.rings.poly.IPolynomial;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.multivar.AMultivariatePolynomial;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Defines #stringify(Stringifiable) method
 */
public interface IStringifier<Element> {
    /**
     * Stringify stringifiable object
     */
    @SuppressWarnings("unchecked")
    default String stringify(Element el) {
        if (el instanceof Stringifiable)
            return ((Stringifiable) el).toString(this);
        else
            return el.toString();
    }

    /**
     * Get stringifier for the specified ring of some underlying elements,
     * should never give null (use dummy() for absent stringifier)
     */
    <UnderlyingElement>
    IStringifier<UnderlyingElement> substringifier(Ring<UnderlyingElement> ring);

    /**
     * Get string binding of corresponding element
     */
    default String getBinding(Element el) {
        return getBinding(el, null);
    }

    /**
     * Get string binding of corresponding element
     *
     * @param defaultStr default string
     */
    default String getBinding(Element el, String defaultStr) {
        return getBindings().getOrDefault(el, defaultStr);
    }

    /**
     * Map of bindings
     */
    Map<Element, String> getBindings();

    //////////////////////////////////////////////////////Factory//////////////////////////////////////////////////////

    /** Dummy stringifier */
    IStringifier DUMMY = new IStringifier() {
        @Override
        public IStringifier substringifier(Ring ring) {
            return this;
        }

        @Override
        public Map getBindings() {
            return Collections.EMPTY_MAP;
        }
    };

    /** Dummy stringifier */
    @SuppressWarnings("unchecked")
    static <E> IStringifier<E> dummy() {
        return (IStringifier<E>) DUMMY;
    }

    /**
     * Simple map-based stringifier
     */
    final class SimpleStringifier<E> implements IStringifier<E> {
        public final Map<E, String> bindings = new HashMap<>();
        public final Map<Ring, IStringifier> substringifiers = new HashMap<>();

        @Override
        @SuppressWarnings("unchecked")
        public <U> IStringifier<U> substringifier(Ring<U> ring) {
            return (IStringifier<U>) substringifiers.getOrDefault(ring, IStringifier.<U>dummy());
        }

        @Override
        public Map<E, String> getBindings() {
            return bindings;
        }
    }

    /** Create simple stringifier */
    static <E> IStringifier<E> mkStringifier(Map<E, String> bindings) {
        SimpleStringifier<E> r = new SimpleStringifier<>();
        r.bindings.putAll(bindings);
        return r;
    }

    /** Create simple stringifier for polynomials with given variables */
    static <Poly extends IPolynomial<Poly>> IStringifier<Poly> mkPolyStringifier(IPolynomialRing<Poly> ring, String... variables) {
        Map<Poly, String> bindings = new HashMap<>();
        for (int i = 0; i < ring.nVariables(); ++i)
            bindings.put(ring.variable(i), variables[i]);
        return mkStringifier(bindings);
    }

    /** Create simple stringifier for polynomials with given variables */
    @SuppressWarnings("unchecked")
    static <Poly extends IPolynomial<Poly>> IStringifier<Poly> mkPolyStringifier(Poly factory, String... variables) {
        Map<Poly, String> bindings = new HashMap<>();
        if (factory instanceof IUnivariatePolynomial)
            bindings.put((Poly) ((IUnivariatePolynomial) factory).createMonomial(1), variables[0]);
        else {
            AMultivariatePolynomial mf = (AMultivariatePolynomial) factory;
            for (int i = 0; i < mf.nVariables; ++i)
                bindings.put((Poly) mf.createMonomial(i, 1), variables[i]);
        }
        return mkStringifier(bindings);
    }

    /**
     * Enclose with math parenthesis if needed (e.g. a+b in (a+b)*c should be enclosed)
     */
    static String encloseMathParenthesisInSumIfNeeded(String cf) {
        if (needParenthesisInSum(cf))
            return "(" + cf + ")";
        else
            return cf;
    }

    /**
     * If required to enclose with math parenthesis (e.g. a+b in (a+b)*c should be enclosed)
     */
    static boolean needParenthesisInSum(String cf) {
        if (cf.startsWith("-") && !hasPlusMinus(1, cf))
            return false;
        return hasPlusMinus(0, cf);
    }

    static boolean hasPlusMinus(int start, String cf) {
        // has +- on a zero bracket level
        int level = 0;
        for (int i = start; i < cf.length(); ++i) {
            char c = cf.charAt(i);
            if (c == '(')
                ++level;
            else if (c == ')')
                --level;
            else if ((c == '+' || c == '-') && level == 0)
                return true;
        }
        return false;
    }

    static boolean hasMulDivPlusMinus(int start, String cf) {
        // has +- on a zero bracket level
        int level = 0;
        for (int i = start; i < cf.length(); ++i) {
            char c = cf.charAt(i);
            if (c == '(')
                ++level;
            else if (c == ')')
                --level;
            else if (level == 0 && (c == '+' || c == '-' || c == '*' || c == '/'))
                return true;
        }
        return false;
    }

    /**
     * Sequence of strings "a", "b", "c" etc.
     *
     * @param nVars number of variable
     */
    static String[] defaultVars(int nVars) {
        if (nVars == 1)
            return new String[]{"x"};
        if (nVars == 2)
            return new String[]{"x", "y"};
        if (nVars == 3)
            return new String[]{"x", "y", "z"};

        String[] vars = new String[nVars];
        for (int i = 1; i <= nVars; i++)
            vars[i - 1] = "x" + i;
        return vars;
    }

    static String defaultVar(int i, int nVars) {
        if (nVars == 1)
            return "x";
        if (nVars == 2)
            return i == 0 ? "x" : "y";
        if (nVars == 3)
            return i == 0 ? "x" : i == 1 ? "y" : "z";

        return "x" + (i + 1);
    }

    static String defaultVar() {
        return defaultVar(0, 1);
    }
}
