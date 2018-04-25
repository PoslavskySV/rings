package cc.redberry.rings.io;

import cc.redberry.rings.Ring;

import java.util.Map;

/**
 * Defines #stringify(Stringifiable) method
 */
public interface IStringifier<Element> {
    /**
     * Stringify stringifiable object
     */
    default String stringify(Stringifiable<Element> el) {
        return el.toString(this);
    }

    /**
     * Get stringifier for the specified ring of some underlying elements
     */
    <UnderlyingElement> IStringifier<UnderlyingElement> substringifier(Ring<UnderlyingElement> ring);

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

    /**
     * Enclose with math parenthesis if needed (e.g. a+b in (a+b)*c should be enclosed)
     */
    static String encloseMathParenthesisIfNeeded(String cf) {
        if (needParenthesis(cf))
            return "(" + cf + ")";
        else
            return cf;
    }

    /**
     * If required to enclose with math parenthesis (e.g. a+b in (a+b)*c should be enclosed)
     */
    static boolean needParenthesis(String cf) {
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
}
