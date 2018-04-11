package cc.redberry.rings.poly.univar;

import cc.redberry.rings.ElementParser;
import cc.redberry.rings.Ring;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Parser for univariate polynomials.
 *
 * @since 1.0
 */
final class Parser {
    private Parser() {}

    public static final Pattern termPattern = Pattern.compile("((?<var>[a-zA-Zz]+)(\\^(?<exponent>[0-9]+))?)?$");

    public static <E> UnivariatePolynomial<E> parse(Ring<E> ring, String string) {
        return parse(ring, ring, string);
    }

    public static <E> UnivariatePolynomial<E> parse(Ring<E> ring, ElementParser<E> eParser, String string) {
        string = string.replace(" ", "");
        string = string.replace("--", "+");
        string = string.replace("-+", "-");
        string = string.replace("+-", "-");
        string = string.replace("++", "+");
        int bracketLevel = 0;
        StringBuilder buffer = new StringBuilder();

        List<E> coefficients = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        String[] varRef = {null};
        for (int i = 0; i < string.length(); i++) {
            char c = string.charAt(i);
            if (c == '(')
                ++bracketLevel;
            if (c == ')')
                --bracketLevel;

            if ((c == '+' || c == '-') && bracketLevel == 0) {
                parseTerm(ring, eParser, buffer.toString(), coefficients, exponents, varRef);
                buffer = new StringBuilder();
            }

            buffer.append(c);
        }

        parseTerm(ring, eParser, buffer.toString(), coefficients, exponents, varRef);
        int degree = exponents.size() == 0 ? 0 : exponents.max();
        E[] data = ring.createZeroesArray(degree + 1);
        for (int i = 0; i < coefficients.size(); i++)
            data[exponents.get(i)] = ring.addMutable(data[exponents.get(i)], coefficients.get(i));

        return UnivariatePolynomial.create(ring, data);
    }

    private static int separatorPosition(String string) {
        int bracketLevel = 0;
        for (int i = string.length() - 1; i >= 0; --i) {
            char c = string.charAt(i);
            if (c == '(')
                ++bracketLevel;
            if (c == ')')
                --bracketLevel;
            if (c == '*' && bracketLevel == 0)
                return i;
        }
        return -1;
    }

    private static <E> void parseTerm(Ring<E> ring, ElementParser<E> eParser, String string, List<E> coefficients, TIntArrayList exponents, String[] varRef) {
        if (string.isEmpty())
            return;

        int p = separatorPosition(string);

        String coefficient, exponent;
        if (p == -1) {
            if (string.contains("(") || string.contains(")") || string.matches("[+-]?\\d+")) {
                coefficient = string;
                exponent = "";
            } else {
                coefficient = "";
                exponent = string;
            }
        } else if (string.substring(p).contains(")")) {
            coefficient = string;
            exponent = "";
        } else {
            coefficient = string.substring(0, p);
            exponent = string.substring(p + 1, string.length());
        }

        if (coefficient.isEmpty()) {
            if (exponent.startsWith("-")) {
                coefficient = "-1";
                exponent = exponent.substring(1);
            } else if (exponent.startsWith("+"))
                exponent = exponent.substring(1);
        }

        coefficients.add(parseCoeff(ring, eParser, coefficient));
        exponents.add(parseExponent(exponent, varRef));
    }

    private static <E> E parseCoeff(Ring<E> ring, ElementParser<E> eParser, String string) {
        boolean negate = false;
        if (string.startsWith("+"))
            string = string.substring(1);
        if (string.startsWith("-")) {
            string = string.substring(1);
            negate = true;
        }
        E result;
        if (string.isEmpty())
            result = ring.getOne();
        else if (string.startsWith("(") && string.endsWith(")"))
            result = eParser.parse(string.substring(1, string.length() - 1));
        else
            result = eParser.parse(string);
        return negate ? ring.negateMutable(result) : result;
    }

    private static int parseExponent(String string, String[] varRef) {
        if (string.isEmpty())
            return 0;
        Matcher matcher = termPattern.matcher(string);
        if (matcher.find()) {
            String var = matcher.group("var");
            String exponent = matcher.group("exponent");
            if (var == null || (varRef[0] != null && !var.equals(varRef[0])))
                throw new RuntimeException("var = " + var + "  ref = " + varRef[0] + " string = " + string);
            if (varRef[0] == null)
                varRef[0] = var;
            if (exponent != null)
                return Integer.parseInt(exponent);
            else
                return 1;
        }
        return 0;
    }
}