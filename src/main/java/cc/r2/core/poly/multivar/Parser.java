package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Integers;
import cc.r2.core.util.ArraysUtil;
import gnu.trove.list.array.TIntArrayList;

import java.util.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class Parser {
    private Parser() {}

    static MultivariatePolynomial<BigInteger> parse(String input, Comparator<DegreeVector> ordering, String... variables) {
        return parse(input, Integers.Integers, ordering, variables);
    }

    static <E> MultivariatePolynomial<E> parse(String input, Domain<E> domain) {
        return parse(input, domain, MonomialOrder.LEX);
    }

    static <E> MultivariatePolynomial<E> parse(String input, Domain<E> domain, Comparator<DegreeVector> ordering, String... variables) {
        List<TMonomialTerm<E>> terms = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        int bracketLevel = 0;
        for (int i = 0; i < input.length(); i++) {
            char c = input.charAt(i);
            if (c == '(')
                ++bracketLevel;
            if (c == ')')
                --bracketLevel;
            if ((c == '+' || c == '-') && bracketLevel == 0 && sb.length() != 0) {
                terms.add(parseMonomial(sb.toString(), domain));
                sb = new StringBuilder();
            }
            sb.append(c);
        }
        if (sb.length() != 0)
            terms.add(parseMonomial(sb.toString(), domain));

        Set<String> allVars = new HashSet<>();
        allVars.addAll(Arrays.asList(variables));
        terms.forEach(t -> allVars.addAll(Arrays.asList(t.variables)));

        List<String> varsList = new ArrayList<>();
        varsList.addAll(Arrays.asList(variables));
        for (String var : allVars) {
            if (!varsList.contains(var))
                varsList.add(var);
        }
        String[] vars = varsList.toArray(new String[allVars.size()]);

        @SuppressWarnings("unchecked")
        Monomial<E>[] mTerms = new Monomial[terms.size()];
        for (int i = 0; i < terms.size(); i++)
            mTerms[i] = terms.get(i).toMonomialTerm(vars);

        return MultivariatePolynomial.create(mTerms[0].exponents.length, domain, ordering, mTerms);
    }

    public static final class TMonomialTerm<E> {
        final String[] variables;
        final E coefficient;
        final int[] exponents;

        public TMonomialTerm(String[] variables, E coefficient, int[] exponents) {
            this.variables = variables;
            this.coefficient = coefficient;
            this.exponents = exponents;
        }

        Monomial<E> toMonomialTerm(String[] map) {
            int[] degrees = new int[map.length];
            for (int i = 0; i < variables.length; i++)
                degrees[ArraysUtil.firstIndexOf(variables[i], map)] = exponents[i];
            return new Monomial<>(degrees, coefficient);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(coefficient);
            if (variables.length > 0) {
                for (int i = 0; i < variables.length; i++) {
                    sb.append("*").append(variables[i]);
                    if (exponents[i] != 1)
                        sb.append("^").append(exponents[i]);
                }
            }
            return sb.toString();
        }
    }

    public static <E> TMonomialTerm<E> parseMonomial(String expression, Domain<E> domain) {
        expression = expression.replace(" ", "");
        E coefficient = domain.getOne();
        if (expression.startsWith("+"))
            expression = expression.substring(1);
        if (expression.startsWith("-")) {
            coefficient = domain.negate(coefficient);
            expression = expression.substring(1);
        }
        String[] elements = splitMultipliers(expression);

        ArrayList<String> variables = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        for (String element : elements) {
            String el = element.trim();
            if (el.startsWith("(") || el.matches("\\d+")) {
                if (el.startsWith("(") && el.endsWith(")"))
                    el = el.substring(1, el.length() - 1);
                try {
                    coefficient = domain.multiply(coefficient, domain.parse(el));
                    continue;
                } catch (Exception e) {}  // not a coefficient
            }

            // variable with exponent
            String[] varExp = el.split("\\^");
            String var = varExp[0].trim();
            if (varExp.length == 2) {
                if (var.startsWith("(") || var.matches("\\d+"))
                    try {
                        coefficient = domain.multiply(coefficient, domain.pow(domain.parse(var), Integer.parseInt(varExp[1].trim())));
                        continue;
                    } catch (Exception e) {}
            }
            variables.add(var);
            exponents.add(varExp.length == 1 ? 1 : Integer.parseInt(varExp[1].trim()));
        }

        return new TMonomialTerm<>(variables.toArray(new String[variables.size()]), coefficient, exponents.toArray());
    }

    private static String[] splitMultipliers(String str) {
        List<String> result = new ArrayList<>();
        int bracketsLevel = 0;
        int prevPosition = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == '(')
                ++bracketsLevel;
            else if (str.charAt(i) == ')')
                --bracketsLevel;
            else if (str.charAt(i) == '*' && bracketsLevel == 0) {
                result.add(str.substring(prevPosition, i));
                prevPosition = i + 1;
            }
        }
        result.add(str.substring(prevPosition, str.length()));
        return result.toArray(new String[result.size()]);
    }
}
