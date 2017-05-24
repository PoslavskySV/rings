package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.poly.Domain;
import cc.r2.core.poly.Integers;
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
        return parse(input, domain, DegreeVector.LEX);
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
        terms.forEach(t -> allVars.addAll(Arrays.asList(t.variables)));
        allVars.addAll(Arrays.asList(variables));

        String[] vars = allVars.toArray(new String[allVars.size()]);
        Arrays.sort(vars);

        @SuppressWarnings("unchecked")
        MonomialTerm<E>[] mTerms = new MonomialTerm[terms.size()];
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

        MonomialTerm<E> toMonomialTerm(String[] map) {
            int[] degrees = new int[map.length];
            for (int i = 0; i < variables.length; i++)
                degrees[Arrays.binarySearch(map, variables[i])] = exponents[i];
            return new MonomialTerm<>(degrees, coefficient);
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
        String[] elements = expression.split("\\*");

        ArrayList<String> variables = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();
        for (String element : elements) {
            String el = element.trim();
            try {
                coefficient = domain.multiply(coefficient, domain.parse(el));
                continue;
            } catch (Exception e) {}  // not a coefficient

            // variable with exponent
            String[] varExp = el.split("\\^");
            if (varExp.length == 2) {
                try {
                    coefficient = domain.multiply(coefficient, domain.pow(domain.parse(varExp[0].trim()), Integer.parseInt(varExp[1].trim())));
                    continue;
                } catch (Exception e) {}
            }
            variables.add(varExp[0].trim());
            exponents.add(varExp.length == 1 ? 1 : Integer.parseInt(varExp[1].trim()));
        }

        return new TMonomialTerm<>(variables.toArray(new String[variables.size()]), coefficient, exponents.toArray());
    }
}
