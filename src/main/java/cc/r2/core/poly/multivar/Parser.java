package cc.r2.core.poly.multivar;

import cc.r2.core.number.BigInteger;
import cc.r2.core.number.BigIntegerArithmetics;
import gnu.trove.list.array.TIntArrayList;

import java.util.*;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
final class Parser {
    private Parser() {}

    static MultivariatePolynomial parse(String input, Comparator<MultivariatePolynomial.DegreeVector> ordering, String... variables) {
        List<ParsedTerm> terms = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < input.length(); i++) {
            char c = input.charAt(i);
            if ((c == '+' || c == '-') && sb.length() != 0) {
                terms.add(parseTerm(sb.toString()));
                sb = new StringBuilder();
            }
            sb.append(c);
        }
        if (sb.length() != 0)
            terms.add(parseTerm(sb.toString()));

        Set<String> allVars = new HashSet<>();
        terms.forEach(t -> allVars.addAll(Arrays.asList(t.variables)));
        allVars.addAll(Arrays.asList(variables));

        String[] vars = allVars.toArray(new String[allVars.size()]);
        Arrays.sort(vars);

        BigInteger[] factor = new BigInteger[terms.size()];
        MultivariatePolynomial.DegreeVector[] vectors = new MultivariatePolynomial.DegreeVector[terms.size()];
        for (int i = 0; i < terms.size(); i++) {
            factor[i] = terms.get(i).factor;
            vectors[i] = terms.get(i).degreeVector(vars);
        }

        return MultivariatePolynomial.create(factor, vectors, ordering);
    }


    static ParsedTerm parseTerm(String input) {
        input = input.trim();

        ArrayList<String> variables = new ArrayList<>();
        TIntArrayList exponents = new TIntArrayList();

        char sign = input.charAt(0);
        int i = (sign == '+' || sign == '-') ? 1 : 0;
        StringBuilder currentBuffer = new StringBuilder();
        for (; i < input.length(); ++i) {
            char c = input.charAt(i);
            if (c == ' ')
                continue;

            if (c == '*' || i == input.length() - 1) {
                if (i == input.length() - 1)
                    currentBuffer.append(c);
                String variable = currentBuffer.toString();
                String[] varExp = variable.split("\\^");
                if (varExp.length == 1) {
                    variables.add(variable);
                    exponents.add(1);
                } else {
                    variables.add(varExp[0]);
                    exponents.add(Integer.parseInt(varExp[1]));
                }
                currentBuffer = new StringBuilder();
            } else
                currentBuffer.append(c);
        }

        BigInteger factor = sign == '-' ? BigInteger.NEGATIVE_ONE : BigInteger.ONE;
        for (i = variables.size() - 1; i >= 0; --i) {
            String s = variables.get(i).trim();
            if (s.matches("^\\d+$")) {
                factor = factor.multiply(BigIntegerArithmetics.pow(new BigInteger(s), exponents.get(i)));
                variables.remove(i);
                exponents.removeAt(i);
            }
        }

        return new ParsedTerm(factor, variables.toArray(new String[variables.size()]), exponents.toArray());
    }

    private static final class ParsedTerm {
        final BigInteger factor;
        final String[] variables;
        final int[] exponents;

        public ParsedTerm(BigInteger factor, String[] variables, int[] exponents) {
            this.factor = factor;
            this.variables = variables;
            this.exponents = exponents;
        }

        public MultivariatePolynomial.DegreeVector degreeVector(String[] map) {
            int[] degrees = new int[map.length];
            for (int i = 0; i < variables.length; i++)
                degrees[Arrays.binarySearch(map, variables[i])] = exponents[i];
            return new MultivariatePolynomial.DegreeVector(degrees);
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(factor);
            for (int i = 0; i < exponents.length; i++)
                sb.append("*").append(variables[i]).append("^").append(exponents[i]);
            return sb.toString();
        }
    }
}
