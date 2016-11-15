package cc.r2.core.number;


public class ArithmeticParser {
//    static final int LPAREN = 0, RPAREN = 1, PLUS = 2, MINUS = 3, MUL = 4, DIV = 5, NUM = 6, SPACE = 7;
//
//    static int parse(char c) {
//        switch (c) {
//            case '(':
//                return LPAREN;
//            case ')':
//                return RPAREN;
//            case '+':
//                return PLUS;
//            case '-':
//                return MINUS;
//            case '*':
//                return MUL;
//            case '/':
//                return DIV;
//            case ' ':
//                return SPACE;
//        }
//        if (Character.isDigit(c))
//            return NUM;
//
//        throw new RuntimeException(String.valueOf(c));
//    }

    public static <R extends RingElement<R>> R parseBrackets(String expression, Ring<R> ring) {
        if (expression.startsWith("(") && expression.endsWith(")")) {
            int pLevel = 0;
            for (int i = 1; i < expression.length() - 1; i++) {
                char c = expression.charAt(i);
                if (c == '(')
                    ++pLevel;
                else if (c == ')')
                    --pLevel;

                if (pLevel < 0)
                    return null;
            }
            if (pLevel != 0)
                return null;
            return parse(expression.substring(1, expression.length() - 1), ring);
        }
        return null;
    }

    public static <R extends RingElement<R>> R parseNum(String expression, Ring<R> ring) {
        if (expression.matches("[0-9]+"))
            return ring.parse(expression);
        return null;
    }

    public static <R extends RingElement<R>> R parseMinus(String expression, Ring<R> ring) {
        if (expression.matches("-[0-9]+"))
            return ring.parse(expression).negate();
        return null;
    }

    public static <R extends RingElement<R>> R parsePlus(String expression, Ring<R> ring) {
        if (!expression.contains("+") && !expression.contains("-"))
            return null;

        R element = ring.getZero();
        StringBuilder buffer = new StringBuilder();

        int pLevel = 0;
        boolean canParse = false;
        for (int i = 0; i < expression.length(); i++) {
            char c = expression.charAt(i);
            if (c == '(')
                ++pLevel;
            if (c == ')')
                --pLevel;
            if (pLevel == 0 && (c == '+' || c == '-')) {
                canParse = true;
            }
            if (pLevel < 0)
                return null;
        }
        if (!canParse)
            return null;

        char lastOperator = '+';
        for (int i = 0; i < expression.length(); i++) {
            char c = expression.charAt(i);
            if (c == '(')
                ++pLevel;
            if (c == ')')
                --pLevel;
            if (pLevel == 0 && (c == '+' || c == '-')) {
                if (lastOperator == '+')
                    element = element.add(parse(buffer.toString(), ring));
                else if (lastOperator == '-')
                    element = element.subtract(parse(buffer.toString(), ring));
                lastOperator = c;
                buffer = new StringBuilder();
            } else buffer.append(c);
        }

        if (pLevel != 0)
            return null;
        if (buffer.length() != 0) {
            if (lastOperator == '+')
                element = element.add(parse(buffer.toString(), ring));
            else if (lastOperator == '-')
                element = element.subtract(parse(buffer.toString(), ring));
        }
        return element;
    }

    public static <R extends RingElement<R>> R parseMul(String expression, Ring<R> ring) {
        if (!expression.contains("*") && !expression.contains("/"))
            return null;

        R element = ring.getOne();
        StringBuilder buffer = new StringBuilder();

        int pLevel = 0;
        boolean canParse = false;
        for (int i = 0; i < expression.length(); i++) {
            char c = expression.charAt(i);
            if (c == '(')
                ++pLevel;
            if (c == ')')
                --pLevel;
            if (pLevel == 0 && (c == '*' || c == '/')) {
                canParse = true;
            }
            if (pLevel < 0)
                return null;
        }
        if (!canParse)
            return null;

        char lastOperator = '*';
        for (int i = 0; i < expression.length(); i++) {
            char c = expression.charAt(i);
            if (c == '(')
                ++pLevel;
            if (c == ')')
                --pLevel;
            if (pLevel == 0 && (c == '*' || c == '/')) {
                if (lastOperator == '*')
                    element = element.multiply(parse(buffer.toString(), ring));
                else if (lastOperator == '/')
                    element = (R) (((FieldElement) element).divide((FieldElement) parse(buffer.toString(), ring)));
                lastOperator = c;
                buffer = new StringBuilder();
            } else buffer.append(c);
        }

        if (pLevel != 0)
            return null;
        if (buffer.length() != 0) {
            if (lastOperator == '*')
                element = element.multiply(parse(buffer.toString(), ring));
            else if (lastOperator == '/')
                element = (R) ((FieldElement) element).divide((FieldElement) parse(buffer.toString(), ring));
        }
        return element;
    }

    public static <R extends RingElement<R>> R parse(String expression, Ring<R> ring) {
        expression = expression.replaceAll(" ", "");
        R result = parseBrackets(expression, ring);
        if (result != null)
            return result;
        result = parseNum(expression, ring);
        if (result != null)
            return result;
        result = parsePlus(expression, ring);
        if (result != null)
            return result;
        result = parseMul(expression, ring);
        if (result != null)
            return result;
        return result;
    }
}
