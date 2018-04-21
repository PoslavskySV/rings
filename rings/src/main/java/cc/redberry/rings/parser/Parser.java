package cc.redberry.rings.parser;

import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.parser.Tokenizer.Token;
import cc.redberry.rings.parser.Tokenizer.TokenType;
import cc.redberry.rings.poly.IPolynomialRing;
import cc.redberry.rings.poly.MultivariateRing;
import cc.redberry.rings.poly.UnivariateRing;
import cc.redberry.rings.poly.multivar.*;
import cc.redberry.rings.poly.univar.IUnivariatePolynomial;
import cc.redberry.rings.poly.univar.UnivariatePolynomial;

import java.util.ArrayDeque;
import java.util.HashMap;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import static cc.redberry.rings.parser.Tokenizer.END;
import static cc.redberry.rings.parser.Tokenizer.TokenType.T_BRACKET_OPEN;
import static cc.redberry.rings.parser.Tokenizer.TokenType.T_VARIABLE;


/**
 * High-level parser of ring elements. Uses simple shunting-yard algorithm.
 *
 * @param <Element> type of resulting elements
 * @param <Term>    underlying polynomial terms
 * @param <Poly>    underlying multivariate polynomials
 */
public final class Parser<
        Element,
        Term extends AMonomial<Term>,
        Poly extends AMultivariatePolynomial<Term, Poly>> {
    /** the base ring */
    private final Ring<Element> baseRing;
    /** map variableName -> Element (if it is a polynomial variable) */
    private final Map<String, Element> eVariables;
    /** auxiliary polynomial ring */
    private final MultivariateRing<Poly> polyRing;
    /** map variableName -> variableIndex (if it is a polynomial variable) */
    private final Map<String, Integer> pVariables;
    /** convert polynomial to base ring elements */
    private final Function<Poly, Element> polyToElement;

    private Parser(Ring<Element> baseRing,
                   Map<String, Element> eVariables,
                   MultivariateRing<Poly> polyRing,
                   Map<String, Integer> pVariables,
                   Function<Poly, Element> polyToElement) {
        this.baseRing = baseRing;
        this.eVariables = eVariables;
        this.polyRing = polyRing;
        this.pVariables = pVariables;
        this.polyToElement = polyToElement;
        // make sure that eVariables contain all pVariables
        if (pVariables != null)
            pVariables.forEach((k, v) -> eVariables.computeIfAbsent(k, __ -> polyToElement.apply(polyRing.variable(v))));
    }

    ////////////////////////////////////////////////////Factory////////////////////////////////////////////////////////

    /**
     * @param baseRing      the base ring
     * @param eVariables    variables bindings (variableString -> base ring element)
     * @param polyRing      auxiliary polynomial ring, to manage intermediate polynomial expressions
     * @param pVariables    polynomial variables bindings (variableString -> polyRing variable index)
     * @param polyToElement convert from auxiliary polynomials to basering
     */
    public static <
            Element,
            Term extends AMonomial<Term>,
            Poly extends AMultivariatePolynomial<Term, Poly>> Parser<Element, Term, Poly>
    mkParser(Ring<Element> baseRing,
             Map<String, Element> eVariables,
             MultivariateRing<Poly> polyRing,
             Map<String, Poly> pVariables,
             Function<Poly, Element> polyToElement) {
        Map<String, Integer> iVariables = new HashMap<>();
        if (pVariables != null)
            for (Map.Entry<String, Poly> v : pVariables.entrySet()) {
                Poly p = v.getValue();
                if (p.isEffectiveUnivariate() && !p.isConstant() && p.degreeSum() == 1)
                    iVariables.put(v.getKey(), v.getValue().univariateVariable());
                else
                    eVariables.put(v.getKey(), polyToElement.apply(p));
            }
        return new Parser<>(baseRing, eVariables, polyRing, iVariables, polyToElement);
    }

    /**
     * Create parser for generic rings
     */
    public static <E> Parser<E, ?, ?> mkGenericParser(Ring<E> baseRing,
                                                      Map<String, E> eVariables) {
        return new Parser<>(baseRing, eVariables, null, null, null);
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Parser<Poly, Term, Poly> mkMultivariateParser(MultivariateRing<Poly> baseRing,
                                                  Map<String, Poly> eVariables) {
        return mkParser(baseRing, eVariables, baseRing, eVariables, Function.identity());
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Parser<Poly, Term, Poly> mkMultivariateParser(MultivariateRing<Poly> baseRing,
                                                  String... variables) {
        Map<String, Poly> pVariables = new HashMap<>();
        for (int i = 0; i < variables.length; ++i)
            pVariables.put(variables[i], baseRing.variable(i));
        return mkMultivariateParser(baseRing, pVariables);
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <E> Parser<MultivariatePolynomial<E>, Monomial<E>, MultivariatePolynomial<E>>
    mkMultivariateParser(MultivariateRing<MultivariatePolynomial<E>> baseRing,
                         Map<String, MultivariatePolynomial<E>> eVariables,
                         Parser<E, ?, ?> elementParser) {
        elementParser.eVariables.forEach((k, v) -> eVariables.put(k, baseRing.factory().createConstant(v)));
        return mkMultivariateParser(baseRing, eVariables);
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <E> Parser<MultivariatePolynomial<E>, Monomial<E>, MultivariatePolynomial<E>>
    mkMultivariateParser(MultivariateRing<MultivariatePolynomial<E>> baseRing,
                         Parser<E, ?, ?> elementParser,
                         String... variables) {
        Map<String, MultivariatePolynomial<E>> eVariables = new HashMap<>();
        for (int i = 0; i < variables.length; ++i)
            eVariables.put(variables[i], baseRing.variable(i));
        return mkMultivariateParser(baseRing, eVariables, elementParser);
    }


    /**
     * Create parser for multivariate polynomials rings
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>>
    Parser<Poly, ?, ?> mkUnivariateParser(IPolynomialRing<Poly> baseRing,
                                          Map<String, Poly> eVariables) {
        MultivariateRing mRing = Rings.MultivariateRing(baseRing.factory().asMultivariate());
        Map<String, AMultivariatePolynomial> pVariables = eVariables.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().asMultivariate()));
        return mkParser(baseRing, eVariables, mRing, pVariables, p -> (Poly) p.asUnivariate());
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>>
    Parser<Poly, ?, ?> mkUnivariateParser(UnivariateRing<Poly> baseRing,
                                          String variable) {
        return mkUnivariateParser(baseRing, new HashMap<String, Poly>() {{put(variable, baseRing.variable(0));}});
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <E> Parser<UnivariatePolynomial<E>, ?, ?>
    mkUnivariateParser(IPolynomialRing<UnivariatePolynomial<E>> baseRing,
                       Map<String, UnivariatePolynomial<E>> eVariables,
                       Parser<E, ?, ?> elementParser) {
        elementParser.eVariables.forEach((k, v) -> eVariables.put(k, baseRing.factory().createConstant(v)));
        return mkUnivariateParser(baseRing, eVariables);
    }

    /**
     * Create parser for multivariate polynomials rings
     */
    public static <E> Parser<UnivariatePolynomial<E>, ?, ?>
    mkUnivariateParser(IPolynomialRing<UnivariatePolynomial<E>> baseRing,
                       Parser<E, ?, ?> elementParser,
                       String variable) {
        HashMap<String, UnivariatePolynomial<E>> eVariables = new HashMap<>();
        eVariables.put(variable, baseRing.variable(0));
        return mkUnivariateParser(baseRing, eVariables, elementParser);
    }

    /**
     * Create parser for nested rings (e.g. fractions over polynomials etc)
     */
    @SuppressWarnings("unchecked")
    public static <E, I>
    Parser<E, ?, ?> mkNestedParser(Ring<E> baseRing,
                                   Map<String, E> eVariables,
                                   Parser<I, ?, ?> innerParser,
                                   Function<I, E> imageFunc) {
        if (baseRing instanceof MultivariateRing && ((MultivariateRing) baseRing).factory() instanceof MultivariatePolynomial)
            return mkMultivariateParser((MultivariateRing) baseRing, (Map) eVariables, innerParser);
        else if (baseRing instanceof UnivariateRing && ((UnivariateRing) baseRing).factory() instanceof UnivariatePolynomial)
            return mkUnivariateParser((UnivariateRing) baseRing, (Map) eVariables, innerParser);

        innerParser.eVariables.forEach((k, v) -> eVariables.put(k, imageFunc.apply(v)));
        return new Parser(
                baseRing,
                eVariables,
                innerParser.polyRing,
                innerParser.pVariables,
                innerParser.polyToElement.andThen(imageFunc));
    }

    ///////////////////////////////////////////////////Implementation///////////////////////////////////////////////////

    /** Parse string */
    public Element parse(String string) {
        return parse(Tokenizer.mkTokenizer(string));
    }

    /**
     * Parse stream of tokens into ring element
     */
    public Element parse(Tokenizer tokenizer) {
        // operators stack
        ArrayDeque<Operator> operators = new ArrayDeque<>();
        // operands stack
        ArrayDeque<IOperand<Poly, Element>> operands = new ArrayDeque<>();

        TokenType previousToken = null;
        Token token;
        while ((token = tokenizer.nextToken()) != END) {
            TokenType tType = token.tokenType;
            // first token is open bracket
            assert previousToken != null || tType == T_BRACKET_OPEN;

            // if this is variable, push it to stack
            if (tType == T_VARIABLE)
                operands.push(mkOperand(token));
            else {
                Operator op = tokenToOp(token.tokenType);
                assert op != null;

                // manage unary operations
                if (Operator.isPlusMinus(op) && previousToken == T_BRACKET_OPEN)
                    operands.push(Zero); // push neutral element

                if (op == Operator.BRACKET_CLOSE) {
                    while (!operators.isEmpty() && operators.peek() != Operator.BRACKET_OPEN)
                        popEvaluate(operators, operands);

                    operators.pop(); // remove opening bracket
                } else {

                    while (canPop(op, operators))
                        popEvaluate(operators, operands);

                    operators.push(op);
                }
            }

            previousToken = tType;
        }

        if (operands.size() > 1 || operators.size() > 1)
            throw new IllegalArgumentException("Can't parse");

        return operands.pop().toElement();
    }

    /** parse operand token */
    private IOperand<Poly, Element> mkOperand(Token operand) {
        if (isInteger(operand.content))
            return new NumberOperand(new BigInteger(operand.content));

        if (operand.tokenType != T_VARIABLE)
            throw new RuntimeException("illegal operand: " + operand);

        // if polynomial
        Integer iVar = pVariables.get(operand.content);
        if (iVar != null)
            return new VarOperand(pVariables.get(operand.content));

        // if base ring element
        Element eVar = eVariables.get(operand.content);
        if (eVar != null)
            return new ElementOperand(baseRing.copy(eVar));

        throw new RuntimeException("illegal operand: " + operand);
    }

    /** whether can pop element from ops stack */
    private boolean canPop(Operator op, ArrayDeque<Operator> opsStack) {
        if (opsStack.isEmpty())
            return false;
        int pOp = op.priority;
        int pOpPrev = opsStack.peek().priority;

        if (pOp < 0 || pOpPrev < 0)
            return false;

        return (op.associativity == Associativity.LeftToRight && pOp >= pOpPrev)
                || (op.associativity == Associativity.RightToLeft && pOp > pOpPrev);
    }

    /** pop two elements from operands stack and apply binary op from ops stack */
    private void popEvaluate(ArrayDeque<Operator> opsStack, ArrayDeque<IOperand<Poly, Element>> exprStack) {
        IOperand<Poly, Element>
                right = exprStack.pop(),
                left = exprStack.pop(),
                result;

        Operator op = opsStack.pop();

        switch (op) {
            case PLUS:
                result = left.plus(right);
                break;
            case MINUS:
                result = left.minus(right);
                break;
            case MULTIPLY:
                result = left.multiply(right);
                break;
            case DIVIDE:
                result = left.divide(right);
                break;
            case POWER:
                if (!(right instanceof Parser.NumberOperand))
                    throw new IllegalArgumentException("Exponents must be positive integers, but got " + right.toElement());
                result = left.pow(((NumberOperand) right).number);
                break;
            default: throw new RuntimeException();
        }

        exprStack.push(result);
    }

    //////////////////////////////////////////////////////Operands//////////////////////////////////////////////////////

    /** optimized operations with operands */
    private interface IOperand<P, E> {
        /** to auxiliary poly */
        P toPoly();

        /** to base ring element */
        E toElement();

        /** whether this operand is already converted to a base ring element */
        default boolean inBaseRing() {
            return this instanceof Parser.ElementOperand;
        }

        IOperand<P, E> plus(IOperand<P, E> oth);

        IOperand<P, E> minus(IOperand<P, E> oth);

        IOperand<P, E> multiply(IOperand<P, E> oth);

        IOperand<P, E> divide(IOperand<P, E> oth);

        IOperand<P, E> pow(BigInteger exponent);
    }

    /** default implementation of operands algebra */
    private abstract class DefaultOperandOps implements IOperand<Poly, Element> {
        @Override
        public Element toElement() {
            return polyToElement.apply(toPoly());
        }

        @Override
        public IOperand<Poly, Element> plus(IOperand<Poly, Element> oth) {
            if (inBaseRing() || oth.inBaseRing())
                return new ElementOperand(baseRing.addMutable(toElement(), oth.toElement()));
            else
                return new PolyOperand(polyRing.addMutable(toPoly(), oth.toPoly()));
        }

        @Override
        public IOperand<Poly, Element> minus(IOperand<Poly, Element> oth) {
            if (inBaseRing() || oth.inBaseRing())
                return new ElementOperand(baseRing.subtractMutable(toElement(), oth.toElement()));
            else
                return new PolyOperand(polyRing.subtractMutable(toPoly(), oth.toPoly()));
        }

        @Override
        public IOperand<Poly, Element> divide(IOperand<Poly, Element> oth) {
            return new ElementOperand(baseRing.divideExactMutable(toElement(), oth.toElement()));
        }

        @Override
        public IOperand<Poly, Element> multiply(IOperand<Poly, Element> oth) {
            if (inBaseRing() || oth.inBaseRing())
                return new ElementOperand(baseRing.multiplyMutable(toElement(), oth.toElement()));
            else
                return new PolyOperand(polyRing.multiplyMutable(toPoly(), oth.toPoly()));
        }

        @Override
        public IOperand<Poly, Element> pow(BigInteger exponent) {
            if (inBaseRing())
                return new ElementOperand(baseRing.pow(toElement(), exponent));
            else
                return new PolyOperand(polyRing.pow(toPoly(), exponent));
        }
    }

    /** zero operand */
    private final NumberOperand Zero = new NumberOperand(BigInteger.ZERO);

    /** A single number */
    private final class NumberOperand extends DefaultOperandOps implements IOperand<Poly, Element> {
        final BigInteger number;

        NumberOperand(BigInteger number) { this.number = number; }

        @Override
        public Poly toPoly() {
            return polyRing.valueOfBigInteger(number);
        }

        @Override
        public Element toElement() {
            return baseRing.valueOfBigInteger(number);
        }

        @Override
        public IOperand<Poly, Element> plus(IOperand<Poly, Element> oth) {
            if (oth instanceof Parser.NumberOperand) {
                return new NumberOperand(number.add(((NumberOperand) oth).number));
            } else
                return super.plus(oth);
        }

        @Override
        public IOperand<Poly, Element> minus(IOperand<Poly, Element> oth) {
            if (oth instanceof Parser.NumberOperand) {
                return new NumberOperand(number.subtract(((NumberOperand) oth).number));
            } else
                return super.minus(oth);
        }

        @Override
        public IOperand<Poly, Element> multiply(IOperand<Poly, Element> oth) {
            if (oth instanceof Parser.NumberOperand) {
                return new NumberOperand(number.multiply(((NumberOperand) oth).number));
            } else
                return oth.multiply(this);
        }

        @Override
        public IOperand<Poly, Element> divide(IOperand<Poly, Element> oth) {
            if (oth instanceof Parser.NumberOperand) {
                BigInteger[] divRem = this.number.divideAndRemainder(((NumberOperand) oth).number);
                if (divRem[1].isZero())
                    return new NumberOperand(divRem[0]);
                else
                    return super.divide(oth);
            } else
                return super.divide(oth);
        }

        @Override
        public IOperand<Poly, Element> pow(BigInteger exponent) {
            return new NumberOperand(Rings.Z.pow(number, exponent));
        }
    }

    /** A single variable */
    private final class VarOperand extends DefaultOperandOps implements IOperand<Poly, Element> {
        final int variable;

        VarOperand(int variable) { this.variable = variable; }

        @Override
        public Poly toPoly() { return polyRing.variable(variable); }

        @Override
        public IOperand<Poly, Element> multiply(IOperand<Poly, Element> oth) {
            if (oth instanceof Parser.NumberOperand) {
                return new MonomialOperand(polyRing.multiplyMutable(toPoly(), oth.toPoly()).lt());
            } else if (oth instanceof Parser.VarOperand) {
                int[] exponents = new int[polyRing.nVariables()];
                exponents[variable] += 1;
                exponents[((VarOperand) oth).variable] += 1;
                return new MonomialOperand(polyRing.factory().monomialAlgebra.create(exponents));
            }
            return oth.multiply(this);
        }

        @Override
        public IOperand<Poly, Element> pow(BigInteger exponent) {
            if (!exponent.isInt())
                return super.pow(exponent);
            int[] exponents = new int[polyRing.nVariables()];
            exponents[variable] += exponent.intValue();
            return new MonomialOperand(polyRing.factory().monomialAlgebra.create(exponents));
        }
    }

    /** A single monomial (x*y^2*z etc) */
    private class MonomialOperand extends DefaultOperandOps implements IOperand<Poly, Element> {
        Term term;

        MonomialOperand(Term term) {
            this.term = term;
        }

        @Override
        public Poly toPoly() {
            return polyRing.factory().create(term);
        }

        @Override
        public IOperand<Poly, Element> multiply(IOperand<Poly, Element> oth) {
            IMonomialAlgebra<Term> monomialAlgebra = polyRing.monomialAlgebra();
            if (oth instanceof Parser.NumberOperand) {
                return new MonomialOperand(monomialAlgebra.multiply(term, ((NumberOperand) oth).number));
            } else if (oth instanceof Parser.VarOperand) {
                int[] exponents = term.exponents;
                exponents[((VarOperand) oth).variable] += 1;
                return new MonomialOperand(term.forceSetDegreeVector(exponents, term.totalDegree + 1));
            } else if (oth instanceof Parser.MonomialOperand) {
                Term othTerm = ((MonomialOperand) oth).term;
                if (((long) othTerm.totalDegree) + term.totalDegree > Short.MAX_VALUE)
                    return super.multiply(oth);
                return new MonomialOperand(monomialAlgebra.multiply(term, othTerm));
            }
            return oth.multiply(this);
        }

        @Override
        public IOperand<Poly, Element> pow(BigInteger exponent) {
            if (exponent.isInt()) {
                int e = exponent.intValue();
                if (((long) term.totalDegree) * e > Short.MAX_VALUE)
                    return super.pow(exponent);
                int[] exponents = term.exponents;
                for (int i = 0; i < exponents.length; ++i)
                    exponents[i] *= e;
                return new MonomialOperand(term.forceSetDegreeVector(exponents, term.totalDegree * e));
            }
            return super.pow(exponent);
        }
    }

    /** A single polynomial */
    private class PolyOperand extends DefaultOperandOps implements IOperand<Poly, Element> {
        final Poly poly;

        PolyOperand(Poly poly) { this.poly = poly; }

        @Override
        public Poly toPoly() { return poly; }
    }

    /** Base ring element */
    private class ElementOperand extends DefaultOperandOps implements IOperand<Poly, Element> {
        final Element element;

        ElementOperand(Element element) {
            this.element = element;
        }

        @Override
        public Poly toPoly() {
            throw new UnsupportedOperationException();
        }

        @Override
        public Element toElement() {
            return element;
        }
    }

    /////////////////////////////////////////////////////Operators//////////////////////////////////////////////////////

    private enum Associativity {
        LeftToRight,
        RightToLeft
    }

    /** Operators */
    private enum Operator {
        // dummy ops
        BRACKET_OPEN(null, -1),
        BRACKET_CLOSE(null, -1),

        // priority = 2
        POWER(Associativity.LeftToRight, 20),

        // priority = 3
        UNARY_PLUS(Associativity.RightToLeft, 30),
        UNARY_MINUS(Associativity.RightToLeft, 30),

        // priority = 5
        MULTIPLY(Associativity.LeftToRight, 50),
        DIVIDE(Associativity.RightToLeft, 51),

        // priority = 6
        PLUS(Associativity.LeftToRight, 60),
        MINUS(Associativity.LeftToRight, 60);

        final Associativity associativity;
        final int priority;

        Operator(Associativity associativity, int priority) {
            this.associativity = associativity;
            this.priority = priority;
        }

        static boolean isPlusMinus(Operator op) {
            return op == PLUS || op == MINUS;
        }

        static Operator toUnaryPlusMinus(Operator op) {
            return op == PLUS ? UNARY_PLUS : op == MINUS ? UNARY_MINUS : null;
        }
    }

    /** convert token to operator */
    private static Operator tokenToOp(TokenType tType) {
        return tokenToOp[tType.ordinal()];
    }

    private static final Operator[] tokenToOp;

    static {
        tokenToOp = new Operator[TokenType.values().length];
        tokenToOp[TokenType.T_BRACKET_OPEN.ordinal()] = Operator.BRACKET_OPEN;
        tokenToOp[TokenType.T_BRACKET_CLOSE.ordinal()] = Operator.BRACKET_CLOSE;
        tokenToOp[TokenType.T_MULTIPLY.ordinal()] = Operator.MULTIPLY;
        tokenToOp[TokenType.T_DIVIDE.ordinal()] = Operator.DIVIDE;
        tokenToOp[TokenType.T_PLUS.ordinal()] = Operator.PLUS;
        tokenToOp[TokenType.T_MINUS.ordinal()] = Operator.MINUS;
        tokenToOp[TokenType.T_EXPONENT.ordinal()] = Operator.POWER;
    }

    private static boolean isInteger(String str) {
        if (str == null) {
            return false;
        }
        int length = str.length();
        if (length == 0) {
            return false;
        }
        int i = 0;
        if (str.charAt(0) == '-') {
            if (length == 1) {
                return false;
            }
            i = 1;
        }
        for (; i < length; i++) {
            char c = str.charAt(i);
            if (c < '0' || c > '9') {
                return false;
            }
        }
        return true;
    }
}
