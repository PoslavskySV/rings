package cc.redberry.rings.io;

import cc.redberry.rings.Rational;
import cc.redberry.rings.Rationals;
import cc.redberry.rings.Ring;
import cc.redberry.rings.Rings;
import cc.redberry.rings.bigint.BigInteger;
import cc.redberry.rings.io.Tokenizer.Token;
import cc.redberry.rings.io.Tokenizer.TokenType;
import cc.redberry.rings.poly.IPolynomial;
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

import static cc.redberry.rings.io.Tokenizer.END;
import static cc.redberry.rings.io.Tokenizer.TokenType.*;

/**
 * High-level parser and stringifier of ring elements. Uses shunting-yard algorithm for parsing.
 *
 * @param <Element> type of resulting elements
 * @param <Term>    underlying polynomial terms
 * @param <Poly>    underlying multivariate polynomials
 */
public class Coder<
        Element,
        Term extends AMonomial<Term>,
        Poly extends AMultivariatePolynomial<Term, Poly>>
        implements IParser<Element>, IStringifier<Element> {
    // parser stuff
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

    // stringifier stuff
    /** toString bindings */
    private final Map<Element, String> bindings;
    /** inner coders */
    private final Map<Ring<?>, Coder<?, ?, ?>> subcoders;

    private Coder(Ring<Element> baseRing,
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

        this.bindings = eVariables.entrySet()
                .stream()
                .collect(Collectors.toMap(Map.Entry::getValue, Map.Entry::getKey));
        this.subcoders = new HashMap<>();
    }

    /** Add string -> element mapping */
    public Coder<Element, Term, Poly> bind(String var, Element el) {
        bindings.put(el, var);
        eVariables.put(var, el);
        return this;
    }

    /** Add stringifier of inner elements */
    public Coder<Element, Term, Poly> withEncoder(Coder<?, ?, ?> subencoder) {
        subcoders.put(subencoder.baseRing, subencoder);
        return this;
    }

    @Override
    @SuppressWarnings("unchecked")
    public <K> IStringifier<K> substringifier(Ring<K> ring) {
        IStringifier<K> s = (IStringifier<K>) subcoders.get(ring);
        if (s == null)
            return IStringifier.dummy();
        else
            return s;
    }

    @Override
    public Map<Element, String> getBindings() {
        return bindings;
    }

    /**
     * Decode element from its string representation (#parse)
     */
    public Element decode(String string) {
        return parse(string);
    }

    /**
     * Encode element to its string representation (#stringify)
     */
    @SuppressWarnings("unchecked")
    public String encode(Element element) {
        return stringify(element);
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
            Poly extends AMultivariatePolynomial<Term, Poly>> Coder<Element, Term, Poly>
    mkCoder(Ring<Element> baseRing,
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
        return new Coder<>(baseRing, eVariables, polyRing, iVariables, polyToElement);
    }

    /**
     * Create coder for generic ring
     *
     * @param ring the ring
     */
    public static <E> Coder<E, ?, ?> mkCoder(Ring<E> ring) {
        return mkCoder(ring, new HashMap<>());
    }

    /**
     * Create coder for generic rings
     *
     * @param ring      the ring
     * @param variables map string_variable -> ring_element
     */
    @SuppressWarnings("unchecked")
    public static <E> Coder<E, ?, ?> mkCoder(Ring<E> ring,
                                             Map<String, E> variables) {
        if (ring instanceof MultivariateRing)
            return mkMultivariateCoder((MultivariateRing) ring, (Map) variables);

        if (ring instanceof IPolynomialRing && ((IPolynomialRing) ring).nVariables() == 1)
            return mkUnivariateCoder((IPolynomialRing) ring, (Map) variables);

        return new Coder<>(ring, variables, null, null, null);
    }

    /**
     * Create coder for multivariate polynomial rings
     *
     * @param ring      the ring
     * @param variables map string_variable -> ring_element
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Coder<Poly, Term, Poly> mkMultivariateCoder(MultivariateRing<Poly> ring,
                                                Map<String, Poly> variables) {
        return mkCoder(ring, variables, ring, variables, Function.identity());
    }

    /**
     * Create coder for multivariate polynomial rings
     *
     * @param ring      the ring
     * @param variables polynomial variables
     */
    public static <Term extends AMonomial<Term>, Poly extends AMultivariatePolynomial<Term, Poly>>
    Coder<Poly, Term, Poly> mkMultivariateCoder(MultivariateRing<Poly> ring,
                                                String... variables) {
        Map<String, Poly> pVariables = new HashMap<>();
        for (int i = 0; i < variables.length; ++i)
            pVariables.put(variables[i], ring.variable(i));
        return mkMultivariateCoder(ring, pVariables);
    }

    /**
     * Create coder for multivariate polynomial rings
     *
     * @param ring      the ring
     * @param cfCoder   coder for coefficient ring elements
     * @param variables map string_variable -> ring_element
     */
    public static <E> Coder<MultivariatePolynomial<E>, Monomial<E>, MultivariatePolynomial<E>>
    mkMultivariateCoder(MultivariateRing<MultivariatePolynomial<E>> ring,
                        Coder<E, ?, ?> cfCoder,
                        Map<String, MultivariatePolynomial<E>> variables) {
        cfCoder.eVariables.forEach((k, v) -> variables.put(k, ring.factory().createConstant(v)));
        return mkMultivariateCoder(ring, variables).withEncoder(cfCoder);
    }

    /**
     * Create parser for multivariate polynomial rings
     *
     * @param ring      the ring
     * @param cfCoder   coder of coefficient ring elements
     * @param variables polynomial variables
     */
    public static <E> Coder<MultivariatePolynomial<E>, Monomial<E>, MultivariatePolynomial<E>>
    mkMultivariateCoder(MultivariateRing<MultivariatePolynomial<E>> ring,
                        Coder<E, ?, ?> cfCoder,
                        String... variables) {
        Map<String, MultivariatePolynomial<E>> eVariables = new HashMap<>();
        for (int i = 0; i < variables.length; ++i)
            eVariables.put(variables[i], ring.variable(i));
        return mkMultivariateCoder(ring, cfCoder, eVariables);
    }


    /**
     * Create coder for univariate polynomial rings
     *
     * @param ring      the ring
     * @param variables map string_variable -> ring_element
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IUnivariatePolynomial<Poly>>
    Coder<Poly, ?, ?> mkUnivariateCoder(IPolynomialRing<Poly> ring,
                                        Map<String, Poly> variables) {
        MultivariateRing mRing = Rings.MultivariateRing(ring.factory().asMultivariate());
        Map<String, AMultivariatePolynomial> pVariables = variables.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, e -> e.getValue().asMultivariate()));
        return mkCoder(ring, variables, mRing, pVariables, p -> (Poly) p.asUnivariate());
    }

    /**
     * Create coder for univariate polynomial rings
     *
     * @param ring     the ring
     * @param variable variable string
     */
    public static <Poly extends IUnivariatePolynomial<Poly>>
    Coder<Poly, ?, ?> mkUnivariateCoder(IPolynomialRing<Poly> ring,
                                        String variable) {
        return mkUnivariateCoder(ring, new HashMap<String, Poly>() {{put(variable, ring.variable(0));}});
    }

    /**
     * Create coder for univariate polynomial rings
     *
     * @param ring      the ring
     * @param cfCoder   coder of coefficient ring elements
     * @param variables map string_variable -> ring_element
     */
    public static <E> Coder<UnivariatePolynomial<E>, ?, ?>
    mkUnivariateCoder(IPolynomialRing<UnivariatePolynomial<E>> ring,
                      Coder<E, ?, ?> cfCoder,
                      Map<String, UnivariatePolynomial<E>> variables) {
        cfCoder.eVariables.forEach((k, v) -> variables.put(k, ring.factory().createConstant(v)));
        return mkUnivariateCoder(ring, variables).withEncoder(cfCoder);
    }

    /**
     * Create coder for univariate polynomial rings
     *
     * @param ring     the ring
     * @param cfCoder  coder of coefficient ring elements
     * @param variable string variable
     */
    public static <E> Coder<UnivariatePolynomial<E>, ?, ?>
    mkUnivariateCoder(IPolynomialRing<UnivariatePolynomial<E>> ring,
                      Coder<E, ?, ?> cfCoder,
                      String variable) {
        HashMap<String, UnivariatePolynomial<E>> eVariables = new HashMap<>();
        eVariables.put(variable, ring.variable(0));
        return mkUnivariateCoder(ring, cfCoder, eVariables);
    }

    /**
     * Create coder for polynomial ring
     *
     * @param ring      the ring
     * @param variables string variables
     */
    @SuppressWarnings("unchecked")
    public static <Poly extends IPolynomial<Poly>> Coder<Poly, ?, ?>
    mkPolynomialCoder(IPolynomialRing<Poly> ring,
                      String... variables) {
        if (ring instanceof MultivariateRing)
            return mkMultivariateCoder((MultivariateRing) ring, variables);
        else
            return mkUnivariateCoder((IPolynomialRing) ring, variables[0]);
    }

    /**
     * Create coder for rational elements
     */
    @SuppressWarnings("unchecked")
    public static <E> Coder<Rational<E>, ?, ?>
    mkRationalsCoder(Rationals<E> ring,
                     Coder<E, ?, ?> elementsCoder) {
        return mkNestedCoder(ring, new HashMap<>(), elementsCoder, e -> new Rational<>(ring.ring, e));
    }

    /**
     * Create coder for nested rings (e.g. fractions over polynomials etc).
     *
     * <p>Example:
     * <pre><code>
     * // GF(17, 3) as polynomials over "t"
     * FiniteField<UnivariatePolynomialZp64> gf = Rings.GF(17, 3);
     * // parser of univariate polynomials over "t" from GF(17, 3)
     * Coder<UnivariatePolynomialZp64, ?, ?> gfParser = Coder.mkUnivariateCoder(gf, mkVars(gf, "t"));
     *
     * // ring GF(17, 3)[x, y, z]
     * MultivariateRing<MultivariatePolynomial<UnivariatePolynomialZp64>> polyRing = Rings.MultivariateRing(3, gf);
     * // parser of multivariate polynomials over GF(17, 3)
     * Coder<MultivariatePolynomial<UnivariatePolynomialZp64>, ?, ?> polyParser = Coder.mkMultivariateCoder(polyRing, gfParser, "x", "y", "z");
     *
     * // field Frac(GF(17, 3)[x, y, z])
     * Rationals<MultivariatePolynomial<UnivariatePolynomialZp64>> fracRing = Rings.Frac(polyRing);
     * // parser of elements in Frac(GF(17, 3)[x, y, z])
     * Coder<Rational<MultivariatePolynomial<UnivariatePolynomialZp64>>, ?, ?> fracParser =
     *              Coder.mkNestedCoder(
     *                   fracRing,                        // the frac field
     *                   new HashMap<>(),                 // variables (no any)
     *                   polyParser,                      // parser of multivariate polynomials
     *                   p -> new Rational<>(polyRing, p) // convert from polys to fractions
     *              );
     * </code></pre>
     *
     * @param ring       the ring
     * @param variables  map string_variable -> ring_element
     * @param innerCoder coder for underlying ring elements
     * @param imageFunc  mapping from @{code I} to @{code E}
     */
    @SuppressWarnings("unchecked")
    public static <E, I>
    Coder<E, ?, ?> mkNestedCoder(Ring<E> ring,
                                 Map<String, E> variables,
                                 Coder<I, ?, ?> innerCoder,
                                 Function<I, E> imageFunc) {
        if (ring instanceof MultivariateRing && ((MultivariateRing) ring).factory() instanceof MultivariatePolynomial)
            return mkMultivariateCoder((MultivariateRing) ring, innerCoder, (Map) variables);
        else if (ring instanceof UnivariateRing && ((UnivariateRing) ring).factory() instanceof UnivariatePolynomial)
            return mkUnivariateCoder((UnivariateRing) ring, innerCoder, (Map) variables);

        innerCoder.eVariables.forEach((k, v) -> variables.put(k, imageFunc.apply(v)));
        return new Coder(
                ring,
                variables,
                innerCoder.polyRing,
                innerCoder.pVariables,
                innerCoder.polyToElement == null ? null : innerCoder.polyToElement.andThen(imageFunc))
                .withEncoder(innerCoder);
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
                if (Operator.isPlusMinus(op) && Operator.isPlusMinus(tokenToOp(previousToken))) {
                    if (op == Operator.PLUS && previousToken == T_MINUS) {
                        // manage --
                        tType = T_MINUS;
                        op = Operator.MINUS;
                    } else if (op == Operator.MINUS && previousToken == T_MINUS) {
                        // manage -+
                        tType = T_PLUS;
                        op = Operator.PLUS;
                    }
                    operands.push(Zero); // push neutral element
                }

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
                if (!(right instanceof Coder.NumberOperand))
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
            return this instanceof Coder.ElementOperand;
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
            if (oth instanceof Coder.NumberOperand) {
                return new NumberOperand(number.add(((NumberOperand) oth).number));
            } else
                return super.plus(oth);
        }

        @Override
        public IOperand<Poly, Element> minus(IOperand<Poly, Element> oth) {
            if (oth instanceof Coder.NumberOperand) {
                return new NumberOperand(number.subtract(((NumberOperand) oth).number));
            } else
                return super.minus(oth);
        }

        @Override
        public IOperand<Poly, Element> multiply(IOperand<Poly, Element> oth) {
            if (oth instanceof Coder.NumberOperand) {
                return new NumberOperand(number.multiply(((NumberOperand) oth).number));
            } else
                return oth.multiply(this);
        }

        @Override
        public IOperand<Poly, Element> divide(IOperand<Poly, Element> oth) {
            if (oth instanceof Coder.NumberOperand) {
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
            if (oth instanceof Coder.NumberOperand) {
                return new MonomialOperand(polyRing.multiplyMutable(toPoly(), oth.toPoly()).lt());
            } else if (oth instanceof Coder.VarOperand) {
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
            if (oth instanceof Coder.NumberOperand) {
                return new MonomialOperand(monomialAlgebra.multiply(term, ((NumberOperand) oth).number));
            } else if (oth instanceof Coder.VarOperand) {
                int[] exponents = term.exponents;
                exponents[((VarOperand) oth).variable] += 1;
                return new MonomialOperand(term.forceSetDegreeVector(exponents, term.totalDegree + 1));
            } else if (oth instanceof Coder.MonomialOperand) {
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
        DIVIDE(Associativity.LeftToRight, 50),

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
            return PLUS == op || MINUS == op;
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
