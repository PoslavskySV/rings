package cc.redberry.rings.parser;

/**
 *
 */
public final class Tokenizer {
    private final CharacterStream stream;

    private Tokenizer(CharacterStream stream) {
        this.stream = stream;
    }

    /** token type */
    public enum TokenType {
        T_VARIABLE,
        T_PLUS,
        T_MINUS,
        T_MULTIPLY,
        T_DIVIDE,
        T_EXPONENT,
        T_SPACE,
        T_BRACKET_OPEN,
        T_BRACKET_CLOSE,
        T_END
    }

    // whether there is a single char in local buffer
    private boolean bufferedCharDefined;
    // local buffer with only one character
    private char bufferedChar;

    /** has more elements to consider */
    private boolean hasNextChar() {
        return bufferedCharDefined || stream.hasNext();
    }

    /** next character */
    private char nextChar() {
        char c;
        if (bufferedCharDefined) {
            bufferedCharDefined = false;
            c = bufferedChar;
        } else
            c = stream.next();

        checkChar(c);
        return c;
    }

    private static void checkChar(char c) {
        switch (c) {
            case '=':
            case '&':
            case '%':
            case '!':
            case '~':
                throw new IllegalArgumentException(String.format("Illegal character %s", c));
        }
    }

    // envelope tokens wit ( )
    private boolean first = true, last = true;

    private Token firstToken() {
        first = false;
        return BRACKET_OPEN;
    }

    private Token lastToken() {
        if (!last)
            return END;
        last = false;
        return BRACKET_CLOSE;
    }

    /**
     * Get the next token from stream
     */
    public Token nextToken() {
        if (first)
            return firstToken();

        if (!hasNextChar())
            return lastToken();

        // peek the first char from the buffer
        char seed = nextChar();
        // check whether this is already a token
        Token token = primitiveToken(seed);

        // get rid of spaces
        if (token == SPACE) {
            do {
                if (!hasNextChar())
                    return lastToken();

                seed = nextChar();
                // check whether this is already a token
                token = primitiveToken(seed);
            } while (token == SPACE);
        }

        if (token != null)
            return token;

        // <-seed char is a start of variable

        // string that contains the beginning of the variable
        String sBegin = stream.currentString();
        int pBegin = stream.indexInCurrentString();

        if (!hasNextChar()) {
            // seed char is the last
            return new Token(TokenType.T_VARIABLE, sBegin.substring(pBegin, pBegin + 1));
        }

        char c = 0;
        int nSpaces = 0;
        boolean needBuffer = false;
        while (hasNextChar()) {
            c = stream.next();
            Token t = primitiveToken(c);
            if (t == SPACE) {// terminating spaces
                ++nSpaces;
                continue;
            }
            if (t == null) {
                if (nSpaces > 0)
                    throw new IllegalArgumentException("spaces in variable name are forbidden");
                continue;
            }
            needBuffer = true;
            break;
        }

        // string that contains the end of the variable
        String sEnd = stream.currentString();
        int pEnd = stream.indexInCurrentString() + 1; // end inclusive

        // put the next character back into the buffer
        if (needBuffer) {
            bufferedChar = c;
            bufferedCharDefined = true;
            pEnd -= 1;
        }

        String variable;
        if (sBegin == sEnd)
            variable = sBegin.substring(pBegin, pEnd - nSpaces);
        else {
            if (nSpaces < pEnd)
                variable = sBegin.substring(pBegin) + sEnd.substring(0, pEnd);
            else
                variable = sBegin.substring(pBegin, sBegin.length() - (pEnd - nSpaces));
        }

        return new Token(TokenType.T_VARIABLE, variable);
    }

    private static Token primitiveToken(char character) {
        switch (character) {
            case '+': return PLUS;
            case '-': return MINUS;
            case '*': return MULTIPLY;
            case '/': return DIVIDE;
            case '^': return EXPONENT;
            case '(': return BRACKET_OPEN;
            case ')': return BRACKET_CLOSE;
            case ' ': return SPACE;
            default: return null;
        }
    }

    public static final Token
            END = new Token(TokenType.T_END, ""),
            PLUS = new Token(TokenType.T_PLUS, "+"),
            MINUS = new Token(TokenType.T_MINUS, "-"),
            MULTIPLY = new Token(TokenType.T_MULTIPLY, "*"),
            DIVIDE = new Token(TokenType.T_DIVIDE, "/"),
            EXPONENT = new Token(TokenType.T_EXPONENT, "^"),
            BRACKET_OPEN = new Token(TokenType.T_BRACKET_OPEN, "("),
            BRACKET_CLOSE = new Token(TokenType.T_BRACKET_CLOSE, ")"),
            SPACE = new Token(TokenType.T_SPACE, " ");

    /** Simple token */
    public static final class Token {
        public final TokenType tokenType;
        public final String content;

        public Token(TokenType tokenType, String content) {
            this.tokenType = tokenType;
            this.content = content;
        }

        @Override
        public String toString() {
            return tokenType + "(" + content + ")";
        }
    }

    public interface CharacterStream {
        boolean hasNext();

        char next();

        String currentString();

        int indexInCurrentString();
    }

    public static CharacterStream mkCharacterStream(String string) {
        return new CharacterStream() {
            int index = 0;
            int currentIndex = 0;

            @Override
            public boolean hasNext() {
                return index < string.length();
            }

            @Override
            public char next() {
                currentIndex = index;
                return string.charAt(index++);
            }

            @Override
            public String currentString() {
                return string;
            }

            @Override
            public int indexInCurrentString() {
                return currentIndex;
            }
        };
    }

    public static Tokenizer mkTokenizer(String string) {
        return new Tokenizer(mkCharacterStream(string));
    }
}
