package cc.redberry.rings.io;

import java.util.ArrayList;
import java.util.List;

/**
 * Simple math expression tokenizer
 */
public final class Tokenizer {
    private final CharacterStream stream;

    /**
     * Create tokenizer of a given char stream
     */
    public Tokenizer(CharacterStream stream) {
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
        T_NEWLINE,
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
            case '\t':
            case '\n':
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

        Token(TokenType tokenType, String content) {
            this.tokenType = tokenType;
            this.content = content;
        }

        @Override
        public String toString() {
            return tokenType + "(" + content + ")";
        }
    }

    /** Stream of chars. Implementations are not synchronized and doesn't support concurrent iteration. */
    public interface CharacterStream {
        /** next char available in this stream */
        boolean hasNext();

        /** next char from this stream */
        char next();

        /** string containing current char */
        String currentString();

        /** index of char in string */
        int indexInCurrentString();

        /**
         * skip all chars preceding the specified char and place caret to the first char after the specified one
         *
         * @return false if stream finished without specified char and true otherwise
         */
        default boolean seek(char c) {
            while (hasNext()) {
                char n = next();
                if (n == c)
                    return true;
            }
            return false;
        }
    }

    private static final class ConcatStreams implements CharacterStream {
        final List<CharacterStream> streams;
        int currentStream;

        ConcatStreams(List<CharacterStream> streams) {
            this(streams, 0);
        }

        ConcatStreams(List<CharacterStream> streams, int currentStream) {
            this.streams = streams;
            this.currentStream = currentStream;
        }

        @Override
        public boolean hasNext() {
            int cs = currentStream;
            for (; cs < streams.size(); ++cs)
                if (streams.get(cs).hasNext())
                    return true;
            return false;
        }

        @Override
        public char next() {
            int cs = currentStream;
            for (; cs < streams.size(); ++cs)
                if (streams.get(cs).hasNext()) {
                    currentStream = cs;
                    return streams.get(cs).next();
                }
            throw new IllegalArgumentException("no more elements in this stream");
        }

        @Override
        public String currentString() {
            return streams.get(currentStream).currentString();
        }

        @Override
        public int indexInCurrentString() {
            return streams.get(currentStream).indexInCurrentString();
        }
    }

    /** Concat char streams */
    public static CharacterStream concat(CharacterStream a, CharacterStream b) {
        List<CharacterStream> streams = new ArrayList<>();
        streams.add(a);
        streams.add(b);
        return new ConcatStreams(streams);
    }

    /**
     * Create character stream from string
     *
     * @param terminateChar if a non-null value specified, stream will terminate on the last char preceding the {@code terminateChar}
     */
    public static CharacterStream mkCharacterStream(String string, Character terminateChar) {
        return new CharacterStream() {
            int index = 0;
            int currentIndex = 0;

            @Override
            public boolean hasNext() {
                return index < string.length()
                        && (terminateChar == null || string.charAt(index) != terminateChar);
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

    /**
     * Create string tokenizer
     *
     * @param terminateChar if a non-null value specified, stream will terminate on the last char preceding the {@code terminateChar}
     */
    public static Tokenizer mkTokenizer(String string, Character terminateChar) {
        return new Tokenizer(mkCharacterStream(string, terminateChar));
    }

    /**
     * Create string tokenizer
     */
    public static Tokenizer mkTokenizer(String string) {
        return new Tokenizer(mkCharacterStream(string, null));
    }
}
