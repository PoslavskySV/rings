package cc.r2.core.poly;

import java.io.*;
import java.util.Base64;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Stanislav Poslavsky
 * @since 1.0
 */
public final class CommonUtils {
    private CommonUtils() {}

    public static void ensureFiniteFieldDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverFiniteField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureFieldDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverField())
                throw new IllegalArgumentException("Polynomial over finite field is expected; " + poly.getClass());
    }

    public static void ensureZDomain(IGeneralPolynomial... polys) {
        for (IGeneralPolynomial poly : polys)
            if (!poly.isOverZ())
                throw new IllegalArgumentException("Polynomial over Z is expected, but got " + poly.getClass());
    }

    public static String compress(Serializable object) {
        try(ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
            GZIPOutputStream zipped = new GZIPOutputStream(outBytes);
            ObjectOutputStream outSer = new ObjectOutputStream(zipped)) {
            outSer.writeObject(object);
            zipped.finish();
            return Base64.getEncoder().encodeToString(outBytes.toByteArray());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @SuppressWarnings("unchecked")
    public static <T> T uncompress(String object) {
        byte[] decoded = Base64.getDecoder().decode(object);
        try(ByteArrayInputStream inBytes = new ByteArrayInputStream(decoded);
            ObjectInputStream inSer = new ObjectInputStream(new GZIPInputStream(inBytes))) {
            return (T) inSer.readObject();
        } catch (IOException|ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
