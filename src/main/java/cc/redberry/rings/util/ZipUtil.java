package cc.redberry.rings.util;

import java.io.*;
import java.util.Base64;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @since 1.0
 */
public final class ZipUtil {
    private ZipUtil() {}

    /**
     * Compress object to a string
     */
    public static String compress(Serializable object) {
        try (ByteArrayOutputStream outBytes = new ByteArrayOutputStream();
             GZIPOutputStream zipped = new GZIPOutputStream(outBytes);
             ObjectOutputStream outSer = new ObjectOutputStream(zipped)) {
            outSer.writeObject(object);
            zipped.finish();
            return Base64.getEncoder().encodeToString(outBytes.toByteArray());
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Decompress object from its string code obtained via {@link ZipUtil#compress(Serializable)}
     */
    @SuppressWarnings("unchecked")
    public static <T> T uncompress(String object) {
        byte[] decoded = Base64.getDecoder().decode(object);
        try (ByteArrayInputStream inBytes = new ByteArrayInputStream(decoded);
             ObjectInputStream inSer = new ObjectInputStream(new GZIPInputStream(inBytes))) {
            return (T) inSer.readObject();
        } catch (IOException | ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }
}
