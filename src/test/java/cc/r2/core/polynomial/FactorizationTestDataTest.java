package cc.r2.core.polynomial;

import org.junit.Assert;
import org.junit.Test;

import java.io.InputStream;

public class FactorizationTestDataTest {

    @Test
    public void test() throws Exception {
        FactorizationTestData fct = FactorizationTestData.decode("181|1,172,85,84,151,122,53,107,117,82,152,133,151,178,1|2,55,1|1,51,166,130,172,1|1,8,89,23,94,1|1,42,110,164,1");
        System.out.println(fct.factorization);
    }

}