package org.broadinstitute.dropseqrna.utils;

import org.testng.annotations.Test;
import static org.testng.Assert.*;

public class ByteArrayWrapperTest {

    @Test
    public void testEquals() {
        byte[] data1 = {1, 2, 3};
        byte[] data2 = {1, 2, 3};
        byte[] data3 = {4, 5, 6};

        ByteArrayWrapper wrapper1 = new ByteArrayWrapper(data1);
        ByteArrayWrapper wrapper2 = new ByteArrayWrapper(data2);
        ByteArrayWrapper wrapper3 = new ByteArrayWrapper(data3);

        assertEquals(wrapper1, wrapper2);
        assertNotEquals(wrapper1, wrapper3);
    }

    @Test
    public void testHashCode() {
        byte[] data1 = {1, 2, 3};
        byte[] data2 = {1, 2, 3};
        byte[] data3 = {4, 5, 6};

        ByteArrayWrapper wrapper1 = new ByteArrayWrapper(data1);
        ByteArrayWrapper wrapper2 = new ByteArrayWrapper(data2);
        ByteArrayWrapper wrapper3 = new ByteArrayWrapper(data3);

        assertEquals(wrapper1.hashCode(), wrapper2.hashCode());
        assertNotEquals(wrapper1.hashCode(), wrapper3.hashCode());
    }

    @Test
    public void testCompareTo() {
        byte[] data1 = {1, 2, 3};
        byte[] data2 = {1, 2, 3};
        byte[] data3 = {4, 5, 6};
        byte[] data4 = {1, 2, 4};

        ByteArrayWrapper wrapper1 = new ByteArrayWrapper(data1);
        ByteArrayWrapper wrapper2 = new ByteArrayWrapper(data2);
        ByteArrayWrapper wrapper3 = new ByteArrayWrapper(data3);
        ByteArrayWrapper wrapper4 = new ByteArrayWrapper(data4);

        assertEquals(wrapper1.compareTo(wrapper2), 0);
        assertTrue(wrapper1.compareTo(wrapper3) < 0);
        assertTrue(wrapper3.compareTo(wrapper1) > 0);
        assertTrue(wrapper1.compareTo(wrapper4) < 0);
    }

    @Test
    public void testGetData() {
        byte[] data = {1, 2, 3};
        ByteArrayWrapper wrapper = new ByteArrayWrapper(data);
        assertEquals(wrapper.getData(), data);
    }
}