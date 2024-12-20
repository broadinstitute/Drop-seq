package org.broadinstitute.dropseqrna.utils;

import java.util.Arrays;

public class ByteArrayWrapper implements Comparable<ByteArrayWrapper> {
    private final byte[] data;

    public ByteArrayWrapper(byte[] data) {
        this.data = data;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || getClass() != obj.getClass()) return false;
        ByteArrayWrapper that = (ByteArrayWrapper) obj;
        return Arrays.equals(this.data, that.data);
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(data);
    }

    @Override
    public int compareTo(ByteArrayWrapper other) {
        return Arrays.compare(this.data, other.data); // Lexicographical comparison
    }

    public byte[] getData() {
        return data;
    }

    public static void main(String[] args) {
        ByteArrayWrapper wrapper1 = new ByteArrayWrapper(new byte[]{1, 2, 3});
        ByteArrayWrapper wrapper2 = new ByteArrayWrapper(new byte[]{1, 2, 4});
        ByteArrayWrapper wrapper3 = new ByteArrayWrapper(new byte[]{1, 2, 3});

        System.out.println(wrapper1.compareTo(wrapper2)); // Output: -1 (wrapper1 < wrapper2)
        System.out.println(wrapper1.compareTo(wrapper3)); // Output: 0 (wrapper1 == wrapper3)
        System.out.println(wrapper2.compareTo(wrapper1)); // Output: 1 (wrapper2 > wrapper1)
    }
}
