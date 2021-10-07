package org.broadinstitute.dropseqrna.sbarro;

@SuppressWarnings("ClassCanBeRecord")
public class IndexedItem<T> {
    private final long index;
    private final T item;

    public IndexedItem(long index, T item) {
        this.index = index;
        this.item = item;
    }

    public long getIndex() {
        return index;
    }

    public T getItem() {
        return item;
    }
}
