package org.broadinstitute.dropseqrna.sbarro;

import java.util.Iterator;

public class IndexedIterator<TItem, TIterator extends Iterator<TItem>>
        implements Iterator<IndexedItem<TItem>>, Iterable<IndexedItem<TItem>> {
    private long index = 0;
    private final long limit;
    private final Iterator<TItem> iterator;

    public IndexedIterator(final TIterator inner, final long skip, final long count) {
        this.iterator = inner;
        this.limit = count < 0 ? -1 : skip + count;

        while (this.index < skip && this.iterator.hasNext()) {
            this.iterator.next();
            this.index++;
        }
    }

    @Override
    public Iterator<IndexedItem<TItem>> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        if (this.limit < 0 || this.limit > this.index) {
            return this.iterator.hasNext();
        } else {
            return false;
        }
    }

    @Override
    public IndexedItem<TItem> next() {
        final IndexedItem<TItem> next = new IndexedItem<>(this.index, this.iterator.next());
        this.index++;
        return next;
    }
}
