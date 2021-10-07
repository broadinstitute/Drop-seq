package org.broadinstitute.dropseqrna.sbarro;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

public class IndexedIterator<TItem, TIterator extends Iterator<TItem> & Iterable<TItem> & Closeable>
        implements Iterator<IndexedItem<TItem>>, Iterable<IndexedItem<TItem>>, Closeable {
    private long index = 0;
    private final long limit;
    private final TIterator inner;

    public IndexedIterator(final TIterator inner, final long skip, final long count) {
        this.inner = inner;
        this.limit = count < 0 ? -1 : skip + count;

        while (this.index < skip && this.inner.hasNext()) {
            this.inner.next();
            this.index++;
        }
    }

    @Override
    public void close() throws IOException {
        this.inner.close();
    }

    @Override
    public Iterator<IndexedItem<TItem>> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        if (this.limit < 0 || this.limit > this.index) {
            return this.inner.hasNext();
        } else {
            return false;
        }
    }

    @Override
    public IndexedItem<TItem> next() {
        final IndexedItem<TItem> next = new IndexedItem<>(this.index, this.inner.next());
        this.index++;
        return next;
    }
}
