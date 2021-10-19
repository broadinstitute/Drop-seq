package org.broadinstitute.dropseqrna.sbarro;

import java.io.Closeable;
import java.io.IOException;
import java.util.Iterator;

public class IndexedIterator<TItem, TIterable extends Iterable<TItem> & Closeable>
        implements Iterator<IndexedItem<TItem>>, Iterable<IndexedItem<TItem>>, Closeable {
    private long index = 0;
    private final long limit;
    private final TIterable inner;
    private final Iterator<TItem> iterator;

    public IndexedIterator(final TIterable inner, final long skip, final long count) {
        this.inner = inner;
        this.iterator = inner.iterator();
        this.limit = count < 0 ? -1 : skip + count;

        while (this.index < skip && this.iterator.hasNext()) {
            this.iterator.next();
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
