/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.utils;

import java.util.Iterator;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IterableOnceIterator;
import htsjdk.samtools.util.PeekableIterator;

public abstract class FilteredIterator<T>
        extends IterableOnceIterator<T>
        implements CloseableIterator<T>,
        PassFailTrackingIteratorI
{

    final PeekableIterator<T> it;
    boolean firstTime = true;
    private final ObjectSink<T> sink;
    private long recordsPass=0;
	private long recordsFail=0;
	
        
    /**
     * If instantiating your FilteredIterator with an ObjectSink, reads that fail the filter will be put in the sink.
     * @param underlyingIterator The iterator.
     * @param filteredReadSink An ObjectSink to place filtered reads in.
     */
    protected FilteredIterator(final Iterator<T> underlyingIterator, final ObjectSink<T> filteredReadSink) {
        it = new PeekableIterator<>(underlyingIterator);
        this.sink=filteredReadSink;
    }
    
    protected FilteredIterator(final Iterator<T> underlyingIterator) {
    	this(underlyingIterator, null);        
    }

    /**
     * @return true if record should be skipped
     */
    public abstract boolean filterOut(final T rec);
    
    /**
     * After the iterator is finished, the subclass logs the number of pass/fail records.  
     * Can be implemented to emit logging messages, but does not log by default
     */
    public void logFilterResults() {
    	// This is a NO-OP unless implemented by the subclass
    }

    // if there's a sink, then filtered reads are added to the sink.
    // 
	private void skipUndesired() {		
        while (it.hasNext()) {
        	T rec = it.peek();
        	boolean filter = filterOut(rec);
        	if (filter) 
        		this.recordsFail++;
        	if (!filter) {
        		this.recordsPass++;
        		return; // break out of skipping records.
        	}
        	// you're filtered, sink if needed.
        	if (this.sink!=null) 
        		sink.add(rec);        	
        	// you're filtered, get the next record.
        	it.next();        	        	
        }			
    }

    private void maybeSkipFirstTime() {
        if (firstTime) {
            skipUndesired();
            firstTime = false;
        }
    }

    @Override
    public void close() {
        it.close();
    }

    @Override
    public boolean hasNext() {
        maybeSkipFirstTime();
        boolean hasNext=it.hasNext();
        if (!hasNext) 
        	logFilterResults();
        return it.hasNext();
    }

    @Override
    public T next() {
        maybeSkipFirstTime();
        final T ret = it.next();
        skipUndesired();
        return ret;
    }

    @Override
    public void remove() {
        it.remove();
    }
    
    /**
     * If you passed in an ObjectSink to filter reads, then you can get it back here.
     * @return
     */
    public ObjectSink<T> getFilteredReadSink () {
    	return this.sink;
    }
    
    @Override
	public long getRecordsPassed() {
		return this.recordsPass;		
	}

	@Override
	public long getRecordsFailed() {
		return this.recordsFail;
	}
}
