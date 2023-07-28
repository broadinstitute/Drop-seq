/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
import java.util.function.Predicate;

import org.broadinstitute.dropseqrna.utils.readiterators.TagValueFilteringIterator;

import htsjdk.samtools.util.Log;

/**
 * FilteredIterator using java.util.function.Predicate
 * Record is included if predicate.test() == true
 */
public class PredicateFilteredIterator<T>
        extends FilteredIterator<T> {
    private final Predicate<T> predicate;
    
	
    /**
     *
     * @param underlyingIterator records to be iterated over
     * @param predicate test() method means include this record
     */
    public PredicateFilteredIterator(Iterator<T> underlyingIterator, Predicate<T> predicate) {
        super(underlyingIterator);
        this.predicate = predicate;
    }

    @Override
    public boolean filterOut(T rec) {
    	return !predicate.test(rec);        
    }

	@Override
	public void logFilterResults() {
		Class<? extends PredicateFilteredIterator> subclass = this.getClass();
		Log log = Log.getInstance(subclass);
		String msg = String.format("Records pass [%d] records fail [%d] ", this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);		
	}


}
