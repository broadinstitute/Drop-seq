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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableOnceIterator;

import java.util.Iterator;

/**
 * Given an iterator of some type, transform each record to another domain object, possibly of the same or a different type.
 * These are 1 to 1 transforms.
 * @author nemesh
 *
 * @param <INPUT> The input iterator object type
 * @param <OUTPUT> The output iterator object type
 */
public abstract class TransformingIterator<INPUT,OUTPUT> extends IterableOnceIterator<OUTPUT> {

	protected final Iterator<INPUT> underlyingIterator;

	public TransformingIterator (final Iterator <INPUT> underlyingIterator) {
		this.underlyingIterator=underlyingIterator;
	}

	@Override
	public boolean hasNext() {
		return this.underlyingIterator.hasNext();
	}

	@Override
	public abstract OUTPUT next();

	@Override
	public void remove() {
		this.underlyingIterator.remove();
	}

	@Override
	public void close() {
		CloserUtil.close(this.underlyingIterator);
		// TODO Auto-generated method stub

	}

}
