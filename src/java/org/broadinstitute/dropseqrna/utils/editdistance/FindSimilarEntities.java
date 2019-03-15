package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.List;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public interface FindSimilarEntities<T extends Comparable<T>, M> {
	public FindSimilarEntitiesResult<T,M> find (final T entity, final List<T> searchSpace, final ObjectCounter<T> counts);	
}
