package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.List;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

/**
 * Interface to plug in various collapse strategies for greedy collapse algorithm in MapBarcodesByEditDistance.collapseBarcodesGeneric
 * @author nemesh
 *
 * @param <T> The object type to be compared for collapsing.  
 * @param <M> A metrics-like reporting object if additional reporting on collapse occurs.
 */
public interface FindSimilarEntities<T extends Comparable<T>, M> {
	/**
	 * For a given entity and a search space of other entities, find the subset that are defined as similar by the implementor of this interface.
	 * @param entity The parent entity
	 * @param searchSpace Other entities, a subset of which may be similar to the parent entity
	 * @param counts The number of observations of entities.  This defines an ordering for entities collapse, from largest to smallest. 
	 * @return The result for this similarity search.  Contains a mapping from the entity to similar other entities in the search space, as well as metrics.
	 */
	public FindSimilarEntitiesResult<T,M> find (final T entity, final List<T> searchSpace, final ObjectCounter<T> counts);	
}
