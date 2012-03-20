/*
 * Created on Apr 3, 2009
 *
 */
package org.reactome.weka;

/**
 * This interface is used to check is a pair of proteins should be treated as
 * positive based on some feature. For example, if a protein pair is contained
 * by a protein-protein interaction data set, this protein pair should be treated
 * as positive. 
 * @author wgm
 *
 */
public interface PositiveChecker {
    /**
     * Check if a protein pair should be regarded as a positive pair based on
     * the undeath feature.
     * @param pair
     * @return
     */
    public boolean isPositive(String pair);
}