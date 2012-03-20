/*
 * Created on Sep 28, 2006
 *
 */
package org.reactome.funcInt;

/**
 * This class is used to map enum to SMALLINT
 * @author guanming
 *
 */
public class HibGeneExpressionType extends IntEnumUserType<GeneExpressionType> {
    
    public HibGeneExpressionType() {
        // Pass all values to the super class.
        super(GeneExpressionType.class,
              GeneExpressionType.values());
    }
    
}
