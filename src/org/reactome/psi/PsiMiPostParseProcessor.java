/*
 * Created on Aug 22, 2006
 *
 */
package org.reactome.psi;

/**
 * This interface is used to handle the job after a PSI-MI file is parsed. For example,
 * save into a database using hibernate or convert it the Reactome data model.
 * @author guanming
 *
 */
public interface PsiMiPostParseProcessor {

    public void postProcess(Object obj, 
                            PsiMiModel model) throws Exception;
    
}
