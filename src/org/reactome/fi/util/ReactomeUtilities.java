/*
 * Created on Jun 9, 2017
 *
 */
package org.reactome.fi.util;

import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;

/**
 * Some extra utility methods collected here.
 * @author gwu
 *
 */
public class ReactomeUtilities {
    
    public static void grepGenesFromEntity(GKInstance pe,
                                           Set<String> genes) throws Exception {
        Set<GKInstance> refEntities = InstanceUtilities.grepReferenceEntitiesForPE(pe);
        for (GKInstance refEntity : refEntities) {
            if (refEntity.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
                String geneName = (String) refEntity.getAttributeValue(ReactomeJavaConstants.geneName);
                if (geneName != null)
                    genes.add(geneName);
            }
        }
    }
    
}
