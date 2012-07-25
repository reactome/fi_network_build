/*
 * Created on Jul 25, 2012
 *
 */
package org.reactome.data;

import java.util.Iterator;
import java.util.List;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;

/**
 * Customized ReactomeAnalyzer in order to control sizes of genes in pathways. In annotated NCI-PID
 * pathways, many pathways have been linked to other pathways. In this class, if a pathway is linked
 * to other pathways, but still contain reactions, other pathways will be excluded from generating
 * interactions. If a pathway contains just sub-pathways, its sub-pathways will be used for generating
 * gene sets, but not this pathway.
 * @author gwu
 *
 */
public class AnnotatedNCIPICAnalyzer extends CPathAnalyzer {
    
    public AnnotatedNCIPICAnalyzer() {
    }

    @Override
    protected List<GKInstance> getTopics() throws Exception {
        List<GKInstance> pathways = super.getTopics();
        for (Iterator<GKInstance> it = pathways.iterator(); it.hasNext();) {
            GKInstance pathway = (GKInstance) it.next();
            if (hasSubPathwaysOnly(pathway))
                it.remove(); // There are no need to process these pathways, since sub-pathways will be processed.
            excludeSubPathway(pathway);
        }
        return pathways;
    }
    
    /**
     * Linked sub-pathways will be excluded in the attribute list for further processes:
     * generating functional intearctions and gene sets.
     * @param pathway
     * @return
     * @throws Exception
     */
    private void excludeSubPathway(GKInstance pathway) throws Exception {
        List<GKInstance> list = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        for (Iterator<GKInstance> it = list.iterator(); it.hasNext();) {
            GKInstance inst = it.next();
            if (inst.getSchemClass().isa(ReactomeJavaConstants.Pathway))
                it.remove();
        }
    }
    
    /**
     * Check if this passed pathway is used as a container, which contains sub-pathways only 
     * without reactions.
     * @param pathway
     * @return
     * @throws Exception
     */
    private boolean hasSubPathwaysOnly(GKInstance pathway) throws Exception {
        List<GKInstance> list = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        if (list == null || list.size() == 0)
            return true;
        for (GKInstance inst : list) {
            if (inst.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent))
                return false;
        }
        return true;
    }
    
}
