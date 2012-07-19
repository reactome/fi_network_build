/*
 * Created on Jul 19, 2012
 *
 */
package org.reactome.data;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;

/**
 * @author gwu
 *
 */
public class EncodeInteractionAnalyzer extends TargetedInteractionAnalyzer {
    
    public EncodeInteractionAnalyzer() {
        
    }

    /**
     * Only ENCODE interactions have extra evidences should be exported as functional
     * interactions.
     */
    @Override
    protected boolean isNeededInteraction(GKInstance interaction) throws Exception {
        String definition = (String) interaction.getAttributeValue(ReactomeJavaConstants.definition);
        if (definition != null && definition.contains("supported by"))
            return true;
        return false;
    }
    
    
}
