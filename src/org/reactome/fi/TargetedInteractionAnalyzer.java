/*
 * Created on Apr 15, 2009
 *
 */
package org.reactome.fi;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;

/**
 * This class is used to fetch FIs from Targetted interactions.
 * @author wgm
 *
 */
public class TargetedInteractionAnalyzer extends CPathAnalyzer {
    
    public TargetedInteractionAnalyzer() {
    }

    @Override
    protected Collection prepareComplexes() throws Exception {
        return new HashSet<GKInstance>();
    }

    @Override
    protected Collection prepareReactions() throws Exception {
        return new HashSet<GKInstance>();
    }
    
    @Override
    protected Set<String> extractInteractionEvents() throws Exception {
        Set<String> interactions = new HashSet<String>();
        Collection eventInteractions = prepareInteractions();
        GKInstance eventInteraction = null;
        GKInstance factor;
        GKInstance target;
        Set<GKInstance> interactorSet = new HashSet<GKInstance>();
        for (Iterator it = eventInteractions.iterator(); it.hasNext();) {
            eventInteraction = (GKInstance) it.next();
            factor = (GKInstance) eventInteraction.getAttributeValue(ReactomeJavaConstants.factor);
            target = (GKInstance) eventInteraction.getAttributeValue(ReactomeJavaConstants.target);
            if (factor == null || target == null)
                continue;
            interactorSet.clear();
            interactorSet.add(factor);
            interactorSet.add(target);
            generateInteractions(interactorSet, interactions, eventInteraction);
        }
        return interactions;
    }
    
    @Override
    protected Collection prepareInteractions() throws Exception {
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) {
            // Load necessary attributes
            MySQLAdaptor dba = (MySQLAdaptor) adaptor;
            GKInstance dataSource = getDataSource();
            Collection interactions = null;
            if (dataSource == null)
                interactions = dba.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
            else
                interactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.TargettedInteraction,
                                                            ReactomeJavaConstants.dataSource,
                                                            "=",
                                                            dataSource);
            SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.TargettedInteraction);
            SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.factor);
            dba.loadInstanceAttributeValues(interactions, att);
            att = cls.getAttribute(ReactomeJavaConstants.target);
            dba.loadInstanceAttributeValues(interactions, att);
            return interactions;
        }
        else {
            XMLFileAdaptor fileAdaptor = (XMLFileAdaptor) adaptor;
            return fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
        }
    }

    @Override
    protected List<GKInstance> getTopics() throws Exception {
        return new ArrayList<GKInstance>();
    }
    
    
}
