/*
 * Created on Sep 5, 2006
 *
 */
package org.reactome.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

public class CPathAnalyzer extends PantherAnalyzer {
    
    public CPathAnalyzer() {
        //dataSourceId = 229095L; // Cancer Cell Map
        //dataSourceId = 424025L; // NCI-Nature Curated Pathways
        //dataSourceId = 435967L; // NCI-Nature imported BiCarta Pathways
        //dataSourceId = 444163L; // KEGG
        //dataSourceId = 379958L; // BIND
        //dataSourceId = 317977L; // HPRD
        dataSourceId = 415130L; // IntAct
    }
    
    public Set<String> extractInteractionSet() throws Exception {
        Set<String> interactions = super.extractInteractionSet();
        Set<String> dataFromInteractions = extractInteractionEvents();
        System.out.println("Total interactions from Interaction Event: " + dataFromInteractions.size());
        interactions.addAll(dataFromInteractions);
        return interactions;
    }
    
    @Override
    public Set<String> extractInteractionSetWithComplexAndSet() throws Exception {
        Set<String> interactions = super.extractInteractionSetWithComplexAndSet();
        Set<String> dataFromInteractions = extractInteractionEventsWithComplexAndSet();
        System.out.println("Total interactions from Interaction Event: " + dataFromInteractions.size());
        interactions.addAll(dataFromInteractions);
        return interactions;
    }
    
    protected Set<String> extractInteractionEventsWithComplexAndSet() throws Exception {
        Set<String> interactions = new HashSet<String>();
        Collection eventInteractions = prepareInteractions();
        GKInstance eventInteraction = null;
        List interactors = null;
        Set<GKInstance> interactorSet = new HashSet<GKInstance>();
        for (Iterator it = eventInteractions.iterator(); it.hasNext();) {
            eventInteraction = (GKInstance) it.next();
            interactors = eventInteraction.getAttributeValuesList(ReactomeJavaConstants.interactor);
            if (interactors == null || interactors.size() == 0)
                continue;
            interactorSet.clear();
            for (Iterator it1 = interactors.iterator(); it1.hasNext();) {
                interactorSet.add((GKInstance)it1.next());
            }
            generateInteractionsWithComplexAndSet(interactorSet, 
                                                  interactions, 
                                                  eventInteraction);
        }
        return interactions;
    }

    protected Set<String> extractInteractionEvents() throws Exception {
        Set<String> interactions = new HashSet<String>();
        Collection eventInteractions = prepareInteractions();
        GKInstance eventInteraction = null;
        List interactors = null;
        Set<GKInstance> interactorSet = new HashSet<GKInstance>();
        for (Iterator it = eventInteractions.iterator(); it.hasNext();) {
            eventInteraction = (GKInstance) it.next();
            // Avoid loading missing interactions
            String interactionType = (String) eventInteraction.getAttributeValue(ReactomeJavaConstants.interactionType);
            if (interactionType != null && interactionType.contains("missing interaction"))
                continue;
            interactors = eventInteraction.getAttributeValuesList(ReactomeJavaConstants.interactor);
            if (interactors == null || interactors.size() == 0)
                continue;
            interactorSet.clear();
            for (Iterator it1 = interactors.iterator(); it1.hasNext();) {
                interactorSet.add((GKInstance)it1.next());
            }
            // Try to get PPIs less than 4: expermentally.
            if (interactorSet.size() > 2)
                continue; 
            generateInteractions(interactorSet, interactions, eventInteraction);
        }
        return interactions;
    }
    
    protected Collection prepareInteractions() throws Exception {
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) {
            // Load necessary attributes
            MySQLAdaptor dba = (MySQLAdaptor) adaptor;
            GKInstance cpath = getDataSource();
            Collection reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Interaction,
                                                                ReactomeJavaConstants.dataSource,
                                                                "=",
                                                                cpath);
            // Load precedingEvent values
            SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Interaction);
            SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.interactor);
            dba.loadInstanceAttributeValues(reactions, att);
            att = cls.getAttribute(ReactomeJavaConstants.interactionType);
            dba.loadInstanceAttributeValues(reactions, att);
            return reactions;
        }
        else {
            XMLFileAdaptor fileAdaptor = (XMLFileAdaptor) adaptor;
            return fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        }
    }
    
    /**
     * hasComponent should be replaced by hasEvent. However, it is not done in the database.
     * @throws Exception
     */
    @Test
    public void fixHasEventProblems() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "CellMap_Fix_hasEvent.rtpj";
        fileAdaptor.setSource(fileName);
        Collection events = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Event);
        System.out.println("Total events: " + events.size());
        Set<String> names = new HashSet<String>();
        Map<String, GKInstance> crossRefToEvent = new HashMap<String, GKInstance>();
        // Find a single mapping
        for (Object event : events) {
            GKInstance inst = (GKInstance) event;
            GKInstance crossRef = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.crossReference);
            names.add(crossRef.getDisplayName());
            crossRefToEvent.put(crossRef.getDisplayName(), inst);
        }
        System.out.println("Total names: " + names.size());
        if (names.size() != events.size())
            throw new IllegalStateException("Keys are not unique!");
        // Get the original project.
        XMLFileAdaptor adaptor1 = new XMLFileAdaptor();
        adaptor1.setSource("datasets/cellmap_may_2006/CellMap.rtpj");
        Collection pathways = adaptor1.fetchInstancesByClass(ReactomeJavaConstants.Pathway);
        Map<String, GKInstance> nameToPathway = new HashMap<String, GKInstance>();
        for (Object obj : pathways) {
            GKInstance pathway = (GKInstance) obj;
            nameToPathway.put(pathway.getDisplayName(), pathway);
        }
        Collection newPathways = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Pathway);
        for (Object obj : newPathways) {
            GKInstance newPathway = (GKInstance) obj;
            GKInstance oldPathway = nameToPathway.get(newPathway.getDisplayName());
            List oldValues = oldPathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
            List newValues = new ArrayList<GKInstance>(oldValues.size());
            for (Object value : oldValues) {
                GKInstance inst = (GKInstance) value;
                GKInstance crossRef = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.crossReference);
                GKInstance newInst = crossRefToEvent.get(crossRef.getDisplayName());
                if (newInst == null)
                    throw new IllegalStateException("Cannot find event for: " + crossRef.getDisplayName());
                newValues.add(newInst);
            }
            newPathway.setIsDirty(true);
            newPathway.setAttributeValue(ReactomeJavaConstants.hasEvent, newValues);
        }
        fileAdaptor.save();
    }
    
}
