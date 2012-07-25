/*
 * Created on Jul 10, 2007
 *
 */
package org.reactome.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This helper class is used to process topic related tasks.
 * @author guanming
 *
 */
@SuppressWarnings("unchecked")
public class ReactomeAnalyzerTopicHelper {
    
    public ReactomeAnalyzerTopicHelper() {
    }
    
    /**
     * Get the numbers of ids used in a passed pathway instance.
     * @param topic
     * @return
     * @throws Exception
     */
    protected Map<String, Integer> grepIDToNumberFromTopic(GKInstance topic) throws Exception {
        List<String> idList = new ArrayList<String>();
        List<GKInstance> participants = grepPathwayParticipantsInList(topic);
        for (GKInstance participant : participants) {
            Set<GKInstance> refPeptides = grepRefPepSeqs(participant);
            for (GKInstance ref : refPeptides) {
                GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                // It should be faster by comparing DB_ID
                if (refDb == null)
                    continue;
                String identifier = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier != null)
                    idList.add(identifier);
            }
        }
        Map<String, Integer> idToNumber = InteractionUtilities.countTermUsageInList(idList);
        return idToNumber;
    }
    
    /**
     * Get protein identifiers in a passed pathway (topic) instance.
     * @param topic
     * @return
     * @throws Exception
     */
    protected Set<String> grepIDsFromTopic(GKInstance topic) throws Exception {
        Set<String> ids = new HashSet<String>();
        // First load all PhysicalEntities involved in Reactions
        Set<GKInstance> participants = grepPathwayParticipants(topic);
        // Grep ReferencePeptideSequence
        for (GKInstance participant : participants) {
            Set<GKInstance> refPeptides = grepRefPepSeqs(participant);
            for (GKInstance ref : refPeptides) {
                GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                // It should be faster by comparing DB_ID
                if (refDb == null || !refDb.getDBID().equals(2L))
                    continue;
                String identifier = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier != null)
                    ids.add(identifier);
            }
        }
        return ids;
    }
    
    /**
     * Get a set of RefPepSeq instances from a passed PhysicalEntity.
     * @param interactor
     * @return
     * @throws Exception
     */
    public Set<GKInstance> grepRefPepSeqs(GKInstance interactor) throws Exception {
        Set<GKInstance> refPepSeq = new HashSet<GKInstance>();
        if (interactor.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
            GKInstance ref = (GKInstance) interactor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (ref != null) {
                if (ref.getSchemClass().isa(ReactomeJavaConstants.ReferencePeptideSequence) ||
                    ref.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct))
                    refPepSeq.add(ref);
            }
        }
        else if (interactor.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            grepRefPepSeqFromInstanceRecursively(interactor, refPepSeq);
        }
        else if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
            grepRefPepSeqFromInstanceRecursively(interactor, refPepSeq);
        }
        return refPepSeq;
    }
    
    private void grepRefPepSeqFromInstanceRecursively(GKInstance complex, 
                                                      Set<GKInstance> refPepSeq) throws Exception {
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(complex);
        Set<GKInstance> next = new HashSet<GKInstance>();
        Set<GKInstance> children = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (GKInstance inst : current) {
                children.clear();
                if (inst.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent)) {
                    List list = inst.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                    if (list != null)
                        children.addAll(list);
                }
                if (inst.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember)) {
                    List list = inst.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                    if (list != null)
                        children.addAll(list);
                }
                // Check for candidate set: added on June 24, 2009
                if (inst.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasCandidate)) {
                    List list = inst.getAttributeValuesList(ReactomeJavaConstants.hasCandidate);
                    if (list != null)
                        children.addAll(list);
                }
                if (children.size() == 0)
                    continue;
                for (Iterator it = children.iterator(); it.hasNext();) {
                    GKInstance tmp = (GKInstance) it.next();
                    if (tmp.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                        next.add(tmp);
                    else if (tmp.getSchemClass().isa(ReactomeJavaConstants.Complex))
                        next.add(tmp);
                    else if (tmp.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
                        GKInstance ref = (GKInstance) tmp.getAttributeValue(ReactomeJavaConstants.referenceEntity);
                        if (ref != null &&
                            (ref.getSchemClass().isa(ReactomeJavaConstants.ReferencePeptideSequence) ||
                             ref.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)))
                            refPepSeq.add(ref);
                    }
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    /**
     * Get event component in a specified pathway instance. Complex is not included.
     * @param pathway
     * @return
     * @throws Exception
     */
    public Set<GKInstance> grepPathwayEventComponents(GKInstance pathway) throws Exception {
        return InstanceUtilities.grepPathwayEventComponents(pathway);
    }
    
    /**
     * Grep reactions, interactions, and complexes into a set. Note: complexes have been
     * loaded into this set.
     * @param pathway
     * @return
     * @throws Exception
     */
    protected Set<GKInstance> grepPathwayComponents(GKInstance pathway) throws Exception {
        // First load all PhysicalEntities involved in Reactions
        Set<GKInstance> components = new HashSet<GKInstance>();
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(pathway);
        Set<GKInstance> next = new HashSet<GKInstance>();
        // To avoid self-containing as in Nature-PID
        Set<GKInstance> checked = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (GKInstance tmp : current) {
                checked.add(tmp);
                if (tmp.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                    components.add(tmp);
                    // Need to get complex
                    Set reactionParticipants = InstanceUtilities.getReactionParticipants(tmp);
                    for (Iterator it = reactionParticipants.iterator(); it.hasNext();) {
                        GKInstance p = (GKInstance) it.next();
                        if (p.getSchemClass().isa(ReactomeJavaConstants.Complex))
                            components.add(p);
                    }
                }
                else if (tmp.getSchemClass().isa(ReactomeJavaConstants.Interaction)) {
                    components.add(tmp);
                    List interactors = tmp.getAttributeValuesList(ReactomeJavaConstants.interactor);
                    if (interactors != null) {
                        for (Iterator it = interactors.iterator(); it.hasNext();) {
                            GKInstance i = (GKInstance) it.next();
                            if (i.getSchemClass().isa(ReactomeJavaConstants.Complex))
                                components.add(i);
                        }
                    }
                }
                else if (tmp.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent)) {
                    List values = tmp.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                    if (values != null)
                        next.addAll(values);
                }
                else if (tmp.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    next.addAll(tmp.getAttributeValuesList(ReactomeJavaConstants.hasEvent));
                else if (tmp.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasSpecialisedForm))
                    next.addAll(tmp.getAttributeValuesList(ReactomeJavaConstants.hasSpecialisedForm));
                else if (tmp.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember))
                    next.addAll(tmp.getAttributeValuesList(ReactomeJavaConstants.hasMember));
            }
            current.clear();
            current.addAll(next);
            next.clear();
            current.removeAll(checked);
        }
        return components;
    }
    
    /**
     * Grep Entity pathway participants in a list so that the usage of each entity can be counted.
     * @param pathway
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    private List<GKInstance> grepPathwayParticipantsInList(GKInstance pathway) throws Exception {
        // First load all PhysicalEntities involved in Reactions
        List<GKInstance> participants = new ArrayList<GKInstance>();
        // Complexes have be pushed into this set too.
        Set<GKInstance> components = grepPathwayEventComponents(pathway);
        for (GKInstance tmp : components) {
            if (tmp.getSchemClass().isa(ReactomeJavaConstants.Reaction)) {
                List inputs = tmp.getAttributeValuesList(ReactomeJavaConstants.input);
                if (inputs != null)
                    participants.addAll(inputs);
                List outputs = tmp.getAttributeValuesList(ReactomeJavaConstants.output);
                if (outputs != null)
                    participants.addAll(outputs);
                List cas = tmp.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
                if (cas != null && cas.size() > 0) {
                    for (Iterator it = cas.iterator(); it.hasNext();) {
                        GKInstance ca = (GKInstance) it.next();
                        List catalysts = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                        if (catalysts != null) {
                            for (Iterator it1 = catalysts.iterator(); it1.hasNext();) {
                                GKInstance catalyst = (GKInstance) it1.next();
                                if (catalyst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                                    participants.add(catalyst);
                            }
                        }
                    }
                }
                Collection regulations = tmp.getReferers(ReactomeJavaConstants.regulatedEntity);
                if (regulations != null && regulations.size() > 0) {
                    for (Iterator it = regulations.iterator(); it.hasNext();) {
                        GKInstance regulation = (GKInstance) it.next();
                        List regulators = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
                        if (regulators != null) {
                            for (Iterator it1 = regulators.iterator(); it1.hasNext();) {
                                GKInstance regulator = (GKInstance) it1.next();
                                if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                                    participants.add(regulator);
                            }
                        }
                    }
                }
            }
            else if (tmp.getSchemClass().isa(ReactomeJavaConstants.Interaction)) {
                List interactors = tmp.getAttributeValuesList(ReactomeJavaConstants.interactor);
                if (interactors != null)
                    participants.addAll(interactors);
            }
        }
        return participants;
    }
    
    /**
     * Get entity components in a specified pathway instance.
     * @param pathway
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    public Set<GKInstance> grepPathwayParticipants(GKInstance pathway) throws Exception {
        return InstanceUtilities.grepPathwayParticipants(pathway);
    }
    
    /**
     * This method is used to generate a list of numbers of pathway event particiapnts and entity participants
     * for a list of pathways, which are fetched from a local file. The file contains a list of pathway ids
     * only: one id per line.
     * @throws Exception
     */
    @Test
    public void tallyPathways() throws Exception {
        String dirName = "/Users/wgm/Documents/gkteam/peter/";
        String fileName = dirName + "pathway_IDs.txt";
        FileUtility fu = new FileUtility();
        Set<String> idSet = fu.loadInteractions(fileName);
        List<GKInstance> instances = new ArrayList<GKInstance>();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "gk_central_031309",
                                            "root",
                                            "macmysql01",
                                            3306);
        for (String id : idSet) {
            Long dbId = Long.parseLong(id);
            GKInstance pathway = dba.fetchInstance(dbId);
            if (pathway == null)
                throw new IllegalStateException(id + " cannot be found!");
            instances.add(pathway);
        }
        InstanceUtilities.sortInstances(instances);
        fileName = dirName + "PathwaysWithNumbers.txt";
        fu.setOutput(fileName);
        fu.printLine("Pathway\tDB_ID\tNumberOfReactions\tNumberOfEntities\tRatioOfReactionToEntity");
        for (GKInstance pathway : instances) {
            Set<GKInstance> eventComponents = grepPathwayEventComponents(pathway);
            Set<GKInstance> entityComponents = grepPathwayParticipants(pathway);
            double ratio = (double) eventComponents.size() / entityComponents.size();
            fu.printLine(pathway.getDisplayName() + "\t" +
                         pathway.getDBID() + "\t" + 
                         eventComponents.size() + "\t" +
                         entityComponents.size() + "\t" + 
                         ratio);
        }
        fu.close();
    }
    
    @Test
    public void checkGenesInPathways() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("GENE_TO_TOPIC");
        FileUtility fu = new FileUtility();
        Map<String, Set<String>> geneToPathways = fu.loadSetMap(fileName);
        Map<String, Set<String>> pathwayToGenes = InteractionUtilities.switchKeyValues(geneToPathways);
//        System.out.println("Total pathways: " + pathwayToGenes.size());
        for (String pathway : pathwayToGenes.keySet()) {
            if (!pathway.endsWith("(N)"))
                continue;
            Set<String> genes = pathwayToGenes.get(pathway);
            System.out.println(pathway + "\t" + genes.size());
        }
    }
    
}
