/*
 * Created on Jul 10, 2007
 *
 */
package org.reactome.data;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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
@SuppressWarnings({"unchecked", "rawtypes"})
public class ReactomeAnalyzerTopicHelper {
	// A flag to control is hadCandidate is used
	private boolean needCandidateMemebers = false;
    
    public ReactomeAnalyzerTopicHelper() {
    }
    
    public boolean isNeedCandidateMemebers() {
		return needCandidateMemebers;
	}

	public void setNeedCandidateMemebers(boolean needCandidateMemebers) {
		this.needCandidateMemebers = needCandidateMemebers;
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
        Set<GKInstance> refSeqs = new HashSet<GKInstance>();
        grepRefSeqs(interactor, refSeqs);
        return refSeqs;
    }
    
    /**
     * This is a big change: as of September 23, 2019, both ReferenceGeneProdcuts (for proteins)
     * and REferenceDNASequence (for genes) are collected so that we can catch FIs between two genes, 
     * and proteins and genes.
     * @param pe
     * @param refSeqs
     * @throws Exception
     */
    private void grepRefSeqs(GKInstance pe, Set<GKInstance> refSeqs) throws Exception {
    	Set<GKInstance> ewases = null;
    	if (needCandidateMemebers) {
    		ewases = InstanceUtilities.getContainedInstances(pe,
                    ReactomeJavaConstants.hasComponent,
                    ReactomeJavaConstants.hasMember,
                    ReactomeJavaConstants.hasCandidate);
    	}
    	else {
    		// As of December 15, 2014, hasCandidate will not be used, which
    		// reduces the total FIs about 12% (from 144733 to 127382).
    		ewases = InstanceUtilities.getContainedInstances(pe,
    				ReactomeJavaConstants.hasComponent,
    				ReactomeJavaConstants.hasMember);
    	}
        ewases.add(pe);
        for (GKInstance ewas : ewases) {
            if (!ewas.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                continue;
            GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refEntity == null)
                continue;
            if (refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct) ||
                refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceDNASequence)) {
                refSeqs.add(refEntity);
            }
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
                Collection regulations = InstanceUtilities.getRegulations(tmp);
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
    
    @Test
    public void countTotalPathways() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("GENE_TO_TOPIC");
        Set<String> pathways = Files.lines(Paths.get(fileName)).map(line -> line.split("\t")[1]).collect(Collectors.toSet());
        System.out.println("Total pathways: " + pathways.size());
    }
    
    private Set<String> generateFIs(Set<GKInstance> interactors,
                                    ReactomeAnalyzer analyzer) throws Exception {
        Set<String> fis = new HashSet<String>();
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        for (int i = 0; i < list.size() - 1; i++) {
            GKInstance int1 = list.get(i);
            Set<GKInstance> refs1 = InstanceUtilities.grepReferenceEntitiesForPE(int1);
            if (refs1.size() == 0)
                continue;
            for (Iterator<GKInstance> it = refs1.iterator(); it.hasNext();) {
                GKInstance ref = it.next();
                GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                if (!refDb.getDisplayName().equals("UniProt"))
                    it.remove();;
            }
            if (refs1.size() == 0)
                continue;
            for (int j = i + 1; j < list.size(); j++) {
                GKInstance int2 = list.get(j);
                Set<GKInstance> refs2 = InstanceUtilities.grepReferenceEntitiesForPE(int2);
                if (refs2.size() == 0)
                    continue;
                for (Iterator<GKInstance> it = refs2.iterator(); it.hasNext();) {
                    GKInstance ref = it.next();
                    GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                    if (!refDb.getDisplayName().equals("UniProt"))
                        it.remove();;
                }
                if (refs2.size() == 0)
                    continue;
                analyzer.generateFIs(refs1, refs2, fis, false);
            }
        }
//        if (fis.size() > 0)
//            return fis;
//        // Try to use complexes
//        for (GKInstance inst : interactors) {
//            if (inst.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
//                Set<GKInstance> complexInteractors = 
//            }
//        }
        return fis;
    }
    
    @Test
    public void checkReactomeFIsCoverage() throws Exception {
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Collection<GKInstance> reactions = analyzer.prepareReactions();
        // Use reactions only
        for (Iterator<GKInstance> it = reactions.iterator(); it.hasNext();) {
            GKInstance inst = it.next();
            if (inst.getSchemClass().isa(ReactomeJavaConstants.Reaction))
                continue;
            it.remove();
        }
        Collection<GKInstance> complexes = analyzer.prepareComplexes();
        System.out.println("Total reactions: " + reactions.size());
        Map<GKInstance, Set<String>> rxtToFIs = new HashMap<GKInstance, Set<String>>();
        int noIntCount = 0;
        int proteinChemicalCount = 0;
        int oneProteinOneChemical = 0;
        int oneFICount = 0;
        int associtations = 0;
        for (GKInstance rxt : reactions) {
            Set<GKInstance> interactors = new HashSet<GKInstance>();
            analyzer.extractInteractorsFromReaction(rxt, interactors);
            Set<String> fis = generateFIs(interactors, analyzer);
            if (fis.size() > 0) {
                rxtToFIs.put(rxt, fis);
                if (fis.size() == 1)
                    oneFICount ++;
                continue;
            }
            Set<GKInstance> proteins = new HashSet<GKInstance>();
            Set<GKInstance> chemicals = new HashSet<GKInstance>();
            for (GKInstance inst : interactors) {
                Set<GKInstance> refs = InstanceUtilities.grepReferenceEntitiesForPE(inst);
                for (GKInstance ref : refs) {
                    GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                    if (refDb.getDisplayName().equals("ChEBI"))
                        chemicals.add(ref);
                    else if (refDb.getDisplayName().equals("UniProt"))
                        proteins.add(ref);
                }
            }
            if (proteins.size() > 0 && chemicals.size() > 0) {
                proteinChemicalCount ++;
                if (proteins.size() == 1 && chemicals.size() == 1)
                    oneProteinOneChemical ++;
            }
            else {
                // Check if it is associtation
                List<GKInstance> outputs = rxt.getAttributeValuesList(ReactomeJavaConstants.output);
                if (outputs.size() == 1) {
                    GKInstance output = outputs.get(0);
                    if (output.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
                        List<GKInstance> inputs = rxt.getAttributeValuesList(ReactomeJavaConstants.input);
                        if (inputs.size() > 1) {
                            associtations ++;
                            continue;
                        }
                    }
                }
                noIntCount ++;
//                System.out.println(rxt);
            }
        }
        System.out.println("Reactions having FIs: " + rxtToFIs.size());
        System.out.println("Reactions having one FI: " + oneFICount);
        System.out.println("Reactions don't have FIs: " + noIntCount);
        System.out.println("Interactions between proteins and chemicals: " + proteinChemicalCount);
        System.out.println("One protein and one chemical: " + oneProteinOneChemical);
        System.out.println("Possible associtations: " + associtations); // This number is added to reactions having FIs to get the total reactions that have FIs.
//        if (true)
//            return;
        
        // Load Proteins in the Reactome FI
        String resultDir = FIConfiguration.getConfiguration().get("RESULT_DIR");
        // Interactions
        String[] files = new String[] {
                FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"),
                FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_TO_HUMAN_PPI_FILE"),
                FIConfiguration.getConfiguration().get("IREFINDEX_FLY_TO_HUMAN_PPI_FILE"),
                FIConfiguration.getConfiguration().get("IREFINDEX_WORM_TO_HUMAN_PPI_FILE"),
                FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_TO_HUMAN_PPI_FILE")
        };
        FileUtility fu = new FileUtility();
        Set<String> interactions = new HashSet<String>();
        for (String file : files) {
            Set<String> tmp = fu.loadInteractions(file);
//            System.out.println("Interactions in " + file + ": " + tmp.size());
            interactions.addAll(tmp);
//            System.out.println("Merged: " + interactions.size());
        }
        System.out.println("Total interactions: " + interactions.size());
        Map<GKInstance, Set<String>> cannotMapped = new HashMap<GKInstance, Set<String>>();
        for (GKInstance rxt : rxtToFIs.keySet()) {
            Set<String> fis = rxtToFIs.get(rxt);
            Set<String> shared = InteractionUtilities.getShared(fis, interactions);
            if (shared.size() == 0 && fis.size() > 1) {
                cannotMapped.put(rxt, fis);
            }
        }
        System.out.println("Reactions cannot map to interactions: " + cannotMapped.size());
        // Check with domain domain interactions
        PfamAnalyzer pfamAnalyzer = new PfamAnalyzer();
        int count = 0;
        int cannotMappedCount = 0;
        for (GKInstance rxt : rxtToFIs.keySet()) {
            Set<String> fis = rxtToFIs.get(rxt);
            for (String fi : fis) {
                boolean isMapped = pfamAnalyzer.checkIfInteracting(fi);
                if (isMapped) {
                    count ++;
                    if (cannotMapped.keySet().contains(rxt))
                        cannotMappedCount ++;
                    break;
                }
            }
        }
        System.out.println("Reactions mapped to domain interactions: " + count);
        System.out.println("No interaction reactions mapped to domain interactions: " + cannotMappedCount);
    }
    
}
