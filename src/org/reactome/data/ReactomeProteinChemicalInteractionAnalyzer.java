/*
 * Created on Jan 30, 2017
 *
 */
package org.reactome.data;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to analyze interactions between proteins and small molecules.
 * @author gwu
 *
 */
@SuppressWarnings("unchecked")
public class ReactomeProteinChemicalInteractionAnalyzer extends ReactomeAnalyzer {
    
    /**
     * Default constructor.
     */
    public ReactomeProteinChemicalInteractionAnalyzer() {
    }
    
    /**
     * Call this method to generate interactions between proteins and chemicals.
     */
    @Test
    public void extractInteractions() throws Exception {
        Collection<GKInstance> reactions = prepareReactions();
        Collection<GKInstance> complexes = prepareComplexes();
//        complexes.clear();
        Map<String, Set<GKInstance>> interactionToSources = extractInteractions(reactions, complexes);
        
        outputInteractions(interactionToSources, complexes.size() == 0 ? false : true);
        
        // The following code is used for test
//        Long pathwayId = 1257604L; // PIP3 activates AKT signaling
//        GKInstance pathway = getMySQLAdaptor().fetchInstance(pathwayId);
//        Set<GKInstance> events = InstanceUtilities.grepPathwayEventComponents(pathway);
//        Set<GKInstance> entities = InstanceUtilities.grepPathwayParticipants(pathway);
//        checkReactions(reactions);
        // Check for complexes
//        checkComplexes(complexes);
    }

    private void outputInteractions(Map<String, Set<GKInstance>> interactionToSources,
                                    boolean needComplex) throws Exception {
//        String dirName = FIConfiguration.getConfiguration().get("RESULT_DIR");
//        String fileName = dirName + File.separator + "FIs_ProteinChemical.txt";
        
        String fileName = "/Users/gwu/git/Ogmios/results/FIs_ProteinChemical_032017.txt";
        fileName = "/Users/gwu/git/Ogmios/results/FIs_ProteinChemicalWithComplex_032017.txt";
        
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        if (needComplex)
            fu.printLine("ChemicalName\tCheEBIId\tGeneName\tUniProtId\tReactionIds\tComplexeIds");
        else 
            fu.printLine("ChemicalName\tCheEBIId\tGeneName\tUniProtId\tReactionIds");
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        StringBuilder builder = new StringBuilder();
        for (String interaction : interactionToSources.keySet()) {
            String[] ids = interaction.split(":");
            GKInstance chemical = dba.fetchInstance(new Long(ids[0]));
            builder.append(chemical.getAttributeValue(ReactomeJavaConstants.name));
            builder.append("\t").append(chemical.getAttributeValue(ReactomeJavaConstants.identifier));
            GKInstance protein = dba.fetchInstance(new Long(ids[1]));
            builder.append("\t").append(protein.getAttributeValue(ReactomeJavaConstants.name));
            builder.append("\t").append(protein.getAttributeValue(ReactomeJavaConstants.identifier));
            Set<GKInstance> sources = interactionToSources.get(interaction);
            // Get reactions first
            builder.append("\t");
            for (GKInstance source : sources) {
                if (source.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) 
                    builder.append(source.getDBID()).append(",");
            }
            if (builder.charAt(builder.length() - 1) == ',')
                builder.deleteCharAt(builder.length() - 1);
            if (needComplex) {
                // Get complexes
                builder.append("\t");
                for (GKInstance source : sources) {
                    if (source.getSchemClass().isa(ReactomeJavaConstants.Complex)) 
                        builder.append(source.getDBID()).append(",");
                }
                if (builder.charAt(builder.length() - 1) == ',')
                    builder.deleteCharAt(builder.length() - 1);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    private Map<String, Set<GKInstance>> extractInteractions(Collection<GKInstance> reactions,
                                                             Collection<GKInstance> complexes) throws Exception {
        // Key: id of protein : id of chemical; Value: Set of sources
        Map<String, Set<GKInstance>> interactionToSources = new HashMap<String, Set<GKInstance>>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        int totalCheckedReactions = 0;
        for (GKInstance rxn : reactions) {
            if (!isChemicalReaction(rxn))
                continue;
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactionToSources, rxn);
            interactors.clear();
            totalCheckedReactions ++;
        }
        System.out.println("Total checked reactions: " + totalCheckedReactions);
        System.out.println("Total interactions from reactions: " + interactionToSources.size());
        int totalCheckedComplexes = 0;
        for (GKInstance complex : complexes) {
            if (!isChemicalComplex(complex))
                continue;
            interactors.clear();
            grepComplexComponents(complex, interactors);
            // One more step to get the small molecule
            List<GKInstance> components = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
            for (GKInstance comp : components) {
                interactors.add(comp); // This will add chemicals back
            }
            generateInteractions(interactors, interactionToSources, complex);
            totalCheckedComplexes ++;
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total checked complexes: " + totalCheckedComplexes);
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total interactions from Reactome: " + interactionToSources.size());
        return interactionToSources;
    }
    
    protected void generateInteractions(Set<GKInstance> interactors, 
                                        Map<String, Set<GKInstance>> interactionToSources,
                                        GKInstance source) throws Exception {
        // Split into two types
        Set<GKInstance> chemicals = new HashSet<GKInstance>();
        Set<GKInstance> proteins = new HashSet<GKInstance>();
        for (GKInstance interactor : interactors) {
            if (isMacromolecule(interactor))
                proteins.add(interactor);
            else if (isChemical(interactor))
                chemicals.add(interactor);
        }
        // Just in case
        if (chemicals.size() == 0 || proteins.size() == 0)
            return;
        for (GKInstance chemical : chemicals) {
            Set<GKInstance> refMolecules = new HashSet<GKInstance>();
            grepReferenceMolecules(chemical, refMolecules);
            if (refMolecules.size() == 0)
                continue;
            for (GKInstance protein : proteins) {
                Set<GKInstance> refGeneProducts = grepRefPepSeqs(protein);
                if (refGeneProducts.size() == 0)
                    continue;
                generateInteractions(refMolecules, 
                                     refGeneProducts, 
                                     interactionToSources,
                                     source);
            }
        }
    }
    
    private void generateInteractions(Set<GKInstance> refMolecules,
                                      Set<GKInstance> refProteins,
                                      Map<String, Set<GKInstance>> interactionToSource,
                                      GKInstance source) throws Exception {
        for (GKInstance refMol : refMolecules) {
            for (GKInstance refProt : refProteins) {
                String interaction = refMol.getDBID() + ":" + refProt.getDBID();
                InteractionUtilities.addElementToSet(interactionToSource, interaction, source);
            }
        }
    }
    
    private void grepReferenceMolecules(GKInstance entity,
                                        Set<GKInstance> refEntities) throws Exception {
        if (entity.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity)) {
            GKInstance refEntity = (GKInstance) entity.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refEntity != null && refEntity.getSchemClass().isa(ReactomeJavaConstants.ReferenceMolecule))
                refEntities.add(refEntity);
        }
        if (entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            Set<GKInstance> members = InstanceUtilities.getContainedInstances(entity,
                                                                              ReactomeJavaConstants.hasMember,
                                                                              ReactomeJavaConstants.hasCandidate);
            for (GKInstance member : members) {
                grepReferenceMolecules(member, refEntities);
            }
        }
    }
    
    private void checkComplexes(Collection<GKInstance> complexes) throws Exception {
        // Extract interactions from reactions first
        int count = 1;
        System.out.println("Count\tDB_ID\tName");
        for (GKInstance complex : complexes) {
            // Don't want to have disease reaction
            GKInstance disease = (GKInstance) complex.getAttributeValue(ReactomeJavaConstants.disease);
            if (disease != null)
                continue; // Escape it
            if (isChemicalComplex(complex)) {
                System.out.println(count + "\t" + 
                                   complex.getDBID() + "\t" +
                                   complex.getDisplayName());
                count ++;
            }
        }
    }
    
    private boolean isChemicalComplex(GKInstance complex) throws Exception {
        // Don't drill down for this type of analysis since we want to focus on
        // high reliable interaction between protein and simple entity.
        List<GKInstance> entities = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        // Search for complexes having both macromolecules and simple entity
        boolean hasChemical = false;
        boolean hasMacromolecule = false;
        for (GKInstance entity : entities) {
            if (isChemical(entity)) {
                hasChemical = true;
            }
            else if (isMacromolecule(entity))
                hasMacromolecule = true;
        }
        return (hasChemical && hasMacromolecule);
    }

    private void checkReactions(Collection<GKInstance> reactions) throws InvalidAttributeException, Exception {
        // Extract interactions from reactions first
        int count = 1;
        System.out.println("Count\tDB_ID\tName\tIsPhosphorylation\tIsDephosphorylation\tIsATPPowered");
        for (GKInstance event : reactions) {
            // Don't want to have disease reaction
            GKInstance disease = (GKInstance) event.getAttributeValue(ReactomeJavaConstants.disease);
            if (disease != null)
                continue; // Escape it
            if (!event.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) 
                continue;
            if (isChemicalReaction(event)) {
                boolean isPhosphorylation = isPhosphorylation(event);
                boolean isDephosphorylation = isDephosphorylation(event);
                boolean isATPPowered = isATPPowered(event);
                if (isATPPowered) { // Then the reaction should not be a (de)phosphorylation
                    isPhosphorylation = false;
                    isDephosphorylation = false; 
                }
                System.out.println(count + "\t" + event.getDBID() + "\t" +
                                   event.getDisplayName() + "\t" + 
                                   (isPhosphorylation ? 1 : 0) + "\t" +
                                   (isDephosphorylation ? 1 : 0) + "\t" + 
                                   (isATPPowered ? 1 : 0));
                count ++;
            }
        }
    }
    
    private boolean isMacromolecueInInput(GKInstance reaction) throws Exception {
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        if (inputs == null || inputs.size() == 0)
            return false;
        for (GKInstance input : inputs) {
            if (isMacromolecule(input))
                return true;
        }
        return false;
    }
    
    private boolean isChemical(GKInstance entity) throws Exception {
        if (entity.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity))
            return true;
        if (entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            Set<GKInstance> members = InstanceUtilities.getContainedInstances(entity,
                                                                              ReactomeJavaConstants.hasMember,
                                                                              ReactomeJavaConstants.hasCandidate);
            for (GKInstance member : members) {
                boolean rtn = isChemical(member);
                if (rtn)
                    return true;
            }
        }
        return false;
    }
    
    private boolean isMacromolecule(GKInstance entity) throws Exception {
        if (entity.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
            entity.getSchemClass().isa(ReactomeJavaConstants.GenomeEncodedEntity) ||
            entity.getSchemClass().isa(ReactomeJavaConstants.Polymer))
            return true;
        if (entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
            Set<GKInstance> members = InstanceUtilities.getContainedInstances(entity,
                                                                              ReactomeJavaConstants.hasMember,
                                                                              ReactomeJavaConstants.hasCandidate);
            for (GKInstance member : members) {
                boolean rtn = isMacromolecule(member);
                if (rtn)
                    return true;
            }
        }
        return false;
    }
    
    private boolean isATPPowered(GKInstance reaction) throws Exception {
        Set<GKInstance> inputATP = new HashSet<GKInstance>();
        Set<GKInstance> outputADP = new HashSet<GKInstance>();
        Set<GKInstance> inputH2O = new HashSet<GKInstance>();
        Set<GKInstance> outputPi = new HashSet<GKInstance>();
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        for (GKInstance input : inputs) {
            if (isATP(input))
                inputATP.add(input);
            else if (isH2O(input))
                inputH2O.add(input);
        }
        List<GKInstance> outputs = reaction.getAttributeValuesList(ReactomeJavaConstants.output);
        for (GKInstance output : outputs) {
            if (isADP(output))
                outputADP.add(output);
            else if (isPi(output))
                outputPi.add(output);
        }
        if (inputATP.size() > 0 && inputH2O.size() > 0 && outputADP.size() > 0 && outputPi.size() > 0 &&
            inputATP.size() == inputH2O.size() && 
            outputADP.size() == outputPi.size() &&
            inputATP.size() == outputADP.size())
            return true;
        return false;
    }
    
    /**
     * A phosphorylation has one or more ATP and one and more ADP.
     * @param reaction
     * @return
     * @throws Exception
     */
    private boolean isPhosphorylation(GKInstance reaction) throws Exception {
        if (!isMacromolecueInInput(reaction))
            return false;
        Set<GKInstance> inputATP = new HashSet<GKInstance>();
        Set<GKInstance> outputADP = new HashSet<GKInstance>();
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        for (GKInstance input : inputs) {
            if (isATP(input))
                inputATP.add(input);
        }
        List<GKInstance> outputs = reaction.getAttributeValuesList(ReactomeJavaConstants.output);
        for (GKInstance output : outputs) {
            if (isADP(output))
                outputADP.add(output);
        }
        if (inputATP.size() > 0 && outputADP.size() > 0 && inputATP.size() == outputADP.size())
            return true;
        return false;
    }
    
    private boolean isDephosphorylation(GKInstance reaction) throws Exception {
        if (!isMacromolecueInInput(reaction))
            return false;
        Set<GKInstance> inputH2O = new HashSet<GKInstance>();
        Set<GKInstance> outputPi = new HashSet<GKInstance>();
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        for (GKInstance input : inputs) {
            if (isH2O(input))
                inputH2O.add(input);
        }
        List<GKInstance> outputs = reaction.getAttributeValuesList(ReactomeJavaConstants.output);
        for (GKInstance output : outputs) {
            if (isPi(output))
                outputPi.add(output);
        }
        if (inputH2O.size() > 0 && outputPi.size() > 0 && inputH2O.size() == outputPi.size())
            return true;
        return false;
    }
    
    private boolean isH2O(GKInstance entity) throws Exception {
        // 114728: water [ChEBI: 15377]
        return isReferredChemical(entity, 114728L);
    }
    
    private boolean isPi(GKInstance entity) throws Exception {
        // 114736: phosphate(3-) [ChEBI:18367]
        return isReferredChemical(entity, 114736L);
    }
    
    private boolean isATP(GKInstance entity) throws Exception {
        // 114729: ATP [ChEBI:15422]
        return isReferredChemical(entity, 114729L);
    }
    
    private boolean isADP(GKInstance entity) throws Exception {
        // 114735: ADP [ChEBI:16761]
        return isReferredChemical(entity, 114735L);
    }
    
    private boolean isReferredChemical(GKInstance entity,
                                       Long dbId) throws Exception {
        if (!entity.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity))
            return false;
        GKInstance refEntity = (GKInstance) entity.getAttributeValue(ReactomeJavaConstants.referenceEntity);
        if (refEntity == null)
            return false;
        return refEntity.getDBID().equals(dbId);
    }
    
    /**
     * Check if a reaction is a chemical related reaction, which is defined as a reaction using
     * a SimpleEntity as its input, catalyst, activator, or inhibitor.
     * @param reaction
     * @return
     * @throws Exception
     */
    private boolean isChemicalReaction(GKInstance reaction) throws Exception {
        Set<GKInstance> set = new HashSet<GKInstance>();
        List<GKInstance> inputs = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        for (GKInstance input : inputs) {
            if (isChemical(input))
                return true;
        }
        List<GKInstance> cas = reaction.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null && cas.size() > 0) {
            for (Iterator it = cas.iterator(); it.hasNext();) {
                GKInstance ca = (GKInstance) it.next();
                GKInstance catalyst = (GKInstance) ca.getAttributeValue(ReactomeJavaConstants.physicalEntity);
                if (catalyst != null && isChemical(catalyst))
                    return true;
            }
        }
        Collection<GKInstance> regulations = reaction.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null && regulations.size() > 0) {
            for (Iterator it = regulations.iterator(); it.hasNext();) {
                GKInstance regulation = (GKInstance) it.next();
                GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
                // Only take physical entity
                if (regulator != null && isChemical(regulator))
                    return true;
            }
        }
        return false;
    }
    
}
