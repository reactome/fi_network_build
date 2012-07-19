package org.reactome.data;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * Converter that converts (TF -> Target files) downloaded from the ENCODE project from
 * Mark Geistein's group (http://archive.gersteinlab.org/proj/encodenets/) to Reactome 
 * Curator Tool projects.
 *
 * @author Adrian Duong
 */
public class EncodeTFTargetToReactomeConverter {
	private XMLFileAdaptor fileAdaptor;
	private Map<String, GKInstance> nameToEntityMap;
	// private Map<String, GKInstance> interactionMap;
	private final GKInstance human;

	public EncodeTFTargetToReactomeConverter() throws Exception { // should this throw an exception?
		fileAdaptor = new XMLFileAdaptor();
		nameToEntityMap = new HashMap<String, GKInstance>();
		// interactionMap = new HashMap<String, GKInstance>();

		human = fileAdaptor.createNewInstance(ReactomeJavaConstants.Species);
		human.setAttributeValue(ReactomeJavaConstants.name, "Homo sapiens");
	}

    @Test
	public void convert() throws Exception {
		convert(FIConfiguration.getConfiguration().get("ENCODE_TFF_FILE"),
				FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE"));
		attachEvidences();
	}
    
    /**
     * Attach some evidences for the converted ENCODE TF/target interactions to be used later on.
     * Currently evidences have been attached in the definition slot.
     * @throws Exception
     */
    private void attachEvidences() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE"));
        Collection<GKInstance> interactions = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
        attachEvidences(interactions);
        fileAdaptor.save();
    }
    
    @Test
    public void testAttachEvidences() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE"));
        Collection<GKInstance> interactions = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
        attachEvidences(interactions);
        fileAdaptor.save("tmp/encodedTf_target_withEvidences.rtpj");
    }
    
    /**
     * Use this method to attach supporting evidneces to ENCODE TargettedInteractions that
     * have been loaded into the database.
     * @throws Exception
     */
    @Test
    public void attachEvidencesInDB() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost", 
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"), 
                                            FIConfiguration.getConfiguration().get("DB_USER"), 
                                            FIConfiguration.getConfiguration().get("DB_PWD"));
        // Get the encode database instance first
        Collection<GKInstance> c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceDatabase, 
                                                                ReactomeJavaConstants._displayName, 
                                                                "=",
                                                                "ENCODE");
        // There should be only one ENCODE DB
        GKInstance encode = c.iterator().next();
        c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.TargettedInteraction,
                                         ReactomeJavaConstants.dataSource, 
                                         "=",
                                         encode);
        System.out.println("Total ENCODE TargettedInteractions: " + c.size());
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.TargettedInteraction);
        dba.loadInstanceAttributeValues(c, cls.getAttribute(ReactomeJavaConstants.definition));
        dba.loadInstanceAttributeValues(c, cls.getAttribute(ReactomeJavaConstants.target));
        dba.loadInstanceAttributeValues(c, cls.getAttribute(ReactomeJavaConstants.factor));
        attachEvidences(c);
        try {
            dba.startTransaction();
            for (GKInstance inst : c) {
                String defintion = (String) inst.getAttributeValue(ReactomeJavaConstants.definition);
                if (defintion != null && defintion.contains("supported by")) {
                    System.out.println("Update " + inst);
                    dba.updateInstanceAttribute(inst, ReactomeJavaConstants.definition);
                }
            }
            dba.commit();
        }
        catch(Exception e) {
            dba.rollback();
            e.printStackTrace();
        }
    }
    
    private void attachEvidences(Collection<GKInstance> interactions) throws Exception {
        // Check for  co-expression
        MicroarrayDataAnalyzer microarray = new MicroarrayDataAnalyzer();
        Set<String> coexpressions = microarray.loadCoExpFromPavlidis();
        coexpressions.addAll(microarray.loadCoExpFromPrietoCarlos());
//        System.out.println("Total gene expression pairs: " + coexpressions.size());
        // Check for BP sharing
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> goTerms = goAnalyzer.loadProteinToGOBPTerms();
        
        for (GKInstance interaction : interactions) {
            GKInstance factor = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.factor);
            GKInstance target = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.target);
            // Note: EWAS instances are used always
            GKInstance factorRef = (GKInstance) factor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            GKInstance targetRef = (GKInstance) target.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (factorRef == null || targetRef == null)
                continue;
            String factorId = (String) factorRef.getAttributeValue(ReactomeJavaConstants.identifier);
            String targetId = (String) targetRef.getAttributeValue(ReactomeJavaConstants.identifier);
            if (factorId == null || targetId == null)
                continue;
            int compare = factorId.compareTo(targetId);
            String intInIds = null;
            if (compare > 0) {
                intInIds = targetId + "\t" + factorId;
            }
            else if (compare < 0)
                intInIds = factorId + "\t" + targetId;
            if (intInIds == null)
                continue;
            boolean isCoexpressed = coexpressions.contains(intInIds);
            boolean isBpShared = goAnalyzer.isTermShared(intInIds, goTerms);
            String definition = (String) interaction.getAttributeValue(ReactomeJavaConstants.definition);
            if (definition == null)
                definition = "";
            else
                definition = definition + "; ";
            if (isCoexpressed && isBpShared) 
                definition = definition + "supported by co-expression and GO BP sharing";
            else if (isCoexpressed)
                definition = definition + "supported by co-expression";
            else if (isBpShared)
                definition = definition + "supported by GO BP sharing";
            if (isCoexpressed || isBpShared)
                interaction.setAttributeValue(ReactomeJavaConstants.definition, definition);
        }
    }
    
    /**
     * Check TF/Target interactions with gene expressions
     * @throws Exception
     */
    @Test
    public void checkWithEvidences() throws Exception {
        MicroarrayDataAnalyzer microarray = new MicroarrayDataAnalyzer();
        Set<String> coexpressions = microarray.loadCoExpFromPavlidis();
        coexpressions.addAll(microarray.loadCoExpFromPrietoCarlos());
//        System.out.println("Total gene expression pairs: " + coexpressions.size());
        // Load the converted project
        XMLFileAdaptor fileAdpator = new XMLFileAdaptor();
        fileAdpator.setSource(FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE"));
        Collection<GKInstance> interactions = fileAdpator.fetchInstancesByClass(ReactomeJavaConstants.TargettedInteraction);
        System.out.println("Total TF/Target interactions: " + interactions.size());
        Set<String> interactionsInIds = new HashSet<String>();
        for (GKInstance interaction : interactions) {
            GKInstance factor = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.factor);
            GKInstance target = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.target);
            // Note: EWAS instances are used always
            GKInstance factorRef = (GKInstance) factor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            GKInstance targetRef = (GKInstance) target.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (factorRef == null || targetRef == null)
                continue;
            String factorId = (String) factorRef.getAttributeValue(ReactomeJavaConstants.identifier);
            String targetId = (String) targetRef.getAttributeValue(ReactomeJavaConstants.identifier);
            if (factorId == null || targetId == null)
                continue;
            int compare = factorId.compareTo(targetId);
            String intInIds = null;
            if (compare > 0) {
                intInIds = targetId + "\t" + factorId;
            }
            else if (compare < 0)
                intInIds = factorId + "\t" + targetId;
            if (intInIds == null)
                continue;
            interactionsInIds.add(intInIds);
        }
        System.out.println("Total interactions in ids: " + interactionsInIds.size());
        Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(interactionsInIds);
        System.out.println("Total proteins: " + proteins.size());
        
        Set<String> coexpShared = InteractionUtilities.getShared(coexpressions, interactionsInIds);
        System.out.println("Shared with co-expression: " + coexpShared.size());
        proteins = InteractionUtilities.grepIDsFromInteractions(coexpShared);
        System.out.println("Total proteins: " + proteins.size());
        
        // Check with GO BP
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> goTerms = goAnalyzer.loadProteinToGOBPTerms();
        Set<String> goBpShared = new HashSet<String>();
        for (String intInIds : interactionsInIds) {
            boolean isShared = goAnalyzer.isTermShared(intInIds, goTerms);
            if (isShared)
                goBpShared.add(intInIds);
        }
        System.out.println("Shared GO BP terms: " + goBpShared.size());
        proteins = InteractionUtilities.grepIDsFromInteractions(goBpShared);
        System.out.println("Total proteins: " + proteins.size());
        
        Set<String> bothShared = InteractionUtilities.getShared(coexpShared, goBpShared);
        System.out.println("Both shared: " + bothShared.size());
        proteins = InteractionUtilities.grepIDsFromInteractions(bothShared);
        System.out.println("Total proteins: " + proteins.size());
        
        // Merged two sharing
        Set<String> merged = new HashSet<String>(goBpShared);
        merged.addAll(coexpShared);
        System.out.println("Merged: " + merged.size());
        proteins = InteractionUtilities.grepIDsFromInteractions(merged);
        System.out.println("Total proteins: " + proteins.size());
    }

	private void convert(String tftargetFileName, String destFileName) throws Exception {
		MySQLAdaptor db = new MySQLAdaptor("localhost",
										   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
										   FIConfiguration.getConfiguration().get("DB_USER"),
										   FIConfiguration.getConfiguration().get("DB_PWD"),
										   3306);

		FileUtility fu = new FileUtility();
		Map<String, Set<String>> tfToTargetsMap = fu.loadSetMap(tftargetFileName);

		// make GKInstances of the interactions

		for (String tfName : tfToTargetsMap.keySet()) {
            GKInstance factor = nameToEntityMap.get(tfName);
            if (factor == null) {
                factor = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
                factor.addAttributeValue(ReactomeJavaConstants.name, tfName);
                factor.setDisplayName(tfName);

                nameToEntityMap.put(tfName, factor);
            }

			Set<String> targetNames = tfToTargetsMap.get(tfName);
			for (String targetTypeName : targetNames) {
                int tabIndex = targetTypeName.indexOf("\t");
                if (tabIndex < 0) {
                    continue;
                }

                String targetType = targetTypeName.substring(0, tabIndex);
                String targetName = targetTypeName.substring(tabIndex+1);

				GKInstance target = nameToEntityMap.get(targetName);
                if (target == null) {
                    target = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
                    target.addAttributeValue(ReactomeJavaConstants.name, targetName);
                    target.setDisplayName(targetName);

                    nameToEntityMap.put(targetName, target);
                }

				String interactionKey = factor.getDBID() + "-" + target.getDBID();
				GKInstance interaction = fileAdaptor.createNewInstance(ReactomeJavaConstants.TargettedInteraction);
				interaction.setAttributeValue(ReactomeJavaConstants.factor, factor);
				interaction.setAttributeValue(ReactomeJavaConstants.target, target);
				interaction.setAttributeValue(ReactomeJavaConstants.species, human);
				interaction.setDisplayName(factor.getDisplayName() + "-" + target.getDisplayName());
                interaction.setAttributeValue(ReactomeJavaConstants.definition,
                                              "ENCODE " + targetType + " TF/target interaction");
			}
		}

        EncodeTFTargetToReactomePostProcessor postProcessor = new EncodeTFTargetToReactomePostProcessor();
        postProcessor.postProcess(db,
                                  fileAdaptor);

		fileAdaptor.save(destFileName);
	}
}