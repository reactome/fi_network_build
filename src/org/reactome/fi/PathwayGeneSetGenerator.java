/*
 * Created on Mar 21, 2012
 *
 */
package org.reactome.fi;

import static org.gk.model.ReactomeJavaConstants.Pathway;
import static org.gk.model.ReactomeJavaConstants.dataSource;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
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

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.MySQLAdaptor.QueryRequest;
import org.gk.render.Renderable;
import org.gk.render.RenderablePathway;
import org.junit.Test;
import org.reactome.data.ProteinIdFilters;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.ReactomeUtilities;
import org.reactome.gsea.ReactomeToMsigDBExport;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.kegg.KeggAnalyzer;

/**
 * This class is used to generate gene sets from pathways. Methods in this class
 * are copied from original TopicAnalyzer in the caBigR3 project.
 * @author gwu
 *
 */
public class PathwayGeneSetGenerator {
    private static Logger logger = Logger.getLogger(PathwayGeneSetGenerator.class);
    private FileUtility fu = new FileUtility();
    private final int MINIMUM_PATHWAY_SIZE = new Integer(FIConfiguration.getConfiguration().get("MINIMUM_PAHTWAY_SIZE"));
    
    public PathwayGeneSetGenerator() {
    }
    
    @Test
    public void generateProteinNameToPathwayMap() throws Exception {
        generateIdToTopicMap();
        generateNameToTopicMap();
        fixKeggPathwayGeneSets();
    }
    
    @Test
    public void generateIdToTopicMap() throws Exception {
        // Get a map from protein ids to topics (pathways)
        Map<String, Set<String>> id2Topics = new HashMap<String, Set<String>>();
        // Use all four databases
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        ProteinIdFilters filters = new ProteinIdFilters();
        for (ReactomeAnalyzer analyzer : analyzers) {
            String source = ReactomeAnalyzer.getSourceLetter(analyzer);
            logger.info("Source: " + source);
            Map<GKInstance, Set<String>> topic2Ids = analyzer.grepIDsFromTopics();
            generateIdToTopicMap(id2Topics, 
                                 filters,
                                 topic2Ids,
                                 source);
        }
        // Output this map
        fu.saveSetMap(id2Topics,
                      FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
//        fu.saveSetMapInSort(id2Topics, "Test.txt");
//        id2Topics = fu.loadSetMap(R3Constants.PROTEIN_ID_TO_TOPIC);
//        fu.saveSetMapInSort(id2Topics, "Test1.txt");
    }
    
    @Test
    public void countProteinsInPathways() throws Exception {
        Map<String, Set<String>> id2Topics = fu.loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
        System.out.println("Total ids in pathways: " + id2Topics.size());
        Map<String, Set<String>> topicToIds = InteractionUtilities.switchKeyValues(id2Topics);
        Set<String> reactomeIds = new HashSet<String>();
        Set<String> keggIds = new HashSet<String>();
        for (String topic : topicToIds.keySet()) {
            if (topic.endsWith("(R)"))
                reactomeIds.addAll(topicToIds.get(topic));
            else if (topic.endsWith("(K)"))
                keggIds.addAll(topicToIds.get(topic));
        }
        System.out.println("Ids in Reactome: " + reactomeIds.size());
        System.out.println("Ids in KEGG: " + keggIds.size());
        Set<String> shared = InteractionUtilities.getShared(reactomeIds, keggIds);
        System.out.println("Shared: " + shared.size());
        
        // Filter based on SwissProt
        Set<String> swissProtIds = new UniProtAnalyzer().loadSwissProtIds();
        int total = swissProtIds.size();
        System.out.println("\nTotal SwissProt: " + total);
        // Check SwissIds in Reactome and KEGG
        id2Topics.keySet().retainAll(swissProtIds);
        System.out.println("Ids from pathways in SwissProt: " + id2Topics.size());
        
        reactomeIds.retainAll(swissProtIds);
        keggIds.retainAll(swissProtIds);
        shared = InteractionUtilities.getShared(reactomeIds, keggIds);
        System.out.println("SwissProt ids in Reactome: " + reactomeIds.size());
        System.out.println("SwissProt ids in KEgg: " + keggIds.size());
        System.out.println("Shared SwissProt ids: " + shared.size());
    }

    private void generateIdToTopicMap(Map<String, Set<String>> id2Topics,
                                      ProteinIdFilters filters,
                                      Map<GKInstance, Set<String>> topic2Ids,
                                      String source) throws Exception {
        for (Iterator<GKInstance> it = topic2Ids.keySet().iterator(); 
             it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topic2Ids.get(topic);
            logger.info(topic + "(" + source + ")");
//            logger.info("Before filtering: " + ids.size());
            ids = filters.filter(ids);
//            logger.info("After filtering: " + ids.size());
            if (ids.size() < MINIMUM_PATHWAY_SIZE) { 
                // Remove these small pathways
                logger.info("Size is too small: " + ids.size());
                continue; 
            }
            for (String id : ids) {
                Set<String> topicNames = id2Topics.get(id);
                if (topicNames == null) {
                    topicNames = new HashSet<String>();
                    id2Topics.put(id, topicNames);
                }
                if (source == null)
                    topicNames.add(topic.getDisplayName());
                else
                    topicNames.add(topic.getDisplayName() + "(" + source + ")");
            }
        }
    }
    
    /**
     * This method is used to generate a map from protein name to topics based on another file
     * from protein UniProt ids to topics.
     * @throws Exception
     */
    @Test
    public void generateNameToTopicMap() throws Exception {
        Map<String, Set<String>> idToTopics = fu.loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
        String fileName = FIConfiguration.getConfiguration().get("GENE_TO_TOPIC");
        
        generateNameToTopicMap(idToTopics, fileName);
    }

    private void generateNameToTopicMap(Map<String, Set<String>> idToTopics,
                                        String fileName) throws Exception, IOException {
        HibernateFIReader hibernateAnalyzer = new HibernateFIReader();
        Map<String, String> idToNames = hibernateAnalyzer.generateAccessionToProteinNames();
        Map<String, Set<String>> nameToTopics = new HashMap<String, Set<String>>();
        for (String id : idToTopics.keySet()) {
            Set<String> topics = idToTopics.get(id);
            String name = idToNames.get(id);
            if (name == null)
                continue;
            Set<String> nameTopics = nameToTopics.get(name);
            if (nameTopics == null) {
                nameTopics = new HashSet<String>();
                nameToTopics.put(name, nameTopics);
            }
            nameTopics.addAll(topics);
        }
//        fu.saveSetMap(nameToTopics, 
//                      R3Constants.RESULT_DIR + "ProteinNameToTopics051109.txt");
//        fu.saveSetMap(nameToTopics, 
//                      R3Constants.RESULT_DIR + "ProteinNameToTopics080410.txt");
        fu.saveSetMap(nameToTopics,
                      fileName);
    }
    
    /**
     * KEGG pathways contain genes having interactions only. However, some genes in KEGG pathways
     * have no interactions with others. These genes should be included in kegg gene sets too.
     * @throws Exception
     */
    @Test
    public void fixKeggPathwayGeneSets() throws Exception {
        String srcFileName = FIConfiguration.getConfiguration().get("GENE_TO_TOPIC");
        File srcFile = new File(srcFileName);
        // Make a temp copy of the source file name
        String tmpFileName = srcFileName + ".tmp";
        File tempFile = new File(tmpFileName);
        srcFile.renameTo(tempFile); // Make a temp file
        KeggAnalyzer keggAnalyzer = new KeggAnalyzer();
        keggAnalyzer.augmentGeneNameToPathwayMap(tempFile.getAbsolutePath(),
                                                 srcFileName);
        tempFile.delete();
    }
    
    /**
     * Some pathway diagrams have sub-pathways listed. This method is used to grep
     * these sub-pathways based on the list of pathways generated in method {@link 
     * generateReactomePathwayListBasedOnDiagrams() generateReactomePathwayListBasedOnDiagrams}.
     * The example is this:  https://reactome.org/PathwayBrowser/#/R-HSA-73857&PATH=R-HSA-74160
     * (see sub-pathway Generic Transcription Pathway) (release 63). The generated list of
     * pathways are used for batch PGM or BN analysis.
     * @throws Exception
     */
    @Test
    public void boostPathwayListForModeling() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        List<Long> dbIds = Files.lines(Paths.get(fileName))
                .map(line -> line.split("\t")[0])
                .map(text -> new Long(text))
                .collect(Collectors.toList());
        logger.info("Total DB_IDs: " + dbIds.size());

        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = (MySQLAdaptor) reactomeAnalyzer.getMySQLAdaptor();

        Set<GKInstance> pathways = new HashSet<>();
        for (Long dbId : dbIds) {
            GKInstance pathway = dba.fetchInstance(dbId);
            boostPathwayListForModeling(pathways, pathway, dba);
        }
        
        // Some normal pathway diagrams have disease pathways drawn as potential downstream 
        // pathways. Using the following statements to exclude them.
        GKInstance disease = dba.fetchInstance(1643685L); // The topmost Disease pathway
        Set<GKInstance> diseasePathways = InstanceUtilities.getContainedEvents(disease);
        diseasePathways.add(disease);
        diseasePathways.retainAll(pathways);
        logger.info("Pathways in disease: " + diseasePathways.size());
        diseasePathways.forEach(System.out::println);
        
        logger.info("Total pathways after boosting: " + pathways.size());
        pathways.removeAll(diseasePathways);
        logger.info("Total pathways after removing disease pathways: " + pathways.size());
        String[] tokens = fileName.split("\\.");
        String newFileName = tokens[0] + "_ForModeling." + tokens[1];
        outputPathways(pathways, newFileName);
    }
    
    @SuppressWarnings("unchecked")
    private void boostPathwayListForModeling(Set<GKInstance> pathways,
                                             GKInstance pathway,
                                             MySQLAdaptor dba) throws Exception {
        if (pathways.contains(pathway))
            return;
        logger.info("Check " + pathway);
        Collection<GKInstance> diagrams = dba.fetchInstanceByAttribute(ReactomeJavaConstants.PathwayDiagram,
                ReactomeJavaConstants.representedPathway, 
                "=",
                pathway);
        if (diagrams == null || diagrams.size() != 1) {
            String message = pathway + " doesn't have diagram or has more than one diagram: " + diagrams;
            logger.error(message);
            // A skip for release 79
            if (pathway.getDBID().equals(177128L) ||
                pathway.getDBID().equals(622312L))
            	return;
            throw new IllegalStateException(message);
        }
        // Need to check its contained pathways
        GKInstance diagram = diagrams.iterator().next();
        DiagramGKBReader reader = new DiagramGKBReader();
        RenderablePathway renderableDiagram = reader.openDiagram(diagram);
        boolean hasEntity = false;
        // Use this set to avoid circular link: PathwayA <-> PathwayB
        Set<GKInstance> containedPathways = new HashSet<>();
        for (Object obj : renderableDiagram.getComponents()) {
            Renderable r = (Renderable) obj;
            if (r.getReactomeId() == null)
                continue;
            GKInstance inst = dba.fetchInstance(r.getReactomeId());
            if (inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                hasEntity = true;
            }
            else if (inst.getSchemClass().isa(ReactomeJavaConstants.Pathway)) {
                containedPathways.add(inst);
            }
        }
        if (hasEntity)
            pathways.add(pathway);
        for (GKInstance containedPathway : containedPathways)
            boostPathwayListForModeling(pathways, containedPathway, dba);
    }
    
    /**
     * Using this method to generate a list of pathways from Reactome. Each pathway should
     * have a corresponding pathway diagrams with entities fully laid-out.
     * @throws Exception
     */
    @Test
    @SuppressWarnings("unchecked")
    public void generateReactomePathwayListBasedOnDiagrams() throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = (MySQLAdaptor) reactomeAnalyzer.getMySQLAdaptor();
        // Start from the top level pathways so that we can exclude pathways under the disease topics
        Collection<GKInstance> frontPages = dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        if (frontPages == null || frontPages.size() != 1) {
            logger.error("No FrontPage or more than one FrontPage instance!");
            throw new IllegalStateException("No FrontPage or more than one FrontPage instance!");
        }
        GKInstance frontPage = (GKInstance) frontPages.iterator().next(); // There should be at least one FrontPageItem
        List<GKInstance> frontPageItems = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        Set<GKInstance> pathways = new HashSet<GKInstance>();
        for (GKInstance item : frontPageItems) {
            // Escape the disease topic since they should be covered by associated normal pathways
            if (item.getDisplayName().equalsIgnoreCase("Disease"))
                continue;
            // Exclude non-human pathways
            GKInstance species = (GKInstance) item.getAttributeValue(ReactomeJavaConstants.species);
            if (species == null ||
                !species.getDBID().equals(ReactomeUtilities.HOMO_SAPIENS_DB_ID))
                continue;
            grepPathwaysWithDiagrams(item, pathways, dba);
        }
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        outputPathways(pathways, fileName);
    }

    private void outputPathways(Set<GKInstance> pathways, String fileName) throws IOException {
        // Want to sort it before output
        List<GKInstance> list = new ArrayList<GKInstance>(pathways);
        InstanceUtilities.sortInstances(list);
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        for (GKInstance pathway : list)
            fu.printLine(pathway.getDBID() + "\t" + pathway.getDisplayName());
        fu.close();
        logger.info("Total Pathways from Reactome: " + list.size());
    }
    
    @SuppressWarnings("unchecked")
    private void grepPathwaysWithDiagrams(GKInstance pathway, 
                                          Set<GKInstance> pathways,
                                          MySQLAdaptor dba) throws Exception {
        if (pathways.contains(pathway))
            return; // Just in case this has been checked before since a pathway can be contained in several different places.
        // Just in case something is not right
        if (!pathway.getSchemClass().isa(ReactomeJavaConstants.Pathway)) {
            logger.error(pathway + " is not a pathway!");
            throw new IllegalArgumentException(pathway + " is not a pathway!");
        }
        Collection<GKInstance> diagrams = dba.fetchInstanceByAttribute(ReactomeJavaConstants.PathwayDiagram,
                                                                       ReactomeJavaConstants.representedPathway, 
                                                                       "=",
                                                                       pathway);
        if (diagrams == null || diagrams.size() != 1) {
            String message = pathway + " doesn't have diagram or has more than one diagram: " + diagrams;
            logger.error(message);
            throw new IllegalStateException(message);
        }
        GKInstance diagram = diagrams.iterator().next();
        DiagramGKBReader reader = new DiagramGKBReader();
        RenderablePathway renderableDiagram = reader.openDiagram(diagram);
        for (Object obj : renderableDiagram.getComponents()) {
            Renderable r = (Renderable) obj;
            if (r.getReactomeId() == null)
                continue;
            GKInstance inst = dba.fetchInstance(r.getReactomeId());
            if (inst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
                pathways.add(pathway);
                break;
            }
        }
        if (pathways.contains(pathway))
            return; 
        // Need to check its subpathways
        List<GKInstance> hasEvent = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        for (GKInstance comp : hasEvent)
            grepPathwaysWithDiagrams(comp, pathways, dba);
    }
    
    private Collection<GKInstance> getPathwaysForExport(ReactomeAnalyzer reactomeAnalyzer,
                                                        Long speciesId) throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) reactomeAnalyzer.getMySQLAdaptor();
        GKInstance species = dba.fetchInstance(speciesId);
        List<QueryRequest> query = new ArrayList<MySQLAdaptor.QueryRequest>();
        query.add(dba.createAttributeQueryRequest(Pathway,
                                                  ReactomeJavaConstants.species, 
                                                  "=",
                                                  species));
        // So that it can be used for the normal release database
        if (dba.getSchema().getClassByName(Pathway).isValidAttribute(dataSource)) {
            query.add(dba.createAttributeQueryRequest(Pathway, 
                                                      dataSource, 
                                                      "IS NULL", 
                                                      null));
        }
        @SuppressWarnings("unchecked")
        Collection<GKInstance> pathways = dba.fetchInstance(query);
        logger.info("Total " + species.getDisplayName() + " pathways in Reactome: " + pathways.size());
        // Did a fix on July 23, 2016 to get non-disease pathways only
        if (species.getDisplayName().equals("Homo sapiens"))
            pathways = getNonDiseasePathways(dba, pathways); // For mouse, we will just use all pathways since disease pathways are not inferred.
                                                             // However, this may change in the future.
        logger.info("After removing pathways in Disease: " + pathways.size());
        return pathways;
    }
    
    /**
     * This method is used to generate gene name to pathway map for all pathways in the Reactome
     * database. This method has been deprecated. Use generateHumanFiles() instead.
     */
    @Test
    @Deprecated
    public void generateReactomeGeneToPathwayMap() throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Long humanDbId = 48887L;
        Collection<GKInstance> pathways = getPathwaysForExport(reactomeAnalyzer, humanDbId);
        logger.info("generateReactomeGeneToPathwayMap(): total " + pathways.size() + " pathways.");
        // Get a map from protein ids to topics (pathways)
        Map<String, Set<String>> id2Topics = new HashMap<String, Set<String>>();
        ProteinIdFilters filters = new ProteinIdFilters();
        Map<GKInstance, Set<String>> topic2Ids = reactomeAnalyzer.grepIDsFromTopics(pathways);
        generateIdToTopicMap(id2Topics, 
                             filters,
                             topic2Ids,
                             null);
        // Output this map
        String fileName = FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_REACTOME_PATHWAYS");
        fu.saveSetMap(id2Topics,
                      fileName);
        
        // Convert UniProt ids to gene names
        fileName = FIConfiguration.getConfiguration().get("GENE_TO_REACTOME_PATHWAYS");
        generateNameToTopicMap(id2Topics, fileName);
    }
    
    /**
     * This method is used to generate mouse gene name to pathway map for all pathways in the Reactome
     * database. Both genes to pathways and the mouse gmt file are generated here.
     */
    @Test
    public void generateMouseFiles() throws Exception {
        Long mouseDbId = 48892L;
        String gmtFileName = FIConfiguration.getConfiguration().get("MOUSE_REACTOME_GMT_FILE_NAME");
        String geneToPathwayFileName = FIConfiguration.getConfiguration().get("MOUSE_GENE_TO_REACTOME_PATHWAYS");
        
        generatePathwayEnrichmentFiles(mouseDbId, gmtFileName, geneToPathwayFileName);
    }
    
    /**
     * This method is used to fulfill a request from an outside Reactome user.
     * @throws Exception
     */
    @Test
    public void generateBovineFiles() throws Exception {
        Long bovineDbId = 48898L; // Bos taurus
        String gmtFileName = "Bovine_Release75_01062020.gmt";
        String geneToPathwayFileName = "Bovine_Release75_Gene2Pathway_01062020.txt";
        
        generatePathwayEnrichmentFiles(bovineDbId, gmtFileName, geneToPathwayFileName);
    }
    
    /**
     * This method is used to generate human gene name to pathway map for all pathways in the Reactome
     * database. Both genes to pathways and the mouse gmt file are generated here.
     */
    @Test
    public void generateHumanFiles() throws Exception {
        Long humanDbId = 48887L;
        String gmtFileName = FIConfiguration.getConfiguration().get("REACTOME_GMT_FILE_NAME");
        String geneToPathwayFileName = FIConfiguration.getConfiguration().get("GENE_TO_REACTOME_PATHWAYS");
        
        generatePathwayEnrichmentFiles(humanDbId, gmtFileName, geneToPathwayFileName);
    }

    private void generatePathwayEnrichmentFiles(Long mouseDbId, String gmtFileName, String geneToPathwayFileName)
            throws Exception, FileNotFoundException, IOException {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Collection<GKInstance> pathways = getPathwaysForExport(reactomeAnalyzer, mouseDbId);
        logger.info("generateMouseReactomeGeneToPathwayMap(): total " + pathways.size() + " pathways.");
        
        ReactomeToMsigDBExport msigDbExport = new ReactomeToMsigDBExport();
        msigDbExport.setSpeciesId(mouseDbId);
        // Use gene names. Reactome release database have gene names filled for mouse ReferenceGeneProducts since release 72.
        msigDbExport.setUseUniProt(false); 
        msigDbExport.setSizeCutoff(MINIMUM_PATHWAY_SIZE);
        msigDbExport.exportInGMT(pathways, new FileOutputStream(gmtFileName));
        
        // Read back the GMT file so that we can generate gene to pathway mapping. This is a little bit weird though.
        // But it works.
        fu.setInput(gmtFileName);
        fu.setOutput(geneToPathwayFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            for (int i = 2; i < tokens.length; i++) {
                // The first token is pathway name and the second is the stable id. All others are gene names.
                fu.printLine(tokens[i] + "\t" + tokens[0]);
            }
        }
        fu.close();
    }
    
    // This is a fix performed on July 23, 2016 to remove pathways under Disease since they have
    // not listed in the hierarchy for Reactome pathway analysis
    // However, we cannot just remove pathways listed under Disease since normal pathways may be listed
    // there too. We don't want to remove them. So we do another away around.
    private Set<GKInstance> getNonDiseasePathways(MySQLAdaptor dba,
                                                  Collection<GKInstance> pathways) throws Exception {
        Set<GKInstance> nonDiseasePathways = new HashSet<GKInstance>();
        Collection<GKInstance> frontPages = dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        for (GKInstance frontPage : frontPages) {
            List<GKInstance> items = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
            for (GKInstance item : items) {
                if (item.getDisplayName().equals("Disease"))
                    continue; // Escape disease
                if (pathways.contains(item)) {
                    Set<GKInstance> events = InstanceUtilities.getContainedEvents(item);
                    for (GKInstance event : events)
                        if (event.getSchemClass().isa(ReactomeJavaConstants.Pathway))
                            nonDiseasePathways.add(event);
                }
            }
        }
        return nonDiseasePathways;
    }
    
}
