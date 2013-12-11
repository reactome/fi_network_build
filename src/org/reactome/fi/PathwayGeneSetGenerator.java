/*
 * Created on Mar 21, 2012
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.security.auth.login.Configuration;

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
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
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
            grepPathwaysWithDiagrams(item, pathways, dba);
        }
        // Want to sort it before output
        List<GKInstance> list = new ArrayList<GKInstance>(pathways);
        InstanceUtilities.sortInstances(list);
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
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
            }
        }
        if (pathways.contains(pathway))
            return; 
        // Need to check its subpathways
        List<GKInstance> hasEvent = pathway.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
        for (GKInstance comp : hasEvent)
            grepPathwaysWithDiagrams(comp, pathways, dba);
    }
    
    /**
     * This method is used to generate gene name to pathway map for all pathways in the Reactome
     * database.
     */
    @Test
    public void generateReactomeGeneToPathwayMap() throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = (MySQLAdaptor) reactomeAnalyzer.getMySQLAdaptor();
        // Got all human pathways
        // This is humna DB_ID
        Long humanDbId = 48887L;
        GKInstance human = dba.fetchInstance(humanDbId);
        List<QueryRequest> query = new ArrayList<MySQLAdaptor.QueryRequest>();
        query.add(dba.createAttributeQueryRequest(ReactomeJavaConstants.Pathway,
                                                  ReactomeJavaConstants.species, 
                                                  "=",
                                                  human));
        query.add(dba.createAttributeQueryRequest(ReactomeJavaConstants.Pathway, 
                                                  ReactomeJavaConstants.dataSource, 
                                                  "IS NULL", 
                                                  null));
        @SuppressWarnings("unchecked")
        Collection<GKInstance> pathways = dba.fetchInstance(query);
        logger.info("Total human pathways in Reactome: " + pathways.size());
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
}
