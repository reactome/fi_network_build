/*
 * Created on Mar 21, 2012
 *
 */
package org.reactome.fi;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
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
    private FileUtility fu = new FileUtility();
    
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
            System.out.println("Source: " + source);
            Map<GKInstance, Set<String>> topic2Ids = analyzer.grepIDsFromTopics();
            for (Iterator<GKInstance> it = topic2Ids.keySet().iterator(); 
                 it.hasNext();) {
                GKInstance topic = it.next();
                Set<String> ids = topic2Ids.get(topic);
                System.out.println(topic + "(" + source + ")");
                System.out.println("Before filtering: " + ids.size());
                ids = filters.filter(ids);
                System.out.println("After filtering: " + ids.size());
                if (ids.size() < 10) {
                    // Remove these small pathways
                    System.out.println("Size is too small: " + ids.size());
                    continue; 
                }
                for (String id : ids) {
                    Set<String> topicNames = id2Topics.get(id);
                    if (topicNames == null) {
                        topicNames = new HashSet<String>();
                        id2Topics.put(id, topicNames);
                    }
                    topicNames.add(topic.getDisplayName() + "(" + source + ")");
                }
            }
        }
        // Output this map
        FileUtility fu = new FileUtility();
        fu.saveSetMap(id2Topics,
                      FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
//        fu.saveSetMapInSort(id2Topics, "Test.txt");
//        id2Topics = fu.loadSetMap(R3Constants.PROTEIN_ID_TO_TOPIC);
//        fu.saveSetMapInSort(id2Topics, "Test1.txt");
    }
    
    /**
     * This method is used to generate a map from protein name to topics based on another file
     * from protein UniProt ids to topics.
     * @throws Exception
     */
    @Test
    public void generateNameToTopicMap() throws Exception {
        Map<String, Set<String>> idToTopics = fu.loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
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
                      FIConfiguration.getConfiguration().get("GENE_TO_TOPIC"));
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
    
}
