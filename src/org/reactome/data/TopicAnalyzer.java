/*
 * Created on Mar 28, 2007
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.util.*;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.fi.FIFileAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.hibernate.HibernateFIReader;

/**
 * This class is used to handle pathway or topic related processing.
 * @author guanming
 *
 */
public class TopicAnalyzer {
    // Want to cache this map since the loading is slow
    private final String CLUSTER_FILE = FIConfiguration.getConfiguration().get("RESULT_DIR") + "PathwayClusters070425.txt";
    private FileUtility fu = new FileUtility();
    
//    @Test
//    public void extractPathwayCentrality() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "Pathways_Reactome_28_plus_i_betweeness_centrality.txt";
//        fu.setInput(fileName);
//        String line = null;
//        int totalLine = 0;
//        List<String> lines = new ArrayList<String>();
//        int totalPathways = 0;
//        while ((line = fu.readLine()) != null) {
//            if (line.startsWith("#")) {
//                if (lines.size() > 1) {
//                    for (String tmp : lines)
//                        System.out.println(tmp);
//                    totalPathways ++;
//                }
//                lines.clear();
//                lines.add(line);
//            }
//            else {
//                String[] tokens = line.split("\t");
//                double value = Double.parseDouble(tokens[2]);
//                if (value >= 0.20) {
//                    lines.add(line);
//                    totalLine ++;
//                }
//            }
//        }
//        fu.close();
//        System.out.println("Total genes: " + totalLine);
//        System.out.println("Total pathways: " + totalPathways);
//    }
//    
//    public Map<String, Map<String, Double>> loadPathwayToGeneCentralities() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "Pathways_Reactome_28_plus_i_betweeness_centrality.txt";
//        fu.setInput(fileName);
//        Map<String, Map<String, Double>> rtn = new HashMap<String, Map<String,Double>>();
//        String line = null;
//        Map<String, Double> geneToCentrality = new HashMap<String, Double>();
//        while ((line = fu.readLine()) != null) {
//            if (line.startsWith("#")) {
//                // Pathway
//                int index = line.lastIndexOf("\t");
//                String pathway = line.substring(1, index);
//                geneToCentrality = new HashMap<String, Double>();
//                rtn.put(pathway, geneToCentrality);
//            }
//            else {
//                String[] tokens = line.split("\t");
//                // Use normalization one
//                geneToCentrality.put(tokens[0], new Double(tokens[2]));
//            }
//        }
//        fu.close();
//        return rtn;
//    }
//    
//    
//    public Map<String, Set<String>> getClusteredTopicToIdMap(int hop) throws IOException {
//        Map<String, Set<String>> idToTopics = generateIdToClusteredTopicsWithHop(hop);
//        Map<String, Set<String>> clusteredTopicToIds = InteractionUtilities.switchKeyValues(idToTopics);
//        return clusteredTopicToIds;
//    }
//    
//    @Test
//    public void generatePathwayClusters() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Hierarchical_PI_CUT_OFF_50.txt";
//        fu.setInput(fileName);
//        String line = null;
//        // This is for output
//        FileUtility outputFu = new FileUtility();
//        outputFu.setOutput(CLUSTER_FILE);
//        List<String> annotations = loadClusterAnnotations();
//        int c = 0;
//        outputFu.printLine("Index\tClusterIndex\tAnnotations\tPathways");
//        StringBuilder builder = new StringBuilder();
//        while ((line = fu.readLine()) != null) {
//            //1 (0.5714285714285714): Antigen processing and presentation(K)(73), antigen processing and presentation(B)(7)
//            int pos1 = line.indexOf('(');
//            int pos2 = line.indexOf(')');
//            String index = line.substring(0, pos1).trim();
//            String clusterIndex = line.substring(pos1 + 1, pos2);
//            String pathways = line.substring(pos2 + 3);
//            String[] tokens = pathways.split("\\),");
//            builder.append(index).append("\t");
//            builder.append(clusterIndex).append("\t");
//            builder.append(annotations.get(c)).append("\t");
//            for (String token : tokens) {
//                token = token.trim();
//                // want to get rid of numbers
//                pos1 = token.lastIndexOf("(");
//                token = token.substring(0, pos1);
//                builder.append(token).append("\t");
//            }
//            outputFu.printLine(builder.toString());
//            builder.setLength(0);
//            c ++;
//        }
//        outputFu.close();
//        fu.close();
//    }
//    
//    public List<String> loadClusterAnnotations() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Hierarchical_PI_CUT_OFF_50_With_Label.txt";
//        fu.setInput(fileName);
//        String line = null;
//        List<String> annotations = new ArrayList<String>();
//        while ((line = fu.readLine()) != null) {
//            // 1 (0.5714285714285714) (Antigen processing and presentation): Antigen ...
//            int index1 = line.indexOf(")");
//            index1 = line.indexOf("(", index1) + 1;
//            int index2 = line.indexOf("):", index1);
//            String annotation = null;
//            if (index2 == -1) {
//                // No annotations. Use the pathways.
//                index1 = line.indexOf(":") + 2;
//                index2 = line.length();
//                annotation = line.substring(index1, index2);
//                annotation = removeSourceAndNumber(annotation);
//            }
//            else
//                annotation = line.substring(index1, index2);
//            annotations.add(annotation);
//        }
//        fu.close();
//        System.out.println("Size of list: " + annotations.size());
//        for (String annot : annotations)
//            System.out.println(annot);
//        Set<String> set = new HashSet<String>(annotations);
//        System.out.println("Size of set: " + set.size());
//        return annotations;
//    }
//    
//    private String removeSourceAndNumber(String line) {
//        List<String> pathways = new ArrayList<String>();
//        // Find the second last index of "("
//        String[] tokens = line.split(",");
//        for (String token : tokens) {
//            token = token.trim();
//            int index = token.lastIndexOf("(");
//            index = token.lastIndexOf('(', index - 1);
//            //System.out.println(token);
//            pathways.add(token.substring(0, index));
//        }
//        if (tokens.length == 1)
//            return pathways.get(0);
//        else {
//            StringBuilder builder = new StringBuilder();
//            for (Iterator<String> it = pathways.iterator(); it.hasNext();) {
//                builder.append(it.next());
//                if (it.hasNext())
//                    builder.append(", ");
//            }
//            return builder.toString();
//        }
//    }
//    
//    @Test
//    public void generateTopicToIdAndNumber() throws Exception {
//        // Topic names to a map of protein ids to numbers
//        Map<String, Map<String, Integer>> topicToIdNumber = new HashMap<String, Map<String, Integer>>();
//        // Use all databases
//        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
//        for (ReactomeAnalyzer analyzer : analyzers) {
//            Map<GKInstance, Map<String, Integer>> topicInstanceToIdNumber = analyzer.grepIDNumberFromTopics();
//            String source = ReactomeAnalyzer.getSourceLetter(analyzer);
//            for (Iterator<GKInstance> it = topicInstanceToIdNumber.keySet().iterator(); 
//                 it.hasNext();) {
//                GKInstance topic = it.next();
//                Map<String, Integer> idToNumber = topicInstanceToIdNumber.get(topic);
//                if (idToNumber.size() == 0)
//                    continue; // In case like this: Focal Adhesion appears twice in KEGG. One is just an empty pathway.
//                              // This should be excluded.
//                String name = topic.getDisplayName() + "(" + source + ")";
//                topicToIdNumber.put(name, topicInstanceToIdNumber.get(topic));
//            }
//        }
//        // Export three columns
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "TopicIDNumber.txt";
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        fu.printLine("Topic\tID\tNumber");
//        for (Iterator<String> it = topicToIdNumber.keySet().iterator(); it.hasNext();) {
//            String name = it.next();
//            Map<String, Integer> map = topicToIdNumber.get(name);
//            for (Iterator<String> it1 = map.keySet().iterator(); it1.hasNext();) {
//                String id = it1.next();
//                Integer number = map.get(id);
//                fu.printLine(name + "\t" + id + "\t" + number);
//            }
//        }
//        fu.close();
//    }
//    
//    public Map<String, Map<String, Integer>> loadTopicToIDNumber() throws IOException {
//        Map<String, Map<String, Integer>> topicToIdNumber = new HashMap<String, Map<String, Integer>>();
//        FileUtility fu = new FileUtility();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "TopicIDNumber.txt";
//        fu.setInput(fileName);
//        // Title line
//        String line = fu.readLine();
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            Map<String, Integer> idToNumber = topicToIdNumber.get(tokens[0]);
//            if (idToNumber == null) {
//                idToNumber = new HashMap<String, Integer>();
//                topicToIdNumber.put(tokens[0],
//                                    idToNumber);
//            }
//            idToNumber.put(tokens[1],
//                           Integer.parseInt(tokens[2]));
//        }
//        fu.close();
//        return topicToIdNumber;
//    }
//    
    /**
     * This method is used to generate a file mapping from protein ids to pathways.
     * @throws Exception
     */
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
        fu.saveSetMap(id2Topics, FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC"));
//        fu.saveSetMapInSort(id2Topics, "Test.txt");
//        id2Topics = fu.loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC);
//        fu.saveSetMapInSort(id2Topics, "Test1.txt");
    }
//    
//    /**
//     * This method is used to generate a map from protein name to topics based on another file
//     * from protein UniProt ids to topics.
//     * @throws Exception
//     */
//    @Test
//    public void generateNameToTopicMap() throws Exception {
//        Map<String, Set<String>> idToTopics = fu.loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC);
//        HibernateFIReader hibernateAnalyzer = new HibernateFIReader();
//        Map<String, String> idToNames = hibernateAnalyzer.generateAccessionToProteinNames();
//        Map<String, Set<String>> nameToTopics = new HashMap<String, Set<String>>();
//        for (String id : idToTopics.keySet()) {
//            Set<String> topics = idToTopics.get(id);
//            String name = idToNames.get(id);
//            if (name == null)
//                continue;
//            Set<String> nameTopics = nameToTopics.get(name);
//            if (nameTopics == null) {
//                nameTopics = new HashSet<String>();
//                nameToTopics.put(name, nameTopics);
//            }
//            nameTopics.addAll(topics);
//        }
////        fu.saveSetMap(nameToTopics, 
////                      FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics051109.txt");
////        fu.saveSetMap(nameToTopics, 
////                      FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics080410.txt");
//        fu.saveSetMap(nameToTopics,
//                      FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics101110.txt");
//    }
//    
//    public Map<String, Set<String>> loadNameToTopicsMap() throws IOException {
//        //return fu.loadSetMap(FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics100608.txt");
//        return fu.loadSetMap(FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics051109.txt");
//    }
//    
//    public Map<String, Set<String>> loadNameToPathwayClustersMap() throws IOException {
//        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Pathways_Reactome_28_plus_i_clusters_cutoff_2.gmt";
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Pathways_Reactome_28_plus_i_betweeness_clusters.gmt";
//        fu.setInput(fileName);
//        String line = fu.readLine();
//        Map<String, Set<String>> nameToClusters = new HashMap<String, Set<String>>();
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            String pathway = tokens[0];
//            for (int i = 2; i < tokens.length; i++) {
//                String gene = tokens[i];
//                Set<String> pathways = nameToClusters.get(gene);
//                if (pathways == null) {
//                    pathways = new HashSet<String>();
//                    nameToClusters.put(gene, pathways);
//                }
//                pathways.add(pathway);
//            }
//        }
//        fu.close();
//        return nameToClusters;
//    }
//    
//    /**
//     * Load pathway clusters to genes. Pathway clusters are generated by method
//     * TCGAOvarianCancerAnalyzer.generatePathwayClusters().
//     * @return
//     * @throws IOException
//     */
//    public Map<String, Set<String>> loadPathwayClusterToGenes() throws IOException {
//        Map<String, Set<String>> pathwayToGenes = new HashMap<String, Set<String>>();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Pathways_Reactome_28_plus_i_clusters_cutoff_2.gmt";
//        fu.setInput(fileName);
//        String line = null;
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            Set<String> genes = new HashSet<String>();
//            for (int i = 2; i < tokens.length; i++)
//                genes.add(tokens[i]);
//            pathwayToGenes.put(tokens[0], genes);
//        }
//        fu.close();
//        return pathwayToGenes;
//    }
//    
//    /**
//     * This method is used to generate a map from topics to gene names contained by
//     * the topics. The species for the pathways can be set in the code.
//     * @throws Exception
//     */
//    @Test
//    public void generateTopicToProteinIdsMap() throws Exception {
//        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "gk_current_ver24",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
//        analyzer.setMySQLAdaptor(dba);
//        ReactomeAnalyzerTopicHelper topicHelper = new ReactomeAnalyzerTopicHelper();
//        GKInstance mouse = dba.fetchInstance(48892L);
//        Collection pathways = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway, 
//                                                           ReactomeJavaConstants.species, 
//                                                           "=",
//                                                           mouse);
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "MousePathwayToProteins.txt";
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        StringBuilder builder = new StringBuilder();
//        for (Iterator it = pathways.iterator(); it.hasNext();) {
//            GKInstance pathway = (GKInstance) it.next();
//            // Need to avoid human pathways since species are a multiple value
//            List species = pathway.getAttributeValuesList(ReactomeJavaConstants.species);
//            if (species.size() > 1)
//                continue;
//            Set<String> proteinIds = new HashSet<String>();
//            Set<GKInstance> participants = analyzer.grepPathwayParticipants(pathway);
//            // Want to get gene names from participants
//            for (GKInstance part : participants) {
//                Set<GKInstance> refPepSeqs = topicHelper.grepRefPepSeqs(part);
//                for (GKInstance refPepSeq : refPepSeqs) {
//                    String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
//                    if (identifier != null)
//                        proteinIds.add(identifier);
//                }
//            }
//            System.out.println(pathway + ": " + proteinIds.size());
//            //if (++c > 20)
//            //    break;
//            builder.append(pathway.getDisplayName());
//            for (String id : proteinIds) {
//                builder.append("\t");
//                builder.append(id);
//            }
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
//    }
//    
//    @Test
//    public void countProteinsInTopics() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinIdToTopic052708.txt";
//        Map<String, Set<String>> protein2Topic = fu.loadSetMap(fileName);
//        Map<String, Integer> proteinInTopics = countProteinsInTopics(protein2Topic);
//        exportProteinNumberMap(fu, 
//                               FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNumberInTopics052708.txt",
//                               proteinInTopics);
//    }
//    
//    @Test
//    public void countProteinsInTopicsWithHops() throws IOException {
//        FileUtility fu = new FileUtility();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinIdToClusteredPathwayWithHop.txt";
//        fu.setOutput(fileName);
//        List<String> topics = null;
//        int total = 0;
//        int preNumber = 0;
//        // Cache all results for printing
//        List<Map<String, Integer>> results = new ArrayList<Map<String, Integer>>();
//        List<Integer> totals = new ArrayList<Integer>();
//        for (int i = 0; i < 10; i++) {
//            Map<String, Set<String>> idToTopics = generateIdToClusteredTopicsWithHop(i);
//            total = idToTopics.size();
//            totals.add(total);
//            if (total == preNumber)
//                break; // Don't need do again
//            Map<String, Integer> topicToNumber = countProteinsInTopics(idToTopics);
//            if (topics == null) {
//                topics = new ArrayList<String>(topicToNumber.keySet());
//                Collections.sort(topics);
//            }
//            results.add(topicToNumber);
//            preNumber = total;
//        }
//        int size = results.size();
//        StringBuilder builder = new StringBuilder();
//        builder.append("Pathway\t");
//        for (int i = 0; i < size; i++) {
//            builder.append("hop" + i);
//            if (i < size - 1)
//                builder.append("\t");
//        }
//        fu.printLine(builder.toString());
//        // Print out totals
//        builder.setLength(0);
//        builder.append("Total\t");
//        for (int i = 0; i < size; i++) {
//            builder.append(totals.get(i));
//            if (i < size - 1)
//                builder.append("\t");
//        }
//        fu.printLine(builder.toString());
//        for (String topic : topics) {
//            int hop = 0;
//            builder.setLength(0);
//            builder.append(topic).append("\t");
//            for (Map<String, Integer> map : results) {
//                Integer number = map.get(topic);
//                builder.append(number);
//                if (hop < results.size() - 1)
//                    builder.append("\t");
//                hop ++;
//            }
//            fu.printLine(builder.toString());
//        }
//        fu.close();
//    }
//
//    public Map<String, Integer> countProteinsInTopics(Map<String, Set<String>> protein2Topic) {
//        Map<String, Integer> topicToProteinNumber = new HashMap<String, Integer>();
//        for (Iterator<String> it = protein2Topic.keySet().iterator(); it.hasNext();) {
//            String id = it.next();
//            Set<String> topics = protein2Topic.get(id);
//            for (String t : topics) {
//                Integer c = topicToProteinNumber.get(t);
//                if (c == null)
//                    topicToProteinNumber.put(t, 1);
//                else {
//                    topicToProteinNumber.put(t, ++c);
//                }
//            }
//        }
//        return topicToProteinNumber;
//    }
//
//    private void exportProteinNumberMap(FileUtility fu,
//                                        String fileName,
//                                        final Map<String, Integer> proteinInTopics) throws IOException {
//        List<String> topics = new ArrayList<String>(proteinInTopics.keySet());
//        // Sort based on protein numbers
//        Collections.sort(topics, new Comparator<String>() {;
//            public int compare(String obj1, String obj2) {
//                Integer c1 = proteinInTopics.get(obj1);
//                Integer c2 = proteinInTopics.get(obj2);
//                return c2.compareTo(c1);
//            }
//        });
//        fu.setOutput(fileName);
//        fu.printLine("Index\tProteinNumber\tPathway");
//        int index = 1;
//        for (String t : topics) {
//            Integer c = proteinInTopics.get(t);
//            fu.printLine(index + "\t" + c + "\t" + t);
//            index ++;
//        }
//        fu.close();
//    }
//    
//    @Test
//    public void generateFIToTopicMap() throws Exception {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIToTopic.txt";
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
//        for (ReactomeAnalyzer analyzer : analyzers) {
//            String sourceLetter = ReactomeAnalyzer.getSourceLetter(analyzer);
//            Map<GKInstance, Set<String>> topicToFIs = analyzer.grepInteractionsForTopics();
//            for (Iterator<GKInstance> it = topicToFIs.keySet().iterator(); it.hasNext();) {
//                GKInstance topic = it.next();
//                Set<String> fis = topicToFIs.get(topic);
//                String topicName = topic.getDisplayName() + "(" + sourceLetter + ")\t";
//                for (String fi : fis)
//                    fu.printLine(topicName + fi);
//            }
//        }
//        fu.close();
//    }
//    
//    public Map<String, Double> loadTopicToFIRatio() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "TopicToFINumber.txt";
//        Map<String, Double> topicToFIRatio = new HashMap<String, Double>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        while ((line = fu.readLine()) != null) {
//            String[] tokens = line.split("\t");
//            topicToFIRatio.put(tokens[0],
//                               new Double(tokens[2]));
//        }
//        fu.close();
//        return topicToFIRatio;
//    }
//    
//    @Test
//    public void generateTopicToFINumber() throws IOException {
//        Map<String, Set<String>> fiToTopic = loadFIToTopics();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "TopicToFINumber.txt";
//        final Map<String, Integer> topicToNumber = new HashMap<String, Integer>();
//        for (Iterator<String> it = fiToTopic.keySet().iterator(); it.hasNext();) {
//            String fi = it.next();
//            Set<String> topics = fiToTopic.get(fi);
//            for (String topic : topics) {
//                Integer c = topicToNumber.get(topic);
//                if (c == null)
//                    topicToNumber.put(topic, 1);
//                else
//                    topicToNumber.put(topic, ++c);
//            }
//        }
//        // Want to output
//        // Sort first
//        List<String> topics = new ArrayList<String>(topicToNumber.keySet());
//        Collections.sort(topics, new Comparator<String>() {
//            public int compare(String topic1, String topic2) {
//                Integer c1 = topicToNumber.get(topic1);
//                Integer c2 = topicToNumber.get(topic2);
//                return c2 - c1;
//            }
//        });
//        FileUtility fu = new FileUtility();
//        fu.setOutput(fileName);
//        for (String topic : topics) {
//            Integer number = topicToNumber.get(topic);
//            double ratio = (double) number / FIConfiguration.getConfiguration().get("TOTAL_FI_INTERACTION;
//            fu.printLine(topic + "\t" + number + "\t" + ratio);
//        }
//        fu.close();
//    }
//    
//    public Map<String, Set<String>> loadFIToTopics() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIToTopic.txt";
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        Map<String, Set<String>> fiToTopics = new HashMap<String, Set<String>>();
//        int index = 0;
//        String line = null;
//        String topic;
//        String fi;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            // Translation(R)  P55884 P62263
//            topic = line.substring(0, index);
//            fi = line.substring(index + 1);
//            Set<String> topics = fiToTopics.get(fi);
//            if (topics == null) {
//                topics = new HashSet<String>();
//                fiToTopics.put(fi, topics);
//            }
//            topics.add(topic);
//        }
//        fu.close();
//        return fiToTopics;
//    }
//    
//    public Map<String, Set<String>> getTopicToIdMap() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC;
//        return getTopicToIdMap(fileName);
//    }
//
//    private Map<String, Set<String>> getTopicToIdMap(String fileName) throws IOException {
//        Map<String, Set<String>> topicToIds = new HashMap<String, Set<String>>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        int index = 0;
//        String id, topic;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            id = line.substring(0, index);
//            topic = line.substring(index + 1);
//            Set<String> set = topicToIds.get(topic);
//            if (set == null) {
//                set = new HashSet<String>();
//                topicToIds.put(topic, set);
//            }
//            set.add(id);
//        }
//        return topicToIds;
//    }
//    
//    public Map<String, Set<String>> getTopicToNamesMap() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ProteinNameToTopics051109.txt";
//        return getTopicToIdMap(fileName);
//    }
//    
//    /**
//     * This method is used to generate a GMT file for all pathways in database reactome_plus_i.
//     * @throws IOException
//     */
//    @Test
//    public void generateGSEAGmtFile() throws IOException {
//        Map<String, Set<String>> topicToNames = getTopicToNamesMap();
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "Pathways_Reactome_28_plus_i.gmt";
//        fu.setOutput(fileName);
//        StringBuilder builder = new StringBuilder();
//        for (String pathway : topicToNames.keySet()) {
//            Set<String> names = topicToNames.get(pathway);
//            builder.append(pathway).append("\tna");
//            for (String name : names)
//                builder.append("\t").append(name);
//            fu.printLine(builder.toString());
//            builder.setLength(0);
//        }
//        fu.close();
//    }
//    
//    /**
//     * Protein ids used in pathways
//     * @return
//     * @throws IOException
//     */
//    public Set<String> getTopicIds() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC;
//        Set<String> ids = new HashSet<String>();
//        FileUtility fu = new FileUtility();
//        fu.setInput(fileName);
//        String line = null;
//        int index = 0;
//        String id;
//        while ((line = fu.readLine()) != null) {
//            index = line.indexOf("\t");
//            id = line.substring(0, index);
//            ids.add(id);
//        }
//        return ids;
//    }
//    
//    public Map<String, Set<String>> getTopicToFIMap() throws IOException {
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIToTopic.txt";
//        FileUtility fu = new FileUtility();
//        return fu.loadSetMap(fileName);
//    }
//    
//    public Map<String, Set<String>> getTopicToFIInNameMap() throws Exception {
//        Map<String, Set<String>> topicToFIs = getTopicToFIMap();
//        // Need to convert to gene names
//        HibernateFIReader fiReader = new HibernateFIReader();
//        Map<String, String> accToName = fiReader.generateAccessionToProteinNames();
//        Map<String, Set<String>> topicToFIInNames = new HashMap<String, Set<String>>();
//        for (String topic : topicToFIs.keySet()) {
//            Set<String> fis = topicToFIs.get(topic);
//            Set<String> fiInNames = new HashSet<String>();
//            for (String fi : fis) {
//                int index = fi.indexOf(" ");
//                String acc1 = fi.substring(0, index);
//                String acc2 = fi.substring(index + 1);
//                String name1 = accToName.get(acc1);
//                String name2 = accToName.get(acc2);
//                if (name1 == null || name2 == null)
//                    continue;
//                int compare = name1.compareTo(name2);
//                if (compare < 0)
//                    fiInNames.add(name1 + "\t" + name2);
//                else if (compare > 0)
//                    fiInNames.add(name2 + "\t" + name1);
//            }
//            topicToFIInNames.put(topic, fiInNames);
//        }
//        return topicToFIInNames;
//    }
//    
//    public Map<String, Set<String>> generateIdToTopicWithHop(int totalHop) throws IOException {
//        Map<String, Set<String>> idToTopics = new FileUtility().loadSetMap(FIConfiguration.getConfiguration().get("PROTEIN_ID_TO_TOPIC);
//        return generateIdToTopics(totalHop, idToTopics);
//    }
//    
//    public Map<String, Set<String>> generateIdToClusteredTopicsWithHop(int totalHop) throws IOException {
//        List<PathwayCluster> clusters = loadPathwayClusters();
//        Map<String, Double> topicToRatio = new HashMap<String, Double>();
//        Map<String, Set<String>> idToClusters = new HashMap<String, Set<String>>();
//        for (PathwayCluster cluster : clusters) {
//            String label = cluster.getLabel();
//            Set<String> proteins = cluster.getProteins();
//            topicToRatio.put(label, proteins.size() / (double) FIConfiguration.getConfiguration().get("TOTAL_PROTEIN_IDS);
//            for (String id : proteins) {
//                Set<String> pathways = idToClusters.get(id);
//                if (pathways == null) {
//                    pathways = new HashSet<String>();
//                    idToClusters.put(id, pathways);
//                }
//                pathways.add(label);
//            }
//        }
//        return generateIdToTopics(totalHop, idToClusters);
//    }
//
//    public Map<String, Set<String>> generateIdToTopics(int totalHop, 
//                                                       Map<String, Set<String>> idToTopics) throws IOException {
//        FIFileAnalyzer fiAnalyzer = new FIFileAnalyzer();
//        Set<String> totalIds = fiAnalyzer.loadInteractionIds();
//        totalIds.removeAll(idToTopics.keySet());
//        Map<String, Set<String>> idToPartners = fiAnalyzer.loadIdToPartners();
//        if (totalHop == -1)
//            totalHop = Integer.MAX_VALUE; // Give it a real big number
//        int hop = 0;
//        while (hop < totalHop && totalIds.size() > 0) {
//            Map<String, Set<String>> hopIdToTopics = annotateIdsByHop(totalIds, 
//                                                                      idToTopics, 
//                                                                      idToPartners);
//            if (hopIdToTopics.size() == 0)
//                break; // Cannnot annotate more via hopping
//            totalIds.removeAll(hopIdToTopics.keySet());
//            idToTopics.putAll(hopIdToTopics);
//            hop ++;
//        }
//        return idToTopics;
//    }
//    
//    public Map<String, Set<String>> annotateIdsByHop(Set<String> idsToBeAnnotated,
//                                                      Map<String, Set<String>> idToTopics,
//                                                      Map<String, Set<String>> idToPartners) {
//        Map<String, Set<String>> hopIdToTopics = new HashMap<String, Set<String>>();
//        for (String id : idsToBeAnnotated) {
//            Set<String> partners = idToPartners.get(id);
//            for (String partner : partners) {
//                Set<String> topics = idToTopics.get(partner);
//                if (topics == null || topics.size() == 0)
//                    continue;
//                Set<String> newTopics = hopIdToTopics.get(id);
//                if (newTopics == null) {
//                    newTopics = new HashSet<String>();
//                    hopIdToTopics.put(id, newTopics);
//                }
//                newTopics.addAll(topics);
//            }
//        }
//        return hopIdToTopics;
//    }
}
