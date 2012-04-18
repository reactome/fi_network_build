/*
 * Created on Apr 25, 2007
 *
 */
package org.reactome.fi;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.tred.TREDAnalyzer;

/**
 * This class is used to handle the interaction file generated from the hibernate adaptor.
 * @author guanming
 *
 */
public class FIFileAnalyzer {
    private final String pathwayFIFile = FIConfiguration.getConfiguration().get("RESULT_DIR") + "PathwayFIs040909.txt";
    private final String predictedFIFile = FIConfiguration.getConfiguration().get("PREDICTED_FI_FILE");
    
    private FileUtility fu;
    
    public FIFileAnalyzer() {
        fu = new FileUtility();
    }
    
//    /**
//     * This method is used to remove the ZNF clique.
//     * @throws IOException
//     */
//    @Test
//    public void generateFIFileWithZNFRemoved() throws IOException {
//        Set<String> fis = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
//        Collection<String> znfClique = new GraphAnalyzer().searchZNFClique(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
////        System.out.println("Genes in the ZNF clique: " + znfClique.size());
////        System.out.println(znfClique);
//        Set<String> znfFIs = InteractionUtilities.getFIs(znfClique, fis);
//        fis.removeAll(znfFIs);
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIsInGene_No_ZNF_042810.txt";
//        fu.saveInteractions(fis, fileName);
//    }
    
    /**
     * This method is used to check ZNF genes in the FI files.
     * @throws IOException
     */
    @Test
    public void checkZNFInFIs() throws IOException {
        Set<String> fis = fu.loadInteractions(FIConfiguration.getConfiguration().get("GENE_FI_FILE_NAME"));
        Set<String> genes = InteractionUtilities.grepIDsFromInteractions(fis);
        Set<String> znfs = new HashSet<String>();
        for (String gene : genes) {
            if (gene.startsWith("ZNF"))
                znfs.add(gene);
        }
        System.out.println("Total Genes: " + genes.size());
        System.out.println("Total ZNFs: " + znfs.size() + " (" + (double)znfs.size() / genes.size() + ")");
        Set<String> touchedZNFFIs = new HashSet<String>();
        Set<String> allZNFFIs = new HashSet<String>();
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String gene1 = fi.substring(0, index);
            String gene2 = fi.substring(index + 1);
            if (gene1.startsWith("ZNF") || gene2.startsWith("ZNF"))
                touchedZNFFIs.add(fi);
            if (gene1.startsWith("ZNF") && gene2.startsWith("ZNF"))
                allZNFFIs.add(fi);
        }
        System.out.println("Total FIs: " + fis.size());
        System.out.println("Total FIs having at least one ZNF: " + touchedZNFFIs.size() + " (" + (double)touchedZNFFIs.size() / fis.size() + ")");
        System.out.println("Total FIs having both ZNFs: " + allZNFFIs.size() + " (" + (double)allZNFFIs.size() / fis.size() + ")");
    }
    
    /**
     * Use this method to save a set of FIs into an order list.
     * @param fis
     * @param fileName
     * @throws IOException
     */
    public void saveFIInOrder(Set<String> fis,
                              String fileName) throws IOException {
        List<String> list = new ArrayList<String>(fis);
        Collections.sort(list);
        fu.setOutput(fileName);
        for (String fi : list)
            fu.printLine(fi);
        fu.close();
    }
    
    /**
     * Load all FIs in UniProt accession numbers.
     * @return
     * @throws IOException
     */
    public Set<String> loadFIs() throws IOException {
        Set<String> fis = new HashSet<String>();
        fis.addAll(loadPathwayFIs());
        fis.addAll(loadTFTargetInteractions());
        fis.addAll(loadPredictedFIs());
        return fis;
        //return fu.loadInteractions(FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME);
    }
    
    private Set<String> loadArtificalFIs() throws IOException {
        Set<String> fis = new HashSet<String>();
        fis.add("A B");
        fis.add("A C");
        fis.add("A D");
        fis.add("A G");
        fis.add("B D");
        fis.add("C D");
        fis.add("C F");
        fis.add("D E");
        fis.add("D G");
        fis.add("D H");
        fis.add("E F");
        fis.add("E H");
        fis.add("G H");
        return fis;
    }
    
    /**
     * This class is used to check the SwissProt coverage in the functional
     * interaction file.
     * @throws IOException
     */
    @Test
    public void checkSwissProtCoverage() throws IOException {
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        // Need to use this map in case some UniProt accession numbers
        // used in pathways are not the first one!
        Map<String, String> map = uniProtAnalyzer.loadSwissProtIDsMap();
        Set<String> mapIds = new HashSet<String>(map.values());
        // There are three files that need to be check
        String[] fileNames = new String[] {
//                "FIInteractions73_021108_Pathway.txt",
//                "FIInteractions73_021108_PPI.txt",
//                "FIInteractions73_021108.txt"
                "FI73_041408.txt"
        };
        for (String name : fileNames) {
            System.out.println("File: " + name);
            Set<String> interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + name);
            System.out.println("Total interactions: " + interactions.size());
            Set<String> totalIds = InteractionUtilities.grepIDsFromInteractions(interactions);
            System.out.println("Total IDs: " + totalIds.size());
            totalIds = removeSpliceIsoform(totalIds);
            System.out.println("Remove isoforms: " + totalIds.size());
            // 25205 is the total identifiers in HPRD and used as the total
            // gene numbers
            System.out.println("Total coverage: " + totalIds.size() / 25205.0);
            // Check in other way
            Set<String> tmpIds = new HashSet<String>();
            for (String id : totalIds) {
                String tmp = map.get(id);
                if (tmp != null)
                    tmpIds.add(tmp);
            }
            System.out.println("Swiss Prot Ids in FIs: " + tmpIds.size());
            System.out.println("Swiss Coverage: " + tmpIds.size() + "/" + 
                               mapIds.size() + "=" + (double)tmpIds.size() / mapIds.size());
            System.out.println();
        }
    }
    
    private Set<String> removeSpliceIsoform(Set<String> ids) {
        Set<String> rtn = new HashSet<String>();
        int index = 0;
        for (String id : ids) {
            index = id.indexOf("-");
            if (index > 0)
                rtn.add(id.substring(0, index));
            else
                rtn.add(id);
        }
        return rtn;
    }
    
    public Set<String> loadInteractionIds() throws IOException {
        Set<String> fis = loadFIs();
        return InteractionUtilities.grepIDsFromInteractions(fis);
    }
    
    @Test
    public void checkHumanIDsInFIs() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtIds = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = uniProtIds.keySet();
        Set<String> fis = loadFIs();
        int index = 0;
        int totalFIs = 0;
        int removeFIs = 0;
        Set<String> idsFromHPRDOrNCBI = new HashSet<String>();
        Set<String> removedIds = new HashSet<String>();
        // Need to get alternative form out
        for (String pair : fis) {
            totalFIs ++;
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            if (id1.matches("^[0-9]+")) {
                idsFromHPRDOrNCBI.add(id1);
            }
            if (id2.matches("^[0-9]+"))
                idsFromHPRDOrNCBI.add(id2);
            index = id1.indexOf("-");
            if (index > 0)
                id1 = id1.substring(0, index);
            index = id2.indexOf("-");
            if (index > 0)
                id2 = id2.substring(0, index);
            // Check if id1 or id2 are numbers only: for NCBI or HPRD
            if (!id1.matches("^[0-9]+") && !uniSet.contains(id1)) {
                removedIds.add(id1);
                removeFIs ++;
                continue;
            }
            if (!id2.matches("^[0-9]+") && !uniSet.contains(id2)) {
                removedIds.add(id2);
                removeFIs ++;
                continue;
            }
            // Otherwise, have to make sure they are from human UniProt IDs
        }
        System.out.println("Total FIs: " + totalFIs);
        System.out.println("Remove FIs: " + removeFIs);
        System.out.println("Total IDs from NCBI or HPRD: " + idsFromHPRDOrNCBI.size() + ": " + idsFromHPRDOrNCBI);
        System.out.println("Removed Ids: " + removedIds.size());
        for (String id : removedIds)
            System.out.println(id);
    }
    
//    public Map<String, Set<String>> loadIdToPartners() throws IOException {
//        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
//        Set<String> fis = loadFIs();
//        return new BreadthFirstSearch().generateIdToPartnersMap(fis);
//    }
    
//    
//    @Test
//    public void analyzeNonPathwayIds() throws IOException {
//        Set<String> totalIds = loadInteractionIds();
//        TopicAnalyzer topicAnalyzer = new TopicAnalyzer();
//        Set<String> topicIds = topicAnalyzer.getTopicIds();
//        System.out.println("Total Ids: " + totalIds.size());
//        System.out.println("Pathway Ids: " + topicIds.size());
//        totalIds.removeAll(topicIds);
//        System.out.println("Non Pathway Ids: " + totalIds.size());
//        Map<String, Set<String>> idToPartners = loadIdToPartners();
//        int hop = 1;
//        while (totalIds.size() > 0) {
//            System.out.println("Hop: " + hop);
//            Set<String> hopAnnotated = checkNonPathwayIds(totalIds, 
//                                                          topicIds, 
//                                                          idToPartners);
//            if (hopAnnotated.size() == 0) {
//                // Cannot be annotated by hop anymore
//                System.out.println("Cannot annotated by hop!");
//                break;
//            }
//            System.out.println("Annotated by hopping: " + hopAnnotated.size());
//            totalIds.removeAll(hopAnnotated);
//            System.out.println("Not Annotated: " + totalIds.size());
//            topicIds.addAll(hopAnnotated);
//            hop ++;
//        }
//    }
    
    private Set<String> checkNonPathwayIds(Set<String> nonAnnotatedIds, 
                                           Set<String> topicIds, 
                                           Map<String, Set<String>> idToPartners) {
        Set<String> annotatedAfterHop = new HashSet<String>();
        for (Iterator<String> it = nonAnnotatedIds.iterator(); it.hasNext();) {
            String id = it.next();
            Set<String> partners = idToPartners.get(id);
            // Check if any partners is annotated
            for (String partner : partners) {
                if (topicIds.contains(partner)) {
                    annotatedAfterHop.add(id);
                    break;
                }
            }
        }
        return annotatedAfterHop;
    }
    
    /**
     * Generate a FI file in upper case
     * @throws IOException
     */
    @Test
    public void upperCaseForFIFile() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGeneUpperCase_111208.txt";
        String inFile = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGene_102908.txt";
        Set<String> fis = fu.loadInteractions(inFile);
        Set<String> newFis = new HashSet<String>();
        for (String fi : fis) {
            int index = fi.indexOf("\t");
            String name1 = fi.substring(0, index).toUpperCase();
            String name2 = fi.substring(index + 1).toUpperCase();
            int compare = name1.compareTo(name2);
            if (compare < 0)
                newFis.add(name1 + "\t" + name2);
            else if (compare > 0)
                newFis.add(name2 + "\t" + name1);
        }
        saveFIInOrder(newFis, fileName);
    }
    
    /**
     * This method is used to check the name case (upper or lower) usage for protein or
     * gene names in a FI file.
     * @throws IOException
     */
    @Test
    public void analyzeCaseInFIInNames() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FI73InGene_102908.txt";
        Set<String> fis = fu.loadInteractions(fileName);
        Set<String> names = InteractionUtilities.grepIDsFromInteractions(fis);
        Map<String, Set<String>> upperToNames = new HashMap<String, Set<String>>();
        for (String name : names) {
            String upper = name.toUpperCase();
            Set<String> set = upperToNames.get(upper);
            if (set == null) {
                set = new HashSet<String>();
                upperToNames.put(upper, set);
            }
            set.add(name);
        }
        // Check more than one case
        System.out.println("Names having more than one cases:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() > 1) {
                System.out.println(upper + ": " + set);
            }
        }
        System.out.println("\n\n");
        System.out.println("Names using lower cases:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() == 1) {
                String name = set.iterator().next();
                if (!name.equals(upper))
                    System.out.println(upper + ": " + name);
            }
        }
        Map<String, Set<String>> proteinToPartners = InteractionUtilities.generateProteinToPartners(fis);
        // Check two cases have the same interactions
        System.out.println("\n\nCheck partners:");
        for (String upper : upperToNames.keySet()) {
            Set<String> set = upperToNames.get(upper);
            if (set.size() > 1) {
                System.out.println(upper);
                Iterator<String> it = set.iterator();
                String name1 = it.next();
                String name2 = it.next();
                Set<String> partners1 = proteinToPartners.get(name1);
                Set<String> partners2 = proteinToPartners.get(name2);
                if (partners1.equals(partners2))
                    System.out.println(name1 + ", " + name2 + " have the same partners!");
                else {
                    Set<String> shared = new HashSet<String>(partners1);
                    shared.retainAll(partners2);
                    partners1.removeAll(shared);
                    partners2.removeAll(shared);
                    System.out.println(name1 + ": " + partners1);
                    System.out.println(name2 + ": " + partners2);
                }
            }
        }
    }
    
    /**
     * This method is used to create files for pathway FIs.
     * @throws Exception
     */
    @Test
    public void dumpPathwayFIs() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        ProteinIdFilters filters = new ProteinIdFilters();
        for (ReactomeAnalyzer analyzer : analyzers) {
            Set<String> fis = analyzer.extractInteractionSet();
            Set<String> normalized = filters.normalizeProteinPairs(fis);
            GKInstance dataSource = analyzer.getDataSource();
            String sourceName = null;
            if (dataSource == null)
                sourceName = "Reactome";
            else
                sourceName = dataSource.getDisplayName();
            System.out.println("Done data source: " + sourceName);
            String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "/FIs_" + sourceName + ".txt";
            fu.saveInteractions(normalized, fileName);
            System.out.println();
        }
    }
    
    /**
     * Merge all pre-dumped pathway FI files into one: PathwayFIs????.txt.
     * @throws Exception
     */
    @Test
    public void generateOnePathwayFIFile() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        Set<String> pathwayFIs = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : analyzers) {
            GKInstance dataSource = analyzer.getDataSource();
            String sourceName = null;
            if (dataSource == null)
                sourceName = "Reactome";
            else
                sourceName = dataSource.getDisplayName();
            String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_" + sourceName + ".txt";
            Set<String> fis = fu.loadInteractions(fileName);
            pathwayFIs.addAll(fis);
        }
        fu.saveInteractions(pathwayFIs, pathwayFIFile);
    }
    
    @Test
    public void checkTotalPathwayFIs() throws Exception {
        Set<String> allFIs = loadPathwayFIsFromFiles();
        System.out.println("Total pathway FIs: " + allFIs.size());
        Set<String> tredFIs = loadTFTargetInteractions();
        System.out.println("TRED FIs: " + tredFIs.size());
        allFIs.addAll(tredFIs);
        System.out.println("After merging: " + allFIs.size());
        ProteinAndInteractionCount counter = new ProteinAndInteractionCount();
        counter.countVsSwissProt(InteractionUtilities.grepIDsFromInteractions(allFIs));
    }
    
    /**
     * Load pathway FIs in UniProt ids.
     * @return
     * @throws IOException
     */
    @Deprecated
    public Set<String> loadPathwayFIs() throws IOException {
        return loadPathwayFIsFromFiles();
//        return fu.loadInteractions(pathwayFIFile);
    }
    
    @Test
    public void generateFIFiles() throws IOException {
        Set<String> fis = loadPathwayAndTFTargetFIs();
        // This file contains both pathways and TF/Target FIs.
        fu.saveInteractions(fis, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Pathway_043009.txt");
        Set<String> predictedFIs = loadPredictedFIs();
        fu.saveInteractions(predictedFIs, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Predicted_043009.txt");
        fis.addAll(predictedFIs);
        fu.saveInteractions(fis, FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_043009.txt");
    }

    public Set<String> loadPathwayAndTFTargetFIs() throws IOException {
        // Pathway FIs
        Set<String> fis = loadPathwayFIsFromFiles();
        Set<String> tredFIs = loadTFTargetInteractions();
        fis.addAll(tredFIs);
        return fis;
    }
    
    /**
     * Load predicted FIs in UniProt ids.
     * @return
     * @throws IOException
     */
    public Set<String> loadPredictedFIs() throws IOException {
        return fu.loadInteractions(predictedFIFile);
    }
    
    /**
     * Load TF/Target interactions.
     * @return
     * @throws IOException
     */
    public Set<String> loadTFTargetInteractions() throws IOException {
        return fu.loadInteractions(FIConfiguration.getConfiguration().get("TRED_FI_FILE"));
        //        return fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + "TREDInteractionsInUniProt.txt");
    }
    
    /**
     * This file is used to check the final FI network size.
     * @throws Exception
     */
    @Test
    public void checkFinalNetworkSize() throws Exception {
        ProteinAndInteractionCount counter = new ProteinAndInteractionCount();
        Set<String> pathwayFIs = loadPathwayFIs();
        System.out.println("Pathway FIs:");
        countFinalNetwork(counter, pathwayFIs);
        Set<String> predictedFIs = loadPredictedFIs();
        System.out.println("Predicted FIs:");
        countFinalNetwork(counter, predictedFIs);
        Set<String> tfTargetInteractions = new TREDAnalyzer().loadTFTargetInteractions();
        System.out.println("TF/Target interactions:");
        countFinalNetwork(counter, tfTargetInteractions);
        // Pathways and predicted
        predictedFIs.addAll(pathwayFIs);
        System.out.println("Pathway + Predicted FIs:");
        countFinalNetwork(counter, predictedFIs);
        // Pathways, predicted and TF/Targets
        predictedFIs.addAll(tfTargetInteractions);
        System.out.println("Pathway + Predicted + TF/Target FIs:");
        countFinalNetwork(counter, predictedFIs);
        // Count pathways and TF/Targets
        tfTargetInteractions.addAll(pathwayFIs);
        System.out.println("Pathway + TF/Target FIs:");
        countFinalNetwork(counter, tfTargetInteractions);
    }

    private void countFinalNetwork(ProteinAndInteractionCount counter,
                                   Set<String> pathwayFIs) throws IOException {
        Set<String> pathwayIds = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        System.out.println("Total FIs: " + pathwayFIs.size() + " (" + pathwayIds.size() + ")");
        counter.countVsSwissProt(pathwayIds);
        System.out.println();
    }

    /**
     * Load a set of pre-generated normalized pathway FIs from files.
     * @return
     * @throws IOException
     */
    public Set<String> loadPathwayFIsFromFiles() throws IOException {
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("REACTOME_FI_FILE"),
                FIConfiguration.getConfiguration().get("KEGG_FI_FILE"),
                FIConfiguration.getConfiguration().get("NCI_PID_FI_FILE"),
                FIConfiguration.getConfiguration().get("NCI_PID_BIOCARTA_FI_FILE"),
                FIConfiguration.getConfiguration().get("PANTHER_FI_FILE")
        };
        Set<String> allFIs = new HashSet<String>();
        for (String fileName : fileNames) {
            Set<String> fis = fu.loadInteractions(fileName);
            allFIs.addAll(fis);
        }
        return allFIs;
    }
}
