/*
 * Created on Jun 30, 2006
 *
 */
package org.reactome.weka;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.Value;

public class WEKADataAnalyzer {
    // Help to generate ARFF file
    private ARFFGenerator arffGenerator;
    
    public WEKADataAnalyzer() {
        arffGenerator = new ARFFGenerator();
    }
    
    @Test
    public void generateLeaveOneOutDataSet() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        List<ReactomeAnalyzer> otherAnalyzers = new ArrayList<ReactomeAnalyzer>();
        // This loops is used to get all non-reactome datasources
        for (ReactomeAnalyzer analyzer : analyzers) {
            GKInstance dataSource = analyzer.getDataSource();
            if (dataSource != null)
                otherAnalyzers.add(analyzer);
        }
        String outputFileName;
        Set<String> interactions = new HashSet<String>();
        Set<String> negativePairs = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : otherAnalyzers) {
            String sourceLetter = ReactomeAnalyzer.getSourceLetter(analyzer);
            // The following code is used to generate the training dataset
            outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "Leave_" + sourceLetter + "_Out121407.arff";
            interactions.clear();
            negativePairs.clear();
            for (ReactomeAnalyzer dataAnalyzer : analyzers) {
                if (dataAnalyzer == analyzer)
                    continue;
                Set<String> fetchedSet = dataAnalyzer.extractInteractionSet();
                interactions.addAll(fetchedSet);
                Set<String> fetchedPairs = dataAnalyzer.generateUniProtPairsFromTopics();
                negativePairs.addAll(fetchedPairs);
            }
            negativePairs.removeAll(interactions);
            System.out.println("Leaving " + sourceLetter + " out...");
            System.out.println("Positive pairs: " + interactions.size() + ", negative pairs: " + negativePairs.size());
            filterToHumanIds(interactions);
            filterToHumanIds(negativePairs);
            System.out.println("After filtering...");
            System.out.println("Positive pairs: " + interactions.size() + ", negative pairs: " + negativePairs.size());
            generateDataSetForWEKA(interactions, 
                                             negativePairs, 
                                             outputFileName);
            // The following code is used to generate the test dataset
            String testOutputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + sourceLetter + "_Test_121407.arff";
            Set<String> testInteractions = analyzer.extractInteractionSet();
            Set<String> testNegativePairs = analyzer.generateUniProtPairsFromTopics();
            testNegativePairs.removeAll(testInteractions);
            System.out.println("Information for test dataset...");
            System.out.println("Before removing training pairs: positive " + testInteractions.size() + ", negative " + testNegativePairs.size());
            testInteractions.removeAll(interactions);
            testNegativePairs.removeAll(negativePairs);
            System.out.println("After removing training pairs: positive " + testInteractions.size() + ", negative " + testNegativePairs.size());
            filterToHumanIds(testInteractions);
            filterToHumanIds(testNegativePairs);
            generateDataSetForWEKA(testInteractions, 
                                             testNegativePairs, 
                                             testOutputFileName);
        }
    }
    
    /**
     * Use a pre-dumped FI file to create arff file in order to increase the performance.
     * @throws Exception
     */
    @Test
    public void generateDatasetFromFile() throws Exception {
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeFIs040309_Neg_10_No_Filter.arff";
        Set<String> interactions = analyzer.loadFIsFromFile();
        // Use random pairs 
        // As estimated, positive pairs is about 1%.
        Set<String> negativePairs = generateRandomPairsFromAllProteins(interactions.size() * 10);
        negativePairs.removeAll(interactions);
        System.out.println("Total positive pairs: " + interactions.size());
        System.out.println("Total negaive pairs: " + negativePairs.size());
        // Check ratio
        double ratio = (double) interactions.size() / negativePairs.size();
        System.out.println("Positive / Negative: " + ratio);
        // No need to do filter
        arffGenerator.setIsForVersion3(true);
        generateDataSetForWEKA(interactions, 
                               negativePairs, 
                               outputFileName);       
    }
    
    @Test
    public void generateDataSetFromCombinedDatabase() throws Exception {
        //ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        //String outputFileName = "results/ReactomeData091806.arff";
        //ReactomeAnalyzer analyzer = new PantherAnalyzer();
        //String outputFileName = "results/PantherDataNoReactome091806_1.arff";
        //ReactomeAnalyzer analyzer = new CPathAnalyzer();
        //String outputFileName = "results/CellMapDataNoReactome091806_1.arff";
        //String outputFileName = "results/BIND090806.arff";
        //ReactomeAnalyzer analyzer = new HPRDAnalyzer();
        //String outputFileName = "results/HPRD2H091106.arff";
        //ReactomeAnalyzer analyzer = new INOHAnalyzer();
        //String outputFileName = "results/INOHDataNoReactome120206.arff";
        // The following code is used to aggreate all interactions from three pathway
        // database together: Reactome, Panther, CellMap. INOH is excluded since its 
        // low precision rate
        // Regulators were not considered during extracting pathway participants
        // for file SixDBInteractions020507.arff, and all files before Feb 9, 2007.
        //String outputFileName = "results/v2/SixDBInteractions020507.arff";
        // Regulators are considers in this file.
        //String outputFileName = "results/v2/SixDBInteractions020907.arff";
        // A mistake is found in ReactomeAnalyzer to grep protein pair from pathways.
        // This mistake affects only pairs from the Reactome database.
        //String outputFileName = "results/v2/SixDBInteractions041807.arff";
        // This arff file uses pairs randomly picked from all proteins
        // as a negative data set. Of couse, interactions are removed from
        // this negative data sets.
        //String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "SixDBInteractionsWithRandomPairs.arff";
        //String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "SixDBInteractions082807.arff";
        // This arff uses negative data sets generated from different cell compartments
        // in file: NoInteractions090506.txt, which was based on GO annotation.
        //String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "SixDBInteractions121307_2.arff";
        // An arff file without considering complex to see if TP can be increased
        //String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions020708_DiffCompart.arff";
        // Just a quick test
        //String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions021009_1.arff";
        String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeInteractions030609.arff";
        Set<String> interactions = new HashSet<String>();
        Set<String> negativePairs = new HashSet<String>();
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        for (ReactomeAnalyzer analyzer : analyzers) {
            //analyzer.setExcludeComplex(true);
            Set<String> fetchedSet = analyzer.extractInteractionSet();
            interactions.addAll(fetchedSet);
            Set<String> fetchedPairs = analyzer.generateUniProtPairsFromTopics();
            negativePairs.addAll(fetchedPairs);
        }
        // This is only for random pairs
        //pairsFromTopics.addAll(generateRandomPairsFromAllProteins());
        // These statements are used to get negative pairs from proteins not
        // in the same compartments
        //negativePairs.addAll(generateNoInteractionPairsFromCompartment());
        negativePairs.removeAll(interactions);
        // Try to keep positive and negative in the same amounts
        //negativePairs = randomPickPairs(negativePairs, interactions.size());
        //removeDataFromReactome(interactions, pairsFromTopics);
        // Try to get the similar ratio no-interaction pairs from cellular compartments
        //GODataAnalyzer goAnalyzer = new GODataAnalyzer();
        //negativePairs = goAnalyzer.generateProteinPairsInNucleusAndExRegion(negativePairs.size());
        //negativePairs = generateRandomPairsFromAllProteins(negativePairs.size());
//        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
//        negativePairs = generateRandomPairsFromProteinIDs(interactions.size() * 10, ids);
        System.out.println("Total negative pairs: " + negativePairs.size());
 //       negativePairs.removeAll(interactions);
        System.out.println("Remove interactions: " + negativePairs.size());
        System.out.println("Positive pairs: " + interactions.size() + ", negative pairs: " + negativePairs.size());
        filterToHumanIds(interactions);
        filterToHumanIds(negativePairs);
       //negativePairs = randomPickPairs(negativePairs, interactions.size() * 24);
        System.out.println("After filtering...");
        System.out.println("Positive pairs: " + interactions.size() + ", negative pairs: " + negativePairs.size());
        generateDataSetForWEKA(interactions, 
                               negativePairs, 
                               outputFileName);       
    }
    
    protected Set<String> randomPickPairs(Set<String> pairs,
                                        int targetSize) {
        System.out.println("Target size: " + targetSize);
        Set<String> rtn = new HashSet<String>();
        List<String> list = new ArrayList<String>(pairs);
        // Cannot do anything better
        if (targetSize >= pairs.size())
            return pairs;
        int index = 0;
        int size = 0;
        while (rtn.size() < targetSize) {
            size = list.size();
            index = (int)(Math.random() * size);
            rtn.add(list.remove(index));
        }
        System.out.println("Return Size: " + rtn.size());
        return rtn;
    }
    
    /**
     * This method will be replaced by ProteinIdFilter.cleanUpVsUniProt() soon.
     * @param pairs
     * @throws IOException
     */
    @Deprecated
    private void filterToHumanIds(Set<String> pairs) throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtIds = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = uniProtIds.keySet();
        for (Iterator<String> it = pairs.iterator(); it.hasNext();) {
            String pair = it.next();
            int index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            if (!uniSet.contains(id1) ||
                !uniSet.contains(id2)) {
                it.remove();
                //System.out.println(id1 + " or " + id2);
            }
        }
    }
    
    private Set<String> generateRandomPairsFromAllProteins(int size) throws IOException {
        Set<String> swissProtIds = new UniProtAnalyzer().loadSwissProtIds();
        return generateRandomPairsFromProteinIDs(size, swissProtIds);
    }

    private Set<String> generateRandomPairsFromProteinIDs(int size, Set<String> swissProtIds) {
        List<String> idList = new ArrayList<String>(swissProtIds);
        Set<String> pairs = new HashSet<String>();
        // Based on the total lines in the file SixDBInteractions041807.arff
        int total = size;
        int comp = 0;
        int index1, index2;
        String id1, id2;
        int totalSize = idList.size();
        while (pairs.size() < total) {
            index1 = (int) (Math.random() * totalSize);
            index2 = (int) (Math.random() * totalSize);
            if (index1 == index2)
                continue;
            id1 = idList.get(index1);
            id2 = idList.get(index2);
            comp = id1.compareTo(id2);
            if (comp < 0)
                pairs.add(id1 + " " + id2);
            else if (comp > 0)
                pairs.add(id2 + " " + id1);
        }
        return pairs;
    }
    
    public void generateDataSetFromNoInteractionDataSetFromLocalization() throws Exception {
        //String outputFileName = "results/NoLocalizationDataNoReactome091806_1.arff";
        String outputFileName = "/Users/wgm/Documents/gkteam/Carl/RafY2H.arff";
        //String outputFileName = "results/ThreePPIsData091806.arff";
        //String outputFileName = "results/HumanPPIsData091806.arff";
        Set<String> posSet = new HashSet<String>();
        //String fileName = "results/interaction/NoInteractions090506.txt";
        //String fileName = "results/interaction/ThreePPIs091406.txt";
        //String fileName = "results/interaction/HumanInteractions091106.txt";
        String fileName = "/Users/wgm/Documents/gkteam/Carl/RafY2H.txt";
        Set<String> negSet = new FileUtility().loadInteractions(fileName);
        System.out.println("Negative Interactions: " + negSet.size());
        System.out.println("Positive Interactions: " + posSet.size());
        //removeDataFromReactome(posSet, negSet);
        generateDataSetForWEKA(posSet, negSet, outputFileName);  
    }
    
    private Set<String> generateNoInteractionPairsFromCompartment() throws Exception {
        String fileName = "results/interaction/NoInteractions090506.txt";
        FileUtility fu = new FileUtility();
        return fu.loadInteractions(fileName);
    }
    
    public void generateDataSetForWEKA(Set<String> positiveInteractions,
                                       Set<String> negativeInteractions,
                                       String outFileName) throws Exception {
        Map<String, Value> values = new HashMap<String, Value>();
        // Initialize values first
        Value value;
        for (String pair : positiveInteractions) {
            value = new Value();
            value.functionalInteraction = Boolean.TRUE;
            values.put(pair, value);
        }
        for (String pair : negativeInteractions) {
            value = new Value();
            value.functionalInteraction = Boolean.FALSE;
            values.put(pair, value);
        }
        arffGenerator.generateDataSet(values);
        arffGenerator.exportDataForWEKA(outFileName, values);
    }
    
    public void generateDataSet(Map<String, Value> values) throws Exception {
        arffGenerator.generateDataSet(values);
    }
    
    private void removeDataFromReactome(Set<String> posSet,
                                        Set<String> negSet) throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        if (posSet.size() > 0) {
            Set<String> reactomeInteractions = reactomeAnalyzer.extractInteractionSet();
            System.out.println("Size of PosSet: " + posSet.size());
            posSet.removeAll(reactomeInteractions);
            System.out.println("Size of PosSet after removing Reactome: " + posSet.size());
        }
        if (negSet.size() > 0) {
            Set<String> reactomePairs = reactomeAnalyzer.generateUniProtPairsFromTopics();
            System.out.println("Size of NegSet: " + negSet.size());
            negSet.removeAll(reactomePairs);
            System.out.println("Size of NegSet after removing Reactome: " + negSet.size());
        }
    }
    
    private void divideIds(List<String> allIds,
                           List<String> ids1, 
                           List<String> ids2, 
                           double ratio) {
        int totalSize = allIds.size();
        int list1Size = (int) (totalSize * ratio);
        // Using no-replacement sampling
        Set<Integer> list1Index = new HashSet<Integer>();
        while (true) {
            list1Index.add((int)(Math.random() * totalSize));
            if (list1Index.size() > list1Size)
                break;
        }
        for (Integer index : list1Index) 
            ids1.add(allIds.get(index));
        ids2.addAll(allIds);
        ids2.removeAll(ids1);
    }
    
}
