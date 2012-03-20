/*
 * Created on Sep 27, 2006
 *
 */
package org.reactome.fi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.Value;
import org.reactome.weka.WEKADataAnalyzer;
import org.reactome.weka.WEKAResultAnalyzer;

import weka.classifiers.Classifier;
import weka.core.Instance;
import weka.core.Instances;

/**
 * This class is used to get statistic data for functional interactions.
 * @author guanming
 *
 */
public class FunctionalInteractionAnalyzer extends TestCase {
    
    //private final String NBC_MODEL_FILE_NAME = "results/NoMFDisc20NaiveBayes091506.model";
    //private final String NBC_MODEL_FILE_NAME = "results/FourDBInteractions120506.model";
    private Classifier classifier;
    private FileUtility fu;
    private ProteinIdFilters idFilter;
    
    public FunctionalInteractionAnalyzer() {
    }
    
    public void setUp() throws Exception {
        fu = new FileUtility();
        classifier = (Classifier) fu.loadObject(FIConfiguration.getConfiguration().get("NBC_MODEL_NAME"));
    }
    
    public void generateFIInteractionsWithComplexAndSet() throws Exception {
        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
        Set<String> all = new HashSet<String>();
        for (ReactomeAnalyzer analyzer : analyzers) {
            Set<String> interactions = analyzer.extractInteractionSetWithComplexAndSet();
            all.addAll(interactions);
        }
        // Output
        String outputFileName = "results/v2/FIInteractionsWithComplexAndSet.txt";
        fu.saveInteractions(all, outputFileName);
    }
    
    public void generateNonReduncyInteractions() throws Exception {
        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "HumanInteractions020507.txt";
        //String fileName = "results/interaction/OrthoInteractions.txt";
        String fileName = "results/interaction/YeastInteractions.txt";
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> interactions = fu.loadInteractions(fileName);
        int index = 0;
        Set<String> newInteractions = new HashSet<String>();
        for (String i : interactions) {
            index = i.indexOf(" ");
            String acc1 = i.substring(0, index);
            String acc2 = i.substring(index + 1);
            if (acIdMap.containsKey(acc1))
                acc1 = acIdMap.get(acc1);
            if (acIdMap.containsKey(acc2))
                acc2 = acIdMap.get(acc2);
            int compare = acc1.compareTo(acc2);
            if (compare < 0)
                newInteractions.add(acc1 + " " + acc2);
            else if (compare > 0)
                newInteractions.add(acc2 + " " + acc1);
        }
        System.out.println("Old number: " + interactions.size());
        System.out.println("New number: " + newInteractions.size());
        //fu.saveInteractions(newInteractions, "HumanInteractions031808.txt");
        //fu.saveInteractions(newInteractions, FIConfiguration.getConfiguration().get("RESULT_DIR + "OrthoInteractions031908.txt");
        fu.saveInteractions(newInteractions, FIConfiguration.getConfiguration().get("RESULT_DIR") + "YeastInteractions031908.txt");
    }
    
    public Set<String> filterRedundencyInteractions(Set<String> interactions,
                                                     Map<String, String> acIdMap) throws Exception {
        System.out.println("Before renduandancy filtering: " + interactions.size());
        Set<String> filtered = new HashSet<String>();
        int index = 0;
        for (String i : interactions) {
            index = i.indexOf(" ");
            String id1 = i.substring(0, index);
            String id2 = i.substring(index + 1);
            String tmpId1 = acIdMap.get(id1);
            if (tmpId1 == null)
                tmpId1 = id1;
            String tmpId2 = acIdMap.get(id2);
            if (tmpId2 == null)
                tmpId2 = id2;
            int compare = tmpId1.compareTo(tmpId2);
            if (compare < 0) {
                filtered.add(tmpId1 + " " + tmpId2);
            }
            else if (compare > 0) {
                filtered.add(tmpId2 + " " + tmpId1);
            }
        }
        System.out.println("After redundancy filtering: " + filtered.size());
        return filtered;
    }
    
    /**
     * A helper method to filter non-human FIs.
     * @param interactions
     * @param uniSet
     * @param uniAnalyzer
     * @deprecated use @See ProteinIdFilters.filterNonHumanIs(Set<String>)
     */
    @Deprecated
    public void filterNonHumanIds(Set<String> interactions,
                                  Set<String> uniSet,
                                  UniProtAnalyzer uniAnalyzer) {
        System.out.println("Before filtering: " + interactions.size());
        for(Iterator<String> it = interactions.iterator(); it.hasNext();) {
            String interaction = it.next();
            if (containsNonHumanId(interaction, uniSet, uniAnalyzer))
                it.remove();
        }
        System.out.println("After filtering: " + interactions.size());
    }
    
    /**
     * A helper method to check if a protein pair has non-human identifiers
     * @param interaction
     * @param uniSet
     * @param uniAnalyzer
     * @return
     */
    private boolean containsNonHumanId(String interaction,
                                       Set<String> uniSet,
                                       UniProtAnalyzer uniAnalyzer) {
        int index = interaction.indexOf(" ");
        String id1 = interaction.substring(0, index);
        String id2 = interaction.substring(index + 1);
        if (!uniAnalyzer.isHumanID(uniSet, id1) ||
            !uniAnalyzer.isHumanID(uniSet, id2)) {
            return true;
        }
        return false;
    }
    
    /**
     * A quick check method to find predicated FI partners for a provided list of protein 
     * identifiers.
     * @throws Exception
     */
    public void checkPredicatedPartners() throws Exception {
        String[] pathwayIntFileNames = new String[] {
                "ReactomeInteractions020507.txt",
                "PantherInteractions020507.txt",
                "INOHInteractions020507.txt",
                "NciNatureCuratedInteractions020507.txt",
                "NciNatureBiCartaInteractions020507.txt",
                "KEGGInteractions020507.txt",
                "CellMapInteractions020507.txt"
        };
        Set<String> pathwayInteractions = new HashSet<String>();
        Set<String> humanInteractions = new HashSet<String>();
        for (String pInt : pathwayIntFileNames) {
            Set<String> set = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + pInt);
            pathwayInteractions.addAll(set);
        }
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = acIdMap.keySet();
        filterNonHumanIds(pathwayInteractions, uniSet, uniAnalyzer);
        pathwayInteractions = filterRedundencyInteractions(pathwayInteractions, acIdMap);
        String[] ppiFileNames = new String[] {
                "BINDInteractions020507.txt",
                "IntActInteractions020507.txt",
                "HPRDInteractions020507.txt",
                "results/interaction/OrthoInteractions.txt",
                "results/interaction/YeastInteractions.txt",
                "ScoredPairsFromGOBP.txt"
        };
        String dirName = FIConfiguration.getConfiguration().get("RESULT_DIR");
        for (String ppi : ppiFileNames) {
            Set<String> set = null;
            if (ppi.contains("/"))
                set = fu.loadInteractions(ppi);
            else
                set = fu.loadInteractions(dirName + ppi);
            humanInteractions.addAll(set);
        }
        filterNonHumanIds(humanInteractions, uniSet, uniAnalyzer);
        humanInteractions = filterRedundencyInteractions(humanInteractions, acIdMap);
        System.out.println("Total pairs: " + humanInteractions.size());
        humanInteractions.removeAll(pathwayInteractions);
        System.out.println("Remove pathway FIs: " + humanInteractions.size());
        Set<String> filterPairs = filterBasedOnScores(humanInteractions,
                                                      new Double(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE")));
        System.out.println("Scored pairs: " + filterPairs.size());
        String[] checkIds = new String[] {
                "P42224",
                "P23458",
                "O60674",
                "P38484",
                "P15260"
        };
        Set<String> partners = new HashSet<String>();
        for (String id : checkIds) {
            for (String pair : filterPairs) {
                if (pair.contains(id)) {
                    partners.add(pair);
                }
            }
            System.out.println(id + ": " + partners.size());
            for (String pair : partners)
                System.out.println("\t" + pair);
            partners.clear();
        }
    }
    
    public void checkRafflaellaIDLst() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "IDListForRaffaella.txt";
        // Need to read in all lists
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        List<String> ids = new ArrayList<String>();
        String line = null;
        while ((line = fu.readLine()) != null)
            ids.add(line);
        fu.close();
        String egfr = "P00533";
        // Need to construct a interaction set
        Set<String> interactions = new HashSet<String>();
        int compare = 0;
        for (String id : ids) {
            compare = id.compareTo(egfr);
            if (compare < 0)
                interactions.add(id + " " + egfr);
            else
                interactions.add(egfr + " " + id);
        }
        // Use Ids
        Map<String, Value> valueMap = generateValues(interactions);
        // Assign values
        new WEKADataAnalyzer().generateDataSet(valueMap);
        // Need an empty dataset
        Instances dataset = new WEKAResultAnalyzer().createDataSet();
        double cutoff = 0.5d;
        String key = null;
        for (String id : ids) {
            compare = id.compareTo(egfr);
            if (compare < 0)
                key = id + " " + egfr;
            else
                key = egfr + " " + id;
            Value value = valueMap.get(key);
            double[] prob = calculateProbabilities(value, dataset);
            if (prob[0] > cutoff)
                value.functionalInteraction = true;
            else
                value.functionalInteraction = false;
            System.out.println(key + ": " + prob[0] + ", " + 
                               value.functionalInteraction + ", " +
                               value.humanInteraction + ", " + 
                               value.orthoInteraction + ", " +
                               value.yeastInteraction + ", " + 
                               value.geneExp + ", " +
                               value.goBPSemSimilarity);
        }
    }
    
    /**
     * This method is used to generate pairs from GO annotated proteins whose NBC scores
     * are higher than a certain of cutoff value.
     * @throws Exception
     */
    public void generateScoredPairsFromGO() throws Exception {
        GODataAnalyzer goAnalyzer = new GODataAnalyzer();
        Set<String> goIds = goAnalyzer.getAnnotatedGOBPProteins();
        System.out.println("Total GO IDs: " + goIds.size());
        // Use SwissProt ids only
        Set<String> goSwissProtIds = new HashSet<String>();
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> swissIds = uniAnalyzer.loadSwissProtIDsMap();
        for (String id : goIds) {
            String mapped = swissIds.get(id);
            if (mapped != null)
                goSwissProtIds.add(mapped);
        }
        System.out.println("Filtered to SwissProt IDs: " + goSwissProtIds.size());
        List<String> idList = new ArrayList<String>(goSwissProtIds);
        int compare = 0;
        Set<String> pairs = new HashSet<String>();
        Set<String> passedPairs = new HashSet<String>();
        int stepSize = 1000000;
        long time1 = System.currentTimeMillis();
        for (int i = 0; i < idList.size() - 1; i++) {
            String id1 = idList.get(i);
            for (int j = i + 1; j < idList.size(); j++) {
                String id2 = idList.get(j); 
                compare = id1.compareTo(id2);
                if (compare < 0)
                    pairs.add(id1 + " " + id2);
                else if (compare > 0)
                    pairs.add(id2 + " " + id1);
            }
            if (pairs.size() >= stepSize) {
                // Want to have a little conserved
                passedPairs.addAll(filterBasedOnScores(pairs, 0.5d));
                pairs.clear();
                long time2 = System.currentTimeMillis();
                System.out.println("*** Time for one loop: " + (time2 - time1));
                time1 = time2;
            }
        }
        if (pairs.size() > 0)
            passedPairs.addAll(filterBasedOnScores(pairs, 0.5d));
        System.out.println("Total passed pairs: " + passedPairs.size());
        fu.outputSet(passedPairs, FIConfiguration.getConfiguration().get("RESULT_DIR") + "ScoredPairsFromGOBP_041008.txt");
    }
    
    /**
     * Use this method to filter a list of protein pairs for a passed cutoff values.
     * @param pairs
     * @param cutoff
     * @return
     * @throws Exception
     */
    public Set<String> filterBasedOnScores(Set<String> pairs,
                                            double cutoff) throws Exception {
        Map<String, Value> valueMap = generateValues(pairs);
        // Assign values
        new WEKADataAnalyzer().generateDataSet(valueMap);
        // Need an empty dataset
        Instances dataset = new WEKAResultAnalyzer().createDataSet();
        Set<String> passedPairs = new HashSet<String>();
        for (Iterator<String> it = valueMap.keySet().iterator(); it.hasNext();) {
            String pair = it.next();
            Value value = valueMap.get(pair);
            classify(value, dataset, cutoff);
            if (value.functionalInteraction) {
                passedPairs.add(pair);
            }
        }
        return passedPairs;
    }

    /**
     * Calculate true and false probabilities for the specified value in the specified dataset.
     * @param value
     * @param dataset
     * @return
     * @throws Exception
     */
    public double[] calculateProbabilities(Value value,
                                           Instances dataset) throws Exception {
        // Construct a new Instance to be fed into the classified
        Instance instance = new Instance(dataset.numAttributes());
        instance.setDataset(dataset);
        instance.setClassMissing();
        if (Boolean.TRUE.equals(value.humanInteraction))
            instance.setValue(1, 0);
        else
            instance.setValue(1, 1);
        if (Boolean.TRUE.equals(value.orthoInteraction))
            instance.setValue(2, 0);
        else
            instance.setValue(2, 1);
        if (Boolean.TRUE.equals(value.yeastInteraction))
            instance.setValue(3, 0);
        else
            instance.setValue(3, 1);
//        //Newly added by xin 
//        if (Boolean.TRUE.equals(value.unambGenewaysPPI))
//            instance.setValue(4,0);
//        else
//            instance.setValue(4,1);
//        if (Boolean.TRUE.equals(value.ambGenewaysPPI))
//                instance.setValue(5,0);
//            else
//                instance.setValue(5,1);
//        if (Boolean.TRUE.equals(value.unambGenewaysFlybaseDmePPI))
//            instance.setValue(6,0);
//        else
//            instance.setValue(6,1);
//        if (Boolean.TRUE.equals(value.f347274))
//            instance.setValue(7,0);
//        else
//            instance.setValue(7,1);
//        if (Boolean.TRUE.equals(value.f364033))
//            instance.setValue(8,0);
//        else
//            instance.setValue(8,1);
//        if (Boolean.TRUE.equals(value.f342698))
//            instance.setValue(9,0);
//        else
//            instance.setValue(9,1);
//        if (Boolean.TRUE.equals(value.f239492))
//            instance.setValue(10,0);
//        else
//            instance.setValue(10,1);
        if (value.geneExp == null)
            instance.setValue(4, 2);
        else if (value.geneExp.equals("pos"))
            instance.setValue(4, 0);
        else if (value.geneExp.equals("neg"))
            instance.setValue(4, 1);
        else
            instance.setValue(4, 2);
        if (value.goBPSemSimilarity != null)
            instance.setValue(5, value.goBPSemSimilarity);
        double[] prob = classifier.distributionForInstance(instance);
        return prob;
    }
    
    protected void classify(Value value, 
                            Instances dataset,
                            double cutoff) throws Exception {
        double[] prob = calculateProbabilities(value, dataset);
        if (prob[0] > cutoff)
            value.functionalInteraction = true;
        else
            value.functionalInteraction = false;
    }
    
    /**
     * Use this help method to generate Value for protein pairs so that these pairs
     * can be scored by a trained NBC.
     * @param proteinPair
     * @return
     */
    protected Map<String, Value> generateValues(Set<String> proteinPair) {
        Map<String, Value> values = new HashMap<String, Value>();
        for (String pair : proteinPair) {
            Value value = new Value();
            values.put(pair, value);
        }
        return values;
    }
    
}
