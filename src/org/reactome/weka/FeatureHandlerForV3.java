/*
 * Created on Apr 6, 2009
 *
 */
package org.reactome.weka;

import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.GODataAnalyzerV2;
import org.reactome.fi.MicroarrayDataAnalyzer;
import org.reactome.fi.PfamAnalyzer;
import org.reactome.fi.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.Value;

/**
 * This class is used to handle features used in NBC.
 * @author wgm
 *
 */
public class FeatureHandlerForV3 {
    private FileUtility fu = new FileUtility();
    
    public FeatureHandlerForV3() {
    }
    
    public List<String> getFeatureList() {
        List<String> features = new ArrayList<String>();
        //features.add("intactHumanPPI");
        //features.add("hprdHumanPPI"); 
        //features.add("biogridHumanPPI");
        features.add("humanInteraction");
        features.add("scePPI");
        features.add("celPPI");
        features.add("dmePPI");
        // Add a new feature in March, 2012 for mouse PPI.
        features.add("mousePPI");
        // As of March, 2012, this feature is not used any more since
        // many PPIs in the databases are extracted from literatures
//        features.add("genewaysPPI");
        features.add("carlosGeneExp");
        features.add("pavlidisGeneExp");
        features.add("pfamDomainInt");
        features.add("goBPSharing");
        //features.add("goMFSharing");
        //features.add("goCCSharing");
        return features;
    }
    
    /**
     * Convert protein pairs to Value objects and assign feature values to these
     * converted Value objects.
     * @throws Exception
     */
    public Map<String, Value> convertPairsToValues(Set<String> pairs, boolean isTrue) throws Exception {
        Map<String, Value> pairToValue = new HashMap<String, Value>();
        for (String pair : pairs) {
            Value value = new Value();
            pairToValue.put(pair, value);
            value.functionalInteraction = isTrue;
        }
        generateDataSetForV3(pairToValue);
        return pairToValue;
    }
    
    /**
     * This method is used to load feature to checker
     * @return
     * @throws IOException
     */
    public Map<String, PositiveChecker> loadFeatureToChecker() throws IOException {
        Map<String, Set<String>> featureToPairs = loadFeatureToPairs();
        Map<String, PositiveChecker> featureToChecker = new HashMap<String, PositiveChecker>();
        for (String feature : featureToPairs.keySet()) {
            final Set<String> featurePairs = featureToPairs.get(feature);
            PositiveChecker checker = new PositiveChecker() {
                public boolean isPositive(String pair) {
                    return featurePairs.contains(pair);
                }
            };
            featureToChecker.put(feature, checker);
        }
        // Need to figure out other no-pair features
        final PfamAnalyzer pfamAnalyzer = new PfamAnalyzer();
        PositiveChecker checker = new PositiveChecker() {
            public boolean isPositive(String pair) {
                try {
                    return pfamAnalyzer.checkIfInteracting(pair);
                }
                catch(Exception e) {}
                return false;
            }
        };
        featureToChecker.put("pfamDomainInt", checker);
        // GO BP
        final GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        final Map<String, Set<String>> proteinToBPTerms = goAnalyzer.loadProteinToGOBPTerms();
        checker = new PositiveChecker() {
            public boolean isPositive(String pair) {
                return goAnalyzer.isTermShared(pair, proteinToBPTerms);
            }
        };
        featureToChecker.put("goBPSharing", checker);
//        final Map<String, Set<String>> proteinToMFTerms = goAnalyzer.loadProteinToGOMFTerms();
//        checker = new PositiveChecker() {
//            public boolean isPositive(String pair) {
//                return goAnalyzer.isTermShared(pair, proteinToMFTerms);
//            }
//        };
//        featureToChecker.put("goMFSharing", checker);
//        final Map<String, Set<String>> proteinToCCTerms = goAnalyzer.loadProteinToGOCCTerms();
//        checker = new PositiveChecker() {
//            public boolean isPositive(String pair) {
//                return goAnalyzer.isTermShared(pair, proteinToCCTerms);
//            }
//        };
//        featureToChecker.put("goCCSharing", checker);
        return featureToChecker;
    }
    
    /**
     * This method is used to load protein pairs for features.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadFeatureToPairs() throws IOException {
        Map<String, Set<String>> featureToPairs = new HashMap<String, Set<String>>();
        // As of March, 2012, all PPIs datasets are downloaded from iRefIndex
        Set<String> humanPPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"));
        featureToPairs.put("humanInteraction", humanPPIs);
        // For PPIs mapped from other species
        Set<String> scePPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_YEAST_FILE"));
        featureToPairs.put("scePPI", scePPIs);
        Set<String> celPPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_WORM_FILE"));
        featureToPairs.put("celPPI", celPPIs);
        Set<String> dmePPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_FLY_FILE"));
        featureToPairs.put("dmePPI", dmePPIs);
        Set<String> mousePPIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_MOUSE_FILE"));
        featureToPairs.put("mousePPI", mousePPIs);
        // For Gene Expression
        MicroarrayDataAnalyzer arrayAnalyzer = new MicroarrayDataAnalyzer();
        Set<String> pavlidisGeneExp = arrayAnalyzer.loadCoExpFromPavlidis();
        featureToPairs.put("pavlidisGeneExp", pavlidisGeneExp);
        Set<String> carlosGeneExp = arrayAnalyzer.loadCoExpFromPrietoCarlos();
        featureToPairs.put("carlosGeneExp", carlosGeneExp);
        return featureToPairs;
    }
    
    @Test
    public void countInFeatures() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissIds = uniAnalyzer.loadSwissProtIds();
        System.out.println("Total Swiss IDs: " + swissIds.size());
        System.out.println("Feature\tPairs\tProtein_IDs\tSwissProt_IDs\tPercentage");
        Map<String, Set<String>> featureToPairs = loadFeatureToPairs();
        for (String feature : featureToPairs.keySet()) {
            Set<String> pairs = featureToPairs.get(feature);
            Set<String> proteinIds = InteractionUtilities.grepIDsFromInteractions(pairs);
            int size = proteinIds.size();
            proteinIds.retainAll(swissIds);
            double percentage = (double) proteinIds.size() / swissIds.size();
            System.out.println(feature + "\t" + pairs.size() + "\t" + size + "\t" + proteinIds.size() + "\t" + percentage);
        }
        // Another two features
        // Need to figure out other no-pair features
        PfamAnalyzer pfamAnalyzer = new PfamAnalyzer();
        Map<String, Set<String>> uni2PfamMap = pfamAnalyzer.getUni2PfamMap();
        int size = uni2PfamMap.size();
        Set<String> proteins = new HashSet<String>(uni2PfamMap.keySet());
        proteins.retainAll(swissIds);
        double percentage = (double) proteins.size() / swissIds.size();
        System.out.println("PfamDomainInt" + "\tNA\t" + size + "\t" + proteins.size() + "\t" + percentage);
        // GO BP
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> proteinToBPTerms = goAnalyzer.loadProteinToGOBPTerms();
        size = proteinToBPTerms.size();
        proteins = new HashSet<String>(proteinToBPTerms.keySet());
        proteins.retainAll(swissIds);
        percentage = (double) proteins.size() / swissIds.size();
        System.out.println("GOBPSharing" + "\tNA\t" + size + "\t" + proteins.size() + "\t" + percentage);
    }
    
    /**
     * Assign feature values to the passed values map.
     * @param pairToValue
     * @throws Exception
     */
    public void generateDataSetForV3(Map<String, Value> pairToValue) throws Exception {
        Map<String, PositiveChecker> featureToChecker = loadFeatureToChecker();
        Map<String, Field> featureToField = Value.convertValueFeatureToField(getFeatureList());
        for (Iterator<String> it = pairToValue.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Value value = pairToValue.get(key);
            for (String feature : featureToChecker.keySet()) {
                PositiveChecker checker = featureToChecker.get(feature);
                Boolean featureValue = checker.isPositive(key);
                Field field = featureToField.get(feature);
                field.set(value, featureValue);
            }
        }
    }
    
}
