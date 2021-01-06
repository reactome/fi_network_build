/*
 * Created on Apr 8, 2009
 *
 */
package org.reactome.fi;

import java.io.IOException;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math.stat.inference.TestUtils;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.data.GODataAnalyzerV2;
import org.reactome.data.PfamAnalyzer;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.MathUtilities;
import org.reactome.fi.util.PositiveChecker;
import org.reactome.fi.util.Value;
import org.reactome.tred.TREDAnalyzer;
import org.reactome.weka.FeatureHandlerForV3;
import org.reactome.weka.NaiveBayesClassifier;

/**
 * This class is used to analyze results from class NaiveBayesClassifier.
 * @author wgm
 *
 */
public class NBCAnalyzer {
    private static final Logger logger = Logger.getLogger(NBCAnalyzer.class);
    // Saved NBC file name
    private final String NBC_FILE_NAME = FIConfiguration.getConfiguration().get("RESULT_DIR") + "/NaiveBayesClassifier_100_Random.ser";
    private int NEGATIVE_TO_POSITIVE_RATIO = 100;
    private FileUtility fu = new FileUtility();
    
    public NBCAnalyzer() {
    }
    
    /**
     * This method uses conditional probabilities from individual features. 
     * @throws Exception
     */
    @Test
    public void calculateNBCFromReactomeBasedOnSingleFeature() throws Exception {
        NaiveBayesClassifier nbc = new NaiveBayesClassifier();
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        nbc.calculatePriorProbability(fis);
        List<String> featureList = featureHandler.getFeatureList();
        nbc.setFeatureList(featureList);
        nbc.calculateConditionalProbabilities(fis, NEGATIVE_TO_POSITIVE_RATIO);
        Map<String, Value> positivePairToValue = featureHandler.convertPairsToValues(fis, true);
        //filterPairWithOneFeatureLeast(positivePairToValue, featureList);
        // For negative
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissIds = uniAnalyzer.loadSwissProtIds();
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(swissIds, 
                                                                           fis.size() * NEGATIVE_TO_POSITIVE_RATIO, 
                                                                           fis);
        Map<String, Value> negativePairToValue = featureHandler.convertPairsToValues(randomPairs, false);
        // Check true positive rate
        double cutoff = 0.5d;
        int tp = nbc.calculatePositiveCount(positivePairToValue, cutoff);
        System.out.println("True positive rate: " + (double) tp / positivePairToValue.size());
        int fp = nbc.calculatePositiveCount(negativePairToValue, cutoff);
        System.out.println("False positive rate: " + (double) fp / negativePairToValue.size());
        nbc.checkNBC();
    }
    
    private void filterFIsBasedOnFeatureProteins(Set<String> fis,
                                                 Set<String> featureProteins) {
        System.out.println("Total FIs: " + fis.size());
        int index = 0;
        for (Iterator<String> it = fis.iterator(); it.hasNext();) {
            String fi = it.next();
            index = fi.indexOf(" ");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            if (!featureProteins.contains(id1) ||
                    !featureProteins.contains(id2))
                it.remove();
        }
        System.out.println("Total FIs after filtering: " + fis.size());
    }
    
    /**
     * This method is used to generate two files as the training data set and test data set.
     * @throws Exception
     */
    @Test
    public void generateDatasetFile() throws Exception {
        NaiveBayesClassifier nbc = new NaiveBayesClassifier();
        // Get the positive data set
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        // Keep a positive copy for filtering negative values.
        Set<String> originalFIs = new HashSet<String>(fis);
        // Feature handler
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        // Do a filter: use only pairs having at least one feature
        filterPairWithOneFeatureLeast(fis, featureToChecker);
        nbc.calculatePriorProbability(fis);
        // For negative
        Set<String> proteinIds = InteractionUtilities.grepIDsFromInteractions(fis);
        // Generate random pairs from proteins in the filtered FIs to limit the negative data set.
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(proteinIds, 
                                                                           fis.size() * NEGATIVE_TO_POSITIVE_RATIO, 
                                                                           originalFIs);
        String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "TrainingDataset.txt";
        generateDatasetFile(fis, 
                            randomPairs, 
                            featureToChecker,
                            outFileName);
        nbc.setFeatureList(featureHandler.getFeatureList());
        nbc.calPosDatasetInSeqWay(fis, featureToChecker);
        nbc.calNegDatasetInSeqWay(randomPairs, featureToChecker);
        // Load a test data set from non-Reactome pathway FIs.
        Set<String> testFIs = loadTestFIs();
        Set<String> originalTestFIs = new HashSet<String>(testFIs);
        System.out.println("test fis before filtering: " + testFIs.size());
        testFIs.removeAll(originalFIs);
        testFIs.removeAll(randomPairs);
        filterPairWithOneFeatureLeast(testFIs, featureToChecker);
        // Same as the negative training data set, create random pairs from filtered test FIs.
        Set<String> testNegativePairs = InteractionUtilities.generateRandomPairs(InteractionUtilities.grepIDsFromInteractions(testFIs), 
                                                                                 testFIs.size() * NEGATIVE_TO_POSITIVE_RATIO,
                                                                                 originalTestFIs);
        // Remove any known FIs from Reactome.
        testNegativePairs.removeAll(originalFIs);
        testNegativePairs.removeAll(randomPairs);
        outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "TestDataset.txt";
        generateDatasetFile(testFIs,
                            testNegativePairs,
                            featureToChecker,
                            outFileName);
        // Check true positive rate
        for (int i = 1; i < 10; i++) {
            double cutoff = 0.1d * i;
            System.out.println("Cutoff: " + cutoff);
            int tp = nbc.calPosCountInSeqWay(testFIs, featureToChecker, cutoff);
            System.out.println("True positive rate: " + (double) tp / testFIs.size());
            int fp = nbc.calPosCountInSeqWay(testNegativePairs, featureToChecker, cutoff);
            System.out.println("False positive rate: " + (double) fp / testNegativePairs.size());
        }
    }
    
    private void generateDatasetFile(Set<String> posPairs,
                                     Set<String> negPairs,
                                     Map<String, PositiveChecker> featureToChecker,
                                     String outFileName) throws IOException {
        fu.setOutput(outFileName);
        List<String> featureList = new ArrayList<String>(featureToChecker.keySet());
        StringBuilder builder = new StringBuilder();
        builder.append("ProteinPair\tFunctionalInteraction");
        for (String feature : featureList) {
            builder.append("\t").append(feature);
        }
        fu.printLine(builder.toString());
        builder.setLength(0);
        boolean value = false;
        for (String pair : posPairs) {
            builder.append(pair).append("\ttrue");
            for (String feature : featureList) {
                PositiveChecker checker = featureToChecker.get(feature);
                value = checker.isPositive(pair);
                builder.append("\t").append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        for (String pair : negPairs) {
            builder.append(pair).append("\tfalse");
            for (String feature : featureList) {
                PositiveChecker checker = featureToChecker.get(feature);
                value = checker.isPositive(pair);
                builder.append("\t").append(value);
            }
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
    
    /**
     * This method is used to check correlation between two data sets based on Reactome Fis.
     * @throws Exception
     */
    @Test
    public void checkCorrelationAmongFeaturse() throws Exception {
        // Get the positive data set
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        
        // Feature handler
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        List<String> featureList = featureHandler.getFeatureList();
        Map<String, Value> positivePairToValue = featureHandler.convertPairsToValues(fis, true);
        Map<String, Set<String>> featureToPairs = featureHandler.loadFeatureToPairs();
        Map<String, Field> featureToField = Value.convertValueFeatureToField(featureList);
        for (int i = 0; i < featureList.size() - 1; i++) {
            String feature1 = featureList.get(i);
            Field field1 = featureToField.get(feature1);
            for (int j = i + 1; j < featureList.size(); j++) {
                String feature2 = featureList.get(j);
                Field field2 = featureToField.get(feature2);
                long c1 = 0;
                long c2 = 0;
                long shared = 0;
                for (String pair : positivePairToValue.keySet()) {
                    Value value = positivePairToValue.get(pair);
                    boolean v1 = (Boolean) field1.get(value);
                    if (v1)
                        c1 ++;
                    boolean v2 = (Boolean) field2.get(value);
                    if (v2)
                        c2 ++;
                    if (v1 && v2)
                        shared ++;
                }
                //                double pvalue = MathUtilities.calculateHypergeometricPValue(fis.size(), 
                //                                                                            c1,
                //                                                                            c1, 
                //                                                                            shared);
                long counts[][] = new long[][] {
                        {shared, c2 - shared},
                        {c1 - shared, fis.size() - c1 - c2 + shared}
                };
                double pvalue = TestUtils.chiSquareTest(counts);
                System.out.println(feature1 + "\t" + feature2 + "\t" + pvalue);
            }
        }
    }
    
    /**
     * This method is used to calculate NBC based on Reactome data set.
     * @throws Exception
     */
    @Test
    public void calculateNBCBasedOnReactome() throws Exception {
        _calculateNBCBasedOnReactome(true, true, 10, null);
    }
    
    @Test
    public void calculateROCPoints() throws Exception {
        String rocFileName = FIConfiguration.getConfiguration().get("ROC_CURVE_FILE");
        _calculateNBCBasedOnReactome(false, 
                                     false,
                                     100, 
                                     rocFileName);
    }
    
    private void _calculateNBCBasedOnReactome(boolean needCheckNBC,
                                              boolean needSave,
                                              int thresholdPoints,
                                              String rocPointFile) throws Exception {
        // The engine
        NaiveBayesClassifier nbc = new NaiveBayesClassifier();
        // Get the positive data set
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        
        // Feature handler
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        List<String> featureList = featureHandler.getFeatureList();
        nbc.setFeatureList(featureList);
        Map<String, Value> positivePairToValue = featureHandler.convertPairsToValues(fis, true);
        // Do a filter: use only pairs having at least one feature
        filterPairWithOneFeatureLeast(positivePairToValue, featureList);
        // For negative
        Set<String> filteredFIs = positivePairToValue.keySet();
        Set<String> proteinIds = InteractionUtilities.grepIDsFromInteractions(filteredFIs);
        nbc.calculatePriorProbability(filteredFIs);
        nbc.calculatePositiveDataset(positivePairToValue.values());
        // Generate random pairs from proteins in the filtered FIs to limit the negative data set.
        Set<String> randomPairs = InteractionUtilities.generateRandomPairs(proteinIds, 
                                                                           fis.size() * NEGATIVE_TO_POSITIVE_RATIO, 
                                                                           fis);
        System.out.println("Total negative: " + randomPairs.size());
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        // The following filter, which is used for the positive data set, should NOT be applied to the
        // negative data set. Otherwise, the odds ratio will be decreased a lot, which also make the trained NBC
        // not useful at all. This may need to be discussed. -- Guanming on April 8, 2009.
        //        filterPairWithOneFeatureLeast(randomPairs, featureToChecker);
        //System.out.println("After filtering: " + randomPairs.size());
        //calculateNegativeDataset(negativePairToValue.values());
        nbc.calNegDatasetInSeqWay(randomPairs, 
                                  featureToChecker);
        // Load a test data set from non-Reactome pathway FIs.
        Set<String> testFIs = loadTestFIs();
        System.out.println("test fis before filtering: " + testFIs.size());
        testFIs.removeAll(fis);
        System.out.println("Removing fis in the training data set: " + testFIs.size());
        // Same as the positive training data set, do a filtering too.
        filterPairWithOneFeatureLeast(testFIs, featureToChecker);
        System.out.println("Total test fis: " + testFIs.size());
        // Same as the negative training data set, create random pairs from filtered test FIs.
        Set<String> testNegativePairs = InteractionUtilities.generateRandomPairs(InteractionUtilities.grepIDsFromInteractions(testFIs), 
                                                                                 testFIs.size() * NEGATIVE_TO_POSITIVE_RATIO,
                                                                                 testFIs);
        // Remove any known FIs from Reactome.
        testNegativePairs.removeAll(fis);
        System.out.println("Total negative fis for AUC: " + testNegativePairs.size());
        // Check true positive rate
        FileUtility rocFU = null;
        if (rocPointFile != null) {
            rocFU = new FileUtility();
            rocFU.setOutput(rocPointFile);
        }
        System.out.println("Cutoff\tFalse_Positive_Rate\tTrue_Positive_Rate");
        if (rocFU != null)
            rocFU.printLine("Cutoff\tFalse_Positive_Rate\tTrue_Positive_Rate");
        int totalPoints = thresholdPoints;
        double step = 1.0 / totalPoints;
        nbc.enableCache();
        for (int i = 0; i < totalPoints + 1; i++) {
            double cutoff = 0.0 + step * i;
            //System.out.println("Cutoff: " + cutoff);
            int tp = nbc.calPosCountInSeqWay(testFIs, featureToChecker, cutoff);
            double tpr = (double) tp / testFIs.size();
            //System.out.println("True positive rate: " + (double) tp / testFIs.size());
            int fp = nbc.calPosCountInSeqWay(testNegativePairs, featureToChecker, cutoff);
            double fpr = (double) fp / testNegativePairs.size();
            System.out.println(cutoff + "\t" + fpr + "\t" + tpr);
            if (rocFU != null)
                rocFU.printLine(cutoff + "\t" + fpr + "\t" + tpr);
            //System.out.println("False positive rate: " + (double) fp / testNegativePairs.size());
            //             The following should NOT be used. Using an independent test data should be reliable than
            //             10-fold cross-validation.
            //            int tp = nbc.calculatePositiveCount(positivePairToValue, cutoff);
            //            System.out.println("True positive rate: " + (double) tp / positivePairToValue.size());
            //            int fp = nbc.calPosCountInSeqWay(randomPairs, featureToChecker, cutoff);
            //            System.out.println("False positive rate: " + (double) fp / randomPairs.size());
        }
        if (rocFU != null)
            rocFU.close();
        //        // Used as control to generate ARFF file.
        //        ARFFGenerator generator = new ARFFGenerator();
        //        generator.setIsForVersion3(true);
        //        Map<String, Value> pairToValue = new HashMap<String, Value>(positivePairToValue);
        //        pairToValue.putAll(negativePairToValue);
        //        generator.exportDataForWEKA(FIConfiguration.getConfiguration().get("RESULT_DIR + "Reactome040809.arff",
        //                                    pairToValue);
        //        generator.exportDataForWEKA(FIConfiguration.getConfiguration().get("RESULT_DIR + "Test040809.arff",
        //                                    testPairToValue);
        if (needCheckNBC)
            nbc.checkNBC();
        if (needSave)
            nbc.saveToFile(NBC_FILE_NAME);
    }
    
    private void filterPairWithOneFeatureLeast(Map<String, Value> pairToValue,
                                               List<String> featureList) throws IllegalAccessException {
        System.out.println("Total pairs: " + pairToValue.size());
        Map<String, Field> featureToField = Value.convertValueFeatureToField(featureList);
        boolean isValid = false;
        for (Iterator<String> it = pairToValue.keySet().iterator(); it.hasNext();) {
            String pair = it.next();
            Value value = pairToValue.get(pair);
            isValid = false;
            for (String feature : featureToField.keySet()) {
                Field field = featureToField.get(feature);
                Boolean featureValue = (Boolean) field.get(value);
                if (featureValue != null && featureValue) {
                    isValid = true;
                    break;
                }
            }
            if (isValid)
                continue;
            it.remove();
        }
        System.out.println("After filtering: " + pairToValue.size());
    }
    
    private void filterPairWithOneFeatureLeast(Set<String> pairs,
                                               Map<String, PositiveChecker> featureToChecker) {
        boolean isValid = false;
        for (Iterator<String> it = pairs.iterator(); it.hasNext();) {
            String pair = it.next();
            isValid = false;
            for (String feature : featureToChecker.keySet()) {
                PositiveChecker checker = featureToChecker.get(feature);
                if (checker.isPositive(pair)) {
                    isValid = true;
                    break;
                }
            }
            if (isValid)
                continue;
            it.remove();
        }
    }
    
    private Set<String> loadTestFIs() throws IOException {
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("KEGG_FI_FILE"),
                FIConfiguration.getConfiguration().get("NCI_PID_FI_FILE"),
                FIConfiguration.getConfiguration().get("NCI_PID_BIOCARTA_FI_FILE"),
                FIConfiguration.getConfiguration().get("PANTHER_FI_FILE")
                // Don't use TRED as a test FIs since TRED interactions are not related to pathways
                // As of March, 2012, pathways from cell map is not used any more because of its low quality!
        };
        Set<String> fis = new HashSet<String>();
        for (String fileName : fileNames) {
            Set<String> fis1 = fu.loadInteractions(fileName);
            fis.addAll(fis1);
        }
        return fis;
    }
    
    public NaiveBayesClassifier loadSavedNBC() throws IOException, ClassNotFoundException {
        NaiveBayesClassifier nbc = (NaiveBayesClassifier) fu.loadObject(NBC_FILE_NAME);
        return nbc;
    }
    
    @Test
    public void testLoadSavedNBC() throws Exception {
        NaiveBayesClassifier classifier = loadSavedNBC();
        System.out.println("Prior probability: " + classifier.getPriorProbability());
        classifier.checkNBC();
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        double cutoff = 0.5d;
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, Value> pairToValue = featureHandler.convertPairsToValues(fis, true);
        filterPairWithOneFeatureLeast(pairToValue, 
                                      featureHandler.getFeatureList());
        int tp = classifier.calculatePositiveCount(pairToValue, cutoff);
        System.out.println("True positive rate: " + (double) tp / pairToValue.size());
    }
    
    /**
     * This method is used to check how big a random set should be. The check is based on
     * standard deviation size.
     * @throws Exception
     */
    @Test
    public void checkSizeOfRandomPairs() throws Exception {
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        Set<String> fis = analyzer.loadFIsFromFile();
        NaiveBayesClassifier nbc = new NaiveBayesClassifier();
        nbc.calculatePriorProbability(fis);
        List<String> featureList = featureHandler.getFeatureList();
        nbc.setFeatureList(featureList);
        Map<String, Value> positivePairToValue = featureHandler.convertPairsToValues(fis, true);
        // Filter out Values without any value
        filterPairWithOneFeatureLeast(positivePairToValue, featureList);
        nbc.calculatePositiveDataset(positivePairToValue.values());
        // For negative: generate from random pairs from SwissProts
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissIds = uniAnalyzer.loadSwissProtIds();
        // Generate a list of size:
        double[] ratios = new double[] {
                1.0d,
                10.0d,
                100.0d
        };
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        for (double ratio : ratios) {
            int size = (int)(fis.size() * ratio);
            System.out.println("Random size: " + size + " (" + ratio + ")");
            DescriptiveStatistics stat = new DescriptiveStatistics();
            for (int i = 0; i < 5; i++) { // Try to run 5 times
                Set<String> randomPairs = InteractionUtilities.generateRandomPairs(swissIds, 
                                                                                   size,
                                                                                   fis);
                nbc.calNegDatasetInSeqWay(randomPairs, featureToChecker);
                // Check true positive rate
                double cutoff = 0.5d;
                int tp = nbc.calculatePositiveCount(positivePairToValue, cutoff);
                System.out.println("True positive rate: " + (double) tp / positivePairToValue.size());
                int fp = nbc.calPosCountInSeqWay(randomPairs, 
                                                 featureToChecker,
                                                 cutoff);
                double fpRate = (double) fp / randomPairs.size();
                stat.addValue(fpRate);
                System.out.println("False positive rate: " + (double) fp / randomPairs.size());
            }
            System.out.println("Average false positive rate: " + stat.getMean() + " +- " + stat.getStandardDeviation());
            System.out.println();
        }
    }
    
    /**
     * Use this method to dump predicated FIs into an external file.
     * @throws Exception
     */
    @Test
    public void generatePredictedFIs() throws Exception {
        NaiveBayesClassifier classifier = loadSavedNBC();
        // Get pairs from all features
        Set<String> allPairs = loadPairForPrediction();
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        //String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "PredictedFIs.txt";
        //String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "PredictedFIsHumanPPIsGeneExp.txt";
        String outFileName = FIConfiguration.getConfiguration().get("PREDICTED_FI_FILE");
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        double cutoff = new Double(FIConfiguration.getConfiguration().get("CUT_OFF_VALUE"));
        logger.info("Cutoff: " + cutoff);
        Set<String> predicted = new HashSet<String>();
        for (String pair : allPairs) {
            Double score = classifier.calculateScore(pair, featureToChecker);
            if (score >= cutoff) {
                predicted.add(pair);
                fu.printLine(pair);
            }
        }
        fu.close();
        Set<String> pathwayFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        countForDBPathwayFIsAndPredFIs(pathwayFIs, predicted);
    }
    
    /**
     * @return
     * @throws Exception
     */
    public Set<String> loadPairForPrediction() throws Exception {
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Set<String> allPairs = loadFeaturePairs(featureHandler);
        Set<String> pathwayFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        allPairs.removeAll(pathwayFIs); // Don't count pathway FIs.
        return allPairs;
    }
    
    /**
     * This method is used to check cutoff values to get enough size of 
     * combined FI network.
     * @throws Exception
     */
    @Test
    public void checkCutoffValueForPredictedFIs() throws Exception {
        NaiveBayesClassifier classifier = loadSavedNBC();
        // Get pairs from all features
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Set<String> allPairs = loadPairForPrediction();
        System.out.println("Total pairs: " + allPairs.size());
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        Set<String> pathwayFIs = new FIFileAnalyzer().loadPathwayAndTFTargetFIs();
        //        double cutoffs[] = new double[]{0.25, 0.35, 0.45, 0.55};
        //        for (double cutoff : cutoffs) {
        for (int i = 1; i < 10; i++) {
            double cutoff = 0.1d * i;
            System.out.println("Cutoff: " + cutoff);
            Set<String> predicted = new HashSet<String>();
            for (String pair : allPairs) {
                Double score = classifier.calculateScore(pair, featureToChecker);
                if (score >= cutoff) {
                    predicted.add(pair);
                }
            }
            // Make copies in case pathways set are changed.
            countForDBPathwayFIsAndPredFIs(new HashSet<String>(pathwayFIs),
                                           predicted);
        }
    }
    
    private Set<String> loadFeaturePairs(FeatureHandlerForV3 featureHandler) throws IOException {
        Map<String, Set<String>> featureToPairs = featureHandler.loadFeatureToPairs();
        Set<String> allPairs = new HashSet<String>();
        for (String feature : featureToPairs.keySet()) {
            //            if (feature.equals("intactHumanPPI") ||
            //                feature.equals("hprdHumanPPI") ||
            //                feature.equals("biogridHumanPPI") ||
            //                feature.equals("pavlidisGeneExp") ||
            //                feature.equals("carlosGeneExp")) {
            Set<String> featurePairs = featureToPairs.get(feature);
            allPairs.addAll(featurePairs);
            //            }
        }
        // Add pairs from BP and Domain interactions. If a pair can be predicted from one of these
        // two data sets, it much be shared.
        Set<String> bpDomainPairs = fu.loadInteractions(FIConfiguration.getConfiguration().get("BP_DOMAIN_SHARED_PAIRS"));
        allPairs.addAll(bpDomainPairs);
        return allPairs;
    }
    
    @Test
    public void checkFeaturePairMatrixForPredictedFIs() throws Exception {
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        Map<String, Set<String>> featureToPairs = featureHandler.loadFeatureToPairs();
        
        Set<String> predictedFIs = new FIFileAnalyzer().loadPredictedFIs();
        System.out.println("Total predicted FIs: " + predictedFIs.size());
        
        Map<String, Set<String>> featureToPredictedFIs = new HashMap<String, Set<String>>();
        for (String feature : featureToChecker.keySet()) {
            PositiveChecker checker = featureToChecker.get(feature);
            Set<String> set = new HashSet<String>();
            for (String pair : predictedFIs) {
                if (checker.isPositive(pair))
                    set.add(pair);
            }
            featureToPredictedFIs.put(feature, set);
        }
        
        String[] features = new String[] {
                "humanInteraction",
                "dmePPI",
                "celPPI",
                "scePPI",
                "pfamDomainInt",
                "pavlidisGeneExp",
                "carlosGeneExp",
                "goBPSharing",
                "genewaysPPI",
        };
        System.out.println("Feature\tInteractions");
        Set<String> pathwayFIs = new FIFileAnalyzer().loadPathwayFIs();
        Set<String> tfTargetFIs = new TREDAnalyzer().loadTFTargetInteractions();
        
        for (String feature : features) {
            Set<String> pairs = featureToPredictedFIs.get(feature);
            Set<String> originalPairs = featureToPairs.get(feature);
            if (originalPairs == null)
                originalPairs = new HashSet<String>();
            // Want to remove pairs from pathwayFIs and tfTargetFIs
            Set<String> copy = new HashSet<String>(originalPairs);
            System.out.println("Before remove: " + copy.size());
            copy.removeAll(pathwayFIs);
            copy.removeAll(tfTargetFIs);
            System.out.println("After remove: " + copy.size());
            System.out.println(feature + "\t" + pairs.size() + "\t" + 
                    copy.size() + "\t" + 
                    (double)pairs.size() / copy.size());
        }
        // Generate a matrix
        StringBuilder builder = new StringBuilder();
        for (String feature : features) {
            builder.append("\t").append(feature);
        }
        System.out.println("\n" + builder.toString());
        builder.setLength(0);
        for (String feature : features) {
            Set<String> pairs = featureToPredictedFIs.get(feature);
            builder.append(feature);
            for (String feature1 : features) {
                Set<String> pairs1 = featureToPredictedFIs.get(feature1);
                if (feature1.equals(feature)) {
                    builder.append("\t");
                    continue;
                }
                Set<String> shared = InteractionUtilities.getShared(pairs, pairs1);
                builder.append("\t").append(shared.size());
            }
            System.out.println(builder.toString());
            builder.setLength(0);
        }
    }
    
    @Test
    public void checkPPIsOverlapWithDomainInt() throws Exception {
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, Set<String>> featureToPairs = featureHandler.loadFeatureToPairs();
        Set<String> humanPPIs = featureToPairs.get("humanInteraction");
        Map<String, PositiveChecker> featureToChecker = featureHandler.loadFeatureToChecker();
        PositiveChecker checker = featureToChecker.get("pfamDomainInt");
        // From PPIs to domain PPIs
        int count = 0;
        for (String pair : humanPPIs) {
            if (checker.isPositive(pair))
                count ++;
        }
        System.out.println("Total human PPIs: " + humanPPIs.size());
        System.out.println("Having domain ints: " + count);
        System.out.println("Percentage: " + (double)count / humanPPIs.size());
        // From domain to PPis: do a random sampling
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(humanPPIs);
        Set<String> randomIds = MathUtilities.randomSampling(ids, 5000);
        count = 0;
        int ppiCount = 0;
        List<String> list = new ArrayList<String>(randomIds);
        Collections.sort(list);
        for (int i = 0; i < list.size() - 1; i ++) {
            String id = list.get(i);
            for (int j = i + 1; j < list.size(); j++) {
                String id1 = list.get(j);
                String pair = id + " " + id1;
                if (checker.isPositive(pair)) {
                    count ++;
                    if (humanPPIs.contains(pair))
                        ppiCount ++;
                }
            }
        }
        System.out.println("Domain ints: " + count);
        System.out.println("PPIs from domain int: " + ppiCount);
        System.out.println("Percentage: " + (double) ppiCount / count);
    }
    
    /**
     * This method is used to check coexpression data set in PPIs.
     * @throws Exception
     */
    @Test
    public void checkCoExpAndCoLocInPPIs() throws Exception {
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, Set<String>> featureToPairs = featureHandler.loadFeatureToPairs();
        Set<String> ppiPairs = new HashSet<String>();
        Set<String> ppiFeatures = new HashSet<String>();
        ppiFeatures.add("humanInteraction");
        ppiFeatures.add("scePPI");
        ppiFeatures.add("celPPI");
        ppiFeatures.add("dmePPI");
        ppiFeatures.add("genewaysPPI");
        // Load Gene Expression
        Set<String> geneExpPairs = new HashSet<String>();
        Set<String> geneExpFeatures = new HashSet<String>();
        geneExpFeatures.add("pavlidisGeneExp");
        geneExpFeatures.add("carlosGeneExp");
        for (String feature : featureToPairs.keySet()) {
            System.out.println(feature + ": " + featureToPairs.get(feature).size());
            if (ppiFeatures.contains(feature)) {
                Set<String> featurePairs = featureToPairs.get(feature);
                ppiPairs.addAll(featurePairs);
            }
            else if (geneExpFeatures.contains(feature)) {
                Set<String> featurePairs = featureToPairs.get(feature);
                geneExpPairs.addAll(featurePairs);
            }
        }
        Set<String> ppiIds = InteractionUtilities.grepIDsFromInteractions(ppiPairs);
        System.out.println("Total ppi pairs: " + ppiPairs.size());
        System.out.println("Total ppi ids: " + ppiIds.size());
        Set<String> geneExpIds = InteractionUtilities.grepIDsFromInteractions(geneExpPairs);
        System.out.println("Total gene expression pairs: " + geneExpPairs.size());
        System.out.println("Total gene expression ids: " + geneExpIds.size());
        // Check shared
        Set<String> shared = InteractionUtilities.getShared(ppiPairs, geneExpPairs);
        System.out.println("Total shared: " + shared.size());
        double percentage = (double) shared.size() / ppiPairs.size();
        System.out.println("Percentage: " + percentage);
        // Check ids contains by gene expression;
        Set<String> ppisInExps = new HashSet<String>();
        int index = 0;
        for (String pair : ppiPairs) {
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            if (geneExpIds.contains(id1) && geneExpIds.contains(id2)) {
                ppisInExps.add(pair);
            }
        }
        System.out.println("PPIs having both ids in gene exp: " + ppisInExps.size());
        shared = InteractionUtilities.getShared(ppisInExps, geneExpPairs);
        System.out.println("Total shared after filtering: " + shared.size());
        percentage = (double) shared.size() / ppisInExps.size();
        System.out.println("Percentage: " + percentage);
        
        // Load predicted FIs
        Set<String> predictedFIs = new FIFileAnalyzer().loadPredictedFIs();
        System.out.println("Total predicted FIs: " + predictedFIs.size());
        Set<String> ppisInPredicted = InteractionUtilities.getShared(predictedFIs, ppiPairs);
        System.out.println("Total predicted FIs having ppi: " + ppisInPredicted.size());
        Set<String> ppisInExpInPredicted = InteractionUtilities.getShared(ppisInExps, predictedFIs);
        System.out.println("Total predicited FIs having ppis for gene exp: " + ppisInExpInPredicted.size());
        shared = InteractionUtilities.getShared(ppisInExpInPredicted, geneExpPairs);
        System.out.println("Shared with gene exp: " + shared.size());
        System.out.println("Percentage for ppi: " + (double) shared.size() / ppisInPredicted.size());
        System.out.println("Percentage for ppi only for gene exp ids: " + (double) shared.size() / ppisInExpInPredicted.size());
        
        // The following is used to check if two proteins are in the same compartments
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> proteinToCC = goAnalyzer.loadProteinToGOCCTerms();
        Set<String> ppisInGO = new HashSet<String>();
        for (String pair : ppiPairs) {
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            if (proteinToCC.keySet().contains(id1) && proteinToCC.keySet().contains(id2)) {
                ppisInGO.add(pair);
            }
        }
        System.out.println("\n\nPPIs having two ids in GO CC: " + ppisInGO.size());
        System.out.println(" from total PPI pairs: " + ppiPairs.size());
        int count = 0;
        for (String pair : ppisInGO) {
            if (goAnalyzer.isTermShared(pair, proteinToCC))
                count ++;
        }
        System.out.println("PPIs having term shared: " + count);
        percentage = (double) count / ppisInGO.size();
        System.out.println("Percentage: " + percentage);
        Set<String> ppiInGOInPredicted = InteractionUtilities.getShared(ppisInPredicted, ppisInGO);
        System.out.println("PPIs having two ids in GO CC in predicted FIs: " + ppiInGOInPredicted.size());
        System.out.println(" from total PPIs in predicted FIs: " + ppisInPredicted.size());
        count = 0;
        for (String pair : ppiInGOInPredicted) {
            if (goAnalyzer.isTermShared(pair, proteinToCC))
                count ++;
        }
        System.out.println("Predicted PPIs having CC term shared: " + count);
        percentage = (double) count / ppiInGOInPredicted.size();
        System.out.println("Percentage: " + percentage);
    }
    
    /**
     * This method is used to generate a file containing shared protein pairs from
     * GO BP term sharing and Domain interactions. A pair in this file can be used
     * for scoring checking.
     * @throws Exception
     */
    @Test
    public void checkSharedBPPairAndDomainPair() throws Exception {
        // Used to check proteins from GO BP annotations
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        Map<String, Set<String>> goToBPTerms = goAnalyzer.loadProteinToGOBPTerms();
        Set<String> goIDs = goToBPTerms.keySet();
        List<String> goIdList = new ArrayList<String>(goIDs);
        Collections.sort(goIdList); // To avoid comparing
        System.out.println("Total GO ids: " + goIdList.size());
        // Check with GO ids
        long time1 = System.currentTimeMillis();
        PfamAnalyzer domainAnalyzer = new PfamAnalyzer();
        Set<String> shared = new HashSet<String>();
        for (int j = 0; j < goIdList.size() - 1; j++) {
            String id1 = goIdList.get(j);
            for (int k = j + 1; k < goIdList.size(); k ++) {
                String id2 = goIdList.get(k);
                String pair = id1 + "\t" + id2;
                if (goAnalyzer.isTermShared(pair, goToBPTerms) &&
                        domainAnalyzer.checkIfInteracting(pair)) {
                    shared.add(pair);
                }
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total shared: " + shared.size());
        fu.saveInteractions(shared, FIConfiguration.getConfiguration().get("BP_DOMAIN_SHARED_PAIRS"));
    }
    
    private void countForDBPathwayFIsAndPredFIs(Set<String> pathwayFIs,
                                                Set<String> predictedFIs) throws Exception {
        Set<String> pathwayIds = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        System.out.printf("FIs from pathways: %d (%d)%n",
                          pathwayFIs.size(),
                          pathwayIds.size());
        ProteinAndInteractionCount counter = new ProteinAndInteractionCount();
        counter.countVsSwissProt(pathwayIds);
        Set<String> predictedIDs = InteractionUtilities.grepIDsFromInteractions(predictedFIs);
        System.out.printf("FIs from prediction: %d (%d)%n",
                          predictedFIs.size(),
                          predictedIDs.size());
        counter.countVsSwissProt(predictedIDs);
        // Merge them
        pathwayIds.addAll(predictedIDs);
        pathwayFIs.addAll(predictedFIs);
        System.out.printf("FIs merged: %d (%d)%n",
                          pathwayFIs.size(),
                          pathwayIds.size());
        counter.countVsSwissProt(pathwayIds);
    }
}
