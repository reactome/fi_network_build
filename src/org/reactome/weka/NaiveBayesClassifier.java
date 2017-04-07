/*
 * Created on Apr 3, 2009
 *
 */
package org.reactome.weka;

import java.io.IOException;
import java.io.Serializable;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.PositiveChecker;
import org.reactome.fi.util.Value;

/**
 * This class is used to NBC related calculation. 
 * @author wgm
 *
 */
public class NaiveBayesClassifier implements Serializable {
    private static final long serialVersionUID = 123456789L;
    // This prior probability is pre-configured
    // 0.001 is based on paper: Alexeyenko, A & Sonnhammer, ELL. Genome Research March 2, 2009.
    private double priorPosProb = 0.001d;
    private List<String> featureList;
    private Map<String, Double> conditionalProbabilities;
    // Cache scores for quick ROC calculation
    private Map<String, Double> pairToScore;
    
    public NaiveBayesClassifier() {
    }
    
    public void enableCache() {
        pairToScore = new HashMap<String, Double>();
    }
    
    /**
     * Set the features for this NBC. In the current used NBC, all features are
     * boolean features.
     * @param featureList
     */
    public void setFeatureList(List<String> featureList) {
        this.featureList = featureList;
    }
    
    public List<String> getFeatureList() {
        return this.featureList;
    }
    
    public Map<String, Double> getConditionalProbabilities() {
        return conditionalProbabilities;
    }
    
    public double getPriorProbability() {
        return this.priorPosProb;
    }
    
    /**
     * Calculate the score for a value. This method should be called after positive and negative
     * probabilities have been tested. Otherwise, an exception will be thrown.
     * @param value
     * @return
     */
    public double calculateScore(Value value) {
        if (conditionalProbabilities == null || conditionalProbabilities.size() == 0)
            throw new  IllegalStateException("NaiveBayesClassifier.calculateScore(): " +
            		                         "conditional probabiltiies have not calculated!");
        Map<String, Field> featureToField = Value.convertValueFeatureToField(featureList);
        // Calculate the positive side
        double posProbs = priorPosProb;
        // calculate the negative probability
        double negProbs = 1.0 - priorPosProb;
        try {
            String posKey = null;
            String negKey = null;
            for (String feature : featureList) {
                // Get this feature value from value
                Field field = featureToField.get(feature);
                Boolean featureValue = (Boolean) field.get(value);
                if (featureValue != null && featureValue.booleanValue()) {
                    posKey = feature + "_true|true";
                    negKey = feature + "_true|false";
                }
                else {
                    posKey = feature + "_false|true";
                    negKey = feature + "_false|false";
                }
                Double prob = conditionalProbabilities.get(posKey);
                posProbs *= prob;
                prob = conditionalProbabilities.get(negKey);
                negProbs *= prob;
            }
        }
        catch(IllegalAccessException e) {
            e.printStackTrace();
        }
        // Combine positive and negative to give a score
        return posProbs / (posProbs + negProbs);
    }
    
    public double calculateScore(String pair,
                                 Map<String, PositiveChecker> featureToChecker) {
        if (conditionalProbabilities == null || conditionalProbabilities.size() == 0)
            throw new  IllegalStateException("NaiveBayesClassifier.calculateScore(): " +
                                             "conditional probabiltiies have not calculated!");
        if (pairToScore != null && pairToScore.containsKey(pair))
            return pairToScore.get(pair);
        // Calculate the positive side
        double posProbs = priorPosProb;
        // calculate the negative probability
        double negProbs = 1.0 - priorPosProb;
        String posKey = null;
        String negKey = null;
        for (String feature : featureList) {
            // Get this feature value from value
            PositiveChecker checker = featureToChecker.get(feature);
            Boolean featureValue = checker.isPositive(pair);
//            System.out.println(feature + " -> " + featureValue);
            if (featureValue != null && featureValue.booleanValue()) {
                posKey = feature + "_true|true";
                negKey = feature + "_true|false";
            }
            else {
                posKey = feature + "_false|true";
                negKey = feature + "_false|false";
            }
            Double prob = conditionalProbabilities.get(posKey);
            posProbs *= prob;
            prob = conditionalProbabilities.get(negKey);
            negProbs *= prob;
        }
        // Combine positive and negative to give a score
        double rtn = posProbs / (posProbs + negProbs);
        if (pairToScore != null)
            pairToScore.put(pair, rtn);
        return rtn;
    }
    
    public void calculateConditionalProbabilities(Set<String> posPairs,
                                                  double ratioOfNegToPos) throws Exception {
        if (conditionalProbabilities == null)
            conditionalProbabilities = new HashMap<String, Double>();
        else
            conditionalProbabilities.clear();
        FeatureHandlerForV3 featureHandler = new FeatureHandlerForV3();
        Map<String, PositiveChecker> featureToPairs = featureHandler.loadFeatureToChecker();
        for (String feature : featureToPairs.keySet()) {
            PositiveChecker checker = featureToPairs.get(feature);
            calculateConditionalProbability(feature,
                                            checker,
                                            new HashSet<String>(posPairs),
                                            ratioOfNegToPos); // Make a copy to avoid modification
        }
    }
    
    private void calculateConditionalProbability(String featureName,
                                                 PositiveChecker posChecker, 
                                                 Set<String> posPairs,
                                                 double ratioOfNegToPos) throws IOException {
        System.out.println("Working on " + featureName);
        int originalPosSize = posPairs.size();
        calConditionalProbability(featureName,
                                  posChecker, 
                                  posPairs,
                                  true);
        // Generate negative pairs
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissIds = uniAnalyzer.loadSwissProtIds();
        Set<String> negPairs = InteractionUtilities.generateRandomPairs(swissIds, 
                                                                        (int) (posPairs.size() * ratioOfNegToPos), 
                                                                        posPairs);
        calConditionalProbability(featureName,
                                  posChecker, 
                                  negPairs,
                                  false);
    }
    
    private void calConditionalProbability(String featureName,
                                           PositiveChecker posChecker,
                                           Set<String> pairs,
                                           boolean isForPos) {
        int pos = 0;
        for (String posPair : pairs) {
            if (posChecker.isPositive(posPair)) {
                pos ++;
            }
        }
        double prob = (double) pos / pairs.size();
        String key = featureName + "_true|" + isForPos;
        conditionalProbabilities.put(key, prob);
        key = featureName + "_false|" + isForPos;
        conditionalProbabilities.put(key, 1- prob);
    }

    /**
     * This method is used to calculate all condition properties for the positive
     * data set.
     * @param positiveValues
     */
    public void calculatePositiveDataset(Collection<Value> positiveValues) {
        calculateDatasets(positiveValues, "true");
    }
    
    /**
     * Calculate probabilities for the passed values.
     * @param values
     * @param event true or false.
     */
    private void calculateDatasets(Collection<Value> values,
                                   String event) {
        if (featureList == null) {
            throw new IllegalStateException("NaiveBayesClassifier.calcualtePositiveDataset(): " +
                                            "featureList has not been assigned!");
        }
        if (conditionalProbabilities == null)
            conditionalProbabilities = new HashMap<String, Double>();
        // Need this map for reflection
        Map<String, Field> featureToField = Value.convertValueFeatureToField(featureList);
        try {
            for (String feature : featureList) {
                Field field = featureToField.get(feature);
                // Count the number of feature true
                int featureTrue = 0;
                for (Value value : values) {
                    // Do some reflection
                    Boolean tmp = (Boolean) field.get(value);
                    if (tmp != null && tmp.booleanValue())
                        featureTrue ++;
                }
                double probability = (double) featureTrue / values.size();
                String key = feature + "_true|" + event;
                conditionalProbabilities.put(key, probability);
                key = feature + "_false|" + event;
                conditionalProbabilities.put(key, 1.0 - probability);
            }
        }
        catch(IllegalAccessException e) {
            System.err.println("NaiveBayesClassifier.calculatePositiveDataset(): " + e);
        }
    }
    
    public void calculateNegativeDataset(Collection<Value> negativeValues) {
        calculateDatasets(negativeValues, "false");
    }
    
    public void calculatePriorProbability(Set<String> fis) {
        Set<String> fiIds = InteractionUtilities.grepIDsFromInteractions(fis);
        int totalIds = fiIds.size();
        int totalPairs = totalIds * (totalIds - 1) / 2;
        priorPosProb = (double) fis.size() / totalPairs;
        System.out.println("Prior probability: " + priorPosProb);
    }
    
    public void calNegDatasetInSeqWay(Set<String> pairs,
                                       Map<String, PositiveChecker> featureToChecker) throws Exception {
        calDatasetInSeqWay(pairs, featureToChecker, false);
    }

    private void calDatasetInSeqWay(Set<String> pairs,
                                    Map<String, PositiveChecker> featureToChecker,
                                    boolean isPositive) {
        if (conditionalProbabilities == null)
            conditionalProbabilities = new HashMap<String, Double>();
        for (String feature : featureToChecker.keySet()) {
            PositiveChecker checker = featureToChecker.get(feature);
            int posCount = 0;
            for (String pair : pairs) {
                boolean featureValue = checker.isPositive(pair);
                if (featureValue)
                    posCount ++;
            }
            String key = feature + "_true|" + isPositive;
            double prob = (double) posCount / pairs.size();
            conditionalProbabilities.put(key, prob);
            key = feature + "_false|" + isPositive;
            prob = 1.0 - prob;
            conditionalProbabilities.put(key, prob);
        }
    }
    
    public void calPosDatasetInSeqWay(Set<String> pairs,
                                      Map<String, PositiveChecker> featureToChecker) throws Exception {
        calDatasetInSeqWay(pairs, featureToChecker, true);
    }
    
    public int calPosCountInSeqWay(Set<String> pairs,
                                    Map<String, PositiveChecker> featureToChecker,
                                    double cutoff) throws Exception {
        int posCount = 0;
        double score;
        for (String pair : pairs) {
            score = calculateScore(pair, featureToChecker);
            if (score >= cutoff)
                posCount ++;
        }
        return posCount;
    }
    
    /**
     * Calculate the count for the positive pairs based on the passed pairToValue map.
     * @param pairToValue
     * @param cutoff
     * @return
     */
    public int calculatePositiveCount(Map<String, Value> pairToValue,
                                      double cutoff) {
        
        int posCount = 0;
        for (String pair : pairToValue.keySet()) {
            Value value = pairToValue.get(pair);
            double score = calculateScore(value);
            if (score >= cutoff)
                posCount ++;
        }
        return posCount;
    }
    
    /**
     * Check the learning parameters.
     * @throws IllegalAccessException
     */
    public void checkNBC() throws IllegalAccessException {
        if (conditionalProbabilities == null || conditionalProbabilities.size() == 0)
            throw new  IllegalStateException("NaiveBayesClassifier.calculateScore(): " +
                                             "conditional probabiltiies have not calculated!");
        // Check if one one feature is true
        Map<String, Value> featureToValue = constructOneFeatureValues(featureList);
        System.out.println("\nOne feature contribution:");
        for (String feature : featureToValue.keySet()) {
            Value value = featureToValue.get(feature);
            double score = calculateScore(value);
            System.out.println(feature + ": " + score);
        }
        // Check learned probabilties
        System.out.println("\nLearning probabilties:");
        List<String> keys = new ArrayList<String>(conditionalProbabilities.keySet());
        Collections.sort(keys);
        for (String key : keys) {
            Double prob = conditionalProbabilities.get(key);
            System.out.println(key + ": " + prob);
        }
    }
    
    private Map<String, Value> constructOneFeatureValues(List<String> featureList) throws IllegalAccessException {
        Map<String, Field> featureToField = Value.convertValueFeatureToField(featureList);
        Map<String, Value> featureToValue = new HashMap<String, Value>();
        for (String feature : featureList) {
            Field field = featureToField.get(feature);
            Value value = new Value();
            field.set(value, Boolean.TRUE);
            featureToValue.put(feature, value);
        }
        return featureToValue;
    }
    
    public void saveToFile(String fileName) throws Exception {
        FileUtility fu = new FileUtility();
        fu.saveObjet(this, fileName);
    }
}
