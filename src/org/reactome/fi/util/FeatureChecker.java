/*
 * Created on Mar 31, 2009
 *
 */
package org.reactome.fi.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * This method is used to check a feature's odds ratio.
 * @author wgm
 */
public class FeatureChecker {
    private Set<String> interactionSet = null;
    
    public FeatureChecker() {
    }
    
    private void loadInteractionSet() throws IOException {
        if (interactionSet == null) {
            String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
            interactionSet = new FileUtility().loadInteractions(fileName);
        }
    }
    
    /**
     * This method is used to check a feature's odd ratio. The feature has been wrapped
     * in the passed positiveChecker parameter.
     * @param posChecker
     * @throws Exception
     */
    public void checkFeatureOddsRatio(PositiveChecker posChecker) throws Exception {
        loadInteractionSet();
        checkFeatureOddsRatio(interactionSet,
                              posChecker);
    }
    
    /**
     * This method is used to check a feature's odd ratio. The data for the feature
     * is provided as a set of protein pairs.
     * @param proteinPairs
     * @throws Exception
     */
    public void checkFeatureOddsRatio(final Set<String> proteinPairs) throws Exception {
        loadInteractionSet();
        PositiveChecker isPositivev = new PositiveChecker() {
            public boolean isPositive(String pair) {
                return proteinPairs.contains(pair);
            }
        };
        System.out.println("Total checked pairs: " + proteinPairs.size());
        checkFeatureOddsRatio(interactionSet,
                              isPositivev);
    }
    
    public Set<String> generateRandomPairs(Set<String> fis) {
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(fis);
        List<String> idList = new ArrayList<String>(ids);
        int total = fis.size();
        Set<String> randomPairs = new HashSet<String>();
        int index1 = 0;
        int index2 = 0;
        Random random = new Random();
        int compare = 0;
        while (randomPairs.size() < total) {
            index1 = random.nextInt(idList.size());
            String id1 = idList.get(index1);
            index2 = random.nextInt(idList.size());
            String id2 = idList.get(index2);
            compare = id1.compareTo(id2);
            String pair = null;
            if (compare < 0)
                pair = id1 + "\t" + id2;
            else if (compare > 0)
                pair = id2 + "\t" + id1;
            if (pair == null || fis.contains(pair))
                continue;
            randomPairs.add(pair);
        }
        return randomPairs;
    }
    
    /**
     * Check a feature's odd ratio.
     * @param interactionSet
     * @param randomPairs
     * @param ppis
     */
    public void checkFeatureOddsRatio(Set<String> interactionSet,
                                      Set<String> randomPairs, 
                                      Set<String> ppis) {
        int tp = 0;
        for (String fi : interactionSet) {
            // No need to check isoforms!!!
            if (ppis.contains(fi))
                tp ++;
        }
        System.out.printf("Total: %d%n", interactionSet.size());
        double propOnPos = (double) tp / interactionSet.size();
        System.out.printf("Mapped to ppi: %d (%f)%n", 
                          tp,
                          propOnPos);
        // Check random pairs
        tp = 0;
        for (String pair : randomPairs) {
            if (ppis.contains(pair))
                tp ++;
        }
        System.out.printf("Total random: %d%n", randomPairs.size());
        double propOnNeg = (double) tp / randomPairs.size();
        System.out.printf("Mapped to random: %d (%f)%n", 
                          tp,
                          propOnNeg);
        double odds = calculateOddsRatio(propOnPos, propOnNeg);
        System.out.println("Odds: " + odds);
        System.out.println();
    }
    
    /**
     * In this method, five random pairs have been generated to create a mean
     * odds ratios.
     * @param interactionSet
     * @param ppis
     */
    private void checkFeatureOddsRatio(Set<String> interactionSet,
                                       PositiveChecker posChecker) {
        int tp = 0;
        for (String fi : interactionSet) {
            // No need to check isoforms!!!
            if (posChecker.isPositive(fi))
                tp ++;
        }
        System.out.printf("Total: %d%n", interactionSet.size());
        double propOnPos = (double) tp / interactionSet.size();
        System.out.printf("Mapped to ppi: %d (%f)%n", 
                          tp,
                          propOnPos);
        // Check random pairs
        // Need to run five times
        List<Double> oddsList = new ArrayList<Double>();
        int testNumber = 10;
        for (int i = 0; i < testNumber; i++) {
            tp = 0;
            Set<String> randomPairs = generateRandomPairs(interactionSet);
            for (String pair : randomPairs) {
                if (posChecker.isPositive(pair))
                    tp ++;
            }
            //System.out.printf("Total random: %d%n", randomPairs.size());
            double propOnNeg = (double) tp / randomPairs.size();
            System.out.printf("Mapped to random: %d (%f)%n", 
                              tp,
                              propOnNeg);
            double odds = calculateOddsRatio(propOnPos, propOnNeg);
            oddsList.add(odds);
        }
        DescriptiveStatistics statistics = new DescriptiveStatistics();
        for (Double value : oddsList) {
            if (value.isNaN())
                continue;
            statistics.addValue(value);
        }
        System.out.println("Average odds ratio: " + statistics.getMean() + 
                           " +- " + statistics.getStandardDeviation() +
                           " (from " + testNumber + " tests)");
        System.out.println();
    }
    
    /**
     * Calculate an odds ratio. In this calculation, it is required the numbers of 
     * positive data points and negative data points be the same.
     * @param propOnPos
     * @param propOnNeg
     * @return
     */
    public double calculateOddsRatio(double propOnPos,
                                     double propOnNeg) {
        return propOnPos * (1-propOnNeg) / ((1 - propOnPos) * propOnNeg);
    }
}
