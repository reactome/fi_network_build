/*
 * Created on Apr 6, 2007
 *
 */
package org.reactome.fi.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.HypergeometricDistribution;
import org.apache.commons.math.distribution.HypergeometricDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math.stat.inference.TestUtils;
import org.junit.Test;

import cern.jet.math.Arithmetic;
import cern.jet.random.Binomial;
import cern.jet.random.ChiSquare;
import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;

public class MathUtilities {
    
    //private static DistributionFactory distFactory;
    private static RandomEngine randomeEngine;
    
    static {
        //distFactory = DistributionFactory.newInstance();
        randomeEngine = new DRand();
    }
    
    /**
     * It seems that there are some bugs in the package for LogRankTest(). The results
     * from this implementation is differnet from R.
     * @param time1
     * @param censor1
     * @param time2
     * @param censor2
     * @return
     */
//    public static double logRankSurvivalTest(List<Double> time1,
//                                             List<Double> censor1,
//                                             List<Double> time2,
//                                             List<Double> censor2) {
//        // Convert an double list to a single array
//        double[] time11 = convertDoubleListToArray(time1);
//        double[] censor11 = convertDoubleListToArray(censor1);
//        double[] time21 = convertDoubleListToArray(time2);
//        double[] censor21 = convertDoubleListToArray(censor2);
//        LogRankTest test = new LogRankTest(time11, censor11, time21, censor21);
//        return test.pValue;
//    }
    
    private static double[] convertDoubleListToArray(List<Double> list) {
        double[] rtn = new double[list.size()];
        for (int i = 0; i < list.size(); i++)
            rtn[i] = list.get(i);
        return rtn;
    }
    
    public static double calculateDistance(List<Integer> v1,
                                           List<Integer> v2) {
        if (v1.size() != v2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        double t = 0.0d;
        for (int i = 0; i < v1.size(); i++) {
            int i1 = v1.get(i);
            int i2 = v2.get(i);
            t += (double) (i1 - i2) * (i1 - i2);
        }
        return Math.sqrt(t);
    }
    
    public static double log2(double value) {
        return Math.log(value) / Math.log(Math.E);
    }
    
    public static double calculateHammingDistance(List<Boolean> vector1,
                                                  List<Boolean> vector2) {
        if (vector1.size() != vector2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        int dist = 0;
        for (int i = 0; i < vector1.size(); i++) {
            Boolean b1 = vector1.get(i);
            Boolean b2 = vector2.get(i);
            if (!b1.equals(b2))
                dist ++;
        }
        return dist;
    }
    
    public static double calculateNetworkDistance(List<Boolean> vector1,
                                                  List<Boolean> vector2,
                                                  List<Set<String>> clusters) {
        if (vector1.size() != vector2.size())
            throw new IllegalArgumentException("The passed two boolean vectors should have the same length!");
        int similarity = 0;
        int length = vector1.size();
        int total = 0;
        for (int i = 0; i < vector1.size(); i++) {
            Set<String> cluster = clusters.get(i);
            int weight = cluster.size();
            total += weight;
            Boolean b1 = vector1.get(i);
            Boolean b2 = vector2.get(i);
            if (b1.equals(b2))
                similarity += weight;
        }
        return total - similarity;
    }
    
    public static double calculateTTest(List<Double> sample1,
                                        List<Double> sample2) throws MathException {
        double[] sampleArray1 = new double[sample1.size()];
        for (int i = 0; i < sample1.size(); i++)
            sampleArray1[i] = sample1.get(i);
        double[] sampleArray2 = new double[sample2.size()];
        for (int i = 0; i < sample2.size(); i++)
            sampleArray2[i] = sample2.get(i);
        return TestUtils.tTest(sampleArray1, sampleArray2);
    }
    
    public static double calculateMean(List<Double> values) {
        double total = 0.0d;
        for (Double value : values)
            total += value;
        return total / values.size();
    }
    
    public static double calculateJaccardIndex(Collection<String> set1,
                                               Collection<String> set2) {
        Set<String> copy1 = new HashSet<String>(set1);
        Set<String> copy2 = new HashSet<String>(set2);
        Set<String> shared = new HashSet<String>(copy1);
        shared.retainAll(copy2);
        copy1.addAll(copy2);
        return (double) shared.size() / copy1.size();
    }
    
    @Test
    public void test() throws Exception {
//        double p = calculateHypergeometricPValue(45, 
//                                                 17, 
//                                                 15, 
//                                                 9);
//        System.out.println("P value: " + p);
        double zvalue = calculateZValue(76, 91, 14, 18);
        double pvalue = calTwoTailStandardNormalPvalue(zvalue);
        System.out.println(pvalue);
    }
    
    public static double calOneTailedBinomPValue(double ratio,
                                         int sampleSize,
                                         int success) {
        // Try to use apache math lib
//        BinomialDistribution binomial = distFactory.createBinomialDistribution(sampleSize, ratio);
//        try {
//            return 1.0d - binomial.cumulativeProbability(success - 1);
//        }
//        catch(MathException e) {
//            e.printStackTrace();
//        }
//        return Double.MAX_VALUE;
        // Try to use cern lib
        //One tailed only!!!
        Binomial cernBinomial = new Binomial(sampleSize, ratio, randomeEngine);
        return 1.0d - cernBinomial.cdf(success - 1);
    }
    
    public static double calculatePValue(double ratio,
                                         int sampleSize,
                                         int success) {
        // Try to use apache math lib
        //BinomialDistribution binomial = distFactory.createBinomialDistribution(sampleSize, ratio);
        //try {
        //return 1.0d - binomial.cumulativeProbability(success - 1);
        //}
        //catch(MathException e) {
        //e.printStackTrace();
        //}
        //return Double.MAX_VALUE;
        // Try to use cern lib
        //One tailed only!!!
        Binomial cernBinomial = new Binomial(sampleSize, ratio, randomeEngine);
        if(success == 0) // To avoid unreasonal value
            success = 1;
        return 1.0d - cernBinomial.cdf(success - 1);
    }

    public static double calOneTailedStandardNormalPvalue(double z){
        Normal stdNormal = new Normal(0,1,randomeEngine);
        return stdNormal.cdf(z);
    }
    
    /**
     * Calculate a up-tailed p-value based on hypergeometic distribution.
     * @param N the total balls in urn
     * @param s the white balls (as success)
     * @param n the sample size
     * @param success the white balls picked up in the sample size
     * @return
     */
    public static double calculateHypergeometricPValue(int N,
                                                       int s,
                                                       int n,
                                                       int success) throws MathException {
        HypergeometricDistributionImpl hyper = new HypergeometricDistributionImpl(N,
                                                                                  s,
                                                                                  n);
        return hyper.upperCumulativeProbability(success);
//        int max = Math.min(s, n);
//        if (success >= max / 2)
//            return hyper.upperCumulativeProbability(success);
//        else
//            return 1.0d - hyper.cumulativeProbability(success - 1);
//        return hyper.cumulativeProbability(success, max);
    }
    
    @Test
    public void testHyper() throws MathException {
        HypergeometricDistribution hyper = new HypergeometricDistributionImpl(9083, 1977, 201);
        double pvalue1 = hyper.cumulativeProbability(54, 201);
        System.out.println(pvalue1);
        double pvalue2 = hyper.cumulativeProbability(53);
        System.out.println(pvalue2);
        System.out.println(pvalue1 + pvalue2);
    }
    
    public Normal normalPDF(){
        return new Normal(0,1,randomeEngine);
    }
    
    public static double calTwoTailStandardNormalPvalue(double z){
        z = Math.abs(z);
        return 2*(1-calOneTailedStandardNormalPvalue(z));
    }
    
    public static double calculateZValue(int success1, int sample1,
                                         int success2, int sample2) {
        double ratio1 = (double) success1 / sample1;
        double ratio2 = (double) success2 / sample2;
        // Need to calculate z value
        double p = (success1 + success2) / (double) (sample1 + sample2);
        double z = (ratio1 - ratio2) / Math.sqrt(p * (1.0 - p) * (1.0 / sample1 + 1.0 / sample2));
        return z;
    }
    
    public static double twoTailGroupZTest(double p0, double p, int sampleSize){
        double z = (p-p0)*1.0/(Math.sqrt((p0*(1-p0))/sampleSize));
        return  calTwoTailStandardNormalPvalue(z);
    }
    
    public static double chiSquare(double df, double x2){
        ChiSquare cs = new ChiSquare(df,randomeEngine);
        return cs.cdf(x2);
    }
    
    public static double calculateDbinom(int success,
                                         int total,
                                         double ratio) {
        return Arithmetic.binomial(total, success) * 
               Math.pow(ratio, success) *
               Math.pow((1- ratio), total - success);
               
    }
    
    public static double calculateEnrichment(double ratio,
                                             int sampleTotal,
                                             int sampleSuccess) {
        double newValue = (double) sampleSuccess / sampleTotal;
        return newValue / ratio;
    }
    
    public static Set<String> randomSampling(Collection<String> set,
                                             int size) {
        Set<String> rtn = new HashSet<String>();
        int total = set.size();
        List<String> list = null;
        if (set instanceof List)
            list = (List<String>) set;
        else
            list = new ArrayList<String>(set);
        int index;
        while (rtn.size() < size) {
            index = (int) (total * Math.random());
            rtn.add(list.get(index));
        }
        return rtn;
    }
    
    /**
     * Use this static method to contruct a PearsonCorrelation object to get Pearson correaltion
     * value and its p-value.
     * @param values1
     * @param values2
     * @return
     */
    public static PearsonsCorrelation constructPearsonCorrelation(List<Double> values1,
                                                                  List<Double> values2) {
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("Two double lists have different lengths: " + 
                                                values1.size() + " and " + values2.size());
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(values1.size(), 2);
        for (int i = 0; i < values1.size(); i++) {
            matrix.addToEntry(i, 0, values1.get(i));
            matrix.addToEntry(i, 1, values2.get(i));
        }
        PearsonsCorrelation correlation = new PearsonsCorrelation(matrix);
        return correlation;
    }
    
    /**
     * This method is used to calculate Pearson correlation coefficient. The provided two
     * double lists should have the same size.
     * @param values1
     * @param values2
     * @return
     */
    public static Double calculatePearsonCorrelation(List<Double> values1,
                                                     List<Double> values2) {
        // Make sure the two lists have the same size
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("The provided two double lists have different sizes!");
        double av1 = 0.0, av2 = 0.0, var1 = 0.0, var2 = 0.0, var12 = 0.0, c;
        
        int size = values1.size();
        for (int i = 0; i < size; i++) {
          av1 += values1.get(i);
          av2 += values2.get(i);
        }
        av1 /= (double) size;
        av2 /= (double) size;
        for (int i = 0; i < size; i++) {
          var1 += (values1.get(i) - av1) * (values1.get(i) - av1);
          var2 += (values2.get(i) - av2) * (values2.get(i) - av2);
          var12 += (values1.get(i) - av1) * (values2.get(i) - av2);
        }
        if (var1 * var2 == 0.0) {
          c = 1.0;
        }
        else {
          c = var12 / Math.sqrt(Math.abs(var1 * var2));
        }
        
        return c;
    }
    
    public static double calculatePearsonCorrelation1(List<Double> values1,
                                                     List<Double> values2) {
        // Make sure the two lists have the same size
        if (values1.size() != values2.size())
            throw new IllegalArgumentException("The provided two double lists have different sizes!");
        double av1 = 0.0, av2 = 0.0, sqav1 = 0.0, sqav2 = 0.0, cross = 0.0d;
        
        int size = values1.size();
        for (int i = 0; i < size; i++) {
          av1 += values1.get(i);
          sqav1 += values1.get(i) * values1.get(i);
          av2 += values2.get(i);
          sqav2 += values2.get(i) * values2.get(i);
          cross += values1.get(i) * values2.get(i);
        }
        av1 /= size;
        sqav1 /= size;
        av2 /= size;
        sqav2 /= size;
        cross /= size;
        return (cross - av1 * av2) / Math.sqrt((sqav1 - av1 * av1) * (sqav2 - av2 * av2));
    }
    
    /**
     * This method should be used instead of randomSampling(Collection<String>, int) since
     * it should be a more standard randomization implementation.
     * @param set
     * @param size
     * @param randomizer
     * @return
     */
    @SuppressWarnings("unchecked")
    public static <T> Set<T> randomSampling(Collection<T> set,
                                             int size,
                                             RandomData randomizer) {
        Object[] objects = randomizer.nextSample(set, size);
        Set<T> rtn = new HashSet<T>();
        for (Object obj : objects)
            rtn.add((T) obj);
        return rtn;
    }
    
    /**
     * Permutate a list of objects.
     * @param list
     * @param randomizer
     * @return
     */
    public static <T> List<T> permutate(List<T> list,
                                        RandomData randomizer) {
        int[] indices = randomizer.nextPermutation(list.size(), list.size());
        List<T> rtn = new ArrayList<T>();
        for (int i : indices) {
            rtn.add(list.get(i));
        }
        return rtn;
    }
    
    @Test
    public void testRandomization() {
        int total = 500;
        RandomDataImpl randomizer = new RandomDataImpl();
        for (int i = 0; i < 10; i++) {
            StringBuilder builder = new StringBuilder();
            for (int j = 0; j < 100; j++) {
                //int random = (int) (total * Math.random());
                int random = randomizer.nextInt(0, total - 1);
                builder.append(random).append(", ");
            }
            System.out.println(builder.toString());
        }
    }
    
}
