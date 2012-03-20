/*
 * Created on Jul 5, 2006
 *
 */
package org.reactome.weka;

import java.io.File;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

import junit.framework.TestCase;

import org.reactome.fi.util.FileUtility;

import weka.classifiers.Classifier;
import weka.classifiers.bayes.BayesNet;
import weka.classifiers.bayes.NaiveBayes;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.Utils;
import weka.core.converters.ArffLoader;
import weka.core.converters.Loader;


public class WEKAResultAnalyzer extends TestCase {
    public WEKAResultAnalyzer() {
    }
    
    /**
     * Count the union of all predicators for Lincoln (Oct 3, 2006)
     * @param dataset
     */
    protected void count(Instances dataset) {
        Enumeration instancesEnum = dataset.enumerateInstances();
        // Calculate number
        int calTc = 0;
        int calFc = 0;
        int totalT = 0;
        int totalF = 0;
        // Check positive dataset first
        int total = 0;
        while (instancesEnum.hasMoreElements()) {
            Instance inst = (Instance) instancesEnum.nextElement();
            double cls = inst.value(0);
            double hi = inst.value(1);
            double oi = inst.value(2);
            double yi = inst.value(3);
            double ge = inst.value(4);
            if (cls == 0.0d) {
                totalT ++;
                if (hi == 0.0d || oi == 0.0d || yi == 0.0d ||
                    ge == 0.0d || ge == 1.0d)
                    calTc ++;
            }
            else {
                totalF ++;
                if (hi == 0.0d || oi == 0.0d || yi == 0.0d ||
                    ge == 0.0d || ge == 1.0d)
                    calFc ++;
            }
            total ++;
        }
        System.out.println("total: " + total + ", " + " totalT: " + totalT + " totalF: " + totalF);
        System.out.println("Positive: " + totalT + ", predicated: " + calTc + ", " + ((double)calTc / totalT));
        System.out.println("Negative: " + totalF + ", predicated: " + calFc + ", " + ((double)calFc / totalF));
    }
    
    public void checkInteractionsRecall() throws Exception {
        FileUtility fu = new FileUtility();
        Classifier naivebayes = (Classifier) fu.loadObject("results/FourDBNaiveBayes091806.model");
        //String fileName = "results/HumanInteractionData091406.arff";
        //String fileName = "results/HumanPPIsData091806.arff";
        //String fileName = "results/ThreePPIsData091806.arff";
        String fileName = "results/ReactomeData091506_1.arff";
        Loader loader = new ArffLoader();
        loader.setSource(new File(fileName));
        Instances dataset = loader.getDataSet();
        //count(dataset);
        //if (true)
        //    return;
        dataset.setClassIndex(0);
        double[] cutoffs = new double[] {
                0.5, 0.6, 0.7, 0.8, 0.9
        };
        System.out.println("Results for " + fileName + ":");
        int c = 0;
        for (double cutoff : cutoffs) {
            Enumeration instancesEnum = dataset.enumerateInstances();
            // Calculate number
            int calTc = 0;
            int calFc = 0;
            int functionInteraction;
            int humanPPI;
            // Check positive dataset first
            double[] props;
            while (instancesEnum.hasMoreElements()) {
                Instance inst = (Instance) instancesEnum.nextElement();
                props = naivebayes.distributionForInstance(inst);
                if (props[0] > cutoff)
                    calTc ++;
                else {
                    calFc ++;
//                    if (c < 100) {
//                        System.out.println(inst + ": " + props[0]);
//                        c ++;
//                    }
                }
            }
            System.out.println("\tcutoff = " + cutoff);
            int total = dataset.numInstances();
            System.out.printf("\tPositive: Original %d; Predicated True %d, False %d%n",
                              total, calTc, calFc);
            float TP = (float) calTc / total;
            System.out.println("\tTotal Instances: " + total);
            System.out.println("\tTP: " + TP);
            System.out.println();
        }
    }
    
    public void checkTestDataForNaiveBayesClassifier() throws Exception {
        FileUtility fu = new FileUtility();
        //Classifier naivebayes = (Classifier) fu.loadObject("results/NaiveBayes072006_2.model");
        Classifier naivebayes = (Classifier) fu.loadObject("results/NaiveBayes091806_1.model");
        //Classifier naivebayes = (Classifier) fu.loadObject("results/NaiveBayes091106_1.model");
        String[] fileNames = new String[] {
                //"results/ReactomeData091806.arff",
                "results/PantherDataNoReactome091806_1.arff",
                "results/CellMapDataNoReactome091806_1.arff",
                "results/INOHDataNoReactome091806_1.arff",
                "results/NoLocalizationDataNoReactome091806_1.arff",
//                "results/PantherDataMixed090806.arff",
//                "results/ReactomeData090606.arff",
//                "results/BIND090806.arff",
//                "results/HPRD2H091106.arff",
//                "results/BINDIntActData091106.arff",
//                "results/HumanInteractionData091406.arff"
        };
        for (String fileName : fileNames) {
            Loader loader = new ArffLoader();
            loader.setSource(new File(fileName));
            Instances dataset = loader.getDataSet();
            dataset.setClassIndex(0);
            Enumeration instancesEnum = dataset.enumerateInstances();
            // Total number
            int c = 0;
            // Calculate number
            int calTc = 0;
            int calFc = 0;
            // Divide the data into two set
            List<Instance> pos = new ArrayList<Instance>();
            List<Instance> neg = new ArrayList<Instance>();
            int functionInteraction;
            int humanPPI;
            while (instancesEnum.hasMoreElements()) {
                Instance instance = (Instance) instancesEnum.nextElement();
                //humanPPI = (int) instance.value(3);
                //if (humanPPI == 1)
                //    continue; // Only humanPPI true is used
                functionInteraction = (int)instance.value(0);
                if (functionInteraction == 0)
                    pos.add(instance);
                else
                    neg.add(instance);
                c ++;
            }
            // Check positive dataset first
            double[] props;
            for (Instance inst : pos) {
                props = naivebayes.distributionForInstance(inst);
                if (props[0] > 0.5)
                    calTc ++;
                else
                    calFc ++;
            }
            System.out.println("Results for " + fileName + ":");
            System.out.printf("Positive: Original %d; Predicated True %d, False %d%n",
                              pos.size(), calTc, calFc);
            float TP = (float) calTc / pos.size();
            calTc = calFc = 0;
            for (Instance inst : neg) {
                props = naivebayes.distributionForInstance(inst);
                if (props[0] > 0.5)
                    calTc ++;
                else
                    calFc ++;
            }
            System.out.printf("Negative: Original %d; Predicated True %d, False %d%n",
                              neg.size(), calTc, calFc);
            System.out.println("Total Instances: " + c);
            float FP = (float) calTc / neg.size();
            System.out.println("TP: " + TP);
            System.out.println("FP: " + FP);
            System.out.println("Precision: " + TP / (TP + FP));
            System.out.println();
        }
    }
    
    public void checkTestDataForBayesNet() throws Exception {
        FileUtility fu = new FileUtility();
        BayesNet bayesNet = (BayesNet) fu.loadObject("results/BayesNet.mod");
        Instances instances = bayesNet.m_Instances;
        System.out.println("total number of attributes: " + instances.numAttributes());
        // Create a new dataset
        FastVector attributes = new FastVector();
        for (int i = 0; i < instances.numAttributes(); i++) {
            if (i < 7)
                attributes.addElement(instances.attribute(i).copy());
            else 
                attributes.addElement(new Attribute(instances.attribute(i).name()));
        }
        Instances testSet = new Instances("test", attributes, 100);
        testSet.setClassIndex(0);
        // Read in instances
        String testFileName = "results/TestingDataWithOrtho.arff";
        String line = null;
        fu.setInput(testFileName);
        String[] tokens = null;
        double[] probs;
        double predict;
        String predictCls;
        int trueCount = 0;
        int falseCount = 0;
        int totalCount = 0;
        while ((line = fu.readLine()) != null) {
            if (line.length() == 0 || line.startsWith("@"))
                continue;
            tokens = line.split(",");
            Instance inst = new Instance(instances.numAttributes());
            //System.out.println(line);
            inst.setDataset(testSet);
            inst.setClassMissing();
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].equals("?"))
                    continue;
                if (i < 7)
                    inst.setValue(i, tokens[i]);
                else {
                    inst.setValue(i, Double.parseDouble(tokens[i]));
                }
            }
            probs = bayesNet.distributionForInstance(inst);
            predict = Utils.maxIndex(probs);
            predictCls = inst.classAttribute().value((int)predict);
            if (predictCls.equals("true"))
                trueCount ++;
            else if (predictCls.equals("false"))
                falseCount ++;
            totalCount ++;
        }
        System.out.printf("Total %d, True %d, False %d%n", totalCount, trueCount, falseCount);
    }
    
    public void analyzeNaiveBayes() throws Exception {
        FileUtility fu = new FileUtility();
        NaiveBayes naiveBayes = (NaiveBayes) fu.loadObject("results/NaiveBayes070306.model");
        System.out.println("Discretizer is used: " + naiveBayes.getUseSupervisedDiscretization());
        // Create a datasets
        Instances dataset = createDataSet();
        Instance instance = new Instance(8);
        instance.setDataset(dataset);
        instance.setClassMissing();
        instance.setValue(1, "true");
        instance.setValue(2, "true");
        instance.setValue(3, "non");
        //instance.setValue()
        double[] props = naiveBayes.distributionForInstance(instance);
        for (double d : props)
            System.out.print(d + " ");
        System.out.println();
    }
    
    public Instances createDataSetForV3() {
        FastVector attributes = new FastVector();
        // Column
        FastVector values = new FastVector();
        values.addElement("true");
        values.addElement("false");
        Attribute att = new Attribute("FunctionalInteraction", values);
        attributes.addElement(att);

        String[] features = new String[] {
                "HumanPPI",// {true, false}
                "ScePPI",// {true, false}
                "CelPPI",// {true, false}
                "DmePPI",// {true, false}
                "PavlidisGeneExp",// {true, false}
                "CarlosGeneExp",// {true, false}
                "GenewaysPPI",// {true, false}
                "pFamDomainInt",// {true, false}
                "GOBPTermShare",// {true, false}
        };
        for (String feature : features) {
            values = new FastVector();
            values.addElement(Boolean.TRUE.toString());
            values.addElement(Boolean.FALSE.toString());
            att = new Attribute(feature, values);
            attributes.addElement(att);
        }
        
        Instances rtn = new Instances("test", attributes, 1);
        rtn.setClassIndex(0);
        return rtn;
    }
    
    public Instances createDataSet() {
        FastVector attributes = new FastVector();
        FastVector values = new FastVector();
        values.addElement("true");
        values.addElement("false");
        Attribute att = new Attribute("FunctionalInteraction", values);
        attributes.addElement(att);
        values = new FastVector();
        values.addElement("true");
        values.addElement("false");
        att = new Attribute("HumanPPI", values);
        attributes.addElement(att);
        values = new FastVector();
        values.addElement("true");
        values.addElement("false");
        att = new Attribute("OrthoPPI", values);
        attributes.addElement(att);
        values = new FastVector();
        values.addElement("true");
        values.addElement("false");
        att = new Attribute("YeastPPI", values);
        attributes.addElement(att);
        values = new FastVector();
        values.addElement("pos");
        values.addElement("neg");
        values.addElement("non");
        att = new Attribute("GeneExp", values);
        attributes.addElement(att);
        att = new Attribute("GOBPSemSimilarity");
        attributes.addElement(att);
//        att = new Attribute("GOMFSemSimilarity");
//        attributes.addElement(att);
//        att = new Attribute("GOBPOccurence");
//        attributes.addElement(att);
//        att = new Attribute("GOBPDepth");
//        attributes.addElement(att);
//        att = new Attribute("GOMFOccurence");
//        attributes.addElement(att);
//        att = new Attribute("GOMFDepth");
//        attributes.addElement(att);
        Instances rtn = new Instances("test", attributes, 1);
        rtn.setClassIndex(0);
        return rtn;
    }
    
}
