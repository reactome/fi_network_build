/*
 * Created on Apr 3, 2009
 *
 */
package org.reactome.weka;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.reactome.fi.GODataAnalyzer;
import org.reactome.fi.GODataAnalyzer.GOGraphNode;
import org.reactome.fi.PfamAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.Value;

/**
 * This class is used to generate WEKA arff file format.
 * @author wgm
 *
 */
public class ARFFGenerator {
    private FileUtility fu = new FileUtility();
    // To switch to new version
    private boolean isForVersion3 = false;
    
    public ARFFGenerator() {
    }
    
    public void setIsForVersion3(boolean isFor3) {
        this.isForVersion3 = isFor3;
    }
    
    private void generateDataSetForV3(Map<String, Value> values) throws Exception {
        FeatureHandlerForV3 handler = new FeatureHandlerForV3();
        handler.generateDataSetForV3(values);
    }
    
    public void generateDataSet(Map<String, Value> values) throws Exception {
        if (isForVersion3) {
            // Check if the new version of this method should be used.
            generateDataSetForV3(values);
            return;
        }
        Value value = null;
        // For PPI datas
        //Set<String> humanInteractions = fu.loadInteractions(RESULT_DIR + "HumanInteractions.txt");
        //Set<String> humanInteractions = fu.loadInteractions(RESULT_DIR + "HumanInteractions090806.txt");
        //Set<String> humanInteractions = fu.loadInteractions(RESULT_DIR + "HPRDInteractions090806.txt");
        //Set<String> humanInteractions = fu.loadInteractions(RESULT_DIR + "HumanInteractions120506.txt");
        //Set<String> humanInteractions = fu.loadInteractions("results/v2/HumanInteractions020507.txt");
        Set<String> humanInteractions = fu.loadInteractions(FIConfiguration.getConfiguration().get("HUMAN_INTERACTION_FILE_NAME"));
        // For Ortho Interaction data
        //Set<String> orthoInteractions = fu.loadInteractions(RESULT_DIR + "OrthoInteractions.txt");
        Set<String> orthoInteractions = fu.loadInteractions(FIConfiguration.getConfiguration().get("ORTHO_INTERAACTION_FILE_NAME"));
        // Datasets mapped from yeast interactions
        //Set<String> yeastInteractions = fu.loadInteractions(RESULT_DIR + "YeastInteractions.txt");
        Set<String> yeastInteractions = fu.loadInteractions(FIConfiguration.getConfiguration().get("YEAST_INTERACTION_FILE_NAME"));
        // For geneways interactions
//        Set<String> unambGenewaysPPI = fu.loadInteractions(RESULT_DIR+"Copy of genewaysInteractions_ambigous.txt");
//        Set<String> ambGenewaysPPI = fu.loadInteractions(RESULT_DIR+"Copy of genewaysInteractions_unambigous_human.txt");
//        Set<String> unambGenewaysFlybaseDmePPI = fu.loadInteractions(RESULT_DIR+"unambigousPairs_flybase_flyHsaUni.txt");
        // Please refer to Value.java for the meanings of these numbers
//        Set<String> f347274 = fu.loadInteractions(RESULT_DIR+"orthomcl/ortho_347274.txt");
//        Set<String> f364033 = fu.loadInteractions(RESULT_DIR+"orthomcl/ortho_364033.txt");
//        Set<String> f342698= fu.loadInteractions(RESULT_DIR+"orthomcl/ortho_342698.txt");
//        Set<String> f239492 = fu.loadInteractions(RESULT_DIR+"orthomcl/ortho_239492.txt");
       /* Set<String> f361401= fu.loadInteractions(RESULT_DIR+"orthomcl/ortho_361401.txt");
        Set<String> f192607 = fu.loadInteractions(RESULT_DIR+"swissprot_192607.txt");
        Set<String> f216051 = fu.loadInteractions(RESULT_DIR+"swissprot_216051.txt");
        Set<String> f310914 = fu.loadInteractions(RESULT_DIR+"swissprot_310914.txt");*/
//        int debugCounter = 0;
        
        // New patterns in order to incorporate
       // Set<String> debugSet = new HashSet<String>();
        // To control the memory footprint, GeneExp and GO datasets are not loaded. See below.
        // for functions
        long time1 = System.currentTimeMillis();
        String key = null;
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            key = it.next();
            value = values.get(key);
            if (humanInteractions.contains(key))
                value.humanInteraction = Boolean.TRUE;
            if (orthoInteractions.contains(key))
                value.orthoInteraction = Boolean.TRUE;
            if (yeastInteractions.contains(key))
                value.yeastInteraction = Boolean.TRUE;
//            if (unambGenewaysPPI.contains(key))
//              value.unambGenewaysPPI = Boolean.TRUE;
//            if (ambGenewaysPPI.contains(key)) {
//              value.ambGenewaysPPI = Boolean.TRUE;
//            }
//            if (unambGenewaysFlybaseDmePPI.contains(key)) {
//              value.unambGenewaysFlybaseDmePPI = Boolean.TRUE;
//            }
//            if(f347274.contains(key)) {
//              value.f347274 = Boolean.TRUE;
//              debugCounter = 347274;
//            }
//            if(f364033.contains(key)) {
//              value.f364033 = Boolean.TRUE;
//              debugCounter = 364033;
//            }
//            if(f342698.contains(key)) {
//              value.f342698 = Boolean.TRUE;
//              debugCounter = 361401;
//            }
//            if(f239492.contains(key)) {
//              value.f239492 = Boolean.TRUE;
//              debugCounter = 363515;
//            }
          /*  if(f310914.contains(key)) {
                value.f310914 = Boolean.TRUE;
                debugCounter = 310914;
            }*/
           /* if(f216051.contains(key)) {
                value.f216051 = Boolean.TRUE;
                debugCounter = 216051;
            }
            if(f192607.contains(key)) {
                value.f192607 = Boolean.TRUE;debugCounter = 192607;
            }
            if(f361401.contains(key)) {
                value.f361401= Boolean.TRUE;debugCounter = 188996;
            }*/
            values.put(key, value);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for PPIs: " + (time2 - time1));
        // Handle gene expression: check loaded gene expression data
        generateGeneExpData(values);
        long time3 = System.currentTimeMillis();
        System.out.println("Time for GeneExp: " + (time3 - time2));
        //generateGOData();
        generateGOSemanticData(values);
        long time4 = System.currentTimeMillis();
        System.out.println("Time for Generate GO: " + (time4 - time3));
        // Check Pfam
        //generatePfamData();
    }
    
    private void generateGOSemanticData(Map<String, Value> values) throws IOException {
        // Now work for GO datasets
        // GO Occurence
        GODataAnalyzer goAnalyzer = new GODataAnalyzer();
        Map<String, Set<String>> prot2BPMap = new HashMap<String, Set<String>>();
        Map<String, Set<String>> prot2MFMap = new HashMap<String, Set<String>>();
        goAnalyzer.generateProtein2GOMap(prot2MFMap, prot2BPMap);
        Map<String, GOGraphNode> id2Node = goAnalyzer.generateGOTree();
        //Map<String, Double> mfId2InfoContent = goAnalyzer.loadInformationContents("results/go/GOMFInfoContent.txt");
        // Map<String, Double> bpId2InfoContent = goAnalyzer.loadInformationContents("results/go/GOBPInfoContent.txt");
        Map<String, Double> bpId2InfoContent = goAnalyzer.loadBPInformationContents();
        String pair;
        int index;
        Value valueObj;
        Double bpSemSim;
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            pair = it.next();
            valueObj = values.get(pair);
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            //valueObj.goMFSemSimilarity = goAnalyzer.calculateSemanticSimilarity(prot2MFMap.get(id1),
            //                                                                    prot2MFMap.get(id2),
            //                                                                    mfId2InfoContent,
            //                                                                    id2Node);
            
            bpSemSim = goAnalyzer.calculateSemanticSimilarity(prot2BPMap.get(id1),
                                                              prot2BPMap.get(id2),
                                                              bpId2InfoContent,
                                                              id2Node);
            if (bpSemSim != null)
                valueObj.goBPSemSimilarity = discretizeGOBP(bpSemSim);
        }
    }
    
    
    private int discretizeGOBP(double value) {
//      '(-inf-0.474079]'   79147
//      '(0.474079-0.948157]'   92448
//      '(0.948157-1.422236]'   73472
//      '(1.422236-1.896315]'   39246
//      '(1.896315-2.370393]'   24713
//      '(2.370393-2.844472]'   11451
//      '(2.844472-3.318551]'   7150
//      '(3.318551-3.792629]'   2714
//      '(3.792629-4.266708]'   7366
//      '(4.266708-4.740786]'   3173
//      '(4.740786-5.214865]'   2390
//      '(5.214865-5.688944]'   1014
//      '(5.688944-6.163022]'   1595
//      '(6.163022-6.637101]'   99
//      '(6.637101-7.11118]'    51
//      '(7.11118-7.585258]'    120
//      '(7.585258-8.059337]'   69
//      '(8.059337-8.533416]'   44
//      '(8.533416-9.007494]'   0
//      '(9.007494-inf)'    4
        double[] breaks = new double[]{
                0.474079, 
                0.948157, 
                1.422236,
                1.896315,
                2.370393,
                2.844472,
                3.318551,
                3.792629,
                4.266708,
                4.740786,
                5.214865,
                5.688944,
                6.163022,
                6.637101,
                7.11118,
                7.585258,
                8.059337,
                8.533416,
                9.007494
        };
        for (int i = breaks.length; i > 0; i--) {
            if (value > breaks[i - 1])
                return i;
        }
        return 0;
    }
    
    private int discretizeGOMP(double value) {
//      '(-inf-1.037472]'   52197
//      '(1.037472-2.074945]'   15006
//      '(2.074945-3.112417]'   4293
//      '(3.112417-4.14989]'    1061
//      '(4.14989-5.187362]'    1022
//      '(5.187362-6.224835]'   150
//      '(6.224835-7.262307]'   29
//      '(7.262307-8.29978]'    64
//      '(8.29978-9.337252]'    1
//      '(9.337252-inf)'    1
      if (value > 9.337252)
          return 9;
      if (value > 8.29978)
          return 8;
      if (value > 7.262307)
          return 7;
      if (value > 6.224835)
          return 6;
      if (value > 5.187362)
          return 5;
      if (value > 4.14989)
          return 4;
      if (value > 3.112417)
          return 3;
      if (value > 2.074945)
          return 2;
      if (value > 1.037472)
          return 1;
      return 0;
    }
    
    private void generatePfamData(Map<String, Value> values) throws IOException {
        PfamAnalyzer pfamAnalyzer = new PfamAnalyzer();
        String pair = null;
        Value valueObj;
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            pair = it.next();
            valueObj = values.get(pair);
            if (pfamAnalyzer.checkIfInteracting(pair))
                valueObj.pfamDomainInt = Boolean.TRUE;
        }
    }
    
    private void generateGeneExpData(Map<String, Value> values) throws IOException {
        Value value;
        //String geneExpDataName = "results/microarray/GeneExpFromPavlidis.txt";
        //String geneExpDataName = "results/microarray/GeneExpWith3FromPavlidis.txt";
        String geneExpDataName = FIConfiguration.getConfiguration().get("GENE_EXP_FILE_NAME");
        // Format: "+/-" + "\t" + id1 + PfamAnalyzer.java" " + id2 (id1 and id2 are sorted)
        fu.setInput(geneExpDataName);
        String line = null;
        int index = 0;
        String pair = null;
        String geneExpValue = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            pair = line.substring(index + 1);
            value = values.get(pair);
            if (value != null) {
                geneExpValue = line.substring(0, index);
                if (geneExpValue.equals("+"))
                    value.geneExp = "pos";
                else if (geneExpValue.equals("-"))
                    value.geneExp = "neg";
            }
        }
        fu.close();
    }
    
    private void generateGOOccurenceData(boolean isForBP,
                                         Map<String, Set<String>> prot2GOMap,
                                         GODataAnalyzer analyzer,
                                         Map<String, Value> values) throws IOException {
        String fileName;
        if (isForBP) 
            fileName = "results/go/GOPOccurence.txt";
        else 
            fileName = "results/go/GOFOccurence.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        Map<String, Integer> occurenceMap = new HashMap<String, Integer>();
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            occurenceMap.put(line.substring(0, index), 
                             Integer.parseInt(line.substring(index + 1)));
        }
        fu.close();
        Value valueObj;
        String pair;
        int occurence;
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            pair = it.next();
            valueObj = values.get(pair);
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            occurence = analyzer.findSharedTermOccurence(prot2GOMap.get(id1),
                                                         prot2GOMap.get(id2),
                                                         occurenceMap);
            if (occurence == Integer.MAX_VALUE)
                continue;
            if (isForBP)
                valueObj.goBPOccurence = occurence;
            else
                valueObj.goMFOccurence = occurence;
        }
    }
    
    private void generateGOData(Map<String, Value> values) throws IOException {
        // Now work for GO datasets
        // GO Occurence
        GODataAnalyzer goAnalyzer = new GODataAnalyzer();
        Map<String, Set<String>> prot2BPMap = new HashMap<String, Set<String>>();
        Map<String, Set<String>> prot2MFMap = new HashMap<String, Set<String>>();
        goAnalyzer.generateProtein2GOMap(prot2MFMap, prot2BPMap);
        generateGOOccurenceData(true, prot2BPMap, goAnalyzer, values);
        generateGOOccurenceData(false, prot2MFMap, goAnalyzer, values);
        // GO Depth
        Map<String, Integer> depthMap = new HashMap<String, Integer>();
        // Load depth values
        FileUtility fu = new FileUtility();
        fu.setInput("results/go/GODepth.txt");
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            depthMap.put(line.substring(0, index), Integer.parseInt(line.substring(index + 1)));
        }
        fu.close();
        generateGODepthData(true, depthMap, prot2BPMap, goAnalyzer, values);
        generateGODepthData(false, depthMap, prot2MFMap, goAnalyzer, values);
    }
    
    private void generateGODepthData(boolean isForBP,
                                     Map<String, Integer> depthMap,
                                     Map<String, Set<String>> prot2GOMap,
                                     GODataAnalyzer analyzer,
                                     Map<String, Value> values) throws IOException {
        int index = 0;
        Value valueObj;
        String pair;
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            pair = it.next();
            valueObj = values.get(pair);
            index = pair.indexOf(" ");
            String id1 = pair.substring(0, index);
            String id2 = pair.substring(index + 1);
            int depth = analyzer.findSharedDepth(prot2GOMap.get(id1),
                    prot2GOMap.get(id2),
                    depthMap);
            if (depth == Integer.MIN_VALUE)
                continue;
            if (isForBP)
                valueObj.goBPDepth = depth;
            else
                valueObj.goMFDepth = depth;
        }
    }    
    
    private boolean isValidValue(Value value) {
        if (value.intactHumanPPI != null ||
            value.hprdHumanPPI != null ||
            value.biogridHumanPPI != null ||
            value.scePPI != null ||
            value.celPPI != null ||
            value.dmePPI != null ||
            value.pavlidisGeneExp != null ||
            value.carlosGeneExp != null ||
            value.genewaysPPI != null ||
            value.pfamDomainInt != null ||
            value.goBPSharing != null)
            return true;
        return false;
    }
    
    /**
     * New version to do exporting.
     * @param fileName
     * @param values
     * @throws IOException
     */
    private void exportDataForWEKAForV3(String fileName,
                                        Map<String, Value> values) throws IOException {
        StringBuilder builder = new StringBuilder();
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        fu.printLine("@relation FunctionTraining");
        fu.printLine("");
        // Have to keep the following sequence
        fu.printLine("@attribute FunctionalInteraction {true, false}");
        //fu.printLine("@attribute IntActHumanPPI {true, false}");
        //fu.printLine("@attribute HPRDPPI {true, false}");
        //fu.printLine("@attribute BioGridPPI {true, false}");
        fu.printLine("@attribute HumanPPI {true, false}");
        fu.printLine("@attribute ScePPI {true, false}");
        fu.printLine("@attribute CelPPI {true, false}");
        fu.printLine("@attribute DmePPI {true, false}");
        fu.printLine("@attribute PavlidisGeneExp {true, false}");
        fu.printLine("@attribute CarlosGeneExp {true, false}");
        fu.printLine("@attribute GenewaysPPI {true, false}");
        fu.printLine("@attribute pFamDomainInt {true, false}");
        fu.printLine("@attribute GOBPTermShare {true, false}");
        fu.printLine("");
        fu.printLine("@data");
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            String pair = it.next();
            Value value = values.get(pair);
            // Exclude nonsense data points
            //if (!isValidValue(value))
            //    continue;
            builder.append(value.functionalInteraction).append(",");
            //exportFeatureInValue(builder, value.intactHumanPPI);
            //exportFeatureInValue(builder, value.hprdHumanPPI);
            //exportFeatureInValue(builder, value.biogridHumanPPI);
            exportFeatureInValue(builder, value.humanInteraction);
            exportFeatureInValue(builder, value.scePPI);
            exportFeatureInValue(builder, value.celPPI);
            exportFeatureInValue(builder, value.dmePPI);
            exportFeatureInValue(builder, value.pavlidisGeneExp);
            exportFeatureInValue(builder, value.carlosGeneExp);
            exportFeatureInValue(builder, value.genewaysPPI);
            exportFeatureInValue(builder, value.pfamDomainInt);
            exportFeatureInValue(builder, value.goBPSharing);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }

    private void exportFeatureInValue(StringBuilder builder, Boolean feature) {
        if (feature == null)
            builder.append("false,");
        else
            builder.append(feature).append(",");
    }
    
    public void exportDataForWEKA(String fileName, 
                                   Map<String, Value> values) throws IOException {
        if (isForVersion3) {
            exportDataForWEKAForV3(fileName, values);
            return;
        }
        StringBuilder builder = new StringBuilder();
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        fu.printLine("@relation FunctionTraining");
        fu.printLine("");
        fu.printLine("@attribute FunctionalInteraction {true, false}");
        fu.printLine("@attribute HumanPPI {true, false}");
        fu.printLine("@attribute OrthoPPI {true, false}");
        fu.printLine("@attribute YeastPPI {true, false}");
        //fu.printLine("@attribute PfamDomainInt {true, false}");
        //MicroarrayDataAnalyzer analyzer = new MicroarrayDataAnalyzer();
        //List<String> gdsList = analyzer.getGDSList();
        //for (String gds : gdsList) {
        //    fu.printLine("@attribute " + gds + " real");
        //}
        //fu.printLine("@attribute GeneExp real");
//        fu.printLine("@attribute unambGenewaysPPI {true, false}");
//        fu.printLine("@attribute ambGenewaysPPI {true, false}");
//        fu.printLine("@attribute unambGenewaysFlybaseDmePPI {true, false}");
//        fu.printLine("@attribute f347274 {true, false}");
//        fu.printLine("@attribute f364033 {true, false}");
//        fu.printLine("@attribute f342698 {true, false}");
//        fu.printLine("@attribute f239492 {true, false}");
       /* fu.printLine("@attribute f361401{true, false}");
        fu.printLine("@attribute f192607 {true, false}");
        fu.printLine("@attribute f216051 {true, false}");
        fu.printLine("@attribute f310914 {true, false}");*/
        fu.printLine("@attribute GeneExp {pos, neg, non}");
        //fu.printLine("@attribute GOBPOccurence real");
        //fu.printLine("@attribute GOBPDepth real");
        //fu.printLine("@attribute GOMFOccurence real");
        //fu.printLine("@attribute GOMFDepth real");
        //fu.printLine("@attribute GOBPSemSimilarity real");
        fu.printLine("@attribute GOBPSemSimilarity {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}");
        //fu.printLine("@attribute GOMFSemSimilarity real");
        fu.printLine("");
        fu.printLine("@data");
        for (Iterator<String> it = values.keySet().iterator(); it.hasNext();) {
            String pair = it.next();
            Value value = values.get(pair);
            // Exclude nosense data points
            if (value.humanInteraction == null &&
                value.orthoInteraction == null &&
                value.yeastInteraction == null &&
                //value.pfamDomainInt == null &&
                value.geneExp == null &&
                value.goBPSemSimilarity == null // &&
//                value.ambGenewaysPPI == null &&
//                value.unambGenewaysPPI == null &&
//                value.unambGenewaysFlybaseDmePPI == null &&
//                value.f347274 == null &&
//                value.f364033 == null &&
//                value.f342698 == null &&
//                value.f239492 == null 
            /*  value.f361401== null */
      /*          value.f192607 == null &&
                value.f216051 == null &&
                value.f310914 == null*/)// &&
            
                //value.goMFSemSimilarity == null)
                continue;
            builder.append(value.functionalInteraction).append(",");
            if (value.humanInteraction == null)
                builder.append("false,");
            else
                builder.append(value.humanInteraction).append(",");
            if (value.orthoInteraction == null)
                builder.append("false,");
            else
                builder.append(value.orthoInteraction).append(",");
            if (value.yeastInteraction == null)
                builder.append("false,");
            else
                builder.append(value.yeastInteraction).append(",");
//            if (value.unambGenewaysPPI == null)
//                builder.append("false,");
//            else
//                builder.append(value.unambGenewaysPPI).append(",");
//            if (value.ambGenewaysPPI == null)
//                builder.append("false,");
//            else
//                builder.append(value.ambGenewaysPPI).append(",");
//            if (value.unambGenewaysFlybaseDmePPI == null)
//                builder.append("false,");
//            else
//                builder.append(value.unambGenewaysFlybaseDmePPI).append(",");
//            if (value.f347274 == null)
//              builder.append("false,");
//            else
//              builder.append(value.f347274).append(",");
//            if (value.f364033 == null)
//              builder.append("false,");
//            else
//              builder.append(value.f364033).append(",");
//            if (value.f342698 == null)
     //         builder.append("false,");
       //     else
         //     builder.append(value.f342698).append(",");
          //  if (value.f239492 == null)
           //   builder.append("false,");
          //  else
          //    builder.append(value.f239492).append(",");
          /*  if (value.f361401== null)
                builder.append("false,");
            else
                builder.append(value.f361401).append(",");*/
        /*    if (value.f192607 == null)
                builder.append("false,");
            else
                builder.append(value.f192607).append(",");
            if (value.f216051 == null)
                builder.append("false,");
            else
                builder.append(value.f216051).append(",");
            if (value.f310914 == null)
                builder.append("false,");
            else
                builder.append(value.f310914).append(",");*/
//            if (value.pfamDomainInt == null)
//                builder.append("false,");
//            else
//                builder.append(value.pfamDomainInt).append(",");
//          List<Float> geneExpValues = value.geneExpList;
//          for (Float f : geneExpValues) {
//          if (f == null)
//          builder.append("?,");
//          else
//          builder.append(f).append(",");
//          }
            if (value.geneExp == null)
                builder.append("non,");
            else
                builder.append(value.geneExp).append(",");
            if (value.goBPSemSimilarity == null)
                builder.append("?,");
            else
                builder.append(value.goBPSemSimilarity).append(",");
//            if (value.goMFSemSimilarity == null)
//                builder.append("?");
//            else
//                //builder.append(value.goMFSemSimilarity);
//                builder.append(discretizeGOMP(value.goMFSemSimilarity));
//            if (value.goBPOccurence == null)
//                builder.append("?,");
//            else
//                builder.append(value.goBPOccurence).append(","); 
//            if (value.goBPDepth == null)
//                builder.append("?,");
//            else
//                builder.append(value.goBPDepth).append(",");
//            if (value.goMFOccurence == null)
//                builder.append("?,");
//            else
//                builder.append(value.goMFOccurence).append(",");
//            if (value.goMFDepth == null)
//                builder.append("?");
//            else
//                builder.append(value.goMFDepth);
            fu.printLine(builder.toString());
            builder.setLength(0);
        }
        fu.close();
    }
}
