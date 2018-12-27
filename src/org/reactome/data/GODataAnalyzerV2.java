/*
 * Created on Mar 31, 2009
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.MathUtilities;
import org.reactome.fi.util.PositiveChecker;

/**
 * This is a new analyzer for GO related stuff.
 * @author wgm
 *
 */
public class GODataAnalyzerV2 {
    private FileUtility fu = new FileUtility();
    private String GOA_FILE_NAME = FIConfiguration.getConfiguration().get("GOA_FILE_NAME");
    private GOTermLoader termLoader;
    
    public GODataAnalyzerV2() {
        termLoader = new GOTermLoader();
        // As of the 2018 version, we will escape annotations loaded directly from Reactome.
        termLoader.setEscapeReactomeAnnotations(true);
        termLoader.setGoaFileName(GOA_FILE_NAME);
    }
    
    public Map<String, Set<String>> loadProteinToGOBPTerms() throws IOException {
        return termLoader.loadProteinToGOBPTerms();
    }
    
    public Map<String, Set<String>> loadProteinToGOMFTerms() throws IOException {
        return termLoader.loadProteinToGOMFTerms();
    }
    
    public Map<String, Set<String>> loadProteinToGOCCTerms() throws IOException {
        return termLoader.loadProteinToGOCCTerms();
    }

    /**
     * Check if a protein pair share any term
     * @param proteinPair
     * @param proteinToTerms
     * @return
     */
    public boolean isTermShared(String proteinPair,
                                Map<String, Set<String>> proteinToTerms) {
        int index = proteinPair.indexOf("\t");
        String id1 = proteinPair.substring(0, index);
        String id2 = proteinPair.substring(index + 1);
        Set<String> set1 = proteinToTerms.get(id1);
        if (set1 == null)
            return false;
        Set<String> set2 = proteinToTerms.get(id2);
        if (set2 == null)
            return false;
        // Check is shared
        for (String term1 : set1) {
            if (set2.contains(term1))
                return true;
        }
        return false;
    }
    
    /**
     * This method is used to check the difference between hub proteins and non-hub proteins
     * regarding GO annotations in a GO slim way. This method is rewritten from Xin's code.
     * @throws Exception
     */
    @Test
    public void analyzeGOSlimEnrichments() throws Exception {
        // Check the whole network
        String intFileName = FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME");
        //String intFileName = "results/v2/FI73_042108.txt";
        // Get hub and non-hub proteins based on the connection degree
        Set<String> interactions = fu.loadInteractions(intFileName);
        final Map<String, Integer> proteinToDegree = InteractionUtilities.generateProteinToDegree(interactions);
        List<String> proteinList = new ArrayList<String>(proteinToDegree.keySet());
        // Sort proteins based on connection degrees
        Collections.sort(proteinList, new Comparator<String>() {
           public int compare(String protein1, String protein2) { 
               Integer degree1 = proteinToDegree.get(protein1);
               Integer degree2 = proteinToDegree.get(protein2);
               return degree2 - degree1;
           }
        });
        // Pick the first 5% as the hub proteins
        int hubSize = (int) (proteinList.size() * 0.05);
        int nonHubSize = proteinList.size() - hubSize;
        Set<String> hubs = new HashSet<String>();
        Set<String> nonHubs = new HashSet<String>();
        for (int i = 0; i < proteinList.size(); i++) {
            String protein = proteinList.get(i);
            if (i < hubSize)
                hubs.add(protein);
            else
                nonHubs.add(protein);
        }
        compareGOSlimEnricuments(hubs, nonHubs);
    }
    
    /**
     * This method is used to compare the GO slim annotations between proteins in the FI network
     * and all proteins in the SwissProt database.
     * @throws Exception
     */
    @Test
    public void compareGOSlimForFINetwork() throws Exception {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissProteins = uniAnalyzer.loadSwissProtIds();
        // Add a filter for proteins having GO annotations only to avoid the bias of functional
        // annotations
        // Map protein to GO terms
        Map<String, Set<String>> proteinToGOBP = loadProteinToGOBPTerms();
        Map<String, Set<String>> proteinToGOMF = loadProteinToGOMFTerms();
        Set<String> allGOProteins = new HashSet<String>();
        allGOProteins.addAll(proteinToGOBP.keySet());
        allGOProteins.addAll(proteinToGOMF.keySet());
        System.out.println("Total GO proteins: " + allGOProteins.size());
        System.out.println("Total swiss protein: " + swissProteins.size());
        swissProteins.retainAll(allGOProteins);
        System.out.println("SwissProt having GO: " + swissProteins.size());
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME"), // Whole network
                FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Pathway_043009.txt", // Pathway FIs
                FIConfiguration.getConfiguration().get("RESULT_DIR") + "FIs_Predicted_043009.txt" // Predicted FIs
        };
        for (String fileName : fileNames) {
            System.out.println("\nFile name: " + fileName);
            // Load proteins from all FI network
            Set<String> fis = fu.loadInteractions(fileName);
            Set<String> fiProteins = InteractionUtilities.grepIDsFromInteractions(fis);
            fiProteins.retainAll(swissProteins);
            compareGOSlimEnricuments(fiProteins, 
                                     swissProteins);
        }
    }

    private void compareGOSlimEnricuments(Set<String> hubs, 
                                          Set<String> nonHubs) throws IOException {
        // Map protein to GO terms
        Map<String, Set<String>> proteinToGOBP = loadProteinToGOBPTerms();
        Map<String, Set<String>> proteinToGOMF = loadProteinToGOMFTerms();
        Map<String, Set<String>> goToSlim = loadGOToSlim();
        // Record the total proteins having a specific terms for both hub and non-hub proteins
        Map<String, Set<String>> hubToGOBPSlim = new HashMap<String, Set<String>>();
        Map<String, Set<String>> hubToGOMFSlim = new HashMap<String, Set<String>>();
        Map<String, Set<String>> nonHubToGOBPSlim = new HashMap<String, Set<String>>();
        Map<String, Set<String>> nonHubToGOMFSlim = new HashMap<String, Set<String>>();
        for (String protein : hubs) {
            Set<String> goBP = proteinToGOBP.get(protein);
            if (goBP != null) {
                Set<String> mappedSlim = mapToSlim(goBP, goToSlim);
                hubToGOBPSlim.put(protein, mappedSlim);
            }
            Set<String> goMF = proteinToGOMF.get(protein);
            if (goMF != null) {
                Set<String> mappedSlim = mapToSlim(goMF, goToSlim);
                hubToGOMFSlim.put(protein, mappedSlim);
            }
        }
        for (String protein : nonHubs) {
            Set<String> goBP = proteinToGOBP.get(protein);
            if (goBP != null) {
                Set<String> mappedSlim = mapToSlim(goBP, goToSlim);
                nonHubToGOBPSlim.put(protein, mappedSlim);
            }
            Set<String> goMF = proteinToGOMF.get(protein);
            if (goMF != null) {
                Set<String> mappedSlim = mapToSlim(goMF, goToSlim);
                nonHubToGOMFSlim.put(protein, mappedSlim);
            }
        }
        // Want to count 
        Map<String, Set<String>> goBPSlimToHub = InteractionUtilities.switchKeyValues(hubToGOBPSlim);
        Map<String, Set<String>> goBPSlimToNonHub = InteractionUtilities.switchKeyValues(nonHubToGOBPSlim);
        Map<String, Set<String>> goMFSlimToHub = InteractionUtilities.switchKeyValues(hubToGOMFSlim);
        Map<String, Set<String>> goMFSlimToNonHub = InteractionUtilities.switchKeyValues(nonHubToGOMFSlim);
        // Want to have term names
        Map<String, String> goIdToTerm = new GODataAnalyzer().loadGOIdToTermMap();
        // Check for BP terms
        System.out.println("BP Slim terms");
        analyzeGOSlimEnrichments(hubs.size(), 
                                 nonHubs.size(),
                                 goBPSlimToHub, 
                                 goBPSlimToNonHub,
                                 goIdToTerm);
        // Check for MF terms
        System.out.println("\nMF slim terms:");
        analyzeGOSlimEnrichments(hubs.size(), 
                                 nonHubs.size(), 
                                 goMFSlimToHub, 
                                 goMFSlimToNonHub,
                                 goIdToTerm);
    }

    private void analyzeGOSlimEnrichments(int hubSize,
                                          int nonHubSize,
                                          Map<String, Set<String>> goSlimToHub,
                                          Map<String, Set<String>> goSlimToNonHub,
                                          Map<String, String> goIdToTerm) {
        final List<GOSlimEnrichmentData> list = new ArrayList<GOSlimEnrichmentData>();
        Set<String> goSlims = new HashSet<String>();
        goSlims.addAll(goSlimToHub.keySet());
        goSlims.addAll(goSlimToNonHub.keySet());
        //System.out.println("Total go slim: " + goSlims.size());
        for (String slim : goSlims) {
            Set<String> hubProteins = goSlimToHub.get(slim);
            int hubSlimSize = 0;
            if (hubProteins != null)
                hubSlimSize = hubProteins.size();
            Set<String> nonHubProteins = goSlimToNonHub.get(slim);
            int nonHubSlimSize = 0;
            if (nonHubProteins != null)
                nonHubSlimSize = nonHubProteins.size();
            double z = MathUtilities.calculateZValue(hubSlimSize, hubSize, 
                                                     nonHubSlimSize, nonHubSize);
            double hubRatio = (double) hubSlimSize / hubSize;
            double nonHubRatio = (double) nonHubSlimSize / nonHubSize;
            double pvalue = MathUtilities.calTwoTailStandardNormalPvalue(z);
            String term = goIdToTerm.get(slim);
            //System.out.println(slim + "\t" + term + "\t" + hubRatio + "\t" + nonHubRatio + "\t" + pvalue);
            GOSlimEnrichmentData data = new GOSlimEnrichmentData();
            data.id = slim;
            data.term = term;
            data.hubRatio = hubRatio;
            data.nonHubRatio = nonHubRatio;
            data.pvalue = pvalue;
            list.add(data);
        }
        Collections.sort(list, new Comparator<GOSlimEnrichmentData>() {
            public int compare(GOSlimEnrichmentData data1, GOSlimEnrichmentData data2) {
                // Check pvalue first
                int reply = data1.pvalue.compareTo(data2.pvalue);
                if (reply == 0) {
                    reply = data2.hubRatio.compareTo(data1.hubRatio);
                }
                return reply;
            }
        });
        System.out.println("ID\tTerm\tFI\tSwissProt\tpvalue");
        for (GOSlimEnrichmentData data : list) {
            System.out.println(data.id + "\t" + data.term + "\t" +
                               data.hubRatio + "\t" + data.nonHubRatio + "\t" +
                               data.pvalue);
        }
    }
    
    private Set<String> mapToSlim(Set<String> goIds, 
                                  Map<String, Set<String>> goToSlim) {
        Set<String> rtn = new HashSet<String>();
        for (String id : goIds) {
            Set<String> slim = goToSlim.get(id);
            if (slim != null)
                rtn.addAll(slim);
        }
        return rtn;
    }
    
    /**
     * This method is used to check GO terms as features based on Odds ratio.
     * @throws Exception
     */
    @Test
    public void testGOFeatures() throws Exception {
        FeatureChecker checker = new FeatureChecker();
        // Get the proteinToTerms
        List<Map<String, Set<String>>> proteinToTermsList = new ArrayList<Map<String, Set<String>>>();
        Map<String, Set<String>> proteinToTerms = loadProteinToGOBPTerms();
        proteinToTermsList.add(proteinToTerms);
        proteinToTerms = loadProteinToGOMFTerms();
        proteinToTermsList.add(proteinToTerms);
        proteinToTerms = loadProteinToGOCCTerms();
        proteinToTermsList.add(proteinToTerms);
//        if (true) {
//            for (int i = 0; i < proteinToTermsList.size(); i++) {
//                Map<String, Set<String>> map = proteinToTermsList.get(i);
//                System.out.println(i + ": " + map.size());
//            }
//            return;
//        }
        int index = 0;
        for (Map<String, Set<String>> proteinToTerms1 : proteinToTermsList) {
            if (index == 0)
                System.out.println("GO BP:");
            else if (index == 1)
                System.out.println("GO MF:");
            else if (index == 2)
                System.out.println("GO CP:");
            index ++;   
            // Required by an anonymous class.
            final Map<String, Set<String>> map = proteinToTerms1;
            PositiveChecker posChecker = new PositiveChecker() {
                public boolean isPositive(String pair) {
                    return isTermShared(pair, map);
                }
            };
            checker.checkFeatureOddsRatio(posChecker);
        }
    }
     
    /**
     * Load GO id to slim ids. A GO term can be mapped to multiple slim ids since a GO id
     * can have multiple parents in the DAG.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadGOToSlim() throws IOException {
        String mapFileName = FIConfiguration.getConfiguration().get("GO_DIR") + "goaslim.map";
        fu.setInput(mapFileName);
        String line = null;
        Map<String, Set<String>> goToSlim = new HashMap<String, Set<String>>();
        // Want to check how many slim terms are
        Set<String> slimTerms = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("!"))
                continue;
            String[] tokens = line.split("\t");
            String goTerm = tokens[0];
            String slimTerm = tokens[1];
            slimTerms.add(slimTerm);
            Set<String> set = goToSlim.get(goTerm);
            if (set == null) {
                set = new HashSet<String>();
                goToSlim.put(goTerm, set);
            }
            set.add(slimTerm);
        }
        fu.close();
        System.out.println("Total slim terms: " + slimTerms.size());
//        // Check the mapping
//        for (String go : goToSlim.keySet()) {
//            Set<String> set = goToSlim.get(go);
//            if (set.size() > 1)
//                System.out.println(go + ": " + set);
//        }
        return goToSlim;
    }
    
    private class GOSlimEnrichmentData {
        String id;
        String term;
        Double hubRatio;
        Double nonHubRatio;
        Double pvalue;
    }
    
}
