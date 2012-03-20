/*
 * Created on Apr 10, 2009
 *
 */
package org.reactome.tred;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.data.HPRDAnalyzer;
import org.reactome.data.ProteinIdFilters;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to do TRED database related stuff.
 * @author wgm
 *
 */
public class TREDAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public TREDAnalyzer() {
    }
    
    /**
     * Merge three RefGene files and add a new speciesCode column
     * @throws IOException
     */
    @Test
    public void processRefGeneFiles() throws IOException {
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("TRED_DIR") + "flat_pub/flat_RefGene_hg15.txt",
                FIConfiguration.getConfiguration().get("TRED_DIR") + "flat_pub/RefGene_mm3.txt",
                FIConfiguration.getConfiguration().get("TRED_DIR") + "flat_pub/RefGene_rn2.txt"
        };
        String outFileName = FIConfiguration.getConfiguration().get("TRED_DIR") + "flat_pub/RefGene.txt";
        String line = null;
        boolean isHeaderDone = false;
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            line = fu.readLine(); // header
            if (!isHeaderDone) {
                line = line + "\tspeciesCode";
                outFu.printLine(line);
            }
            // Get the species code from file name
            int index = fileName.lastIndexOf("_");
            int index1 = fileName.indexOf(".");
            String specicesCode = fileName.substring(index + 1, index1);
            while ((line = fu.readLine()) != null) {
                line = line + "\t" + specicesCode;
                outFu.printLine(line);
            }
            fu.close();
        }
        outFu.close();
    }
    
    public Set<String> loadTFTargetInteractions() throws IOException {
        return fu.loadInteractions(FIConfiguration.getConfiguration().get("TRED_FI_FILE"));
//        return fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + "TREDInteractionsInUniProt.txt");
    }
    
    /**
     * Used to count genes in the TF/Target interctions.
     * @throws IOException
     */
    @Test
    public void counts() throws IOException {
        String[] fileNames = new String[] {
                "TREDTFTargetInteractions.txt",
                "TREDTFTargetInteractions_All.txt",
                "TREDInteractionsInUniProt.txt",
                "TREDInteractionsInUniProt_All.txt"
        };
        for (String fileName : fileNames) {
            Set<String> interactions = null;
            Set<String> genes = null;
            if (fileName.contains("UniProt")) {
                interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
                genes = InteractionUtilities.grepIDsFromInteractions(interactions);
            }
            else {
                //Set<String> interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + fileName);
                Map<String, Set<String>> tfTargets = fu.loadSetMap(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
                genes = new HashSet<String>();
                interactions = new HashSet<String>();
                for (String tf : tfTargets.keySet()) {
                    genes.add(tf);
                    Set<String> targets = tfTargets.get(tf);
                    genes.addAll(targets);
                    for (String target : targets) {
                        int compare = tf.compareTo(target);
                        if (compare < 0)
                            interactions.add(tf + "\t" + target);
                        else if (compare > 0)
                            interactions.add(target + "\t" + tf);
                    }
                }
                System.out.println("Total TF: " + tfTargets.size());
            }
            System.out.println(fileName);
            System.out.println("Total interactions: " + interactions.size());
            System.out.println("Total genes: " + genes.size());
            System.out.println();
        }
    }
    
    public Map<String, String> mapGeneNameToUniProtId(Set<String> geneNames) throws Exception {
        Map<String, String> geneNameToId = loadGeneNameToIdMap();
        Map<String, String> geneIdToUniProt = new HPRDAnalyzer().loadGeneIdToUniProtAccMap();
        Map<String, String> geneNameToUniProt = new UniProtAnalyzer().loadGeneNameToUniProt();
        // Some manual mapping for TFs
        geneNameToUniProt.put("DP-1", "Q14186");
        // This factor cannot be mapped to anywhere!
        //geneNameToUniProt.put("E2F+p107", );
        Map<String, String> nameToUniProt = new HashMap<String, String>();
        Map<String, Set<String>> primaryNameToAllNames = new TREDHiberanteReader().loadPrimaryNameToAllNames();
        for (String geneName : geneNames) {
            String uniProtId = mapGeneNameToUniProt(geneName, 
                                                    primaryNameToAllNames, 
                                                    geneNameToId, 
                                                    geneIdToUniProt, 
                                                    geneNameToUniProt);
            if (uniProtId != null)
                nameToUniProt.put(geneName, uniProtId);
        }
        return nameToUniProt;
    }
    
    /**
     * Use this method to convert TF/Target interactions in names to interactions in UniProt
     * accession numbers.
     * @throws IOException
     */
    @Test
    public void createInteractionsInUniProtIds() throws Exception {
        String[] fileNames = new String[] {
                "TREDTFTargetInteractions.txt",
                "TREDTFTargetInteractions_All.txt"
        };
        Set<String> notMappedTFs = new HashSet<String>();
        Set<String> notMappedNames = new HashSet<String>();
        ProteinIdFilters filters = new ProteinIdFilters();
        for (String fileName : fileNames) {
            Set<String> interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
            Set<String> geneNames = InteractionUtilities.grepIDsFromInteractions(interactions);
            Map<String, String> geneNameToUniProt = mapGeneNameToUniProtId(geneNames);
            Map<String, Set<String>> tfTargets = fu.loadSetMap(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
            interactions.clear(); // To hold interactions in UniProts
            for (String tf : tfTargets.keySet()) {
                String uniProtAcc = geneNameToUniProt.get(tf);
                if (uniProtAcc == null) {
                    notMappedNames.add(tf);
                    notMappedTFs.add(tf);
                    continue;
                }
                Set<String> targets = tfTargets.get(tf);
                for (String target : targets) {
                    String targetUniProtAcc = geneNameToUniProt.get(target);
                    if (targetUniProtAcc == null) {
                        notMappedNames.add(target);
                        continue;
                    }
                    int compare = uniProtAcc.compareTo(targetUniProtAcc);
                    if (compare < 0)
                        interactions.add(uniProtAcc + " " + targetUniProtAcc);
                    else if (compare > 0)
                        interactions.add(targetUniProtAcc + " " + uniProtAcc);
                }
            }
            // Add a filtering
            interactions = filters.normalizeProteinPairs(interactions);
            System.out.println("total interactions: " + interactions.size());
            String outFileName = null;
            if (fileName.contains("All"))
                outFileName = "TREDInteractionsInUniProt_All.txt";
            else
                outFileName = "TREDInteractionsInUniProt.txt";
            fu.saveInteractions(interactions, 
                                FIConfiguration.getConfiguration().get("RESULT_DIR") + outFileName);
        }
        System.out.println("Not mapped names: " + notMappedNames.size());
        for (String name : notMappedNames)
            System.out.println(name);
        System.out.println("Not mapped TFs: " + notMappedTFs.size());
        for (String name : notMappedTFs)
            System.out.println(name);
    }
    
    private String mapGeneNameToUniProt(String geneName,
                                        Map<String, Set<String>> geneNameToAllNames,
                                        Map<String, String> geneNameToId,
                                        Map<String, String> geneIdToUniProt,
                                        Map<String, String> geneNameToUniProt) {
        String geneId = geneNameToId.get(geneName);
        String uniProtAcc = null;
        if (geneId != null)
            uniProtAcc = geneIdToUniProt.get(geneId);
        if (uniProtAcc == null)
            uniProtAcc = geneNameToUniProt.get(geneName);
        if (uniProtAcc == null) {
            uniProtAcc = mapGeneNameToUniProtViaAllNames(geneName,
                                                         geneNameToAllNames, 
                                                         geneNameToUniProt);
        }
        return uniProtAcc;
    }
    
    private String mapGeneNameToUniProtViaAllNames(String geneName,
                                                   Map<String, Set<String>> geneNameToAllNames,
                                                   Map<String, String> geneNameToUniProt) {
        Set<String> allNames = geneNameToAllNames.get(geneName);
        if (allNames == null)
            return null;
        for (String name : allNames) {
            String uniprot = geneNameToUniProt.get(name);
            if (uniprot != null)
                return uniprot;
        }
        return null;
    }
    
    private Map<String, String> loadGeneNameToIdMap() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("TRED_DIR") + "GeneName2UniGeneID.txt";
        return fu.importMap(fileName);
    }
    
}
