/*
 * Created on Sep 18, 2006
 *
 */
package org.reactome.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.PositiveChecker;

public class PfamAnalyzer {
    private String pFamDirName = FIConfiguration.getConfiguration().get("PFAM_DIR_NAME");
    private String uniprotDirName = FIConfiguration.getConfiguration().get("UNIPROT_DIR");
    private FileUtility fu;
    // Cache these mappings
    private Map<String, Set<String>> uni2PfamMap;
    private Set<String> intPfam;
    
    public PfamAnalyzer() {
        fu = new FileUtility();
    }
    
    public String getpFamDirName() {
        return pFamDirName;
    }

    public void setpFamDirName(String pFamDirName) {
        this.pFamDirName = pFamDirName;
    }

    public String getUniprotDirName() {
        return uniprotDirName;
    }

    public void setUniprotDirName(String uniprotDirName) {
        this.uniprotDirName = uniprotDirName;
    }

    /**
     * This method is used to check domain-domain interaction as a feature
     * using Odds ratio.
     * @throws Exception
     */
    @Test
    public void testPfamFeature() throws Exception {
        FeatureChecker checker = new FeatureChecker();
        PositiveChecker posChecker = new PositiveChecker() {
            public boolean isPositive(String pair) {
                try {
                    return checkIfInteracting(pair);
                }
                catch(IOException e) {
                    e.printStackTrace();
                }
                return false;
            }
        };
        checker.checkFeatureOddsRatio(posChecker);
    }
    
    public boolean checkIfInteracting(String uniProtPair) throws IOException {
        int index = uniProtPair.indexOf("\t");
        String id1 = uniProtPair.substring(0, index);
        String id2 = uniProtPair.substring(index + 1);
        return checkIfInteracting(id1, id2);
    }
    
    public boolean checkIfInteracting(String uniProt1,
                                      String uniProt2) throws IOException {
        if (uni2PfamMap == null)
            uni2PfamMap = loadUni2PfamMap();
        if (intPfam == null)
            intPfam = loadIntPfam();
        Set<String> pfamSet1 = uni2PfamMap.get(uniProt1);
        if (pfamSet1 == null)
            return false;
        Set<String> pfamSet2 = uni2PfamMap.get(uniProt2);
        if (pfamSet2 == null)
            return false;
        String key = null;
        for (String pfam1 : pfamSet1) {
            for (String pfam2 : pfamSet2) {
                if (pfam1.compareTo(pfam2) < 0) 
                    key = pfam1 + "\t" + pfam2;
                else
                    key = pfam2 + "\t" + pfam1;
                if (intPfam.contains(key))
                    return true;
            }
        }
        return false;
    }
    
    public Map<String, Set<String>> getUni2PfamMap() throws IOException {
        if (uni2PfamMap == null)
            uni2PfamMap = loadUni2PfamMap();
        return uni2PfamMap;
    }
    
    private Map<String, Set<String>> loadUni2PfamMap() throws IOException {
        //String fileName = FileNameManager.getManager().getFileName("Uni2Pfam.txt");
        String fileName = uniprotDirName + "Uni2Pfam.txt";
        return fu.loadSetMap(fileName);
    }
    
    public Set<String> loadIntPfam() throws IOException {
        String fileName = pFamDirName + "IntPFamIDs.txt";
        return fu.loadInteractions(fileName);
    }
    
    public void countIntPFam() throws IOException {
        Set<String> intPfam = loadIntPfam();
        Set<String> pfam = InteractionUtilities.grepIDsFromInteractions(intPfam);
        System.out.println("Total PFam used in iFam: " + pfam.size());
    }
    
    @Test
    public void convertIntToPfamIDs() throws IOException {
        //String intFileName = PFAM_DIR_NAME + "int_pfamAs.txt";
        String intFileName = pFamDirName + "pfamA_interactions.txt";
        String destFileName = pFamDirName + "IntPFamIDs.txt";
        // As of 2015, the above file is what we need actually. So
        // We just make a copy of this file to the required file name
        Files.copy(FileSystems.getDefault().getPath(intFileName), 
                   FileSystems.getDefault().getPath(destFileName), 
                   StandardCopyOption.REPLACE_EXISTING);
        
//        Map<String, String> db2Pfam = getDBId2PFamId();
//        fu.setInput(intFileName);
//        String[] tokens = null;
//        String line = null;
//        String pfam1, pfam2;
//        Set<String> intPFamSet = new HashSet<String>();
//        while ((line = fu.readLine()) != null) {
//            tokens = line.split("\t");
//            pfam1 = db2Pfam.get(tokens[0]);
//            pfam2 = db2Pfam.get(tokens[1]);
//            if (pfam1 == null || pfam2 == null) {
//                // This is very strange:
//                // This line is wrong: '3417'   '3417'
//                // throw new IllegalStateException(line + " has unmapped ids.");
//                System.out.println("Cannot be mapped: " + line);
//                continue;
//            }
////            pfam1 = removeQuote(pfam1);
////            pfam2 = removeQuote(pfam2);
//            // Need to sort since some int_ids are duplicated in file intFilename
//            if (pfam1.compareTo(pfam2) < 0)
//                intPFamSet.add(pfam1 + "\t" + pfam2);
//            else
//                intPFamSet.add(pfam2 + "\t" + pfam1);
//        }
//        fu.close();
//        fu.saveInteractions(intPFamSet, destFileName);
    }
    
    private Map<String, String> getDBId2PFamId() throws IOException {
        String pfamFileName = pFamDirName + "pfamA.txt";
        Map<String, String> db2Pfam = new HashMap<String, String>();
        String line = null;
        fu.setInput(pfamFileName);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            db2Pfam.put(tokens[0],
                        tokens[1]);
        }
        fu.close();
        return db2Pfam;
    }
    
    private String removeQuote(String txt) {
        return txt.substring(1, txt.length() - 1);
    }
    
    
    @Test
    public void checkCoverage() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Set<String> swissIds = uniAnalyzer.loadSwissProtIds();
        System.out.println("Total Swiss IDs: " + swissIds.size());
        Map<String, Set<String>> uni2PfamMap = getUni2PfamMap();
        int size = uni2PfamMap.size();
        Set<String> proteins = new HashSet<String>(uni2PfamMap.keySet());
        proteins.retainAll(swissIds);
        double percentage = (double) proteins.size() / swissIds.size();
        System.out.println("PfamDomainInt" + "\t" + size + "\t" + proteins.size() + "\t" + percentage);
        // Load Proteins in the Reactome FI
        String resultDir = FIConfiguration.getConfiguration().get("RESULT_DIR");
        Set<String> reactomeFIs = new FileUtility().loadInteractions(resultDir + File.separator + "FIs_Reactome.txt");
        Set<String> reactomeFIsProteins = InteractionUtilities.grepIDsFromInteractions(reactomeFIs);
        System.out.println("Proteins in ReactomeFIs: " + reactomeFIsProteins.size());
        Set<String> inSwissProt = new HashSet<String>(reactomeFIsProteins);
        inSwissProt.retainAll(swissIds);
        System.out.println("\tIn SwissProt: " + inSwissProt.size() + " (" + (double)inSwissProt.size() / reactomeFIsProteins.size() + ")");
        Set<String> inpFAM = new HashSet<String>(reactomeFIsProteins);
        inpFAM.retainAll(uni2PfamMap.keySet());
        System.out.println("\tIn pFAM: " + inpFAM.size() + " (" + (double) inpFAM.size() / reactomeFIsProteins.size() + ")");
        inpFAM = new HashSet<String>(reactomeFIsProteins);
        inpFAM.retainAll(proteins);
        System.out.println("\tFilter to SwissProt: " + inpFAM.size() + " (" + (double) inpFAM.size()/reactomeFIsProteins.size() + ")");
    }
    
    
}
