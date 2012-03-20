/*
 * Created on Jan 18, 2012
 *
 */
package org.reactome.panther;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to map from Panther pathway component to UniProt access numbers.
 * @author gwu
 *
 */
public class PantherIdToUniProtMapper {
    
    public PantherIdToUniProtMapper() {
        
    }
    
    private List<String> getReliableConfidenceCodes() {
        String[] codes = new String[] {
                "IGI",
                "IPI",
                "IDA",
                "IEP",
                "TAS",
                "IC",
                "IMP",
                "RCA"
        };
        return Arrays.asList(codes);
    }
    
    /**
     * Load the PathwayComponent to UniProt id mappings from the mapping file downloaded from the 
     * pantherdb web site.
     * @param uniIdMap the human UniProt ids. This map is used to filter to human ids only, and also
     * used to filter redundant mapping.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadMapping(Map<String, String> uniIdMap) throws IOException {
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        String line = null;
        FileReader fileReader = new FileReader(FIConfiguration.getConfiguration().get("PANTHER_MAPPING_FILE"));
        BufferedReader reader = new BufferedReader(fileReader);
        // Used to filter no-reliable annotation
        List<String> reliableConfidenceCodes = getReliableConfidenceCodes();
        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\t");
            String pantherId = tokens[2];
            String uniIdToken = tokens[4];
            String uniId = getHumanUniProtId(uniIdToken);
            if (!uniIdMap.containsKey(uniId))
                continue; // Filter to human UniProt ids only
            String confidenceCode = tokens[6];
            if (!reliableConfidenceCodes.contains(confidenceCode))
                continue; // Filter to reliable annotations only
            String evidence = tokens[7];
            if (evidence == null || evidence.length() == 0)
                continue; // Filter to  mapping having PubMed or OMIM support.
            // Get the non-redundant
            uniId = uniIdMap.get(uniId);
            Set<String> ids = map.get(pantherId);
            if (ids == null) {
                ids = new HashSet<String>();
                map.put(pantherId, ids);
            }
            ids.add(uniId);
        }
        reader.close();
        fileReader.close();
        return map;
    }
    
    /**
     * This method is used for version 3.0.1. Previous version doesn't need to use this method
     * since UniProt accession number is listed there directly.
     * @param token
     * @return
     */
    private String getHumanUniProtId(String token) {
        if (!token.startsWith("HUMAN"))
            return null; // Not for human
        String key = "UniProtKB=";
        if (!token.contains(key))
            return null; // No UniProt available
        int index = token.indexOf(key);
        return token.substring(index + key.length());
    }
    
    public void testPantherToUniProtMap() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniIDMap = uniAnalyzer.loadUniProtIDsMap();
        Map<String, Set<String>> panther2UniMap = loadMapping(uniIDMap);
        // Get all mapped UniIDs
        Set<String> ids = new HashSet<String>();
        for (Iterator<String> it = panther2UniMap.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            ids.addAll(panther2UniMap.get(id));
        }
        System.out.println("Total Mapped UniProt Ids: " + ids.size());
    }
    
    /**
     * This method is used to check the usage of the evidence codes.
     * @throws IOException
     */
    @Test
    public void checkEvidenceCodeUsage() throws IOException {
        Map<String, List<EvidenceUsage>> componentToIdEv = new HashMap<String, List<EvidenceUsage>>();
        String line = null;
        FileUtility fu = new FileUtility();
        fu.setInput(FIConfiguration.getConfiguration().get("PANTHER_MAPPING_FILE"));
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String component = tokens[2];
            String id = tokens[4];
            String ec = tokens[6];
            List<EvidenceUsage> list = componentToIdEv.get(component);
            if (list == null) {
                list = new ArrayList<EvidenceUsage>();
                componentToIdEv.put(component, list);
            }
            EvidenceUsage usage = new EvidenceUsage();
            usage.evidenceCode = ec;
            usage.uniProtId = id;
            usage.evidence = tokens[7];
            usage.evidenceType = tokens[8];
            list.add(usage);
        }
        System.out.println("Total pathway components: " + componentToIdEv.size());
        // Check the total UniProt ids used
        Set<String> uniProtIds = new HashSet<String>();
        for (String component : componentToIdEv.keySet()) {
            List<EvidenceUsage> list = componentToIdEv.get(component);
            for (EvidenceUsage usage : list) {
                uniProtIds.add(usage.uniProtId);
            }
        }
        System.out.println("Total uniport ids: " + uniProtIds.size());
        // Check the usage of evidence codes
        // If we use the experimental evidence plus curators assigment
        uniProtIds.clear();
        // The following confidence should be reliable
        List<String> reliableCodes = getReliableConfidenceCodes();
        // Check if there is any reliable mapping exisiting for a pathway component
        boolean isFound = false;
        int totalNoMapped = 0;
        int totalNoEvidence = 0;
        for (String component : componentToIdEv.keySet()) {
            List<EvidenceUsage> list = componentToIdEv.get(component);
            isFound = false;
            for (EvidenceUsage usage : list) {
                if (reliableCodes.contains(usage.evidenceCode)) {
                    if (usage.evidence.length() == 0) {
//                       System.out.println("Has null evidence: " + 
//                                          component + ", " + 
//                                          usage.uniProtId + ", " + 
//                                          usage.evidenceCode); 
                       totalNoEvidence ++;
                       continue;
                    }
                    isFound = true;
                    uniProtIds.add(usage.uniProtId);
                }
            }
            if (!isFound) {
                //System.out.println(component + " cannot find a reliable mapping!");
                totalNoMapped ++;
            }
        }
        System.out.println("Total uniprot ids after reliable mapping filtering: " + uniProtIds.size());
        System.out.println("Total not mapping: " + totalNoMapped);
        System.out.println("Total no evidence: " + totalNoEvidence);
    }
    
    @Test
    public void checkEvidenceCodes() throws IOException {
        Set<String> evidenceCodes = new HashSet<String>();
        String line = null;
        FileUtility fu = new FileUtility();
        fu.setInput(FIConfiguration.getConfiguration().get("PANTHER_MAPPING_FILE"));
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String ec = tokens[6];
            if (ec.length() > 0)
                evidenceCodes.add(ec);
        }
        System.out.println("Total evidence codes: " + evidenceCodes.size());
        for (String ev : evidenceCodes)
            System.out.println(ev);
    }
    
    private class EvidenceUsage {
        String uniProtId;
        String evidenceCode;
        String evidence;
        String evidenceType;
    }
    
}
