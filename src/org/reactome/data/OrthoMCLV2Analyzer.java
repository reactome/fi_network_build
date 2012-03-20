/*
 * Created on Mar 26, 2009
 *
 */
package org.reactome.data;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
import org.reactome.fi.util.ChecksumCalculator;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;


/**
 * This class is used to analyze OrthoMCL version 2 related data files.
 * @author wgm
 *
 */
public class OrthoMCLV2Analyzer extends OrthologyAnalyzer {
    // To calculate Checksum
    private ChecksumCalculator checksum;
    private FileUtility fu;
    
    public OrthoMCLV2Analyzer() {
        checksum = new ChecksumCalculator();
        fu = new FileUtility();
    }
    
    @Test
    public void testReadSeqFile() throws IOException {
        Map<String, Set<String>> groups = readSeqFile();
        // Check the output
        int c = 0;
        for (String group : groups.keySet()) {
            Set<String> members = groups.get(group);
            System.out.println(group);
            for (String member : members) {
                System.out.println("\t" + member);
            }
            c ++;
            if (c > 10)
                break;
        }
    }
    
    @Test
    public void checkHumanProteins() throws IOException {
        Map<String, Set<String>> groups = readSeqFile();
        List<String> humanProteins = new ArrayList<String>();
        for (String group : groups.keySet()) {
            Set<String> members = groups.get(group);
            for (String member : members) {
                if (member.startsWith("hsa:"))
                    humanProteins.add(member);
            }
        }
        System.out.println("Total human proteins in list: " + humanProteins.size());
        Set<String> set = new HashSet<String>(humanProteins);
        System.out.println("Total human proteins in set: " + set.size());
    }
    
    /**
     * This method will generate the OrthoMCL group file in checksum. Since we are going
     * to use this file to map proteins from yeast, celegans and drosophia to human, only
     * groups contains human proteins are in this file.
     * @throws IOException
     */
    @Test
    public void generateGroupInChecksumFile() throws IOException {
        Map<String, Set<String>> groups = readSeqFile();
        System.out.println("Total groups: " + groups.size());
        // Do some filtering
        Set<String> toBeDeleted = new HashSet<String>();
        for (String group : groups.keySet()) {
            Set<String> set = groups.get(group);
            if (set.size() < 2) {
                toBeDeleted.add(group);
            }
            else if (!isValidGroup(set)) {
                toBeDeleted.add(group);
            }
        }
        System.out.println("Total deleted groups: " + toBeDeleted.size());
        for (String group : toBeDeleted)
            groups.remove(group);
        // Output
        String outFileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "OrthoMCLHumanChecksum.txt";
        fu.saveSetMapInSort(groups, 
                            outFileName);
    }
    
    public Map<String, Set<String>> loadYeastToHumanMapInCheckcum() throws IOException {
        return loadToHumanMapInChecksum("sce");
    }
    
    public Map<String, Set<String>> loadFlyToHumanMapInChecksum() throws IOException {
        return loadToHumanMapInChecksum("dme");
    }
    
    public Map<String, Set<String>> loadWormToHumanMapInChecksum() throws IOException {
        return loadToHumanMapInChecksum("cel");
    }
    
    private Map<String, Set<String>> loadToHumanMapInChecksum(String otherSpecies) throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "OrthoMCLHumanChecksum.txt";
        Map<String, Set<String>> groups = fu.loadSetMap(fileName);
        Map<String, Set<String>> rtn = new HashMap<String, Set<String>>();
        List<String> humanCSList = new ArrayList<String>();
        List<String> otherCSList = new ArrayList<String>();
        for (String group : groups.keySet()) {
            Set<String> csSet = groups.get(group);
            if (!hasSpeciesInCSet(csSet, otherSpecies))
                continue;
            // Generate map from other species to human
            humanCSList.clear();
            otherCSList.clear();
            splitChecksumsForSpecies(humanCSList, 
                                     otherCSList, 
                                     otherSpecies, 
                                     csSet);
            // Create mapping
            for (int i = 0; i < otherCSList.size(); i++) {
                String otherCS = otherCSList.get(i);
                Set<String> set = rtn.get(otherCS);
                if (set == null) {
                    set = new HashSet<String>();
                    rtn.put(otherCS, set);
                }
                set.addAll(humanCSList);
            }
        }
        return rtn;
    }
    
    private void splitChecksumsForSpecies(List<String> humanCSList,
                                          List<String> otherCSList,
                                          String otherSpecies,
                                          Set<String> csSet) {
        for (String cs : csSet) {
            if (cs.startsWith("hsa")) {
                // This is human
                humanCSList.add(cs.substring(4));
            }
            else if (cs.startsWith(otherSpecies)) {
                // For other species
                otherCSList.add(cs.substring(4));
            }
        }
    }
    
    private boolean hasSpeciesInCSet(Set<String> csSet,
                                    String species) {
        for (String cs : csSet) {
            if (cs.startsWith(species))
                return true;
        }
        return false;
    }
    
    private boolean isValidGroup(Set<String> group) {
        // Get all species
        Set<String> species = new HashSet<String>();
        for (String cs : group) {
            species.add(cs.substring(0, 3));
        }
        if (species.contains("hsa") &&
            species.size() > 1)
            return true;
        return false;
    }
    
    private List<String> getSpecies() {
        List<String> list = new ArrayList<String>();
        list.add("hsa");
        list.add("sce");
        list.add("dme");
        list.add("cel");
        return list;
    }
    
    /**
     * This method is used to generate human checksums to annotation lines from the fast files.
     * @throws IOException
     */
    @Test
    public void generateHumanChecksumToAnnotationFile() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "seqs_orthomcl-2.fasta";
        String line = null;
        fu.setInput(fileName);
        Map<String, String> checksumToAnnotation = new HashMap<String, String>();
        StringBuilder builder = new StringBuilder();
        boolean isNeeded = false;
        String annotation = null;
        String species = null;
        String group = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith(">")) {
                if (builder.length() > 0) {
                    String seq = builder.toString();
                    if (seq.endsWith("*")) // Remove the end *. Don't know why it is there!!!
                        seq = seq.substring(0, seq.length() - 1);
                    String cs = checksum.calculateChecksum(seq);
                    checksumToAnnotation.put(cs, annotation);
                    builder.setLength(0);
                }
                isNeeded = false;
                species = line.substring(1, 4);
                if (!species.equals("hsa")) {
                    continue;
                }
                group = parseGroup(line);
                if (group.equals("no_group")) {
                    continue;
                }
                isNeeded = true;
                annotation = line;
            }
            else if (isNeeded){
                // Sequence line
                builder.append(line);
            }
        }
        fu.close();
        String outFileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "human_checksum_annotation.txt";
        fu.exportMap(checksumToAnnotation, 
                     outFileName);
    }
    
    public Map<String, Set<String>> readSeqFile() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "seqs_orthomcl-2.fasta";
        Map<String, Set<String>> groups = new HashMap<String, Set<String>>();
        List<String> speciesList = getSpecies();
        String line = null;
        fu.setInput(fileName);
        StringBuilder builder = new StringBuilder();
        boolean isNeeded = false;
        String group = null;
        String species = null;
        int totalSequence = 0;
        Map<String, Set<String>> checksumMap = new HashMap<String, Set<String>>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith(">")) {
                if (builder.length() > 0) {
                    //System.out.println(builder.toString());
                    Set<String> set = groups.get(group);
                    if (set == null) {
                        set = new HashSet<String>();
                        groups.put(group, set);
                    }
                    String seq = builder.toString();
                    if (seq.endsWith("*")) // Remove the end *. Don't know why it is there!!!
                        seq = seq.substring(0, seq.length() - 1);
                    String cs = checksum.calculateChecksum(seq);
                    set.add(species + ":" + cs);
                    builder.setLength(0);
                    totalSequence ++;
                    Set<String> csSet = checksumMap.get(cs);
                    if (csSet == null) {
                        csSet = new HashSet<String>();
                        checksumMap.put(cs, csSet);
                    }
                    csSet.add(seq);
                }
                isNeeded = false;
                species = line.substring(1, 4);
                if (!speciesList.contains(species)) {
                    continue;
                }
                group = parseGroup(line);
                if (group.equals("no_group")) {
                    continue;
                }
                isNeeded = true;
            }
            else if (isNeeded){
                // Sequence line
                builder.append(line);
            }
        }
        fu.close();
        System.out.println("Total sequence: " + totalSequence);
        System.out.println("Total checksum: " + checksumMap.size());
        // Check sequence in checksum
        // The following check makes sure there is one-to-one map between
        // checksum and aa sequence
        for (String cs : checksumMap.keySet()) {
            Set<String> csSet = checksumMap.get(cs);
            if (csSet.size() == 1)
                continue;
            System.out.println(cs);
            for (String seq : csSet)
                System.out.println("\t" + seq);
        }
        return groups;
    }
    
    private String parseGroup(String line) {
        int index1 = line.indexOf("|", 5);
        int index2 = line.indexOf("|", index1 + 1);
        return line.substring(index1 + 1, index2).trim();
    }
    
    /**
     * Load the human checksum to uniprot ids generated from the seq file.
     * @return
     * @throws IOException
     */
    private Map<String, String> extractHumanChecksumToUniProt() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "human_checksum_annotation.txt";
        Map<String, String> checksumToAnnotation = fu.importMap(fileName);
        Map<String, String> checksumToId = new HashMap<String, String>();
        int index = 0;
        for (String cs : checksumToAnnotation.keySet()) {
            String annotation = checksumToAnnotation.get(cs);
            //System.out.println(annotation);
            // ...... [Source:Uniprot/SPTREMBL;Acc:Q3KQX2]
            index = annotation.indexOf("Uniprot");
            if (index < 0)
                continue;
            index = annotation.indexOf("Acc:", index);
            // Some annotations have been truncated
            if (index + 4 >= annotation.length())
                continue;
            String id = annotation.substring(index + 4, annotation.length() - 1);
            checksumToId.put(cs, id);
        }
        return checksumToId;
    }
    
    @Test
    public void generateHumanENSPToUniProtIdMap() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "human_checksum_annotation.txt";
        Map<String, String> checksumToAnnotation = fu.importMap(fileName);
        Map<String, String> ensemblToUni = new HashMap<String, String>();
        int len = "ENSP00000265529".length();
        int index = 0;
        for (String line : checksumToAnnotation.values()) {
            String ensembl = line.substring(5, 5 + len);
            index = line.indexOf("Uniprot");
            String id = null;
            if (index < 0) {
                if (line.contains("[Source:"))
                    continue; // From other non-UniProt source
                id = queryIdFromURL(ensembl);
                if (id != null)
                    ensemblToUni.put(ensembl, id);
                continue;
            }
            index = line.indexOf("Acc:", index);
            // Some annotations have been truncated
            if (index + 4 >= line.length()) {
                id = queryIdFromURL(ensembl);
                if (id != null)
                    ensemblToUni.put(ensembl, id);
                continue;
            }
            id = line.substring(index + 4, line.length() - 1);
            ensemblToUni.put(ensembl, id);
        }
        fu.exportMap(ensemblToUni,
                     FIConfiguration.getConfiguration().get("ORTHO_MCL_DIR") + "ensembl2Uniprot.txt");
    }
    
    private String queryIdFromURL(String ensembl) throws Exception {
        String local = "http://dec2006.archive.ensembl.org/Homo_sapiens/protview?peptide=" + ensembl;
        System.out.println(local);
        URL url = new URL(local);
        InputStream is = url.openStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader reader = new BufferedReader(isr);
        String line = null;
        StringBuilder builder = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            builder.append(line).append("\n");
        }
        Pattern pattern = Pattern.compile(">Source: Uniprot/(SWISSPROT|SPTREMBL)(.+)</a");
        Matcher matcher = pattern.matcher(builder.toString());
        String id = null;
        if (matcher.find()) {
            id = matcher.group(2);
            System.out.println(id);
        }
        reader.close();
        is.close();
        return id;
    }
}
