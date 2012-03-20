/*
 * Created on Mar 31, 2009
 *
 */
package org.reactome.psi.data;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.junit.Test;
import org.reactome.fi.EnsemblAnalyzer;
import org.reactome.fi.OrthoMCLV2Analyzer;
import org.reactome.fi.ProteinIdFilters;
import org.reactome.fi.ProteinSequenceHandler;
import org.reactome.fi.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.weka.FeatureChecker;

/**
 * This class is used to handle orthology mapping for PPIs from different taxons.
 * @author wgm
 *
 */
public class PsiMiOrthologyAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public PsiMiOrthologyAnalyzer() {
    }
    
    /**
     * Normalize mapped PPIs.
     * @throws Exception
     */
    @Test
    public void normalizePPIs() throws Exception {
        String[] intFileNames = new String[] {
                "humanPPIsFromWormInUniProt.txt",
                "humanPPIsFromFlyInUniProt.txt",   
                "humanPPIsFromYeastInUniProt.txt",
        };
        ProteinIdFilters filters = new ProteinIdFilters();
        for (String fileName : intFileNames) {
            Set<String> ppis = fu.loadInteractions(FIConfiguration.getConfiguration().get("INTACT_DIR") + fileName);
            Set<String> normalized = filters.normalizeProteinPairs(ppis);
            int index = fileName.indexOf(".");
            String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName.substring(0, index) + "_Norm.txt";
            fu.outputSet(normalized, outFileName);
            System.out.printf("Before normalization: %d, after: %d%n",
                              ppis.size(),
                              normalized.size());
        }
    }
    
    private void normalizePPIs(String inFileName,
                               String outFileName) throws Exception {
        ProteinIdFilters filters = new ProteinIdFilters();
        Set<String> ppis = fu.loadInteractions(inFileName);
        Set<String> normalized = filters.normalizeProteinPairs(ppis);
        fu.outputSet(normalized, outFileName);
        System.out.printf("Before normalization: %d, after: %d%n",
                          ppis.size(),
                          normalized.size());
    }

    
    /**
     * This method is used to check the utility of the mapped human PPIs
     * based on extracted FIs from Reactome database.
     * @throws Exception
     */
    @Test
    public void testPPIsCoverage() throws Exception {
        String[] intFileNames = new String[] {
//                "humanPPIsFromCaeel.txt",       
//                "humanPPIsFromWormInUniProt.txt",
//                "humanPPIsFromDrome.txt",      
//                "humanPPIsFromFlyInUniProt.txt",   
//                "humanPPIsFromYeast.txt",
//                "humanPPIsFromYeastInUniProt.txt",
                FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"),
                FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_YEAST_FILE"),
                FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_WORM_FILE"),
                FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_FLY_FILE"),
                FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_MOUSE_FILE")
        };
        FeatureChecker checker = new FeatureChecker();
        for (String ppiName : intFileNames) {
            System.out.println("File: " + ppiName);
            Set<String> ppis = fu.loadInteractions(ppiName);
            checker.checkFeatureOddsRatio(ppis);
        }
    }

    @Test
    public void checkPPIsInChecksum() throws IOException {
        ProteinSequenceHandler handler = new ProteinSequenceHandler();
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        Map<String, String> idToChecksum = seqHandler.getUniProtIdToChecksum();
        remapPPIs(idToChecksum);
    }
    
    /**
     * This method is used to generate PPIs in SwissProt only.
     * @throws IOException
     */
    @Test
    public void checkPPIsInSwissProt() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> swissIds = uniAnalyzer.loadSwissProtIDsMap();
        remapPPIs(swissIds);
    }

    private void remapPPIs(Map<String, String> newMap) throws IOException {
        String[] intFileNames = new String[] {
                "humanPPIsFromCaeel.txt",       
                "humanPPIsFromWormInUniProt.txt",
                "humanPPIsFromDrome.txt",      
                "humanPPIsFromYeast.txt",
                "humanPPIsFromFlyInUniProt.txt",   
                "humanPPIsFromYeastInUniProt.txt",
        };
        for (String fileName : intFileNames) {
            System.out.println("File name: " + fileName);
            String inFile = FIConfiguration.getConfiguration().get("INTACT_DIR") + fileName;
            Set<String> ppis = fu.loadInteractions(inFile);
            Set<String> swissPPIs = new HashSet<String>();
            int index = 0;
            for (String ppi : ppis) {
                index = ppi.indexOf(" ");
                String id1 = ppi.substring(0, index);
                String mapped1 = newMap.get(id1);
                if (mapped1 == null)
                    continue;
                String id2 = ppi.substring(index + 1);
                String mapped2 = newMap.get(id2);
                if (mapped2 == null)
                    continue;
                int compare = mapped1.compareTo(mapped2);
                if (compare < 0)
                    swissPPIs.add(mapped1 + " " + mapped2);
                else if (compare > 0)
                    swissPPIs.add(mapped2 + " " + mapped1);
            }
            System.out.println("PPIs in UniProts: " + ppis.size());
            System.out.println("PPIs after new mapping: " + swissPPIs.size());
            System.out.println();
        }
    }
    
    /**
     * This method is used to map PPIs from mouse to human based on Ensembl Compara.
     * @throws IOException
     */
    @Test
    public void generateHumanPPIsFromMouseInUniProt() throws Exception {
        EnsemblAnalyzer analyzer = new EnsemblAnalyzer();
        Map<String, Set<String>> mouseToHuman = analyzer.loadMouseToHumanMapInUniProt();
        System.out.println("Total mouse proteins: " + mouseToHuman.size());
        generateHumanPPIsFromOtherSpeciesInUniprot(FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_PPI_FILE"), 
                                                   FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_TO_HUMAN_PPI_FILE"),
                                                   mouseToHuman);
        normalizePPIs(FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_TO_HUMAN_PPI_FILE"),
                      FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_MOUSE_FILE"));
    }
    
    /**
     * This method is used to map PPIs from yeast to human based on Ensembl Compara.
     * @throws IOException
     */
    @Test
    public void generateHumanPPIsFromYeastInUniProt() throws Exception {
        EnsemblAnalyzer analyzer = new EnsemblAnalyzer();
        Map<String, Set<String>> yeastToHuman = analyzer.loadYeastToHumanMapInUniProt();
        System.out.println("Total yeast proteins: " + yeastToHuman.size());
        generateHumanPPIsFromOtherSpeciesInUniprot(FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_PPI_FILE"), 
                                                   FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_TO_HUMAN_PPI_FILE"),
                                                   yeastToHuman);
        normalizePPIs(FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_TO_HUMAN_PPI_FILE"),
                      FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_YEAST_FILE"));
    }
    
    @Test
    public void generateHumanPPIsFromFlyInUniProt() throws Exception {
        EnsemblAnalyzer analyzer = new EnsemblAnalyzer();
        Map<String, Set<String>> flyToHuman = analyzer.loadFlyToHumanMapInUniProt();
        generateHumanPPIsFromOtherSpeciesInUniprot(FIConfiguration.getConfiguration().get("IREFINDEX_FLY_PPI_FILE"), 
                                                   FIConfiguration.getConfiguration().get("IREFINDEX_FLY_TO_HUMAN_PPI_FILE"),
                                                   flyToHuman);
        normalizePPIs(FIConfiguration.getConfiguration().get("IREFINDEX_FLY_TO_HUMAN_PPI_FILE"),
                      FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_FLY_FILE"));
    }
    
    @Test
    public void generateHumanPPIsFromWormInUniProt() throws Exception {
        EnsemblAnalyzer analyzer = new EnsemblAnalyzer();
        Map<String, Set<String>> flyToHuman = analyzer.loadWormToHumanMapInUniProt();
        generateHumanPPIsFromOtherSpeciesInUniprot(FIConfiguration.getConfiguration().get("IREFINDEX_WORM_PPI_FILE"), 
                                                   FIConfiguration.getConfiguration().get("IREFINDEX_WORM_TO_HUMAN_PPI_FILE"),
                                                   flyToHuman);
        normalizePPIs(FIConfiguration.getConfiguration().get("IREFINDEX_WORM_TO_HUMAN_PPI_FILE"),
                      FIConfiguration.getConfiguration().get("HUMAN_PPIS_FROM_WORM_FILE"));
    }
    
    /**
     * This method is used to map PPIs from yeast to human based on OrthoMCL v2 groups.
     * @throws IOException
     */
    @Test
    public void generateHumanPPIsFromYeast() throws Exception {
        String inFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "yeast_interactions_in_checksum.txt";
        String outFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "humanPPIsFromYeast.txt";
        OrthoMCLV2Analyzer orthoMcl = new OrthoMCLV2Analyzer();
        Map<String, Set<String>> yeastToHuman = orthoMcl.loadYeastToHumanMapInCheckcum();
        generateHumanPPIsFromOtherSpecies(inFileName, 
                                          outFileName,
                                          yeastToHuman);
    }
    
    @Test
    public void generateHumanPPIsFromFly() throws IOException {
        String inFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "drome_interactions_in_checksum.txt";
        String outFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "humanPPIsFromDrome.txt";
        OrthoMCLV2Analyzer orthoMcl = new OrthoMCLV2Analyzer();
        Map<String, Set<String>> flyToHumanMap = orthoMcl.loadFlyToHumanMapInChecksum();
        generateHumanPPIsFromOtherSpecies(inFileName,
                                          outFileName, 
                                          flyToHumanMap);
    }
    
    @Test
    public void generateHumanPPIsFromWorm() throws IOException {
        String inFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "caeel_interactions_in_checksum.txt";
        String outFileName = FIConfiguration.getConfiguration().get("INTACT_DIR") + "humanPPIsFromCaeel.txt";
        OrthoMCLV2Analyzer orthoMcl = new OrthoMCLV2Analyzer();
        Map<String, Set<String>> wormToHumanMap = orthoMcl.loadWormToHumanMapInChecksum();
        generateHumanPPIsFromOtherSpecies(inFileName, 
                                          outFileName,
                                          wormToHumanMap);
    }
    
    private void generateHumanPPIsFromOtherSpeciesInUniprot(String inFileName,
                                                            String outFileName,
                                                            Map<String, Set<String>> toHumanMap) throws IOException {
        // Need to map to unique ids
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> idMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> otherPPIs = fu.loadInteractions(inFileName);
        Set<String> humanPPIs = new HashSet<String>();
        int index = 0;
        int compare = 0;
        for (String otherPPI : otherPPIs) {
            index = otherPPI.indexOf("\t");
            String yeastId1 = otherPPI.substring(0, index);
            Set<String> humanIds1 = toHumanMap.get(yeastId1);
            if (humanIds1 == null)
                continue;
            String yeastId2 = otherPPI.substring(index + 1);
            Set<String> humanIds2 = toHumanMap.get(yeastId2);
            if (humanIds2 == null)
                continue;
            // Create human PPIs
            for (String humanId1 : humanIds1) {
                String mapped1 = idMap.get(humanId1);
                if (mapped1 == null)
                    continue;
                for (String humanId2 : humanIds2) {
                    String mapped2 = idMap.get(humanId2);
                    if (mapped2 == null)
                        continue;
                    compare = mapped1.compareTo(mapped2);
                    if (compare < 0)
                        humanPPIs.add(mapped1 + "\t" + mapped2);
                    else if (compare > 0)
                        humanPPIs.add(mapped2 + "\t" + mapped1);
                }
            }
        }
        // Total human PPIs in checksums
        System.out.printf("Mapped human PPIs from other species: %d -> %d%n",
                          otherPPIs.size(),
                          humanPPIs.size());
        fu.outputSet(humanPPIs, 
                     outFileName);
    }

    private void generateHumanPPIsFromOtherSpecies(String inFileName,
                                                   String outFileName,
                                                   Map<String, Set<String>> toHumanMap) throws IOException {
        Set<String> yeastPPIs = fu.loadInteractions(inFileName);
        Set<String> humanPPIsInChecksum = new HashSet<String>();
        int index = 0;
        int compare = 0;
        for (String yeastPPI : yeastPPIs) {
            index = yeastPPI.indexOf(" ");
            String yeastCs1 = yeastPPI.substring(0, index);
            Set<String> humanCSList1 = toHumanMap.get(yeastCs1);
            if (humanCSList1 == null)
                continue;
            String yeastCs2 = yeastPPI.substring(index + 1);
            Set<String> humanCSList2 = toHumanMap.get(yeastCs2);
            if (humanCSList2 == null)
                continue;
            // Create human PPIs
            for (String humanCS1 : humanCSList1) {
                for (String humanCS2 : humanCSList2) {
                    compare = humanCS1.compareTo(humanCS2);
                    if (compare < 0)
                        humanPPIsInChecksum.add(humanCS1 + " " + humanCS2);
                    else if (compare > 0)
                        humanPPIsInChecksum.add(humanCS2 + " " + humanCS1);
                }
            }
        }
        // Total human PPIs in checksums
        System.out.println("Total human PPIs in checksums: " + humanPPIsInChecksum.size());
        //fu.outputSet(humanPPIsInChecksum, 
        //             outFileName);
        // Want to convert checksums to UniProt ids
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        Map<String, String> checksumToId = seqHandler.getChecksumToUniprotId();
        //Map<String, String> checksumToId = extractHumanChecksumToUniProt();
        Set<String> humanPPIsInUniProt = new HashSet<String>();
        Set<String> notMappedChecksums = new HashSet<String>();
        for (String ppi : humanPPIsInChecksum) {
            index = ppi.indexOf(" ");
            String cs1 = ppi.substring(0, index);
            String id1 = checksumToId.get(cs1);
            if (id1 == null) {
                notMappedChecksums.add(cs1);
                continue;
            }
            String cs2 = ppi.substring(index + 1);
            String id2 = checksumToId.get(cs2);
            if (id2 == null) {
                notMappedChecksums.add(cs2);
                continue;
            }
            compare = id1.compareTo(id2);
            if (compare < 0)
                humanPPIsInUniProt.add(id1 + " " + id2);
            else if (compare > 0)
                humanPPIsInUniProt.add(id2 + " " + id1);
        }
        System.out.println("Total human PPIs in UniProt ids: " + humanPPIsInUniProt.size());
        fu.outputSet(humanPPIsInUniProt, 
                     outFileName);
        System.out.println("Total not mapped checksums: " + notMappedChecksums.size());
        Set<String> totalChecksums = InteractionUtilities.grepIDsFromInteractions(humanPPIsInChecksum);
        System.out.println("Total human checksums: " + totalChecksums.size());
        // Print ten to see why these cannot be mapped
        int c = 0;
        for (String cs : notMappedChecksums) {
            System.out.println(cs);
            c ++;
            if (c > 10)
                break;
        }
    }
    
}
