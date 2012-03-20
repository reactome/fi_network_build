/*
 * Created on Apr 16, 2008
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.data.GODataAnalyzer;
import org.reactome.data.MicroarrayDataAnalyzer;
import org.reactome.data.ProteinSequenceHandler;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.fi.util.Value;
import org.reactome.weka.WEKADataAnalyzer;
import org.reactome.weka.WEKAResultAnalyzer;

import weka.core.Instances;

/**
 * This class is a group of methods that are used to count numbers of proteins and protein
 * pairs from any files or FIs extracted from pathways directly.
 * This class is refactored from FuctionalInteractionAnalyzer.
 * @author wgm
 */
public class ProteinAndInteractionCount {
    private FileUtility fu;
    private FunctionalInteractionAnalyzer fiAnalyzer;
    // Cache these maps for counting
    private Map<String, String> swissProtACIDMap;
    private int totalSwiss;
    private Map<String, String> uniProtACIDMap;
    private int totalUniProt;
    
    public ProteinAndInteractionCount() throws Exception {
        fu = new FileUtility();
        fiAnalyzer = new FunctionalInteractionAnalyzer();
        //fiAnalyzer.setUp();
    }
    
    public void checkReduncyInInteractions() throws Exception {
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("RESULT_DIR") + "BINDInteractions020507.txt",
                FIConfiguration.getConfiguration().get("RESULT_DIR") + "IntActInteractions020507.txt",
                FIConfiguration.getConfiguration().get("RESULT_DIR") + "HPRDInteractions020507.txt",
                "results/interaction/OrthoInteractions.txt",
                "results/interaction/YeastInteractions.txt"
        };
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = acIdMap.keySet();
        Set<String> pathwayInteractions = new HashSet<String>();
        for (String pInt : fileNames) {
            Set<String> interactions = fu.loadInteractions(pInt);
            System.out.println(pInt);
            fiAnalyzer.filterNonHumanIds(interactions,
                              uniSet, 
                              uniAnalyzer);
            report(pathwayInteractions, interactions);
            // Filter interactions
            interactions = fiAnalyzer.filterRedundencyInteractions(interactions, acIdMap);
            System.out.println("After redundancy filtering...");
            report(pathwayInteractions, interactions);
        }
    }
    
    /**
     * Check a list of PPIs only
     * @throws Exception
     */
    @Test
    public void checkPPIsBasedOnSequence() throws Exception {
//        String[] fileNames = new String[] {
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "HumanPPIs_intact.txt",
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "HumanPPIs_HPRD.txt",
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "HumanPPIs_BioGrid.txt",
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "humanPPIsFromYeastInUniProt_Norm.txt",
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "humanPPIsFromWormInUniProt_Norm.txt",
//                FIConfiguration.getConfiguration().get("RESULT_DIR + "humanPPIsFromFlyInUniProt_Norm.txt",
//                FIConfiguration.getConfiguration().get("GENEWAYS_DIR + "HumanInteractionsInUniProt.txt"
//        };
        String[] fileNames = new String[] {
                "FIs_BioCarta - Imported by PID.txt",
                "FIs_KEGG.txt",
                "FIs_The Cancer Cell Map.txt",
                "FIs_Pathway Interaction Database.txt",    
                "FIs_pantherdb.txt",
                "FIs_Reactome.txt"
        };
        // Used to consolidate sequences
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        Set<String> emptySet = new HashSet<String>();
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = acIdMap.keySet();
        Map<String, Set<String>> nameToPairs = new HashMap<String, Set<String>>();
        for (String fileName : fileNames) {
            Set<String> interactions = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
            nameToPairs.put(fileName, interactions);
        }
        // Add two gene expression
        MicroarrayDataAnalyzer arrayAnalyzer = new MicroarrayDataAnalyzer();
        Set<String> pavlidisSet = arrayAnalyzer.loadCoExpFromPavlidis();
        nameToPairs.put("Pavlidis Co-Exp", pavlidisSet);
        Set<String> carlosSet = arrayAnalyzer.loadCoExpFromPrietoCarlos();
        nameToPairs.put("Carlos Co-Exp", carlosSet);
        for (String name : nameToPairs.keySet()) {
            System.out.println(name);
            Set<String> interactions = nameToPairs.get(name);
            fiAnalyzer.filterNonHumanIds(interactions,
                                         uniSet, 
                                         uniAnalyzer);
            report(emptySet, interactions);
            System.out.println("Filtering redunduncies...");
            interactions = fiAnalyzer.filterRedundencyInteractions(interactions, acIdMap);
            report(emptySet, interactions);
            System.out.println("Consolidating based on sequences...");
            interactions = seqHandler.consolidateInteractionsUseChecksum(interactions);
            report(emptySet, interactions);
        }
    }
    
    /**
     * All FI files should have been pre-processed to remove any redundancies.
     * @throws Exception
     */
    @Test
    public void countProteinsAndFIsFromPathwayDBs() throws Exception {
        String[] fileNames = new String[] {
                "FIs_Reactome.txt",
                "FIs_pantherdb.txt",
                "FIs_The Cancer Cell Map.txt",
                "FIs_Pathway Interaction Database.txt",
                "FIs_BioCarta - Imported by PID.txt",
                "FIs_KEGG.txt",
                "TREDInteractionsInUniProt.txt"
        };
        Set<String> swissProtIds = new UniProtAnalyzer().loadSwissProtIds();
        int total = swissProtIds.size();
        System.out.println("total swiss prot: " + total);
        System.out.println("\nDatabase\tProteins\tSwissProt\tcoverage\tInteractions");
        for (String fileName : fileNames) {
            Set<String> fis = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
            Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(fis);
            int proteinTotal = proteins.size();
            proteins.retainAll(swissProtIds);
            int swissprotTotal = proteins.size();
            double percentage = (double)swissprotTotal / total;
            System.out.println(fileName + "\t" + proteinTotal + "\t" + 
                               swissprotTotal + "\t" + percentage + "\t" +
                               fis.size());
        }
    }
       
    /**
     * Use this method to check the numbers of proteins and interactions from
     * FI files. This method is used to count these numbers for individual files.
     * @throws Exception
     */
    @Test
    public void countProteinAndFIFromFiles() throws Exception {
        String dirName = FIConfiguration.getConfiguration().get("RESULT_DIR");
        String[] interactionFileNames = new String[] {
                "ReactomeInteractions020507.txt",
                "PantherInteractions020507.txt",
                "CellMapInteractions020507.txt",
                "INOHInteractions020507.txt",
                "NciNatureCuratedInteractions020507.txt",
                "NciNatureBiCartaInteractions020507.txt",
                "KEGGInteractions020507.txt",
                "BINDInteractions020507.txt",
                "IntActInteractions020507.txt",
                "HPRDInteractions020507.txt",
                FIConfiguration.getConfiguration().get("ORTHO_INTERAACTION_FILE_NAME"),
                FIConfiguration.getConfiguration().get("YEAST_INTERACTION_FILE_NAME")
        };
        //Set<String> humanInteractions = new HashSet<String>();
        Set<String> emptySet = new HashSet<String>();
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = acIdMap.keySet();
        // Used to consolidate sequences
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        
        for (String pInt : interactionFileNames) {
            Set<String> interactions = null;
            if (pInt.contains("/"))
                interactions = fu.loadInteractions(pInt);
            else
                interactions = fu.loadInteractions(dirName + pInt);
            System.out.println(pInt);
            fiAnalyzer.filterNonHumanIds(interactions,
                              uniSet, 
                              uniAnalyzer);
            report(emptySet, interactions);
            System.out.println("Filtering redunduncies...");
            interactions = fiAnalyzer.filterRedundencyInteractions(interactions, acIdMap);
            report(emptySet, interactions);
            System.out.println("Consolidating based on sequences...");
            interactions = seqHandler.consolidateInteractionsUseChecksum(interactions);
            report(emptySet, interactions);
        }
        // The following statements are used to calculate numbers for all pathway FIs
        Set<String> pathwayFIs = new HashSet<String>();
        for (int i = 0; i < 7; i++) {
            Set<String> interactions = fu.loadInteractions(dirName + interactionFileNames[i]);
            pathwayFIs.addAll(interactions);
        }
        fiAnalyzer.filterNonHumanIds(pathwayFIs, uniSet, uniAnalyzer);
        System.out.println("Results for all pathway FIs:");
        report(pathwayFIs, emptySet);
        pathwayFIs = fiAnalyzer.filterRedundencyInteractions(pathwayFIs, acIdMap);
        System.out.println("Result for all pathway FIs after filtering redundancy:");
        pathwayFIs = seqHandler.consolidateInteractionsUseChecksum(pathwayFIs);
        report(pathwayFIs, emptySet);
        // For other data types
        // From GO annotations
        GODataAnalyzer goAnalyzer = new GODataAnalyzer();
        Set<String> goIds = goAnalyzer.getAnnotatedGOBPProteins();
        System.out.println("Total GO IDs: " + goIds.size());
        countVsSwissProt(goIds);
        // From Gene Expression
        MicroarrayDataAnalyzer geneAnalyzer = new MicroarrayDataAnalyzer();
        Set<String> geneIds = geneAnalyzer.getGeneExpPairWiseData();
        System.out.println("Gene Exp: ");
        report(emptySet, geneIds);
        fiAnalyzer.filterNonHumanIds(geneIds, uniSet, uniAnalyzer);
        System.out.println("Filter non-human ids:");
        report(emptySet, geneIds);
        System.out.println("Filter redundancy:");
        geneIds = fiAnalyzer.filterRedundencyInteractions(geneIds, acIdMap);
        report(emptySet, geneIds);
        System.out.println("Consolidating based on sequences...");
        geneIds = seqHandler.consolidateInteractionsUseChecksum(geneIds);
        report(emptySet, geneIds);
    }
    
    /**
     * This method is used to check PPI recall
     * @throws Exception
     */
    @Test
    public void checkPPIRecall() throws Exception {
        String dirName = FIConfiguration.getConfiguration().get("RESULT_DIR");
        String[] pathwayIntFileNames = new String[] {
                "ReactomeInteractions020507.txt",
                "PantherInteractions020507.txt",
                "INOHInteractions020507.txt",
                "NciNatureCuratedInteractions020507.txt",
                "NciNatureBiCartaInteractions020507.txt",
                "KEGGInteractions020507.txt",
                "CellMapInteractions020507.txt"
        };
        String[] ppiFileNames = new String[] {
                "BINDInteractions020507.txt",
                "IntActInteractions020507.txt",
                "HPRDInteractions020507.txt",
        };
        Set<String> pathwayInteractions = new HashSet<String>();
        Set<String> humanInteractions = new HashSet<String>();
        for (String pInt : pathwayIntFileNames) {
            Set<String> set = fu.loadInteractions(dirName + pInt);
            pathwayInteractions.addAll(set);
        }
        for (String ppi : ppiFileNames) {
            Set<String> set = null;
            set = fu.loadInteractions(dirName + ppi);
            humanInteractions.addAll(set);
        }
        // UniProtAnalyzer is used to filter
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = uniIdMap.keySet();
        // For sequence handler
        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        // Do some filters for pathways
        fiAnalyzer.filterNonHumanIds(pathwayInteractions, uniSet, uniAnalyzer);
        pathwayInteractions = fiAnalyzer.filterRedundencyInteractions(pathwayInteractions, uniIdMap);
        pathwayInteractions = seqHandler.consolidateInteractionsUseChecksum(pathwayInteractions);
        // Do some filters for PPIs
        fiAnalyzer.filterNonHumanIds(humanInteractions, uniSet, uniAnalyzer);
        humanInteractions = fiAnalyzer.filterRedundencyInteractions(humanInteractions, uniIdMap);
        humanInteractions = seqHandler.consolidateInteractionsUseChecksum(humanInteractions);
        // Starting counting...
        countPathwayAndPPIUniProts(pathwayInteractions, humanInteractions);
    }
    
    /**
     * Use this method to count the numbers of proteins and interactions for FIs extracted
     * from pathways and PPIs or other pair-wise relationships.
     * @throws Exception
     */
    @Test
    public void countPathwayAndHumanPPIUniProtsFromFiles() throws Exception {
//        BINDInteractions020507.txt              KEGGInteractions020507.txt
//        CellMapInteractions020507.txt           NciNatureBiCartaInteractions020507.txt
//        HPRDInteractions020507.txt              NciNatureCuratedInteractions020507.txt
//        INOHInteractions020507.txt              PantherInteractions020507.txt
//        IntActInteractions020507.txt            ReactomeInteractions020507.txt
        String dirName = "results/v2/";
        String[] pathwayIntFileNames = new String[] {
//                "ReactomeInteractions020507.txt",
//                "PantherInteractions020507.txt",
//                "INOHInteractions020507.txt",
//                "NciNatureCuratedInteractions020507.txt",
//                "NciNatureBiCartaInteractions020507.txt",
//                "KEGGInteractions020507.txt",
//                "CellMapInteractions020507.txt"
                "FIs_Reactome.txt",
                "FIs_pantherdb.txt",
                "FIs_The Cancer Cell Map.txt",
                "FIs_Pathway Interaction Database.txt",
                "FIs_BioCarta - Imported by PID.txt",
                "FIs_KEGG.txt"
        };
        String[] ppiFileNames = new String[] {
//                "BINDInteractions020507.txt",
//                "IntActInteractions020507.txt",
//                "HPRDInteractions020507.txt",
//                "results/interaction/OrthoInteractions.txt",
//                "results/interaction/YeastInteractions.txt",
//                "ScoredPairsFromGOBP.txt"
        };
        Set<String> pathwayInteractions = new HashSet<String>();
        Set<String> humanInteractions = new HashSet<String>();
        for (String pInt : pathwayIntFileNames) {
            Set<String> set = fu.loadInteractions(dirName + pInt);
            pathwayInteractions.addAll(set);
        }
        for (String ppi : ppiFileNames) {
            Set<String> set = null;
            if (ppi.contains("/"))
                set = fu.loadInteractions(ppi);
            else
                set = fu.loadInteractions(dirName + ppi);
            humanInteractions.addAll(set);
        }
        MicroarrayDataAnalyzer geneAnalyzer = new MicroarrayDataAnalyzer();
        Set<String> geneIds = geneAnalyzer.getGeneExpPairWiseData();
        humanInteractions.addAll(geneIds);
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> uniSet = uniIdMap.keySet();
        fiAnalyzer.filterNonHumanIds(pathwayInteractions, uniSet, uniAnalyzer);
        pathwayInteractions = fiAnalyzer.filterRedundencyInteractions(pathwayInteractions, uniIdMap);
        fiAnalyzer.filterNonHumanIds(humanInteractions, uniSet, uniAnalyzer);
        humanInteractions = fiAnalyzer.filterRedundencyInteractions(humanInteractions, uniIdMap);
        countPathwayAndPPIUniProts(pathwayInteractions, humanInteractions);
    }
    
    /**
     * A helper method to count numbers of proteins and interactions based on a trained
     * NBC.
     * @param pathwayInteractions
     * @param ppiInteractions
     * @throws Exception
     */
    private void countPathwayAndPPIUniProts(Set<String> pathwayInteractions,
                                            Set<String> ppiInteractions) throws Exception {
        // Check human interactions
        // Combined ids
        report(pathwayInteractions, ppiInteractions);
        double[] cutoffs = new double[] {
                0.5d, 0.6d, 0.7d, 0.8d, 0.9d,
                //0.68, 0.72, 0.73, 0.74, 0.75
                //0.53, 0.55, 0.57
                0.73
        };
        // Use Ids
        Map<String, Value> valueMap = fiAnalyzer.generateValues(ppiInteractions);
        System.out.println("Value map size: " + valueMap.size());
        // Assign values
        new WEKADataAnalyzer().generateDataSet(valueMap);
        // Need an empty dataset
        Instances dataset = new WEKAResultAnalyzer().createDataSet();
        // Fill in data
        int trueCounter = 0;
        int falseCounter = 0;
        Set<String> truePair = new HashSet<String>();
        for (double cutoff : cutoffs) {
            trueCounter = falseCounter = 0;
            truePair.clear();
            for (Iterator<String> it = valueMap.keySet().iterator(); it.hasNext();) {
                String pair = it.next();
                Value value = valueMap.get(pair);
                fiAnalyzer.classify(value, dataset, cutoff);
                if (value.functionalInteraction) {
                    trueCounter ++;
                    truePair.add(pair);
                }
                else
                    falseCounter ++;
            }
            System.out.println("Cutoff: " + cutoff);
            System.out.printf("    True: %d, False: %d, total: %d, recall: %f%n",
                              trueCounter, falseCounter, valueMap.size(), (double)trueCounter / valueMap.size());
            report(pathwayInteractions, truePair);
//            if (cutoff == 0.6d) {
//                Set<String> total = new HashSet<String>(pathwayInteractions);
//                total.addAll(truePair);
//                new FileUtility().saveInteractions(total, "results/FIInteractions120506_06.txt");
//            }
        }
    }
    
    /**
     * A simple method to print out numbers of proteins and interactions in the passed
     * interaction files.
     * @param reactomeInteractions
     * @param humanInteractions
     * @throws IOException
     */
    private void report(Set<String> reactomeInteractions,
                        Set<String> humanInteractions) throws IOException {
        // Combined ids
        Set<String> combinedIds = new HashSet<String>();
        Set<String> reactomeIds = InteractionUtilities.grepIDsFromInteractions(reactomeInteractions);
        Set<String> humanIntIds = InteractionUtilities.grepIDsFromInteractions(humanInteractions);
        combinedIds.addAll(reactomeIds);
        combinedIds.addAll(humanIntIds);
        System.out.println("Total Combined IDs: " + combinedIds.size());
        System.out.println("    ids from Reactome: " + reactomeIds.size());
        System.out.println("    ids from human interactions: " + humanIntIds.size());
        countVsSwissProt(combinedIds);
        countVsUniProt(combinedIds);
        // combined interactions
        Set<String> combinedInteractions = new HashSet<String>();
        combinedInteractions.addAll(reactomeInteractions);
        combinedInteractions.addAll(humanInteractions);
        System.out.println("Total Combined interactions: " + combinedInteractions.size());
        System.out.println("    interactions from Reactome: " + reactomeInteractions.size());
        System.out.println("    interactions from humanPPI: " + humanInteractions.size());
        System.out.println();
    }
    
    public void countVsSwissProt(Set<String> ids) throws IOException {
        if (swissProtACIDMap == null) {
            swissProtACIDMap = new UniProtAnalyzer().loadSwissProtIDsMap();
            totalSwiss = new HashSet<String>(swissProtACIDMap.values()).size();
        }
        Set<String> swissIds = new HashSet<String>();
        for (String id : ids) {
            String swissId = swissProtACIDMap.get(id);
            if (swissId != null)
                swissIds.add(swissId);
            else {
                // Try to use alternative form
                int index = id.indexOf("-");
                if (index > 0)
                    id = id.substring(0, index);
                swissId = swissProtACIDMap.get(id);
                if (swissId != null)
                    swissIds.add(swissId);
            }
        }
        System.out.println("    SwissProt IDs: " + swissIds.size());
        System.out.println("    SwissProt Coverage: " + (double)swissIds.size() / totalSwiss);
    }
    
    private void countVsUniProt(Set<String> ids) throws IOException {
        if (uniProtACIDMap == null) {
            uniProtACIDMap = new UniProtAnalyzer().loadUniProtIDsMap();
            totalUniProt = new HashSet<String>(uniProtACIDMap.values()).size();
        }
        Set<String> uniIds = new HashSet<String>();
        for (String id : ids) {
            String swissId = uniProtACIDMap.get(id);
            if (swissId != null)
                uniIds.add(swissId);
            else {
                // Try to use alternative form
                int index = id.indexOf("-");
                if (index > 0)
                    id = id.substring(0, index);
                swissId = swissProtACIDMap.get(id);
                if (swissId != null)
                    uniIds.add(swissId);
            }
        }
        System.out.println("    UniProt IDs: " + uniIds.size());
        System.out.println("    UniProt Coverage: " + (double)uniIds.size() / totalUniProt);
    }
    
    /**
     * This method is used to calculate protein and PPI overlapping among three PPI
     * databases.
     * @throws Exception
     */
    @Test
    public void countPPIOverlapping() throws Exception {
        String[] fileNames = new String[] {
//                "BINDInteractions020507.txt",
//                "IntActInteractions020507.txt",
//                "HPRDInteractions020507.txt"
                "HumanPPIs_intact.txt",
                "HumanPPIs_HPRD.txt",
                "HumanPPIs_BioGrid.txt"
        };
        Map<String, Set<String>> fileToPPIs = new HashMap<String, Set<String>>();
        Map<String, Set<String>> fileToIds = new HashMap<String, Set<String>>();
        // Load PPIs and IDs
//        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
//        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
//        Set<String> uniSet = acIdMap.keySet();
//        // For sequence filtering
//        ProteinSequenceHandler seqHandler = new ProteinSequenceHandler();
        for (String file : fileNames) {
            System.out.println("File Name: " + file);
            Set<String> ppis = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + file);
            // Need to filter out no-human IDs
            //fiAnalyzer.filterNonHumanIds(ppis, uniSet, uniAnalyzer);
            // Need to filter out redundancy ids
            //ppis = fiAnalyzer.filterRedundencyInteractions(ppis, acIdMap);
            //ppis = seqHandler.consolidateInteractionsUseChecksum(ppis);
            Set<String> ids = InteractionUtilities.grepIDsFromInteractions(ppis);
            fileToIds.put(file, ids);
            fileToPPIs.put(file, ppis);
        }
        // Output ids and files
        System.out.println("IDs:");
        calculateOverlaps(fileToIds);
        System.out.println("PPIs:");
        calculateOverlaps(fileToPPIs);
    }

    private void calculateOverlaps(Map<String, Set<String>> fileToIds) {
        for (Iterator<String> it = fileToIds.keySet().iterator(); it.hasNext();) {
            String fileName = it.next();
            Set<String> ids = fileToIds.get(fileName);
            System.out.println(fileName + ": " + ids.size());
        }
        List<Set<String>> ids = new ArrayList<Set<String>>(fileToIds.values());
        for (int i = 0; i < ids.size() - 1; i++) {
            Set<String> ids1 = ids.get(i);
            String fileName = getKey(ids1, fileToIds);
            System.out.println(fileName);
            for (int j = i + 1; j < ids.size(); j++) {
                Set<String> ids2 = ids.get(j);
                String fileName2 = getKey(ids2, fileToIds);
                Set<String> tmp = new HashSet<String>(ids2);
                tmp.retainAll(ids1);
                System.out.println("\t" + fileName2 + ": " + tmp.size());
            }
        }
        // All three overlapps
        Set<String> tmp = new HashSet<String>(ids.get(0));
        for (int i = 1; i < ids.size(); i++) {
            Set<String> ids1 = ids.get(i);
            tmp.retainAll(ids1);
        }
        System.out.println("All overlapping: " + tmp.size());
    }
    
    private String getKey(Set<String> set, 
                          Map<String, Set<String>> map) {
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Set<String> value = map.get(key);
            if (value == set)
                return key;
        }
        return null;
    }
    
    /**
     * This method is used to count the total UniProt identifiers before the converted pathways
     * are stored into the database.
     * @throws Exception
     */
    @Test
    public void checkUniProtNumbersInConvertedDBs() throws Exception {
        // For a non-redundant count against SwissProt ids
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> swissProtIds = uniAnalyzer.loadSwissProtIDsMap();
        // Reactome database
        MySQLAdaptor reactomeDb = new MySQLAdaptor("localhost",
                                                   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                                   FIConfiguration.getConfiguration().get("DB_USER"),
                                                   FIConfiguration.getConfiguration().get("DB_PWD"),
                                                   3306);
        GKInstance human = reactomeDb.fetchInstance(48887L);
        // Want to have human proteins only
        Collection<?> ewases = reactomeDb.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                                    ReactomeJavaConstants.species, 
                                                                    "=",
                                                                    human);
        SchemaClass cls = reactomeDb.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.referenceEntity);
        reactomeDb.loadInstanceAttributeValues(ewases, att);
        Set<String> reactomeIds = new HashSet<String>();
        for (Iterator<?> it = ewases.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refEntity == null)
                continue;
            GKInstance species = (GKInstance) refEntity.getAttributeValue(ReactomeJavaConstants.species);
            if (species == human) {
                String identifier = (String) refEntity.getAttributeValue(ReactomeJavaConstants.identifier);
                reactomeIds.add(identifier);
            }
        }
        System.out.println("Total ids from Reactome: " + reactomeIds.size() + "\n");
        // Ids from KEGG
        //      String keggFileName = FIConfiguration.getConfiguration().get("KEGG_DIR + "kegg.rtpj";
        String keggFileName = FIConfiguration.getConfiguration().get("KEGG_CONVERTED_FILE");
        countUniProtIds(keggFileName, 
                        "KEGG", 
                        reactomeIds,
                        swissProtIds);
        // Ids from curated Nature-PID
        //    String pidFileName = FIConfiguration.getConfiguration().get("NATURE_PID_DIR + "NCI-Nature_Curated.rtpj";
        String pidFileName = FIConfiguration.getConfiguration().get("NATURE_PID_CURATED_CONVERTED");
        countUniProtIds(pidFileName, 
                        "Nature-PID", 
                        reactomeIds,
                        swissProtIds);
        
        // Ids from BioCarta from PID
        //    String biocartaFileName = FIConfiguration.getConfiguration().get("NATURE_PID_DIR + "BioCarta.rtpj";
        String biocartaFileName = FIConfiguration.getConfiguration().get("NATURE_PID_BIOCARTA_CONVERTED");
        countUniProtIds(biocartaFileName, 
                        "Biocarta-PID", 
                        reactomeIds,
                        swissProtIds);
        String tredFileName = FIConfiguration.getConfiguration().get("TRED_CONVERTED_FILE");
        countUniProtIds(tredFileName,
                        "TRED",
                        reactomeIds,
                        swissProtIds);      
        // Ids from Panther
        //        String pantherFileName = FIConfiguration.getConfiguration().get("PANTHER_DIR + "Panther_2_5.rtpj";
        String pantherFileName = FIConfiguration.getConfiguration().get("PANTHER_CONVERTED_FILE");
        countUniProtIds(pantherFileName, 
                        "Panther",
                        reactomeIds,
                        swissProtIds);
        
        // Check with PPIs
        Set<String> ppis = fu.loadInteractions(FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"));
        Set<String> totalIds = InteractionUtilities.grepIDsFromInteractions(ppis);
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtIdMap = uniProtAnalyzer.loadUniProtIDsMap();
        Map<String, String> swissProtIdMap = uniProtAnalyzer.loadSwissProtIDsMap();
        System.out.println("\nTotal ids from iRefIndex PPIs: " + totalIds.size());
        reactomeIds.addAll(totalIds);
        System.out.println("Total ids after merging: " + reactomeIds.size());
        reactomeIds.retainAll(swissProtIdMap.keySet());
        int totalSwissProtId = new HashSet<String>(swissProtIdMap.values()).size();
        double percentage = reactomeIds.size() / (double) totalSwissProtId;
        System.out.println("Total SwissProt ids after merging: " + reactomeIds.size() + 
                           " (" + percentage + ")");
        
        //        // Ids from CellMap
        //        String cellmapFileName = "datasets/cellmap_may_2006/CellMap.rtpj";
        //        countUniProtIds(cellmapFileName, 
        //                        "CellMap",
        //                        reactomeIds,
        //                        swissProtIds);     
        ////        // Ids from INOH
        ////        String inohFileName = FIConfiguration.getConfiguration().get("INOH_DIR + "INOH.rtpj";
        ////        countUniProtIds(inohFileName,
        ////                        "INOH",
        ////                        reactomeIds,
        ////                        swissProtIds);
        //        // Check with ids from GeneExpression
        //        String geneExpFile = "datasets/microarray/PrietoCarlos/union60InUniProt.txt";
        //        Set<String> idsFromExp = new FileUtility().loadInteractionTerms(geneExpFile, "\t");
        //        // Do a simple clean-up: remove spliceform
        //        Set<String> idCopy = new HashSet<String>();
        //        countUniProtIds(idsFromExp, 
        //                        "Array Data Set",
        //                        reactomeIds,
        //                        swissProtIds);
        //        // Check ids from IntAct proteins
        //        // Note: Some of ids in IntAct are from non-human species (e.g. Rat)
        //        String intActFileName = FIConfiguration.getConfiguration().get("INTACT_DIR + "IntAct.rtpj";
        //        countUniProtIds(intActFileName, 
        //                        "IntAct", 
        //                        reactomeIds, 
        //                        swissProtIds);
        //        String biogridFileName = FIConfiguration.getConfiguration().get("BIOGRID_DIR + "BioGrid.rtpj";
        //        countUniProtIds(biogridFileName,
        //                        "BioGrid",
        //                        reactomeIds,
        //                        swissProtIds);
        //        String hprdFileName = FIConfiguration.getConfiguration().get("HPRD_DIR + "HPRD.rtpj";
        //        countUniProtIds(hprdFileName,
        //                        "HPRD",
        //                        reactomeIds,
        //                        swissProtIds);
    }
    
    private void countUniProtIds(String projectFileName,
                                 String databaseName,
                                 Set<String> totalIds,
                                 Map<String, String> swissProtIdMap) throws Exception {
        Set<String> ids = loadUniProtIdsFromConvertedPathways(projectFileName);
        countUniProtIds(ids, databaseName, totalIds, swissProtIdMap);
    }


    private void countUniProtIds(Set<String> ids, 
                                 String databaseName,
                                 Set<String> totalIds,
                                 Map<String, String> swissProtIdMap) {
        System.out.println("Total ids from " + databaseName + ": "  + ids.size());
        // merge nature-pid ids to others
        // Want to get a percentage
        int originalSize = totalIds.size();
        totalIds.addAll(ids);
        int currentSize = totalIds.size();
        int increase = currentSize - originalSize;
        double percent = (double) increase / ids.size();
        System.out.println("Total ids after merging: " + currentSize + " (" + percent + ", " + increase + ")");
        // Count against SwissProt
        Set<String> swissProtIds = new HashSet<String>();
        for (String id : totalIds) {
            String mapped = swissProtIdMap.get(id);
            if (mapped != null)
                swissProtIds.add(mapped);
        }
        // Print out swissprot ids
        int totalSwissProtId = new HashSet<String>(swissProtIdMap.values()).size();
        double swissProtPercent = (double) swissProtIds.size() / totalSwissProtId;
        System.out.println("Total SwissProt ids after merging: " + swissProtIds.size() + " (" + swissProtPercent + ")");
        System.out.println();
    }
    
    private Set<String> loadUniProtIdsFromConvertedPathways(String fileName) throws Exception {
        Set<String> uniProtIds = new HashSet<String>();
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(fileName);
        // Some ReferenceGeneProduct instances don't have species, and some of them
        // are shell instances
        Collection<?> c = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferenceGeneProduct);
        for (Iterator<?> it = c.iterator(); it.hasNext();) {
            GKInstance inst = (GKInstance) it.next();
            // Just escape shell instances: there are not too many of them.
            if (inst.isShell())
                continue;
            String identifier = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier != null)
                uniProtIds.add(identifier);
        }
        // Make sure all UniProtIds are validated
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtIdMap = uniAnalyzer.loadUniProtIDsMap();
        Set<String> rtn = new HashSet<String>();
        for (String id : uniProtIds) {
            String tmp = uniProtIdMap.get(id);
            if (tmp != null)
                rtn.add(tmp);
        }
        return rtn;
    }
    
    /**
     * This method is used to count FIs from pathways in DB, and FIs predicted dumped in a
     * file.
     * @throws Exception
     */
    @Test
    public void countForDBPathwayFIsAndFilePredFIs() throws Exception {
//        List<ReactomeAnalyzer> analyzers = ReactomeAnalyzer.getPathwayDBAnalyzers();
//        ProteinIdFilters filters = new ProteinIdFilters();
//        Set<String> pathwayFIs = new HashSet<String>();
//        for (ReactomeAnalyzer analyzer : analyzers) {
//            Set<String> fis = analyzer.extractInteractionSet();
//            Set<String> filtered = filters.cleanUpVsUniProt(fis);
//            GKInstance dataSource = analyzer.getDataSource();
//            String source = null;
//            if (dataSource == null)
//                source = "Reactome";
//            else
//                source = dataSource.getDisplayName();
//            System.out.println("FIs from " + source + ": " + filtered.size());
//            pathwayFIs.addAll(filtered);
//        }
//        fu.saveInteractions(pathwayFIs, FIConfiguration.getConfiguration().get("RESULT_DIR + "PathwayFIs.txt");
        Set<String> pathwayFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + "PathwayFIs.txt");
        Set<String> pathwayIds = InteractionUtilities.grepIDsFromInteractions(pathwayFIs);
        System.out.printf("FIs from pathways: %d (%d)%n",
                          pathwayFIs.size(),
                          pathwayIds.size());
        countVsSwissProt(pathwayIds);
        //Set<String> predictedFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + "PredictedFIs.txt");
        Set<String> predictedFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + "PredictedFIsHumanPPIsGeneExp.txt");
        //Set<String> predictedFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR + "PredictedFIsCutoff09.txt");
        Set<String> predictedIDs = InteractionUtilities.grepIDsFromInteractions(predictedFIs);
        System.out.printf("FIs from prediction: %d (%d)%n",
                          predictedFIs.size(),
                          predictedIDs.size());
        countVsSwissProt(predictedIDs);
        // Merge them
        pathwayIds.addAll(predictedIDs);
        pathwayFIs.addAll(predictedFIs);
        System.out.printf("FIs merged: %d (%d)%n",
                          pathwayFIs.size(),
                          pathwayIds.size());
        countVsSwissProt(pathwayIds);
    }
    
    /**
     * This method is used to check total UniProt ids and SwissProt ids in a 
     * specified database. Alternative splicing forms are not counted.
     * @throws Exception
     */
    @Test
    public void checkUniProtIdsInDatabase() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_28_plus_i",
                                            "root",
                                            "macmysql01",
                                            3306);
        Collection refPepSeqs = dba.fetchInstancesByClass(ReactomeJavaConstants.ReferencePeptideSequence);
        SchemaClass refpepSeqcls = dba.getSchema().getClassByName(ReactomeJavaConstants.ReferencePeptideSequence);
        SchemaAttribute att = refpepSeqcls.getAttribute(ReactomeJavaConstants.identifier);
        dba.loadInstanceAttributeValues(refPepSeqs, att);
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> idMap = uniAnalyzer.loadUniProtIDsMap();
        Map<String, String> spMap = uniAnalyzer.loadSwissProtIDsMap();
        Set<String> uniIds = new HashSet<String>();
        Set<String> spIds = new HashSet<String>();
        for (Iterator it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = (GKInstance) it.next();
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            String mapped = idMap.get(identifier);
            if (mapped != null)
                uniIds.add(mapped);
            mapped = spMap.get(identifier);
            if (mapped != null)
                spIds.add(mapped);
        }
        double coverage = (double) uniIds.size() / new HashSet<String>(idMap.values()).size();
        System.out.println("Total UniProt Ids: " + uniIds.size() + " (" + coverage + ")");
        coverage = (double) spIds.size() / new HashSet<String>(spMap.values()).size();
        System.out.println("Total SwissProt Ids: " + spIds.size() + " (" + coverage + ")");
    }
    
    @Test
    public void countProteinsInChecksumFile() throws IOException {
        String dirName = FIConfiguration.getConfiguration().get("INTACT_DIR");
        File dir = new File(dirName);
        File[] fileList = dir.listFiles();
        for (File file : fileList) {
            String fileName = file.getName();
            if (fileName.endsWith("checksum.txt")) {
                Set<String> interactions = fu.loadInteractions(file.getAbsolutePath());
                Set<String> proteins = InteractionUtilities.grepIDsFromInteractions(interactions);
                System.out.println(fileName + ": " + proteins.size() + ", " + interactions.size());
            }
        }
    }
}
