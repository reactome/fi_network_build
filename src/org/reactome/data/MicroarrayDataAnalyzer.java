/*
 * Created on May 2, 2006
 *
 */
package org.reactome.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to analyze microarray data. This is a subclass of TestCase so that all
 * public methods can be invoked easily in Eclipse.
 * @author guanming
 *
 */
public class MicroarrayDataAnalyzer {
    private final String DATASET_DIR = "datasets/";
    private final String ARRAY_DIR = DATASET_DIR + "microarray/";
    private final String GPL_DIR = ARRAY_DIR + "GPL/";
    private final String GDS_DIR = ARRAY_DIR + "GDS/";
    private final String RESULT_DIR = "results/microarray/";
    
    private Map<String, String> id2ValueLine;
    private FileUtility fu;
    
    public MicroarrayDataAnalyzer() {  
        fu = new FileUtility();
    }
    
    /**
     * This method is used to check overlapping between two gene expresion 
     * data sets.
     * @throws IOException
     */
    @Test
    public void checkOverlapping() throws IOException {
        Set<String> prietoCarlosSet = loadCoExpFromPrietoCarlos();
        int pcSize = prietoCarlosSet.size();
        System.out.println("PrietoCarlos: " + prietoCarlosSet.size());
        Set<String> pavlidisSet = loadCoExpFromPavlidis();
        int pvSize = pavlidisSet.size();
        System.out.println("Pavlidis: " + pavlidisSet.size());
        prietoCarlosSet.retainAll(pavlidisSet);
        System.out.println("Shared pairs: " + prietoCarlosSet.size());
        System.out.println("Percentage: " + (double)prietoCarlosSet.size() / pcSize + 
                           ", " + (double)prietoCarlosSet.size() / pvSize);
    }
    
    /**
     * This method is used to check co-exp from Pavilidis.
     */
    @Test
    public void checkCoExpFromPavlidis() throws Exception {
        Set<String> coExp = loadCoExpFromPavlidis();
        System.out.println("Total Co-expression: " + coExp.size());
        FeatureChecker checker = new FeatureChecker();
        checker.checkFeatureOddsRatio(coExp);
    }

    /**
     * This method is used to load co-exp genes from Pavlidis dataset. Positive and
     * negative pairs have been merged in this file.
     * @return
     * @throws IOException
     */
    public Set<String> loadCoExpFromPavlidis() throws IOException {
        String fileName = FIConfiguration.getConfiguration().get("LEE_GENE_EXP_FILE");
        return fu.loadInteractions(fileName);
    }
    
    @Test
    public void replaceSpaceWithTab() throws IOException {
//        fu.setInput(FIConfiguration.getConfiguration().get("LEE_GENE_EXP_FILE);
        fu.setInput("results/v3/PavlidisCoExp_Norm.txt");
        String line = null;
        String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "LeeGeneExp.txt";
        fu.setOutput(outFileName);
        while ((line = fu.readLine()) != null) {
            line = line.replace(" ", "\t");
            fu.printLine(line);
        }
        fu.close();
//        fu.setInput(FIConfiguration.getConfiguration().get("PRIETO_GENE_EXP_FILE);
        fu.setInput("results/v3/CarlosCoExp_Norm.txt");
        outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "PrietoGeneExp.txt";
        fu.setOutput(outFileName);
        while ((line = fu.readLine()) != null) {
            line = line.replace(" ", "\t");
            fu.printLine(line);
        }
        fu.close();
    }
    
    /**
     * This method is used to check the utility of gene expression feature in 
     * NBC training based on odds ratio.
     * @throws Exception
     */
    @Test
    public void checkCoExpFromPrietoCarlos() throws Exception {
        Set<String> coExpGenes = loadCoExpFromPrietoCarlos();
        System.out.println("Co-expressed genes: " + coExpGenes.size());
        FeatureChecker checker = new FeatureChecker();
        checker.checkFeatureOddsRatio(coExpGenes);
    }
    
    /**
     * Load the co-expressed protein pairs from Prieto Caralos's paper.
     * @return
     * @throws IOException
     */
    public Set<String> loadCoExpFromPrietoCarlos() throws IOException {
        //String fileName = ARRAY_DIR + "PrietoCarlos/union60InUniProt.txt";
//        String fileName = "results/v3/CarlosCoExp_Norm.txt";
        String fileName = FIConfiguration.getConfiguration().get("PRIETO_GENE_EXP_FILE");
        return fu.loadInteractions(fileName);
    }
    
    @Test
    public void normalizeLeeGeneExp() throws Exception {
        // This file contains UniProt pairs mapped from link-data.txt using refseq-hs-annots.txt
        // All pairs have been supported by at least 3 experiments.
        String inFileName = FIConfiguration.getConfiguration().get("LEE_GENE_EXP_FILE_SOURCE");
//        System.out.println("inFileName: " + inFileName);
        fu.setInput(inFileName);
        String line = null;
        Set<String> pairs = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("( |\t)");
            pairs.add(tokens[1] + "\t" + tokens[2]);
        }
        fu.close();
        ProteinIdFilters filters = new ProteinIdFilters();
        System.out.println("Pairs before normalizing: " + pairs.size());
        Set<String> normalized = filters.normalizeProteinPairs(pairs);
        System.out.println("Pairs after normalizing: " + normalized.size());
        fu.saveInteractions(normalized, FIConfiguration.getConfiguration().get("LEE_GENE_EXP_FILE"));
    }
    
    /**
     * This method is used to normalize two co-exp file.
     * @throws Exception
     */
    @Test
    public void normalizeCoExpPairs() throws Exception {
        ProteinIdFilters filters = new ProteinIdFilters();
        Set<String> pairs = loadCoExpFromPavlidis();
        System.out.println("Before normalization: " + pairs.size());
        Set<String> normalized = filters.normalizeProteinPairs(pairs);
        System.out.println("After: " + normalized.size());
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "PavlidisCoExp_Norm.txt";
//        fu.outputSet(normalized, fileName);
        pairs = loadCoExpFromPrietoCarlos();
        System.out.println("Before normalization: " + pairs.size());
        normalized = filters.normalizeProteinPairs(pairs);
        System.out.println("After: " + normalized.size());
//        fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "CarlosCoExp_Norm.txt";
//        fu.outputSet(normalized, fileName);
    }
    
    /**
     * Use this method to generate a UniProt pairs based on the download gene pairs
     * from the Prieto Carlos expression data file.
     * @throws IOException
     */
    @Test
    public void generatePrietoCarlosGeneExpFile() throws Exception {
        String srcFileName = "datasets/microarray/PrietoCarlos/union60.txt";
        Set<String> genePairs = fu.loadInteractions(srcFileName);
        System.out.println("Total gene pairs: " + genePairs.size());
        long time1 = System.currentTimeMillis();
        Map<String, Set<String>> geneNameToUniAcces = new UniProtAnalyzer().generateGeneNameToUniAccesssMap(true);
//        Map<String, Set<String>> geneNameToUniAcces = new UniProtAnalyzer().loadGeneNameToUniProtAcces();
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for loading: " + (time2 - time1));
        System.out.println("Size of geneNameToUnAcces: " + geneNameToUniAcces.size());
        Set<String> uniPairs = new HashSet<String>();
        int unmapped = 0;
        for (String genePair : genePairs) {
            String[] tokens = genePair.split("\t");
            Set<String> access1 = geneNameToUniAcces.get(tokens[0]);
            Set<String> access2 = geneNameToUniAcces.get(tokens[1]);
            // It is possible that no mapping is available
            if (access1 == null || access2 == null) {
                unmapped ++;
                continue;
            }
            for (String acc1 : access1) {
                for (String acc2 : access2) {
                    if (acc1.equals(acc2))
                        continue; // Don't want to have self interaction
                    if (acc1.compareTo(acc2) < 0)
                        uniPairs.add(acc1 + "\t" + acc2);
                    else
                        uniPairs.add(acc2 + "\t" + acc1);
                }
            }
//            System.out.println(genePair + " -> " + uniPairs.size());
        }
        System.out.println("Total unmapped: " + unmapped);
        System.out.println("Total UniProt pairs before filtering: " + uniPairs.size());
        ProteinIdFilters filters = new ProteinIdFilters();
        Set<String> normalized = filters.normalizeProteinPairs(uniPairs);
        System.out.println("After filtering: " + normalized.size());
        fu.saveInteractions(normalized, 
                            FIConfiguration.getConfiguration().get("PRIETO_GENE_EXP_FILE"));
    }
    
    /**
     * This method is used to convert a gene pair to a UniProt pair for a specific file.
     * @throws IOException
     */
    @Test
    @Deprecated
//    public void convertGeneNameToUniProtIdentifiers() throws Exception {
//        String dirName = "datasets/microarray/PrietoCarlos/";
//        String inFileName = dirName + "union60.txt";
//        String outFileName = dirName + "union60InUniProt.txt";
//        FileUtility fu = new FileUtility();
//        Set<String> fisInNames = fu.loadInteractions(inFileName);
//        System.out.println("Total Fis in gene names: " + fisInNames.size());
//        Set<String> names = InteractionUtilities.grepIDsFromInteractions(fisInNames);
//        UCSCDataAnalyzer ucscAnalyzer = new UCSCDataAnalyzer();
//        // Use the following method, one gene name is mapped to one UniProt name only.
//        // However, this may not be the case.
//        Map<String, String> nameToUniIds = ucscAnalyzer.mapToUniProt(names);
//        Set<String> fisInUniIds = new HashSet<String>();
//        int index = 0;
//        String name1 = null;
//        String name2 = null;
//        // Check how many names cannot be mapped
//        Set<String> notMappedNames = new HashSet<String>();
//        for (String fi : fisInNames) {
//            index = fi.indexOf("\t");
//            name1 = fi.substring(0, index);
//            name2 = fi.substring(index + 1);
//            String id1 = nameToUniIds.get(name1);
//            if (id1 == null) {
//                notMappedNames.add(name1);
//                continue;
//            }
//            String id2 = nameToUniIds.get(name2);
//            if (id2 == null) {
//                notMappedNames.add(name2);
//                continue;
//            }
//            // Do a permutation
//            
//            int compare = id1.compareTo(id2);
//            if (compare < 0) {
//                fisInUniIds.add(id1 + " " + id2);
//            }
//            else if (compare > 0) {
//                fisInUniIds.add(id2 + " " + id1);
//            }
//        }
//        System.out.println("Total FIs in UniProt ids: " + fisInUniIds.size());
//        fu.outputSet(fisInUniIds, outFileName);
//        System.out.println("Names cannot be mapped: " + notMappedNames.size());
//        for (String name : notMappedNames)
//            System.out.println(name);
//    }
    
    public void rnaGiToUniAccMapGenerate() throws IOException {
        String gene2acc = DATASET_DIR + "ncbi/gene2accession";
        String uniprot2gi = DATASET_DIR + "iproclass/iproclass.tb";
        FileUtility fu = new FileUtility();
        fu.setInput(gene2acc);
        String line = null;
        Map<String, String> protGiTornaGi = new HashMap<String, String>();
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (!tokens[0].equals("9606"))
                continue; // Only human is needed
            if (tokens[4].equals("-") ||
                tokens[6].equals("-"))
                continue; // No values
            protGiTornaGi.put(tokens[6], tokens[4]);
        }
        System.out.println("Total RNA to Protein GI Map: " + protGiTornaGi.size());
        fu.close();
        fu.setInput(uniprot2gi);
        Map<String, String> rnaToUniProt = new HashMap<String, String>();
        String[] gis = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            if (tokens[4].length() == 0)
                continue; // Nothing
            gis = tokens[4].split("; ");
            // Check if gis is in the first map
            for (String gi : gis) {
                if (protGiTornaGi.containsKey(gi)) {
                    rnaToUniProt.put(protGiTornaGi.get(gi), tokens[0]);
                }
            }
        }
        fu.close();
        System.out.println("Total RNA to UniProt: " + rnaToUniProt.size());
        // Output the map
        String outFile = "results/rnaGi2UniProt.txt";
        fu.exportMap(rnaToUniProt, outFile);
    }
    
    public void processGPLForFiles() throws IOException {
        File[] gplFiles = new File(GPL_DIR).listFiles();
        FileUtility fu = new FileUtility();
        Map<String, String> gi2uni = fu.importMap("results/rnaGi2UniProt.txt");
        int total = 0;
        int uni = 0;
        Map<String, String> id2uni = new HashMap<String, String>();
        String[] tokens = null;
        String line = null;
        for (File file : gplFiles) {
            fu.setInput(file.getAbsolutePath());
            total = 0;
            uni = 0;
            id2uni.clear();
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("#") ||
                    line.startsWith("!") ||
                    line.startsWith("^")) 
                    continue; // All are comment lines
                total ++;
                tokens = line.split("\t");
                // Get gi
                if (tokens.length > 7 && tokens[6].length() > 0) {
                    String uniAc = gi2uni.get(tokens[6]);
                    if (uniAc != null) {
                        uni ++;
                        id2uni.put(tokens[0], uniAc);
                    }
                }
            }
            fu.close();
            System.out.printf("%s (UniProt): %d (%f)%n", 
                              file.getName(), 
                              uni, 
                              uni / (float)total);
        }
    }
    
    public void processGDS() throws IOException {
        String gdsPickup = RESULT_DIR + "PickupGDSSimple.txt";
        Map<String, String> gds2gpl = new HashMap<String, String>();
        List<String> gdsList = new ArrayList<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(gdsPickup);
        String line = fu.readLine(); // escape the first line
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Take the second 
            gds2gpl.put(tokens[0], tokens[2]);
            gdsList.add(tokens[0]);
        }
        fu.close();
        // Load all map first
        Map<String, Map<String, String>> gplMap = loadGPLMap(gds2gpl.values());
        // Start process
        int sampleSize = 0;
        StringBuilder builder = new StringBuilder();
        FileUtility outFu = new FileUtility();
        for (String gds : gdsList) {
            String gdsFile = GDS_DIR + gds + ".soft";
            String outFile = RESULT_DIR + gds + ".txt";
            outFu.setOutput(outFile);
            Map<String, String> id2uni = gplMap.get(gds2gpl.get(gds));
            if (id2uni == null)
                throw new IllegalStateException("No map found for '" + gds + "'");
            fu.setInput(gdsFile);
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("#") ||
                    line.startsWith("!") ||
                    line.startsWith("^")) 
                        continue; // All are comment lines
                if (line.startsWith("ID_REF")) {
                    tokens = line.split("\t");
                    sampleSize = tokens.length - 2; // Don't count ID_REF and IDENTIFIER
                    continue;
                }
                tokens = line.split("\t");
                if (tokens.length - 1 < sampleSize)
                    throw new IllegalStateException("Sample size is not right: " + line + "(" + gds + ")");
                // The first token is probe_id
                String uni = id2uni.get(tokens[0]);
                if (uni == null)
                    continue;
                // Prepare output
                builder.setLength(0);
                builder.append(uni);
                for (int i = 2; i < tokens.length; i++)
                    builder.append("\t").append(tokens[i]);
                outFu.printLine(builder.toString());
            }
            outFu.close();
            fu.close();
            System.out.println("Done..." + gds);
        }
        averageGDS();
    }
    
    public void averageGDS() throws IOException {
        String gdsPickup = RESULT_DIR + "PickupGDSSimple.txt";
        List<String> gdsList = new ArrayList<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(gdsPickup);
        String line = fu.readLine(); // escape the first line
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            gdsList.add(line.substring(0, index));
        }
        fu.close();
        String id = null;
        List<String> lines = null;
        for (String gds : gdsList) {
            // Used to hold data line by line keyed by UniProt IDs
            Map<String, List<String>> id2Lines = new HashMap<String, List<String>>();
            String fileName = RESULT_DIR + gds + ".txt";
            fu.setInput(fileName);
            while ((line = fu.readLine()) != null) {
                index = line.indexOf("\t");
                id = line.substring(0, index);
                lines = id2Lines.get(id);
                if (lines == null) {
                    lines = new ArrayList<String>();
                    id2Lines.put(id, lines);
                }
                lines.add(line);
            }
            fu.close();
            // Check if average is needed
            average(gds, id2Lines);
            System.out.println("Done..." + gds);
            // Delete the original file
            File file = new File(fileName);
            file.delete();
        }
    }
    
    private void average(String gds, Map<String, List<String>> id2Lines) throws IOException {
        String outName = RESULT_DIR + gds + ".ave.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(outName);
        for (Iterator<String> it = id2Lines.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            List<String> lines = id2Lines.get(id);
            if (lines.size() == 1) {
                fu.printLine(lines.get(0));
            }
            else {
                fu.printLine(average(id, lines));
            }
        }
        fu.close();
    }
    
    private String average(String id, List<String> lines) throws IOException {
        StringBuilder builder = new StringBuilder();
        int size = lines.size();
        // Generate data set
        Float[][] values = new Float[size][];
        int firstIndex = 0;
        for (String line : lines) {
            String[] tokens = line.split("\t");
            values[firstIndex] = new Float[tokens.length - 1]; // Escape the first column
            for (int i = 1; i < tokens.length; i ++) {
                if (tokens[i].equals("NULL"))
                    continue; // Just NULL
                values[firstIndex][i - 1] = new Float(tokens[i]);
            }
            firstIndex ++;
        }
        // Average
        builder.append(id);
        int sampleSize = values[0].length;
        float total = 0.0f;
        int sizeTotal = 0;
        for (int i = 0; i < sampleSize; i++) {
            total = 0.0f;
            sizeTotal = 0;
            for (int j = 0; j < size; j++) {
                if (values[j][i] != null) {
                    sizeTotal ++;
                    total += values[j][i];
                }
            }
            if (sizeTotal > 0)
                builder.append("\t").append(total / sizeTotal);
            else
                builder.append("\t").append("NULL");
        }
        return builder.toString();
    }
    
    private Map<String, Map<String, String>> loadGPLMap(Collection<String> gpls) throws IOException {
        Map<String, Map<String, String>> gplMap = new HashMap<String, Map<String, String>>();
        FileUtility fu = new FileUtility();
        for (String gpl : gpls) {
            String fileName = RESULT_DIR + gpl + ".map";
            Map<String, String> map = fu.importMap(fileName);
            gplMap.put(gpl, map);
        }
        return gplMap;
    }
    
    public void processGPL() throws IOException {
        String gdsPickup = RESULT_DIR + "PickupGDSSimple.txt";
        Set<String> platforms = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(gdsPickup);
        String line = fu.readLine(); // escape the first line
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Take the second 
            platforms.add(tokens[2]);
        }
        fu.close();
        System.out.printf("Platforms: %d %s%n", platforms.size(), platforms);
        // Load map from rna gi to uniprot ac
        Map<String, String> gi2uni = fu.importMap("results/rnaGi2UniProt.txt");
        // Check the total
        int total = 0;
        int uni = 0;
        Map<String, String> id2uni = new HashMap<String, String>();
        // Test
        //platforms.clear();
        //platforms.add("GPL201");
        int commaIndex;
        for (String gpl : platforms) {
            if (gpl.equals("GPL1074")) {
                id2uni.clear();
                handleGPL1074(id2uni);
                continue;
            }
            fu.setInput(GPL_DIR + gpl + ".annot");
            total = 0;
            uni = 0;
            id2uni.clear();
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("#") ||
                    line.startsWith("!") ||
                    line.startsWith("^")) 
                    continue; // All are comment lines
                total ++;
                tokens = line.split("\t");
                // Get gi
                if (tokens.length > 7 && tokens[6].length() > 0) {
                    // might contain ",". If it mapped to multiple GI, exclude it since
                    // it makes values for uniprot unreliable.
                    commaIndex = tokens[6].indexOf(",");
                    if (commaIndex > 0)
                        continue; // escape
                    String uniAc = gi2uni.get(tokens[6]);
                    if (uniAc != null) {
                        uni ++;
                        id2uni.put(tokens[0], uniAc);
                    }
                }
            }
            fu.close();
            System.out.printf("%s (UniProt): %d (%f)%n", 
                              gpl, 
                              uni, 
                              uni / (float)total);
            if (id2uni.size() > 0)
                fu.exportMap(id2uni, RESULT_DIR + gpl + ".map");
        }
    }
    
    private void handleGPL1074(Map<String, String> map) throws IOException {
        String file = DATASET_DIR + "microarray/norvartis/gnf1b-anntable.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(file);
        String line = fu.readLine(); // Escape the first line
        String[] tokens = null;
        int total = 0;
        int hit = 0;
        while ((line = fu.readLine()) != null) {
            total ++;
            tokens = line.split("\t");
            if (tokens.length > 10 &&
                tokens[9].length() > 0 &&
                tokens[9].indexOf(";") < 0) {
                hit ++;
                map.put(tokens[3], tokens[9]);
            }
        }
        System.out.printf("GPL1074 (uniprot): %d (%f)%n",
                          hit,
                          hit / (float)total);
        fu.exportMap(map, RESULT_DIR + "GPL1074.map");
    }
    
    public void generateInteractionIndex() throws IOException {
        String fileName = "results/GeneExpFromPavlidis.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int c = 0;
        Set<String> ids = new HashSet<String>();
        long time1 = System.currentTimeMillis();
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String sub = line.substring(index + 1);
            index = sub.indexOf(" ");
            ids.add(sub.substring(0, index));
            ids.add(sub.substring(index + 1));
        }
        fu.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for loading: " + (time2 - time1));
        List<String> idList = new ArrayList<String>(ids);
        Collections.sort(idList);
        long time3 = System.currentTimeMillis();
        System.out.println("Time for sorting: " + (time3 - time2));
        System.out.println("Total Memory: " + Runtime.getRuntime().totalMemory());
        fu = new FileUtility();
        fu.setOutput("results/GeneExpIDs.txt");
        for (String id : idList)
            fu.printLine(id);
        fu.close();
        Map<Integer, List<Integer>> posMap = new HashMap<Integer, List<Integer>>();
        Map<Integer, List<Integer>> negMap = new HashMap<Integer, List<Integer>>();
        fu = new FileUtility();
        fu.setInput(fileName);
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String value = line.substring(0, index);
            String sub = line.substring(index + 1);
            index = sub.indexOf(" ");
            String id1 = sub.substring(0, index);
            String id2 = sub.substring(index + 1);
            int partner1 = Collections.binarySearch(idList, id1);
            int partner2 = Collections.binarySearch(idList, id2);
            if (value.equals("+")) {
                // Positive
                List<Integer> list = posMap.get(partner1);
                if (list == null) {
                    list = new ArrayList<Integer>();
                    posMap.put(partner1, list);
                }
                list.add(partner2);
            }
            else {
                // Negative
                List<Integer> list = negMap.get(partner1);
                if (list == null) {
                    list = new ArrayList<Integer>();
                    negMap.put(partner1, list);
                }
                list.add(partner2);
            }
        }
        fu.close();
        long time4 = System.currentTimeMillis();
        System.out.println("Time for creatimg map: " + (time4 - time3));
        System.out.println("Total Memory after creating map: " + Runtime.getRuntime().totalMemory());
        // Output interaction index
        fu.setOutput("results/GeneExpPosCol.txt");
        outputIndexMap("results/GeneExpPosCol.txt", posMap);
        outputIndexMap("results/GeneExpNegCol.txt", negMap);
    }
    
    private void outputIndexMap(String fileName, Map<Integer, List<Integer>> valueMap) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        StringBuilder builder = new StringBuilder();
        for (Iterator<Integer> it = valueMap.keySet().iterator(); it.hasNext();) {
            Integer index1 = it.next();
            List<Integer> list = valueMap.get(index1);
            builder.setLength(0);
            builder.append(index1).append(" ");
            for (Integer i : list)
                builder.append(i).append(" ");
            fu.printLine(builder.toString());
        }
        fu.close();
    }
    
    public void countExpForLinkData() throws IOException {
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/microarray/Pavlidis/link-data.txt";
        Map<Integer, Integer> countMap = new HashMap<Integer, Integer>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String[] tokens = null;
        String line;
        line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            Integer index = Integer.parseInt(tokens[1]);
            if (countMap.get(index) == null)
                countMap.put(index, 1);
            else {
                countMap.put(index, countMap.get(index) + 1);
            }
        }
        for (Iterator<Integer> it = countMap.keySet().iterator(); it.hasNext();) {
            Integer key = it.next();
            Integer value = countMap.get(key);
            System.out.println(key + " " + value);
        }
        fu.close();
    }
    
    public void processGeneExpLinkData() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/microarray/Pavlidis/";
        String mapFileName = dirName + "refseq-hs-annots.txt";
        String srcFileName = dirName + "link-data.txt";
        // Create a map from ids used in src to uniprot accession numbers
        Map<String, String> idToUni = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        fu.setInput(mapFileName);
        String[] tokens = null;
        int c = 0;
        // Escape the first line
        String line = fu.readLine();
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Take the first one as id and the third one as SwissProt id
            if (tokens[2].length() == 0) {
                c++;
                continue;
            }
            idToUni.put(tokens[0], tokens[2]);
        }
        fu.close();
        System.out.println("Not Mapped: " + c);
        //String outputFileName = "results/GeneExpFromPavlidis.txt";
        // This file contains links confirmed by 3 or more than 3 experiments
        String outputFileName = "results/microarray/GeneExpWith3FromPavlidis.txt";
        FileUtility output = new FileUtility();
        output.setOutput(outputFileName);
        fu.setInput(srcFileName);
        line = fu.readLine();
        String uni1, uni2;
        int comp = 0;
        c = 0;
        int expNumber;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            expNumber = Integer.parseInt(tokens[1]);
            if (expNumber < 3)
                continue; // Take links confirmed by more than 3
            uni1 = idToUni.get(tokens[2]);
            uni2 = idToUni.get(tokens[3]);
            if (uni1 == null || uni2 == null)
                continue; // Make sure both of them are present
            comp = uni1.compareTo(uni2);
            if (comp < 0)
                output.printLine(tokens[0] + "\t" + uni1 + " " + uni2);
            else if (comp > 0)
                output.printLine(tokens[0] + "\t" + uni2 + " " + uni1);
            c ++;
        }
        fu.close();
        output.close();
        System.out.println("Total pair: " + c);
    }
    
    /**
     * Use this method to load a text file into MySQL database
     * @throws Exception
     */
    public void loadDataIntoDB() throws Exception {
        Class.forName("com.mysql.jdbc.Driver");
        Connection connection = DriverManager.getConnection("jdbc:mysql://localhost/BNDataSource?" +
                "user=root&password=macmysql01");
        //String insertSql = "INSERT INTO GeneExpLink (uniprot1, uniprot2, value) VALUES (?, ?, ?)";
        String insertSql = "INSERT INTO GeneExpPair (pair, value) VALUES (?, ?)";
        PreparedStatement stat = connection.prepareStatement(insertSql);
        FileUtility source = new FileUtility();
        source.setInput("results/GeneExpFromPavlidis.txt");
        String line = null;
        String[] tokens;
        while ((line = source.readLine()) != null) {
            tokens = line.split("\t");
            stat.setString(1, tokens[1]);
            stat.setString(2, tokens[0]);
            //stat.setString(2, tokens[2]);
            //stat.setString(3, tokens[0]);
            stat.execute();
        }
        stat.close();
        connection.close();
    }
    
    /**
     * Generate a map file from EMBL ids to UniProt accession numbers.
     * @throws IOException
     */
    public void generateEMBLToUniMap() throws IOException {
        long time1 = System.currentTimeMillis();
        // Only SwissProt part is handled for the pilot project
        String uniFileName = "/Users/wgm/Documents/caBIG_R3/datasets/UniProt/uniprot_sprot_human.dat";
        FileReader fileReader = new FileReader(uniFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        // For output
        String outputFileName = "results/embl2uni.txt";
        FileWriter fileWriter = new FileWriter(outputFileName);
        PrintWriter printWriter = new PrintWriter(fileWriter);
        String line = null;
        String ac = null;
        String emblString = null;
        String[] emblIds = null;
        int index = 0;
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("AC")) {
                index = line.indexOf(";");
                ac = line.substring(5, index);
            }
            else if (line.startsWith("DR   EMBL")) {
                emblString = line.substring(11); // Escape "; " after EMBL
                emblIds = emblString.split("; ");
                // Sometimes there is no "-"
                int c = emblIds.length;
                int lastC = 0;
                // Have to find "-" since more tokens might exist
                for (int i = 0; i < c; i++) {
                    if (emblIds[i].equals("-") ||
                        emblIds[i].equals("JOINED")) {
                        lastC = i;
                        break;
                    }
                }
                if (lastC == 0)
                    lastC = c - 1;
                for (int i = 0; i < lastC; i++) {
                    printWriter.println(emblIds[i] + "\t" + ac);
                }
            }
        }
        // Need to close all writers and readers.
        reader.close();
        fileReader.close();
        printWriter.close();
        fileWriter.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateEMBLToUniMap: " + (time2 - time1));
    } 
    
    public void pickupGDSFiles() throws IOException {
        // GPL in this list should be escaped
        String[] gpls = new String[]{
                "GPL781","GPL169","GPL178","GPL181","GPL182",
                "GPL183","GPL4"
        };
        List<String> escapedList = Arrays.asList(gpls);
        String inFile = RESULT_DIR + "PickupGDS.txt";
        String outFile = RESULT_DIR + "PickupGDSSimple.txt";
        FileUtility inFu = new FileUtility();
        inFu.setInput(inFile);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFile);
        String line = null;
        StringBuilder builder = new StringBuilder();
        outFu.printLine("Accession\tSummary\tParent Platform\t" +
                "Reference Series\tType\tSubsets\t" +
        "Subset Total\tSamples");
        Pattern pattern = Pattern.compile("^(\\d+): (GDS\\d+)");
        Pattern numberPattern = Pattern.compile("(\\d+)");
        int index = 0;
        String platform = null;
        while ((line = inFu.readLine()) != null) {
            Matcher matcher = pattern.matcher(line);
            if (matcher.find()) {
                if (builder.length() > 0) {
                    if (!escapedList.contains(platform)) 
                        outFu.printLine(builder.toString());
                    //System.out.println(builder.toString());
                    builder.setLength(0);
                }
                //System.out.println(line);
                builder.append(matcher.group(2));
            }
            if (line.startsWith("Summary:")) {
                index = line.indexOf(":");
                builder.append("\t").append(line.substring(index + 2));
                while (true) {
                    line = inFu.readLine();
                    if (line.startsWith("Parent Platform"))
                        break;
                    builder.append(" ").append(line);
                }
                index = line.indexOf(":");
                platform = line.substring(index + 2);
                builder.append("\t").append(platform);
            }
            if (line.startsWith("Reference Series") ||
                line.startsWith("Type")) {
                index = line.indexOf(":");
                builder.append("\t").append(line.substring(index + 2));
            }
            if (line.startsWith("Subsets")) {
                index = line.indexOf(":");
                builder.append("\t").append(line.substring(index + 2));
                // Extract the total number of subsets
                matcher = numberPattern.matcher(line);
                int start = 0;
                int total = 0;
                while (matcher.find(start)) {
                    String num = matcher.group(1);
                    total += Integer.parseInt(num);
                    start = matcher.end();
                }
                builder.append("\t").append(total);
            }
            if (line.startsWith("Samples")) {
                // Find the first number
                matcher = numberPattern.matcher(line);
                if (matcher.find()) {
                    builder.append("\t").append(matcher.group(1));
                }
            }
        }
        // might be lost
        if (builder.length() > 0 && !escapedList.contains(platform))
            outFu.printLine(builder.toString());
        inFu.close();
        outFu.close();
    }
    
    public void convertData() throws IOException {
        File dir = new File(RESULT_DIR);
        File[] files = dir.listFiles();
        List<String[]> values = new ArrayList<String[]>();
        FileUtility fu = new FileUtility();
        String line = null;
        String[] tokens = null;
        for (File file : files) {
            String fileName = file.getName();
            if (fileName.endsWith(".ave.txt")) {
                fu.setInput(file.getAbsolutePath());
                while ((line = fu.readLine()) != null) {
                    tokens = line.split("\t");
                    String[] tmp = new String[tokens.length - 1];
                    System.arraycopy(tokens, 1, tmp, 0, tokens.length - 1);
                    values.add(tmp);
                }
                fu.close();
                String outFile = RESULT_DIR + fileName + ".val";
                fu.setOutput(outFile);
                fu.printLine("Value\tSample");
                for (String[] tmp : values) {
                    for (int i = 0; i < tmp.length; i++) {
                        fu.printLine(tmp[i] + "\tV" + i);
                    }
                }
                fu.close();
                values.clear();
            }
        }
    }
    
    public void countAllData() throws IOException {
        String file = RESULT_DIR + "AllData.txt";
        int c = 0;
        String line = null;
        FileUtility fu = new FileUtility();
        fu.setInput(file);
        int total = 1600;
        int cc20 = 0;
        int cc10 = 0;
        int cc15 = 0;
        int totalUni = 0;
        String outFile = RESULT_DIR + "AllData15.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFile);
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            c = -1;
            for (String token : tokens) {
                if (token.equals("NULL"))
                    c++;
            }
            float ratio = c / (float) total;
            if (ratio < 0.2f)
                cc20++;
            if (ratio < 0.15f)
                cc15++;
            else
                outFu.printLine(line);
            if (ratio < 0.1f)
                cc10 ++;
            totalUni ++;
        }
        fu.close();
        outFu.close();
        System.out.println("Total UniProt: " + totalUni);
        System.out.println("Total Sample: " + total);
        System.out.println("Less than 20%: " + cc20);
        System.out.println("Less than 15%: " + cc15);
        System.out.println("Less than 10%: " + cc10);
    }
    
    public void aggregate() throws IOException {
        String outputFile = RESULT_DIR + "AllData.txt";
        List<String> gdsList = getGDSList();
        Map<String, Integer> sampleSize = getSampleSize();
        Map<String, String> allData = new HashMap<String, String>();
        String fileName = null;
        FileUtility fu = new FileUtility();
        String line = null;
        int index = 0;
        Map<String, String> data = new HashMap<String, String>();
        int sizeInAll = 0;
        int c = 1;
        for (String gds : gdsList) {
            fileName = RESULT_DIR + gds + ".ave.txt";
            fu.setInput(fileName);
            data.clear();
            while ((line = fu.readLine()) != null) {
                index = line.indexOf("\t");
                String id = line.substring(0, index);
                data.put(id, line.substring(index + 1));
            }
            Set<String> notMerged = new HashSet<String>(allData.keySet());
            notMerged.removeAll(data.keySet());
            // postLine is needed usually
            String postLine = "";
            int sSize = sampleSize.get(gds);
            for(int i = 0; i < sSize; i++) {
                postLine += "NULL";
                if (i < sSize - 1)
                    postLine += "\t";
            }
            String preLine = null;
            for (Iterator<String> it = data.keySet().iterator(); it.hasNext();) {
                String id = it.next();
                line = data.get(id);
                String dataLine = allData.get(id);
                if (dataLine == null) {
                    // Need to generate preLine
                    if (preLine == null) {
                        preLine = "";
                        for (int i = 0; i < sizeInAll; i++) {
                            preLine += "NULL";
                            if (i < sizeInAll - 1)
                                preLine += "\t";
                        }
                    }
                    dataLine = preLine;
                }
                if (dataLine.length() > 0)
                    dataLine += ("\t" + line);
                else
                    dataLine = line; // Work for the first file
                // Need to reinsert
                allData.put(id, dataLine);
            }
            for (String id : notMerged) {
                String dataLine = allData.get(id);
                allData.put(id, dataLine + "\t" + postLine);
            }
            fu.close();
            sizeInAll += sSize;
            System.out.println("Done..." + gds + " " + c++);
        }
        fu.exportMap(allData, outputFile);
    }
    
    private Map<String, Integer> getSampleSize() throws IOException {
        String gdsPickup = RESULT_DIR + "PickupGDSSimple.txt";
        Map<String, Integer> sampleSize = new HashMap<String, Integer>();
        FileUtility fu = new FileUtility();
        fu.setInput(gdsPickup);
        String line = fu.readLine(); // escape the first line
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            sampleSize.put(tokens[0], Integer.parseInt(tokens[7]));
        }
        fu.close();
        return sampleSize;
    }
    
    public List<String> getGDSList() throws IOException {
        String gdsPickup = RESULT_DIR + "PickupGDSSimple.txt";
        List<String> gdsList = new ArrayList<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(gdsPickup);
        String line = fu.readLine(); // escape the first line
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Take count only
            int index = tokens[4].indexOf(",");
            String count = tokens[4].substring(index + 1).trim();
            if (count.equals("count"))
                gdsList.add(tokens[0]);
        }
        fu.close();
        return gdsList;
    }
    
    public void loadData(String fileName) throws IOException {
        id2ValueLine = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id2ValueLine.put(line.substring(0, index), line.substring(index + 1));
        }        
    }
    
    public Float calculateCorrelation(String id1, String id2) {
        if (id2ValueLine == null)
            throw new IllegalStateException("Id2 to value map has not been loaded yet.");
        String line1 = id2ValueLine.get(id1);
        if (line1 == null)
            return null;
        String line2 = id2ValueLine.get(id2);
        if (line2 == null)
            return null;
        Float[] values1 = parseValueLine(line1);
        Float[] values2 = parseValueLine(line2);
        // Make sure two arrays has the same pairwise values
        for (int i = 0; i < values1.length; i++) {
            if (values1[i] == null && values2[i] != null)
                values2[i] = null;
            if (values1[i] != null && values2[i] == null)
                values1[i] = null;
        }
        List<Float> valuesList1 = new ArrayList<Float>();
        for (Float f : values1) {
            if (f != null)
                valuesList1.add(f);
        }
        if (valuesList1.size() == 0)
            return null; // Cannot calculate
        List<Float> valuesList2 = new ArrayList<Float>();
        for (Float f : values2) {
            if (f != null)
                valuesList2.add(f);
        }
        // Get the total
        float mean1 = mean(valuesList1);
        float mean2 = mean(valuesList2);
        float num = 0.0f;
        float denom1 = 0.0f;
        float denom2 = 0.0f;
        float temp1, temp2;
        for (int i = 0; i < valuesList1.size(); i++) {
            temp1 = valuesList1.get(i) - mean1;
            temp2 = valuesList2.get(i) - mean2;
            num += temp1 * temp2;
            denom1 += temp1 * temp1;
            denom2 += temp2 * temp2;
        }
        return num / (float)Math.sqrt(denom1 * denom2);
    }
       
    private float mean(List<Float> values) {
        float total = 0.0f;
        int c = 0;
        for (Float f : values) {
            if (f == null)
                continue;
            total += f;
            c++;
        }
        return total / c;
    }
    
    private Float[] parseValueLine(String line) {
        String[] tokens = line.split("\t");
        Float[] values = new Float[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
            if (tokens[i].equals("NULL"))
                continue;
            values[i] = Float.parseFloat(tokens[i]);
        }
        return values;
    }
    
    public void anova() throws IOException {
        List<String> gdsList = getGDSList();
        List<Double> values = new ArrayList<Double>();
        int totalSample = 0;
        double total = 0.0f;
        double ssdw = 0.0f;
        List<Double> means = new ArrayList<Double>();
        List<Integer> sampleSizeList = new ArrayList<Integer>();
        FileUtility fu = new FileUtility();
        String[] tokens = null;
        String line = null;
        for (String gds : gdsList) {
            String fileName = RESULT_DIR + gds + ".ave.txt";
            fu.setInput(fileName);
            values.clear();
            while ((line = fu.readLine()) != null) {
                tokens = line.split("\t");
                for (int i = 1; i < tokens.length; i++) {
                    if (tokens[i].equals("NULL"))
                        continue;
                    Double v = Double.parseDouble(tokens[i]);
                    if (v.equals(Double.NaN))
                        System.out.println("Not a number: " + tokens[i]);
                    values.add(Double.parseDouble(tokens[i]));
                }
            }
            double total1 = 0.0f;
            for (double f : values) 
                total1 += f;
            double mean = total1 / values.size();
            means.add(mean);
            for (double f : values)
                ssdw += (mean - f) * (mean - f);
            fu.close();
            totalSample += values.size();
            sampleSizeList.add(values.size());
            total += total1;
        }
        // Calculate ssdb
        double ssdb = 0.0f;
        double totalMean = total / totalSample;
        int index = 0;
        for (double mean : means) {
            int sampleSize = sampleSizeList.get(index ++);
            ssdb += sampleSize * (mean - totalMean) * (mean - totalMean);
        }
        System.out.println("SSDW: " + ssdw);
        System.out.println("SSDB: " + ssdb);
        System.out.println("Total Sample: " + means.size());
        double f = (ssdb / (means.size() - 1)) / (ssdw / (totalSample - means.size()));
        System.out.println("DF1: " + (means.size() - 1));
        System.out.println("DF2: " + (totalSample - means.size()));
        System.out.println("F Value: " + f);
    }
    
    public void generateRandomGenePair() throws IOException {
        loadData(RESULT_DIR + "GDS824.ave.txt");
        int[] index1 = new int[100];
        int[] index2 = new int[100];
        int total = id2ValueLine.size();
        int sample = 100;
        for (int i = 0; i < sample; i++) {
            index1[i] = (int) (Math.random() * total);
            index2[i] = (int) (Math.random() * total);
        }
        List<String> ids = new ArrayList<String>(id2ValueLine.keySet());
        FileUtility fu = new FileUtility();
        fu.setOutput(RESULT_DIR + "RandomeCorBasedGDS824.txt");
        Set<String> touched = new HashSet<String>();
        int c = 0;
        for (int i : index1) {
            String id1 = ids.get(i);
            for (int j : index2) {
                String id2 = ids.get(j);
                if (id1.equals(id2)) {
                    continue;
                }
                int comp = id1.compareTo(id2);
                String key = null;
                if (comp < 0)
                    key = id1 + " " + id2;
                else
                    key = id2 + " " + id1;
                if (touched.contains(key))
                    continue;
                touched.add(key);
                Float cor = calculateCorrelation(id1, id2);
                if (cor == null)
                    continue;
                fu.printLine(cor.toString());
                c ++;
            }
        }
        fu.close();
        System.out.println("Total: " + c);
    }
    
    public void countMap() throws IOException {
        File dir = new File(RESULT_DIR);
        File[] list = dir.listFiles();
        FileUtility fu = new FileUtility();
        String line = null;
        int c = 0;
        for (File f : list) {
            String fileName = f.getName();
            if (fileName.endsWith(".map")) {
                fu.setInput(f.getAbsolutePath());
                c = 0;
                while ((line = fu.readLine()) != null)
                    c ++;
                fu.close();
                System.out.println(fileName + ": " + c);
            }
        }
    }
    
    /**
     * This method is used to extract positive correlated Reactome IDs
     * @throws IOException
     */
    public void checkReactomeIDsInGeneExp() throws Exception {
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Set<String> pairs = reactomeAnalyzer.generateUniProtPairsFromTopics();
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(pairs);
        String geneExpDataName = "results/microarray/GeneExpFromPavlidis.txt";
        // Format: "+/-" + "\t" + id1 + " " + id2 (id1 and id2 are sorted)
        FileUtility fu = new FileUtility();
        fu.setInput(geneExpDataName);
        String line = null;
        int index = 0;
        String pair = null;
        String geneExpValue = null;
        String id1, id2;
        FileUtility outputFu = new FileUtility();
        outputFu.setOutput("results/microarray/ReactomePos.txt");
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            pair = line.substring(index + 1);
            geneExpValue = line.substring(0, index);
            if (!geneExpValue.equals("+"))
                continue;
            index = pair.indexOf(" ");
            id1 = pair.substring(0, index);
            id2 = pair.substring(index + 1);
            if (ids.contains(id1) && ids.contains(id2)) {
                outputFu.printLine(pair);
            }
        }
        fu.close();
        outputFu.close();
    }
    
    public void generateReactomeIDFromPos() throws IOException {
        FileUtility fu = new FileUtility();
        Set<String> interactions = fu.loadInteractions("results/microarray/ReactomePos.txt");
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(interactions);
        fu.setOutput("results/microarray/ReactomePosIDs.txt");
        for (String id : ids)
            fu.printLine(id);
        fu.close();
    }
    
    @Test
    public void removeRendundancyInPairs() throws IOException {
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniAnalyzer.loadUniProtIDsMap();
        String geneExpDataName = "results/microarray/GeneExpWith3FromPavlidis.txt";
        Set<String> set = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(geneExpDataName);
        String line = null;
        int index = 0;
        String pair = null;
        String geneExpValue = null;
        Set<String> lines = new HashSet<String>();
        int total = 0;
        String value = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            value = line.substring(0, index);
            pair = line.substring(index + 1);
            index = pair.indexOf(" ");
            String acc1 = pair.substring(0, index);
            if (acIdMap.containsKey(acc1))
                acc1 = acIdMap.get(acc1);
            String acc2 = pair.substring(index + 1);
            if (acIdMap.containsKey(acc2))
                acc2 = acIdMap.get(acc2);
            int compare = acc1.compareTo(acc2);
            if (compare < 0) {
                lines.add(value + "\t" + acc1 + " " + acc2);
            }
            else if (compare > 0) {
                lines.add(value + "\t" + acc2 + " " + acc1);
            }
            total ++;
        }
        fu.close();
        System.out.println("Original lines: " + total);
        System.out.println("Output lines: " + lines.size());
        String outputFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "GeneExpWith3FromPavlidis_041108.txt";
        fu.setOutput(outputFileName);
        for (String line1 : lines)
            fu.printLine(line1);
        fu.close();
    }
    
    public Set<String> getGeneExpPairWiseData() throws IOException {
        //String geneExpDataName = "results/microarray/GeneExpWith3FromPavlidis.txt";
        String geneExpDataName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "GeneExpWith3FromPavlidis_041108.txt";
        Set<String> set = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(geneExpDataName);
        String line = null;
        int index = 0;
        String pair = null;
        String geneExpValue = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            pair = line.substring(index + 1);
            set.add(pair);
        }
        fu.close();
        return set;
    }
   
}
