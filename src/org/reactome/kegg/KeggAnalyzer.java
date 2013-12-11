/*
 * Created on Jan 30, 2007
 *
 */
package org.reactome.kegg;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.CharBuffer;
import java.nio.channels.FileChannel;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.junit.Test;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to process KEGG data
 * @author guanming
 *
 */
public class KeggAnalyzer{
    //private final String DIR_NAME = "/Users/wgm/Documents/caBIG_R3/datasets/KEGG/";
    private final String DIR_NAME = FIConfiguration.getConfiguration().get("KEGG_DIR");
    
    public KeggAnalyzer() {
    }
    
    @Test
    public void checkNonMetabolicPathways() throws Exception {
        String dirName = DIR_NAME + "KGML/hsa/";
        File dir = new File(dirName);
        for (File file : dir.listFiles()) {
            String fileName = file.getName();
            if (fileName.endsWith(".xml")) {
                SAXBuilder builder = new SAXBuilder();
                Document doc = builder.build(file);
                Element root = doc.getRootElement();
                String org = root.getAttributeValue("org");
                String number = root.getAttributeValue("number");
                String title = root.getAttributeValue("title");
                System.out.println(org + number + "\t" + title);
            }
        }
    }
    
    public Map<String, String> loadKeggToUniProtMaps() throws Exception {
        // Used to map to a non-redundant UniProt id
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtMap = uniAnalyzer.loadUniProtIDsMap();
        //hsa:1   up:P04217
        //hsa:10  up:P11245
        Map<String, String> map = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        //String fileName = DIR_NAME + "hsa_genes_uniprot.list";
        String fileName = FIConfiguration.getConfiguration().get("KEGG_ID_TO_UNIPROT_MAP_FILE"); // For KEGG release downloaded on March 12, 2009
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String kegg = line.substring(0, index);
            String up = line.substring(index + 1);
            String uniProtId = up.substring(3);
            String mapped = uniProtMap.get(uniProtId);
            if (mapped == null)  {// Just in case
                System.err.println(uniProtId + " cannot be found in the download UniProt files!");
                mapped = uniProtId;
            }
            map.put(kegg, mapped);
        }
        fu.close();
        return map;
    }
    
    /**
     * A KEGG id may be mapped to multiple UniProt ids. The above method loadKeggToUniProtMaps() can
     * map one KEGG id to one UniProt accession number. Use this method to get all mappings.
     * @return
     * @throws Exception
     */
    public Map<String, Set<String>> loadKeggToAllUniProtIdsMap() throws Exception {
        // Used to map to a non-redundant UniProt id
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniProtMap = uniAnalyzer.loadUniProtIDsMap();
        //hsa:1   up:P04217
        //hsa:10  up:P11245
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        FileUtility fu = new FileUtility();
        String fileName = FIConfiguration.getConfiguration().get("KEGG_ID_TO_UNIPROT_MAP_FILE"); // For KEGG release downloaded on March 12, 2009
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String kegg = line.substring(0, index);
            String up = line.substring(index + 1);
            String uniProtId = up.substring(3);
            String mapped = uniProtMap.get(uniProtId);
            if (mapped == null)  {// Just in case
                System.err.println(uniProtId + " cannot be found in the download UniProt files!");
                mapped = uniProtId;
            }
            InteractionUtilities.addElementToSet(map, kegg, mapped);
        }
        fu.close();
        return map;
    }
    
    @Test
    public void testLoadKeggToAllUniProtIdsMap() throws Exception {
        Map<String, Set<String>> keggIdToUniProts = loadKeggToAllUniProtIdsMap();
        int multipleCounts = 0;
        for (String keggId : keggIdToUniProts.keySet()) {
            Set<String> uniProtIds = keggIdToUniProts.get(keggId);
            if (uniProtIds.size() > 1)
                multipleCounts ++;
        }
        System.out.println("Total KEGG IDs: " + keggIdToUniProts.size());
        System.out.println("Multiple counts: " + multipleCounts);
    }
    
    @Test
    public void testGetNonMetabolismPathways() throws Exception {
        List<String> pathways = getNonMetabolismPathways();
        System.out.println("Total non-metabolic pathways: " + pathways.size());
        for (String pathway : pathways)
            System.out.println(pathway);
    }
    
    /**
     * Parse the index.html file to get the list of all non-metabolism pathways.
     * Only non-metabolism pathways will be used for data analysis.
     * @return
     * @throws Exception
     */
    public List<String> getNonMetabolismPathways() throws Exception {
        List<String> list = new ArrayList<String>();
        String fileName = DIR_NAME + "hsa/index.html";
        FileInputStream fis = new FileInputStream(fileName);
        FileChannel channel = fis.getChannel();
        ByteBuffer buffer = ByteBuffer.allocate((int)channel.size()).order(ByteOrder.LITTLE_ENDIAN);
        channel.read(buffer);
        buffer.rewind();
        Pattern pattern = Pattern.compile("<BIG><B>(.+)</B></BIG>");
        //Matcher matcher = pattern.matcher("This is a test<BIG><B>Pathway Names</B></BIG>");
        Charset charset = Charset.forName("UTF-8");
        CharBuffer charBuffer = charset.decode(buffer);
        Matcher matcher = pattern.matcher(charBuffer);
        int start = 0;
        boolean needGet = false;
        while (matcher.find(start)) {
            String name = matcher.group(1);
            //System.out.println(name);
            if (name.contains("Sorting and Degradation")) { // Don't include pathways from Transcription and Translation, most of
                                                            // which are just big structures.
                getPathways(matcher.end(), list, charBuffer);
                break;
            }
            start = matcher.end();
        }
        return list;
    }
    
    private void getPathways(int start,
                             List<String> pathways,
                             CharBuffer charBuffer) throws Exception {
        Pattern pattern = Pattern.compile("<B>hsa(\\d+)</B>");
        Matcher matcher = pattern.matcher(charBuffer);
        while (matcher.find(start)) {
            start = matcher.end();
            String name = matcher.group(1);
            if (name.equals("01430"))
                continue; // Escape this pathway: Cell junction, which is basically a cell junction structure drawing.
            pathways.add("hsa" + name);
        }
    }
    
    /**
     * This method is used to generate a map from pathway to genes based on downloaded KEGG file
     * has.list and map_title.tab.
     * @return
     * @throws IOException
     */
    public Map<String, Set<String>> loadPathwayToGeneNamesMapFromHsaList() throws IOException {
        String dirName = FIConfiguration.getConfiguration().get("KEGG_DIR");
        String listFile = dirName + File.separator + "hsa" + File.separator + "hsa.list";
        String titleFile = dirName + File.separator + "map_title.tab";
        Map<String, Set<String>> pathwayToGenes = new HashMap<String, Set<String>>();
        FileUtility fu = new FileUtility();
        Map<String, String> pathwayIdToName = new HashMap<String, String>();
        // Load names
        fu.setInput(titleFile);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            pathwayIdToName.put("path:hsa" + tokens[0], tokens[1] + "(K)");
        }
        fu.close();
        // Load genes
        fu.setInput(listFile);
        int index;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            if (tokens.length < 3 || tokens[2].trim().length() == 0) // For some bug from KEGG that generates an empty string token
                continue; // Small molecules and other types
            String pathwayName = pathwayIdToName.get(tokens[0]);
            index = tokens[2].indexOf(" ");
            String gene = tokens[2].substring(4, index);
            InteractionUtilities.addElementToSet(pathwayToGenes, pathwayName, gene);
        }
        fu.close();
        return pathwayToGenes;
    }
    
    public void augmentGeneNameToPathwayMap(String srcFileName,
                                            String targetFileName) throws IOException {
        Map<String, Set<String>> pathwayToNames = loadPathwayToGeneNamesMapFromHsaList();
        Map<String, Set<String>> nameToPathways = InteractionUtilities.switchKeyValues(pathwayToNames);
        FileUtility fu = new FileUtility();
        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ProteinNameToTopics101110.txt";
        fu.setInput(srcFileName);
        //String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ProteinNameToTopics_all_KEGG_101110.txt";
        fu.setOutput(targetFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.endsWith("(K)"))
                continue; // Escape all KEGG pathways
            fu.printLine(line);
        }
        for (String name : nameToPathways.keySet()) {
            Set<String> pathways = nameToPathways.get(name);
            for (String pathway : pathways)
                fu.printLine(name + "\t" + pathway);
        }
        fu.close();
    }
    
    @Test
    public void testGeneratePathwayToGeneNamesMapFromHasList() throws IOException {
        Map<String, Set<String>> pathwayToGenes = loadPathwayToGeneNamesMapFromHsaList();
        Set<String> totalGenes = InteractionUtilities.grepAllGenes(pathwayToGenes);
        System.out.println("Total genes: " + totalGenes.size());
        List<String> pathways = new ArrayList<String>(pathwayToGenes.keySet());
        Collections.sort(pathways);
        for (String pathway : pathways) {
            Set<String> genes = pathwayToGenes.get(pathway);
            System.out.println(pathway + ": " + genes.size() + " " + genes);
        }
    }
    
}
