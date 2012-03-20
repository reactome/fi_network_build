/*
 * Created on Jul 26, 2006
 *
 */
package org.reactome.data;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.funcInt.Protein;

/**
 * This class is used to process data from HPRD.
 * @author guanming
 *
 */
public class HPRDAnalyzer extends CPathAnalyzer {
    public static final String HPRD_DIR = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/";
    
    public HPRDAnalyzer() {
        //dataSourceId = 317977L;
        dataSourceId = 816565L;
    }
    
    /**
     * Nothing should be returned.
     */
    public Set<String> generateUniProtPairsFromTopics() throws Exception {
        return new HashSet<String>();
    }
    
    public void extractExperimentMethod() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        GKInstance ds = getDataSource();
        Collection summations = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Summation, 
                                                             ReactomeJavaConstants.dataSource, 
                                                             "=",
                                                             ds);
        Set<String> methods = new HashSet<String>();
        for (Iterator it = summations.iterator(); it.hasNext();) {
            GKInstance summation = (GKInstance) it.next();
            methods.add(summation.getDisplayName());
        }
        System.out.println("Methods: " + methods);
    }
    
    protected Collection prepareInteractions() throws Exception {
        // Load necessary attributes
        MySQLAdaptor dba = (MySQLAdaptor)getMySQLAdaptor();
        GKInstance hprd = getDataSource();
        Collection interactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Interaction,
                                                            ReactomeJavaConstants.dataSource,
                                                            "=",
                                                            hprd);
        // Just want to have 2H dataset
        List<GKInstance> y2hInteractions = new ArrayList<GKInstance>();
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            GKInstance interaction = (GKInstance) it.next();
            GKInstance summation = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.summation);
            if (summation == null)
                continue;
            if (summation.getDisplayName().equals("Experiment: . Detection method: 2H"))
                y2hInteractions.add(interaction);
        }
        // Load precedingEvent values
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Interaction);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.interactor);
        dba.loadInstanceAttributeValues(y2hInteractions, att);
        return y2hInteractions;
    }
    
    public Map<String, Set<Integer>> loadInteractionToPmidMap() throws Exception {
        // Load necessary attributes
        MySQLAdaptor dba = (MySQLAdaptor)getMySQLAdaptor();
        GKInstance hprd = getDataSource();
//        Collection interactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Interaction,
//                                                               ReactomeJavaConstants.dataSource,
//                                                               "=",
//                                                               hprd);
        Collection interactions = dba.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        // Try to speed up a little bit
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Interaction);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.interactor);
        dba.loadInstanceAttributeValues(interactions, att);
        Map<String, Set<Integer>> interactionToPmids = new HashMap<String, Set<Integer>>();
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            GKInstance instance = (GKInstance) it.next(); // An interaction
            // Get interaction in gene name
            List interactors = instance.getAttributeValuesList(ReactomeJavaConstants.interactor);
            List<String> geneNames = extractGeneNames(interactors);
            if (geneNames.size() < 2) {
                System.err.println(instance + " has less than two genes!");
                continue;
            }
            Set<Integer> pmids = extractPmids(instance);
            for (int i = 0; i < geneNames.size() - 1; i++) {
                String gene1 = geneNames.get(i);
                for (int j = i + 1; j < geneNames.size(); j++) {
                    String gene2 = geneNames.get(j);
                    String fi = gene1 + "\t" + gene2;
                    Set<Integer> set = interactionToPmids.get(fi);
                    if (set == null)
                        interactionToPmids.put(fi, pmids);
                    else
                        set.addAll(pmids);
                }
            }
        }
        return interactionToPmids;
    }
    
    @Test
    public void dumpInteractionToPMIDs() throws Exception {
        Map<String, Set<Integer>> interactionToPmids = loadInteractionToPmidMap();
        FileUtility fu = new FileUtility();
        //String fileName = R3Constants.RESULT_DIR + "HPRDInteraction2PMIDS.txt";
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "Interaction2PMIDs.txt";
        fu.setOutput(fileName);
        for (String fi : interactionToPmids.keySet()) {
            Set<Integer> pmids = interactionToPmids.get(fi);
            for (Integer pmid : pmids)
                fu.printLine(fi + "\t" + pmid);
        }
        fu.close();
    }
    
    private Set<Integer> extractPmids(GKInstance interaction) throws Exception {
        Set<Integer> rtn = new HashSet<Integer>();
        // PMIDs in summations
        List summations = interaction.getAttributeValuesList(ReactomeJavaConstants.summation);
        if (summations != null && summations.size() > 0) {
            for (Iterator it = summations.iterator(); it.hasNext();) {
                GKInstance summation = (GKInstance) it.next();
                GKInstance lit = (GKInstance) summation.getAttributeValue(ReactomeJavaConstants.literatureReference);
                if (lit == null)
                    continue;
                Integer pmid = (Integer) lit.getAttributeValue(ReactomeJavaConstants.pubMedIdentifier);
                if (pmid != null)
                    rtn.add(pmid);
            }
        }
        return rtn;
    }
    
    private List<String> extractGeneNames(List interactors) throws Exception {
        List<String> names = new ArrayList<String>();
        for (Iterator it = interactors.iterator(); it.hasNext();) {
            GKInstance interactor = (GKInstance) it.next();
            String geneName = null;
            if (interactor.getSchemClass().isValidAttribute(ReactomeJavaConstants.referenceEntity)) {
                GKInstance re = (GKInstance) interactor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
                if (re != null)
                    geneName = (String) re.getAttributeValue(ReactomeJavaConstants.geneName);
            }
            if (geneName == null) {
                //System.err.println(re + " has no gene name!");
                geneName = interactor.getDisplayName();
            }
            names.add(geneName);
        }
        Collections.sort(names);
        return names;
    }
    
    public void map() throws IOException {
        FileUtility fu = new FileUtility();
        Map<String, String> hprdToUniProt = fu.importMap("results/HPRD2UniProtViaGeneSymbols.txt");
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/Test.txt";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String uniprotId = hprdToUniProt.get(tokens[0]);
            System.out.printf("%s %s%n", tokens[0], uniprotId);
        }
        fu.close();
    }
    
    public void mapToUniProtViaGeneSymbols() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/";
        // map from gene symbols to uniprot ids
        String fileName = dirName + "GeneSymbolToUniProtFromHGNC.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        Map<String, String> geneToUniProt = new HashMap<String, String>();
        String line = fu.readLine(); // header
        while ((line = fu.readLine()) != null) {
            // HGNC ID\tApproved Symbol\tApproved Name\tUniProt ID (mapped data)
            String[] tokens = line.split("\t");
            if (tokens.length < 4)
                continue; // cannot be mapped
            if (tokens[3] != null && tokens[3].length() > 0)
                geneToUniProt.put(tokens[1], tokens[3]);
        }
        fu.close();
        System.out.println("Total Gene Symbols to UniProts: " + geneToUniProt.size());
        // map HPRD to gene symbols
        dirName = dirName + "/HPRD_XML_060106";
        File dir = new File(dirName);
        Map<String, String> hprdToGene = new HashMap<String, String>();
        File[] files = dir.listFiles();
        int index = 0;
        String hprdID = null;
        String geneSymbol = null;
        for (File file : files) {
            fileName = file.getName();
            index = fileName.indexOf(".");
            hprdID = fileName.substring(0, index);
            fu.setInput(file.getAbsolutePath());
            while ((line = fu.readLine()) != null) {
                if (line.startsWith("    <gene_symbol>")) {
                    index = line.indexOf("</");
                    geneSymbol = line.substring(17, index);
                    hprdToGene.put(hprdID, geneSymbol);
                    break;
                }
            }
            fu.close();
        }
        System.out.println("Total HPRD to gene symbols: " + hprdToGene.size());
        // Map from HPRD to UniProts
        Map<String, String> hprdToUniProt = new HashMap<String, String>();
        String uniProt = null;
        for (Iterator<String> it = hprdToGene.keySet().iterator(); it.hasNext();) {
            hprdID = it.next();
            geneSymbol = hprdToGene.get(hprdID);
            uniProt = geneToUniProt.get(geneSymbol);
            if (uniProt != null)
                hprdToUniProt.put(hprdID, uniProt);
        }
        System.out.println("Total HPRD to UniProt: " + hprdToUniProt.size());
        fu.exportMap(hprdToUniProt, "results/HPRD2UniProtViaGeneSymbols.txt");
    }
        
    public void extractUniProtIDs() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD_XML_060106";
        File dir = new File(dirName);
        Map<String, String> map = new HashMap<String, String>();
        File[] files = dir.listFiles();
        String fileName = null;
        int index = 0;
        String hprdID = null;
        String uniID = null;
        FileUtility fu = new FileUtility();
        String line = null;
        for (File file : files) {
            fileName = file.getName();
            index = fileName.indexOf(".");
            hprdID = fileName.substring(0, index);
            fu.setInput(file.getAbsolutePath());
            while ((line = fu.readLine()) != null) {
                if (line.equals("    <EXTERNAL_LINKS>")) {
                    while (true) {
                        line = fu.readLine();
                        if (line.equals("    </EXTERNAL_LINKS>"))
                            break;
                        line = line.trim();
                        if (line.startsWith("<SwissProt>")) {
                            index = line.indexOf("</");
                            uniID = line.substring(11, index);
                            if (!uniID.equals("None"))
                                map.put(hprdID, uniID);
                            break;
                        }
                    }
                }
            }
            fu.close();
        }
        fu.exportMap(map, "results/HPRD2UniProt.txt");
    }
    
    public List<Protein> loadHPRDWithNames() throws IOException {
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD2Names060106.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        List<Protein> proteins = new ArrayList<Protein>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            Protein protein = new Protein();
            protein.setPrimaryDbName("HPRD");
            protein.setPrimaryAccession(tokens[0]);
            if (!tokens[1].equals("null"))
                protein.setName(tokens[1]);
            if (!tokens[2].equals("null"))
                protein.setShortName(tokens[2]);
            proteins.add(protein);
        }
        fu.close();
        return proteins;
    }
    
    @Test
    public void generateHPRDToNames() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD_XML_060106";
        File dir = new File(dirName);
        File[] files = dir.listFiles();
        String filename = null;
        int index = 0;
        int index1 = 0;
        String hprdID = null;
        String geneName = null;
        FileUtility fu = new FileUtility();
        String line = null;
        String title = null;
        // For output
        String output = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD2Names060106.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(output);
        for (File file : files) {
            if (title != null) {
                outFu.printLine(hprdID + "\t" + title + "\t" + geneName);
            }
            filename = file.getName();
            if (!filename.endsWith(".xml"))
                continue;
            //System.out.println("FileName: " + filename);
            index = filename.indexOf(".");
            hprdID = filename.substring(0, index);
            fu.setInput(file.getAbsolutePath());
            title = geneName = null; // reset
            while ((line = fu.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("<title>")) {
                    title = getAttributeValue(line);
                    if (title.startsWith("Hypothetical protein ")) {
                        title = title.substring("Hypothetical protein ".length());
                    }
                }
                else if (line.startsWith("<gene_symbol>")) {
                    geneName = getAttributeValue(line);
                }
                else if (line.startsWith("<gene_map_locus>"))
                    break;
            }
            fu.close();
        }
        outFu.close();
    }
    
    @Test
    public void generateGeneNameToHPRDMap() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD_XML_060106";
        File dir = new File(dirName);
        File[] files = dir.listFiles();
        String filename = null;
        int index = 0;
        int index1 = 0;
        String hprdID = null;
        String geneName = null;
        FileUtility fu = new FileUtility();
        String line = null;
        String title = null;
        // For output
        String output = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/GeneName2HPRD060106.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(output);
        for (File file : files) {
            filename = file.getName();
            if (!filename.endsWith(".xml"))
                continue;
            //System.out.println("FileName: " + filename);
            index = filename.indexOf(".");
            hprdID = filename.substring(0, index);
            fu.setInput(file.getAbsolutePath());
            while ((line = fu.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("<title>")) {
                    title = getAttributeValue(line);
                    if (title.startsWith("Hypothetical protein ")) {
                        title = title.substring("Hypothetical protein ".length());
                    }
                    outFu.printLine(title + "\t" + hprdID);
                }
                else if (line.startsWith("<alt_title>") ||
                         line.startsWith("<gene_symbol>")) {
                    geneName = getAttributeValue(line);
                    if (geneName != null)
                        outFu.printLine(geneName + "\t" + hprdID);
                }
                else if (line.startsWith("<gene_map_locus>"))
                    break;
            }
            fu.close();
        }
        outFu.close();
    }
    
    private String getAttributeValue(String line) {
        int index1 = line.indexOf(">");
        int index2 = line.indexOf("<", index1 + 1);
        if (index2 < 0)
            return null;
        return line.substring(index1 + 1, index2);
    }
    
    public Map<String, String> loadGeneToHPRDMap() throws IOException {
        Map<String, String> hprdToGeneName = new FileUtility().importMap("/Users/wgm/Documents/caBIG_R3/datasets/HPRD/GeneName2HPRD060106.txt");
        return hprdToGeneName;
    }
    
    public Map<String, String> loadHprdIdToUniProtMap() throws IOException {
        Map<String, String> map = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        String fileName = FIConfiguration.getConfiguration().get("HPRD_DIR") + "HPRD_ID_MAPPINGS.txt";
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String hprd = tokens[0];
            String uniprot = tokens[6];
            if (uniprot.equals("-"))
                continue; // Should be escaped
            map.put(hprd, uniprot);
        }
        return map;
    }
    
    @Test
    public void tetsLoadHprdToUniProtMap() throws Exception {
        Map<String, String> hprdToUniprot = loadHprdIdToUniProtMap();
        System.out.println("Total mapped HPRDs: " + hprdToUniprot.size());
    }

    /**
     * Use the id mapping from HPRD to map Gene ids to UniProt accessions.
     * @return
     * @throws Exception
     */
    public Map<String, String> loadGeneIdToUniProtAccMap() throws IOException {
        Map<String, String> geneIdToUniProt = new HashMap<String, String>();
        String fileName = FIConfiguration.getConfiguration().get("HPRD_DIR") + "HPRD_ID_MAPPINGS.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String geneId = tokens[4];
            String uniProtId = tokens[6];
            if (uniProtId.equals("-"))
                continue;
            geneIdToUniProt.put(geneId, uniProtId);
        }
        fu.close();
        return geneIdToUniProt;
    }
    
    @Test
    public void testLoadGeneIdToUniProtAccMap() throws IOException {
        Map<String, String> geneIdToUniProt = loadGeneIdToUniProtAccMap();
        System.out.println("Size: " + geneIdToUniProt.size());
    }
}
