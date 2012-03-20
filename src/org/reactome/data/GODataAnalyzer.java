/*
 * Created on Apr 28, 2006
 *
 */
package org.reactome.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import junit.framework.TestCase;

import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to handle GO related data.
 * @author guanming
 *
 */
public class GODataAnalyzer extends TestCase {
    private final String RESULT_DIR = "results/go/";
    
    public GODataAnalyzer() {
    }   
    
    public Map<String, Set<String>> convertIdsToNames(Map<String, Set<String>> geneToIds) throws IOException {
        Map<String, String> idToName = loadGOIdToTermMap();
        Map<String, Set<String>> geneToNames = new HashMap<String, Set<String>>();
        for (String gene : geneToIds.keySet()) {
            Set<String> ids = geneToIds.get(gene);
            Set<String> names = new HashSet<String>(ids.size());
            for (String id : ids) {
                String name = idToName.get(id);
                names.add(name);
            }
            geneToNames.put(gene, names);
        }
        return geneToNames;
    }
    
    public Map<String, String> loadGOIdToTermMap() throws IOException {
        FileUtility fu = new FileUtility();
        //String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/GO/GO.terms_and_ids.txt";
        String fileName = FIConfiguration.getConfiguration().get("GO_DIR") + "GO.terms_and_ids.txt";
        fu.setInput(fileName);
        String line = null;
        Map<String, String> id2Term = new HashMap<String, String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("!"))
                continue;
            String[] tokens = line.split("\t");
            id2Term.put(tokens[0], tokens[1]);
        }
        return id2Term;
    }
    
    public Map<String, List<String>> loadProtein2LocationTerms() throws IOException {
        String fileName = "results/go/goa_c.txt";
        Map<String, List<String>> protein2Locations = new HashMap<String, List<String>>();
        FileUtility fileUtility = new FileUtility();
        fileUtility.setInput(fileName);
        String line = null;
        while ((line = fileUtility.readLine()) != null) {
            String[] tokens = line.split("\\t");
            List<String> terms = protein2Locations.get(tokens[0]);
            if (terms == null) {
                terms = new ArrayList<String>();
                protein2Locations.put(tokens[0], terms);
            }
            terms.add(tokens[3]);
        }
        fileUtility.close();
        return protein2Locations;
    }
    
    public Map<String, List<String>> loadProtein2LocationIds() throws IOException {
        String fileName = "results/go/goa_c.txt";
        Map<String, List<String>> protein2Locations = new HashMap<String, List<String>>();
        FileUtility fileUtility = new FileUtility();
        fileUtility.setInput(fileName);
        String line = null;
        while ((line = fileUtility.readLine()) != null) {
            String[] tokens = line.split("\\t");
            List<String> ids = protein2Locations.get(tokens[0]);
            if (ids == null) {
                ids = new ArrayList<String>();
                protein2Locations.put(tokens[0], ids);
            }
            ids.add(tokens[2]);
        }
        fileUtility.close();
        return protein2Locations;
    }
    
    protected boolean isColocalized(List<String> localizations1, List<String> localizations2) {
        for (String loc : localizations1) {
            if (localizations2.contains(loc))
                return true;
        }
        return false;
    }
    
    public void analyzeCellularLocations() throws IOException {
        Map<String, List<String>> protein2Locations = loadProtein2LocationTerms();
        // check how many proteins are outside of cell
        //int cellSurfaceNumber = 0;
        List<String> proteinsInNucleus = new ArrayList<String>();
        for (Iterator<String> it = protein2Locations.keySet().iterator(); it.hasNext();) {
            String protein = it.next();
            List<String> terms = protein2Locations.get(protein);
//            if (terms.contains("cell surface")) {
//                cellSurfaceNumber ++;
//                System.out.println(protein + ": " + terms);
//            }
            if (terms.size() == 1) {
                String term = terms.get(0);
                if (term.equals("nucleus")) {
                    proteinsInNucleus.add(protein);
                }
            }
        }
        //System.out.println("cell surface: " + cellSurfaceNumber);
        System.out.println("nucleusNumber: " + proteinsInNucleus.size());
        // Want to pick 100 proteins in nucleus randomly
//        List<String> proteinsInNucleusForTraining = new ArrayList<String>();
//        Set<Integer> touchedNumber = new HashSet<Integer>();
//        int generatedIndex = 0;
//        int size = proteinsInNucleus.size();
//        // Pick 2400 from the list of nucleous proteins. The left should be
//        // used for testing.
//        for (int i = 0; i < 2400; i++) {
//            while (true) {
//                generatedIndex = (int) (Math.random() * size);
//                if (!touchedNumber.contains(generatedIndex)) {
//                    touchedNumber.add(generatedIndex);
//                    break;
//                }
//            }
//            proteinsInNucleusForTraining.add(proteinsInNucleus.get(generatedIndex));
//        }
//        System.out.println("Proteins in Nucleus for Training: " + proteinsInNucleusForTraining.size());
//        List<String> surfaceProteins = loadSurfaceProteins();
//        String fileName = "results/interaction/NoInteractions090506.txt";
//        generateProteinPairsNotInteracting(proteinsInNucleus, surfaceProteins, fileName);
//        generateProteinPairsNotInteracting(proteinsInNucleusForTraining, surfaceProteins, "results/NoInteractionsForTrain.txt");
//        proteinsInNucleus.removeAll(proteinsInNucleusForTraining);
//        System.out.println("Proteins in Nucleus for Testing: " + proteinsInNucleus.size());
//        generateProteinPairsNotInteracting(proteinsInNucleus, surfaceProteins, "results/NoInteractionsForTest.txt");
    }
    
    private void generateProteinPairsNotInteracting(List<String> proteinsInNucleus,
                                                    List<String> proteinsInSurface,
                                                    String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        int comp;
        for (String pInN : proteinsInNucleus) {
            for (String pInS : proteinsInSurface) {
                comp = pInN.compareTo(pInS);
                if (comp < 0)
                    fu.printLine(pInN + " " + pInS);
                else if (comp > 0)
                    fu.printLine(pInS + " " + pInN);
            }
        }
        fu.close();
    }
    
    private List<String> loadSurfaceProteins() throws IOException {
        List<String> proteins = new ArrayList<String>();
        String fileName = "results/go/ProteinsInCellSurface.txt";
        FileReader fileReader = new FileReader(fileName);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = null;
        int index = 0;
        String protein = null;
        while ((line = reader.readLine()) != null) {
            index = line.indexOf(":");
            protein = line.substring(0, index);
            proteins.add(protein);
        }
        return proteins;
    }
    
    public void generateShareDepth() throws IOException {
        Map<String, Integer> goDepthMap = loadMap("results/GODepth.txt");
        // Need to create twp maps from Proteins to GO terms
        Map<String, Set<String>> protein2GOF = new HashMap<String, Set<String>>();
        Map<String, Set<String>> protein2GOP = new HashMap<String, Set<String>>();
        generateProtein2GOMap(protein2GOF, protein2GOP);
        System.out.println("Protein2GOMap is done: " + protein2GOF.size() + 
                "\t" + protein2GOP.size());
        generateShareDepth(protein2GOF, goDepthMap, "results/ShareGOFDepth.txt");
        generateShareDepth(protein2GOP, goDepthMap, "results/ShareGOPDepth.txt");
    }
    
    private void generateShareDepth(Map<String, Set<String>> protein2GO,
                                    Map<String, Integer> goDepthMap,
                                    String outputFileName) throws IOException {
        long time1 = System.currentTimeMillis();
        FileWriter fileWriter = new FileWriter(outputFileName);
        PrintWriter writer = new PrintWriter(fileWriter);
        List<String> proteinIds = new ArrayList<String>(protein2GO.keySet());
        int size = proteinIds.size();
        for (int i = 0; i < size; i++) {
            String id1 = proteinIds.get(i);
            Set<String> terms1 = protein2GO.get(id1);
            for (int j = i + 1; j < size; j++) {
                String id2 = proteinIds.get(j);
                // Find share term
                Set<String> terms2 = protein2GO.get(id2);
                int sharedOccurence = findSharedDepth(terms1, terms2, goDepthMap);
                if (sharedOccurence != Integer.MIN_VALUE) {
                    writer.println(id1 + "\t" + id2 + "\t" + sharedOccurence);
                }
            }
        }
        writer.close();
        fileWriter.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateShareUsage: " + (time2 - time1));
    }
    
    public int findSharedDepth(Set<String> terms1, 
                                  Set<String> terms2, 
                                  Map<String, Integer> goDepthMap) {
        int min = Integer.MIN_VALUE;
        if (terms1 == null || terms2 == null)
            return min;
        int max = min;
        int value;
        for (String term : terms1) {
            if (terms2.contains(term)) {
                value = goDepthMap.get(term);
                if (value > max)
                    max = value;
            }
        }
        return max;
    }
    
    private Set<String> getProteinsInOneCompartment(String goTermId,
                                                    Map<String, List<String>> protein2Locations) throws IOException {
        Map<String, GOGraphNode> idToNodeMap = generateGOTree();
        // extracelluar region: GO:0005576
        //GOGraphNode node = idToNodeMap.get("GO:0005576");
        // nucleus: GO:0005634
        //GOGraphNode node = idToNodeMap.get("GO:0005634");
        GOGraphNode node = idToNodeMap.get(goTermId);
        Set<GOGraphNode> targetNodes = new HashSet<GOGraphNode>();
        grepChildren(targetNodes, node);
        Set<String> surfaceIds = new HashSet<String>();
        for (GOGraphNode graphNode : targetNodes) {
            surfaceIds.add(graphNode.getId());
        }
        Set<String> proteinsInTargetCompartment = new HashSet<String>();
        for (Iterator<String> it = protein2Locations.keySet().iterator(); it.hasNext();) {
            String protein = it.next();
            List<String> locations = protein2Locations.get(protein);
            if (locations != null && locations.size() == 1) {
                String location = locations.get(0);
                if (surfaceIds.contains(location)) {
                    proteinsInTargetCompartment.add(protein);
                }
            }
        }
        //System.out.println("Total proteins: " + proteinsInTargetCompartment.size());
        return proteinsInTargetCompartment;
    }
    
    public Set<String> generateProteinPairsInNucleusAndExRegion(int totalSize) throws Exception {
        Map<String, List<String>> protein2Locations = loadProtein2LocationIds();
        // Get proteins in nucleus only: nucleus: GO:0005634
        Set<String> proteinsInNucleus = getProteinsInOneCompartment("GO:0005634", 
                                                                    protein2Locations);
        // Get proteins in extracelluar region only: GO:0005576
        Set<String> proteinsInExtra = getProteinsInOneCompartment("GO:0005576",
                                                                  protein2Locations);
        Set<String> rtnSet = new HashSet<String>();
        List<String> list1 = new ArrayList<String>(proteinsInNucleus);
        int size1 = list1.size();
        List<String> list2 = new ArrayList<String>(proteinsInExtra);
        int size2 = list2.size();
        while (rtnSet.size() < totalSize) {
            int index1 = (int) (size1 * Math.random());
            int index2 = (int) (size2 * Math.random());
            String protein1 = list1.get(index1);
            String protein2 = list2.get(index2);
            int compare = protein1.compareTo(protein2);
            if (compare < 0)
                rtnSet.add(protein1 + " " + protein2);
            else if (compare > 0) // In case they are the same
                rtnSet.add(protein2 + " " + protein1); 
        }
        return rtnSet;
    }
    
    private void grepChildren(Set<GOGraphNode> children,
                              GOGraphNode node) {
        children.add(node);
        if (node.getChildren() == null || node.getChildren().size() == 0)
            return;
        for (GOGraphNode child : node.getChildren()) {
            grepChildren(children, child);
        }
    }
                              
    
    public Map<String, GOGraphNode> generateGOTree() throws IOException {
//        String goFileName = FileNameManager.getManager().getFileName("gene_ontology.obo");
        String goFileName = "D:\\documents\\Stein_lab\\Reactome\\Data\\gene_ontology.obo";
        Map<String, GOGraphNode> idToNodeMap = new HashMap<String, GOGraphNode>();
        FileReader fileReader = new FileReader(goFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = null;
        String id = null;
        String parentId = null;
        GOGraphNode node = null;
        GOGraphNode parentNode = null;
        int index;
        String ns = null;
        String isObsolete = null;
        List<GOGraphNode> mfOboNodes = new ArrayList<GOGraphNode>();
        List<GOGraphNode> bpOboNodes = new ArrayList<GOGraphNode>();
        List<GOGraphNode> ccOboNodes = new ArrayList<GOGraphNode>();
        while ((line = reader.readLine()) != null) {
            isObsolete = null;
            if (line.startsWith("id:")) {
                id = line.substring(4);
                node = getGraphNode(id, idToNodeMap);
            }
            else if (line.startsWith("is_a:")) {
                index = line.indexOf("!");
                parentId = line.substring(6, index - 1);
                parentNode = getGraphNode(parentId, idToNodeMap);
                node.addParent(parentNode);
            }
            else if (line.startsWith("relationship: part_of")) {
                index = line.indexOf("!");
                parentId = line.substring(22, index - 1);
                parentNode = getGraphNode(parentId, idToNodeMap);
                node.addParent(parentNode);
            }
            else if (line.startsWith("namespace:")) {
                ns = line.substring(11);
            }
            else if (line.startsWith("is_obsolete:")) {
                isObsolete = line.substring(13);
            }
            if (isObsolete != null && isObsolete.equals("true")) {
                if (ns.equals("molecular_function"))
                    mfOboNodes.add(node);
                else if (ns.equals("biological_process"))
                    bpOboNodes.add(node);
                else if (ns.equals("cellular_component"))
                    ccOboNodes.add(node);
            }
        }
        reader.close();
        fileReader.close();
        // Attach these orphan nodes to the top nodes
        String mfTerm = "GO:0003674";
        GOGraphNode mfRootNode = idToNodeMap.get(mfTerm);
        for (GOGraphNode node1 : mfOboNodes) {
            if (node1.getParents() == null || 
                node1.getParents().size() == 0) {
                // an obosolote node might be a child of another obosolete node
                node1.addParent(mfRootNode);
            }
        }
        String bpTerm = "GO:0008150";
        GOGraphNode bpRootNode = idToNodeMap.get(bpTerm);
        for (GOGraphNode node1 : bpOboNodes) {
            if (node1.getParents() == null ||
               node1.getParents().size() == 0)
               node1.addParent(bpRootNode);
        }
        String ccTerm = "GO:0005575";
        GOGraphNode ccRootNode = idToNodeMap.get(ccTerm);
        for (GOGraphNode node1 : ccOboNodes) {
            if (node1.getParents() == null ||
                node1.getParents().size() == 0)
                node1.addParent(ccRootNode);
        }
        // Need to rehash children 
        Set<GOGraphNode> parents;
        for (Iterator<String> it = idToNodeMap.keySet().iterator(); it.hasNext();) {
            id = it.next();
            node = idToNodeMap.get(id);
            parents = node.getParents();
            if (parents == null || parents.size() == 0)
                continue;
            for (GOGraphNode parent : parents) {
                parent.addChild(node);
            }
        }
        System.out.println("Total GO: " + idToNodeMap.size());
        return idToNodeMap;
    }
    
    public void calculateDepth() throws IOException {
        Map<String, GOGraphNode> idToNodeMap = generateGOTree();
        Map<String, Integer> goDepthMap = new HashMap<String, Integer>();
        for (Iterator<String> it = idToNodeMap.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            GOGraphNode node = idToNodeMap.get(id);
            List<Integer> paths = new ArrayList<Integer>();
            int depth = calculateDepth(node);
            goDepthMap.put(id, depth);
        }
        outputMap(goDepthMap, "results/GODepth.txt");
    }
    
    private int calculateDepth(GOGraphNode node) {
        Set<GOGraphNode> current = new HashSet<GOGraphNode>();
        current.add(node);
        Set<GOGraphNode> next = new HashSet<GOGraphNode>();
        int depth = -1;
        boolean isAtRoot = false;
        // Width first search
        while (!isAtRoot && current.size() > 0) {           
            for (GOGraphNode goNode : current) {
                Set<GOGraphNode> parents = goNode.getParents();
                if (parents == null || parents.size() == 0) {
                    isAtRoot = true;
                    break;
                }
                next.addAll(parents);
            }
            depth ++;
            current.clear();
            current.addAll(next);
            next.clear();
        }
        return depth;
    }
    
    private GOGraphNode getGraphNode(String id, Map<String, GOGraphNode> idToNodeMap) {
        GOGraphNode node = idToNodeMap.get(id);
        if (node == null) {
            node = new GOGraphNode();
            node.setId(id);
            idToNodeMap.put(id, node);
        }
        return node;
    }
    
    public void generateShareUsage() throws IOException {
        Map<String, Integer> fValueMap = new HashMap<String, Integer>();
        Map<String, Integer> pValueMap = new HashMap<String, Integer>();
        createUsageFrequences(fValueMap, pValueMap);
        // Need to create twp maps from Proteins to GO terms
        Map<String, Set<String>> protein2GOF = new HashMap<String, Set<String>>();
        Map<String, Set<String>> protein2GOP = new HashMap<String, Set<String>>();
        generateProtein2GOMap(protein2GOF, protein2GOP);
        System.out.println("Protein2GOMap is done: " + protein2GOF.size() + 
                "\t" + protein2GOP.size());
        // generate pairwise mapping
        generateShareUsage(protein2GOF, fValueMap, "results/ShareGOFOccurence.txt");
        generateShareUsage(protein2GOP, pValueMap, "results/ShareGOPOccurence.txt");
    }
    
    private void generateShareUsage(Map<String, Set<String>> protein2GO, 
                                    Map<String, Integer> valueMap,
                                    String outputFileName) throws IOException {
        long time1 = System.currentTimeMillis();
        FileWriter fileWriter = new FileWriter(outputFileName);
        PrintWriter writer = new PrintWriter(fileWriter);
        List<String> proteinIds = new ArrayList<String>(protein2GO.keySet());
        int size = proteinIds.size();
        for (int i = 0; i < size; i++) {
            String id1 = proteinIds.get(i);
            Set<String> terms1 = protein2GO.get(id1);
            for (int j = i + 1; j < size; j++) {
                String id2 = proteinIds.get(j);
                // Find share term
                Set<String> terms2 = protein2GO.get(id2);
                int sharedOccurence = findSharedTermOccurence(terms1, terms2, valueMap);
                if (sharedOccurence != Integer.MAX_VALUE) {
                    writer.println(id1 + "\t" + id2 + "\t" + sharedOccurence);
                }
            }
        }
        writer.close();
        fileWriter.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateShareUsage: " + (time2 - time1));
    }
    
    public int findSharedTermOccurence(Set<String> terms1, 
                                        Set<String> terms2, 
                                        Map<String, Integer> valueMap) {
        int max = Integer.MAX_VALUE; 
        if (terms1 == null || terms2 == null)
            return max;
        int min = max;
        int value;
        for (String term : terms1) {
            if (terms2.contains(term)) {
                value = valueMap.get(term);
                if (value < min)
                    min = value;
            }
        }
        return min;
    }
    
    public void countProteins() throws IOException {
        Map<String, Set<String>> protein2GOF = new HashMap<String, Set<String>>();
        Map<String, Set<String>> protein2GOP = new HashMap<String, Set<String>>();
        generateProtein2GOMap(protein2GOF, protein2GOP);
        System.out.println("Total Proteins in MF: " + protein2GOF.size());
        System.out.println("Total Proteins in BP: " + protein2GOP.size());
    }
    
    public Set<String> getAnnotatedGOBPProteins() throws IOException {
        Map<String, Set<String>> protein2GOF = new HashMap<String, Set<String>>();
        Map<String, Set<String>> protein2GOP = new HashMap<String, Set<String>>();
        generateProtein2GOMap(protein2GOF, protein2GOP);
        System.out.println("Total Proteins in MF: " + protein2GOF.size());
        System.out.println("Total Proteins in BP: " + protein2GOP.size());
        return protein2GOP.keySet();
    }
    
    public void generateProtein2GOCCMap(Map<String, Set<String>> protein2GOC) throws IOException {
        // These two terms should be escaped since there are no meanings at all
        String sourceFileName = "D:\\documents\\Stein_lab\\Reactome\\Data\\gene_association.goa_human";
//        String sourceFileName = FileNameManager.getManager().getFileName("gene_association.goa_human");
        FileUtility fu = new FileUtility();
        fu.setInput(sourceFileName);
        String line = null;
        String[] tokens = null;
        String proteinId = null;
        Set<String> goIds = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\\t");
            if (!tokens[0].equals("UniProt"))
                continue; // Work with UniProt only
            proteinId = tokens[1];
            if (tokens[8].equals("C")) {
                goIds = protein2GOC.get(proteinId);
                if (goIds == null) {
                    goIds = new HashSet<String>();
                    protein2GOC.put(proteinId, goIds);
                }
                goIds.add(tokens[4]);
            }
        }
        fu.close();
    }
    
    public void generateProtein2GOMap(Map<String, Set<String>> protein2GOF, 
                                      Map<String, Set<String>> protein2GOP) throws IOException {
        // These two terms should be escaped since there are no meanings at all
        String BP_ESCAPE_TERM = "GO:0000004"; // biological process unknown
        String MF_ESCAPE_TERM = "GO:0005554"; // molecular function unknown
        //String sourceFileName = "D:\\documents\\Stein_lab\\Reactome\\Data\\gene_association.goa_human";
        // String sourceFileName = FileNameManager.getManager().getFileName("gene_association.goa_human");
        String sourceFileName = "datasets" + File.separator + "GO" + File.separator + "gene_association.goa_human";
        FileUtility fu = new FileUtility();
        fu.setInput(sourceFileName);
        String line = null;
        String[] tokens = null;
        String proteinId = null;
        Set<String> goIds = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\\t");
                     //   if (!tokens[0].equals("UniProt"))
            //This modification is done according to the goa_human file
            // obtained at dec 28,2007
            //if (!tokens[0].equals("UniProtKB"))
            if (!tokens[0].equals("UniProt"))
                continue; // Work with UniProt only
            proteinId = tokens[1];
            if (tokens[8].equals("F")) {
                if (tokens[4].equals(MF_ESCAPE_TERM))
                    continue;
                goIds = protein2GOF.get(proteinId);
                if (goIds == null) {
                    goIds = new HashSet<String>();
                    protein2GOF.put(proteinId, goIds);
                }
                goIds.add(tokens[4]);
            }
            else if (tokens[8].equals("P")) {
                if (tokens[4].equals(BP_ESCAPE_TERM))
                    continue;
                goIds = protein2GOP.get(proteinId);
                if (goIds == null) {
                    goIds = new HashSet<String>();
                    protein2GOP.put(proteinId, goIds);
                }
                goIds.add(tokens[4]);
            }
        }
        fu.close();
    }
    
    private void createUsageFrequences(Map<String, Integer> fValueMap, Map<String, Integer> pValueMap) throws IOException {
        String sourceFileName = "/Users/wgm/Documents/caBIG_R3/datasets/GO/gene_association.goa_human";
        Map<String, Set<String>> fMap = new HashMap<String, Set<String>>();
        Map<String, Set<String>> pMap = new HashMap<String, Set<String>>();
        FileReader fileReader = new FileReader(sourceFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = null;
        String[] tokens = null;
        String goId = null;
        Set<String> geneIds = null;
        Set<String> allGeneIds = new HashSet<String>();
        while ((line = reader.readLine()) != null) {
            tokens = line.split("\\t");
            goId = tokens[4];
            allGeneIds.add(tokens[1]);
            if (tokens[8].equals("F")) {
                geneIds = fMap.get(goId);
                if (geneIds == null) {
                    geneIds = new HashSet<String>();
                    fMap.put(goId, geneIds);
                }
                geneIds.add(tokens[1]);
            }
            else if (tokens[8].equals("P")) {
                geneIds = pMap.get(goId);
                if (geneIds == null) {
                    geneIds = new HashSet<String>();
                    pMap.put(goId, geneIds);
                }
                geneIds.add(tokens[1]);
            }
        }
        reader.close();
        fileReader.close();
        createValueMap(fValueMap, fMap);
        createValueMap(pValueMap, pMap);
    }
    
    private void createValueMap(final Map<String, Integer> valueMap, Map<String, Set<String>> map) {
        Set<String> geneIds = null;
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String id = it.next();
            geneIds = map.get(id);
            //value = (5.0f + (float)Math.log10((float)geneIds.size() / totalGeneNumber)) / 5.0f;
            valueMap.put(id, geneIds.size());
        }
    }
    
    public void generateUsageFrequence() throws IOException {
        Map<String, Integer> fValueMap = new HashMap<String, Integer>();
        Map<String, Integer> pValueMap = new HashMap<String, Integer>();
        createUsageFrequences(fValueMap, pValueMap);
        outputMap(fValueMap, "results/GOFOccurence.txt");
        outputMap(pValueMap, "results/GOPOccurence.txt");
    }
    
    private void outputMap(final Map<String, Integer> valueMap, 
                           String fileName) throws IOException {
        // Want to sort based on the values
        List<String> goIds = new ArrayList<String>(valueMap.keySet());
        Collections.sort(goIds, new Comparator<String>() {
            public int compare(String id1, String id2) {
                Integer value1 = valueMap.get(id1);
                Integer value2 = valueMap.get(id2);
                return value1.compareTo(value2);
            }
        });
        StringBuilder builder = new StringBuilder();
        for (String id : goIds) {
            builder.append(id).append("\t");
            builder.append(valueMap.get(id));
            builder.append("\n");
        }
        FileWriter fileWriter = new FileWriter(fileName);
        BufferedWriter writer = new BufferedWriter(fileWriter);
        writer.write(builder.toString());
        writer.close();
        fileWriter.close();
    }
    
    private Map<String, Integer> loadMap(String fileName) throws IOException {
        Map<String, Integer> map = new HashMap<String, Integer>();
        FileReader fileReader = new FileReader(fileName);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line = null;
        int index = 0;
        while ((line = bufferedReader.readLine()) != null) {
            index = line.indexOf("\t");
            map.put(line.substring(0, index), new Integer(line.substring(index + 1)));
        }
        return map;
    }
    
    /**
     * The following is the way how to discretize GO occurence data:
     * 1). Get all occurence values and sort them.
     * 2). Choose the total number of categories as 10
     * 3). Get the total number of the values (e.g. 158)
     * 4). Get how many values should be in one category (e.g. 16)
     * 5). Create categories by defing bin edges (e.g. 2-17, inclusive)
     * 6). Assign categories to all gene pairs by creating a new file.
     * No discretization is done for GO depth data. 
     * @throws IOException
     */
    public void discretizeOccurence() throws IOException {
        int NUMBER_OF_GROUPS = 10;
        String fileName = "ShareGOPOccurence.txt"; 
        //String fileName = "ShareGOFOccurence.txt";
        Map<String, Integer> occurences = new HashMap<String, Integer>();
        String line = null;
        int index = 0;
        String tmp;
        Integer occurence;
        FileReader fileReader = new FileReader("results/" + fileName);
        BufferedReader reader = new BufferedReader(fileReader);
        while ((line = reader.readLine()) != null) {
            index = line.lastIndexOf("\t");
            tmp = line.substring(index + 1);
            occurence = occurences.get(tmp);
            if (occurence == null) {
                occurences.put(tmp, 1);
            }
            else {
                occurences.put(tmp, occurence + 1);
            }
        }
        reader.close();
        fileReader.close();
        int[] values = new int[occurences.size()];
        index = 0;
        for (Iterator<String> it = occurences.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            values[index ++] = Integer.parseInt(key);
        }
        Arrays.sort(values);
        // Create categories
        List<int[]> groups = new ArrayList<int[]>();
        int valueInEachGroup = (int) Math.ceil((double)values.length / NUMBER_OF_GROUPS);
        System.out.println("ValuesInEachGroup: " + valueInEachGroup);
        int rightBound;
        for (int i = 0; i < values.length; i += valueInEachGroup) {
            int[] bounds = new int[2];
            bounds[0] = values[i];
            rightBound = i + valueInEachGroup - 1;
            if (rightBound > values.length - 1)
                rightBound = values.length - 1;
            bounds[1] = values[rightBound];
            groups.add(bounds);
        }
        // Check the bounds
        for (int[] bounds : groups) {
            System.out.printf("%d %d%n", bounds[0], bounds[1]);
        }
        System.out.println("Starting discreting...");
        discretizeOccurences(groups, fileName);
    }
    
    private void discretizeOccurences(List<int[]> groups, String fileName) throws IOException {
         // Input 
         FileReader fileReader = new FileReader("results/" + fileName);
         BufferedReader reader = new BufferedReader(fileReader);
         // Output
         FileWriter fileWriter = new FileWriter("results/" + fileName + ".dis");
         PrintWriter printWriter = new PrintWriter(fileWriter);
         String line = null;
         int index = 0;
         int value;
         while ((line = reader.readLine()) != null) {
             index = line.lastIndexOf("\t");
             value = Integer.parseInt(line.substring(index + 1));
             // Check which group this value is in
             for (int i = 0; i < groups.size(); i++) {
                 int[] bounds = groups.get(i);
                 if (value >= bounds[0] && value <= bounds[1]) {
                     printWriter.println(line.substring(0, index + 1) + i);
                 }
             }
         }
         printWriter.close();
         fileWriter.close();
         reader.close();
         fileReader.close();
    }
    
    public void analyzeOccurence() throws IOException {
        String[] fileNames = new String[] {
                //"ShareGOPOccurence.txt",
                "ShareGOFOccurence.txt"
        };
        Map<String, Integer> occurences = new HashMap<String, Integer>();
        String line = null;
        int index = 0;
        String tmp;
        Integer occurence;
        for (String name : fileNames) {
            FileReader fileReader = new FileReader("results/" + name);
            BufferedReader reader = new BufferedReader(fileReader);
            while ((line = reader.readLine()) != null) {
                index = line.lastIndexOf("\t");
                tmp = line.substring(index + 1);
                occurence = occurences.get(tmp);
                if (occurence == null) {
                    occurences.put(tmp, 1);
                }
                else {
                    occurences.put(tmp, occurence + 1);
                }
            }
        }
        List<String> sorted = new ArrayList<String>(occurences.keySet());
        Collections.sort(sorted, new Comparator<String>() {
            public int compare(String s1, String s2) {
                int i1 = Integer.parseInt(s1);
                int i2 = Integer.parseInt(s2);
                return i1 - i2;
            }
        });
        for (String s : sorted) {
            System.out.printf("%s %d%n", s, occurences.get(s));
        }
    }
    
    public void analyzeDepth() throws IOException {
        String[] fileNames = new String[] {
                "ShareGOPDepth.txt",
                "ShareGOFDepth.txt"
        };
        Set<String> values = new HashSet<String>();
        String line = null;
        int index = 0;
        String tmp;
        for (String name : fileNames) {
            FileReader fileReader = new FileReader("results/" + name);
            BufferedReader reader = new BufferedReader(fileReader);
            while ((line = reader.readLine()) != null) {
                index = line.lastIndexOf("\t");
                tmp = line.substring(index + 1);
                values.add(tmp);
            }
        }
        List<Integer> sorted = new ArrayList<Integer>();
        for (String s : values)
            sorted.add(new Integer(s));
        Collections.sort(sorted);
        System.out.println("Total Values: " + sorted.size());
        for (int s : sorted)
            System.out.println(s);
    }
    
    public void calculateInformationContents() throws IOException {
        // Load all occurences first
        Map<String, Integer> occurenceMap = new HashMap<String, Integer>();
        String fileName = RESULT_DIR + "GOFOccurence.txt";
        String outputFileName = RESULT_DIR + "GOMFInfoContent.txt";
        String rootId = "GO:0003674"; // molecular_function
        //String fileName = RESULT_DIR + "GOPOccurence.txt";
        //String outputFileName = RESULT_DIR + "GOBPInfoContent.txt";
        //String rootId = "GO:0008150"; // biological processes
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            String id = line.substring(0, index);
            Integer occurence = Integer.parseInt(line.substring(index + 1));
            occurenceMap.put(id, occurence);
        }
        fu.close();
        // Have to re-calculate the occurence based on the DAG
        Map<String, GOGraphNode> id2Node = generateGOTree();
        // Now it is possible to re-calculate the occurences
        final Map<String, Integer> newOccMap = new HashMap<String, Integer>();
        Integer occurence = null;
        Set<GOGraphNode> children;
        String id;
        GOGraphNode node;
        for (Iterator<String> it = id2Node.keySet().iterator(); it.hasNext();) {
            id = it.next();
            occurence = occurenceMap.get(id);
            node = id2Node.get(id);
            children = node.getChildren();
            if (children != null && children.size() > 0) {
                occurence = calculateOccurenceWithDescendents(node, occurenceMap);
            }
            if (occurence != null && occurence.intValue() != 0)
                newOccMap.put(id, occurence);
        }
        // Calculate information contents
        Map<String, Double> id2InfoContent = new HashMap<String, Double>();
        int maximum = newOccMap.get(rootId);
        double infoContent;
        for (Iterator<String> it = newOccMap.keySet().iterator(); it.hasNext();) {
            id = it.next();
            occurence = newOccMap.get(id);
            infoContent = (double) occurence / maximum;
            id2InfoContent.put(id, infoContent);
        }
        fu.setOutput(outputFileName);
        for (Iterator<String> it = id2InfoContent.keySet().iterator(); it.hasNext();) {
            id = it.next();
            infoContent = id2InfoContent.get(id);
            fu.printLine(id + "\t" + infoContent);
        }
        fu.close();
    }
    
    public Map<String, Double> loadInformationContents(String fileName) throws IOException {
        Map<String, Double> id2InfoContent = new HashMap<String, Double>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id2InfoContent.put(line.substring(0, index),
                               Double.valueOf(line.substring(index + 1)));
        }
        return id2InfoContent;
    }
    
    public Map<String, Double> loadBPInformationContents() throws IOException {
        String fileName = "datasets" + File.separator + "GO" + File.separator + "GOBPInfoContent.txt";
        return loadInformationContents(fileName);
    }
    
    /**
     * This help method is used to calculate the occurence by considering the passed node's descendents.
     * @param node
     * @param occurenceMap
     * @return
     */
    private int calculateOccurenceWithDescendents(GOGraphNode node, Map<String, Integer> occurenceMap) {
        int occurence = 0;
        Set<GOGraphNode> current = new HashSet<GOGraphNode>();
        current.add(node);
        Set<GOGraphNode> next = new HashSet<GOGraphNode>();
        Integer originalOccurence;
        while (current.size() > 0) {
            for (GOGraphNode currentNode : current) {
                originalOccurence = occurenceMap.get(currentNode.getId());
                if (originalOccurence == null)
                    originalOccurence = 0; 
                occurence += originalOccurence;
                if (currentNode.getChildren() != null)
                    next.addAll(currentNode.getChildren());
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
        return occurence;
    }
    
    // Used to cache values
    private Map<String, Double> cachedValues = new HashMap<String, Double>();
    
    private double calculateSemanticSimilarity(String id1, 
                                               String id2,
                                               Map<String, Double> id2InfoContent,
                                               Map<String, GOGraphNode> id2Node) {
        int compare = id1.compareTo(id2);
        String key = null;
        if (compare < 0)
            key = id1 + " " + id2;
        else
            key = id2 + " " + id1;
        Double rtn = cachedValues.get(key);
        if (rtn != null)
            return rtn;
        //Double rtn;
        if (compare == 0) { // both terms are the same
            rtn = (-Math.log(id2InfoContent.get(id1)));
        }
        else {
            Set<String> ancestors1 = grepAncestors(id2Node.get(id1));
            Set<String> ancestors2 = grepAncestors(id2Node.get(id2));
            ancestors2.retainAll(ancestors1);
            // Find the minimum
            double min = Double.MAX_VALUE;
            for (String id : ancestors2) {
                Double tmp = id2InfoContent.get(id);
                if (tmp == null) {
                     continue;
                }
                if (tmp < min)
                    min = tmp;
            }
            // Use the minimum
            // Min should not be possible to be Double.MAX_VALUE
            if (min == Double.MAX_VALUE)
                //rtn = 0.0d;
                throw new IllegalStateException(id1 + " " + id2 + ": have empty information content");
            else
                rtn = (-Math.log(min));
        }
        cachedValues.put(key, rtn);
        return rtn;
    }
    
    public Double calculateSemanticSimilarity(Set<String> terms1,
                                              Set<String> terms2,
                                              Map<String, Double> id2InfoContent,
                                              Map<String, GOGraphNode> id2Node) {
        if (terms1 == null || terms1.size() == 0)
            return null;
        if (terms2 == null || terms2.size() == 0)
            return null;
        // For counting
        int c = 0;
        double totalSimilarity = 0.0d;
        double similarity = 0.0d;
        double min;
        for (String term1 : terms1) {
            for (String term2 : terms2) {
                c ++;
                similarity = calculateSemanticSimilarity(term1, term2, id2InfoContent, id2Node);
                totalSimilarity += similarity;
            }
        }
        return totalSimilarity / c;
    }
    
    // Used to cache ancestors to improve the performance
    private Map<String, Set<String>> id2Ancestor = new HashMap<String, Set<String>>();
    /**
     * A helper method to grep all ancestors for the passed GOGraphNode. The id of the
     * passed GOGraphNode is included in the returnred set.
     * @param node
     * @return
     */
    private Set<String> grepAncestors(GOGraphNode node) {
        // Check if it is cached
        Set<String> ids = id2Ancestor.get(node.getId());
        if (ids != null)
            return new HashSet<String>(ids);
        ids = new HashSet<String>();
        Set<GOGraphNode> current = new HashSet<GOGraphNode>();
        current.add(node);
        Set<GOGraphNode> next = new HashSet<GOGraphNode>();
        while (current.size() > 0) {
            for (GOGraphNode tmp : current) {
                ids.add(tmp.getId());
                if (tmp.getParents() != null)
                    next.addAll(tmp.getParents());
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
        id2Ancestor.put(node.getId(), ids);
        return new HashSet<String>(ids);
    }
    
    public static class GOGraphNode {
        Set<GOGraphNode> parents;
        Set<GOGraphNode> children;
        String id;
        
        public GOGraphNode() {
            
        }
        
        public void setId(String id) {
            this.id = id;
        }
        
        public String getId() {
            return this.id;
        }
        
        public void addParent(GOGraphNode node) {
            if (parents == null)
                parents = new HashSet<GOGraphNode>();
            parents.add(node);
        }
        
        public Set<GOGraphNode> getParents() {
            return this.parents;
        }
        
        public void addChild(GOGraphNode node) {
            if (children == null)
                children = new HashSet<GOGraphNode>();
            children.add(node);
        }
        
        public Set<GOGraphNode> getChildren() {
            return children;
        }
    }
   
}
