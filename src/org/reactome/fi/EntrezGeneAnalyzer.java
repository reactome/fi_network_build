package org.reactome.fi;

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

import org.junit.Test;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 *  Entrez gene is the hub. This class is written to process the data files obtained
 *  from it.
 * @author feng Oct 18, 2007
 *
 */
public class EntrezGeneAnalyzer {
	private String dirName = "datasets/ncbi/080409/";
	private FileUtility fu = new FileUtility();
	// Cache the mapped synonym to gene names
	private Map<String, String> synonymToName;
	private Set<String> allGenes;
    
    public EntrezGeneAnalyzer() {
	}

    /**
     * This method is used to load all standard gene names.
     * @return
     */
    public Set<String> getAllGenes() throws IOException {
        if (synonymToName == null)
            synonymToName = getGeneSynonymToNameMap();
        return new HashSet<String>(synonymToName.values());
    }
    
    /**
	 * 
	 * ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
	 * 
	 * a sample line is like this:
	 * 9	1246500	repA1	pLeuDn_01	-	-	-	-	putative replication-associated protein	protein-coding	-	-	-	-	20071004
	 */
	public ArrayList<Map<String,String>> loadGeneInfoFile() throws Exception{
		ArrayList<Map<String,String>> resultMap = new ArrayList<Map<String,String>>();
		Map<String,String> name2IDMap = new HashMap<String,String>();
		Map<String,String> locus2IDMap = new HashMap<String,String>();
		Map<String,String> id2NameMap = new HashMap<String,String>();
		FileUtility fu = new FileUtility();
		Set<String> lines = fu.loadInteractions("data/EntrezGene/gene_info");
		for(String line : lines){
			String[] tokens = line.split("\t");
			//name2IDMap.put(tokens[2],tokens[1]);
			//locus2IDMap.put(tokens[3], tokens[1]);
			id2NameMap.put(tokens[1], tokens[11]);
		}
		resultMap.add(name2IDMap);
		resultMap.add(locus2IDMap);
		resultMap.add(id2NameMap);
		return resultMap;
	}
	
	@Test
	public void convertGeneIdsToSymbols() throws Exception {
	    Map<String, String> rtn = loadGeneIdToSymbol();
	    String[] ids = new String[] {
	            "2335",
	            "2736",
	            "3452",
	            "3838",
	            "4609",
	            "4846",
	            "672",
	            "51701",
	            "558",
	            "80324",
	            "7157"
	    };
	    for (String id : ids) {
	        String symbol = rtn.get(id);
	        System.out.println(id + "\t" + symbol);
	    }
	    for (String id : ids) 
	        System.out.println(rtn.get(id));
	}
	
	public Map<String, String> loadGeneIdToSymbol() throws IOException {
	    Map<String, String> rtn = new HashMap<String, String>();
	    String fileName = dirName + "Homo_sapiens.gene_info";
	    fu.setInput(fileName);
	    String line = fu.readLine();
	    while ((line = fu.readLine()) != null) {
	        String[] tokens = line.split("\t");
	        String geneId = tokens[1];
	        String symbol = tokens[2];
	        rtn.put(geneId, symbol);
	    }
	    fu.close();
	    return rtn;
	}
	
	/**
	 * This method is used to check if one gene id is corresponding to one gene name.
	 * @throws IOException
	 */
	@Test
	public void checkGeneIdsAndNames() throws IOException {
	    String fileName = dirName + "Homo_sapiens.gene_info";
	    fu.setInput(fileName);
	    String line = fu.readLine();
	    List<String> geneIds = new ArrayList<String>();
	    List<String> geneNames = new ArrayList<String>();
	    while ((line = fu.readLine()) != null) {
	        String[] tokens = line.split("\t");
	        String geneId = tokens[1];
	        if (!geneId.equals("-"))
	            geneIds.add(geneId);
	        String geneName = tokens[2];
	        if (!geneName.equals("-"))
	            geneNames.add(geneName);
	    }
	    System.out.println("Gene Ids in list: " + geneIds.size());
	    System.out.println("Gene Names in list:  " + geneNames.size());
	    // Change to set
	    Set<String> geneIdSet = new HashSet<String>(geneIds);
	    Set<String> geneNameSet = new HashSet<String>(geneNames);
	    System.out.println("Gene ids in set: " + geneIdSet.size());
	    System.out.println("Gene names in set: " + geneNameSet.size());
	    // Find what projects have been repeated
	    Map<String, Integer> geneNameToNumber = InteractionUtilities.countTermUsageInList(geneNames);
	    for (Iterator<String> it = geneNameToNumber.keySet().iterator(); it.hasNext();) {
	        String geneName = it.next();
	        Integer count = geneNameToNumber.get(geneName);
	        if (count > 1)
	            System.out.println(geneName);
	    }
	}
	
	@Test
	public void testGetGeneSynonymToNameMap() throws IOException {
	    Map<String, String> map = getGeneSynonymToNameMap();
	    System.out.println("Size of map: " + map.size());
	    String inFileName = "/Users/wgm/Documents/wgm/OHSU/TargetGenesList.txt";
	    Set<String> totalGenes = fu.loadInteractions(inFileName);
	    System.out.println("Total genes: " + totalGenes.size());
	    List<String> geneList = new ArrayList<String>(totalGenes);
	    Collections.sort(geneList);
	    int count = 0;
	    List<String> mapped = new ArrayList<String>();
	    int notTheSameCount = 0;
	    for (String gene : geneList) {
	        String name = map.get(gene);
	        if (name == null) {
	            name = gene;
	            count ++;
	        }
	        //System.out.println(gene + "\t" + name);
	        mapped.add(name);
	        if (!name.equals(gene)) {
	            //System.out.println("Not the same!");
	            System.out.println(gene + "\t" + name);
	            notTheSameCount ++;
	        }
	    }
	    System.out.println("Cannot mapped: " + count);
	    System.out.println("Not the same: " + notTheSameCount);
	    Collections.sort(mapped);
	    for (String name : mapped)
	        System.out.println(name);
	}
	
    /**
     * This method is used to map all gene names to standard names based
     * on NCBI gene info file.
     * @return
     */
    public Set<String> normalizeGeneNames(Collection<String> geneNames) throws IOException {
        if (synonymToName == null) {
            synonymToName = getGeneSynonymToNameMap();
            allGenes = new HashSet<String>(synonymToName.values());
        }
        Set<String> rtn = new HashSet<String>();
        for (String name : geneNames) {
            String tmp = synonymToName.get(name);
//            System.out.println(name + "\t" + tmp);
            if (tmp == null) {
                rtn.add(name);
                //System.out.println("cannot mapped: " + name);
            }
            else
                rtn.add(tmp);
        }
        return rtn;
    }
    
    /**
     * This method is used to map all names used in FIs to standard names based on
     * NCBI gene info file.
     * @param fis
     * @return
     * @throws IOException
     */
    public Set<String> normalizeFIInNames(Set<String> fis) throws IOException {
        if (synonymToName == null) {
            synonymToName = getGeneSynonymToNameMap();
            allGenes = new HashSet<String>(synonymToName.values());
        }
        int index = 0;
        Set<String> rtn = new HashSet<String>();
        String name1, name2;
        String mapped1, mapped2;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            name1 = fi.substring(0, index);
            name2 = fi.substring(index + 1);
            mapped1 = getMappedName(name1,
                                    synonymToName,
                                    allGenes);
            mapped2 = getMappedName(name2, 
                                    synonymToName,
                                    allGenes);
            int compare = mapped1.compareTo(mapped2);
            if (compare < 0)
                rtn.add(mapped1 + "\t" + mapped2);
            else if (compare > 0)
                rtn.add(mapped2 + "\t" + mapped1);
        }
        return rtn;
    }
    
    private String getMappedName(String name,
                                 Map<String, String> synonymToName,
                                 Set<String> allGenes) {
        name = name.toUpperCase();
        // the name itself is the standard symbol. See "PIK3CD" and "CLSTN1" for an example
        if (allGenes.contains(name))
            return name;
        String mapped = synonymToName.get(name);
        if (mapped != null)
            return mapped;
        return name;
    }
	
	/**
	 * This method is used to generate a map from a gene synonym to a gene name. In this map,
	 * gene name is also listed too. In other words, a gene name is mapped to gene name. This 
	 * map returns to human genes only. 
	 * @return
	 * @throws IOException
	 */
	public Map<String, String> getGeneSynonymToNameMap() throws IOException {
	    String fileName = dirName + "Homo_sapiens.gene_info";
	    fu.setInput(fileName);
	    String line = fu.readLine();
	    Map<String, String> synonymToName = new HashMap<String, String>();
	    while ((line = fu.readLine()) != null) {
	        String[] tokens = line.split("\t");
	        String name = tokens[2];
	        synonymToName.put(name, name);
	        String synonyms = tokens[4];
	        if (synonyms.equals("-"))
	            continue;
	        tokens = synonyms.split("\\|");
	        for (String synonym : tokens)
	            synonymToName.put(synonym, 
	                              name);
	    }
	    fu.close();
	    return synonymToName;
	}
	
	public Map<String, String> loadRNAAccessionToGeneId() throws IOException {
	    String srcFileName = "datasets/ncbi/080311/human_gene2accession.txt";
	    fu.setInput(srcFileName);
	    Map<String, String> accToGeneId = new HashMap<String, String>();
	    String line = fu.readLine();
	    int index = 0;
	    int duplicated = 0;
	    while ((line = fu.readLine()) != null) {
	        String[] tokens = line.split("\t");
	        String acc = tokens[3];
	        if (acc.equals("-"))
	            continue;
	        String geneId = tokens[1];
	        if (geneId.equals("-"))
	            continue;
	        // Remove version number
	        index = acc.indexOf(".");
	        acc = acc.substring(0, index);
	        if (accToGeneId.containsKey(acc)) {
	            String tmp = accToGeneId.get(acc);
	            if (!tmp.equals(geneId)) {
//	                System.out.println("Duplicated accession with differenet geneId: " + line);
	                duplicated ++;
	            }
	        }
	        accToGeneId.put(acc, geneId);
	    }
	    System.out.println("Duplicated accessions: " + duplicated);
	    return accToGeneId;
	}
	
	@Test
	public void testLoadRNAAccessionToGeneId() throws IOException {
	    Map<String, String> accToGeneId = loadRNAAccessionToGeneId();
	    System.out.println("Total accessions: " + accToGeneId.size());
	    Map<String, String> geneIdToSymbol = loadGeneIdToSymbol();
	    System.out.println("Total gene ids to symbols: " + geneIdToSymbol.size());
	}
	
	/**
	 * Use this method to load RNA accession number to gene symbols.
	 * @return
	 * @throws IOException
	 */
	public Map<String, String> loadRNAAccToGeneSymbol() throws IOException {
	    Map<String, String> accToGeneId = loadRNAAccessionToGeneId();
        System.out.println("Total accessions: " + accToGeneId.size());
        Map<String, String> geneIdToSymbol = loadGeneIdToSymbol();
        System.out.println("Total gene ids to symbols: " + geneIdToSymbol.size());
        Map<String, String> accToSymbol = new HashMap<String, String>();
        for (String acc : accToGeneId.keySet()) {
            String geneId = accToGeneId.get(acc);
            String symbol = geneIdToSymbol.get(geneId);
            if (symbol != null)
                accToSymbol.put(acc, symbol);
        }
        return accToSymbol;
	}
	
	
	/**
	 * This method is used to grep human genes to accession information from the original
	 * download gene2accession file.
	 * @throws IOException
	 */
	@Test
	public void grepHumanGeneToAccession() throws IOException {
	    String dirName = "datasets/ncbi/080311/";
	    String srcFileName = dirName + "gene2accession";
	    String destFilename = dirName + "human_gene2accession.txt";
	    fu.setInput(srcFileName);
	    fu.setOutput(destFilename);
	    String line = fu.readLine(); // header line
	    fu.printLine(line);
	    while ((line = fu.readLine()) != null) {
	        if (line.startsWith("9606")) // Check if a line is for human
	            fu.printLine(line);
	    }
	    fu.close();
	}
	
}
