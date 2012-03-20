/*
 * Created on Jun 21, 2006
 *
 */
package org.reactome.fi;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import org.junit.Test;
import org.reactome.fi.util.FileUtility;

public class OrthologyAnalyzer {
    private final String ORTHOMCL_SRC_DIR = "data/orthomcl/";
    private final String INPARANOID_SRC_DIR = "/Users/guanming/Documents/caBIG_R3/datasets/inparanoid/";
    private final String RESULT_DIR = "results/orthomcl/";

    public OrthologyAnalyzer() {   
    }
    
    public void generatePairwiseMapFromInParanoid() throws IOException {
        String srcFileName = INPARANOID_SRC_DIR + "sqltable.modSC-ensHS.txt";
        //1   3010    modSC   1.000   S000001208  100%
        //1   3010    ensHS   1.000   ENSP00000304350 100%
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        String line = null;
        String[] tokens = null;
        List<String> scIds = new ArrayList<String>();
        List<String> hsIds = new ArrayList<String>();
        String clusterId = "0";
        //Output
        String outputFileName = RESULT_DIR + "sce2hsaInParanoid.txt";
        // Have to generate the output without keeping in a map
        // since one Sce protein can map to several different hsa proteins and vice versa.
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outputFileName);
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            // Check cluster id
            if (!tokens[0].equals(clusterId)) {
                // Need to generate pairwise 
                for (String sc : scIds) {
                    for (String hs : hsIds) {
                        outFu.printLine(sc + "\t" + hs);
                    }
                }
                // New cluster
                clusterId = tokens[0];
                scIds.clear();
                hsIds.clear();
            }
            if (tokens[2].equals("modSC"))
                scIds.add(tokens[4]);
            else if (tokens[2].equals("ensHS"))
                hsIds.add(tokens[4]);
        }
        fu.close();
        outFu.close();
    }
    
    private Map<String, List<String>> loadRedundantMap(String fileName) throws IOException {
        Map<String, List<String>> map = new HashMap<String, List<String>>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        int index = 0;
        String line;
        String id1, id2;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id1 = line.substring(0, index);
            id2 = line.substring(index + 1);
            List<String> list = map.get(id1);
            if (list == null) {
                list = new ArrayList<String>();
                map.put(id1, list);
            }
            list.add(id2);
        }
        fu.close();
        return map;
    }
    
    public void convertInteractionForHuman() throws IOException {
        FileUtility fu = new FileUtility();
        String[] fileNames = new String[]{
                RESULT_DIR + "yeastNonY2HInteraction.txt",
                RESULT_DIR + "yeastY2HInteraction.txt"
                //RESULT_DIR + "flyNonY2HInteraction.txt",
                //RESULT_DIR + "flyY2HInteraction.txt"
                //RESULT_DIR + "wormNonY2HInteraction.txt",
                //RESULT_DIR + "wormY2HInteraction.txt"
        };
        Map<String, List<String>> toHsa = loadRedundantMap(RESULT_DIR + "sce2hsaOrthoMCLUpdated.txt");
        //Map<String, List<String>> toHsa = loadOrthoMap(RESULT_DIR + "sce2hsaOrthoMCL.txt");
        //Map<String, List<String>> toHsa = loadRedundantMap(RESULT_DIR + "dme2hsaOrthoMCL.txt");
        //Map<String, List<String>> toHsa = loadOrthoMap(RESULT_DIR + "cel2hsaOrthoMCL.txt");
        String line;
        int index;
        FileUtility outFu = new FileUtility();
        String partner1 = null;
        List<String> list1 = null;
        String partner2 = null;
        List<String> list2 = null;
        int c = 0;
        int compareResult = 0;
        Set<String> set = new HashSet<String>();
        for (String in : fileNames) {
            set.clear();
            // Create a output Fu
            index = in.lastIndexOf("/");
            int index1 = in.lastIndexOf(".");
            String outFileName = RESULT_DIR + in.substring(index, index1) + "Hsa.txt";
            outFu.setOutput(outFileName);
            fu.setInput(in);
            while ((line = fu.readLine()) != null) {
                index = line.indexOf(" ");
                partner1 = line.substring(0, index);
                list1 = toHsa.get(partner1);
                if (list1 == null || list1.size() == 0)
                    continue;
                partner2 = line.substring(index + 1);
                list2 = toHsa.get(partner2);
                if (list2 == null || list2.size() == 0)
                    continue;
                for (String s1 : list1) {
                    for (String s2 : list2) {
                        compareResult = s1.compareTo(s2);
                        if (compareResult < 0)
                            set.add(s1 + " " + s2);
                        else
                            set.add(s2 + " " + s1);
                        c ++;
                    }
                }
            }
            for (String s : set)
                outFu.printLine(s);
            outFu.close();
            fu.close();
            System.out.println(in + ": " + c + ", " + set.size());
        }
    }
    
    private Map<String, String> loadFlybaseMap() throws IOException {
        String inFile = ORTHOMCL_SRC_DIR + "FBgn.acode";
        FileUtility fu = new FileUtility();
        fu.setInput(inFile);
        String line = null;
        //ID|FBgn0025724
        //SYN|CG6699
        String flybaseId = null;
        String geneName = null;
        Map<String, String> map = new HashMap<String, String>();
        boolean inIn = false;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("GENR")) {
                line = fu.readLine(); // should be "{";
                while (true) {
                    line = fu.readLine();
                    if (line.startsWith("}"))
                        break;
                    if (line.startsWith("ID|FBgn")) {
                        flybaseId = line.substring(3);
                    }
                    else if (line.startsWith("SYN|CG")) {
                        int index = line.indexOf(" ", 4);
                        if (index < 0)
                            geneName = line.substring(4);
                        else
                            geneName = line.substring(4, index);
                        map.put(geneName, flybaseId);
                        while (true) {
                            line = fu.readLine();
                            if (line.startsWith("|CG")) {
                                index = line.indexOf(" ", 4);
                                if (index < 0)
                                    geneName = line.substring(1);
                                else
                                    geneName = line.substring(1, index);
                                map.put(geneName, flybaseId);
                            }
                            else
                                break;
                        }
                    }
                    // In case this line is from the above else branch
                    if (line.startsWith("NAM|CG")) {
                        geneName = line.substring(4);
                        map.put(geneName, flybaseId);
                    }
                    //map.put(geneName, flybaseId);
                }
            }
        }
        System.out.println("Total Mapping: " + map.size());
        return map;
    }
    
    
    @Test
    public void generatePairwiseMap() throws IOException {
     //   Map<String, String> flybaseMap = loadFlybaseMap();
//        FileUtility fu1 = new FileUtility();
//        fu1.exportMap(flybaseMap, "results/flybasemap.txt");
//        if (true)
//            return;
        String mclIdFile = ORTHOMCL_SRC_DIR + "BAE_geneid_anno.txt";
        Map<String, String> ortho2id = new HashMap<String, String>();
        FileUtility fu = new FileUtility();
        fu.setInput(mclIdFile);
        String line = null;
        int index = 0;
        int index1 = 0;
        String ortho;
        String id;
        int c = 0;
        int total = 0;
        Set<String> missed = new HashSet<String>();
        System.out.println("sdfsdf");
        while ((line = fu.readLine()) != null) {
          /*  if (line.startsWith("sce")) {
                //sce2    YAL002W VPS8 SGDID:S0000002, Chr I from 143709-147533, Verified ORF
                index = line.indexOf("\t");
                ortho = line.substring(0, index);
                index = line.indexOf("SGDID:");
                index1 = line.indexOf(",", index);
                id = line.substring(index + 6, index1);
                ortho2id.put(ortho, id);
            }
        
            else if (line.startsWith("dme")) {
                //dme100  CG4280-PB pep:known chromosome:BDGP3.2.1:2L:449919:454688:-1 gene:CG4280 transcript:CG4280-RB
                index = line.indexOf("\t");
                ortho = line.substring(0, index);
                index = line.indexOf("gene:");
                index1 = line.indexOf(" ", index);
                id = line.substring(index + 5, index1);
                // Use Fng##### id
                String flybaseId = flybaseMap.get(id);
                if (flybaseId == null) {
                    c++;
                    missed.add(id);
                    continue;
                }
                //if (flybaseId == null)
                //    throw new IllegalStateException("Cannot find flybase id for: " + id);
                ortho2id.put(ortho, flybaseId);
                total ++;
            }*/
           if (line.startsWith("hsa")) {
        	   //hsa33867	ENSP00000351823 Y chromosome spermatogenesis candidate protein (RBM) pseudogene 
               StringTokenizer st = new StringTokenizer(line);
               ortho = "";
               id = "";
                while(st.hasMoreTokens()){
                	String token = st.nextToken();
                	if(token.startsWith("hsa"))
                		ortho = token;
                	if(token.startsWith("ENSP"))
                		id = token;
                	if(ortho.startsWith("hsa") && id.startsWith("ENSP")) break;
                }
                ortho2id.put(ortho, id);
            }
           //Added by xin for the Zhong's cel conservation analysis oct 23,2007
           else if (line.startsWith("cel")) {
               //cel100  B0035.1a 
        	   //The data file has been changed!
        	   StringTokenizer st = new StringTokenizer(line);
               ortho = "";
               id = "";
                while(st.hasMoreTokens()){
                	String token = st.nextToken();
                	if(token.startsWith("cel"))
                		ortho = token;
                	if(token.contains("."))
                		id = token;
                	if(ortho.startsWith("hsa") && id.contains(".")) break;
                }
                ortho2id.put(ortho, id);
           }
        /*   else if(line.startsWith("mmu")) {
         	   //mmu3	ENSMUSP00000027032 Oxygen-regulated protein 1 (Retinitis pigmentosa RP1 protein homolog). 
        	   index = line.indexOf("\t");
               ortho = line.substring(0, index);
               id = line.substring(index+1,index+19);
               ortho2id.put(ortho, id);  
           }
           else if(line.startsWith("gga")) {
         	   //gga143	ENSGALP00000013235 Ensembl Family: SERINE/THREONINE KINASE 23 EC_
        	   index = line.indexOf("\t");
               ortho = line.substring(0, index);
               id = line.substring(index+1,index+19);
               ortho2id.put(ortho, id);  
           }
           else if(line.startsWith("rno")) {
         	   //rno132	ENSRNOP00000044949 similar to 60S RIBOSOMAL PROTEIN L27A 
        	   index = line.indexOf("\t");
               ortho = line.substring(0, index);
               id = line.substring(index+1,index+19);
               ortho2id.put(ortho, id);  
           }
           else if(line.startsWith("sce")) {
         	   //sce6	YAL007C ERP2, Protein that forms a heterotrimeric complex with Erp1p
        	   index = line.indexOf("\t");
               ortho = line.substring(0, index);
               if(index+ 8 >= line.length()) continue;
               id = line.substring(index+1,index+8);
               ortho2id.put(ortho, id);  
           }
           else if(line.startsWith("dme")) {
         	   //dme18	CG2674-PE S-adenosylmethionine synthetase (EC 2.5.1.6) (Methionine adenosyltransferase) (AdoMet synthetase). 
        	   index = line.indexOf("\t");
               ortho = line.substring(0, index);
               if(index+ 10 >= line.length()) continue;
               id = line.substring(index+1,index+10);
               ortho2id.put(ortho, id);  
           }*/
        }
        fu.close();
        System.out.println("Total: " + ortho2id.size());
        System.out.println("Missed mapping: " + c);
        for (String s : missed)
            System.out.println(s);
        generatePairwiseMap(ortho2id);
    }

    private void generatePairwiseMap(Map<String, String> ortho2Id) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setInput(ORTHOMCL_SRC_DIR + "humanOrthoMCL_10232007.txt");
        // These files will be generated
        String sce2hsaFile = ORTHOMCL_SRC_DIR + "sce2hsaOrthoMCL.txt";
        FileUtility sce2hsaFu = new FileUtility();
        sce2hsaFu.setOutput(sce2hsaFile);
        List<String> sceIDs = new ArrayList<String>();
        String cel2hsaFile = ORTHOMCL_SRC_DIR + "cel2hsaOrthoMCL.txt";
        FileUtility cel2hsaFu = new FileUtility();
        cel2hsaFu.setOutput(cel2hsaFile);
        List<String> celIDs = new ArrayList<String>();
        String dme2hsaFile = ORTHOMCL_SRC_DIR + "dme2hsaOrthoMCL.txt";
        FileUtility dme2hsaFu = new FileUtility();
        dme2hsaFu.setOutput(dme2hsaFile);
        List<String> dmeIDs = new ArrayList<String>();
        String mmu2hsaFile = ORTHOMCL_SRC_DIR + "mmu2hsaOrthoMCL.txt";
        FileUtility mmu2hsaFu = new FileUtility();
        mmu2hsaFu.setOutput(mmu2hsaFile);
        Set<String> mmuIDs = new HashSet<String>();
        String rno2hsaFile = ORTHOMCL_SRC_DIR + "rno2hsaOrthoMCL.txt";
        FileUtility rno2hsaFu = new FileUtility();
        rno2hsaFu.setOutput(rno2hsaFile);
        Set<String> rnoIDs = new HashSet<String>();
        String gga2hsaFile = ORTHOMCL_SRC_DIR + "gga2hsaOrthoMCL.txt";
        FileUtility gga2hsaFu = new FileUtility();
        gga2hsaFu.setOutput(gga2hsaFile);
        Set<String> ggaIDs = new HashSet<String>();
        
        // For human
        Set<String> hsaIDs = new HashSet<String>();
        
        String[] tokens = null;
        String line = null;
        String orthoId;
        String id;
        Set<String> sce2hsaSet = new HashSet<String>();
        Set<String> dme2hsaSet = new HashSet<String>();
        Set<String> cel2hsaSet = new HashSet<String>();
        Set<String> mmu2hsaSet = new HashSet<String>();
        Set<String> rno2hsaSet = new HashSet<String>();
        Set<String> gga2hsaSet = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            hsaIDs.clear();
            sceIDs.clear();
            celIDs.clear();
            dmeIDs.clear();
            mmuIDs.clear();
            rnoIDs.clear();
            ggaIDs.clear();
            tokens = line.split("\t");
            for (int i = 1; i < tokens.length; i++) {
                orthoId = tokens[i];
                id = ortho2Id.get(orthoId);
                if (id == null)
                    continue;
                if (orthoId.startsWith("hsa"))
                    hsaIDs.add(id);
                else if (orthoId.startsWith("sce"))
                    sceIDs.add(id);
                else if (orthoId.startsWith("cel"))
                    celIDs.add(id);
                else if (orthoId.startsWith("dme"))
                    dmeIDs.add(id);
                else if (orthoId.startsWith("mmu"))
                    mmuIDs.add(id);
                else if (orthoId.startsWith("rno"))
                    rnoIDs.add(id);
                else if (orthoId.startsWith("gga"))
                    ggaIDs.add(id);
            }
            if (hsaIDs.size() == 0)
                continue; // Nothing to generated
            for (String hsaId : hsaIDs) {
                for (String sceId : sceIDs) 
                    sce2hsaSet.add(sceId + "\t" + hsaId);
                for (String celId : celIDs) 
                    cel2hsaSet.add(celId + "\t" + hsaId);
                for (String dmeId : dmeIDs)
                    dme2hsaSet.add(dmeId + "\t" + hsaId);
                for (String mmuId : mmuIDs) 
                    mmu2hsaSet.add(mmuId + "\t" + hsaId);
                for (String rnoId : rnoIDs) 
                    rno2hsaSet.add(rnoId + "\t" + hsaId);
                for (String ggaId : ggaIDs)
                    gga2hsaSet.add(ggaId + "\t" + hsaId);
            }
        }
        for (String pair : sce2hsaSet)
            sce2hsaFu.printLine(pair);
        sce2hsaFu.close();
        for (String pair : cel2hsaSet)
            cel2hsaFu.printLine(pair);
        cel2hsaFu.close();
        for (String pair : dme2hsaSet)
            dme2hsaFu.printLine(pair);
        dme2hsaFu.close();
        for (String pair : mmu2hsaSet)
            mmu2hsaFu.printLine(pair);
        mmu2hsaFu.close();
        for (String pair : rno2hsaSet)
            rno2hsaFu.printLine(pair);
        rno2hsaFu.close();
        for (String pair : gga2hsaSet)
            gga2hsaFu.printLine(pair);
        gga2hsaFu.close();
        fu.close();
    }
    
    @Test
    public void extractOrthoMCLForHuman() throws IOException {
        String mclFile = ORTHOMCL_SRC_DIR + "all_orthomcl.out";
        StringBuilder builder = new StringBuilder();
        FileUtility fu = new FileUtility();
        fu.setInput(mclFile);
        String line = null;
        String[] tokens = null;
        List<String> entries = new ArrayList<String>();
        boolean hasHuman = false;
        int index = 0;
        String clusterId = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf(":");
            clusterId = line.substring(0, index);
            line = line.substring(index + 1);
            tokens = line.split(" ");
            entries.clear();
            hasHuman = false;
            // The first token is cluster id
            // check if human (hsa) is in this cluster
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].startsWith("hsa")) {
                    hasHuman = true;
                    break;
                }
            }
            /*
             * Modified by xin @ aug.16.2007
             */
            if (!hasHuman)
                continue; // Don't do anything
            for (int i = 1; i < tokens.length; i++) {
                if (tokens[i].startsWith("hsa"))
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("sce")) // yeast
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("cel")) // C. elegans
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("dme")) // Fly
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("mmu")) // Mus musculus
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("rno")) // Rattus norvegicus
                    entries.add(tokens[i]);
                else if (tokens[i].startsWith("gga")) // Gallus gallus
                    entries.add(tokens[i]);
            }
            // Have to make sure there are at least some genes from other genomes
            boolean nonHuman = false;
            for (String entry : entries) {
                if (!entry.startsWith("hsa")) {
                    nonHuman = true;
                    break;
                }
            }
            if (!nonHuman)
                continue;
            // Convert into a String line
            builder.append(clusterId);
            for (String entry : entries) {
                index = entry.indexOf("(");
                builder.append("\t").append(entry.substring(0, index));
            }
            builder.append("\n");
        }
        fu.close();
        String outFileName = ORTHOMCL_SRC_DIR + "humanOrthoMCL_10232007.txt";
        fu.setOutput(outFileName);
        fu.printLine(builder.toString());
        fu.close();
    }
    
    public void convertEnsemlToUniProt() throws IOException {
        String[] fileNames = new String[] {
                //"flyNonY2HInteractionHsa.txt", 
                //"flyY2HInteractionHsa.txt",
                //"wormY2HInteractionHsa.txt",
                //"wormNonY2HInteractionHsa.txt",          
                "yeastNonY2HInteractionHsa.txt",
                "yeastY2HInteractionHsa.txt"
        };
        Map<String, List<String>> ensembl2uni = loadRedundantMap("results/ensembl2uni.txt");
        int index = 0;
        FileUtility fu = new FileUtility();
        Set<String> pairs = new HashSet<String>();
        String line;
        String partner1;
        String partner2;
        List<String> list1;
        List<String> list2;
        for (String fileName : fileNames) {
            pairs.clear();
            fu.setInput(RESULT_DIR + fileName);
            while ((line = fu.readLine()) != null) {
                index = line.indexOf(" ");
                partner1 = line.substring(0, index);
                list1 = ensembl2uni.get(partner1);
                if (list1 == null || list1.size() == 0)
                    continue;
                partner2 = line.substring(index + 1);
                list2 = ensembl2uni.get(partner2);
                if (list2 == null || list2.size() == 0)
                    continue;
                for (String s1 : list1) {
                    for (String s2 : list2) {
                        pairs.add(createPair(s1, s2));
                    }
                }
            }
            fu.close();
            if (pairs.size() == 0)
                continue;
            FileUtility outFu = new FileUtility();
            index = fileName.lastIndexOf(".");
            outFu.setOutput(RESULT_DIR + fileName.substring(0, index) + "Uni.txt");
            for (String pair : pairs)
                outFu.printLine(pair);
            outFu.close();
            System.out.println(fileName + " in uni: " + pairs.size());
        }
    }
    
    private String createPair(String partner1, String partner2) {
        int compare = partner1.compareTo(partner2);
        if (compare < 0)
            return partner1 + " " + partner2;
        else
            return partner2 + " " + partner1;
    }
    
    /**
     * OrthoMCL uses the old 6 digits, which should be updated to 8 digits. This method is used
     * to update ids by padding two 0 at the left end.
     * @throws IOException
     */
    public void updateOrthoPairForYeast() throws IOException {
        //S0005372  ENSG00000105968
        Set<String> maps = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(RESULT_DIR + "sce2hsaOrthoMCL.txt");
        FileUtility outFu = new FileUtility();
        outFu.setOutput(RESULT_DIR + "sce2hsaOrthoMCLUpdated.txt");
        String line;
        while ((line = fu.readLine()) != null) {
            outFu.printLine(line.substring(0, 1) + "00" + line.substring(1));
        }
        fu.close();
        outFu.close();
    }
    /*
     * These method will generate files that can map the ensembl id to
     * uniprot accessions directly using the files extrated by the ensembl
     * biomart viewer:
     * http://www.ensembl.org/biomart/martview/
     * a typical line in the mart file:
     * ENSMUSG00000053211		ENSMUSP00000069364		P20662	Q8C638
     * 
     * Both the Uniprot accession number and unified Uniprot accession number will be used.
     */
    private Map<String,String> generateUni2EnsemblMap(String mapFile) throws Exception {
    	FileUtility fu = new FileUtility();
    	fu.setInput(mapFile);
    	Map<String,String> uni2EnsMap = new HashMap<String,String>();
    	String line;
    	while((line = fu.readLine()) != null){
    		StringTokenizer st = new StringTokenizer(line);
    		int stSize = st.countTokens();
    		String tokenSG = new String();
    		String tokenSP = new String();
    		String tokenUA = new String();
    		String tokenUUA = new String();
    		//skip the first string
    		if(stSize == 4) {
	    		tokenSG = st.nextToken();
	    		tokenSP = st.nextToken();
	    		tokenUA = st.nextToken();
	    		tokenUUA = st.nextToken();
    		}
    		else if(stSize == 3) {
    			tokenSG = st.nextToken();
	    		tokenSP = st.nextToken();
	    		tokenUA = st.nextToken();
    		}
    			uni2EnsMap.put(tokenUA, tokenSP);
        		uni2EnsMap.put(tokenUUA, tokenSP);
    	}
    	return uni2EnsMap;
    }
    private Map<String,Set<String>> generateEnsembl2UniMap(String mapFile) throws Exception {
    	FileUtility fu = new FileUtility();
    	fu.setInput(mapFile);
    	Map<String,Set<String>> ens2UniMap = new HashMap<String,Set<String>>();
    	String line;
    	while((line = fu.readLine()) != null){
    		StringTokenizer st = new StringTokenizer(line);
    		int stSize = st.countTokens();
    		String tokenSG = new String();
    		String tokenSP = new String();
    		String tokenUA = new String();
    		String tokenUUA = new String();
    		//skip the first string
    		if(stSize == 4) {
	    		tokenSG = st.nextToken();
	    		tokenSP = st.nextToken();
	    		tokenUA = st.nextToken();
	    		tokenUUA = st.nextToken();
    		}
    		else if(stSize == 3) {
    			tokenSG = st.nextToken();
	    		tokenSP = st.nextToken();
	    		tokenUUA = st.nextToken();
    		}
    		Set<String> tmp;
    		if((tmp = ens2UniMap.get(tokenSP)) == null) {
    			tmp = new HashSet<String>();
    			ens2UniMap.put(tokenSP, tmp);
    		}
    		tmp.add(tokenUA);
    		tmp.add(tokenUUA);
    	}
    	return ens2UniMap;
    }
    
    @Test
    public void mapUniInteractions() throws Exception {
    	
    	FileUtility fu = new FileUtility();
    	//364033 MMU 
    	Set<String> orthoInteractions = new HashSet<String>();
    	Set<String> genewaysInteractions = fu.loadInteractions("data/swissprot_361401.txt");
    	Map<String,String> uni2EnsMap = generateUni2EnsemblMap("data/embl/uni2dme.txt");
    	Map<String,Set<String>> ENSMUSP2ENSP = fu.loadSetMap("results/orthomcl/dme2hsaOrthoMCL.txt", "\t",true);
    	Map<String,Set<String>> ENSP2Uni = generateEnsembl2UniMap("data/embl/homo2uni.txt");
    	int debugCounter = 0;

    	for(String interaction:genewaysInteractions) {
//    		StringTokenizer st = new StringTokenizer(interaction);
//    		String term1 = uni2EnsMap.get(st.nextToken());//term1 == ENSMUSP*********;
//    		String term2 = uni2EnsMap.get(st.nextToken());//term2 == ENSMUSP*********;
//    		if(term1 != null && term2 != null) {
//    			Set<String> ENSPTerm1Set = ENSMUSP2ENSP.get(term1);//ENSPTerm1 == ENSP**********
//    			Set<String> ENSPTerm2Set = ENSMUSP2ENSP.get(term2);
//    			//P01581
//    			if(ENSPTerm1Set != null && ENSPTerm2Set != null) {
//    				if(ENSPTerm1Set.size() >= 1 && ENSPTerm2Set.size() >= 1){
//    					//Currently we allow unambigous mapping
//    					Set<String> enspInteractionSet = new HashSet<String>();
//    					tool.buildReactionPairs(enspInteractionSet,ENSPTerm1Set, ENSPTerm2Set);
//    					for(String enspInteraction : enspInteractionSet) {
//    						String tmpTerm1 = tool.getTermFromPairs(enspInteraction, " ", 1);
//    						String tmpTerm2 = tool.getTermFromPairs(enspInteraction, " ", 2);
//    						Set<String> finalUni1Set = ENSP2Uni.get(tmpTerm1);
//        					Set<String> finalUni2Set = ENSP2Uni.get(tmpTerm2);
//        					if(finalUni1Set != null && finalUni2Set != null) {
//        						if(finalUni1Set.size() >= 1 && finalUni2Set.size() >=1 ) {
//        							tool.buildReactionPairs(orthoInteractions, finalUni1Set, finalUni2Set);
//        						}
//        					}
//    					}
//    				}
//    			}
//    		}
    	}
    	fu.saveInteractions(orthoInteractions, "results/orthomcl/dmeOrthoInteractions.txt");
    }
    	
}
