/*
 * Created on Feb 27, 2012
 *
 */
package org.reactome.data;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

/**
 * This class is used to analyze the MITTab file format downloaded from iRefIndex:
 * http://irefindex.uio.no/wiki/?title=README_MITAB2.6_for_iRefIndex&redirect=no. iRefIndex
 * should be used the data source for all protein-protein interaction related work. It is 
 * hoped that this data source can be updated regularly. Currently seemingly it is updated
 * pretty regularly. But we may need to pay attention to it!
 * Note: All PPIs reported in the iIndexRef files are positive interactions.
 * @TODO: Check if iRefIndex database can be updated regularly. 
 * @author gwu
 *
 */
public class IRefIndexMITTabAnalyzer {
    private FileUtility fu = new FileUtility();
    
    public IRefIndexMITTabAnalyzer() {
    }
    
    /**
     * PPIs extracted from mouse (taxon id: 10090).
     * @throws IOException
     */
    @Test
    public void loadMousePPIs() throws IOException {
        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_FILE"),
                    FIConfiguration.getConfiguration().get("IREFINDEX_MOUSE_PPI_FILE"),
                    false);
    }
    
    /**
     * PPIs from yeast (S288c, not baker's yeast) are extracted into a file.
     * @throws IOException
     */
    @Test
    public void loadYeastPPIs() throws IOException {
        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_FILE"),
                    FIConfiguration.getConfiguration().get("IREFINDEX_YEAST_PPI_FILE"),
                    false);
        // Only 63 PPIs can be extracted from the following file
//        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_DIR + "4932.mitab.10182011.txt",
//                    FIConfiguration.getConfiguration().get("IREFINDEX_DIR + "4932PPIsInUniProt.txt",
//                    false);
    }
    
    @Test
    public void loadFlyPPIs() throws IOException {
        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_FLY_FILE"),
                    FIConfiguration.getConfiguration().get("IREFINDEX_FLY_PPI_FILE"),
                    false);
    }
    
    @Test
    public void loadWormPPIs() throws IOException {
        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_WORM_FILE"),
                    FIConfiguration.getConfiguration().get("IREFINDEX_WORM_PPI_FILE"),
                    false);
    }
    
    @Test
    public void loadHumanPPIs() throws IOException {
        extractPPIs(FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_FILE"), 
                    FIConfiguration.getConfiguration().get("IREFINDEX_HUMAN_PPI_FILE"),
                    true);
    }
    
    private void extractPPIs(String srcFileName,
                             String outFileName,
                             boolean needToFilterTohuman) throws IOException {
        long time1 = System.currentTimeMillis();
        fu.setInput(srcFileName);
        String line = fu.readLine(); // MITTAB header
        int compare = 0;
        Set<String> ppis = new HashSet<String>();
        while ((line = fu.readLine()) != null) {
            String[] tokens = line.split("\t");
            String protein1 = tokens[0];
            String uniprotId1 = getUniProtId(protein1);
            String protein2 = tokens[1];
            String uniprotId2 = getUniProtId(protein2);
            if (uniprotId1 == null || uniprotId2 == null)
                continue; // Want to handle UniProt ids only
            compare = uniprotId1.compareTo(uniprotId2);
            if (compare < 0)
                ppis.add(uniprotId1 + "\t" + uniprotId2);
            else if (compare > 0)
                ppis.add(uniprotId2 + "\t" + uniprotId1); // Ignore self-interaction
        }
        fu.close();
        System.out.println("Total PPIs: " + ppis.size());
        if (needToFilterTohuman) {
            Set<String> totalIds = InteractionUtilities.grepIDsFromInteractions(ppis);
            checkWithUniProtIds(totalIds);
            UniProtAnalyzer uniProtHelper = new UniProtAnalyzer();
            Map<String, String> uniProtIdMap = uniProtHelper.loadUniProtIDsMap();
            Set<String> uniProtIds = uniProtIdMap.keySet();
            for (Iterator<String> it = ppis.iterator(); it.hasNext();) {
                String ppi = it.next();
                String[] tokens = ppi.split("\t");
                if (!uniProtIds.contains(tokens[0]) ||
                    !uniProtIds.contains(tokens[1]))
                    it.remove();
            }
            totalIds = InteractionUtilities.grepIDsFromInteractions(ppis);
            checkWithUniProtIds(totalIds);
        }
        fu.saveInteractions(ppis,
                            outFileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time: " + (time2 - time1));
    }
    
    private void checkWithUniProtIds(Set<String> totalIds) throws IOException {
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, String> swissProtMap = uniProtAnalyzer.loadSwissProtIDsMap();
        Map<String, String> uniProtMap = uniProtAnalyzer.loadUniProtIDsMap();
        System.out.println("Total ids: " + totalIds.size());
        Set<String> shared = InteractionUtilities.getShared(totalIds, swissProtMap.keySet());
        System.out.println(" in SwissProt: " + shared.size());
        shared = InteractionUtilities.getShared(totalIds, uniProtMap.keySet());
        System.out.println(" in UniProt: " + shared.size());
        totalIds.removeAll(shared);
//        int count = 0;
//        for (String id : totalIds) {
//            System.out.println(id);
//            count ++;
//            if (count == 100)
//                break;
//        }
    }
    
    private String getUniProtId(String protein) {
        if (!protein.startsWith("uniprotkb"))
            return null;
        int index = protein.indexOf(":");
        return protein.substring(index + 1);
    }
    
}
