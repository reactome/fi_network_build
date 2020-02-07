/*
 * Created on May 14, 2008
 *
 */
package org.reactome.fi;

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
import org.hibernate.Query;
import org.hibernate.SessionFactory;
import org.hibernate.classic.Session;
import org.junit.Test;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.ReactomeAnalyzerTopicHelper;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.InteractionUtilities;
import org.reactome.r3.util.R3Constants;

/**
 * This class is used to analyze the growth of FIs between a period of time 
 * @author wgm
 */
public class ReactomeFIGrowthAnalyzer {
    // cache these values to increate the performance
    private UniProtAnalyzer uniProtAnalyzer;
    private Map<String, String> acIdMap;
    private Set<String> uniSet;
    // Make methods call easy
    private Set<String> allFIs;
    private Set<String> pathwayFIs;
    private Set<String> predicatedFIs;
    
    public ReactomeFIGrowthAnalyzer() {
    }
    
    /**
     * Call this method to load FIs from pre-generated FI files from the FI
     * network.
     * @throws IOException
     */
    private void loadFIs() throws IOException {
        // Load fis
        String fileName = R3Constants.RESULT_DIR + "FI73_042108.txt";
        FileUtility fu = new FileUtility();
        allFIs = fu.loadInteractions(fileName);
        System.out.println("Total network FIs: " + allFIs.size());
        fileName = R3Constants.RESULT_DIR + "FI73_Predicated_042108.txt";
        predicatedFIs = fu.loadInteractions(fileName);
        System.out.println("Network FIs from predicated: " + predicatedFIs.size());
        fileName = R3Constants.RESULT_DIR + "FI73_Pathway_042108.txt";
        pathwayFIs = fu.loadInteractions(fileName);
        System.out.println("Network FIs from pathways: " + pathwayFIs.size());
    }
        
    /**
     * This method is used to calculate how many FIs in the FI network
     * for a pathway have been curated in the new Reactome database.
     * @throws Exception
     */
    @Test
    public void checkIndividualPathwayFIIncrease() throws Exception {
        Long[] dbIds = getIndividualTestPathways();
        // Load fis
        loadFIs();
        MySQLAdaptor oldDba = getOldDba();
        MySQLAdaptor newDba = getNewDba();
        ReactomeAnalyzer oldAnalyzer = new ReactomeAnalyzer();
        oldAnalyzer.setMySQLAdaptor(oldDba);
        ReactomeAnalyzer newAnalyzer = new ReactomeAnalyzer();
        newAnalyzer.setMySQLAdaptor(newDba);
        for (Long dbId : dbIds) {
            GKInstance oldPathway = oldDba.fetchInstance(dbId);
            GKInstance newPathway = newDba.fetchInstance(dbId);
            Set<String> oldIds = oldAnalyzer.grepIDsFromTopic(oldPathway);
            oldIds = filterIds(oldIds);
            Set<String> oldFIs = oldAnalyzer.grepInteractionsForTopic(oldPathway);
            oldFIs = filterFIs(oldFIs);
            Set<String> newFIs = newAnalyzer.grepInteractionsForTopic(newPathway);
            newFIs = filterFIs(newFIs);
            newFIs.removeAll(oldFIs);
            System.out.println(oldPathway.getDisplayName());
            System.out.println("Total new FIs: " + newFIs.size());
            // Want to check how many new FIs can be added to the old pahtway
            Set<String> addableFIsFromAll = new HashSet<String>();
            Set<String> addableFIsFromPathway = new HashSet<String>();
            Set<String> addableFIsFromPredicated = new HashSet<String>();
            checkAddableFIs(allFIs, oldIds, oldFIs, addableFIsFromAll);
            checkAddableFIs(pathwayFIs, oldIds, oldFIs, addableFIsFromPathway);
            checkAddableFIs(predicatedFIs, oldIds, oldFIs, addableFIsFromPredicated);
            // We want to compare the addable FIs and newly added FIs
            // First print out information
            System.out.println("\tAddable FIs from all:");
            printOutFIInfo(newFIs, addableFIsFromAll);
            System.out.println("\tAddable FIs from pathway:");
            printOutFIInfo(newFIs, addableFIsFromPathway);
            System.out.println("\tAddable FIs from predicated:");
            printOutFIInfo(newFIs, addableFIsFromPredicated);
            System.out.println();
        }
    }

    private void printOutFIInfo(Set<String> newFIs,
                                Set<String> addableFIs) {
        int size = addableFIs.size();
        System.out.println("\t\t" + size);
        addableFIs.retainAll(newFIs);
        int size1 = addableFIs.size();
        System.out.println("\t\tIn new pathway: " + size1 + 
                           " (" + size1 / (double)size + ")");
    }
    
    private void checkAddableFIs(Set<String> fis,
                                 Set<String> ids,
                                 Set<String> existingFis,
                                 Set<String> addableFIs) {
        for (String fi : fis) {
            if (existingFis.contains(fi))
                continue; // Don't add FIs in the old pathway
            // Check this fi can be addable
            for (String id : ids) {
                if (fi.contains(id)) {
                    addableFIs.add(fi);
                    break;
                }
            }
        }
    }
    
    private Long[] getIndividualTestPathways() {
        // The following pathways will be checked
//      [Pathway:109582] Hemostasis
//      [Pathway:109581] Apoptosis
//      [Pathway:170877] Oncogenic and Tumor suppressor pathways\
//      [Pathway:168256] Immune System signaling\
//      [Pathway:74160] Gene Expression
//      [Pathway:166520] Nerve Growth Factor (NGF) signaling
        Long[] dbIds = new Long[] {
                109582L,
                109581L,
                170877L,
                168256L,
                74160L,
                166520L
        };
        return dbIds;
    }
    
    /**
     * This method is used to check FI coverage by filtering some big complexes first
     * Here ribosome in gene expression is tested to see if the coverage will be alternated.
     * @throws Exception
     */
    @Test
    public void checkIndividualPathwayFICoverageWithFilter() throws Exception {
        loadFIs();
        Long pathwayId = 74160L;
        Long ribosomeId = 72500L;
        MySQLAdaptor newDba = getNewDba();
        ReactomeAnalyzer newAnalyzer = new ReactomeAnalyzer();
        newAnalyzer.setMySQLAdaptor(newDba);
        MySQLAdaptor oldDba = getOldDba();
        ReactomeAnalyzer oldAnalyzer = new ReactomeAnalyzer();
        Set<String> oldFullFIs = getOldFISet();
        oldFullFIs = filterFIs(oldFullFIs);
        oldAnalyzer.setMySQLAdaptor(oldDba);
        // Try to get all ribosome protein ids
        GKInstance oldRibosome = oldDba.fetchInstance(ribosomeId);
        GKInstance newRibosome = newDba.fetchInstance(ribosomeId);
        Set<String> oldProteinsInRibosome = grepIdsFromComplex(oldRibosome);
        Set<String> newProteinsInRibosome = grepIdsFromComplex(newRibosome);
        System.out.println("Ids in old ribosome: " + oldProteinsInRibosome.size());
        System.out.println("Ids in new ribosome: " + newProteinsInRibosome.size());
        Set<String> oldFIs = fetchFIsForPathway(pathwayId, oldAnalyzer);
        Set<String> newFIs = fetchFIsForPathway(pathwayId, newAnalyzer);
        Set<String> newCopy = new HashSet<String>(newFIs);
        Set<String> oldCopy = new HashSet<String>(oldFIs);
        // Filter all FIs coming from ribosome
        filterFIsForComplex(newCopy, newProteinsInRibosome);
        filterFIsForComplex(oldCopy, oldProteinsInRibosome);
        
        Set<String> oldCopy1 = new HashSet<String>(oldCopy);
        oldCopy.removeAll(newCopy);
        newCopy.removeAll(oldCopy1);
        System.out.println("\tNewly added FIs: " + newCopy.size());
        System.out.println("\tDeleted FIs: " + oldCopy.size());
        // check coverage now
        Set<String> newlyFIs = new HashSet<String>(newCopy);
        newlyFIs.retainAll(allFIs);
        System.out.println("\tNew FIs covered by all FIs: " + newlyFIs.size() + 
                           " (" + newlyFIs.size() / (double) newCopy.size() + ")");
        newlyFIs = new HashSet<String>(newCopy);
        newlyFIs.retainAll(pathwayFIs);
        System.out.println("\tNew FIs covered by pathway FIs: " + newlyFIs.size() + 
                           " (" + newlyFIs.size() / (double) newCopy.size() + ")");
        newlyFIs = new HashSet<String>(newCopy);
        newlyFIs.retainAll(predicatedFIs);
        System.out.println("\tNew FIs covered by predicated FIs: " + newlyFIs.size() + 
                           " (" + newlyFIs.size() / (double)newCopy.size() + ")");
    }
    
    /**
     * This helper method is used to filter FIs that brought by complexes.
     * @param fis
     * @param complexIds
     */
    private void filterFIsForComplex(Set<String> fis,
                                     Set<String> complexIds) {
        System.out.println("Before complex filtering: " + fis.size());
        int index = 0;
        for (Iterator<String> it = fis.iterator(); it.hasNext();) {
            String fi = it.next();
            index = fi.indexOf(" ");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            if (complexIds.contains(id1) ||
                complexIds.contains(id2))
                it.remove();
        }
        System.out.println("After complex filtering: " + fis.size());
    }
    
    /**
     * This method is used to grep protein ids from a complex.
     * @param complex
     * @return
     * @throws Exception
     */
    private Set<String> grepIdsFromComplex(GKInstance complex) throws Exception {
        Set<GKInstance> refPepSeqs = new ReactomeAnalyzerTopicHelper().grepRefPepSeqs(complex);
        Set<String> ids = new HashSet<String>();
        for (GKInstance refPepSeq : refPepSeqs) {
            String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            if (id != null)
                ids.add(id);
        }
        return ids;
    }
    
    private Set<String> fetchFIsForPathway(Long pathwayId,
                                           ReactomeAnalyzer analyzer) throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) analyzer.getMySQLAdaptor();
        GKInstance pathway = dba.fetchInstance(pathwayId);
        Set<String> fis = analyzer.grepInteractionsForTopic(pathway);
        return filterFIs(fis);
    }
    
    /**
     * This method is used to check if newly added FIs to some pathways
     * have been predicated by the FI network.
     * @throws Exception
     */
    @Test
    public void checkIndividualPathwayFICoverage() throws Exception {
        Long[] dbIds = getIndividualTestPathways();
        // Load fis
        loadFIs();
        MySQLAdaptor newDba = getNewDba();
        ReactomeAnalyzer newAnalyzer = new ReactomeAnalyzer();
        newAnalyzer.setMySQLAdaptor(newDba);
        MySQLAdaptor oldDba = getOldDba();
        ReactomeAnalyzer oldAnalyzer = new ReactomeAnalyzer();
        Set<String> oldFullFIs = getOldFISet();
        oldFullFIs = filterFIs(oldFullFIs);
        oldAnalyzer.setMySQLAdaptor(oldDba);
        for (Long dbId : dbIds) {
            GKInstance oldInstance = oldDba.fetchInstance(dbId);
            GKInstance newInstance = newDba.fetchInstance(dbId);
            Set<String> oldFIs = oldAnalyzer.grepInteractionsForTopic(oldInstance);
            oldFIs = filterFIs(oldFIs);
            Set<String> newFIs = newAnalyzer.grepInteractionsForTopic(newInstance);
            newFIs = filterFIs(newFIs);
            // Get newly added ids
            Set<String> oldCopy = new HashSet<String>(oldFIs);
            Set<String> newCopy = new HashSet<String>(newFIs);
            oldCopy.removeAll(newFIs); // Deleted FIs in the new database
            newCopy.removeAll(oldFIs); // Added FIs in the new database
            //newCopy.removeAll(oldFullFIs); // Remove any FIs annotated in the old Reactome database
            System.out.println(oldInstance.getDisplayName());
            System.out.println("\tNewly added FIs: " + newCopy.size());
            System.out.println("\tDeleted FIs: " + oldCopy.size());
            // check coverage now
            Set<String> newlyFIs = new HashSet<String>(newCopy);
            newlyFIs.retainAll(allFIs);
            System.out.println("\tNew FIs covered by all FIs: " + newlyFIs.size() + 
                               " (" + newlyFIs.size() / (double) newCopy.size() + ")");
            newlyFIs = new HashSet<String>(newCopy);
            newlyFIs.retainAll(pathwayFIs);
            System.out.println("\tNew FIs covered by pathway FIs: " + newlyFIs.size() + 
                               " (" + newlyFIs.size() / (double) newCopy.size() + ")");
            newlyFIs = new HashSet<String>(newCopy);
            newlyFIs.retainAll(predicatedFIs);
            System.out.println("\tNew FIs covered by predicated FIs: " + newlyFIs.size() + 
                               " (" + newlyFIs.size() / (double)newCopy.size() + ")");
            System.out.println();
        }
    }
    
    private MySQLAdaptor getNewDba() throws Exception {
        MySQLAdaptor newDba = new MySQLAdaptor("localhost",
                                               "gk_central_051208",
                                               "root",
                                               "macmysql01",
                                               3306);
        return newDba;
    }
    
    private MySQLAdaptor getOldDba() throws Exception {
        MySQLAdaptor oldDba = new MySQLAdaptor("localhost",
                                               "gk_central_101606",
                                               "root",
                                               "macmysql01",
                                               3306);
        return oldDba;
    }
    
    /**
     * This method is used to check pathway growth based on FIs.
     * @throws Exception
     */
    @Test
    public void checkPathwayGrowthBasedOnFIs() throws Exception {
        MySQLAdaptor newDba = getNewDba();
        Map<Long, Set<String>> newPathwayToFIs = getPathwayToFIs(newDba);
        MySQLAdaptor oldDba = getOldDba();
        Map<Long, Set<String>> oldPathwayToFIs = getPathwayToFIs(oldDba);
        // Compare
        for (Iterator<Long> it = newPathwayToFIs.keySet().iterator(); it.hasNext();) {
            Long dbId = it.next();
            Set<String> newIds = newPathwayToFIs.get(dbId);
            Set<String> oldIds = oldPathwayToFIs.get(dbId);
            if (oldIds == null)
                continue;
            GKInstance pathway = oldDba.fetchInstance(dbId);
            System.out.println(pathway + "\t" + oldIds.size() + "\t" + newIds.size() + "\t" + (newIds.size() - oldIds.size()));
        }
    }
    
    /**
     * This method is used to compare the size of pathways listed in the FrontPageItem
     * instances between two snapshots of gk_central databases.
     * @throws Exception
     */
    @Test
    public void checkPathwayGrowthBasedOnIds() throws Exception {
        MySQLAdaptor newDba = getNewDba();
        Map<Long, Set<String>> newModuleToProteinIds = getPathwayToProteinIds(newDba);
        MySQLAdaptor oldDba = getOldDba();
        Map<Long, Set<String>> oldModuleToProteinsIds = getPathwayToProteinIds(oldDba);
        // Compare
        for (Iterator<Long> it = newModuleToProteinIds.keySet().iterator(); it.hasNext();) {
            Long dbId = it.next();
            Set<String> newIds = newModuleToProteinIds.get(dbId);
            Set<String> oldIds = oldModuleToProteinsIds.get(dbId);
            if (oldIds == null)
                continue;
            GKInstance pathway = oldDba.fetchInstance(dbId);
            System.out.println(pathway + "\t" + oldIds.size() + "\t" + newIds.size() + "\t" + (newIds.size() - oldIds.size()));
        }
    }
    
    /**
     * This helper method is used to get a map from a pathway to a set of FIs.
     * @param dba
     * @return
     * @throws Exception
     */
    private Map<Long, Set<String>> getPathwayToFIs(MySQLAdaptor dba) throws Exception {
        List<GKInstance> pathways = getPathwayList(dba);
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        analyzer.setMySQLAdaptor(dba);
        Map<Long, Set<String>> pathwayToFIs = new HashMap<Long, Set<String>>();
        for (GKInstance pathway : pathways) {
            Set<String> fis = analyzer.grepInteractionsForTopic(pathway);
            pathwayToFIs.put(pathway.getDBID(), fis);
        }
        return pathwayToFIs;
    }
    
    /**
     * This helper method is used to generate a map from a pathway to protein ids.
     * @param dba
     * @return
     * @throws Exception
     */
    private Map<Long, Set<String>> getPathwayToProteinIds(MySQLAdaptor dba) throws Exception {
        List<GKInstance> pathways = getPathwayList(dba);
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        analyzer.setMySQLAdaptor(dba);
        // Need to use DB_ID since _displayName may have been changed esp. in Signaling
        // pathways.
        Map<Long, Set<String>> pathwayToIds = new HashMap<Long, Set<String>>();
        for (GKInstance pathway : pathways) {
            Set<String> ids = analyzer.grepIDsFromTopic(pathway);
            pathwayToIds.put(pathway.getDBID(), 
                             ids);
        }
        return pathwayToIds;
    }
    
    /**
     * This helper method is used to generate a pathway list for comparision.
     * @param dba
     * @return
     * @throws Exception
     */
    private List<GKInstance> getPathwayList(MySQLAdaptor dba) throws Exception {
        // Need to get the frontpageItem
        Collection collection = dba.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        GKInstance frontPageItem = (GKInstance) collection.iterator().next();
        List modules = frontPageItem.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        // Want to break up Signaling Pathway module into individual pathways
        List<GKInstance> pathways = new ArrayList<GKInstance>();
        for (Iterator it = modules.iterator(); it.hasNext();) {
            GKInstance module = (GKInstance) it.next();
            if (module.getDisplayName().equals("Signaling Pathways")) {
                Collection list = null;
                if (module.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    list = module.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
                else if (module.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent))
                    list = module.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                pathways.addAll(list);
            }
            else
                pathways.add(module);
        }
        return pathways;
    }
    
    /**
     * This method is used to check how many newly added FIs
     * have fallen into our FI networks.
     * @throws Exception
     */
    @Test
    public void checkFIsCoverage() throws Exception {
        Set<String> newlyAddedFIs = getNewlyAddedFIs();
        System.out.println("Newly added FIs: " + newlyAddedFIs.size());
        // Load fis
        loadFIs();
        // Check how many in allFIs
        allFIs.retainAll(newlyAddedFIs);
        System.out.println("Covered in all fis: " + allFIs.size());
        predicatedFIs.retainAll(newlyAddedFIs);
        System.out.println("Covered in predicated: " + predicatedFIs.size());
        pathwayFIs.retainAll(newlyAddedFIs);
        System.out.println("Covered in pathways: " + pathwayFIs.size());
    }
    
    /**
     * This method is used to pull out newly added human FIs from a new Reactome
     * database.
     * @return
     * @throws Exception
     */
    private Set<String> getNewlyAddedFIs() throws Exception {
        Set<String> fiSet1 = getPathwayFIsFromHiberante();
        Set<String> fiSet2 = getNewFISet();
        fiSet2 = filterFIs(fiSet2);
        fiSet2.removeAll(fiSet1);
        return fiSet2;
    }
    
    /**
     * This method is used to check the growth of the FIs between two timepoints
     * of gk_central.
     * @throws Exception
     */
    @Test
    public void checkGrowth() throws Exception {
        //Set<String> fiSet1 = getOldFISet();
        Set<String> fiSet1 = getPathwayFIsFromHiberante();
        fiSet1 = filterFIs(fiSet1);
        Set<String> proteinIds1 = InteractionUtilities.grepIDsFromInteractions(fiSet1);
        System.out.println("Reactome on 10/16/06: " + fiSet1.size());
        System.out.println("\tProteins: " + proteinIds1.size());
        System.out.println();
        
        Set<String> fiSet2 = getNewFISet();
        fiSet2 = filterFIs(fiSet2);
        Set<String> proteinIds2 = InteractionUtilities.grepIDsFromInteractions(fiSet2);
        System.out.println("Reactome on 05/12/08: " + fiSet2.size());
        System.out.println("\tProteins: " + proteinIds2.size());
        
        // Check how many new FIs at the newly released FIs
        Set<String> copy1 = new HashSet<String>(fiSet1);
        Set<String> copy2 = new HashSet<String>(fiSet2);
        copy1.removeAll(fiSet2);
        copy2.removeAll(fiSet1);
        System.out.println("\nNewly added FIs: " + copy2.size());
        System.out.println("Deleted FIs: " + copy1.size());
        // Want to a quick check regarding deleted FIs
        int c = 0;
        for (String fi : copy1) {
            System.out.println(fi);
            c ++;
            if (c > 10)
                break;
        }
    }
    
    private Set<String> filterFIs(Set<String> fis) throws Exception {
        if (uniProtAnalyzer == null) {
            uniProtAnalyzer = new UniProtAnalyzer();
            acIdMap = uniProtAnalyzer.loadUniProtIDsMap();
            uniSet = acIdMap.keySet();
        }
        FunctionalInteractionAnalyzer fiAnalyzer = new FunctionalInteractionAnalyzer();
        fiAnalyzer.filterNonHumanIds(fis, uniSet, uniProtAnalyzer);
        return fiAnalyzer.filterRedundencyInteractions(fis, 
                                                       acIdMap);
    }
    
    private Set<String> filterIds(Set<String> ids) throws Exception {
        if (uniProtAnalyzer == null) {
            uniProtAnalyzer = new UniProtAnalyzer();
            acIdMap = uniProtAnalyzer.loadUniProtIDsMap();
            uniSet = acIdMap.keySet();
        }
        System.out.println("Before id filtering: " + ids.size());
        Set<String> rtn = new HashSet<String>();
        for (String id : ids) {
            String mapped = acIdMap.get(id);
            if (mapped != null)
                rtn.add(mapped);
        }
        System.out.println("After id filtering: " + rtn.size());
        return rtn;
    }
    
    private Set<String> getNewFISet() throws Exception {
        ReactomeAnalyzer analyzer2 = new ReactomeAnalyzer();
        MySQLAdaptor dba2 = getNewDba();
        analyzer2.setMySQLAdaptor(dba2);
        Set<String> fiSet2 = analyzer2.extractInteractionSet();
        return fiSet2;
    }
    
    private Set<String> getOldFISet() throws Exception {
        ReactomeAnalyzer analyzer1 = new ReactomeAnalyzer();
        MySQLAdaptor dba1 = getOldDba();
        analyzer1.setMySQLAdaptor(dba1);
        Set<String> fiSet1 = analyzer1.extractInteractionSet();
        return fiSet1;
    }
    
    /**
     * This helper method is used to get a list of FIs from pathways via Hiberante API.
     * @return
     * @throws Exception
     */
    private Set<String> getPathwayFIsFromHiberante() throws Exception {
        HibernateFIReader hiFiAnalyzer = new HibernateFIReader();
        hiFiAnalyzer.initSession();
        SessionFactory sf = hiFiAnalyzer.getSessionFactory();
        Session session = sf.openSession();
        Query query = session.createQuery("FROM Interaction as i WHERE i.evidence is null");
        List list = query.list();
        System.out.println("Total interactions from pathways: " + list.size());
        boolean isFromReactome = false;
        Set<String> reactomeFIs = new HashSet<String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Interaction interaction = (Interaction) it.next();
            Set<ReactomeSource> sources = interaction.getReactomeSources();
            isFromReactome = false;
            for (ReactomeSource source : sources) {
                if (source.getDataSource().equals("Reactome")) {
                    isFromReactome = true;
                    break;
                }
            }
            if (!isFromReactome)
                continue;
            Protein protein1 = interaction.getFirstProtein();
            Protein protein2 = interaction.getSecondProtein();
            String fi = protein1.getPrimaryAccession() + " " +
                        protein2.getPrimaryAccession();
            reactomeFIs.add(fi);
        }
        return reactomeFIs;
    }
}
