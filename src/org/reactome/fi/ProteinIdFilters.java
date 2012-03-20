/*
 * Created on May 27, 2008
 *
 */
package org.reactome.fi;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.reactome.fi.ProteinSequenceHandler.Sequence;


/**
 * This class groups a list of human protein ids. Basically human protein ids
 * will be filtered based if an id is a human protein id, or based on redundancies.
 * Also a list of ids can be consolidated based on sequences. The filters used here
 * should be the same used to generate the FI database so that the user of this class
 * doesn't need to call the database to do filtering.
 * @author wgm
 *
 */
public class ProteinIdFilters {
    // Cached to inrease performance
    private UniProtAnalyzer uniProtAnalyzer;
    private Map<String, String> acIdMap;
    private Set<String> uniIdSet;
    // For sequence based consolidation
    private ProteinSequenceHandler seqHandler;
    private Map<String, Sequence> idToSequence;
    private Map<String, String> checksumToId;
    
    public ProteinIdFilters() {
    }
    
    /**
     * Call this method to do several filters together.
     * @param ids
     * @return
     * @throws Exception
     */
    public Set<String> filter(Collection<String> ids) throws Exception {
        Set<String> filtered = filterNonHumanIds(ids);
        filtered = removeRedundancy(filtered);
        filtered = consolidateBasedOnSequence(filtered);
        return filtered;
    }
    
    /**
     * Clean up a set of protein pairs to human only
     * @param idMap
     * @param ppis
     * @return a new set of PPIs after filter. Note: a new version of Set<String> has to be
     * used by the client.
     */
    public Set<String> cleanUpVsUniProt(Set<String> ppis) throws Exception {
        iniUniProtAnalyzer();
        Set<String> rtn = new HashSet<String>();
        int index = 0;
        int compare;
        for (String ppi : ppis) {
            index = ppi.indexOf("\t");
            String id1 = ppi.substring(0, index);
            String id2 = ppi.substring(index + 1);
            String mapped1 = acIdMap.get(id1);
            String mapped2 = acIdMap.get(id2);
            if (mapped1 == null || mapped2 == null) {
//                System.out.println(id1 + " or " + id2 + " cannot be mapped!");
                continue; // Use UniProt ids only
            }
            compare = mapped1.compareTo(mapped2);
            if (compare < 0)
                rtn.add(mapped1 + "\t" + mapped2);
            else if (compare > 0)
                rtn.add(mapped2 + "\t" + mapped1);
        }
        return rtn;
    }
    
    /**
     * This filter is used to remove any non-human protein ids.
     * @param ids
     * @return
     * @throws Exception
     */
    public Set<String> filterNonHumanIds(Collection<String> ids) throws Exception {
        iniUniProtAnalyzer();
        Set<String> filtered = new HashSet<String>();
        for (String id : ids) {
            if (uniProtAnalyzer.isHumanID(uniIdSet, id)) {
                filtered.add(id);
            }
        }
        return filtered;
    }
    
    private void iniUniProtAnalyzer() throws Exception {
        if (uniProtAnalyzer == null) {
            uniProtAnalyzer = new UniProtAnalyzer();
            acIdMap = uniProtAnalyzer.loadUniProtIDsMap();
            uniIdSet = acIdMap.keySet();
        }
    }
    
    /**
     * Use this method to remove redundancy protein ids from Trembl.
     * @param ids
     * @return
     * @throws Exception
     */
    public Set<String> removeRedundancy(Collection<String> ids) throws Exception {
        iniUniProtAnalyzer();
        Set<String> filtered = new HashSet<String>();
        for (String id : ids) {
            String mapped = acIdMap.get(id);
            if (mapped != null)
                filtered.add(mapped);
            else
                filtered.add(id);
        }
        return filtered;
    }
    
    /**
     * A helper method to filter non-human FIs.
     * @param interactions
     * @param uniSet
     * @param uniAnalyzer
     */
    public void filterNonHumanIdsForPairs(Set<String> interactions) throws Exception {
        iniUniProtAnalyzer();
        for(Iterator<String> it = interactions.iterator(); it.hasNext();) {
            String interaction = it.next();
            if (containsNonHumanId(interaction))
                it.remove();
        }
        System.out.println("After filtering: " + interactions.size());
    }
    
    /**
     * A helper method to check if a protein pair has non-human identifiers
     * @param interaction
     * @param uniSet
     * @param uniAnalyzer
     * @return
     */
    private boolean containsNonHumanId(String interaction) {
        int index = interaction.indexOf("\t");
        String id1 = interaction.substring(0, index);
        String id2 = interaction.substring(index + 1);
        if (!uniProtAnalyzer.isHumanID(uniIdSet, id1) ||
            !uniProtAnalyzer.isHumanID(uniIdSet, id2)) {
//            System.out.println("Not human proteins: " + interaction);
            return true;
        }
        return false;
    }
    
    public Set<String> filterRedundencyInteractions(Set<String> interactions) throws Exception {
        iniUniProtAnalyzer();
        Set<String> filtered = new HashSet<String>();
        int index = 0;
        for (String i : interactions) {
            index = i.indexOf("\t");
            String id1 = i.substring(0, index);
            String id2 = i.substring(index + 1);
            String tmpId1 = acIdMap.get(id1);
            if (tmpId1 == null)
                tmpId1 = id1;
            String tmpId2 = acIdMap.get(id2);
            if (tmpId2 == null)
                tmpId2 = id2;
            int compare = tmpId1.compareTo(tmpId2);
            if (compare < 0) {
                filtered.add(tmpId1 + "\t" + tmpId2);
            }
            else if (compare > 0) {
                filtered.add(tmpId2 + "\t" + tmpId1);
            }
        }
        return filtered;
    }
    
    /**
     * This method is used to normalize protein pair sets to remove any non-human ids,
     * redudancy based on id mappings and sequence checksums.
     * @param pairs
     * @return
     * @throws Exception
     */
    public Set<String> normalizeProteinPairs(Set<String> pairs) throws Exception {
        filterNonHumanIdsForPairs(pairs);
        pairs = filterRedundencyInteractions(pairs);
        initSequenecHandler();
        pairs = seqHandler.consolidateInteractionsUseChecksum(pairs);
        return pairs;
    }
    
    /**
     * Use this method to consolidate a list of protein ids based on amino acid sequences.
     * @param ids
     * @throws Exception
     */
    public Set<String> consolidateBasedOnSequence(Collection<String> ids) throws Exception {
        initSequenecHandler();
        Set<String> filtered = new HashSet<String>();
        int index = 0;
        for (String id : ids) {
            // default as UniProt
            if (!id.contains(":"))
                id = "UniProt:" + id;
            Sequence seq = idToSequence.get(id);
            if (seq == null) {
                throw new IllegalStateException(id + ": no sequence available!");
            }
            String mappedId = checksumToId.get(seq.getChecksum());
            if (mappedId == null) {
                throw new IllegalStateException(id + ": no sequence available!");
            }
            // Remove the database part
            index = mappedId.indexOf(":");
            if (index > 0)
                mappedId = mappedId.substring(index + 1);
            filtered.add(mappedId);
        }
        return filtered;
    }
    
    private void initSequenecHandler() throws Exception {
        if (seqHandler != null)
            return;
        seqHandler = new ProteinSequenceHandler();
        idToSequence = seqHandler.getAllSequences();
        checksumToId = seqHandler.generateUniqueChecksumToAccessionMap();
    }
    
}
