/*
 * Created on Jun 30, 2006
 *
 */
package org.reactome.fi.util;

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
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;

/**
 * A utility class for common interaction related parsing.
 * @author guanming
 */
public class InteractionUtilities {
    
    public static <E> String joinStringElements(String delimit,
                                                Collection<E> elements) {
        StringBuilder builder = new StringBuilder();
        for (Iterator<E> it = elements.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append(delimit);
        }
        return builder.toString();
    }
    
    public static boolean isShared(Collection<String> c1,
                                   Collection<String> c2) {
        for (String name1 : c1) {
            if (c2.contains(name1))
                return true;
        }
        return false;
    }
    
    public static Set<String> grepAllGenes(Map<String, Set<String>> sampleToGenes) {
        Set<String> genes = new HashSet<String>();
        for (Set<String> set : sampleToGenes.values())
            genes.addAll(set);
        return genes;
    }
    
    public static Set<String> getShared(Collection<String> c1, Collection<String> c2) {
        Set<String> shared = new HashSet<String>(c1);
        shared.retainAll(c2);
        return shared;
    }
    
    public static String generateFIFromGene(String gene1, String gene2) {
        int compare = gene1.compareTo(gene2);
        if (compare < 0)
            return gene1 + "\t" + gene2;
        else
            return gene2 + "\t" + gene1;
    }
    
    /**
     * This method is used to add an element to a keyed map.
     * @param keyToSet
     * @param key
     * @param element
     */
    public static <K, V> void addElementToSet(Map<K, Set<V>> keyToSet,
                                              K key,
                                              V element) {
        Set<V> set = keyToSet.get(key);
        if (set == null) {
            set = new HashSet<V>();
            keyToSet.put(key, set);
        }
        set.add(element);
    }

    public static Set<String> generateRandomPairs(Set<String> ids,
                                                  int size,
                                                  Set<String> excluded) {
        List<String> idList = new ArrayList<String>(ids);
        Set<String> pairs = new HashSet<String>();
        int total = size;
        int comp = 0;
        int index1, index2;
        String id1, id2;
        int totalSize = idList.size();
        String pair = null;
        while (pairs.size() < total) {
            index1 = (int) (Math.random() * totalSize);
            index2 = (int) (Math.random() * totalSize);
            if (index1 == index2)
                continue;
            id1 = idList.get(index1);
            id2 = idList.get(index2);
            comp = id1.compareTo(id2);
            if (comp < 0)
                pair = id1 + "\t" + id2;
            else if (comp > 0)
                pair = id2 + "\t" + id1;
            if (excluded.contains(pair))
                continue;
            pairs.add(pair);
        }
        return pairs;
    }
    
    public static List<String> convertArrayToList(Object[] objs) {
        List<String> list = new ArrayList<String>();
        for (Object obj : objs)
            list.add(obj.toString());
        return list;
    }
    
    public static Set<String> grepIDsFromInteractions(Set<String> interactions) {
        Set<String> ids = new HashSet<String>();
        int index = 0;
        for (String pair : interactions) {
            // Sometimes tab is used for delimiting: e.g. FIs in names
            if (pair.contains("\t"))
                index = pair.indexOf("\t");
            else
                index = pair.indexOf(" ");
            ids.add(pair.substring(0, index));
            ids.add(pair.substring(index + 1));
        }
        return ids;
    }
    
    public static Set<String> grepFIsContains(String id, Set<String> fis) {
        Set<String> rtn = new HashSet<String>();
        int index = 0;
        String id1, id2;
        for (String fi : fis) {
            //if (fi.contains(id)) // This check is not right. For example, gene METAP1 will be treated as true
            // for MET
            //    rtn.add(fi);
            index = fi.indexOf("\t");
            id1 = fi.substring(0, index);
            id2 = fi.substring(index + 1);
            if (id1.equals(id) ||
                id2.equals(id))
                rtn.add(fi);
        }
        return rtn;
    }
    
    public static Set<String> grepFIsContains(Collection<String> ids, Set<String> fis) {
        Set<String> rtn = new HashSet<String>();
        int index = 0;
        String id1 = null;
        String id2 = null;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            id1 = fi.substring(0, index);
            id2 = fi.substring(index + 1);
            if (ids.contains(id1) ||
                ids.contains(id2))
                rtn.add(fi);
        }
        return rtn;
    }
    
    public static Set<String> getFIs(Collection<String> genes, 
                                     Set<String> fis) {
        Set<String> rtn = new HashSet<String>();
        List<String> list = new ArrayList<String>(genes);
        int compare = 0;
        for (int i = 0; i < list.size() - 1; i++) {
            String gene1 = list.get(i);
            for (int j = i + 1; j < list.size(); j++) {
                String gene2 = list.get(j);
                compare = gene1.compareTo(gene2);
                String key = null;
                if (compare < 0)
                    key = gene1 + "\t" + gene2;
                else
                    key = gene2 + "\t" + gene1;
                if (fis.contains(key))
                    rtn.add(key);
                //else
                //    System.out.println(key);
            }
        }
        return rtn;
    }
    
    public static Set<String> grepIDsFromFuncInt(List<Interaction> interactions) {
        Set<String> ids = new HashSet<String>();
        for (Interaction i : interactions) {
            Protein first = i.getFirstProtein();
            Protein second = i.getSecondProtein();
            ids.add(first.getPrimaryAccession());
            ids.add(second.getPrimaryAccession());
        }
        return ids;
    }
    
    // This method may be moved to other places
    public static Map<String, Set<String>> switchKeyValues(Map<String, Set<String>> map) {
        Map<String, Set<String>> rtn = new HashMap<String, Set<String>>();
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Set<String> values = map.get(key);
            for (String value : values) {
                Set<String> keys = rtn.get(value);
                if (keys == null) {
                    keys = new HashSet<String>();
                    rtn.put(value, keys);
                }
                keys.add(key);
            }
        }
        return rtn;
    }
    
    public static <T> void generateCytoscapeAttributeFile(Map<String, T> geneToValue,
                                                          String attName,
                                                          Class<T> cls,
                                                          String fileName) throws IOException {
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        fu.printLine(attName + " (class=" + cls.getName() + ")");
        for (String gene : geneToValue.keySet())
            fu.printLine(gene + "=" + geneToValue.get(gene));
        fu.close();
    }
    
    public static Map<String, String> swapKeyValue(Map<String, String> map) {
        Map<String, String> rtn = new HashMap<String, String>();
        for (String key : map.keySet()) {
            String value = map.get(key);
            rtn.put(value, key);
        }
        return rtn;
    }
    
    public static Map<String, Integer> countTermUsageInList(List<String> list) {
        Map<String, Integer> termToNumber = new HashMap<String, Integer>();
        for (String term : list) {
            Integer c = termToNumber.get(term);
            if (c == null) 
                termToNumber.put(term, 1);
            else
                termToNumber.put(term, ++c);
        }
        return termToNumber;
    }
    
    public static Map<String, Integer> generateProteinToDegree(Set<String> fis) {
        Map<String, Set<String>> proteinToPartners = generateProteinToPartners(fis);
        Map<String, Integer> proteinToDegree = new HashMap<String, Integer>();
        for (String protein : proteinToPartners.keySet()) {
            Set<String> partners = proteinToPartners.get(protein);
            proteinToDegree.put(protein, partners.size());
        }
        return proteinToDegree;
    }
    
    public static Map<String, Set<String>> generateProteinToPartners(Set<String> fis) {
        Map<String, Set<String>> proteinToPartners = new HashMap<String, Set<String>>();
        int index = 0;
        for (String fi : fis) {
            index = fi.indexOf("\t");
            if (index == -1)
                index = fi.indexOf(" ");
            String id1 = fi.substring(0, index);
            String id2 = fi.substring(index + 1);
            Set<String> partners1 = proteinToPartners.get(id1);
            if (partners1 == null) {
                partners1 = new HashSet<String>();
                proteinToPartners.put(id1, partners1);
            }
            partners1.add(id2);
            Set<String> partners2 = proteinToPartners.get(id2);
            if (partners2 == null) {
                partners2 = new HashSet<String>();
                proteinToPartners.put(id2, partners2);
            }
            partners2.add(id1);
        }
        return proteinToPartners;
    }
    
    public static String getPathwayDBSourceLetter(GKInstance dataSource) throws Exception {
        if (dataSource == null)
            return "R"; // For Reactome
        String displayName = dataSource.getDisplayName();
        if (displayName.equals("pantherdb"))
            return "P";
        if (displayName.equals("INOH"))
            return "I";
        if (displayName.equals("The Cancer Cell Map"))
            return "C";
        if (displayName.equals("Pathway Interaction Database"))
            return "N";
        if (displayName.equals("BioCarta - Imported by PID"))
            return "B";
        if (displayName.equals("KEGG"))
            return "K";
        return "R"; // default
    }
    
    public static String removeQuotationMarks(String id) {
        int index1 = id.indexOf("\"");
        int index2 = id.lastIndexOf("\"");
        if (index1 == index2) {
            return id.substring(index1 + 1);
        }
        else
            return id.substring(index1 + 1, index2);
    }
    
}
