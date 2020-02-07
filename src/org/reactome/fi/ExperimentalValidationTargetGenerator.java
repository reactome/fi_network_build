/*
 * Created on Mar 20, 2008
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.data.EnsemblAnalyzer;
import org.reactome.data.HPRDAnalyzer;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.HibernateUtil;

/**
 * This method is used to generate a list of FIs for experimental validation. Most of methods
 * in this class are targeted for IFNg pathway interactions.
 * @author wgm
 *
 */
@SuppressWarnings("unchecked")
public class ExperimentalValidationTargetGenerator {

    private FileUtility fu;
    private SessionFactory sessionFactory;
    
    public ExperimentalValidationTargetGenerator() {
        fu = new FileUtility();
    }
    
    private MySQLAdaptor getDBA() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_plus_i_v2",
                                            "root",
                                            "macmysql01",
                                            3306);
        return dba;
    }
    
    private void initSession() throws Exception {
        String configFileName = "resources/funcIntHibernate.cfg.xml";
        File configFile = new File(configFileName);
        sessionFactory = HibernateUtil.getSessionFactory(configFile);
    }
    
    private Set<String> getPathwayIds() throws Exception {
        ReactomeAnalyzer analyzer = new ReactomeAnalyzer();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            "reactome_plus_i_v2",
                                            "root",
                                            "macmysql01",
                                            3306);
        // The following lists pathways related to JAK/STAT1 and IFN gamma pathway
        Long[] pathwayIds = new Long[] {
                // Positive regulation of (Transcription of SOCS by STAT dimer) in JAK STAT pathway [Canonical] (I)
                230409L, 
                230398L,
                191873L,
                438494L,
                230364L,
                191870L,
                425015L
        };
        Set<String> ids = new HashSet<String>();
        for (Long pathwayId : pathwayIds) {
            GKInstance pathway = dba.fetchInstance(pathwayId);
            Set<String> proteinIds = analyzer.grepIDsFromTopic(pathway);
            System.out.println("Total protein in pathway: " + proteinIds.size());
            ids.addAll(proteinIds);
        }
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, String> acIdMap = uniProtAnalyzer.loadUniProtIDsMap();
        Set<String> totalSwissIds = uniProtAnalyzer.loadSwissProtIds();
        Set<String> mappedIds = new HashSet<String>();
        for (String id : ids) {
            if (acIdMap.containsKey(id)) {
                String mapped = acIdMap.get(id);
                mappedIds.add(mapped);
            }
            else
                mappedIds.add(id);
        }
        System.out.println("Before map: " + ids.size());
        System.out.println("After map: " + mappedIds.size());
        return mappedIds;
    }
    
    private Protein fetchProtein(Session session,
                                 String id) throws Exception {
        String queryStr = "SELECT p from Protein p where p.accession = ?";
        Query query = session.createQuery(queryStr);
        query.setParameter(0, id);
        List result = query.list();
        if (result == null || result.size() == 0)
            return null;
        return (Protein) result.iterator().next();
    }
    
    @Test
    public void generateNegativeControlForPathway() throws Exception {
        Set<String> mappedIds = getPathwayIds();
        // Try to use 20 ids
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Set<String> totalSwissIds = uniProtAnalyzer.loadSwissProtIds();
        Set<String> randomIds = randomPickIds(totalSwissIds,
                                              mappedIds,
                                              20);
        // Check these two ids information
        // We need a hibernate session
        initSession();
        Session session = sessionFactory.openSession();
        Map<String, Protein> idToProtein = uniProtAnalyzer.generateUniAccToProteinMap();
        Map<String, String> idToRefSeq = uniProtAnalyzer.loadUniProtToRefSeqMap();
        // Only keep proteins that have DNA array data
//        System.out.println("Filtering based on gene expression...");
//        Map<String, Integer> geneNameToExp = generateGeneExpData();
//        Map<String, String> idToName = new HashMap<String, String>();
//        for (String id : randomIds){
//            Protein protein = fetchProtein(session, id);
//            if (protein != null)
//                idToName.put(id, protein.getLabel());
//        }
//        filterBasedOnGeneExp(randomIds,
//                             idToName,
//                             geneNameToExp);
        EnsemblAnalyzer ensembler = new EnsemblAnalyzer();
        Map<String, Integer> idToGeneNumber = ensembler.getGeneCopyNumbers(new ArrayList<String>(randomIds));
        // Filter based on member size
        Map<String, Integer> idToInteractionNumber = new HashMap<String, Integer>();
        for (String id : randomIds) {
            Integer geneCopy = idToGeneNumber.get(id);
            // Don't want any gene number > 4
            if (geneCopy == null || geneCopy > 4)
                continue;
            Set<Interaction> interactions = fetchInteractions(session, 
                                                              id, 
                                                              0.0d,
                                                              false);
            // Don't want to have this protein interacting with any proteins in the pathway
            if (isInteractingWithPathway(mappedIds, interactions))
                continue;
            // Don't want any proteins which can interact with more than ten proteins.
            Protein protein = idToProtein.get(id);
            String refSeq = idToRefSeq.get(id);
//            Integer geneExp = getGeneExp(id, 
//                                         idToName,
//                                         geneNameToExp);
            System.out.println(id + "\t" + refSeq + "\t" + geneCopy + "\t" + interactions.size() + "\t" + 
                               protein.getShortName() + "\t" + protein.getName());
        }
        session.close();
    }
    
    private boolean isInteractingWithPathway(Set<String> pathwayIds,
                                             Set<Interaction> interactions) {
        for (Interaction interaction : interactions) { 
            String id1 = interaction.getFirstProtein().getPrimaryAccession();
            String id2 = interaction.getSecondProtein().getPrimaryAccession();
            if (pathwayIds.contains(id1) || pathwayIds.contains(id2))
                return true;
        }
        return false;
    }
    
    private Set<String> randomPickIds(Set<String> totalIds,
                                      Set<String> excludedIds,
                                      int size) {
        Set<String> rtn = new HashSet<String>();
        List<String> list = new ArrayList<String>(totalIds);
        int index = 0;
        String id = null;
        while (rtn.size() < size) {
            index = (int) (Math.random() * list.size());
            id = list.get(index);
            if (excludedIds.contains(id))
                continue;
            rtn.add(id);
        }
        return rtn;
    }
    
    /*
     * This method is used to pull out FI partners for passed set of UniProt ids. The returned ids
     * should be in UniProt ids.
     */
    public Set<String> generateNewInteractionPartners(Set<String> targetIds) throws Exception {
        MySQLAdaptor dba = getDBA();
        boolean usePredicatedFIsOnly = false;
        double cutoff = 0.73d;
        initSession();
        Session session = sessionFactory.openSession();
        EnsemblAnalyzer ensemblAnalyzer = new EnsemblAnalyzer();
        Map<String, Integer> geneNameToExp = generateGeneExpData();
        Set<String> interactors = new HashSet<String>();
        for (String id : targetIds) {
            Map<String, Double> idToScore = new HashMap<String, Double>();
            // Want to print out the interactors:
            int c = 0;
            Set<String> partners = new HashSet<String>();
            Map<String, Integer> idToSetNumber = new HashMap<String, Integer>();
            String partnerId = null;
            Protein partner = null;
            Map<String, String> idToName = new HashMap<String, String>();
            Set<Interaction> interactions = fetchInteractions(session, 
                                                              id,
                                                              cutoff, 
                                                              usePredicatedFIsOnly);
            for (Interaction i : interactions) {
                String id1 = i.getFirstProtein().getPrimaryAccession();
                String id2 = i.getSecondProtein().getPrimaryAccession();
                if (id1.equals(id)) {
                    partnerId = id2;
                    partner = i.getSecondProtein();
                }
                else if (id2.equals(id)) {
                    partnerId = id1;
                    partner = i.getFirstProtein();
                }
                if (partner == null)
                    continue;
                partners.add(partnerId);
                idToScore.put(partnerId, i.getEvidence() == null ? 1.0 : i.getEvidence().getProbability());
                int setMemNumber = checkMemberNumbersInSet(dba, i, partnerId);
                idToSetNumber.put(partnerId, setMemNumber);
                idToName.put(partnerId, partner.getLabel());
            }
            System.out.println("Total partners: " + partners.size());
            Map<String, Integer> idToGeneNumber = ensemblAnalyzer.getGeneCopyNumbers(new ArrayList<String>(partners));
            System.out.println("Done gene copy number fetching...");
            filterBasedOnNumbers(idToGeneNumber, 
                                 partners, 
                                 8);
            // Filter based on member size
            System.out.println("Filtering based on entity set number...");
            filterBasedOnNumbers(idToSetNumber, 
                                 partners, 
                                 8);
            Map<String, Integer> idToInteractionNumber = new HashMap<String, Integer>();
            for (String par : partners) {
                // Want to check how many interactions these partners can have
                interactions = fetchInteractions(session, 
                                                par,
                                                0.0,
                                                usePredicatedFIsOnly); // want to check all in this case
                idToInteractionNumber.put(par, interactions.size());
            }
            System.out.println("Filtering based on interaction number...");
            // 5% percent should be 169. Here we expand the hub to bigger number
            filterBasedOnNumbers(idToInteractionNumber, partners, 200);
            // Only keep proteins that have DNA array data
            System.out.println("Filtering based on gene expression...");
            filterBasedOnGeneExp(partners,
                                 idToName,
                                 geneNameToExp);
            interactors.addAll(partners);
        }
        session.close();
        return interactors;
    }
                                                      
    @Test
    @SuppressWarnings("unchecked")
    public void generateNewInteractionPartners() throws Exception {  
        MySQLAdaptor dba = getDBA();
        boolean usePredicatedFIsOnly = false;
        String[] queryIds = new String[] {
                "P42224", // STAT1
                "P38484", // Interferon gamma receptor beta chain
                "P15260", // Interferon gamma receptor alpha chain
                "P23458", // JAK1
                "O60674"  // JAK2
        };
        double cutOff = 0.73;
        Set<String> pathwayIds = getPathwayIds();
        initSession();
        Session session = sessionFactory.openSession();
        EnsemblAnalyzer ensemblAnalyzer = new EnsemblAnalyzer();
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, Protein> idToProtein = uniProtAnalyzer.generateUniAccToProteinMap();
        Map<String, String> uniProtToRefSeq = uniProtAnalyzer.loadUniProtToRefSeqMap();
        Map<String, String> hprdToUniProt = fu.importMap(HPRDAnalyzer.HPRD_DIR + "HPRD2UniProtViaGeneSymbols.txt");
        Map<String, Integer> geneNameToExp = generateGeneExpData();
        for (String id : queryIds) {
            Set<Interaction> interactions = fetchInteractions(session, 
                                                              id,
                                                              0.73,
                                                              usePredicatedFIsOnly);
            Map<String, Double> idToScore = new HashMap<String, Double>();
            // Want to print out the interactors:
            int c = 0;
            Set<String> partners = new HashSet<String>();
            Map<String, Integer> idToSetNumber = new HashMap<String, Integer>();
            String partnerId = null;
            Protein partner = null;
            Map<String, String> idToName = new HashMap<String, String>();
            for (Interaction i : interactions) {
                String id1 = i.getFirstProtein().getPrimaryAccession();
                String id2 = i.getSecondProtein().getPrimaryAccession();
                if (id1.equals(id)) {
                    if (hprdToUniProt.containsKey(id2))
                        id2 = hprdToUniProt.get(id2);
                    partnerId = id2;
                    partner = i.getSecondProtein();
                }
                else {
                    if (hprdToUniProt.containsKey(id1))
                        id1 = hprdToUniProt.get(id1);
                    partnerId = id1;
                    partner = i.getFirstProtein();
                }
                partners.add(partnerId);
                idToScore.put(partnerId, i.getEvidence() == null ? 1.0 : i.getEvidence().getProbability());
                int setMemNumber = checkMemberNumbersInSet(dba, i, partnerId);
                idToSetNumber.put(partnerId, setMemNumber);
                idToName.put(partnerId, partner.getLabel());
            }
            partners.removeAll(pathwayIds); // Don't want to pick anything from its own pathway.
            System.out.println("Total partners: " + partners.size());
            Map<String, Integer> idToGeneNumber = ensemblAnalyzer.getGeneCopyNumbers(new ArrayList<String>(partners));
            System.out.println("Done gene copy number fetching...");
            filterBasedOnNumbers(idToGeneNumber, 
                                 partners, 
                                 8);
            // Filter based on member size
            System.out.println("Filtering based on entity set number...");
            filterBasedOnNumbers(idToSetNumber, 
                                 partners, 
                                 8);
            Map<String, Integer> idToInteractionNumber = new HashMap<String, Integer>();
            for (String par : partners) {
                // Want to check how many interactions these partners can have
                interactions = fetchInteractions(session, 
                                                par,
                                                0.0,
                                                usePredicatedFIsOnly); // want to check all in this case
                idToInteractionNumber.put(par, interactions.size());
            }
            System.out.println("Filtering based on interaction number...");
            // 5% percent should be 169. Here we expand the hub to bigger number
            filterBasedOnNumbers(idToInteractionNumber, partners, 200);
            // Only keep proteins that have DNA array data
            System.out.println("Filtering based on gene expression...");
            filterBasedOnGeneExp(partners,
                                 idToName,
                                 geneNameToExp);
            //Map<String, Integer> idToGeneNumber = new HashMap<String, Integer>();
            StringBuilder builder = new StringBuilder();
            builder.append(id).append(": ").append(partners.size()).append("\n");
            for (String par : partners) {
                // Want to check how many interactions these partners can have
                interactions = fetchInteractions(session, 
                                                par,
                                                0.0,
                                                usePredicatedFIsOnly); // want to check all in this case
                Integer geneNumber = idToGeneNumber.get(par);
                Integer setMemNumber = idToSetNumber.get(par);
                Protein protein = idToProtein.get(par);
                String refSeq = uniProtToRefSeq.get(par);
                if (refSeq != null && refSeq.length() == 0)
                    refSeq = null;
                double score = idToScore.get(par);
                Integer geneExp = getGeneExp(par,
                                             idToName, 
                                             geneNameToExp);
                builder.append("\t" + par + "\t" + refSeq + "\t" + score + "\t" + setMemNumber + "\t"  + geneExp + "\t" +  geneNumber + "\t" + interactions.size());
                if (protein != null)
                    builder.append("\t" + protein.getShortName() + "\t" + protein.getName());
                builder.append("\n");
            }
            System.out.println(builder.toString());
        }
        session.close();
    }
    
    private Integer getGeneExp(String id,
                               Map<String, String> idToName,
                               Map<String, Integer> nameToExp) {
        String name = idToName.get(id);
        return nameToExp.get(name);
    }
    
    private void filterBasedOnGeneExp(Set<String> partners,
                                      Map<String, String> idToName,
                                      Map<String, Integer> geneNameToExp) throws Exception {
        System.out.println("Before gene exp filter: " + partners.size());
        for (Iterator<String> it = partners.iterator(); it.hasNext();) {
            String id = it.next();
            String name = idToName.get(id);
            if (!geneNameToExp.containsKey(name))
                it.remove();
        }
        System.out.println("After gene exp filter: " + partners.size());
    }
    
    private int checkMemberNumbersInSet(MySQLAdaptor dba,
                                        Interaction interaction,
                                        String id) throws Exception {
        Set<ReactomeSource> sources = interaction.getReactomeSources();
        Set<GKInstance> instances = new HashSet<GKInstance>();
        for (ReactomeSource src : sources) {
            long dbId = src.getReactomeId();
            GKInstance instance = dba.fetchInstance(dbId);
            instances.add(instance);
        }
        // Grep all DefinedSet from GKInstance
        Set<GKInstance> entitySets = new HashSet<GKInstance>();
        for (GKInstance instance : instances) {
            getEntitySetFromInstance(instance, entitySets);
        }
        int number = 0; // Pick up the largest number
        for (GKInstance entitySet : entitySets) {
            Set<GKInstance> ewases = getEWASFromEntitySet(entitySet);
            // Check identifier
            for (GKInstance ewas : ewases) {
                GKInstance refEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
                if (refEntity == null)
                    continue;
                String identifier = (String) refEntity.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier.endsWith(id)) {
                    if (number < ewases.size()) {
                        number = ewases.size();
                    }
                    break;
                }
            }
        }
        return number;
    }
    
    private void filterBasedOnNumbers(Map<String, Integer> idToNumber,
                                      Set<String> partners,
                                      int cutOff) {
        System.out.println("Before filtering based on number: " + partners.size());
        for (Iterator<String> it = partners.iterator(); it.hasNext();) {
            String id = it.next();
            Integer number = idToNumber.get(id);
            if (number > cutOff)
                it.remove();
        }
        System.out.println("After filtering based on number: " + partners.size());
    }
    
    @SuppressWarnings("unchecked") 
    private Set<GKInstance> getEWASFromEntitySet(GKInstance entitySet) throws Exception {
        Set<GKInstance> ewasSet = new HashSet<GKInstance>();
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(entitySet);
        Set<GKInstance> next = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (GKInstance set : current) {
                if (set.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember)) {
                    List values = set.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                    for (Iterator it = values.iterator(); it.hasNext();) {
                        GKInstance value = (GKInstance) it.next();
                        if (value.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                            ewasSet.add(value);
                        else if (value.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                            next.add(value);
                    }
                }
                if (set.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasCandidate)) {
                    List values = set.getAttributeValuesList(ReactomeJavaConstants.hasCandidate);
                    for (Iterator it = values.iterator(); it.hasNext();) {
                        GKInstance value = (GKInstance) it.next();
                        if (value.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                            ewasSet.add(value);
                        else if (value.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                            next.add(value);
                    }
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
        return ewasSet;
    }
                                                 
    
    @SuppressWarnings("unchecked")
    private void getEntitySetFromInstance(GKInstance instance,
                                          Set<GKInstance> entitySet) throws Exception {
        Set<GKInstance> current = new HashSet<GKInstance>();
        if (instance.getSchemClass().isa(ReactomeJavaConstants.Reaction)) {
            Set participants = InstanceUtilities.getReactionParticipants(instance);
            for (Iterator it = participants.iterator(); it.hasNext();) {
                current.add((GKInstance) it.next());
            }
        }
        else if (instance.getSchemClass().isa(ReactomeJavaConstants.Complex))
            current.add(instance);
        Set<GKInstance> next = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (GKInstance tmp : current) {
                if (tmp.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                    entitySet.add(tmp);
                else if (tmp.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
                    List components = tmp.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                    next.addAll(components);
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    private Set<Interaction> fetchInteractions(Session session,
                                               String id,
                                               double cutOff,
                                               boolean usePredicatedFIsOnly) throws Exception {
        String queryStr1 = "SELECT i from Interaction i where i.firstProtein.primaryDbReference.accession = ?";
        String queryStr2 = "SELECT i from Interaction i where i.secondProtein.primaryDbReference.accession = ?";
        Set<Interaction> interactions = new HashSet<Interaction>();
        Query query = session.createQuery(queryStr1);
        query.setParameter(0, id);
        List list = query.list();
        for (Object obj : list) {
            interactions.add((Interaction)obj);
        }
        query = session.createQuery(queryStr2);
        query.setParameter(0, id);
        list = query.list();
        for (Object obj : list) {
            interactions.add((Interaction)obj);
        }
        // Do some filters
        if (usePredicatedFIsOnly) {
            for (Iterator<Interaction> it = interactions.iterator(); it.hasNext();) {
                Interaction in = it.next();
                if (in.getEvidence() == null)
                    it.remove();
                if (in.getEvidence().getProbability() < cutOff)
                    it.remove();
            }
        }
        else {
            for (Iterator<Interaction> it = interactions.iterator(); it.hasNext();) {
                Interaction in = it.next();
                if (in.getEvidence() == null)
                    continue;
                if (in.getEvidence().getProbability() < cutOff)
                    it.remove();
            }
        }
        return interactions;
    }
    
    /**
     * This method is used to generate expression data.
     * @throws IOException
     */
    public Map<String, Integer> generateGeneExpData() throws IOException {
        String dirName = "datasets/Rod_IFNG/";
        String[] fileNames = new String[] {
                "sht007.txt",
                "sht010.txt",
                "shx062.txt"
        };
        Map<String, List<Integer>> nameToValues = new HashMap<String, List<Integer>>();
        FileUtility fu = new FileUtility();
        String line = null;
        for (String fileName : fileNames) {
            fu.setInput(dirName + fileName);
            line = fu.readLine();
            while ((line = fu.readLine()) != null) {
                String[] tokens = line.split("\t");
                String gene = tokens[3];
                // Color1: CH1D_MEAN
                String value1 = tokens[17]; 
                // Color2: CH2DN_MEAN
                String value2 = tokens[28];
                List<Integer> values = nameToValues.get(gene);
                if (values == null) {
                    values = new ArrayList<Integer>();
                    nameToValues.put(gene, values);
                }
                if (value1.length() > 0)
                    values.add(new Integer(value1));
                if (value2.length() > 0)
                    values.add(new Integer(value2));
            }
            fu.close();
        }
        // calculate average values
        Map<String, Integer> nameToAverage = new HashMap<String, Integer>();
        for (Iterator<String> it = nameToValues.keySet().iterator(); it.hasNext();) {
            String name = it.next();
            List<Integer> values = nameToValues.get(name);
            int total = 0;
            for (Integer value : values) 
                total += value;
            int average = total / values.size();
            nameToAverage.put(name, average);
        }
        return nameToAverage;
    }
}
