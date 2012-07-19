/*
 * Created on Sep 29, 2006
 *
 */
package org.reactome.data;

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
import org.reactome.data.ProteinSequenceHandler.Sequence;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.Protein;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.funcInt.ReactomeSourceType;

/**
 * As a wrapper for ReactomeAnalyzer to convert pair-wise interaction to
 * org.reactome.funInt.Interaction.
 * @author guanming
 */
public class ReactomeFuncInteractionExtractor {
    // wrapped ReactomeAnalyzer
    private ReactomeAnalyzer analyzer;
    // Cache these converted values
    private Map<String, Protein> accessionToProteinMap;
    private Map<String, Interaction> pairToInteractionMap;
    private Map<Long, ReactomeSource> idToSourceMap;
    private UniProtAnalyzer uniProtAnalyzer; 
    private Map<String, String> acIdMap; 
    private Set<String> uniSet;
    // Protein names extracted from databases.
    // Used to assign names to proteins
    private Map<String, Protein> dbAccessionToProteinInfoMap;
    // For sequence related information
    private ProteinSequenceHandler seqHandler;
    
    public ReactomeFuncInteractionExtractor() {
        resetMaps();
        seqHandler = new ProteinSequenceHandler();
    }
    
    public void setReactomeAnalyzer(ReactomeAnalyzer analyzer) {
        this.analyzer = analyzer;
    }
    
    private void resetMaps() {
        if (accessionToProteinMap == null)
            accessionToProteinMap = new HashMap<String, Protein>();
        else
            accessionToProteinMap.clear();
        if (pairToInteractionMap == null)
            pairToInteractionMap = new HashMap<String, Interaction>();
        else
            pairToInteractionMap.clear();
        if (idToSourceMap == null)
            idToSourceMap = new HashMap<Long, ReactomeSource>();
        else
            idToSourceMap.clear();
    }
    
    public List<Interaction> getExtractedInteractions() {
        List<Interaction> funIntList = new ArrayList<Interaction>(pairToInteractionMap.values());
        return funIntList;
    }
    
    public void extractFuncInteractions() throws Exception {
       // resetMaps();
        // Extract from gk_central reactions
        Collection reactions = analyzer.prepareReactions();
        Collection complexes = analyzer.prepareComplexes();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        Set<String> interactions = new HashSet<String>();
        long time1 = System.currentTimeMillis();
        int c = 0;
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            interactors.clear();
            interactions.clear();
            analyzer.extractInteractorsFromReaction(rxn, interactors);
            analyzer.generateInteractionsWithDBNames(interactors, interactions, rxn);
            pushPairToFuncInteractions(interactions, 
                                       rxn);
        }
        GKInstance complex = null;
        c = 0;
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            complex = (GKInstance) it.next();
            //System.out.println("Complex: " + c++ + " " + complex.getDBID());
            interactors.clear();
            interactions.clear();
            analyzer.grepComplexComponents(complex, interactors);
            analyzer.generateInteractionsWithDBNames(interactors, interactions, complex);
            pushPairToFuncInteractions(interactions, 
                                       complex);
        }
        // Something special for Interactions
        if (analyzer instanceof CPathAnalyzer) {
            CPathAnalyzer cpathAnalyzer = (CPathAnalyzer) analyzer;
            Collection eventInteractions = cpathAnalyzer.prepareInteractions();
            GKInstance eventInteraction = null;
            c = 0;
            for (Iterator it = eventInteractions.iterator(); it.hasNext();) {
                eventInteraction = (GKInstance) it.next();
                interactors.clear();
                if (eventInteraction.getSchemClass().isa(ReactomeJavaConstants.Interaction)) {
                    String interactionType = (String) eventInteraction.getAttributeValue(ReactomeJavaConstants.interactionType);
                    if (interactionType != null && interactionType.contains("missing interaction"))
                        continue;
                    List interactorList = eventInteraction.getAttributeValuesList(ReactomeJavaConstants.interactor);
                    if (interactorList != null && interactorList.size() > 0) {
                        for (Iterator it1 = interactorList.iterator(); it1.hasNext();) {
                            interactors.add((GKInstance)it1.next());
                        }
                    }
                }
                else if (eventInteraction.getSchemClass().isa(ReactomeJavaConstants.TargettedInteraction)) {
                    TargetedInteractionAnalyzer tiAnalyzer = (TargetedInteractionAnalyzer) analyzer;
                    if (!tiAnalyzer.isNeededInteraction(eventInteraction))
                        continue;
                    GKInstance factor = (GKInstance) eventInteraction.getAttributeValue(ReactomeJavaConstants.factor);
                    GKInstance target = (GKInstance) eventInteraction.getAttributeValue(ReactomeJavaConstants.target);
                    if (factor != null)
                        interactors.add(factor);
                    if (target != null)
                        interactors.add(target);
                }
                interactions.clear();
                cpathAnalyzer.generateInteractionsWithDBNames(interactors, 
                                                              interactions, 
                                                              eventInteraction);
                if (interactions.size() == 0) {
                    System.out.println("Cannot extract interactions: " + eventInteraction);
                    c ++;
                }
                pushPairToFuncInteractions(interactions, 
                                           eventInteraction);
            }
            System.out.println("Total empty interactions: " + c);
        }
        //List<Interaction> funIntList = new ArrayList<Interaction>(pairToInteractionMap.values());
        //return funIntList;
    }
    
    /**
     * Convert a string based accession pairs to Protein based Interaction objects.
     * @param interactions
     * @param source
     * @param pairToInteractionMap
     * @param accessionToProteinMap
     */
    private void pushPairToFuncInteractions(Set<String> interactions,
                                            GKInstance source) throws Exception {
        int index = 0;
        for (String pair : interactions) {
            Interaction funInt = getInteractionFromPair(pair);
            if (funInt == null)
                continue;
            ReactomeSource reactomeSrc = getReactomeSource(source);
            funInt.addReactomeSource(reactomeSrc);
        }
    }
    
    private ReactomeSource getReactomeSource(GKInstance instance) throws Exception {
        ReactomeSource src = idToSourceMap.get(instance.getDBID());
        if (src == null) {
            src = new ReactomeSource();
            src.setReactomeId(instance.getDBID());
            GKInstance dataSource = getDataSource(instance);
            if (dataSource == null)
                src.setDataSource("Reactome"); // default
            else
                src.setDataSource(dataSource.getDisplayName());
            if (instance.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent))
                src.setSourceType(ReactomeSourceType.REACTION);
            else if (instance.getSchemClass().isa(ReactomeJavaConstants.Complex))
                src.setSourceType(ReactomeSourceType.COMPLEX);
            else if (instance.getSchemClass().isa(ReactomeJavaConstants.Interaction))
                src.setSourceType(ReactomeSourceType.INTERACTION);
            else if (instance.getSchemClass().isa(ReactomeJavaConstants.TargettedInteraction))
                src.setSourceType(ReactomeSourceType.TARGETED_INTERACTION);
            idToSourceMap.put(instance.getDBID(), src);
        }
        return src;
    }
    
    private GKInstance getDataSource(GKInstance instance) throws Exception {
        GKInstance dataSource = null;
        if (instance.getSchemClass().isValidAttribute(ReactomeJavaConstants.dataSource)) {
            dataSource = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.dataSource);
        }
        return dataSource;
    }
    
    private Interaction getInteractionFromPair(String pair) throws Exception {
        pair = validatePair(pair);
        if (pair == null)
            return null;
        Interaction funInt = pairToInteractionMap.get(pair);
        if (funInt != null)
            return funInt;
        // Need to create one
        // Have to get first protein and second protein
        int index = pair.indexOf("\t"); // Use "\t" since dbName might contain " "
        // The order has been made when interaction pair is extracted and validated
        String firstAccession = pair.substring(0, index);
        Protein firstProtein = getProtein(firstAccession);
        String secondAccession = pair.substring(index + 1);
        Protein secondProtein = getProtein(secondAccession);
        funInt = new Interaction();
        funInt.setFirstProtein(firstProtein);
        funInt.setSecondProtein(secondProtein);
        pairToInteractionMap.put(pair, funInt);
        return funInt;
    }
    
    private String validatePair(String pair) throws Exception {
        if (uniProtAnalyzer == null) {
            uniProtAnalyzer = new UniProtAnalyzer();
            acIdMap = uniProtAnalyzer.loadUniProtIDsMap();
            uniSet = acIdMap.keySet();
        }
        // Have to make sure ids used in pair are from human
        int index = pair.indexOf("\t");
        String firstAccessionWithDbName = pair.substring(0, index);
        String secondAccessionWithDbName = pair.substring(index + 1);
        // Some errors in the database
        index = firstAccessionWithDbName.indexOf(":");
        String firstAccession = firstAccessionWithDbName.substring(index + 1);
        if (firstAccession.equals("-1"))
            return null;
        String firstDbName = firstAccessionWithDbName.substring(0, index);
        index = secondAccessionWithDbName.indexOf(":");
        String secondAccession = secondAccessionWithDbName.substring(index + 1);
        if (secondAccession.equals("-1"))
            return null;
        String secondDbName = secondAccessionWithDbName.substring(0, index);
        if (firstDbName.startsWith("UniProt")) {
            if (!uniProtAnalyzer.isHumanID(uniSet, firstAccession))
                return null;
            if (acIdMap.containsKey(firstAccession)) // Check to avoid getting null for alternative form
                firstAccession = acIdMap.get(firstAccession);
        }
        if (secondDbName.startsWith("UniProt")) {
            if (!uniProtAnalyzer.isHumanID(uniSet, secondAccession))
                return null;
            if (acIdMap.containsKey(secondAccession))
                secondAccession = acIdMap.get(secondAccession);
        }
        // Have to make sure all proteins having aa sequences
        // Some Entrez Proteins from BIND are not human proteins
        // Need to recreate dbaccession in case it is remapped from acIdMap
        firstAccessionWithDbName = firstDbName + ":" + firstAccession;
        Sequence sequence1 = seqHandler.getSequence(firstAccessionWithDbName);
        if (sequence1 == null)
            return null;
        secondAccessionWithDbName = secondDbName +":" + secondAccession;
        Sequence sequence2 = seqHandler.getSequence(secondAccessionWithDbName);
        if (sequence2 == null)
            return null;
        // Need to merge accessions from different sources but having the same AA sequences
        String seqDbAcc1 = seqHandler.getDbAccFromChecksum(sequence1.getChecksum());
        String seqDbAcc2 = seqHandler.getDbAccFromChecksum(sequence2.getChecksum());
        String[] parsed = parseDbAcc(seqDbAcc1);
        firstDbName = parsed[0];
        firstAccession = parsed[1];
        parsed = parseDbAcc(seqDbAcc2);
        secondDbName = parsed[0];
        secondAccession = parsed[1];
        // Parse dbName and accession
        String rtn = null;
        // Have to make sure the firstAccession is less than the secondAccession
        int compare = firstAccession.compareTo(secondAccession);
        if (compare < 0)
            rtn = firstDbName + ":" + firstAccession + "\t" + secondDbName + ":" + secondAccession;
        else if (compare > 0)
            rtn = secondDbName + ":" + secondAccession + "\t" + firstDbName + ":" + firstAccession;
        return rtn;
    }
    
    private String[] parseDbAcc(String dbAcc) {
        int index = dbAcc.indexOf(":");
        String[] rtn = new String[] {
                dbAcc.substring(0, index),
                dbAcc.substring(index + 1)
        };
        return rtn;
    }
    
    /**
     * Get a protein based on an accession number.
     * @param dbAccession
     * @return
     * @throws Exception
     */
    public Protein getProtein(String dbAccession) throws Exception {
        Protein protein = accessionToProteinMap.get(dbAccession);
        if (protein != null)
            return protein;
        protein = new Protein();
        // Need to parse dbname and accession
        int index = dbAccession.indexOf(":");
        String dbName = dbAccession.substring(0, index);
        String accession = dbAccession.substring(index + 1);
        protein.setPrimaryAccession(accession);
        protein.setPrimaryDbName(dbName);
        addNameToProtein(protein);
        addSequenceToProtein(protein, dbAccession);
        accessionToProteinMap.put(dbAccession, protein);
        return protein;
    }
    
    private void addSequenceToProtein(Protein protein,
                                      String dbAccession) throws Exception {
        Sequence sequence = seqHandler.getSequence(dbAccession);
        if (sequence == null) 
            throw new IllegalStateException("Cannot find sequence for " + dbAccession);
        protein.setCheckSum(sequence.getChecksum());
        protein.setSequence(sequence.getSequence());
    }
    
    /**
     * TODO: Protein names "-", "Hypothetical protein" should be nullified.
     * @param protein
     * @throws Exception
     */
    public void addNameToProtein(Protein protein) throws Exception {
        Map<String, Protein> accToProtein = getDbAccToProteinMap();
        String dbName = protein.getPrimaryDbName();
        String key = generateProteinKey(protein);
        if (accToProtein.containsKey(key)) {
            Protein infoProtein = accToProtein.get(key);
            // Some names from BIND are actually protein description.
            if (infoProtein.getName().length() > 100)
                protein.setName(infoProtein.getName().substring(0, 100));
            else
                protein.setName(infoProtein.getName());
            protein.setShortName(infoProtein.getShortName());
        }
    }
    
    protected String generateProteinKey(Protein protein) {
        if (!protein.getPrimaryDbName().equals("UniProt"))
            return protein.getPrimaryDbName() + ":" + protein.getPrimaryAccession();
        String acc = protein.getPrimaryAccession();
        // To avoid alternative forms
        int index = acc.indexOf("-");
        if (index > 0)
            acc = acc.substring(0, index);
        return protein.getPrimaryDbName() + ":" + acc;
    }
    
    protected Map<String, Protein> getDbAccToProteinMap() throws Exception {
        if (dbAccessionToProteinInfoMap != null)
            return dbAccessionToProteinInfoMap;
        // Load the UniProt part
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, Protein> uniToProtein = uniAnalyzer.generateUniAccToProteinMap();
        dbAccessionToProteinInfoMap = new HashMap<String, Protein>();
        for (Iterator<String> it = uniToProtein.keySet().iterator(); it.hasNext();) {
            String acc = it.next();
            Protein protein = uniToProtein.get(acc);
            dbAccessionToProteinInfoMap.put(protein.getPrimaryDbName() + ":" + protein.getPrimaryAccession(),
                                            protein);
        }
        // The following mapping will not be used in the new version as of April 15, 2009.
//        // Load the Entrez Protein part from reactome_plus_i
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_plus_i_v2",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
//        Long entrezProteinId = new Long(379957L);
//        GKInstance entrezProteinDb = dba.fetchInstance(entrezProteinId);
//        Collection entrezProteins = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferencePeptideSequence, 
//                                                                 ReactomeJavaConstants.referenceDatabase,
//                                                                 "=",
//                                                                 entrezProteinDb);
//        for (Iterator it = entrezProteins.iterator(); it.hasNext();) {
//            GKInstance entrezProtein = (GKInstance) it.next();
//            List names = entrezProtein.getAttributeValuesList(ReactomeJavaConstants.name);
//            Protein protein = new Protein();
//            protein.setPrimaryDbName("Entrez Protein");
//            protein.setPrimaryAccession((String)entrezProtein.getAttributeValue(ReactomeJavaConstants.identifier));
//            if (names != null) {
//                if (names.size() == 1) {
//                    protein.setName(names.get(0).toString());
//                    protein.setShortName(names.get(0).toString());
//                }
//                else {
//                    // Find the shortest and longest name
//                    String shortName = (String) names.get(0);
//                    String longName = (String) names.get(0);
//                    for (int i = 1; i < names.size(); i++) {
//                        String name = (String) names.get(i);
//                        if (name.length() > longName.length())
//                            longName = name;
//                        else if (name.length() < shortName.length())
//                            shortName = name;
//                    }
//                    protein.setName(longName);
//                    protein.setShortName(shortName);
//                }
//            }
//            dbAccessionToProteinInfoMap.put(protein.getPrimaryDbName() + ":" + protein.getPrimaryAccession(),
//                                            protein);
//        }
//        // Load HPRD part
//        HPRDAnalyzer hprdAnalyzer = new HPRDAnalyzer();
//        List<Protein> proteins = hprdAnalyzer.loadHPRDWithNames();
//        for (Protein protein : proteins) {
//            dbAccessionToProteinInfoMap.put(protein.getPrimaryDbName() + ":" + protein.getPrimaryAccession(),
//                                            protein);
//        }
        return dbAccessionToProteinInfoMap;
    }    
}
