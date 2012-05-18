/*
 * Created on May 3, 2006
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.util.*;

import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;

@SuppressWarnings("unchecked")
public class ReactomeAnalyzer {
    // Used to control if complex should be used
    protected boolean excludeComplex = false;
    protected PersistenceAdaptor dba;
    protected Long dataSourceId;
    // A helper class
    private ReactomeAnalyzerTopicHelper topicHelper;
    
    public ReactomeAnalyzer() {
        topicHelper = new ReactomeAnalyzerTopicHelper();
    }
    
    /**
     * Make sure dataSourceIds are correct for the version the database.
     * @return
     */
    public static List<ReactomeAnalyzer> getPathwayDBAnalyzers() {
        List<ReactomeAnalyzer> analyzers = new ArrayList<ReactomeAnalyzer>();
        analyzers.add(new ReactomeAnalyzer());
        PantherAnalyzer pantherAnalyzer = new PantherAnalyzer();
        Long pantherId = new Long(FIConfiguration.getConfiguration().get("PANTHER_DB_ID"));
        pantherAnalyzer.setDataSourceId(pantherId);
        analyzers.add(pantherAnalyzer);
        // INOH is not used in version 3.
        //analyzers.add(new INOHAnalyzer());
        Long[] dataSourceIds = new Long[] {
                new Long(FIConfiguration.getConfiguration().get("NCI_NATURE_CURATED_DB_ID")),
                new Long(FIConfiguration.getConfiguration().get("NCI_NATURE_BIOCARTA_DB_ID")),
                new Long(FIConfiguration.getConfiguration().get("KEGG_DB_ID")),
        };
        for (Long dataSourceId : dataSourceIds) {
            ReactomeAnalyzer tmp = new CPathAnalyzer();
            tmp.setDataSourceId(dataSourceId);
            analyzers.add(tmp);
        }
        // Add targeted interactions (TF/Target from TRED)
        TargetedInteractionAnalyzer tredAnalyzer = new TargetedInteractionAnalyzer();
        Long tredId = new Long(FIConfiguration.getConfiguration().get("TRED_DB_ID"));
        tredAnalyzer.setDataSourceId(tredId);
        analyzers.add(tredAnalyzer);

        TargetedInteractionAnalyzer encodeAnalyzer = new TargetedInteractionAnalyzer();
        encodeAnalyzer.setDataSourceId(new Long(FIConfiguration.getConfiguration().get("ENCODE_TFF_ID")));
        analyzers.add(encodeAnalyzer);
        return analyzers;
    }
    
    public static String getSourceLetter(ReactomeAnalyzer analyzer) throws Exception {
        GKInstance dataSource = analyzer.getDataSource();
        return InteractionUtilities.getPathwayDBSourceLetter(dataSource);
    }
    
    public void setExcludeComplex(boolean value) {
        this.excludeComplex = value;
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    protected PersistenceAdaptor getMySQLAdaptor() throws Exception {
        if (dba == null) {
//            dba = new MySQLAdaptor("localhost",
//                                   "gk_central_101606",//"panther_from_david", 
//                                   //"test_reactome_plus_i_v2",
//                                   "root",
//                                   "macmysql01",
//                                   3306);
//            dba = new MySQLAdaptor("localhost",
//                                   "reactome_28_plus_i",//"panther_from_david", 
//                                   "root",
//                                   "macmysql01",
//                                   3306);
            dba = new MySQLAdaptor("localhost",
                                   FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"), 
                                   FIConfiguration.getConfiguration().get("DB_USER"),
                                   FIConfiguration.getConfiguration().get("DB_PWD"),
                                   3306);
//            dba = new MySQLAdaptor("localhost",
//                                   "gk_central_031309",//"panther_from_david", 
//                                   "root",
//                                   "macmysql01",
//                                   3306);
        }
        return dba;
    }
    
    public void setDataSourceId(Long id) {
        this.dataSourceId = id;
    }
    
    public Set<String> generateUniProtPairsFromTopics() throws Exception {
        Set<String> pairs = new HashSet<String>();
        Map<GKInstance, Set<String>> topics2Ids = grepIDsFromTopics();
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            List<String> idList = new ArrayList<String>(ids);
            Collections.sort(idList);
            for (int i = 0; i < idList.size() - 1; i++) {
                String id1 = idList.get(i);
                for (int j = i + 1; j < idList.size(); j++) {
                    String id2 = idList.get(j);
                    pairs.add(id1 + " " + id2);
                }
            }
        }
        System.out.println("Total Pairs: " + pairs.size());
        return pairs;
    }
    
    /**
     * Using this method to generate interaction list for each topic.
     * @return
     * @throws Exception
     */
    @SuppressWarnings("unchecked")
    public Map<GKInstance, Set<String>> grepInteractionsForTopics() throws Exception {
        Map<GKInstance, Set<String>> topicToInteraction = new HashMap<GKInstance, Set<String>>();
        List<GKInstance> topics = getTopics();
        prepareReactions();
        prepareComplexes();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        for (GKInstance topic : topics) {
            Set<String> interactions = grepInteractionsForTopic(topic);
            topicToInteraction.put(topic, interactions);
        }
        return topicToInteraction;
    }
    
    /**
     * Use this method to load a list of FIs for a specified pathway (here is used as topic).
     * @param topic
     * @return
     * @throws Exception
     */
    public Set<String> grepInteractionsForTopic(GKInstance topic) throws Exception {
        Set<GKInstance> components = topicHelper.grepPathwayComponents(topic);
        Set<String> interactions = new HashSet<String>();
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        for (GKInstance comp : components) {
            if (comp.getSchemClass().isa(ReactomeJavaConstants.Reaction)) {
                interactors.clear();
                extractInteractorsFromReaction(comp, interactors);
                generateInteractions(interactors, interactions, comp);
            }
            else if (comp.getSchemClass().isa(ReactomeJavaConstants.Interaction)) {
                List list = comp.getAttributeValuesList(ReactomeJavaConstants.interactor);
                if (list != null) {
                    interactors.clear();
                    interactors.addAll(list);
                    generateInteractions(interactors, interactions, comp);
                }
            }
            else if (comp.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
                interactors.clear();
                grepComplexComponents(comp, interactors);
                generateInteractions(interactors, interactions, comp);
            }
        }
        return interactions;
    }
    
    public Set<String> generateUniProtIdsFromTopics() throws Exception {
        Map<GKInstance, Set<String>> topics2Ids = grepIDsFromTopics();
        Set<String> uniProtIds = new HashSet<String>();
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            uniProtIds.addAll(ids);
        }
        return uniProtIds;
    }
    
    protected List<GKInstance> getTopics() throws Exception {
        MySQLAdaptor releasedDBA = (MySQLAdaptor) getMySQLAdaptor();
        // The following code should be used for generating NBC training data set (May 11, 2009)
//        Collection frontPages = releasedDBA.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
//        List<GKInstance> topics = new ArrayList<GKInstance>();
//        for (Iterator it = frontPages.iterator(); it.hasNext();) {
//            GKInstance frontPage = (GKInstance) it.next();
//            List items = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
//            for (Iterator it1 = items.iterator(); it1.hasNext();)
//                topics.add((GKInstance)it1.next());
//        }
        // As of April 19, 2007, a list of Reactome pathways is fetched from a semi-manually
        // create file, ReactomePathways.txt.
        List<GKInstance> topics = new ArrayList<GKInstance>();
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        Long id = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            id = Long.parseLong(line.substring(0, index));
            GKInstance pathway = releasedDBA.fetchInstance(id);
            topics.add(pathway);
        }
        return topics;
    }
    
    public Map<GKInstance, Map<String, Integer>> grepIDNumberFromTopics() throws Exception {
        List<GKInstance> topics = getTopics();
        Map<GKInstance, Map<String, Integer>> topicToIDNumber = new HashMap<GKInstance, Map<String, Integer>>();
        for (GKInstance topic : topics) {
            Map<String, Integer> id2Number = topicHelper.grepIDToNumberFromTopic(topic);
            topicToIDNumber.put(topic, id2Number);
        }
        return topicToIDNumber;
    }
    
    public Map<GKInstance, Set<String>> grepIDsFromTopics() throws Exception {
        long time1 = System.currentTimeMillis();
        List<GKInstance> topics = getTopics();
        // Try to get ids for each topics
        Map<GKInstance, Set<String>> topics2Ids = new HashMap<GKInstance, Set<String>>();
        for (GKInstance topic : topics) {
            System.out.println("Topic: " + topic);
            Set<String> ids = grepIDsFromTopic(topic);
            topics2Ids.put(topic, ids);
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for grepping IDs for topics: " + (time2 - time1));
        // Print out the id numbers in each topic.
        for (Iterator<GKInstance> it = topics2Ids.keySet().iterator(); it.hasNext();) {
            GKInstance topic = it.next();
            Set<String> ids = topics2Ids.get(topic);
            System.out.println(topic.getDisplayName() + ": " + ids.size());
        }
        return topics2Ids;
    }
    
    @SuppressWarnings("unchecked")
    protected Set<GKInstance> grepPathwayParticipants(GKInstance pathway) throws Exception {
        return topicHelper.grepPathwayParticipants(pathway);
    }
    
    @SuppressWarnings("unchecked")
    public Set<String> grepIDsFromTopic(GKInstance topic) throws Exception {
        Set<String> ids = new HashSet<String>();
        // First load all PhysicalEntities involved in Reactions
        Set<GKInstance> participants = grepPathwayParticipants(topic);
        // Grep ReferencePeptideSequence
        for (GKInstance participant : participants) {
            if (excludeComplex && participant.getSchemClass().isa(ReactomeJavaConstants.Complex))
                continue;
            Set<GKInstance> refPeptides = grepRefPepSeqs(participant);
            for (GKInstance ref : refPeptides) {
                GKInstance refDb = (GKInstance) ref.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                // It should be faster by comparing DB_ID
                if (refDb == null || !refDb.getDBID().equals(2L))
                    continue;
                String identifier = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
                if (identifier != null)
                    ids.add(identifier);
            }
        }
        return ids;
    }
    
    private Collection loadReactions() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // Get homo sapiens
        GKInstance homosapiens = dba.fetchInstance(48887L);
        // Load all reactions for analyzed: human reactions only
        Collection reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                            ReactomeJavaConstants.species,
                                                            "=",
                                                            homosapiens);
        // Load precedingEvent values
        SchemaClass schema = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        SchemaAttribute att = schema.getAttribute(ReactomeJavaConstants.precedingEvent);
        dba.loadInstanceAttributeValues(reactions, att);
        return reactions;
    }
    
    @SuppressWarnings("unchecked")
    protected void extractInteractorsFromReaction(GKInstance rxn, 
                                                  Set<GKInstance> interactors) throws Exception {
        List input = rxn.getAttributeValuesList(ReactomeJavaConstants.input);
        if (input != null)
            interactors.addAll(input);
        // Get catalyst
        List cas = rxn.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            for (Iterator it1 = cas.iterator(); it1.hasNext();) {
                GKInstance ca = (GKInstance) it1.next();
                List catalysts = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                if (catalysts != null)
                    interactors.addAll(catalysts);
            }
        }
        // Check regulators
        Collection regulations = rxn.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null) {
            for (Iterator it1 = regulations.iterator(); it1.hasNext();) {
                GKInstance regulation = (GKInstance) it1.next();
                List regulators = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
                for (Iterator it2 = regulators.iterator(); it2.hasNext();) {
                    GKInstance regulator = (GKInstance) it2.next();
                    if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                        interactors.add(regulator);
                }
            }
        }
    }
    
    /**
     * This method is used to dump FIs extracted from the Reactome database to 
     * an external file.
     * @throws Exception
     */
    @Test
    public void extractInteractions() throws Exception {
        Set<String> interactionSet = extractInteractionSet();
        FileUtility fu = new FileUtility();
        //String fileName = "results/v2/ReactomeInteractions020507.txt";
        // 28 for the release number
        // Do a filters
        ProteinIdFilters filters = new ProteinIdFilters();
        System.out.println("Before filtering: " + interactionSet.size());
        Set<String> filtered = filters.cleanUpVsUniProt(interactionSet);
        System.out.println("After filtering: " + filtered.size());
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions28_051711.txt";
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions36.txt";
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
        fu.saveInteractions(filtered, fileName);
        // Check what have been removed
        // Note: interactions filtered out are extracted from reactions involed other non-human
        // species (e.g. virus, bacteria, etc).
//        interactionSet.removeAll(filtered);
//        int count = 0;
//        System.out.println("Total removed: " + interactionSet.size());
//        for (String interaction : interactionSet) {
//            System.out.println(interaction);
//            count ++;
//            if (count == 10)
//                break;
//        }
    }
    
    @Test
    public void checkNumberChangesBetweenTwoReleases() throws Exception {
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeOnlyInteractions28_051711.txt";
        FileUtility fu = new FileUtility();
        Set<String> fis1 = fu.loadInteractions(fileName);
        Set<String> fiProteins1 = InteractionUtilities.grepIDsFromInteractions(fis1);
        fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "ReactomeOnlyInteractions36.txt";
        Set<String> fis2 = fu.loadInteractions(fileName);
        Set<String> fiProteins2 = InteractionUtilities.grepIDsFromInteractions(fis2);
        System.out.println("Release 28 FIs: " + fis1.size());
        System.out.println("           Proteins: " + fiProteins1.size());
        System.out.println("Release 36 FIs: " + fis2.size());
        System.out.println("           Proteins: " + fiProteins2.size());
        // Check how many new reactions have been predicted before
        Set<String> totalFIs = fu.loadInteractions(FIConfiguration.getConfiguration().get("INTERACTION_FILE_NAME"));
        Set<String> totalFIProteins = InteractionUtilities.grepIDsFromInteractions(totalFIs);
        System.out.println("Total FIs: " + totalFIs.size());
        Set<String> newFIs = new HashSet<String>(fis2);
        newFIs.removeAll(fis1);
        Set<String> shared = InteractionUtilities.getShared(totalFIs, newFIs);
        Set<String> sharedProteins = InteractionUtilities.grepIDsFromInteractions(shared);
        System.out.println("Total new FIs: " + newFIs.size());
        System.out.println("    Shared FIs: " + shared.size());
        System.out.println("    Proteins: " + sharedProteins.size());
        Set<String> newProteins = new HashSet<String>(fiProteins2);
        newProteins.removeAll(fiProteins1);
        shared = InteractionUtilities.getShared(totalFIProteins, newProteins);
        System.out.println("Total new proteins: " + newProteins.size());
        System.out.println("Shared new Proteins: " + shared.size());
//        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
//        Collection<?> c = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct, 
//                                                                     ReactomeJavaConstants.species, 
//                                                                     "=",
//                                                                     48887L);
//        dba.loadInstanceAttributeValues(c, new String[]{ReactomeJavaConstants.identifier});
//        Set<String> idsInDb = new HashSet<String>();
//        for (Object obj : c) {
//            GKInstance inst = (GKInstance) obj;
//            String id = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
//            idsInDb.add(id);
//        }
//        System.out.println("Total ids in db: " + idsInDb.size());
//        Set<String> notInFIs = new HashSet<String>(idsInDb);
//        notInFIs.removeAll(fiProteins2);
//        System.out.println("Not in FIs: " + notInFIs.size());
//        int count = 0;
//        for (String id : notInFIs) {
//            System.out.println(id);
//            count ++;
//            if (count > 100)
//                break;
//        }
    }
    
    /**
     * Use this method to load a list of pre-generated FIs from the Reactome database.
     * @return
     * @throws IOException
     */
    public Set<String> loadFIsFromFile() throws IOException {
        FileUtility fu = new FileUtility();
        //String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "ReactomeInteractions28.txt";
//        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR + "FIs_Reactome.txt";
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_FI_FILE");
        return fu.loadInteractions(fileName);
    }
    
    public Set<String> extractInteractionSetWithComplexAndSet() throws Exception {
        // Extract from gk_central reactions
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        int c = 0;
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn);
            interactors.clear();
        }
        System.out.println("Total interactions from " + getDataSource() + ": "
                           + interactions.size());
        return interactions;
    }
    
    protected void generateInteractionsWithComplexAndSet(Set<GKInstance> interactors,
                                                         Set<String> interactions,
                                                         GKInstance reaction) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithComplexAndSet(interactor1, 
                                                      interactor2, 
                                                      interactions);
            }
        }
    }
    
    private void generateInteractionsWithComplexAndSet(GKInstance interactor1,
                                                       GKInstance interactor2,
                                                       Set<String> interactions) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs2.size() == 0)
            return;
        String int1 = convertRefPepSeqToString(refPepSeqs1, interactor1);
        if (int1.length() == 0)
            return;
        String int2 = convertRefPepSeqToString(refPepSeqs2, interactor2);
        if (int2.length() == 0)
            return;
        int compare = int1.compareTo(int2);
        if (compare < 0)
            interactions.add(int1 + " " + int2);
        else
            interactions.add(int2 + " " + int1);
    }
              
    private String convertRefPepSeqToString(Set<GKInstance> refPepSeqs,
                                            GKInstance interactor) throws Exception {
        List<String> ids = new ArrayList<String>();
        for (Iterator<GKInstance> it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = it.next();
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier == null)
                continue; // maybe
            ids.add(identifier);
        }
        Collections.sort(ids);
        String delimit = "?"; // default: should not be used
        if (interactor.getSchemClass().isa(ReactomeJavaConstants.Complex))
            delimit = ",";
        else if (interactor.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
            delimit = "|";
        StringBuilder builder = new StringBuilder();
        for (Iterator<String> it = ids.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append(delimit);
        }
        return builder.toString();
    }
    
    public Set<String> extractInteractionSet() throws Exception {
        // Extract from gk_central reactions
        Collection reactions = prepareReactions();
        Collection complexes = prepareComplexes();
        Set<String> interactions = new HashSet<String>();
        GKInstance rxn = null;
        Set<GKInstance> interactors = new HashSet<GKInstance>();
        long time1 = System.currentTimeMillis();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            //System.out.println("Reaction: " + c++);
            extractInteractorsFromReaction(rxn, interactors);
            generateInteractions(interactors, interactions, rxn);
            interactors.clear();
        }
        System.out.println("Total interactions from reactions: " + interactions.size());
        if (!excludeComplex) {
            GKInstance complex = null;
            for (Iterator it = complexes.iterator(); it.hasNext();) {
                complex = (GKInstance) it.next();
                //System.out.println("Complex: " + c++ + " " + complex.getDBID());
                interactors.clear();
                grepComplexComponents(complex, interactors);
                // No need
                //if (interactors.size() > 10)
                //    continue; // cutoff set manually
                generateInteractions(interactors, interactions, complex);
            }
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Time for looping: " + (time2 - time1));
        System.out.println("Total interactions from Reactome: " + interactions.size());
        return interactions;
    }
    
    protected void grepComplexComponents(GKInstance complex, Set<GKInstance> interactors) throws Exception {
        Set<GKInstance> current = new HashSet<GKInstance>();
        current.add(complex);
        Set<GKInstance> next = new HashSet<GKInstance>();
        while (current.size() > 0) {
            for (Iterator it = current.iterator(); it.hasNext();) {
                GKInstance tmp = (GKInstance) it.next();
                List components = tmp.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                if (components == null || components.size() == 0)
                    continue;
                for (Iterator it1 = components.iterator(); it1.hasNext();) {
                    GKInstance tmp1 = (GKInstance) it1.next();
                    if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.EntitySet))
                        interactors.add(tmp1);
                    else if (tmp1.getSchemClass().isa(ReactomeJavaConstants.Complex))
                        next.add(tmp1);
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    protected void generateInteractions(Set<GKInstance> interactors, 
                                        Set<String> interactions,
                                        GKInstance source) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractions(interactor1, interactor2, interactions, source);
            }
        }
    }
    
    protected void generateInteractionsWithDBNames(Set<GKInstance> interactors,
                                                   Set<String> interactions,
                                                   GKInstance source) throws Exception {
        List<GKInstance> list = new ArrayList<GKInstance>(interactors);
        int size = list.size();
        for (int i = 0; i < size - 1; i++) {
            GKInstance interactor1 = list.get(i);
            for (int j = i + 1; j < size; j++) {
                GKInstance interactor2 = list.get(j);
                generateInteractionsWithDBNames(interactor1, interactor2, interactions, source);
            }
        }        
    }
    
    private void generateInteractionsWithDBNames(GKInstance interactor1,
                                                 GKInstance interactor2,
                                                 Set<String> interactions,
                                                 GKInstance source) throws Exception {
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs1.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            String dbName1 = null;
            GKInstance db1 = (GKInstance) ref1.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            if (db1 != null)
                dbName1 = db1.getDisplayName();
            else
                dbName1 = "unknown";
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                String dbName2 = null;
                GKInstance db2 = (GKInstance) ref2.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                if (db2 != null)
                    dbName2 = db2.getDisplayName();
                else
                    dbName2 = "unknown";
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = dbName1 + ":" + uid1 + "\t" + dbName2 + ":" + uid2;
                else if (comp > 0)
                    pair = dbName2 + ":" + uid2 + "\t" + dbName1 + ":" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                }
            }
        }
    }
    
    private void generateInteractions(GKInstance interactor1, 
                                      GKInstance interactor2,
                                      Set<String> interactions,
                                      GKInstance source) throws Exception {
        if (excludeComplex) {
            if (interactor1.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
                interactor2.getSchemClass().isa(ReactomeJavaConstants.Complex))
                return;
        }
        Set<GKInstance> refPepSeqs1 = grepRefPepSeqs(interactor1);
        if (refPepSeqs1.size() == 0)
            return;
        Set<GKInstance> refPepSeqs2 = grepRefPepSeqs(interactor2);
        if (refPepSeqs1.size() == 0)
            return;
        // Permutate members in these two sets
        int comp = 0;
        String pair = null;
        for (GKInstance ref1 : refPepSeqs1) {
            String uid1 = (String) ref1.getAttributeValue(ReactomeJavaConstants.identifier);
            if (uid1 == null)
                continue;
            for (GKInstance ref2 : refPepSeqs2) {
                String uid2 = (String) ref2.getAttributeValue(ReactomeJavaConstants.identifier);
                if (uid2 == null)
                    continue;
                comp = uid1.compareTo(uid2);
                if (comp < 0)
                    pair = uid1 + "\t" + uid2;
                else if (comp > 0)
                    pair = uid2 + "\t" + uid1;
                if (pair != null) {
                    interactions.add(pair); //exclude self interaction
                    // Used for debugging
                    //if (pair.equals("O95405 P01270"))
                    //    System.out.println(pair + " < " + source.getDBID());
                }
            }
        }
    }
    
    protected Set<GKInstance> grepRefPepSeqs(GKInstance interactor) throws Exception {
        return topicHelper.grepRefPepSeqs(interactor);
    }
    
    protected Collection prepareComplexes() throws Exception {
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection complexes = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
                                                            ReactomeJavaConstants.species,
                                                            "=",
                                                            homosapiens);
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Complex);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.hasComponent);
        dba.loadInstanceAttributeValues(complexes, att);
        return complexes;
    }
    
    protected Collection prepareReactions() throws Exception {
        // Load necessary attributes
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        // Load all reactions for analyzed
        GKInstance homosapiens = dba.fetchInstance(48887L);
        Collection reactions = null;
        SchemaClass reactionCls = null;
        // Adjust for new schema
        if (dba.getSchema().isValidClass(ReactomeJavaConstants.ReactionlikeEvent)) {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.ReactionlikeEvent,
                                                    ReactomeJavaConstants.species,
                                                    "=",
                                                    homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.ReactionlikeEvent);
        }
        else {
            reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                     ReactomeJavaConstants.species,
                                                     "=",
                                                     homosapiens);
            reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
        }
        Collection cas = dba.fetchInstancesByClass(ReactomeJavaConstants.CatalystActivity);
        Collection regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
        Collection entities = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                           ReactomeJavaConstants.species,
                                                           "=",
                                                           homosapiens);
        SchemaAttribute att = reactionCls.getAttribute(ReactomeJavaConstants.input);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.output);
        dba.loadInstanceAttributeValues(reactions, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.catalystActivity);
        dba.loadInstanceAttributeValues(reactions, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.CatalystActivity);
        att = reactionCls.getAttribute(ReactomeJavaConstants.physicalEntity);
        dba.loadInstanceAttributeValues(cas, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.Regulation);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulatedEntity);
        dba.loadInstanceAttributeValues(regulations, att);
        att = reactionCls.getAttribute(ReactomeJavaConstants.regulator);
        dba.loadInstanceAttributeValues(regulations, att);
        reactionCls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        att = reactionCls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(entities, att);
        return reactions;
    }
    
    public void findReactionGraphComponents() throws Exception {
        Collection reactions = loadReactions();
        final Map<GKInstance, Set<GKInstance>> componentMap = new HashMap<GKInstance, Set<GKInstance>>();
        GKInstance rxn = null;
        Set<String> preTypes = new HashSet<String>();
        for (Iterator it = reactions.iterator(); it.hasNext();) {
            rxn = (GKInstance) it.next();
            Set<GKInstance> component = componentMap.get(rxn);
            if (component == null) {
                component = new HashSet<GKInstance>();
                component.add(rxn);
                componentMap.put(rxn, component);
            }
            List preceding = rxn.getAttributeValuesList(ReactomeJavaConstants.precedingEvent);
            if (preceding != null) {
                for (Iterator it1 = preceding.iterator(); it1.hasNext();) {
                    GKInstance preRxn = (GKInstance) it1.next();
                    preTypes.add(preRxn.getSchemClass().getName());
                    component.add(preRxn);
                    componentMap.put(preRxn, component);
                }
            }
        }
        Set<Set<GKInstance>> components = new HashSet<Set<GKInstance>>(componentMap.values());
        List<Set<GKInstance>> list = new ArrayList<Set<GKInstance>>(components);
        Collections.sort(list, new Comparator<Set<GKInstance>>() {
            public int compare(Set<GKInstance> set1, Set<GKInstance> set2) {
                return set2.size() - set1.size();
            }
        });
        System.out.printf("Total Components: %d -> %d%n", componentMap.size(), components.size());
        // Print out the first ten
        for (int i = 0; i < 10; i++) {
            Set<GKInstance> set = list.get(i);
            System.out.println(i + ": " + set.size());
        }
        System.out.println("Preceding Type: " + preTypes);
    }
    
    public void checkNoInteractionInReactome() throws IOException {
        String noIntFileName = "results/NoInteractionsForTrain.txt";
        String reactomeIntFileName = "results/ReactomeInteractions.txt";
        Set<String> noIntSet = new HashSet<String>();
        Set<String> rIntSet = new HashSet<String>();
        FileUtility fu = new FileUtility();
        fu.setInput(noIntFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            noIntSet.add(line);
        }
        fu.close();
        fu.setInput(reactomeIntFileName);
        while ((line = fu.readLine()) != null) {
            rIntSet.add(line);
        }
        fu.close();
        System.out.println("Reactome Interactions: " + rIntSet.size());
        int size1 = noIntSet.size();
        System.out.println("No Interactions: " + size1);
        noIntSet.removeAll(rIntSet);
        int size2 = noIntSet.size();
        System.out.println("Reactome in No Interactions: " + (size1 - size2));
    }
    
    public void countTotalUsedUniProtIDs() throws Exception {
        // Fetch human only
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor();
        Collection ews = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                      ReactomeJavaConstants.species,
                                                      "=",
                                                      48887);
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.referenceEntity);
        dba.loadInstanceAttributeValues(ews, att);
        Set<String> ids = new HashSet<String>();
        for (Iterator it = ews.iterator(); it.hasNext();) {
            GKInstance ew = (GKInstance) it.next();
            GKInstance ref = (GKInstance) ew.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (ref == null)
                continue;
            String id = (String) ref.getAttributeValue(ReactomeJavaConstants.identifier);
            if (id != null)
                ids.add(id);
        }
        System.out.println("Total IDs used in gk_central: " + ids.size());
    }

    public GKInstance getDataSource() throws Exception {
        if (dataSourceId == null)
            return null;
        // GKInstance for dataSource pantherdb
        GKInstance dataSource = null;
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) 
            dataSource = ((MySQLAdaptor)adaptor).fetchInstance(dataSourceId);
        else
            dataSource = ((XMLFileAdaptor)adaptor).fetchInstance(dataSourceId);
        return dataSource;
    }
    
    /**
     * @throws Exception
     */
    @Test
    public void countProteinsInPathways() throws Exception {
        // Top level pathways
        MySQLAdaptor dba = (MySQLAdaptor) getMySQLAdaptor ();
        Collection events = dba.fetchInstancesByClass(ReactomeJavaConstants.Event);
        List topLevelPathways = InstanceUtilities.grepTopLevelEvents(events);
        List<GKInstance> topics = getTopics();
        // Check proteins in each topics
        Set<String> totals = new HashSet<String>();
        for (GKInstance topic : topics) {
            Set<String> ids = grepIDsFromTopic(topic);
            System.out.println(topic.getDisplayName() + "(" + topic.getDBID() + "): " + ids.size());
            if (!topLevelPathways.contains(topic))
                System.out.println("\tNot a toplevel!");
            totals.addAll(ids);
            if (ids.size() > 200) {
                // check the subpathways
                List components = null;
                if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent)) 
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasSpecialisedForm))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasSpecialisedForm);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember))
                    components = topic.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                if (components == null || components.size() == 0)
                    continue;
                for (Iterator it = components.iterator(); it.hasNext();) {
                    GKInstance comp = (GKInstance) it.next();
                    Set<String> subIds = grepIDsFromTopic(comp);
                    System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + "): " + subIds.size());
                }
            }
        }
        System.out.println("Total Ids: " + totals.size());
        // Special handling on Signaling Pathways
        GKInstance signalingPathway = dba.fetchInstance(162582L);
        List components = signalingPathway.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        for (Iterator it = components.iterator(); it.hasNext();) {
            GKInstance comp = (GKInstance) it.next();
            Collection referrers = comp.getReferers(ReactomeJavaConstants.hasComponent);
            if (referrers.size() > 1) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
            referrers = comp.getReferers(ReactomeJavaConstants.hasMember);
            if (referrers != null && referrers.size() > 0) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
            referrers = comp.getReferers(ReactomeJavaConstants.hasSpecialisedForm);
            if (referrers != null && referrers.size() > 0) {
                System.out.println("\t" + comp.getDisplayName() + "(" + comp.getDBID() + ") contained by others");
                continue;
            }
        }
    }
    
    /**
     * Generate a list of pathways from the Reactome pathway hierarchy tree so that they can be used
     * in pathway enrichment analysis. The generation is done in two steps:
     * 1). Check item in the FrontPageItem list. If a pathway item has less than 200 proteins, that pathway
     * is listed in the pathway list.
     * 2). Otherwise, a front page pathway's sub-pathways are listed in the list regardless of their sizes if 
     * the size of subpathways are less than 300. If a sub-pathway has size > 300, its sub-pathways are listed.
     * 3). Note: The cutoff values 200 and 300 is chosen rather arbitrary. The reason why 300 is chosen for the
     * second level is based on assumption that the lower level pathways should be more like functional units
     * than upper level pathways.
     * @throws Exception
     */
    @Test
    public void generateListOfPathways() throws Exception {
        Set<GKInstance> pathwaySet = new HashSet<GKInstance>();
        MySQLAdaptor releasedDBA = (MySQLAdaptor) getMySQLAdaptor();
        Collection frontPages = releasedDBA.fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        List<GKInstance> bigTopics = new ArrayList<GKInstance>();
        for (Iterator it = frontPages.iterator(); it.hasNext();) {
            GKInstance frontPage = (GKInstance) it.next();
            List items = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
            for (Iterator it1 = items.iterator(); it1.hasNext();) {
                GKInstance topic = (GKInstance) it1.next();
                // Make sure this is a human pathway
                GKInstance species = (GKInstance) topic.getAttributeValue(ReactomeJavaConstants.species);
                if (!species.getDBID().equals(48887L))
                    continue;
                Set<String> ids = grepIDsFromTopic(topic);
                if (ids.size() > 200) // This is arbitrary
                    bigTopics.add(topic);
                else
                    pathwaySet.add(topic);
            }
        }
        Set<GKInstance> next = new HashSet<GKInstance>();
        for (int i = 0; i < 2; i++) { // Run two levels only
            // Have to split the big topics
            for (GKInstance topic : bigTopics) {
                List comps = null;
                if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasSpecialisedForm))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasSpecialisedForm);
                else if (topic.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasMember))
                    comps = topic.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                if (comps == null || comps.size() == 0) {
                    pathwaySet.add(topic);
                    continue;
                }
                // If there is any reaction in the pathway, it should not split any more
                boolean isAdded = false;
                for (Object obj : comps) {
                    GKInstance subEvent = (GKInstance) obj;
                    if (subEvent.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
                        pathwaySet.add(topic);
                        isAdded = true;
                        break;
                    }
                }
                if (isAdded)
                    continue;
                for (Iterator<?> it = comps.iterator(); it.hasNext();) {
                    GKInstance sub = (GKInstance) it.next();
                    if (i == 1)
                        pathwaySet.add(sub);
                    else {
                        //if (sub.getDBID().equals(163359L)) 
                        //    continue; // Escape Glucagon signaling in metabolic regulation(163359)
                        //              // This pathway has been included by Integration of pathways involved in energy metabolism(163685)
                        // Check sub-pathway size
                        Set<String> ids = grepIDsFromTopic(sub);
                        if (ids.size() > 300) 
                            next.add(sub);
                        else
                            pathwaySet.add(sub);
                    }
                }
            }
            bigTopics.clear();
            bigTopics.addAll(next);
        }
        // Want to sort it before output
        List<GKInstance> list = new ArrayList<GKInstance>(pathwaySet);
        InstanceUtilities.sortInstances(list);
        String fileName = FIConfiguration.getConfiguration().get("REACTOME_PATHWAYS");
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        for (GKInstance pathway : list)
            fu.printLine(pathway.getDBID() + "\t" + pathway.getDisplayName());
        fu.close();
        System.out.println("Total Pathways: " + list.size());
    }
}
