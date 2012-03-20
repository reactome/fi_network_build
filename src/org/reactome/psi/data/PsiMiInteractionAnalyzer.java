/*
 * Created on Apr 27, 2006
 *
 */
package org.reactome.psi.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.HibernateUtil;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.psi.DbReference;
import org.reactome.psi.Experiment;
import org.reactome.psi.Interaction;
import org.reactome.psi.Interactor;
import org.reactome.psi.OpenCV;
import org.reactome.psi.Participant;
import org.reactome.psi.Xref;

/**
 * This class is used to analyze a PSI-MI database.
 * @author guanming
 *
 */
public class PsiMiInteractionAnalyzer extends TestCase {
    private final String RESULT_DIR = "results/interaction/";
    private final String INTERACTION_XML_FILE = RESULT_DIR + "HPRDInteractions.xml";
    
    private SessionFactory sessionFactory;

    public PsiMiInteractionAnalyzer() {
        String configFileName = "resources/hibernate.cfg.xml";
        File configFile = new File(configFileName);
        sessionFactory = HibernateUtil.getSessionFactory(configFile);
    }
    
    public PsiMiInteractionAnalyzer(SessionFactory sf) {
        this.sessionFactory = sf;
    }
    
    public void setSessionFactory(SessionFactory sf) {
        this.sessionFactory = sf;
    }
    
    public void extractInteractions(String fileName) throws Exception {
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from " + Interaction.class.getName());
        int start = 0;
        int max = 2000;
        // Use pagination query to minimize the memory usage.
        query.setFirstResult(start);
        query.setMaxResults(max);
        List list = query.list();
        Interaction interaction = null;
        Map<String, NTInteraction> interactions = new HashMap<String, NTInteraction>();
        while (list.size() > 0) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                interaction = (Interaction) it.next();
                List<Participant> participants = interaction.getParticipantList();
                String key = generateInteractionKey(participants);
                if (key == null)
                    continue;
                NTInteraction ntInteraction = interactions.get(key);
                if (ntInteraction != null) {
                    // Check if there are any new experiments
                    // Check experiments
                    List<Experiment> experiments = interaction.getExperimentList();
                    if (experiments != null) {
                        for (Experiment exp : experiments) {
                            OpenCV openCV = exp.getInteractionDetectionMethod();
                            if (openCV != null) {
                                ntInteraction.addExperimentType(openCV.getNames().getShortLabel());
                            }
                        }
                    }
                    // Check if there are any new interactionTypes
                    OpenCV interactionType = interaction.getInteractionType();
                    if (interactionType != null)
                        ntInteraction.addInteractionType(interactionType.getNames().getShortLabel());
                }
                else {
                    // Create a new interaction
                    ntInteraction = distill(interaction);
                    interactions.put(key, ntInteraction);
                }
            }
            session.clear(); // To empty cache to keep the memory usage small.
            start += list.size();
            query.setFirstResult(start);
            list = query.list();
            System.out.println("Current interactions: " + interactions.size());
        }
        session.close();
        exportToXML(interactions, fileName);
    }
    
    private void exportToXML(Map<String, NTInteraction> interactions, String fileName) 
                 throws IOException {
        StringBuilder builder = new StringBuilder();
        builder.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        builder.append("<interactions>");
        Set<String> keys = interactions.keySet();
        for (String key : keys) {
            NTInteraction interaction = interactions.get(key);
            builder.append("\n");
            interaction.exportToXML(builder);
        }
        builder.append("\n</interactions>");
        FileWriter fileWriter = new FileWriter(fileName);
        BufferedWriter bw = new BufferedWriter(fileWriter);
        bw.write(builder.toString());
        bw.close();
        fileWriter.close();
    }
    
    /**
     * When this method is called, an non-empty key should be created already, which
     * means that there are at least two UniProtKB ids existing.
     * @param interaction
     * @return
     */
    private NTInteraction distill(Interaction interaction) {
        NTInteraction ntInteraction = new NTInteraction();
        ntInteraction.setDbId(interaction.getDbId());
        List<Participant> participants = interaction.getParticipantList();
        for (Participant p : participants) {
            Interactor interactor = p.getInteractor();
            NTInteractor ntInteractor = distill(interactor);
            if (ntInteractor != null)
                ntInteraction.addInteractor(ntInteractor);
        }
        // Check experiments
        List<Experiment> experiments = interaction.getExperimentList();
        if (experiments != null) {
            for (Experiment exp : experiments) {
                OpenCV openCV = exp.getInteractionDetectionMethod();
                if (openCV != null) {
                    ntInteraction.addExperimentType(openCV.getNames().getShortLabel());
                }
            }
        }
        // Check interactionType
        OpenCV interactionType = interaction.getInteractionType();
        if (interactionType != null)
            ntInteraction.addInteractionType(interactionType.getNames().getShortLabel());
        return ntInteraction;
    }
    
    private NTInteractor distill(Interactor interactor) {
        DbReference uniXref = getUniProtKBRef(interactor.getXref());
        if (uniXref == null)
            return null; // Don't need to create.
        NTInteractor ntInteractor = new NTInteractor();
        ntInteractor.setDbId(interactor.getDbId());
        ntInteractor.setId(uniXref.getId());
        ntInteractor.setName(interactor.getNames().getShortLabel());
        // Extract type
        OpenCV opencv = interactor.getInteractorType();
        if (opencv != null)
            ntInteractor.setType(NTInteractorType.valueOf(opencv.getNames().getShortLabel()));
        return ntInteractor;
    }
    
    private String generateInteractionKey(List<Participant> participants) {
        Set<DbReference> uniXrefList = new HashSet<DbReference>();
        for (Participant p : participants) {
            Interactor interactor = p.getInteractor();
            // Get xref
            Xref xref = interactor.getXref();
            if (xref == null) {
                // Most of them are small molecules
                //System.out.println("Interactor has not xref: " + 
                  //      interactor.getDbId() + " " + 
                    //    interactor.getNames().getShortLabel());
                continue;
            }
            DbReference uniRef = getUniProtKBRef(xref);
            if (uniRef != null) {  
                uniXrefList.add(uniRef);
            }
        }
        List<String> ids = new ArrayList<String>();
        for (DbReference ref : uniXrefList) {
            ids.add(ref.getId());
        }
        if (ids.size() < 2)
            return null; // Might be possible
        Collections.sort(ids);
        StringBuilder builder = new StringBuilder();
        for (Iterator it = ids.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append("->");
        }
        return builder.toString();
    }
    
    private DbReference getUniProtKBRef(Xref xref) {
        DbReference dbRef = xref.getPrimaryRef();
        if (dbRef.getDb().startsWith("uniprot"))
            return dbRef;
        List<DbReference> dbRefList = xref.getSecondaryRefList();
        if (dbRefList == null || dbRefList.size() == 0)
            return null;
        for (DbReference ref : dbRefList) {
            if (ref.getDb().startsWith("uniprot"))
                return ref;
        }
        return null;
    }
    
    public void analyzeExperimentNumbers() throws Exception {
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from " + Interaction.class.getName());
        int start = 0;
        int max = 2000;
        // Use pagination query to minimize the memory usage.
        query.setFirstResult(start);
        query.setMaxResults(max);
        List list = query.list();
        Interaction interaction = null;
        Map<String, Integer> expNumbers = new HashMap<String, Integer>();
        while (list.size() > 0) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                interaction = (Interaction) it.next();
                List<Participant> participants = interaction.getParticipantList();
                String key = generateInteractionKey(participants);
                if (key == null)
                    continue;
                List<Experiment> experiments = interaction.getExperimentList();
                if (experiments != null) {
                    Integer tmp = expNumbers.get(key);
                    if (tmp == null)
                        expNumbers.put(key, experiments.size());
                    else
                        expNumbers.put(key, tmp + experiments.size());
                }
            }
            session.clear(); // To empty cache to keep the memory usage small.
            start += list.size();
            query.setFirstResult(start);
            list = query.list();
        }
        session.close();
        analyzeExpNumberMap(expNumbers);
    }
    
    private void analyzeExpNumberMap(Map<String, Integer> map) {
        Map<Integer, Integer> numberMap = new HashMap<Integer, Integer>();
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String key = it.next();
            Integer c = map.get(key);
            Integer tmp = numberMap.get(c);
            if (tmp == null)
                numberMap.put(c, 1);
            else
                numberMap.put(c, tmp + 1);
        }
        for (Iterator<Integer> it = numberMap.keySet().iterator(); it.hasNext();) {
            Integer key = it.next();
            Integer value = numberMap.get(key);
            System.out.println(key + " " + value);
        }
    }
    
    public void analyzeExperimentTypes() throws Exception {
        String interactionFileName = INTERACTION_XML_FILE;
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        System.out.println("Total interactions: " + intElms.size());
        Set<String> totalExpNames = new HashSet<String>();
        int c = 0;
        Set<String> expNamesSet = new HashSet<String>();
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            Element expElm = elm.getChild("experimentTypes");
            String[] expNames = expElm.getText().split(",");
            for (String tmp : expNames)
                expNamesSet.add(tmp);
            if (expNames.length > 1)
                c ++;
        }
        System.out.println("Interactions with more than two experiments: " + c);
        System.out.println("Experiments: " + expNamesSet);
    }
    
    public void extractIDsFromInteractions() throws Exception {
        String interactionFileName = "results/interactions.xml";
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        System.out.println("total interactions: " + intElms.size());
        Set<String> ids = new HashSet<String>();
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            Element expElm = elm.getChild("experimentTypes");
            String[] expNames = expElm.getText().split(",");
            Element actorsElm = elm.getChild("interactors");
            List actorsList = actorsElm.getChildren("interactor");
            Element tmp = (Element) actorsList.get(0);
            String id1 = tmp.getAttributeValue("id");
            tmp = (Element) actorsList.get(1);
            String id2 = tmp.getAttributeValue("id");
            ids.add(id1);
            ids.add(id2);
        }
        System.out.println("Total IDs: " + ids.size());
    }  
    
    public void tally() throws Exception {
        String[] fileNames = new String[] {
            "HumanInteractions.txt",
            "OrthoInteractions.txt",
            "YeastInteractions.txt"
        };
        FileUtility fu = new FileUtility();
        Set<String> interactions = null;
        Set<String> ids = null;
        for (String fileName : fileNames) {
            System.out.println("FileName: " + fileName);
            interactions = fu.loadInteractions(RESULT_DIR + fileName);
            System.out.println("Total Interactions: " + interactions.size());
            ids = InteractionUtilities.grepIDsFromInteractions(interactions);
            System.out.println("Total ids: " + ids.size());
        }
    }
    
    public void loadHumanInteractionsIDs() throws Exception {
        Set<String> y2hIds = new HashSet<String>();
        Set<String> nonY2hIds = new HashSet<String>();
        String interactionFileName = "results/interaction/interactions.xml";
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        System.out.println("total interactions: " + intElms.size());
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            Element expElm = elm.getChild("experimentTypes");
            String[] expNames = expElm.getText().split(",");
            Element actorsElm = elm.getChild("interactors");
            List actorsList = actorsElm.getChildren("interactor");
            Element tmp = (Element) actorsList.get(0);
            String id1 = tmp.getAttributeValue("id");
            tmp = (Element) actorsList.get(1);
            String id2 = tmp.getAttributeValue("id");
            for (String expName : expNames) {
                if (expName.startsWith("two hybrid")) {
                    y2hIds.add(id1);
                    y2hIds.add(id2);
                }
                else {
                    nonY2hIds.add(id1);
                    nonY2hIds.add(id2);
                }
            }
        }
        System.out.println("Total Y2H: " + y2hIds.size());
        System.out.println("Total Non Y2H: " + nonY2hIds.size());
        nonY2hIds.addAll(y2hIds);
        System.out.println("Total IDs: " + nonY2hIds.size());
    }
}