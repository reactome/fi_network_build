/*
 * Created on Jun 26, 2006
 *
 */
package org.reactome.psi.data;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.events.StartElement;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.InteractionUtilities;
import org.reactome.psi.PsiMiModel;

/**
 * Focus on Human PsiMi Interactions only.
 * @author guanming
 */
public class HumanPsiMiInteractionAnalyzer extends PsiMiInteractionExtractor {
    private final String DATASET_DIR = "/Users/wgm/Documents/caBIG_R3/datasets/";
    private final String INTERACTION_XML_FILE = "results/interaction/HPRDInteractions.xml";
    // Used for human only
    private Map<String, Set<String>> id2ProteinGis;
    
    public void generateNonInteractionData() throws Exception {
        Set<String> nonY2H = loadNonY2HInteraction();
        Set<String> ids = InteractionUtilities.grepIDsFromInteractions(nonY2H);
        FileUtility fu = new FileUtility();
        fu.setOutput("results/interaction/HumanNonY2HNoInteractions.txt");
        List<String> idList = new ArrayList<String>(ids);
        String id1, id2;
        String pair = null;
        for (int i = 0; i < idList.size() - 1; i++) {
            id1 = idList.get(i);
            for (int j = i + 1; j < idList.size(); j++) {
                id2 = idList.get(j);
                pair = createProteinPair(id1, id2);
                if (!nonY2H.contains(pair))
                    fu.printLine(pair);
            }
        }
        fu.close();
    }
    
    public void analyzeOverlapping() throws Exception {
        Set<String> y2hInteractions = new HashSet<String>();
        Set<String> nonY2hInteractions = new HashSet<String>();
        loadInteractions(y2hInteractions, nonY2hInteractions);
        // Output
//        FileUtility fu = new FileUtility();
//        fu.saveInteractions(y2hInteractions, "results/interaction/HPRDHumanY2H.txt");
//        fu.saveInteractions(nonY2hInteractions, "results/interaction/HPRDHumanNonY2H.txt");
//        Set<String> merged = new HashSet<String>();
//        merged.addAll(y2hInteractions);
//        merged.addAll(nonY2hInteractions);
//        fu.saveInteractions(merged, "results/interaction/HPRDHumanInteractions.txt");
        y2hInteractions.retainAll(nonY2hInteractions);
        System.out.println("Overlapping interactions: " + y2hInteractions.size());
        Set<String> y2hIds = InteractionUtilities.grepIDsFromInteractions(y2hInteractions);
        Set<String> nonY2hIds = InteractionUtilities.grepIDsFromInteractions(nonY2hInteractions);
        System.out.println("IDs in y2h: " + y2hIds.size());
        System.out.println("IDs in nonY2h: " + nonY2hIds.size());
        nonY2hIds.retainAll(y2hIds);
        System.out.println("Overlapping Ids: " + nonY2hIds.size());
    }
    
    public void analyzeExperiments() throws Exception {
        String interactionFileName = "results/interaction/interactions.xml";
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        System.out.println("total interactions: " + intElms.size());
        // Check for experiments
        Map<Integer, Integer> expNumbers = new HashMap<Integer, Integer>();
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            Element expElm = elm.getChild("experimentTypes");
            String[] expNames = expElm.getText().split(",");
            Integer intNumber = expNumbers.get(expNames.length);
            if (intNumber == null)
                expNumbers.put(expNames.length, 1);
            else
                expNumbers.put(expNames.length, intNumber + 1);
        }
        for (Iterator<Integer> it = expNumbers.keySet().iterator(); it.hasNext();) {
            Integer expNumber = it.next();
            Integer intNumber = expNumbers.get(expNumber);
            System.out.println(expNumber + " " + intNumber);
        }
    }
    
    private void loadInteractions(Set<String> y2hInteractions, Set<String> nonY2HInteractions)
        throws Exception {
        String interactionFileName = "results/interaction/HPRDInteractions.xml";
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
            String prtPair = createProteinPair(id1, id2);
            for (String expName : expNames) {
                if (expName.startsWith("two hybrid") ||
                    expName.startsWith("2H"))
                    y2hInteractions.add(prtPair);
                else
                    nonY2HInteractions.add(prtPair);
            }
        }
        System.out.println("Total Y2H: " + y2hInteractions.size());
        System.out.println("Total Non Y2H: " + nonY2HInteractions.size());
    }
    
    private String createProteinPair(String id1, String id2) {
        int compare = id1.compareTo(id2);
        if (compare < 0)
            return id1 + " " + id2;
        return id2 + " " + id1;
    }
    
    public Set<String> loadHumanInteractions() throws Exception {
        String interactionFileName = INTERACTION_XML_FILE;
        return loadHumanInteractions(interactionFileName);
    }

    public Set<String> loadHumanInteractions(String interactionFileName) throws JDOMException, IOException {
        Set<String> interactions = new HashSet<String>();
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        int compare = 0;
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            //Element expElm = elm.getChild("experimentTypes");
            //String[] expNames = expElm.getText().split(",");
            Element actorsElm = elm.getChild("interactors");
            List actorsList = actorsElm.getChildren("interactor");
            Element tmp = (Element) actorsList.get(0);
            String id1 = tmp.getAttributeValue("id");
            tmp = (Element) actorsList.get(1);
            String id2 = tmp.getAttributeValue("id");
            compare = id1.compareTo(id2);
            if (compare < 0)
                interactions.add(id1 + " " + id2);
            else if (compare > 0)
                interactions.add(id2 + " " + id1);
            // Exclude self interactions
        }
        return interactions;
    }
    
    public Set<String> loadNonY2HInteraction() throws Exception {
        Set<String> interactions = new HashSet<String>();
        String interactionFileName = "results/interaction/interactions.xml";
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build(interactionFileName);
        Element root = document.getRootElement();
        List intElms = root.getChildren("interaction");
        int compare = 0;
        boolean nonY2H = false;
        for (Iterator it = intElms.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            Element expElm = elm.getChild("experimentTypes");
            String[] expNames = expElm.getText().split(",");
            nonY2H = false;
            for (String expName : expNames) {
                if (!expName.startsWith("two hybrid")) {
                    nonY2H = true;
                }
            }
            if (!nonY2H)
                continue;
            Element actorsElm = elm.getChild("interactors");
            List actorsList = actorsElm.getChildren("interactor");
            Element tmp = (Element) actorsList.get(0);
            String id1 = tmp.getAttributeValue("id");
            tmp = (Element) actorsList.get(1);
            String id2 = tmp.getAttributeValue("id");
            compare = id1.compareTo(id2);
            if (compare < 0)
                interactions.add(id1 + " " + id2);
            else if (compare > 0)
                interactions.add(id2 + " " + id1);
        }
        return interactions;
    }
    
    public void parse() throws Exception {
        long time1 = System.currentTimeMillis();
        String primaryDBName = "PROTEIN GI";
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/BioGrid/" +
        "GRID-ORGANISM-Homo_sapiens-2.0.19.psi25.xml";
        String y2hOutFileName = "results/interaction/humanY2HInteraction.txt";
        String nonY2HOutFileName = "results/interaction/humanNonY2HInteraction.txt";
        id2ProteinGis = new HashMap<String, Set<String>>();
        parse(primaryDBName, fileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for parsing: " + (time2 - time1));
        System.out.println("Total interactors: " + id2ProteinGis.size());
        outputInteractions(interactions, y2hOutFileName, nonY2HOutFileName);
    }
    
    public Map<String, Set<String>> generateGI2UniMap() throws IOException {
        String mapFileName = DATASET_DIR + "iproclass/iproclass.tb";
        FileUtility fu = new FileUtility();
        fu.setInput(mapFileName);
        String line = null;
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        String[] tokens = null;
        while ((line = fu.readLine()) != null) {
            tokens = line.split("\t");
            String acc = tokens[0];
            String id = tokens[1];
            if (!id.endsWith("_HUMAN"))
                continue;
            String gis = tokens[4];
            if (gis.length() == 0)
                continue;
            tokens = gis.split("; ");
            for (String gi : tokens) {
                Set<String> set = map.get(gi);
                if (set == null) {
                    set = new HashSet<String>();
                    map.put(gi, set);
                }
                set.add(acc);
            }
        }
        fu.close();
        System.out.println("Total GI to UniProt Map: " + map.size());
        return map;
    }
    
    public void checkGI2UniMap() throws IOException {
        Map<String, Set<String>> map = generateGI2UniMap();
        System.out.println("Total GI Protein: " + map.size());
        // Check how many GIs have more than one UniProts
        int c = 0;
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String gi = it.next();
            Set<String> uniIds = map.get(gi);
            if (uniIds.size() > 1)
                c ++;
        }
        System.out.println("Number of GI having more than one UniProts: " + c);
        // Save it
        String fileName = "results/GI2UniProt.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(fileName);
        for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
            String gi = it.next();
            Set<String> uniIds = map.get(gi);
            for (String uni : uniIds) {
                fu.printLine(gi + "\t" + uni);
            }
        }
        fu.close();
    }
    
    public Map<String, Set<String>> loadGI2UniMap() throws IOException {
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        String fileName = "results/GI2UniProt.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        int index = 0;
        String giId = null;
        String uniId = null;
        Set<String> uniIdSet = null;
        while ((line = fu.readLine()) != null) {
            index = line.indexOf("\t");
            giId = line.substring(0, index);
            uniId = line.substring(index + 1);
            uniIdSet = map.get(giId);
            if (uniIdSet == null) {
                uniIdSet = new HashSet<String>();
                map.put(giId, uniIdSet);
            }
            uniIdSet.add(uniId);
        }
        return map;
    }
    
    protected void outputInteractions(List<Interaction> interactions, 
            String y2hOutFileName,
            String nonY2HOutFileName) throws IOException {
        // Need to replace GI with UniProtein Accs
        Map<String, Set<String>> gi2uni = generateGI2UniMap();
        // Check participants list
        boolean isNotEmpty = true;
        int c = 0;
        HumanInteraction hi;
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            hi = (HumanInteraction) it.next();
            if (hi.participantList == null)
                continue;
            List<Set<String>> uniPartiList = new ArrayList<Set<String>>();
            for (Set<String> set : hi.participantList) {
                Set<String> uniSet = new HashSet<String>();
                for (String gi : set) {
                    Set<String> uni = gi2uni.get(gi);
                    if (uni == null)
                        continue;
                    uniSet.addAll(uni);
                }
                if (uniSet.size() == 0) {
                    isNotEmpty = false;
                    break;
                }
                uniPartiList.add(uniSet);
            }
            if (isNotEmpty)
                hi.participantList = uniPartiList;
            else {
                it.remove();
                c ++;
            }
        }
        System.out.println("Cannot mapped to UniProtein interactions: " + c);
        // Need to make sure participantList is not empty
        FileUtility y2hFu = new FileUtility();
        y2hFu.setOutput(y2hOutFileName);
        FileUtility nony2hFu = new FileUtility();
        nony2hFu.setOutput(nonY2HOutFileName);
        Set<String> nony2hSet = new HashSet<String>();
        Set<String> y2hSet = new HashSet<String>();
        for (Interaction i : interactions) {
            hi = (HumanInteraction) i;
            Set<String> partnerList1 = hi.participantList.get(0);
            Set<String> partnerList2 = hi.participantList.get(1);
            for (String partner1 : partnerList1) {
                for (String partner2 : partnerList2) {
                    int compare = partner1.compareTo(partner2);
                    String pair = null;
                    if (compare < 0)
                        pair = partner1 + " " + partner2;
                    else
                        pair = partner2 + " " + partner1;
                    if (i.expNames == null)
                        nony2hSet.add(pair);
                    else if (i.expNames.size() == 1) {
                        String expName = i.expNames.get(0);
                        if (expName.equals("Two-hybrid"))
                            y2hSet.add(pair);
                        else
                            nony2hSet.add(pair);
                    }
                    else {
                        if (i.expNames.contains("Two-hybrid"))
                            y2hSet.add(pair);
                        nony2hSet.add(pair);
                    }
                }
            }
        }
        c = 0;
        for (String pair : y2hSet) {
            y2hFu.printLine(pair);
            c ++;
        }
        System.out.println("Y2H: " + c);
        c = 0;
        for (String pair : nony2hSet) {
            nony2hFu.printLine(pair);
            c ++;
        }
        System.out.println("NonY2H: " + c);
        y2hFu.close();
        nony2hFu.close();
    }
    
    protected void analyzeInteractions(List<Interaction> interactions) {
        int c = 0;
        for (Interaction interaction : interactions) {
            c ++;
            System.out.println(c + ": " + interaction);
            if (c > 20)
                break;
        }
        int c0 = 0;
        int c3 = 0;
        for (Iterator<Interaction> it = interactions.iterator(); it.hasNext();) {
            HumanInteraction tmp = (HumanInteraction) it.next();
            if (tmp.participantList == null ||
                    tmp.participantList.size() < 2) {
                it.remove();
                c0 ++;
            }
            else if (tmp.participantList.size() > 2)
                c3 ++;
        }
        System.out.println("Interaction after cleanup: " + interactions.size());
        System.out.println("Interaction containing no interactors: " + c0);
        System.out.println("Interaction containing more than 2 interactors: " + c3);
    }
    
    protected void handleStartElement(StartElement startElement, 
            String primaryDBName) throws Exception {
        String elmName = startElement.getName().getLocalPart();
        if (elmName.equals(PsiMiModel.experimentList)) {
            isInExpList = true;
        }
        else if (elmName.equals(PsiMiModel.interactorList)) {
            isInInteractorList = true;
        }
        else if (elmName.equals(PsiMiModel.interactionList)) {
            isInInteractionList = true;
        }
        else if (isInInteractionList && elmName.equals(PsiMiModel.interaction)) {
            isInInteraction = true;
            interaction = new HumanInteraction();
            String id = getAttributeValue(startElement, PsiMiModel.id);
            interaction.id = id;
        }
        else if (elmName.equals(PsiMiModel.interactionDetectionMethod)) {
            isInCV = true;
        }
        else if (isInInteractorList && elmName.equals(PsiMiModel.interactor)) {
            isInInteractor = true;
            id = getAttributeValue(startElement, PsiMiModel.id);
        }
        else if (isInInteractor && elmName.equals(PsiMiModel.xref))
            isInxref = true;
        else if (isInInteraction && isInExpList && elmName.equals(PsiMiModel.experimentRef)) {
            String ref = getAttributeValue(startElement, "ref");
            String name = expId2Name.get(ref);
            if (name == null) {
                throw new IllegalStateException("Experiment cannot be found for " +
                        "interaction: " + ref + " in " + interaction.id);
            }
            interaction.addExp(name);
        }
        else if (elmName.equals(PsiMiModel.participantList)) {
            isInParticipantList = true;
        }
        else if (isInParticipantList && elmName.equals(PsiMiModel.interactorRef)) {
            String ref = getAttributeValue(startElement, "ref");
            Set<String> accession = id2ProteinGis.get(ref);
            if (accession != null)
                ((HumanInteraction)interaction).addParticipant(accession);
        }
        else if (elmName.equals(PsiMiModel.interactionType)) {
            isInInteractionType = true;
            isInCV = true;
        }
        else if (isInInteractor && isInxref && elmName.equals(PsiMiModel.secondaryRef)) {
            String db = getAttributeValue(startElement, PsiMiModel.db);
            String acc = getAttributeValue(startElement, PsiMiModel.id);
            if (db.equals(primaryDBName)) {
                Set<String> gis = id2ProteinGis.get(id);
                if (gis == null) {
                    gis = new HashSet<String>();
                    id2ProteinGis.put(id, gis);
                }
                gis.add(acc);
            }
        }
//      Use xref for dme and sce
        else if (isInExpList && elmName.equals(PsiMiModel.experimentDescription)) {
            id = getAttributeValue(startElement, PsiMiModel.id);
        }
        else if (isInCV && elmName.equals(PsiMiModel.shortLabel)) {
            needText = true;
        }
    }
    
    public void generateHumanInteractions() throws Exception {
        String[] fileNames = new String[] {
                "results/v2/BINDInteractions020507.txt",
                "results/v2/IntActInteractions020507.txt",
                "results/v2/HPRDInteractions020507.txt",
                //"results/interaction/OrthoInteractions.txt",
                //"results/interaction/YeastInteractions.txt"
        };
        FileUtility fu = new FileUtility();
        Set<String> interactions = new HashSet<String>();
        for (String fileName : fileNames) {
            Set<String> set = fu.loadInteractions(fileName);
            interactions.addAll(set);
        }
        System.out.println("Total: " + interactions.size());
        //fu.saveInteractions(interactions, 
        //                    "results/interaction/HumanInteractions090806.txt");
        // HumanInteractions120506.txt: merged three interaction sets, BINDInteractions090806.txt,
        // IntActInteractions090806.txt, and HPRDInteractions120506.txt
        //fu.saveInteractions(interactions, 
        //                    "results/interaction/HumanInteractions120506.txt");
        //fu.saveInteractions(interactions, 
        //                    "results/interaction/ThreePPIs091406.txt");    
        fu.saveInteractions(interactions,
                            "results/v2/HumanInteractions020507.txt");
    }
    
    public void cleanHPRDInteractions() throws IOException {
        String fileName = "results/interaction/HPRDInteractions090806.txt";
        FileUtility fu = new FileUtility();
        Map<String, String> hprdToUniprotViaGene = fu.importMap("results/HPRD2UniProtViaGeneSymbols.txt");
        Set<String> interactions = fu.loadInteractions(fileName);
        Set<String> newInteractions = new HashSet<String>();
        int index = 0;
        String id1, id2;
        int order;
        int c = 0;
        for (String i : interactions) {
            index = i.indexOf(" ");
            id1 = i.substring(0, index);
            id2 = i.substring(index + 1);
            String mapped1 = hprdToUniprotViaGene.get(id1);
            String mapped2 = hprdToUniprotViaGene.get(id2);
            if (mapped1 != null || mapped2 != null) {
                if (mapped1 != null)
                    id1 = mapped1;
                if (mapped2 != null)
                    id2 = mapped2;
                order = id1.compareTo(id2);
                if (order < 0)
                    newInteractions.add(id1 + " " + id2);
                else
                    newInteractions.add(id2 + " " + id1);
                c ++;
            }
            else
                newInteractions.add(i);
        }
        System.out.println("New Mapped HPRD Interactions: " + c);
        fu.saveInteractions(newInteractions, "results/interaction/HPRDInteractions120506.txt");
    }
    
    private class HumanInteraction extends Interaction {
        
        private List<Set<String>> participantList;
        
        public void addParticipant(Set<String> ids) {
            if (participantList == null)
                participantList = new ArrayList<Set<String>>();
            participantList.add(ids);
        }
        
    }
}
