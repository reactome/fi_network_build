/*
 * Created on Jun 19, 2006
 *
 */
package org.reactome.psi.data;

import static javax.xml.stream.XMLStreamConstants.CHARACTERS;
import static javax.xml.stream.XMLStreamConstants.END_DOCUMENT;
import static javax.xml.stream.XMLStreamConstants.END_ELEMENT;
import static javax.xml.stream.XMLStreamConstants.START_DOCUMENT;
import static javax.xml.stream.XMLStreamConstants.START_ELEMENT;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.stream.XMLEventReader;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.events.Attribute;
import javax.xml.stream.events.Characters;
import javax.xml.stream.events.EndElement;
import javax.xml.stream.events.StartElement;
import javax.xml.stream.events.XMLEvent;

import junit.framework.TestCase;

import org.reactome.fi.util.ChecksumCalculator;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.psi.PsiMiModel;

/**
 * This class is used to parse interaction information from PSI-MI files directly.
 * 
 * @author guanming
 *
 */
public class PsiMiInteractionExtractor extends TestCase {
    protected Map<String, String> expId2Name;
    protected Map<String, String> interactorId2Accession;
    protected Map<String, String> interactorId2Sequence;
    protected List<Interaction> interactions;
    // These flags are used to indicate the iterator positions
    protected boolean isInExpList;
    protected boolean isInExp;
    protected boolean isInInteractorList;
    protected boolean isInInteractorRef;
    protected boolean isInInteractionList;
    protected boolean isInCV;
    protected boolean needText;
    protected boolean isInInteractor;
    protected boolean isInxref;
    protected boolean isInParticipantList;
    protected boolean isInInteractionType;
    protected boolean isInInteraction;
    protected boolean isInSequence;
    // Temp values
    protected String id;
    protected Interaction interaction;
    // Use to control some different values
    private boolean useSecondaryRef;
    // To calculate sequence check sum
    private ChecksumCalculator checksum;

    
    public PsiMiInteractionExtractor() {
        expId2Name = new HashMap<String, String>();
        interactorId2Accession = new HashMap<String, String>();
        interactions = new ArrayList<Interaction>();
        interactorId2Sequence = new HashMap<String, String>();
        checksum = new ChecksumCalculator();
    }
    
    private void reset() {
        expId2Name.clear();
        interactorId2Accession.clear();
        interactions.clear();
        interactorId2Sequence.clear();
        id = null;
        interaction = null;
    }
    
    /**
     * This method can be used to process a list of PSI-MI file
     * @throws Exception
     */
    public void batchParse() throws Exception {
        // For Fly
        String speciesName = "caeel";
        long time1 = System.currentTimeMillis();
        String dirName = FIConfiguration.getConfiguration().get("INTACT_DIR");
        String primaryDBName = "uniprotkb";
        useSecondaryRef = false; 
        List<String> fileNames = fetchFileName(speciesName, 
                                               dirName);
        Set<String> interactionsInCheckSum = new HashSet<String>();
        for (String fileName : fileNames) {
            parse(primaryDBName, fileName);
            extractInteractionsInChecksum(interactionsInCheckSum);
            reset();
            //break;
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for parsing: " + (time2 - time1));
        // For output
        FileUtility fu = new FileUtility();
        fu.outputSet(interactionsInCheckSum, dirName + speciesName + "_interactions_in_checksum.txt");
    }
    
    public void batchParsePPIInUniProts() throws Exception {
        String[] speciesNames = new String[] {
                "caeel", // C. elegans
                "drome", // D. meanologa
                "yeast" // S. cerevias
        };
        for (String speciesName : speciesNames) {
            long time1 = System.currentTimeMillis();
            String dirName = FIConfiguration.getConfiguration().get("INTACT_DIR");
            String primaryDBName = "uniprotkb";
            useSecondaryRef = false; 
            List<String> fileNames = fetchFileName(speciesName, 
                                                   dirName);
            Set<String> ppiInUniProt = new HashSet<String>();
            for (String fileName : fileNames) {
                parse(primaryDBName, fileName);
                extractPPIsInUniProt(ppiInUniProt);
                reset();
                //break;
            }
            long time2 = System.currentTimeMillis();
            System.out.println("Total time for parsing: " + (time2 - time1));
            // For output
            FileUtility fu = new FileUtility();
            fu.outputSet(ppiInUniProt, 
                         dirName + speciesName + "_interactions_in_uniprot.txt");
        }
    }
    
    private void extractPPIsInUniProt(Set<String> ppisInUniProt) {
        for (Interaction interaction : interactions) {
            List<String> participants = interaction.participants;
            if (participants == null || participants.size() < 2)
                continue; // Self-interction
            // Convert interactions to uniprot
            for (int i = 0; i < participants.size() - 1; i++) {
                String id1 = participants.get(i);
                for (int j = 1; j < participants.size(); j++) {
                    String id2 = participants.get(j);
                    int compare = id1.compareTo(id2);
                    if (compare < 0) 
                        ppisInUniProt.add(id1 + " " + id2);
                    else if (compare > 0)
                        ppisInUniProt.add(id2 + " " + id1);
                }
            }
        }
    }
    
    private void extractInteractionsInChecksum(Set<String> interactionsInCheckSum) {
        for (Interaction interaction : interactions) {
            List<String> sequences = interaction.sequences;
            if (sequences.size() < 2)
                continue; // Self-interction
            // Convert interactions to checksum
            List<String> checksums = new ArrayList<String>();
            for (String sequence : sequences) {
                checksums.add(checksum.calculateChecksum(sequence));
            }
            for (int i = 0; i < checksums.size() - 1; i++) {
                String checksum1 = checksums.get(i);
                for (int j = 1; j < checksums.size(); j++) {
                    String checksum2 = checksums.get(j);
                    int compare = checksum1.compareTo(checksum2);
                    if (compare < 0) 
                        interactionsInCheckSum.add(checksum1 + " " + checksum2);
                    else if (compare > 0)
                        interactionsInCheckSum.add(checksum2 + " " + checksum1);
                }
            }
        }
    }
    
    private List<String> fetchFileName(String startWith,
                                       String dirName) {
        List<String> fileNames = new ArrayList<String>();
        File dir = new File(dirName);
        File[] files = dir.listFiles();
        for (File file : files) {
            String fileName = file.getName();
            if (fileName.startsWith(startWith) &&
                fileName.endsWith(".xml") &&
                !fileName.endsWith("_negative.xml")) {
                fileNames.add(file.getAbsolutePath());
            }
        }
        return fileNames;
    }
    
    /**
     * Use this method to generate PPIs for other species from BioGrid data sets.
     * @throws Exception
     */
    public void parse() throws Exception {
        long time1 = System.currentTimeMillis();
        String biogridDir = FIConfiguration.getConfiguration().get("BIOGRID_DIR");
        String outputDir = biogridDir;
//      String primaryDBName = "SGD";
//      String fileName = biogridDir + "BIOGRID-ORGANISM-Saccharomyces_cerevisiae-2.0.50.psi25.xml";
//      String y2hOutFileName = outputDir + "yeastY2HInteraction2.0.50.txt";
//      String nonY2HOutFileName = outputDir + "yeastNonY2HInteraction2.0.50.txt";
//      boolean isForCel = false;
//      useSecondaryRef = false; 
//        String primaryDBName = "WormBase";
//        String fileName = biogridDir +
//                          "BIOGRID-ORGANISM-Caenorhabditis_elegans-2.0.32.psi25.xml";
//        String y2hOutFileName = outputDir + "wormY2HInteraction2.0.32.txt";
//        String nonY2HOutFileName = outputDir + "wormNonY2HInteraction2.0.32.txt";
//        boolean isForCel = true;
//        useSecondaryRef = true;
//        String primaryDBName = "FLYBASE";
//        String fileName = biogridDir + "BIOGRID-ORGANISM-Drosophila_melanogaster-2.0.50.psi25.xml";
//        String y2hOutFileName = outputDir + "flyY2HInteraction2.0.50.txt";
//        String nonY2HOutFileName = outputDir + "flyNonY2HInteraction2.0.50.txt"; 
//        useSecondaryRef = false; 
        String primaryDBName = "FLYBASE";
        String fileName = biogridDir + "BIOGRID-ORGANISM-Drosophila_melanogaster-2.0.50.psi25.xml";
        String y2hOutFileName = outputDir + "flyY2HInteraction2.0.50.txt";
        String nonY2HOutFileName = outputDir + "flyNonY2HInteraction2.0.50.txt"; 
        useSecondaryRef = false; 
        parse(primaryDBName, fileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Total time for parsing: " + (time2 - time1));
        outputInteractions(interactions, y2hOutFileName, nonY2HOutFileName);
    }

    protected void parse(String primaryDBName, String fileName) throws Exception {
        XMLInputFactory inputFactory = XMLInputFactory.newInstance();
        inputFactory.setProperty(XMLInputFactory.IS_COALESCING, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_REPLACING_ENTITY_REFERENCES, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_NAMESPACE_AWARE, Boolean.TRUE);
        inputFactory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, Boolean.TRUE);
        FileInputStream fis = new FileInputStream(fileName);
        XMLEventReader eventReader = inputFactory.createXMLEventReader(fis);
        while (eventReader.hasNext()) {
            XMLEvent xmlEvent = eventReader.nextEvent();
            handleXMLEvent(xmlEvent, primaryDBName);
        }
        System.out.println("Total experiments: " + expId2Name.size());
        System.out.println("Total interactors: " + interactorId2Accession.size());
        System.out.println("Total interactions: " + interactions.size());
        //analyzeInteractions(interactions);
    }
    
    protected void outputInteractions(List<Interaction> interactions, 
                                    String y2hOutFileName,
                                    String nonY2HOutFileName) throws IOException {
        FileUtility y2hFu = new FileUtility();
        y2hFu.setOutput(y2hOutFileName);
        FileUtility nony2hFu = new FileUtility();
        nony2hFu.setOutput(nonY2HOutFileName);
        Set<String> nony2hSet = new HashSet<String>();
        Set<String> y2hSet = new HashSet<String>();
        for (Interaction i : interactions) {
            Collections.sort(i.participants);
            String partner1 = i.participants.get(0);
            String partner2 = i.participants.get(1);
            String pair = partner1 + " " + partner2;
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
        int c = 0;
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
            // Check sequences
            List<String> sequences = interaction.sequences;
            for (String seq : sequences) 
                System.out.println(seq);
            for (String seq : sequences) {
                String code = checksum.calculateChecksum(seq);
                System.out.println(code);
            }
            if (c > 10)
                break;
        }
        int c0 = 0;
        int c3 = 0;
        for (Iterator<Interaction> it = interactions.iterator(); it.hasNext();) {
            Interaction tmp = it.next();
            if (tmp.participants == null ||
                tmp.participants.size() < 2) {
                it.remove();
                c0 ++;
            }
            else if (tmp.participants.size() > 2)
                c3 ++;
        }
        System.out.println("Interaction after cleanup: " + interactions.size());
        System.out.println("Interaction containing no interactors: " + c0);
        System.out.println("Interaction containing more than 2 interactors: " + c3);
    }
    
    private void handleXMLEvent(XMLEvent event, String primaryDBName) throws Exception {
        int eventType = event.getEventType();
        switch (eventType) {
            case START_DOCUMENT:
                System.out.println("Starting parsing the document...");
                break;
            case START_ELEMENT:
                StartElement startElement = event.asStartElement();
                handleStartElement(startElement, primaryDBName);
                break;               
            case END_ELEMENT :
                handleEndElement(event.asEndElement());
                break;
            case CHARACTERS :
                handleCharacters(event.asCharacters());
                break;
            case END_DOCUMENT :
                System.out.println("Ending parsing the document.");
                break;
        }
    }
    
    private void handleCharacters(Characters characters) throws Exception {
        if (!needText)
            return;
        String text = characters.getData().trim();
        if (isInExp) {
            String name = expId2Name.get(text);
            if (name == null) {
                throw new IllegalStateException("Experiment cannot be found for " +
                                                "interaction: " + text + " in " + interaction.id);
            }
            interaction.addExp(name);
        }
        else if (isInInteractorRef) {
            String accession = interactorId2Accession.get(text);
            if (accession != null)
                interaction.addParticipant(accession);
            String sequence  = interactorId2Sequence.get(text);
            if (sequence != null)
                interaction.addSequence(sequence);
        }
        if (isInExpList)
            expId2Name.put(id, text);
        if (isInInteractionType)
            interaction.type = text;
        if (isInSequence) {
            interactorId2Sequence.put(id, text);
        }
    }
    
    private void handleEndElement(EndElement endElement) throws Exception {
        String elmName = endElement.getName().getLocalPart();
        if (elmName.equals(PsiMiModel.experimentList))
            isInExpList = false;
        else if (elmName.equals(PsiMiModel.experimentRef)) {
            isInExp = false;
            needText = false;
        }
        else if (elmName.equals(PsiMiModel.interactorList))
            isInInteractorList = false;
        else if (elmName.equals(PsiMiModel.interactorRef)) {
            isInInteractorRef = false;
            needText = false;
        }
        else if (elmName.equals(PsiMiModel.interactionList))
            isInInteractionList = false;
        else if (elmName.equals(PsiMiModel.interaction)) {
            isInInteraction = false;
            interactions.add(interaction);
            interaction = null;
        }
        else if (elmName.equals(PsiMiModel.interactionDetectionMethod))
            isInCV = false;
        else if (elmName.equals(PsiMiModel.shortLabel))
            needText = false;
        else if (elmName.equals(PsiMiModel.interactor))
            isInInteractor = false;
        else if (elmName.equals(PsiMiModel.xref))
            isInxref = false;
        else if (elmName.equals(PsiMiModel.participantList))
            isInParticipantList = false;
        else if (elmName.equals(PsiMiModel.interactionType)) {
            isInInteractionType = false;
            isInCV = false;
        }
        else if (elmName.equals(PsiMiModel.sequence)) {
            isInSequence = false;
            needText = false;
        }
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
            interaction = new Interaction();
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
            if (ref == null) {
                // Try another way: <experimentRef>id</experimentRef>
                isInExp = true;
                needText = true;
            }
            else {
                String name = expId2Name.get(ref);
                if (name == null) {
                    throw new IllegalStateException("Experiment cannot be found for " +
                                                    "interaction: " + ref + " in " + interaction.id);
                }
                interaction.addExp(name);
            }
        }
        else if (elmName.equals(PsiMiModel.participantList)) {
            isInParticipantList = true;
        }
        else if (isInParticipantList && elmName.equals(PsiMiModel.interactorRef)) {
            String ref = getAttributeValue(startElement, "ref");
            if (ref == null) {
                // Need to try another way: <interactorRef>id</interactorRef>
                isInInteractorRef = true;
                needText = true;
            }
            else {
                String accession = interactorId2Accession.get(ref);
                if (accession != null)
                    interaction.addParticipant(accession);
            }
        }
        else if (elmName.equals(PsiMiModel.interactionType)) {
            isInInteractionType = true;
            isInCV = true;
        }
        // Use xref for dme and sce
        else if (!useSecondaryRef && isInxref && elmName.equals(PsiMiModel.primaryRef)) {
            String db = getAttributeValue(startElement, PsiMiModel.db);
            String acc = getAttributeValue(startElement, PsiMiModel.id);
            if (db.equals(primaryDBName))
                interactorId2Accession.put(id, acc);
        }
        else if (useSecondaryRef && isInxref && elmName.equals(PsiMiModel.secondaryRef)) {
            String db = getAttributeValue(startElement, PsiMiModel.db);
            String acc = getAttributeValue(startElement, PsiMiModel.id);
            if (db.equals(primaryDBName))
                interactorId2Accession.put(id, acc);
        }
        else if (isInExpList && elmName.equals(PsiMiModel.experimentDescription)) {
            id = getAttributeValue(startElement, PsiMiModel.id);
        }
        else if (isInCV && elmName.equals(PsiMiModel.shortLabel)) {
            needText = true;
        }
        else if (elmName.equals(PsiMiModel.sequence)) {
            isInSequence = true;
            needText = true;
        }
    }
    
    protected String getAttributeValue(StartElement element, String attName) {
        for (Iterator it = element.getAttributes(); it.hasNext();) {
            Attribute att = (Attribute) it.next();
            String attName1 = att.getName().getLocalPart();
            if (attName1.equals(attName))
                return att.getValue();
        }
        return null;
    }
    
    /**
     * This method is used to merge all interactions together.
     */
    public void mergeInteractions() throws IOException {
        String[] fileNames = new String[] {
//                "flyNonY2HInteractionHsaUni.txt",
//                "wormY2HInteractionHsaUni.txt",
  //              "flyY2HInteractionHsaUni.txt"
                "yeastNonY2HInteractionHsaUni.txt",
                "yeastY2HInteractionHsaUni.txt"
        };
        FileUtility fu = new FileUtility();
        Set<String> merged = new HashSet<String>();
        for (String fileName : fileNames) {
            Set<String> loaded = fu.loadInteractions("results/interaction/" + fileName);
            merged.addAll(loaded);
        }
        System.out.println("Merged Orthologous Interactions: " + merged.size());
        //fu.saveInteractions(merged, "results/interaction/OrthoInteractions.txt");
        fu.saveInteractions(merged, "results/interaction/YeastInteractions.txt");
//        fileNames = new String[] {
//                "HumanNonY2H.txt",
//                "HumanY2H.txt"
//        };
//        merged.clear();
//        for (String fileName : fileNames) {
//            Set<String> loaded = fu.loadInteractions("results/interaction/" + fileName);
//            merged.addAll(loaded);
//        }
//        System.out.println("Merged Human Interactions: " + merged.size());
//        fu.saveInteractions(merged, "results/interaction/HumanInteractions.txt");
    }
    
    
    
    /**
     * An inner class to hold information for interaction
     */
    protected class Interaction {
        String id;
        List<String> participants;
        List<String> expNames;
        List<String> sequences;
        String type;
        
        Interaction() {
        }
        
        public void addSequence(String sequence) {
            if (sequences == null)
                sequences = new ArrayList<String>();
            sequences.add(sequence);
        }
        
        public void addParticipant(String id) {
            if (participants == null)
                participants = new ArrayList<String>();
            participants.add(id);
        }
        
        public void addExp(String expName) {
            if (expNames == null)
                expNames = new ArrayList<String>();
            expNames.add(expName);
        }
        
        public String toString() {
            StringBuilder builder = new StringBuilder();
            builder.append(id).append("\t");
            if (participants != null) {
                builder.append("    participant: ");
                for (String s : participants)
                    builder.append(s).append("\t");
            }
            if (expNames != null) {
                builder.append("    expNames: ");
                for (String s : expNames)
                    builder.append(s).append("\t");
            }
            builder.append("    type: ").append(type);
            return builder.toString();
        }
    }
    
    public void grepUniProts() throws IOException {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/";
        String[] fileNames = new String[] {
                dirName + "HPRD/PSI-MI/HPRD_SINGLE_PSIMI_060106.xml",
                //dirName + "BioGrid/GRID-ORGANISM-Homo_sapiens-2.0.19.psi25.xml"
        };
        FileUtility fu = new FileUtility();
        String line = null;
        Set<String> uniProtIds = new HashSet<String>();
        int index1, index2;
        String id;
        for (String fileName : fileNames) {
            fu.setInput(fileName);
            uniProtIds.clear();
            while ((line = fu.readLine()) != null) {
                if (line.trim().startsWith("<interactorList>")) {
                    while (true) {
                        line = fu.readLine();
                        if (line.trim().startsWith("<secondaryRef")) {
                            if ((line.indexOf("db=\"uniprot\"") > 0) ||
                                    (line.indexOf("db=\"TREMBL\"") > 0) ||
                                    (line.indexOf("db=\"SWISSPROT\"") > 0)) {
                                index1 = line.indexOf("id=");
                                index2 = line.indexOf("\"", index1 + 5);
                                id = line.substring(index1 + 4, index2);
                                uniProtIds.add(id);
                            }
                        }
//                      if (line.trim().startsWith("<primaryRef")) {
//                      // want to get primary db name
//                      index1 = line.indexOf("db=");
//                      index2 = line.indexOf("\"", index1 + 4);
//                      id = line.substring(index1 + 4, index2);
//                      uniProtIds.add(id);
//                      }
                        if (line.trim().startsWith("</interactorList>"))
                            break;
                    }
                }
            }
            fu.close();
            System.out.println(fileName + ": " + uniProtIds.size());
        }
        String mapFileName = dirName + "UniProt/ACIDMap.txt";
        Map<String, String> acidMap = fu.importMap(mapFileName);
        Set<String> mappedSet = new HashSet<String>();
        String mapped = null;
        for (String tmp : uniProtIds) {
            mapped = acidMap.get(tmp);
            if (mapped == null) {
                System.out.println(tmp + " cannot be mapped!");
            }
            else
                mappedSet.add(mapped);
        }
        System.out.println("Mapped: " + mappedSet.size());
        Set<String> idsInPathways = fu.loadSet("results/UniProtInPathways.txt");
        System.out.println("Total Ids in Pathways: " + idsInPathways.size());
        mappedSet.addAll(idsInPathways);
        System.out.println("Combined: " + mappedSet.size());
    }
    
}
