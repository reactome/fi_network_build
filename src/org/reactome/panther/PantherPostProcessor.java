/*
 * Created on Mar 14, 2006
 *
 */
package org.reactome.panther;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.gk.database.util.LiteratureReferenceAttributeAutoFiller;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.slicing.SlicingEngine;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.PostProcessTemplate;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;

/**
 * This class is used to add UniProt instances to EntityWithAccessinedSequence instances
 * that are converted from Panther. The annotation is based on the mapping file provided
 * by the Panther web site, SequenceAssociationPathway1.13.
 * @author guanming
 *
 */
public class PantherPostProcessor extends PostProcessTemplate {
    //private GKInstance uniProtDBInstance;
    //private Map<String, GKInstance> idInstanceMap;
    
    public PantherPostProcessor() {
        //idInstanceMap = new HashMap<String, GKInstance>();
    }
    
    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor, 
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        super.attachDataSource("pantherdb", 
                               "http://www.pantherdb.org/",
                               dbAdaptor,
                               fileAdaptor);
    }

    @Override
    protected void processEWAS(MySQLAdaptor dbAdaptor, 
                               XMLFileAdaptor fileAdaptor) throws Exception {
        annotate(fileAdaptor,
                 dbAdaptor);
    }
    
    

    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor, 
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        super.processEntityCompartment(ReactomeJavaConstants.name, 
                                       dbAdaptor, 
                                       fileAdaptor);
    }
    
    /**
     * Use to process Panther only tasks.
     * @param dbAdaptor
     * @param fileAdaptor
     * @throws Exception
     */
    public void otherProcesses(MySQLAdaptor dbAdaptor,
                              XMLFileAdaptor fileAdaptor) throws Exception {
        processSimpleEntities(dbAdaptor, fileAdaptor);
        // As of December 23, 2020, disable literature reference processing since it is too slow after switching to
        // EUtil for auto-processing PMCIDs. 
//        PostProcessHelper.processLiteratureReferences(dbAdaptor, fileAdaptor);
        attachSpecies(dbAdaptor, fileAdaptor);
    }
    
    private void attachSpecies(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance humanSpecies = fileAdaptor.fetchInstance(new Long(48887L));
        if (humanSpecies == null) {
            GKInstance dbHuman = dbAdaptor.fetchInstance(new Long(48887L));
            humanSpecies = PostProcessHelper.downloadDBInstance(dbHuman, fileAdaptor);
        }
        if (humanSpecies == null)
            throw new IllegalStateException("PantherPostProcessor.attachSpecies(): cannot find homo sapiens.");
        String[] clses = new String[] {
                ReactomeJavaConstants.Event,
                ReactomeJavaConstants.Complex,
                ReactomeJavaConstants.EntitySet,
                ReactomeJavaConstants.GenomeEncodedEntity,
                ReactomeJavaConstants.Polymer
        };
        for (String cls : clses) {
            Collection c = fileAdaptor.fetchInstancesByClass(cls);
            if (c == null || c.size() == 0)
                continue;
            for (Iterator it = c.iterator(); it.hasNext();) {
                GKInstance instance = (GKInstance) it.next();
                if (instance.getDBID() > -1)
                    continue; // Instances from gk_central. Ignore them
                GKInstance species = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.species);
                if (species != null)
                    continue; // don't need to add
                instance.setAttributeValue(ReactomeJavaConstants.species,
                                           humanSpecies);
            }
        }
    }
    
    /**
     * This method is created to meet the requirement for Panther converting. It should be
     * pushed up the PostProcessTemplate in the future so that converting from other 
     * data sources can be done too.
     * @param dbAdaptor
     * @param fileAdaptor
     * @throws Exception
     */
    public void processSimpleEntities(MySQLAdaptor dbAdaptor,
                                      XMLFileAdaptor fileAdaptor) throws Exception {
        SmallMoleculeHandler handler = new SmallMoleculeHandler();
        handler.processSimpleEntities(dbAdaptor, fileAdaptor);
    }
    
    private GKInstance createFlagInstanceEdit(String flag,
                                              XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance ie = fileAdaptor.createNewInstance(ReactomeJavaConstants.InstanceEdit);
        ie.setAttributeValue(ReactomeJavaConstants.note, flag);
        ie.setDisplayName(flag);
        return ie;
    }

    @SuppressWarnings("unchecked")
    public void dumpToDB(String rtpjFileName, MySQLAdaptor targetDBA) throws Exception {
        XMLFileAdaptor fileLoader = new XMLFileAdaptor();
        fileLoader.setSource(rtpjFileName);
        // Get all instances
        Collection instances = fileLoader.fetchInstancesByClass(ReactomeJavaConstants.DatabaseObject);
        // Find the largest DB_ID
        long largest = 0;
        GKInstance instance = null;
        for (Iterator it = instances.iterator(); it.hasNext();) {
            instance = (GKInstance) it.next();
            if (instance.getDBID() > largest)
                largest = instance.getDBID();
        }
        System.out.println("Find the largest DB_ID: " + largest);
        // Want to convert negative to positive id
        // Use as a marker for panther
        GKInstance pantherIE = createFlagInstanceEdit("panther", fileLoader);
        instances.add(pantherIE);
        GKInstance proxyIE = createFlagInstanceEdit("proxy", fileLoader);
        instances.add(proxyIE);
        for (Iterator it = instances.iterator(); it.hasNext();) {
            instance = (GKInstance) it.next();
            if (instance.isShell()) {
                instance.setIsShell(false);
                instance.setAttributeValue(ReactomeJavaConstants.created, 
                                           proxyIE);
                continue;
            }
            if (instance.getDBID() > 0)
                continue;
            instance.setAttributeValue(ReactomeJavaConstants.created,
                                       pantherIE);
            instance.setDBID(++largest);
        }
        System.out.println("The largest DB_ID after reassigning DB_ID: " + largest);
        SlicingEngine slingEngine = new SlicingEngine();
        boolean isTnSupported = targetDBA.supportsTransactions();
        if (isTnSupported)
            targetDBA.startTransaction();
        try {
            for (Iterator it = instances.iterator(); it.hasNext();) {
                instance = (GKInstance) it.next();
                slingEngine.storeInstance(instance, targetDBA);
            }
            if (isTnSupported)
                targetDBA.commit();
        }
        catch (Exception e) {
            targetDBA.rollback();
            System.err.println("dumpToDB(): " + e);
            e.printStackTrace();
            throw e;
        }
    }

    public void annotate(String prjFileName,
                         MySQLAdaptor dbAdaptor) throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(prjFileName);
        annotate(fileAdaptor, dbAdaptor);
    }
    
    /**
     * Three things will be done in this method:
     * 1). Loading the mappping from Panther ID to UniProt ID
     * 2). Find ReferencePeptideSequence from the provided for mapped UniProt IDs.
     * 3). If one Panther ID is mapped to more than one UniProt ID, switch the original
     * type to DefinedSet, create EWAS for each UniProt ID, add these newly created EWAS
     * as members of DefinedSet.
     * @param fileAdaptor
     * @param dbAdaptor
     * @throws Exception
     * @throws IOException
     * @throws InvalidAttributeException
     */
    public void annotate(XMLFileAdaptor fileAdaptor,
                         MySQLAdaptor dbAdaptor) throws Exception {
        // Only EWAS instances are mapped
        Collection ewasCollection = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        if (ewasCollection == null || ewasCollection.size() == 0)
            return; // Nothing to do
        // Create a map between Panther ID and EWAS
        // One Panther ID can be mapped to multiple EWAS
        Map<String, Set<GKInstance>> pantherId2EWAS = new HashMap<String, Set<GKInstance>>();
        GKInstance gkInstance = null;
        GKInstance pantherIdInstance = null;
        for (Iterator it = ewasCollection.iterator(); it.hasNext();) {
            gkInstance = (GKInstance) it.next();
            pantherIdInstance = (GKInstance) gkInstance.getAttributeValue(ReactomeJavaConstants.crossReference);
            if (pantherIdInstance == null)
                continue; // It is possible an EWAS is not specified an ID
            GKInstance db = (GKInstance) pantherIdInstance.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            String dbName = (String) db.getAttributeValue(ReactomeJavaConstants.name);
            if (!dbName.equals("pantherdb"))
                continue;
            String identifier = (String) pantherIdInstance.getAttributeValue(ReactomeJavaConstants.identifier);
            Set<GKInstance> set = pantherId2EWAS.get(identifier);
            if (set == null) {
                set = new HashSet<GKInstance>();
                pantherId2EWAS.put(identifier, set);
            }
            set.add(gkInstance);
        }
        // Load the mapping from Panther IDs to UniProt IDs
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniIDMap = uniAnalyzer.loadUniProtIDsMap();
        PantherIdToUniProtMapper mapper = new PantherIdToUniProtMapper();
        Map<String, Set<String>> panther2UniMap = mapper.loadMapping(uniIDMap);
        // Map UniProt Ids to Reactome ReferencePeptideSequence instances
        Map<String, Set<GKInstance>> panther2RefPepSeqMap = new HashMap<String, Set<GKInstance>>();
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        Map<String, GKInstance> newRefPepSeqs = new HashMap<String, GKInstance>();
        //System.out.println("Size of pantherId2EWAS: " + pantherId2EWAS.size());
        for (Iterator<String> it = pantherId2EWAS.keySet().iterator(); it.hasNext();) {
            String pantherId = it.next();
            Set<String> uniIds = panther2UniMap.get(pantherId);
            // uniIds might be null
            if (uniIds == null) {
                //System.out.println("Panther ID, " + pantherId + ", cannot be mapped!");
                continue;
            }
            Set<GKInstance> set = new HashSet<GKInstance>();
            for (String uniId : uniIds) {
                GKInstance localRefPepSeq = null;
                GKInstance refPepSeq = fetchUniProtInstance(dbAdaptor, uniId);
                if (refPepSeq == null) {
                    localRefPepSeq = newRefPepSeqs.get(uniId);
                    if (localRefPepSeq == null) {
                        localRefPepSeq = createRefPepSeq(uniId, fileAdaptor);
                        newRefPepSeqs.put(uniId, localRefPepSeq);
                    }
                }
                else {
                    // Download refPepSeq to the local project
                    localRefPepSeq = PostProcessHelper.downloadDBInstance(refPepSeq, fileAdaptor);
                }
                set.add(localRefPepSeq);
            }
            panther2RefPepSeqMap.put(pantherId, set);
        }
        // Filled up newly created ReferencePeptideSequence instances.
        // Too many new ReferencePeptideSequence. Don't bother to fetch other information.
        //ReferencePeptideSequenceAutoFiller autoFiller = new ReferencePeptideSequenceAutoFiller();
        //autoFiller.setPersistenceAdaptor(fileAdaptor);
        //System.out.println("New Reference Peptide Sequence: " + newRefPepSeqs.size());
        // Add the minimum information
        fillNewRefPepSeqs(newRefPepSeqs.values(), fileAdaptor);
        //for (GKInstance refPepSeq : newRefPepSeqs.values()) {
        //    autoFiller.process(refPepSeq, null);
        //}
        // Attach the mapped ReferencePeptideSequence to EWAS
        for (Iterator<String> it = pantherId2EWAS.keySet().iterator(); it.hasNext();) {
            String pantherId = it.next();
            Set<GKInstance> ewasSet = pantherId2EWAS.get(pantherId);
            Set<GKInstance> refPepSeqSet = panther2RefPepSeqMap.get(pantherId);
            if (refPepSeqSet == null)
                continue;
            for (GKInstance ewas : ewasSet)
                annotateEwas(ewas, refPepSeqSet, fileAdaptor);
        }
    }
    
    private GKInstance fetchUniProtInstance(MySQLAdaptor dbAdaptor,
                                           String uniId) throws Exception {
        Collection refPepSeqCollection = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                            ReactomeJavaConstants.identifier,
                                                                            "=",
                                                                            uniId);
        if (refPepSeqCollection == null || refPepSeqCollection.size() == 0) {
            return null;
        }
        GKInstance refPepSeq = (GKInstance) refPepSeqCollection.iterator().next();
        return refPepSeq;
    }
    
    private void fillNewRefPepSeqs(Collection<GKInstance> newRefPepSeqs,
                                   XMLFileAdaptor fileAdaptor) throws Exception {
        // UniProt
        GKInstance uniProt = fileAdaptor.fetchInstance(2L);
        assert (uniProt != null);
        // Homo sapiens
        GKInstance human = fileAdaptor.fetchInstance(48887L);
        assert (human != null);
        String identifier = null;
        // These should be solved automatically after ReferencePeptideSequences are updated 
        // automatically
        // Do a query for the UniProt web site
//        ReferencePeptideSequenceAutoFiller autoFiller = new ReferencePeptideSequenceAutoFiller();
//        autoFiller.setPersistenceAdaptor(fileAdaptor);
//        long time1 = System.currentTimeMillis();
//        for (GKInstance refPepSeq : newRefPepSeqs) {
//            autoFiller.process(refPepSeq);
//            InstanceDisplayNameGenerator.setDisplayName(refPepSeq);
//        }
//        long time2 = System.currentTimeMillis();
        for (GKInstance refPepSeq : newRefPepSeqs) {
            identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            refPepSeq.addAttributeValue(ReactomeJavaConstants.name,
                                        identifier);
            refPepSeq.setAttributeValue(ReactomeJavaConstants.species,
                                        human);
            refPepSeq.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                        uniProt);
            InstanceDisplayNameGenerator.setDisplayName(refPepSeq);
        }
    }
    
    private GKInstance createRefPepSeq(String uniId,
                                       XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance refPepSeq = fileAdaptor.createNewInstance(ReactomeJavaConstants.ReferenceGeneProduct);
        refPepSeq.setAttributeValue(ReactomeJavaConstants.identifier,
                                    uniId);
        return refPepSeq;
    }
    
    private void annotateEwas(GKInstance ewas,
                              Set<GKInstance> refPepSeqSet,
                              XMLFileAdaptor fileAdaptor) throws Exception {
        if (refPepSeqSet.size() == 1) {
            // Just attach to EWAS
            GKInstance refPepSeq = refPepSeqSet.iterator().next();
            ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity,
                                   refPepSeq);
        }
        else {
            // It should have more than 1. The minimum is 1.
            // Switch the type of EWAS to DefinedSet
            PostProcessHelper.switchEWASToSet(ewas, refPepSeqSet, fileAdaptor);
        }
    }
    
    public void fillUpLiteratureReferences() throws Exception {
        String inFile = "/Users/wgm/Documents/gkteam/bernard/msb4100057-s1.rtpj";
        String outFile = "/Users/wgm/Documents/gkteam/bernard/msb4100057-s1-v2.rtpj";
        // Load the project file
        XMLFileAdaptor adaptor = new XMLFileAdaptor();
        // This method call will load the project
        adaptor.setSource(inFile);
        Collection litRefs = adaptor.fetchInstancesByClass(ReactomeJavaConstants.LiteratureReference);
        if (litRefs == null || litRefs.size() == 0)
            return;
        LiteratureReferenceAttributeAutoFiller autoFiller = new LiteratureReferenceAttributeAutoFiller();
        autoFiller.setPersistenceAdaptor(adaptor);
        GKInstance litRef = null;
        for (Iterator it = litRefs.iterator(); it.hasNext();) {
            litRef = (GKInstance) it.next();
            autoFiller.process(litRef, null);
            InstanceDisplayNameGenerator.setDisplayName(litRef);
        }
        adaptor.save(outFile);
    }    
}
