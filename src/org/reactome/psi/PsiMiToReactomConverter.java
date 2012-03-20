/*
 * Created on Aug 21, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.schema.InvalidAttributeValueException;
import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;

/**
 * Extract protein protein interaction from PSI-MI to Reactome data model.
 * This is a shallow extraction, not all of information in PSI-MI is extracted.
 * @author guanming
 *
 */
public class PsiMiToReactomConverter implements PsiMiPostParseProcessor {
    static Logger logger = Logger.getLogger(PsiMiToReactomConverter.class);
    
    private SessionFactory sessionFactory;
    private XMLFileAdaptor fileAdaptor;
    private MySQLAdaptor dbAdaptor;
    private PsiMiToReactomePostProcessor postProcessor;
    // Cache these GKInstances to improve the performance by avoiding calling
    // XMLFileAdaptor.fetchInstanceByAttribute()
    private Map<String, GKInstance> dbName2RefDB = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> dbNameId2RefSeq = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> refSeqIdDef2EWAS = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> pmid2Lit = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> journal2Lit = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> litText2Summation = new HashMap<String, GKInstance>();
    private Map<String, GKInstance> dbNameId2DI = new HashMap<String, GKInstance>();

    public PsiMiToReactomConverter() {
        //postProcessor = new PsiMiToReactomePostProcessor();
    }
    
    public void setPostProcessor(PsiMiToReactomePostProcessor processor) {
        this.postProcessor = processor;
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dbAdaptor) {
        this.dbAdaptor = dbAdaptor;
    }
    
    public void setSessionFactory(SessionFactory sf) {
        this.sessionFactory = sf;
    }
    
    public SessionFactory getSessionFactory() {
        return this.sessionFactory;
    }
    
    public void setReactomeAdaptor(XMLFileAdaptor fileAdaptor) {
        this.fileAdaptor = fileAdaptor;
    }
    
    public XMLFileAdaptor getReactomeAdaptor() {
        return this.fileAdaptor;
    }
    
    public void postProcess(Object obj, 
                            PsiMiModel model) throws Exception {
        // Will handle Entry object only
        if (!(obj instanceof Entry))
            return;
        logger.info("Starting converting...");
        long time1 = System.currentTimeMillis();
        Entry entry = (Entry) obj;
        List<Interaction> interactions = entry.getInteractionList();
        for (Interaction interaction : interactions)
            convert(interaction);
        postProcess();
        long time2 = System.currentTimeMillis();
        logger.info("Converting is done!");
        logger.info("Concerting time: " + (time2 - time1));
    }

    /**
     * Convert PSI-MI interactions to Reactome interactions.
     * @throws Exception
     */
    public void convert() throws Exception {
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from " + Interaction.class.getName());
        int start = 0;
        int max = 2000;
        max = 1000;
        // Use pagination query to minimize the memory usage.
        query.setFirstResult(start);
        query.setMaxResults(max);
        List list = query.list();
        Interaction interaction = null;
        GKInstance gkInteraction = null;
        int total = 0;
        while (list.size() > 0) {
            total += list.size();
            for (Iterator it = list.iterator(); it.hasNext();) {
                interaction = (Interaction) it.next();
                gkInteraction = convert(interaction);
            }
            session.clear(); // To empty cache to keep the memory usage small.
            start += list.size();
            query.setFirstResult(start);
            list = query.list();
            break; // JUST FOR TESTING
        }
        session.close();
        postProcess();
        System.out.println("Total Interactions: " + total);
    }
    
    public void postProcess() throws Exception {
        if (postProcessor != null)
            postProcessor.postProcess(dbAdaptor, fileAdaptor);
    }
    
    private GKInstance convert(Interaction interaction) throws Exception {
        GKInstance gkInteraction = fileAdaptor.createNewInstance(ReactomeJavaConstants.Interaction);
        extractNames(interaction.getNames(), gkInteraction);
        extractType(interaction, gkInteraction);
        extractParticipants(interaction, gkInteraction);
        extractXref(interaction.getXref(), gkInteraction);
        convertExperiments(interaction, gkInteraction);
        return gkInteraction;
    }
    
    private void convertExperiments(Interaction interaction,
                                    GKInstance gkInteraction) throws Exception {
        List<Experiment> experiments = interaction.getExperimentList();
        if (experiments == null || experiments.size() == 0)
            return;
        for (Experiment exp : experiments) {
            GKInstance gkSummation = convertExperiment(exp);
            gkInteraction.addAttributeValue(ReactomeJavaConstants.summation,
                                            gkSummation);
        }
    }
    
    private GKInstance convertExperiment(Experiment exp) throws Exception {
        // Get bibxref
        Bibref bibref = exp.getBibref();
        GKInstance gkLitRef = null;
        if (bibref != null) {
            Xref xref = bibref.getXref();
            DbReference dbref = xref.getPrimaryRef();
            gkLitRef = convertDbRefToLiteratureRef(dbref);
        }
        StringBuilder builder = new StringBuilder();
        builder.append("Experiment: ");
        Names names = exp.getNames();
        if (names != null)
            builder.append(names.getShortLabel());
        OpenCV method = exp.getInteractionDetectionMethod();
        if (method != null) {
            builder.append(". Detection method: ").append(method.getNames().getShortLabel());
        }
        String text = builder.toString();
        String key = gkLitRef.getDBID() + ";" + builder.toString();
        GKInstance gkInstance = litText2Summation.get(key);
        if (gkInstance != null)
            return gkInstance;
        gkInstance = fileAdaptor.createNewInstance(ReactomeJavaConstants.Summation);
        gkInstance.setAttributeValue(ReactomeJavaConstants.text, text);
        if (gkLitRef != null)
            gkInstance.addAttributeValue(ReactomeJavaConstants.literatureReference, gkLitRef);
        litText2Summation.put(key, gkInstance);
        return gkInstance;
    }
    
    private GKInstance convertDbRefToLiteratureRef(DbReference ref) throws Exception {
        if (ref.getDb().equalsIgnoreCase("pubmed")) {
            GKInstance rtn = pmid2Lit.get(ref.getId());
            if (rtn != null)
                return rtn;
            Integer id = null;
            try {
                id = Integer.valueOf(ref.getId());
            }
            catch(Exception e) {
                System.out.println(e + ": " + id);
            }
            GKInstance gkInstance = fileAdaptor.createNewInstance(ReactomeJavaConstants.LiteratureReference);
            if (id != null)
                gkInstance.setAttributeValue(ReactomeJavaConstants.pubMedIdentifier, id);
            pmid2Lit.put(ref.getId(), gkInstance);
            return gkInstance;
        }
        else {
            // Maybe from other source
            String journal = ref.getDb() + ": " + ref.getId();
            GKInstance rtn = journal2Lit.get(journal);
            if (rtn != null)
                return rtn;
            GKInstance gkInstance = fileAdaptor.createNewInstance(ReactomeJavaConstants.LiteratureReference);
            gkInstance.setAttributeValue(ReactomeJavaConstants.journal, 
                                         journal);
            journal2Lit.put(journal, gkInstance);
            return gkInstance;
        }
    }
    
    private void extractXref(Xref xref, 
                             GKInstance gkInteraction) throws Exception {
        if (xref == null)
            return;
        DbReference primaryRef = xref.getPrimaryRef();
        GKInstance gkDbIdentifier = getDBIdentifier(primaryRef);
        gkInteraction.addAttributeValue(ReactomeJavaConstants.crossReference,
                                        gkDbIdentifier);
        // Check if there are any secondayRef
        List<DbReference> secondaryRefList = xref.getSecondaryRefList();
        if (secondaryRefList == null)
            return;
        for (DbReference dbRef : secondaryRefList) {
            gkDbIdentifier = getDBIdentifier(dbRef);;
            gkInteraction.addAttributeValue(ReactomeJavaConstants.crossReference,
                                            gkDbIdentifier);
        }
    }
    
    private GKInstance getDBIdentifier(DbReference dbReference) throws Exception {
        String dbName = dbReference.getDb();
        String id = dbReference.getId();
        String key = dbName + ":" + id;
        GKInstance rtn = dbNameId2DI.get(key);
        if (rtn != null)
            return rtn;
        rtn = fileAdaptor.createNewInstance(ReactomeJavaConstants.DatabaseIdentifier);
        rtn.setAttributeValue(ReactomeJavaConstants.identifier,
                              id);
        if (dbReference.getDb() != null) {
            GKInstance refDb = getReferenceDatabase(dbReference.getDb());
            rtn.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                  refDb);
        }
        return rtn;
    }
    
    private GKInstance getReferenceDatabase(String dbName) throws Exception {
        GKInstance gkInstance = dbName2RefDB.get(dbName);
        if (gkInstance != null)
            return gkInstance;
        gkInstance = fileAdaptor.createNewInstance(ReactomeJavaConstants.ReferenceDatabase);
        gkInstance.setAttributeValue(ReactomeJavaConstants.name,
                                     dbName);
        gkInstance.setDisplayName(dbName);
        dbName2RefDB.put(dbName, gkInstance);
        return gkInstance;
    }
    
    private void extractParticipants(Interaction interaction,
                                     GKInstance gkInteraction) throws Exception {
        List<Participant> participants = interaction.getParticipantList();
        if (participants == null || participants.size() == 0)
            return;
        for (Participant participant : participants) {
            GKInstance gkInstance = convert(participant);
            gkInteraction.addAttributeValue(ReactomeJavaConstants.interactor,
                                            gkInstance);
        }
    }
    
    @SuppressWarnings("unchecked")
    private GKInstance convert(Participant participant) throws Exception {
        // Biological is used as a defintion
        String def = convertBiologicalRole(participant);
        Interactor interactor = participant.getInteractor();
        // Search for gkParticipant based on refEntity. Just ignore other values.
        GKInstance refEntity = null;
        //if (interactor != null) {
            refEntity = convert(interactor);
            String key = refEntity.getDBID() + ";" + def;
            GKInstance gkInstance = refSeqIdDef2EWAS.get(key);
            if (gkInstance != null)
                return gkInstance;
        GKInstance gkParticipant = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        if (def != null)
            gkParticipant.addAttributeValue(ReactomeJavaConstants.definition,
                                            def);
        gkParticipant.setAttributeValue(ReactomeJavaConstants.referenceEntity,
                                        refEntity);
        // In case nothing is displayed 
        List names = refEntity.getAttributeValuesList(ReactomeJavaConstants.name);
        if (names != null && names.size() > 0)
            gkParticipant.setAttributeValueNoCheck(ReactomeJavaConstants.name, 
                                                   new ArrayList(names));
        extractXref(participant.getXref(), gkParticipant);
        refSeqIdDef2EWAS.put(key, gkParticipant);
        return gkParticipant;
    }
    
    private GKInstance convert(Interactor interactor) throws Exception {
        GKInstance gkRefDb = null;
        String identifier = null;
        Xref xref = interactor.getXref();
        if (xref != null) {
            DbReference dbRef = getUniProtRef(xref);
            if (dbRef != null) {
                gkRefDb = getReferenceDatabase("UniProt");
            }
            else {
                dbRef = xref.getPrimaryRef();
                gkRefDb = getReferenceDatabase(dbRef.getDb());
            }
            identifier = dbRef.getId();
        }
        GKInstance refEntity = null;
        // Note: gkRefDb might be null
        String key = gkRefDb + ":" + identifier;
        if (identifier != null && gkRefDb != null) {
            refEntity = dbNameId2RefSeq.get(key);
            if (refEntity != null)
                return refEntity;
        }
        OpenCV type = interactor.getInteractorType();
        String typeName = type.getNames().getShortLabel();
        String gkClsName = null;
        if (typeName.equals("protein") || typeName.equals("peptide")) // || typeName.equals("unknown participant")) // Just a hack for Manuel's file
            gkClsName = ReactomeJavaConstants.ReferenceGeneProduct;
        else if (typeName.equals("deoxyribonucleic acid"))
            gkClsName = ReactomeJavaConstants.ReferenceDNASequence;
        else if (typeName.equals("ribonucleic acid") ||
                typeName.equals("rna"))
            gkClsName = ReactomeJavaConstants.ReferenceRNASequence;
        else // Use ReferenceSequence for the time being even for SmallMolecules, Complexes, etc.
            gkClsName = ReactomeJavaConstants.ReferenceSequence;
        refEntity = fileAdaptor.createNewInstance(gkClsName);
        extractNames(interactor.getNames(), refEntity);
        // Get UniProt if it exists
        refEntity.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                    gkRefDb);
        refEntity.setAttributeValue(ReactomeJavaConstants.identifier, 
                                    identifier);
        dbNameId2RefSeq.put(key, refEntity);
        return refEntity;
    }
    
    private DbReference getUniProtRef(Xref xref) {
        DbReference primaryRef = xref.getPrimaryRef();
        if (primaryRef.getDb().startsWith("uniprot")) {
            return primaryRef;
        }
        List<DbReference> secondaryRefList = xref.getSecondaryRefList();
        if (secondaryRefList != null && secondaryRefList.size() > 0) {
            for (DbReference ref : secondaryRefList)
                if (ref.getDb().startsWith("uniprot"))
                    return ref;
        }
        return null;
    }
    
    private String convertBiologicalRole(Participant participant) throws Exception {
        OpenCV role = participant.getBiologicalRole();
        if (role != null) {
            // Put in the defintion slot
            String def = role.getNames().getShortLabel();
            if (def.equals("unspecified"))
                return null; // Want to escape unspecified. Meaningless!
            return def;
        }
        return null;
    }
    
    private void extractType(Interaction interaction, GKInstance gkInteraction) throws InvalidAttributeException, InvalidAttributeValueException {
        OpenCV type = interaction.getInteractionType();
        if (type != null) {
            // Just extract name
            String typeName = type.getNames().getShortLabel();
            gkInteraction.setAttributeValue(ReactomeJavaConstants.interactionType,
                                            typeName);
        }
    }
    
    private void extractNames(Names names, 
                              GKInstance gkInteraction) throws Exception {
        // Get names
        if (names != null) {
            String name = names.getShortLabel();
            gkInteraction.addAttributeValue(ReactomeJavaConstants.name, 
                                            name);
            name = names.getFullName();
            gkInteraction.addAttributeValue(ReactomeJavaConstants.name, 
                                            name);
            List<String> aliases = names.getAlias();
            if (aliases != null)
                for (String name1 : aliases)
                    gkInteraction.addAttributeValue(ReactomeJavaConstants.name, 
                                                    name1);
        }
    }

}
