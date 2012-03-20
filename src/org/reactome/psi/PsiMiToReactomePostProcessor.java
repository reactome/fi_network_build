/*
 * Created on Aug 22, 2006
 *
 */
package org.reactome.psi;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.gk.database.SynchronizationManager;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.PersistenceManager;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.GKSchemaAttribute;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.PostProcessTemplate;
import org.reactome.data.UniProtAnalyzer;


public class PsiMiToReactomePostProcessor extends PostProcessTemplate {
    private final static Logger logger = Logger.getLogger(PsiMiToReactomePostProcessor.class);
    private String dataSourceName;
    private String dataSourceUrl;
    
    public PsiMiToReactomePostProcessor() {
    }
    
    public String getDataSourceName() {
        return dataSourceName;
    }

    public void setDataSourceName(String dataSourceName) {
        this.dataSourceName = dataSourceName;
    }

    public String getDataSourceUrl() {
        return dataSourceUrl;
    }

    public void setDataSourceUrl(String dataSourceUrl) {
        this.dataSourceUrl = dataSourceUrl;
    }

    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor, 
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        if (dataSourceName != null && dataSourceUrl != null) {
            attachDataSource(dataSourceName,
                             dataSourceUrl,
                             dbAdaptor,
                             fileAdaptor);
            return;
        }
        // Start from Interaction. The data source can be fetched from interaction
        Collection interactionCollection = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        GKInstance interaction = null;
        for (Iterator it = interactionCollection.iterator(); it.hasNext();) {
            interaction = (GKInstance) it.next();
            attachDataSourceToInteraction(interaction);
        }
    }
    
    private void attachDataSourceToInteraction(GKInstance interaction) throws Exception {
        // Get ReferenceDatabase from crossReference
        GKInstance crossReference = (GKInstance) interaction.getAttributeValue(ReactomeJavaConstants.crossReference);
        if (crossReference == null)
            return;
        GKInstance refDB = (GKInstance) crossReference.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
        if (refDB == null)
            return;
        addDataSource(crossReference, refDB);
        addDataSource(interaction, refDB);
        // Check interactors
        List interactors = interaction.getAttributeValuesList(ReactomeJavaConstants.interactor);
        for (Iterator it = interactors.iterator(); it.hasNext();) {
            GKInstance interactor = (GKInstance) it.next();
            addDataSource(interactor, refDB);
            GKInstance refEntity = (GKInstance) interactor.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            addDataSource(refEntity, refDB);
        }
        // Check summations
        List summations = interaction.getAttributeValuesList(ReactomeJavaConstants.summation);
        if (summations != null) {
            for (Iterator it = summations.iterator(); it.hasNext();) {
                GKInstance summation = (GKInstance) it.next();
                addDataSource(summation, refDB);
                List ltRefs = summation.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
                if (ltRefs != null) {
                    for (Iterator it1 = ltRefs.iterator(); it1.hasNext();) {
                        GKInstance lit = (GKInstance) it1.next();
                        addDataSource(lit, refDB);
                    }
                }
            }
        }
    }
    
    
    private void addDataSource(GKInstance gkInstance, GKInstance refDB) throws Exception {
        // For the time being, add dataSource once only even though it can come from more than
        // one data source. This issue will be figured out in the future.
        GKInstance ds = (GKInstance) gkInstance.getAttributeValue(ReactomeJavaConstants.dataSource);
        if (ds == null)
            gkInstance.addAttributeValue(ReactomeJavaConstants.dataSource,
                                         refDB);
    }

    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor, 
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        // No need to do anything here.
    }

    @Override
    protected void processEWAS(MySQLAdaptor dbAdaptor, 
                               XMLFileAdaptor fileAdaptor) throws Exception {
        Collection ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        SynchronizationManager manager = SynchronizationManager.getManager();
        PersistenceManager.getManager().setActiveFileAdaptor(fileAdaptor);
        // Want to remove any redundant ids
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> idMap = uniAnalyzer.loadUniProtIDsMap();
        // Need to make a copy in case the local set changed
        List<GKInstance> list = new ArrayList<GKInstance>(ewases);
        logger.info("Process ReferencePeptideSequences...");
        for (GKInstance ewas : list) {
            GKInstance referenceEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (referenceEntity == null)
                continue;
            String identifier = (String) referenceEntity.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier == null)
                continue;
            String mapped = idMap.get(identifier);
            if (mapped == null)
                mapped = identifier;
            GKInstance localRefPepSeq = queryAppropriateLocal(mapped, fileAdaptor);
            if (localRefPepSeq == null || localRefPepSeq.getDBID() < 0) {
                // Need to download
                GKInstance dbRefPepSeq = PostProcessHelper.queryRefPepSeq(mapped,
                                                                          dbAdaptor);
                if (dbRefPepSeq != null) {
                    if (localRefPepSeq == null)
                        localRefPepSeq = PostProcessHelper.downloadDBInstance(dbRefPepSeq, 
                                                                              fileAdaptor);
                    else
                        PostProcessHelper.updateFromDB(localRefPepSeq, 
                                                       dbRefPepSeq,
                                                       manager);
                }
            }
            if (localRefPepSeq != null) // In case the newly mapped identifier cannot be found
                ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity,
                                       localRefPepSeq);
        }
        // Clean up any RefPepSeqs
        logger.info("Clean up ReferencePeptideSequences...");
        cleanUpRefPepSeqs(fileAdaptor);
    }
    
    private GKInstance queryAppropriateLocal(String uniProtId,
                                             XMLFileAdaptor localAdaptor) throws Exception {
        Collection c = null;
        if (uniProtId.contains("-"))
            c = localAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceIsoform, 
                                                     ReactomeJavaConstants.variantIdentifier,
                                                     "=",
                                                     uniProtId);
        else
            c = localAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                 ReactomeJavaConstants.identifier,
                                                 "=",
                                                 uniProtId);
        GKInstance refPepSeq = null;
        if (c != null && c.size() > 0) {
            // Check if a db instance existed in the local project
            for (Iterator it = c.iterator(); it.hasNext();) {
                GKInstance tmp = (GKInstance) it.next();
                if (tmp.getDBID() > 0) {
                    refPepSeq = tmp;
                    break;
                }
            }
            // If a download refPepSeq cannot be found. Just use any local
            if (refPepSeq == null)
                refPepSeq = (GKInstance) c.iterator().next();
        }
        return refPepSeq;
    }
    
    protected void useDBInstanceAsShell(GKInstance db, 
                                        GKInstance local,
                                        XMLFileAdaptor fileAdaptor) throws Exception {
        Long oldId = local.getDBID();
        // Delete all attributes in the local
        for (Iterator it = local.getSchemClass().getAttributes().iterator(); it.hasNext();) {
            GKSchemaAttribute att = (GKSchemaAttribute) it.next();
            local.setAttributeValueNoCheck(att, null);
        }
        local.setDBID(db.getDBID());
        fileAdaptor.dbIDUpdated(oldId, local);
        local.setDisplayName(db.getDisplayName());
        local.setIsShell(true);
        local.setIsDirty(false);
    }
    
    protected GKInstance fetchReferenceSequence(String identifier, MySQLAdaptor dbAdaptor) throws Exception {
        Collection c = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceSequence, 
                                                          ReactomeJavaConstants.identifier, 
                                                          "=", 
                                                          identifier);
        if (c != null && c.size() > 0)
            return (GKInstance) c.iterator().next();
        return null;
    }
}
