/*
 * Created on Jul 26, 2006
 *
 */
package org.reactome.b2rPostProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;
import org.gk.database.SynchronizationManager;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.biopax.BioPAXToReactomePostProcessor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.PostProcessTemplate;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to post process the mapped GKInstances from CellMap database.
 * ReferencePeptideSequence is mapped to UniProt, names, species and crossReference 
 * are copied to EWAS that refer to ReferencePeptideSequences.
 * @author guanming
 *
 */
public class HPRDBToRPostProcessor extends PostProcessTemplate implements BioPAXToReactomePostProcessor {
    private Logger logger = Logger.getLogger(HPRDBToRPostProcessor.class);
    
    public HPRDBToRPostProcessor() {
    }
    
//    public void postProcess(MySQLAdaptor dbAdaptor,
//                            XMLFileAdaptor fileAdaptor) throws Exception {
//        // mergePhysicalEntities(fileAdaptor, 
//        //                     ReactomeJavaConstants.EntityWithAccessionedSequence);
//    }
    
    protected void processEWAS(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        // This is a kind of cheating: For projects converted from BioPAX, ReferencePeptideSequence
        // should be processed.
        processReferencePeptideSequences(dbAdaptor, fileAdaptor);
    }
    
    protected void attachDataSource(MySQLAdaptor dbAdaptor,
                                   XMLFileAdaptor fileAdaptor) throws Exception {
        String dbName = "The Cancer Cell Map";
        String url = "http://www.cellmap.org";
        attachDataSource(dbName, url, dbAdaptor, fileAdaptor);
    }
    
//    private void mergePhysicalEntities(XMLFileAdaptor fileAdaptor,
//                                       String clsName) throws Exception {
//        logger.info("Start merging PhysicalEntities...");
//        Collection list = fileAdaptor.fetchInstancesByClass(clsName);
//        if (list == null || list.size() == 0)
//            return;
//        SchemaClass cls = fileAdaptor.getSchema().getClassByName(clsName);
//        if (!cls.isValidAttribute(ReactomeJavaConstants.referenceEntity))
//            return; // Don't need to do it.
//        // Hold map from ReferenceEntity to GKInstances
//        Map<GKInstance, List<GKInstance>> ref2Entities = new HashMap<GKInstance, List<GKInstance>>();
//        GKInstance gkInstance = null;
//        GKInstance ref = null;
//        logger.info("    creating map for merging...");
//        // The following loop is used to extract referrers to referenceEntity.
//        for (Iterator it = list.iterator(); it.hasNext();) {
//            gkInstance = (GKInstance) it.next();
//            ref = (GKInstance) gkInstance.getAttributeValue(ReactomeJavaConstants.referenceEntity);
//            if (ref != null) {
//                List<GKInstance> tmpList = ref2Entities.get(ref);
//                if (tmpList == null) {
//                    tmpList = new ArrayList<GKInstance>();
//                    ref2Entities.put(ref, tmpList);
//                }
//                tmpList.add(gkInstance);
//            }
//        }
//        logger.info("   merging actually starting...");
//        logger.info("   size of the merging map: " + ref2Entities.size());
//        int c = 0;
//        // Merge based on ReferenceEntity
//        for (Iterator<GKInstance> it = ref2Entities.keySet().iterator(); it.hasNext();) {
//            logger.info("    merging " + c++);
//            ref = it.next();
//            List<GKInstance> tmpList = ref2Entities.get(ref);
//            if (tmpList.size() < 2)
//                continue; // Don't need merging
//            merge(tmpList, fileAdaptor);
//        }
//    }
    
//    private void merge(List<GKInstance> instanceList, XMLFileAdaptor fileAdaptor) throws Exception {
//        GKInstance instance1;
//        GKInstance instance2;
//        // Note: instanceList.size will not change. However, null will be
//        // placed into the list.
//        for (int i = 0; i < instanceList.size() - 1; i++) {
//            instance1 = instanceList.get(i);
//            if (instance1 == null)
//                continue;
//            for (int j = i + 1; j < instanceList.size(); j++) {
//                instance2 = instanceList.get(j);
//                if (instance2 == null)
//                    continue;
//                if (merge(instance1, instance2)) {
//                    instanceList.set(j, null);
//                    fileAdaptor.deleteInstance(instance2);
//                }
//            }
//        }
//    }
    
    /**
     * target will be kept and source should be merged away.
     * @param target
     * @param source
     * @throws Exception
     */
//    private boolean merge(GKInstance target, GKInstance source) throws Exception {
//        // Check if target and source can be merged based on defining attributes
//        GKSchemaClass cls = (GKSchemaClass) target.getSchemClass();
//        Collection definingAtts = cls.getDefiningAttributes();
//        SchemaAttribute att = null;
//        List attValues1 = null;
//        List attValues2 = null;
//        for (Iterator it = definingAtts.iterator(); it.hasNext();) {
//            att = (SchemaAttribute) it.next();
//            attValues1 = target.getAttributeValuesList(att);
//            attValues2 = source.getAttributeValuesList(att);
//            if(!InstanceUtilities.compareAttValues(attValues1, attValues2, att))
//                return false;
//        }
//        // Source can be merged away
//        XMLFileAdaptor fileAdaptor = (XMLFileAdaptor) source.getDbAdaptor();
//        Collection referrers = fileAdaptor.getReferers(source);
//        for (Iterator it = referrers.iterator(); it.hasNext();) {
//            GKInstance referrer = (GKInstance) it.next();
//            InstanceUtilities.replaceReference(referrer, source, target);
//        }
//        return true;
//    }

    protected void processReferencePeptideSequences(MySQLAdaptor dbAdaptor, XMLFileAdaptor fileAdaptor) 
                                                  throws IOException, Exception {
        logger.info("Starting process ReferencePeptideSequence...");
        FileUtility fu = new FileUtility();
        Map<String, String> hprd2UniMap = fu.importMap("resources/HPRD2UniProt.txt");
        Map<String, GKInstance> hprd2LocalGKMap = mapHPRD2LocalGK(fileAdaptor);
        Map<GKInstance, GKInstance> local2DBMap = mapLocal2DB(hprd2LocalGKMap,
                                                              hprd2UniMap,
                                                              dbAdaptor);
        copyProperties(local2DBMap, fileAdaptor);
    }
    
    protected void copyProperties(Map<GKInstance, GKInstance> local2DBMap,
                                XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance localGK = null;
        GKInstance dbGK = null;
        SynchronizationManager manager = SynchronizationManager.getManager();
        for (Iterator<GKInstance> it = local2DBMap.keySet().iterator(); it.hasNext();) {
            localGK = it.next();
            dbGK = local2DBMap.get(localGK);
            copyRefPepPropertiesToReferrers(localGK);
            PostProcessHelper.updateFromDB(localGK, dbGK, manager);
        }
    }
    
    private void copyRefPepPropertiesToReferrers(GKInstance refPep) throws Exception {
        Collection referrers = refPep.getReferers(ReactomeJavaConstants.referenceEntity);
        if (referrers == null || referrers.size() == 0)
            return;
        GKInstance referrer = null;
        for (Iterator it = referrers.iterator(); it.hasNext();) {
            referrer = (GKInstance) it.next();
            if (!referrer.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence))
                continue;
            copyAttributeValue(ReactomeJavaConstants.species, refPep, referrer);
            copyAttributeValue(ReactomeJavaConstants.crossReference, refPep, referrer);
            // Don't copy names. Maybe names in PhysicalEntities are used to indicate
            // something as in NCI-Pathways, e.g., ATM, ATM_Active. These names can be
            // used only in PhysicalEntities. (Jan 26, 2007).
            //copyAttributeValue(ReactomeJavaConstants.name, refPep, referrer);
            InstanceDisplayNameGenerator.setDisplayName(referrer);
        }
    }
    
    @SuppressWarnings("unchecked")
    protected void copyAttributeValue(String attName,
                                    GKInstance source, 
                                    GKInstance target) throws Exception {
        List list = source.getAttributeValuesList(attName);
        if (list != null && list.size() > 0)
            target.setAttributeValueNoCheck(attName, new ArrayList(list));
        else
            target.setAttributeValue(attName, null);
    }
    
    private Map<GKInstance, GKInstance> mapLocal2DB(Map<String, GKInstance> hprd2LocalGKMap,
                                                    Map<String, String> hprd2UniMap,
                                                    MySQLAdaptor dbAdaptor) throws Exception {
        Map<GKInstance, GKInstance> local2DB = new HashMap<GKInstance, GKInstance>();
        String hprdId = null;
        String uniId = null;
        GKInstance localGK = null;
        GKInstance dbGK = null;
        int c = 0;
        for (Iterator<String> it = hprd2LocalGKMap.keySet().iterator(); it.hasNext();) {
            hprdId = it.next();
            uniId = hprd2UniMap.get(hprdId);
            if (uniId == null) {
                System.out.println("HPRD:" + hprdId + " cannot be mapped to UniProt!");
                continue;
            }
            localGK = hprd2LocalGKMap.get(hprdId);
            // Query database for UniProt Identifier
            Collection list = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferencePeptideSequence,
                                                                 ReactomeJavaConstants.identifier,
                                                                 "=",
                                                                 uniId);
            if (list == null || list.size() == 0) {
                c ++;
                continue;
            }
            dbGK = (GKInstance) list.iterator().next();
            local2DB.put(localGK, dbGK);
        }
        System.out.println("HPRD cannot be mapped: " + c);
        return local2DB;
    }
    
    private Map<String, GKInstance> mapHPRD2LocalGK(XMLFileAdaptor fileAdaptor) throws Exception {
        Map<String, GKInstance> hprd2Local = new HashMap<String, GKInstance>();
        Collection list = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferencePeptideSequence);
        GKInstance instance = null;
        for (Iterator it = list.iterator(); it.hasNext();) {
            instance = (GKInstance) it.next();
            // Find HPRD for this instance
            List crossReference = instance.getAttributeValuesList(ReactomeJavaConstants.crossReference);
            if (crossReference == null || crossReference.size() == 0)
                continue;
            for (Iterator it1 = crossReference.iterator(); it1.hasNext();) {
                GKInstance cr = (GKInstance) it1.next();
                // DisplayName might not be correct yet
                GKInstance referenceDB = (GKInstance) cr.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
                if (referenceDB == null)
                    continue;
                String name = (String) referenceDB.getAttributeValue(ReactomeJavaConstants.name);
                if (name.equals("HPRD")) {
                    // Get id
                    String id = (String) cr.getAttributeValue(ReactomeJavaConstants.identifier);
                    hprd2Local.put(id, instance);
                }
            }
        }
        System.out.println("Found HPRD Xref: " + hprd2Local.size());
        return hprd2Local;
    }

    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor, 
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        processEntityCompartment(ReactomeJavaConstants.accession,
                                 dbAdaptor, 
                                 fileAdaptor);
    }
}
