/*
 * Created on Aug 24, 2006
 *
 */
package org.reactome.psi;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.database.SynchronizationManager;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.PersistenceManager;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.GKSchemaClass;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.util.FileUtility;

/**
 * BIND PSI_MI data have not UniProt mapping, use this class to do this mapping. 
 * If a Entrez Protein can mapping to more than one UniProt identifiers, a DefinedSet
 * will be used.
 * @author guanming
 *
 */
public class BINDPsiMiToReactomePostProcessor extends PsiMiToReactomePostProcessor {
    private Map<String, GKInstance> uniProt2RefPepSeq = new HashMap<String, GKInstance>();
    
    public BINDPsiMiToReactomePostProcessor() {    
        setDataSourceName("BIND");
        setDataSourceUrl("http://www.bind.ca");
    }
    
    private Map<String, Set<String>> loadGI2UniMap() throws IOException {
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

    @Override
    protected void processEWAS(MySQLAdaptor dbAdaptor, 
                               XMLFileAdaptor fileAdaptor) throws Exception {
        Map<String, Set<String>> gi2UniMap = loadGI2UniMap();
        // Basically used to load ReferenceSequences from the database
        Collection collection = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferencePeptideSequence);
        GKInstance gkInstance = null;
        GKInstance dbInstance = null;
        String identifier = null;
        SynchronizationManager manager = SynchronizationManager.getManager();
        PersistenceManager.getManager().setActiveFileAdaptor(fileAdaptor);
        // Counting
        int notMapped = 0;
        int moreThanOneMapped = 0;
        int mapped = 0;
        for (Iterator it = collection.iterator(); it.hasNext();) {
            gkInstance = (GKInstance) it.next();
            identifier = (String) gkInstance.getAttributeValue(ReactomeJavaConstants.identifier);
            if (identifier == null) // Maybe possible for some instances
                continue; 
            Set<String> uniIds = gi2UniMap.get(identifier);
            if (uniIds == null) {
                notMapped ++;
                continue;
            }
            else if (uniIds.size() > 1) {
                moreThanOneMapped ++;
                handleMultipleMap(uniIds, gkInstance, dbAdaptor, fileAdaptor);
            }
            else {
                mapped ++;
                identifier = uniIds.iterator().next();
                handleSingleMap(identifier, gkInstance, dbAdaptor, fileAdaptor);
            }
        }
        System.out.println("Not mapped: " + notMapped);
        System.out.println("Single Mapped: " + mapped);
        System.out.println("More Than One Mapped: " + moreThanOneMapped);
    }
    
    private void handleMultipleMap(Set<String> uniIdSet,
                                   GKInstance localRefPepSeq,
                                   MySQLAdaptor dbAdaptor,
                                   XMLFileAdaptor fileAdaptor) throws Exception {
        // EWAS should be converted as DefinedSet
        Collection referrers = localRefPepSeq.getReferers(ReactomeJavaConstants.referenceEntity);
        if (referrers == null || referrers.size() == 0)
            return; // Don't care if no instances use this
        // Convert to List to use the index
        List<String> uniIdList = new ArrayList<String>(uniIdSet);
        // Use the existing one for the first mapping
        List<GKInstance> refPepSeqList = new ArrayList<GKInstance>();
        String uniId = uniIdList.get(0);
        // Other GI might be mapped to this UniId too even though this GI
        // can be mapped to more than one UniId. So getReferrers() should
        // be called before the following calling.
        handleSingleMap(uniId, localRefPepSeq, dbAdaptor, fileAdaptor);
        // localRefPepSeq might be changed. Use map to fetch it.
        refPepSeqList.add(uniProt2RefPepSeq.get(uniId));
        // Generate enough localRefPepSeq to be used for mapping
        for (int i = 1; i < uniIdList.size(); i++) {
            GKInstance newRefPepSeq = fetchRefPepSeqForUniID(uniIdList.get(i), 
                                                             dbAdaptor, 
                                                             fileAdaptor);
            refPepSeqList.add(newRefPepSeq);
        }
        GKSchemaClass definedSetCls = (GKSchemaClass) fileAdaptor.getSchema().getClassByName(ReactomeJavaConstants.DefinedSet);
        for (Iterator it = referrers.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            // Convert this EWAS to DefinedSet
            fileAdaptor.switchType(ewas, definedSetCls);
            // Genreate EWAS based on RefPepSeq
            for (GKInstance refPepSeq : refPepSeqList) {
                GKInstance newEwas = createEWASFromSet(ewas, fileAdaptor);
                newEwas.setAttributeValue(ReactomeJavaConstants.referenceEntity,
                                          refPepSeq);
            }
        }
    }
    
    @SuppressWarnings("unchecked")
    private GKInstance createEWASFromSet(GKInstance set,
                                         XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance newEwas = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        // Copy properties
        List names = set.getAttributeValuesList(ReactomeJavaConstants.name);
        if (names != null)
            newEwas.setAttributeValueNoCheck(ReactomeJavaConstants.name, 
                                             new ArrayList(names));
        set.addAttributeValue(ReactomeJavaConstants.hasMember, newEwas);
        return newEwas;
    }
    
    private GKInstance fetchRefPepSeqForUniID(String uniProtId,
                                              MySQLAdaptor dbAdaptor,
                                              XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance rtn = uniProt2RefPepSeq.get(uniProtId);
        if (rtn != null)
            return rtn;
        GKInstance dbInstance = fetchReferenceSequence(uniProtId, dbAdaptor);
        if (dbInstance != null) {
            rtn = downloadAsShell(dbInstance, fileAdaptor);
            uniProt2RefPepSeq.put(uniProtId, rtn);
        } 
        else {
            rtn = fileAdaptor.createNewInstance(ReactomeJavaConstants.ReferencePeptideSequence);
            GKInstance uniProtDB = PostProcessHelper.getUniProtInstance(dbAdaptor, fileAdaptor);
            rtn.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                  uniProtDB);
            rtn.setAttributeValue(ReactomeJavaConstants.identifier,
                                  uniProtId);
            uniProt2RefPepSeq.put(uniProtId, rtn);
        }
        return rtn;
    }
    
    private GKInstance downloadAsShell(GKInstance dbInstance,
                                       XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance localShell = fileAdaptor.createNewInstance(dbInstance.getSchemClass().getName());
        localShell.setIsDirty(false);
        localShell.setIsShell(true);
        localShell.setDBID(dbInstance.getDBID());
        localShell.setDisplayName(dbInstance.getDisplayName());
        return localShell;
    }
    
    private void handleSingleMap(String uniProtId, 
                                 GKInstance localRefPepSeq,
                                 MySQLAdaptor dbAdaptor,
                                 XMLFileAdaptor fileAdaptor) throws Exception {
        // More than one GIID might be mapped to the same UniProt. In this case,
        // LocalRefPepSeq should be merged.
        GKInstance localMapped = uniProt2RefPepSeq.get(uniProtId);
        if (localMapped != null) {
            // Merge localRefPepSeq to localMapped
            merge(localMapped, localRefPepSeq, fileAdaptor);
        }
        else {
            GKInstance dbInstance = fetchReferenceSequence(uniProtId, dbAdaptor);
            if (dbInstance != null) {
                useDBInstanceAsShell(dbInstance, localRefPepSeq, fileAdaptor);
                uniProt2RefPepSeq.put(uniProtId, localRefPepSeq);
            }
            else {// This will force two or more mapped by the same GI id to be merged.
                // Want to use UniProt
                GKInstance uniProtDB = PostProcessHelper.getUniProtInstance(dbAdaptor, fileAdaptor);
                localRefPepSeq.setAttributeValue(ReactomeJavaConstants.referenceDatabase, 
                                                 uniProtDB);
                localRefPepSeq.setAttributeValue(ReactomeJavaConstants.identifier, 
                                                 uniProtId);
                uniProt2RefPepSeq.put(uniProtId, localRefPepSeq); 
            }
        }
    }
    
    @SuppressWarnings("unchecked")
    private void merge(GKInstance targetRefPepSeq,
                       GKInstance sourceRefPepSeq,
                       XMLFileAdaptor fileAdaptor) throws Exception {
        PostProcessHelper.mergeReferenceEntity(targetRefPepSeq, sourceRefPepSeq, fileAdaptor);
    }
}
