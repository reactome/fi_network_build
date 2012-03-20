/*
 * Created on Mar 24, 2009
 *
 */
package org.reactome.psi;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.UniProtAnalyzer;

/**
 * This class is used to process BioGrid PSI-MI post processing.
 * @author wgm
 *
 */
public class BioGridPsiMiToReactomePostProcessor extends HPRDPsiMiToReactomePostProcessor {
    
    public BioGridPsiMiToReactomePostProcessor() {
        setDataSourceName("BioGrid");
        setDataSourceUrl("http://www.thebiogrid.org");
    }

    @Override
    protected void mapToUniProtRefPepSeqs(XMLFileAdaptor fileAdaptor,
                                          MySQLAdaptor dbAdaptor) throws Exception {
        // Load the mapping file
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, Set<String>> mimToUpIds = uniAnalyzer.loadMimToUniProt();
        GKInstance uniProtDb = PostProcessHelper.getUniProtInstance(dbAdaptor, fileAdaptor);
        Collection c = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferencePeptideSequence);
        GKInstance refPepSeq = null;
        for (Iterator it = c.iterator(); it.hasNext();) {
            refPepSeq = (GKInstance) it.next();
            GKInstance db = (GKInstance) refPepSeq.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            String dbName = (String) db.getAttributeValue(ReactomeJavaConstants.name);
            if (!(dbName.equals("MIM")))
                continue;
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            Set<String> upIds = mimToUpIds.get(identifier);
            if (upIds == null || upIds.size() == 0)
                continue;
            if (upIds.size() == 1) {
                String uniProtId = upIds.iterator().next();
                replaceByUniProt(refPepSeq,
                                 uniProtId,
                                 uniProtDb);
            }
            else {
                handleMultipleMapping(refPepSeq,
                                      upIds,
                                      fileAdaptor,
                                      dbAdaptor);
            }
        }
    }
    
    private void handleMultipleMapping(GKInstance refPepSeq,
                                       Set<String> upIds,
                                       XMLFileAdaptor fileAdaptor,
                                       MySQLAdaptor dbAdaptor) throws Exception {
        Collection referrers = refPepSeq.getReferers(ReactomeJavaConstants.referenceEntity);
        if (referrers == null || referrers.size() == 0) {
            // Just delete it
            fileAdaptor.deleteInstance(refPepSeq);
            return;
        }
        // Create a list of ReferencePeptideSequences from upIds
        Set<GKInstance> refPepSeqs = new HashSet<GKInstance>();
        for (String id : upIds) {
            GKInstance refPepSeq1 = PostProcessHelper.queryRefPepSeq(id, fileAdaptor);
            if (refPepSeq1 == null)
                refPepSeq1 = PostProcessHelper.createLocalRefPepSeq(id, dbAdaptor, fileAdaptor);
            refPepSeqs.add(refPepSeq1);
        }
        // Need to convert referrers to a set
        for (Iterator it = referrers.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            PostProcessHelper.switchEWASToSet(ewas, 
                                              refPepSeqs, 
                                              fileAdaptor);
        }
    }
                                 
}
