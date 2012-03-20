/*
 * Created on Feb 1, 2007
 *
 */
package org.reactome.psi;

import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.HPRDAnalyzer;

public class HPRDPsiMiToReactomePostProcessor extends PsiMiToReactomePostProcessor {
    
    public HPRDPsiMiToReactomePostProcessor() {
        setDataSourceName("HPRD");
        setDataSourceUrl("http://www.hprd.org");
    }

    @Override
    protected void processEWAS(MySQLAdaptor dbAdaptor, XMLFileAdaptor fileAdaptor) throws Exception {
        mapToUniProtRefPepSeqs(fileAdaptor, dbAdaptor);
        super.processEWAS(dbAdaptor, fileAdaptor);
    }
    
    protected void mapToUniProtRefPepSeqs(XMLFileAdaptor fileAdaptor,
                                          MySQLAdaptor dbAdaptor) throws Exception {
        // Load the mapping file
        //FileUtility fu = new FileUtility();
        //String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/HPRD2UniProtViaGeneSymbols.txt";
        //Map<String, String> hprdToUp = fu.importMap(fileName);
        Map<String, String> hprdToUp = new HPRDAnalyzer().loadHprdIdToUniProtMap();
        GKInstance uniProtDb = PostProcessHelper.getUniProtInstance(dbAdaptor, fileAdaptor);
        Collection c = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferencePeptideSequence);
        GKInstance refPepSeq = null;
        for (Iterator it = c.iterator(); it.hasNext();) {
            refPepSeq = (GKInstance) it.next();
            GKInstance db = (GKInstance) refPepSeq.getAttributeValue(ReactomeJavaConstants.referenceDatabase);
            String dbName = (String) db.getAttributeValue(ReactomeJavaConstants.name);
            if (!(dbName.equals("HPRD")))
                continue;
            String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
            String uniProtId = hprdToUp.get(identifier);
            if (uniProtId == null)
                continue;
            replaceByUniProt(refPepSeq,
                             uniProtId,
                             uniProtDb);
        }
    }
    
    protected void replaceByUniProt(GKInstance refPepSeq,
                                    String uniProtId,
                                    GKInstance uniProtDb) throws Exception {
        refPepSeq.setAttributeValue(ReactomeJavaConstants.referenceDatabase,
                                    uniProtDb);
        refPepSeq.setAttributeValue(ReactomeJavaConstants.identifier, 
                                    uniProtId);
    }
}
