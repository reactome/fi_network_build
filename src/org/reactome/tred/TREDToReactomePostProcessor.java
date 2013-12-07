/*
 * Created on Apr 15, 2009
 *
 */
package org.reactome.tred;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.PostProcessTemplate;


/**
 * To handle the post-process for converted CuratorTool project.
 * @author wgm
 *
 */
public class TREDToReactomePostProcessor extends PostProcessTemplate {
    private static final long HOMO_SAPIENS_ID = 48887L;

    public TREDToReactomePostProcessor() {
    }
    
    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor,
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        super.attachDataSource("TRED", 
                               "http://rulai.cshl.edu/cgi-bin/TRED/tred.cgi?process=home",
                               dbAdaptor, 
                               fileAdaptor);
    }
    
    @Override
    @SuppressWarnings("unchecked")
    protected void processEWAS(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        Set<String> unmapped = new HashSet<String>();
        // Preload these instances to avoid slow SQL queries for quick performance
        GKInstance human = dbAdaptor.fetchInstance(HOMO_SAPIENS_ID);
        Collection<GKInstance> refGeneProducts = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct, 
                                                                                    ReactomeJavaConstants.species,
                                                                                    "=",
                                                                                    human);
        dbAdaptor.loadInstanceAttributeValues(refGeneProducts, new String[]{ReactomeJavaConstants.geneName,
                                                                            ReactomeJavaConstants.identifier});
        // Create a map for quick mapping
        Map<String, GKInstance> geneNameToRefGeneProd = new HashMap<String, GKInstance>();
        for (GKInstance refGeneProd : refGeneProducts) {
            List<String> geneNames = refGeneProd.getAttributeValuesList(ReactomeJavaConstants.geneName);
            for (String geneName : geneNames) {
                geneNameToRefGeneProd.put(geneName, refGeneProd);
            }
        }
        // Fetch RefPepSeq from the database and link it to EWAS based on gene names.
        Collection<GKInstance> ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        for (GKInstance ewas : ewases) {
            String geneName = ewas.getDisplayName();
            // Just do a simple search
            GKInstance refGeneProd = geneNameToRefGeneProd.get(geneName);
            if (refGeneProd == null) {
                unmapped.add(geneName);
                logger.warn(geneName + " for EWAS " + ewas.getDBID() + " cannot be mapped to UniProt!");
                continue;
            }
            //TODO: Currently we just pick up one ReferenceGeneProduct only though it is known that a gene name
            // can be mapped to multiple UniProt accessions!!!
            // Want to get ReferenceGeneProduct only instead of its isoforms
            String uniProtId = (String) refGeneProd.getAttributeValue(ReactomeJavaConstants.identifier);
            GKInstance refPepSeq = PostProcessHelper.getRefPepSeq(uniProtId, 
                                                                  dbAdaptor, 
                                                                  fileAdaptor);
            ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity, 
                                   refPepSeq);
        }
        logger.info("Total unmapped names: " + unmapped.size());
    }
    
    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor,
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        // There is no compartment information.
    }
    
}
