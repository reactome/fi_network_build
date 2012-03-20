/*
 * Created on Apr 15, 2009
 *
 */
package org.reactome.tred;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
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
    protected void processEWAS(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        Set<String> unmapped = new HashSet<String>();
        // Fetch RefPepSeq from the database and link it to EWAS based on gene names.
        Collection<?> ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        for (Iterator<?> it = ewases.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            String geneName = ewas.getDisplayName();
            Collection<?> refGeneProds = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                            ReactomeJavaConstants.geneName, 
                                                                            "=",
                                                                            geneName);
            if (refGeneProds == null || refGeneProds.size() == 0) {
                unmapped.add(geneName);
                logger.warn(geneName + " for EWAS " + ewas.getDBID() + " cannot be mapped to UniProt!");
                continue;
            }
            GKInstance refGeneProd = null;
            // Have to make sure GKInstance should be from human only. There are many different species in the database now
            for (Object obj : refGeneProds) {
                GKInstance refGeneProd1 = (GKInstance) obj;
                GKInstance species = (GKInstance) refGeneProd1.getAttributeValue(ReactomeJavaConstants.species);
                if (species == null || !species.getDBID().equals(48887L)) // 48887 is  homo sapiens id
                    continue;
                refGeneProd = refGeneProd1;
                break;
            }
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
        System.out.println("Total unmapped names: " + unmapped.size());
    }
    
    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor,
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        // There is no compartment information.
    }
    
}
