/*
 * Created on Jan 16, 2012
 *
 */
package org.reactome.b2rPostProcessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;

/**
 * @author gwu
 *
 */
public class NciPIDBToRPostProcessor extends HPRDBToRPostProcessor {
    public Logger logger = Logger.getLogger(NciPIDBToRPostProcessor.class);
    private String dataSourceName  = "Pathway Interaction Database";
    
    public NciPIDBToRPostProcessor() {
    }
    
    @Override
    public void postProcess(MySQLAdaptor dbAdaptor, XMLFileAdaptor fileAdaptor) throws Exception {
        super.postProcess(dbAdaptor, fileAdaptor);
        generatePrecedingEventProperties(fileAdaptor);
    }
    
    private void generatePrecedingEventProperties(XMLFileAdaptor fileAdaptor) throws Exception {
        Collection pathways = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Pathway);
        for (Iterator it = pathways.iterator(); it.hasNext();) {
            // Here treat all pathways are top level pathways
            GKInstance topPathway = (GKInstance) it.next();
            PostProcessHelper.generatePrecedingProperties(topPathway);
        }
    }
    
    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor, 
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        String dbName = getDataSourceName();
        String url = "http://pid.nci.nih.gov/";
        attachDataSource(dbName, url, dbAdaptor, fileAdaptor);
    }
    
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    
    public String getDataSourceName() {
        return this.dataSourceName;
    }
    
    @Override
    protected void processReferencePeptideSequences(MySQLAdaptor dbAdaptor, 
                                                    XMLFileAdaptor fileAdaptor) throws IOException, Exception {
        logger.info("Starting process ReferencePeptideSequence...");
        //  Map<GKInstance, GKInstance> local2DBMap = mapLocalToDB(fileAdaptor,
        //                                                         dbAdaptor);
        //  copyProperties(local2DBMap, fileAdaptor);
        processEWASBasedOnRefPepSeqs(fileAdaptor, dbAdaptor);
        cleanUpRefPepSeqs(fileAdaptor);
    }
    
    /**
     * In NCI-Nature pathways, some EWASs are actually protein families that refer to multiple
     * UniProt identifiers. In case like this, EWAS should be converted into DefinedSet, and link to
     * individual RefSeqSeqs. In this method implementation, all RefPepSeqs except that cannot be mapped
     * to UniProt have been regenerated to make the work simplier.
     * @throws Exception
     */
    private void processEWASBasedOnRefPepSeqs(XMLFileAdaptor fileAdaptor,
                                              MySQLAdaptor dbAdaptor) throws Exception {
        // Get the map from the original RefPepSeq to new RefPepSeq
        Map<GKInstance, Set<GKInstance>> originalRefPepToNewRefPepSeqs = remapRefPepSeqs(fileAdaptor, dbAdaptor);
        Collection ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        if (ewases == null || ewases.size() == 0)
            return; // Nothing needed to be done
        // If an EWAS can be mapped to more than one UniProt ids, then this EWAS should be remapped.
        for (Iterator it = ewases.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            //      if (ewas.getDisplayName().equals("pid_x_100031") ||
            //          ewas.getDisplayName().equals("pid_x_100031_26")) {
            //          System.out.println(ewas);
            //          GKInstance modifiedResidue = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.hasModifiedResidue);
            //          System.out.println("has modifiedResidue: " + modifiedResidue);
            //      }
            GKInstance refPepSeq = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            Set<GKInstance> mappedRefPepSeqs = originalRefPepToNewRefPepSeqs.get(refPepSeq);
            if (mappedRefPepSeqs == null)
                continue; // Don't do anything. refPepSeq may be other type of ReferenceSequence (e.g. ReferenceRNASequence).
            // If there is only one mappedRefPepSeqs. Just use it
            if (mappedRefPepSeqs.size() == 1) {
                GKInstance mapped = mappedRefPepSeqs.iterator().next();
                ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity, mapped);
            }
            else {
                // Need to create new DefinedSet
                PostProcessHelper.switchEWASToSet(ewas, 
                                                  mappedRefPepSeqs,
                                                  fileAdaptor);
            }
        }
        // Copy necessary attributes (esp. name) from RefPepSeq to EWAS
        // Refetch this list
        ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        for (Iterator it = ewases.iterator(); it.hasNext();) {
            GKInstance ewas = (GKInstance) it.next();
            // check if species is there. If not, try to copy from the RefPepSeq plus name
            GKInstance species = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.species);
            if (species == null) {
                GKInstance referenceEntity = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
                if (referenceEntity != null) {
                    copyAttributeValue(ReactomeJavaConstants.species,
                                       referenceEntity,
                                       ewas);
                    // Want to keep the original name as the reference to the old ids. But this original name
                    // should be listed at the end of the name list
                    List ewasNames = ewas.getAttributeValuesList(ReactomeJavaConstants.name);
                    List refEntityNames = referenceEntity.getAttributeValuesList(ReactomeJavaConstants.name);
                    if (refEntityNames != null && refEntityNames.size() > 0) {
                        // Make a copy
                        List<String> nameCopy = new ArrayList<String>(refEntityNames);
                        for (int i = 0;i < ewasNames.size(); i++) {
                            String name = (String) ewasNames.get(i);
                            if (nameCopy.contains(name))
                                continue;
                            nameCopy.add(name);
                        }
                        ewas.setAttributeValue(ReactomeJavaConstants.name, nameCopy);
                    }
                }
                InstanceDisplayNameGenerator.setDisplayName(ewas);
            }
        }
    }
    
    /**
     * Re-map RefPepSeq instances based on crossReference values used in the first round converting.
     * Some of RefPepSeqs can be mapped to multiple UniProt ids. In those cases, an original RefPepSeq
     * can be mapped to multiple new RefPepSeqs.
     * @param fileAdpator
     * @param dbAdaptor
     * @return
     * @throws Exception
     */
    private Map<GKInstance, Set<GKInstance>> remapRefPepSeqs(XMLFileAdaptor fileAdaptor,
                                                             MySQLAdaptor dbAdaptor) throws Exception {
        Collection refPepSeqs = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.ReferenceGeneProduct);
        // Used to map UniProt ids to the latest version
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        Map<String, String> idMap = uniProtAnalyzer.loadUniProtIDsMap();
        Map<GKInstance, Set<GKInstance>> originalRefPepToNewRefPeps = new HashMap<GKInstance, Set<GKInstance>>();
        int totalEntrezCases = 0;
        for (Iterator it = refPepSeqs.iterator(); it.hasNext();) {
            GKInstance refPepSeq = (GKInstance) it.next();
            Set<GKInstance> set = new HashSet<GKInstance>();
            originalRefPepToNewRefPeps.put(refPepSeq, set);
            List crossRefs = refPepSeq.getAttributeValuesList(ReactomeJavaConstants.crossReference);
            if (crossRefs == null || crossRefs.size() == 0) {
                set.add(refPepSeq);
                continue;
            }
            // The simple case
            for (int i = 0; i < crossRefs.size(); i++) {
                GKInstance crossRef = (GKInstance) crossRefs.get(i);
                String identifier = (String) crossRef.getAttributeValue(ReactomeJavaConstants.identifier);
                String uniProtId = null;
                if (crossRef.getDisplayName().startsWith("EntrezGene")) {
                    uniProtId = mapEntrezToUniProt(identifier);
                    totalEntrezCases ++;
                }
                else if (crossRef.getDisplayName().startsWith("UniProt"))
                    uniProtId = identifier;
                if (uniProtId == null) {
                    set.add(refPepSeq); // This may be added multiple times.
                }
                else {
                    // Just in case not the latest UniProt id is used.
                    String nonRedundant = idMap.get(uniProtId);
                    if (nonRedundant == null)
                        nonRedundant = uniProtId; 
                    GKInstance mapped = PostProcessHelper.getRefPepSeq(nonRedundant,
                                                                       dbAdaptor, 
                                                                       fileAdaptor);
                    set.add(mapped);
                }
            }
        }
        System.out.println("Total Entrez cases: " + totalEntrezCases);
        return originalRefPepToNewRefPeps;
    }
    
    private String mapEntrezToUniProt(String entrez) throws IOException {
        return PostProcessHelper.mapEntrezToUniProt(FIConfiguration.getConfiguration().get("ENTREZ_TO_UNIPROT_MAP_FILE_NAME"), entrez);
    }
    
}
