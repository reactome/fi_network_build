/*
 * Created on Jan 31, 2007
 *
 */
package org.reactome.kegg;

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
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.convert.common.PostProcessTemplate;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.InteractionUtilities;

public class KEGGPostProcess extends PostProcessTemplate {
    private static Logger logger = Logger.getLogger(KEGGPostProcess.class);
    
    
    @Override
    public void postProcess(MySQLAdaptor dbAdaptor, XMLFileAdaptor fileAdaptor) throws Exception {
        super.postProcess(dbAdaptor, fileAdaptor);
        logger.info("Starting attach species...");
        attachSpecies(dbAdaptor,
                      fileAdaptor);
    }
    
    @Override
    protected void setDisplayNames(XMLFileAdaptor fileAdaptor) throws Exception {
        // Call this method first
        super.setDisplayNames(fileAdaptor);
        Collection complexes = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Complex);
        StringBuilder builder = new StringBuilder();
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            GKInstance complex = (GKInstance) it.next();
            List components = complex.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
            if (components == null || components.size() == 0)
                continue;
            builder.setLength(0);
            for (Iterator it1 = components.iterator(); it1.hasNext();) {
                GKInstance comp = (GKInstance) it1.next();
                builder.append(comp.getDisplayName());
                if (it1.hasNext())
                    builder.append("+");
            }
            complex.setDisplayName(builder.toString());
        }
        // Reset interactions
        Collection interactions = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        for (Iterator it = interactions.iterator(); it.hasNext();) {
            GKInstance inter = (GKInstance) it.next();
            List interactors = inter.getAttributeValuesList(ReactomeJavaConstants.interactor);
            builder.setLength(0);
            for (Iterator it1 = interactors.iterator(); it1.hasNext();) {
                GKInstance interactor = (GKInstance) it1.next();
                builder.append(interactor.getDisplayName());
                if (it1.hasNext())
                    builder.append("-");
            }
            inter.setDisplayName(builder.toString());
        }
    }

    /**
     * Use this method to attach human as species to Complex, EWAS, Interaction and Pathways.
     * @param dbAdaptor
     * @param fileAdaptor
     * @throws Exception
     */
    private void attachSpecies(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        GKInstance human = PostProcessHelper.getHumanInstance(dbAdaptor, fileAdaptor);
        Collection events = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Event);
        for (Iterator it = events.iterator(); it.hasNext();) {
            GKInstance instance = (GKInstance) it.next();
            instance.setAttributeValue(ReactomeJavaConstants.species,
                                       human);
        }
        Collection ewass = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        for (Iterator it = ewass.iterator(); it.hasNext();) {
            GKInstance instance = (GKInstance) it.next();
            instance.setAttributeValue(ReactomeJavaConstants.species,
                                       human);
        }
        Collection complexes = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Complex);
        for (Iterator it = complexes.iterator(); it.hasNext();) {
            GKInstance complex = (GKInstance) it.next();
            complex.setAttributeValue(ReactomeJavaConstants.species,
                                      human);
        }
    }

    @Override
    protected void attachDataSource(MySQLAdaptor dbAdaptor,
                                    XMLFileAdaptor fileAdaptor) throws Exception {
        String dbName = "KEGG";
        String url = "http://www.genome.jp/kegg/";
        attachDataSource(dbName, 
                         url, 
                         dbAdaptor, 
                         fileAdaptor);
    }
    
    /**
     * Find ReferenceGeneProducts based on crossReferences. Several steps have been applied following:
     * 1). A converted KEGG EWAS should have a cross-reference that points to a KEGG id. Get this KEGG ID
     * 2). One KEGG id can be mapped to more than one UniProt id. Create or map to ReferenceGeneProduct for
     * each UniProt id.
     * 3). If an EWAS can be mapped to more than one UniProt based on 2), convert it to DefinedSet.
     */
    @Override
    @SuppressWarnings("unchecked")
    protected void processEWAS(MySQLAdaptor dbAdaptor,
                               XMLFileAdaptor fileAdaptor) throws Exception {
        Collection<?> ewases = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
        // Get all KEGG ids: one KEGG id may be used in multple EWASes
        Map<String, Set<GKInstance>> keggIdToEWASes = new HashMap<String, Set<GKInstance>>();
        for (Object obj : ewases) {
            GKInstance ewas = (GKInstance) obj;
            GKInstance crossRef = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.crossReference);
            if (crossRef == null)
                continue;
            String keggId = (String) crossRef.getAttributeValue(ReactomeJavaConstants.identifier);
            if (keggId == null)
                continue;
            InteractionUtilities.addElementToSet(keggIdToEWASes, keggId, ewas);
        }
        // Map kegg ids to ReferenceGeneProducts
        KeggAnalyzer analyzer = new KeggAnalyzer();
        Map<String, Set<String>> keggIdToUniProtIMap = analyzer.loadKeggToAllUniProtIdsMap();
        for (String keggId : keggIdToEWASes.keySet()) {
            Set<String> uniProtIds = keggIdToUniProtIMap.get(keggId);
            // The follow check has been used for 2009 version
//            if (uniProtIds == null || uniProtIds.size() == 0) {
//                String uniProtId = mapKeggIdToUniProtViaGeneId(keggId);
//                if (uniProtId != null) {
//                    uniProtIds = new HashSet<String>();
//                    uniProtIds.add(uniProtId);
//                }
//            }
            if (uniProtIds == null || uniProtIds.size() == 0){
                logger.warn("KEGG ID cannot be mapped to UniProt: " + keggId);
                continue; // Don't do anything here since there are not many unmapped kegg ids.
            }
            Set<GKInstance> refGeneProds = new HashSet<GKInstance>();
            for (String uniProtId : uniProtIds) {
                GKInstance refGeneProd = PostProcessHelper.getRefPepSeq(uniProtId, dbAdaptor, fileAdaptor);
                if (refGeneProd == null) {
                    logger.error("Cannot get ReferenceGeneProduct: " + uniProtId);
                    continue;
                }
                refGeneProds.add(refGeneProd);
            }
            if (refGeneProds.size() == 0) {
                logger.error("Cannot get ReferenceGeneProduct for KEGG id: " + keggId);
                continue;
            }
            Set<GKInstance> ewasSet = keggIdToEWASes.get(keggId);
            for (GKInstance ewas : ewasSet) {
                if (refGeneProds.size() == 1) {
                    ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity, 
                                           refGeneProds.iterator().next());
                    addNameToEWAS(ewas);
                }
                else {
                    PostProcessHelper.switchEWASToSet(ewas, refGeneProds, fileAdaptor);
                    // Make sure names are not empty in the members
                    // ewas is a DefinedSet now
                    List<GKInstance> members = (List<GKInstance>) ewas.getAttributeValuesList(ReactomeJavaConstants.hasMember);
                    for (GKInstance member : members) 
                        addNameToEWAS(member);
                }
            }
        }
        
//        KeggAnalyzer analyzer = new KeggAnalyzer();
//        Map<String, String> keggIdToUniprot = analyzer.loadKeggToUniProtMaps();
////        Map<String, Set<String>> keggIdToUniProtIds = analyzer.loadKeggToAllUniProtIdsMap();
//        Map<String, GKInstance> keggIdToRefPepSeq = new HashMap<String, GKInstance>();
//        Collection collection = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.EntityWithAccessionedSequence);
//        for (Iterator it = collection.iterator(); it.hasNext();) {
//            GKInstance ewas = (GKInstance) it.next();
//            GKInstance keggRef = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.crossReference);
//            if (keggRef == null)
//                continue; // Cannot map
//            Set<GKInstance> refPepSeqs = new HashSet<GKInstance>();
//            GKInstance refPepSeq = mapKeggRef(keggRef,
//                                              keggIdToUniprot,
//                                              keggIdToRefPepSeq,
//                                              fileAdaptor,
//                                              dbAdaptor);
//            if (refPepSeq == null) {
//                logger.warn(keggRef + " cannot be mapped to uniprot!");
//                continue; // Maybe cannot be mapped
//            }
//            ewas.setAttributeValue(ReactomeJavaConstants.referenceEntity, 
//                                   refPepSeq);
//            // Want to have names
//            if (ewas.getAttributeValue(ReactomeJavaConstants.name) == null) {
//                List names = refPepSeq.getAttributeValuesList(ReactomeJavaConstants.geneName);
//                if (names != null && names.size() > 0)
//                    ewas.setAttributeValueNoCheck(ReactomeJavaConstants.name, 
//                                                  new ArrayList(names));
//                else {
//                    // Just use identifier to avoid unknown names
//                    String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
//                    ewas.addAttributeValue(ReactomeJavaConstants.name,
//                                           identifier);
//                }
//            }
//        }
    }
    
    @SuppressWarnings("unchecked")
    private void addNameToEWAS(GKInstance ewas) throws Exception {
        // Want to have names
        if (ewas.getAttributeValue(ReactomeJavaConstants.name) == null) {
            GKInstance refPepSeq = (GKInstance) ewas.getAttributeValue(ReactomeJavaConstants.referenceEntity);
            if (refPepSeq == null)
                return;
            List<?> names = refPepSeq.getAttributeValuesList(ReactomeJavaConstants.geneName);
            if (names != null && names.size() > 0)
                ewas.setAttributeValueNoCheck(ReactomeJavaConstants.name, 
                                              new ArrayList(names));
            else {
                // Just use identifier to avoid unknown names
                String identifier = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                ewas.addAttributeValue(ReactomeJavaConstants.name,
                                       identifier);
            }
        }
    }
    
//    private GKInstance mapKeggRef(GKInstance keggRef,
//                                  Map<String, String> keggIdToUMap,
//                                  Map<String, GKInstance> keggIdToUPRefPMap,
//                                  XMLFileAdaptor fileAdaptor,
//                                  MySQLAdaptor dbAdaptor) throws Exception {
//        String keggId = (String) keggRef.getAttributeValue(ReactomeJavaConstants.identifier);
//        if (keggIdToUPRefPMap.containsKey(keggId)) {
//            return keggIdToUPRefPMap.get(keggId); // Might be null
//        }
//        GKInstance refPepSeq = null;
//        String uniProtId = keggIdToUMap.get(keggId);
//        if (uniProtId == null) {
//            // On the KEGG download in January, 2012, there are 26 KEGG ids cannot be mapped
//            // to UniProts. Check to ID mapping in UniProt, find three mappings. 
//            // So ignore these ids. The unmapped ids are listed in file UnMappedKeggIds.txt.
//            // Try to get map via GeneId
//            // Mapping map may be created. However, only one is used.
//            logger.info("Cannot map to UniProt, try via Entrez: " + keggId);
//            uniProtId = mapKeggIdToUniProtViaGeneId(keggId);
//        }
//        if (uniProtId != null) 
//            // Have to create a new one
//            refPepSeq = PostProcessHelper.getRefPepSeq(uniProtId, dbAdaptor, fileAdaptor);
//        keggIdToUPRefPMap.put(uniProtId,
//                              refPepSeq); // refPepSeq might be null
//        return refPepSeq;
//    }
    
    private String mapKeggIdToUniProtViaGeneId(String keggId) throws IOException {
        String geneId = keggId.substring(4);
        String fileName = FIConfiguration.getConfiguration().get("ENTREZ_TO_UNIPROT_MAP_FILE_NAME");
        return PostProcessHelper.mapEntrezToUniProt(fileName, geneId);
    }
    
    @Override
    protected void processEntityCompartment(MySQLAdaptor dbAdaptor,
                                            XMLFileAdaptor fileAdaptor) throws Exception {
        // No compartment information is available in the KEGG database.
        // Do nothing here.
    }
    
}
