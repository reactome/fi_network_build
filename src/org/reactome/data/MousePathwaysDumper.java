package org.reactome.data;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.DiagramGKBReader;
import org.gk.persistence.DiagramGKBWriter;
import org.gk.persistence.MySQLAdaptor;
import org.gk.render.Renderable;
import org.gk.render.RenderablePathway;
import org.gk.schema.SchemaAttribute;
import org.gk.slicing.ProjectBasedSlicingEngine;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to migrate mouse pathways from a release database to reactome_plus_i database so that we can do pathway analysis
 * for mouse. The mouse pathways are predicted based on the orthologous script. The methods used for mirgration is based on slicing tool.
 * Note: To avoid DB_IDs conflicts, DB_IDs for some mouse related instances have been changed.
 * @author wug
 *
 */
@SuppressWarnings("unchecked")
public class MousePathwaysDumper extends ProjectBasedSlicingEngine {
    private final Long MOUSE_DB_ID = 48892L;
    private static final Logger logger = Logger.getLogger(MousePathwaysDumper.class);
    private MySQLAdaptor sourceDBA;
    private MySQLAdaptor targetDBA;
    private Map<GKInstance, RenderablePathway> pdToDiagram;
    
    public MousePathwaysDumper() {
    }
    
    private MySQLAdaptor getSourceDBA() throws Exception {
        if (sourceDBA != null)
            return sourceDBA;
        sourceDBA = new MySQLAdaptor("localhost",
                               FIConfiguration.getConfiguration().get("MOUSE_SOURCE_DB_NAME"), 
                               FIConfiguration.getConfiguration().get("DB_USER"),
                               FIConfiguration.getConfiguration().get("DB_PWD"));
        logger.info("Initializing sourceDBA: " + sourceDBA.getDBName());
        return sourceDBA;
    }
    
    private MySQLAdaptor getTargetDBA() throws Exception {
        if (targetDBA != null)
            return targetDBA;
        targetDBA = new MySQLAdaptor("localhost",
                                     FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"), 
                                     FIConfiguration.getConfiguration().get("DB_USER"),
                                     FIConfiguration.getConfiguration().get("DB_PWD"));
        logger.info("Initializing targetDBA: " + targetDBA.getDBName());
        return targetDBA;
    }
    
    /**
     * Basically this is a simplified slice method. However, DB_IDs of some mouse instances have been reset to avoid conflicts.
     * @throws Exception
     */
    @Test
    public void dump() throws Exception {
        // Need to set up the source DBA first
        setSource(getSourceDBA());
        eventMap = extractEvents();
        extractReferences();
        extractRegulations();
        // To be used by pathway digrams
        Map<Long, GKInstance> topHuman2mouse = getTopLevelHumanToMouseMap();
        extractPathwayDiagrams();
        removeTargetInstances();
        if(_dumpInstances()) {
            updateHumanTopLevelPathways(topHuman2mouse);
            fixMouseGeneNames();
        }
    }
    
    /**
     * This most likely is a one-time thing. Previously we don't have the geneName slot filled for mouse
     * ReferenceGeneProduct instances. Since version 72 (?), we have done this. This method is to use
     * version 73 data to fix version 71 problem.
     * @throws Exception
     */
    @Test
    public void fixMouseGeneNames() throws Exception {
        // We cannot use a release database with gene names filled since some of mouse genes
        // may not be in the current database.
        logger.info("Loading id to gene name mapping from an external file...");
        // Note: The following file used contains both reviewed and unreviewed entries because 
        // there are many unreviewed entries have been used in our inference.
        String fileName = "datasets/UniProt/mouse_2020_09/uniprot-proteome-mouse.tab";
        Map<String, List<String>> idToGenes = new HashMap<>();
        try (Stream<String> lines = Files.lines(Paths.get(fileName))) {
            lines.forEach(line -> {
                String[] tokens = line.split("\t");
                String[] geneNames = tokens[4].split(" ");
                idToGenes.put(tokens[0], Stream.of(geneNames).collect(Collectors.toList()));
            });
        }
        logger.info("Total mapping size: " + idToGenes.size());
        // Load all mouse ReferenceGeneProduct instances
        MySQLAdaptor targetDBA = getTargetDBA();
        GKInstance mouse = targetDBA.fetchInstance(MOUSE_DB_ID);
        Collection<GKInstance> mouseRefGeneProduct = targetDBA.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                 ReactomeJavaConstants.species,
                                                                 "=",
                                                                 mouse);
        targetDBA.loadInstanceAttributeValues(mouseRefGeneProduct,
                                              new String[] {ReactomeJavaConstants.identifier});
        logger.info("Mouse ReferenceGeneProducts in the target database: " + mouseRefGeneProduct.size());
        logger.info("Filling gene names...");
        // Try to use transaction
        boolean isTnSupported = targetDBA.supportsTransactions();
        if (isTnSupported)
            targetDBA.startTransaction();
        try {
            for (GKInstance inst : mouseRefGeneProduct) {
                String identifier = (String) inst.getAttributeValue(ReactomeJavaConstants.identifier);
                List<String> geneNames = idToGenes.get(identifier);
                if (geneNames == null) {
                    logger.warn("Cannot find gene names for " + identifier);
                    continue;
                }
                inst.setAttributeValue(ReactomeJavaConstants.geneName,
                                       geneNames);
                targetDBA.updateInstanceAttribute(inst, ReactomeJavaConstants.geneName);
            }
            if (isTnSupported)
                targetDBA.commit();
        }
        catch (Exception e) {
            if (isTnSupported)
                targetDBA.rollback();
            logger.error("SlicingEngine.fixMouseGeneNames(): " + e, e);
        }
        logger.info("Done fixMouseGeneNames.");
    }
    
    protected GKInstance fetchDiagramForPathway(GKInstance pathway,
                                                MySQLAdaptor dba) throws Exception {
        Collection collection = dba.fetchInstanceByAttribute(ReactomeJavaConstants.PathwayDiagram,
                                                             ReactomeJavaConstants.representedPathway,
                                                             "=",
                                                             pathway);    
        if (collection == null || collection.size() == 0)
            return null;
        // Pick up the first one only
        return (GKInstance) collection.iterator().next();
    }
    
    @Override
    protected void extractPathwayDiagrams() throws Exception {
        MySQLAdaptor sourceDBA = getSourceDBA();
        Collection<GKInstance> pathwayDiagrams = sourceDBA.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram);
        sourceDBA.loadInstanceAttributeValues(pathwayDiagrams, new String[]{ReactomeJavaConstants.representedPathway});
        List<GKInstance> mousePDs = new ArrayList<>();
        for (GKInstance pd : pathwayDiagrams) {
            List<GKInstance> repPathways = pd.getAttributeValuesList(ReactomeJavaConstants.representedPathway);
            for (GKInstance repPathway : repPathways) {
                if (eventMap.containsKey(repPathway.getDBID())) {
                    mousePDs.add(pd);
                    break;
                }
            }
        }
        logger.info("Total mouse pathway diagrams: " + mousePDs.size());
        DiagramGKBReader diagramReader = new DiagramGKBReader();
        pdToDiagram = new HashMap<>();
        for (GKInstance pd : mousePDs) {
            // Add this diagram to the slice map
            extractReferencesToInstance(pd);
            RenderablePathway diagram = diagramReader.openDiagram(pd);
            for (Renderable r : (List<Renderable>) diagram.getComponents()) {
                if (r.getReactomeId() == null)
                    continue;
                GKInstance inst = sourceDBA.fetchInstance(r.getReactomeId());
                r.setInstance(inst);
            }
            pdToDiagram.put(pd, diagram);
        }
        logger.info("extractPathwayDiagrams(): " + sliceMap.size());
    }
    
    @Override
    protected void extractReferences() throws Exception {
        // Check all references in the events
        logger.info("Starting extractReference...");
        for (Long dbId : (Collection<Long>) eventMap.keySet()) {
            GKInstance instance = (GKInstance) eventMap.get(dbId);
            extractReferencesToInstance(instance);
        }
        logger.info("extractReferences(): " + sliceMap.size() + " instances");
    }
    
    @Override
    protected void extractReferencesToInstance(GKInstance instance) throws Exception {
        Set<GKInstance> current = new HashSet<>();
        Set<GKInstance> next = new HashSet<>();
        current.add(instance);
        while (current.size() > 0) {
            for (GKInstance tmp : current) {
                if (checkedIDs.contains(tmp.getDBID()))
                    continue; // It has been checked
                checkedIDs.add(tmp.getDBID());
                // Check if an event should be in a slice
                if (tmp.getSchemClass().isa("Event") &&
                   !eventMap.containsKey(tmp.getDBID())) 
                   continue;
                pushToMap(tmp, sliceMap);
                // Don't exclude the following classes. Otherwise DB_IDs will be messed up.
//                if (tmp.getSchemClass().isa(ReactomeJavaConstants.ReferenceEntity) ||
//                    tmp.getSchemClass().isa(ReactomeJavaConstants.InstanceEdit) ||
//                    tmp.getSchemClass().isa(ReactomeJavaConstants.Person))
//                    continue; // Don't need detailed information for these instances
                getSourceDBA().fastLoadInstanceAttributeValues(tmp);
                for (SchemaAttribute att : (Collection<SchemaAttribute>) tmp.getSchemaAttributes()) {
                    if (!att.isInstanceTypeAttribute() ||
                        att.getName().equals(ReactomeJavaConstants.inferredFrom) ||
                        att.getName().equals(ReactomeJavaConstants.orthologousEvent))
                        continue;
                    List<GKInstance> values = tmp.getAttributeValuesList(att);
                    if (values == null || values.size() == 0)
                        continue;
                    for (GKInstance reference : values) {
                        if (checkedIDs.contains(reference.getDBID()))
                            continue;
                        next.add(reference);
                    }
                }
            }
            current.clear();
            current.addAll(next);
            next.clear();
        }
    }
    
    private void updateHumanTopLevelPathways(Map<Long, GKInstance> human2mouse) throws Exception {
        logger.info("Starting update human top level pathways...");
        MySQLAdaptor targetDBA = getTargetDBA();
        // Try to use transaction
        boolean isTnSupported = targetDBA.supportsTransactions();
        if (isTnSupported)
            targetDBA.startTransaction();
        try {
            for (Long dbId : human2mouse.keySet()) {
                GKInstance humanEvent = targetDBA.fetchInstance(dbId);
                GKInstance mouseEvent = human2mouse.get(dbId);
                GKInstance targetMouseEvent = targetDBA.fetchInstance(mouseEvent.getDBID());
                humanEvent.getAttributeValuesList(ReactomeJavaConstants.orthologousEvent); // Load it first
                humanEvent.addAttributeValue(ReactomeJavaConstants.orthologousEvent,
                                             targetMouseEvent);
                targetDBA.updateInstanceAttribute(humanEvent, ReactomeJavaConstants.orthologousEvent);
            }
            if (isTnSupported)
                targetDBA.commit();
        }
        catch (Exception e) {
            if (isTnSupported)
                targetDBA.rollback();
            logger.error("SlicingEngine.dumpInstances(): " + e, e);
        }
        logger.info("Done updateHumanTopLevelPathways.");
    }
    
    @Override
    public Map<Long, GKInstance> extractEvents() throws Exception {
        // This is a simple version to get all mouse events
        MySQLAdaptor sourceDBA = getSourceDBA();
        GKInstance mouse = sourceDBA.fetchInstance(MOUSE_DB_ID);
        Collection<GKInstance> mouseEvents = sourceDBA.fetchInstanceByAttribute(ReactomeJavaConstants.Event,
                                                                                ReactomeJavaConstants.species,
                                                                                "=",
                                                                                mouse);
        Map<Long, GKInstance> eventMap = mouseEvents.stream().collect(Collectors.toMap(e -> e.getDBID(),
                                                                                       Function.identity()));
        logger.info("Total mouse events: " + eventMap.size());
        return eventMap;      
    }
    
    @Override
    protected boolean shouldEventInSlice(GKInstance event) throws Exception {
        return true; // Always true since we don't specify _doRelease for predicted mouse events.
    }
    
    /**
     * Some of instances should be in the target database already. Remove these instances.
     * @throws Exception
     */
    private void removeTargetInstances() throws Exception {
        logger.info("Starting remove target instances...");
        MySQLAdaptor targetDBA = getTargetDBA();
        List<Long> existed = new ArrayList<>();
        for (Long dbId : sliceMap.keySet()) {
            GKInstance instance = (GKInstance) sliceMap.get(dbId);
            GKInstance targetInst = targetDBA.fetchInstance(instance.getDBID());
            if (targetInst == null)
                continue; // We can use this DB_ID
            if (!instance.getSchemClass().getName().equals(targetInst.getSchemClass().getName())) {
                continue; 
            }
            // This is a weak check to make sure they have the same _displayName
            if (instance.getDisplayName() != null && 
                targetInst.getDisplayName() != null &&
                !instance.getDisplayName().equals(targetInst.getDisplayName())) {
                continue;
            }
            // Make sure they have the same species
            if (instance.getSchemClass().isValidAttribute(ReactomeJavaConstants.species)) {
                GKInstance species = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.species);
                GKInstance targetSpecies = (GKInstance) targetInst.getAttributeValue(ReactomeJavaConstants.species);
                if (species != null && targetSpecies != null && !species.getDBID().equals(targetSpecies.getDBID())) {
                    continue;
                }
            }
            existed.add(dbId);
        }
        logger.info("Total existed instances for mouse pathways: " + existed.size());
        sliceMap.keySet().removeAll(existed);
    }
    
    // To dump the instanaces into the database, the original source database
    // should be modified to add new classes and attributes added for the reactome_plus
    // _i database by source ~/resource/SchemaModification.sql.
    private boolean _dumpInstances() throws Exception {
        logger.info("dumpInstances()...");
        logger.info("Total instances to be dumped: " + sliceMap.size());
        long time1 = System.currentTimeMillis();
        boolean success = false;
        // Try to use transaction
        boolean isTnSupported = targetDBA.supportsTransactions();
        if (isTnSupported)
            targetDBA.startTransaction();
        try {
            // Treat all instances as local ones for handling DB_IDs.
            for (Long dbId : sliceMap.keySet()) {
                GKInstance instance = (GKInstance) sliceMap.get(dbId);
                instance.setDBID(-instance.getDBID());
            }
            targetDBA.storeLocalInstances(new ArrayList<>(sliceMap.values()));
            // Do a fix for pathway diagram
            DiagramGKBWriter writer = new DiagramGKBWriter();
            writer.setPersistenceAdaptor(targetDBA);
            for (GKInstance pd : pdToDiagram.keySet()) {
                RenderablePathway diagram = pdToDiagram.get(pd);
                for (Renderable r : (List<Renderable>) diagram.getComponents()) {
                    if (r.getInstance() != null)
                        r.setReactomeId(r.getInstance().getDBID());
                }
                // Don't forget to update reactomeDigaramId. Otherwise, the old, wrong DB_IDs
                // will be saved.
                diagram.setReactomeDiagramId(pd.getDBID());
                String xml = writer.generateXMLString(diagram);
                // Get the newly saved instance
                GKInstance targetPd = targetDBA.fetchInstance(pd.getDBID()); // DBID should be updated.
                targetPd.setAttributeValue(ReactomeJavaConstants.storedATXML, xml);
                targetDBA.updateInstanceAttribute(targetPd,
                                                  ReactomeJavaConstants.storedATXML);
            }
            if (isTnSupported) 
                targetDBA.commit();
            success = true;
        }
        catch (Exception e) {
            if (isTnSupported)
                targetDBA.rollback();
            logger.error("SlicingEngine.dumpInstances(): " + e, e);
        }
        long time2 = System.currentTimeMillis();
        logger.info("Time for dumpInstances(): " + (time2 - time1));
        return success;
    }
    
    /**
     * Get the map from human top level pathways to mouse top level pathways so that we will be able to modify 
     * human top level pathways for navigation.
     */
    private Map<Long, GKInstance> getTopLevelHumanToMouseMap() throws Exception {
        Collection<GKInstance> frontPages = getSourceDBA().fetchInstancesByClass(ReactomeJavaConstants.FrontPage);
        if (frontPages.size() == 0)
            throw new IllegalStateException("Cannot find any FrontPage instance.");
        if (frontPages.size() > 1)
            throw new IllegalStateException("More than one FrontPage instances have been found!");
        Map<Long, GKInstance> human2mouse = new HashMap<>();
        // There should be one instance
        GKInstance frontPage = frontPages.stream().findAny().get();
        List<GKInstance> frontPageItems = frontPage.getAttributeValuesList(ReactomeJavaConstants.frontPageItem);
        String mouse = "Mus musculus";
        for (GKInstance inst : frontPageItems) {
            // Check predicted pathways
            List<GKInstance> orEvents = inst.getAttributeValuesList(ReactomeJavaConstants.orthologousEvent);
            for (GKInstance orEvent : orEvents) {
                // For inferred pathway, only one species needs to be checked.
                GKInstance species = (GKInstance) orEvent.getAttributeValue(ReactomeJavaConstants.species);
                if (species.getDisplayName().equals(mouse)) {
                    human2mouse.put(inst.getDBID(),
                                    orEvent);
                    break;
                }
            }
        }
        logger.info("Total top-level mouse pathways: " + human2mouse.size());
        return human2mouse;
    }
    
}
