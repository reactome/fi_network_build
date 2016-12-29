/*
 * Created on Jan 30, 2007
 *
 */
package org.reactome.kegg;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.gk.database.SynchronizationManager;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.PersistenceManager;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.junit.Before;
import org.junit.Test;
import org.reactome.convert.common.PostProcessHelper;
import org.reactome.fi.util.FIConfiguration;


/**
 * This class is used to convert KEGG protein network to Reactome. KGML data is used
 * for such conversion.
 * @author guanming
 *
 */
public class KeggToReactomeConverter {
    private static Logger logger = Logger.getLogger(KeggToReactomeConverter.class);
    
    private XMLFileAdaptor fileAdaptor;
    // cached value
    private GKInstance keggDb;
    // Want to cache these values to avoid duplication 
    private Map<String, GKInstance> keggIdToInstance;
    // Used to map the instances 
    private MySQLAdaptor dbAdaptor;
    
    public KeggToReactomeConverter() throws Exception {
        fileAdaptor = new XMLFileAdaptor();
        keggIdToInstance = new HashMap<String, GKInstance>();
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dbAdaptor) {
        this.dbAdaptor = dbAdaptor;
    }
    
    public void convert(List<String> fileNames) throws Exception {
        for (String fileName : fileNames) {
            logger.info("Convert " + fileName + "...");
            convertBeforePost(fileName);
        }
        postProcess();
    }
    
    public void convert(String fileName) throws Exception {
        convertBeforePost(fileName);
        postProcess();
    }

    private void convertBeforePost(String fileName) throws Exception {
        // Load the file first
        Document doc = loadKeggPathway(fileName);
        Element root = doc.getRootElement();
        GKInstance pathway = createPathway(root);
        // Get all entries"
        Map<String, GKInstance> entryIdToInstance = convertEntries(root);
        // Second pass to convert interactions
        convertRelations(root,
                         entryIdToInstance,
                         pathway);
        // Post process: map to UniProt
    }
    
    /**
     * Post process is used to map EWAS to RefPepSeq in the database
     * @throws Exception
     */
    private void postProcess() throws Exception {
        KEGGPostProcess processor = new KEGGPostProcess();
        processor.postProcess(dbAdaptor, 
                              fileAdaptor);
    }
    
    private GKInstance createPathway(Element pathwayElm) throws Exception {
        String name = pathwayElm.getAttributeValue("name");
        GKInstance pathway = keggIdToInstance.get(name);
        if (pathway == null)
            pathway = createPathwayFromName(name);
        String title = pathwayElm.getAttributeValue("title");
        pathway.addAttributeValue(ReactomeJavaConstants.name,
                                  title);
        pathway.setDisplayName(title);
        return pathway;
    }
    
    private void convertRelations(Element root,
                                  Map<String, GKInstance> idToInstance,
                                  GKInstance pathway) throws Exception {
        List list = root.getChildren("relation");
        Element relationElm = null;
        for (Iterator it = list.iterator(); it.hasNext();) {
            relationElm = (Element) it.next();
//            <relation entry1="13" entry2="10" type="GErel">
//                <subtype name="expression" value="-->"/>
//                <subtype name="indirect" value="..>"/>
//            </relation>
            String type = relationElm.getAttributeValue("type");
            if (type.equals("maplink"))
                continue; //TODO: cannot handle maplink right now
            String id1 = relationElm.getAttributeValue("entry1");
            GKInstance interactor1 = idToInstance.get(id1);
            if (interactor1 == null)
                throw new IllegalStateException("convertRelations(): " + id1 + " is not mapped");
            String id2 = relationElm.getAttributeValue("entry2");
            GKInstance interactor2 = idToInstance.get(id2);
            if (interactor2 == null)
                throw new IllegalStateException("convertRelations(): " + id2 + " is not mapped");
            // As of December 6, 2013, relations between pathways have been encoded in KGML.
            // For example, see http://www.genome.jp/kegg-bin/show_pathway?hsa05161.
            // Currently we don't support such kind of relations
            if (interactor1.getSchemClass().isa(ReactomeJavaConstants.Event) ||
                interactor2.getSchemClass().isa(ReactomeJavaConstants.Event)) {
                logger.warn("Relation is related to event: " + interactor1 + " and " + interactor2);
                continue;
            }
            // extract type
            type = extractType(relationElm);
            // Use lower case: some bugs in the pathways
            String id = interactor1.getDBID() + " " + type.toLowerCase() + " " + interactor2.getDBID();
            // Same interaction might appear several times. E.g. hsa04730.xml.
            GKInstance interaction = keggIdToInstance.get(id);
            if (interaction != null) {
                // Check if new schema is used
                if (pathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent))
                    pathway.addAttributeValue(ReactomeJavaConstants.hasComponent,
                                              interaction);
                else if (pathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                    pathway.addAttributeValue(ReactomeJavaConstants.hasEvent,
                                              interaction);
                continue;
            }
            interaction = fileAdaptor.createNewInstance(ReactomeJavaConstants.Interaction);
            keggIdToInstance.put(id, interaction); 
            interaction.addAttributeValue(ReactomeJavaConstants.interactor,
                                          interactor1);
            interaction.addAttributeValue(ReactomeJavaConstants.interactor,
                                          interactor2);
            
            interaction.addAttributeValue(ReactomeJavaConstants.interactionType, 
                                          type);
            // Check if new schema is used
            if (pathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasComponent))
                pathway.addAttributeValue(ReactomeJavaConstants.hasComponent,
                                          interaction);
            else if (pathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent))
                pathway.addAttributeValue(ReactomeJavaConstants.hasEvent,
                                          interaction);
        }
    }
    
    private String extractType(Element relationElm) {
        StringBuilder builder = new StringBuilder();
        String type = relationElm.getAttributeValue("type");
        builder.append(type);
        // Check subtypes
        List subTypes = relationElm.getChildren("subtype");
        if (subTypes.size() > 0) {
            builder.append(": ");
            for (Iterator it = subTypes.iterator(); it.hasNext();) {
                Element elm = (Element) it.next();
                builder.append(elm.getAttributeValue("name"));
                if (it.hasNext())
                    builder.append(", ");
            }
        }
        return builder.toString();
    }
    
    private GKInstance getKEGGDB() throws Exception {
        if (keggDb != null)
            return keggDb;
        keggDb = fileAdaptor.createNewInstance(ReactomeJavaConstants.ReferenceDatabase);
        keggDb.addAttributeValue(ReactomeJavaConstants.name, "KEGG");
        keggDb.addAttributeValue(ReactomeJavaConstants.url, "http://www.genome.jp/kegg/");
        keggDb.setDisplayName("KEGG");
        return keggDb;
    }
    
    private Map<String, GKInstance> convertEntries(Element root) throws Exception {
        List entries = root.getChildren("entry");
        Element entryElm = null;
        GKInstance instance = null;
        Map<String, GKInstance> entryIdToInstance = new HashMap<String, GKInstance>();
        // First pass to convert to instances
        // Note: a GKInstance might be repeated several times in the KEGG pathways.
        // e.g.GNA13... in hsa04730.xml. Use GKInstance to Element map.
        List<Element> groupElms = new ArrayList<Element>();
        for (Iterator it = entries.iterator(); it.hasNext();) {
            entryElm = (Element) it.next();
            // Defer the group 
            String type = entryElm.getAttributeValue("type");
            if (type.equals("group")) {
                groupElms.add(entryElm);
                continue;
            }
            instance = convertNonGroupEntry(entryElm);
            if (instance == null)
                continue;
            String id = entryElm.getAttributeValue("id");
            entryIdToInstance.put(id,
                                  instance);
        }
        // Handle group elements
        for (Element groupElm : groupElms) {
            instance = convertGroupEntry(groupElm, entryIdToInstance);
            String id = groupElm.getAttributeValue("id");
            entryIdToInstance.put(id,
                                  instance);
        }
        return entryIdToInstance;
    }
    
    private GKInstance convertGroupEntry(Element groupElm,
                                         Map<String, GKInstance> entryIdToInstance) throws Exception {
        List componentList = groupElm.getChildren("component");
        List<String> compIds = new ArrayList<String>();
        for (Iterator it = componentList.iterator(); it.hasNext();) {
            Element comp = (Element) it.next();
            String id = comp.getAttributeValue("id");
            compIds.add(id);
        }
        // Generate a key
        StringBuilder builder = new StringBuilder();
        for (String id : compIds) {
            GKInstance comp = entryIdToInstance.get(id);
            // ortholog might be used, which cannot be mapped right now
            if (comp == null) {
                logger.warn("Complex component " + id + " is not mapped!");
                continue;
            }
            builder.append(comp.getDBID());
            builder.append("+");
        }
        String key = builder.toString();
        GKInstance complex = keggIdToInstance.get(key);
        if (complex != null)
            return complex;
        complex = fileAdaptor.createNewInstance(ReactomeJavaConstants.Complex);
        keggIdToInstance.put(key, complex);
        // Add components
        for (String id : compIds) {
            GKInstance comp = entryIdToInstance.get(id);
            complex.addAttributeValue(ReactomeJavaConstants.hasComponent,
                                      comp);
        }
        return complex;
    }
    
    private void handleDefinedSet(GKInstance set,
                                  Element entryElm) throws Exception {
        String names = entryElm.getAttributeValue("name");
        String[] tokens = names.split(" ");
        for (String name : tokens) {
            GKInstance ewas = createEwasFromName(name);
            set.addAttributeValue(ReactomeJavaConstants.hasMember, 
                                  ewas);
        }
    }
    
    private GKInstance createEwasFromName(String name) throws Exception {
        GKInstance ewas = keggIdToInstance.get(name);
        if (ewas != null)
            return ewas;
        ewas = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        keggIdToInstance.put(name, ewas);
        // Need to create a crossReference for kegg id
        GKInstance reference = createRefIdFromName(name);
        ewas.addAttributeValue(ReactomeJavaConstants.crossReference, 
                               reference);
        return ewas;
    }
    
    private GKInstance createPathwayFromName(String name) throws Exception {
        GKInstance pathway = keggIdToInstance.get(name);
        if (pathway != null)
            return pathway;
        pathway = fileAdaptor.createNewInstance(ReactomeJavaConstants.Pathway);
        keggIdToInstance.put(name, pathway);
        // Need to create a crossReference for kegg id
        GKInstance reference = createRefIdFromName(name);
        pathway.addAttributeValue(ReactomeJavaConstants.crossReference, 
                                  reference);
        return pathway;
    }
    
    private GKInstance createRefIdFromName(String name) throws Exception {
        GKInstance reference = fileAdaptor.createNewInstance(ReactomeJavaConstants.DatabaseIdentifier);
        reference.setAttributeValue(ReactomeJavaConstants.identifier, name);
        reference.setAttributeValue(ReactomeJavaConstants.referenceDatabase, 
                                    getKEGGDB());
        InstanceDisplayNameGenerator.setDisplayName(reference);
        return reference;
    }
    
    private GKInstance convertNonGroupEntry(Element entryElm) throws Exception {
        String type = entryElm.getAttributeValue("type");
        String name = entryElm.getAttributeValue("name");
        // Try to get the converted one
        GKInstance instance = null;
        if (name != null && name.length() > 0) {
            instance = keggIdToInstance.get(name);
            if (instance != null)
                return instance;
        }
//      ortholog
//      enzyme
//      gene    the node is a gene product (mostly a protein)
//      group   the node is a complex of gene products (mostly a protein complex)
//      compound
//      map
//      brite   see this type in the 2016 version of KGML. This is a link to ATC, which is for drug identification.
//              for the time being. It will be converted into a simple small molecules. 
//      other   this type should be removed in KGML 0.70 in 2010. However,  4 entries are found in pathway hsa03320 (
//              PPAR Signaling Pathway). They are groups of drugs. Just map them into other type in Reactome.        
        if (type.equals("enzyme") ||
            type.equals("gene")) {
            String[] tokens = name.split(" ");
            if (tokens.length == 1) {
                instance = createEwasFromName(name);
            }
            else if (tokens.length > 1) {
                // DefinedSet is needed
                instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.DefinedSet);
                handleDefinedSet(instance, entryElm);
            }
        }
        //else if (type.equals("group")) {
        //    instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.Complex);
        //}
        else if (type.equals("compound") || type.equals("brite"))
            instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.SimpleEntity);
        else if (type.equals("map")) {
            instance = createPathwayFromName(name);
        }
        else if (type.equals("other")) // Used an other type in 2015 version of KEGG: hsa03320.
            instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.OtherEntity); 
        //TODO: need to map "Ortholog"
        // Try to get the name from graphics
        if (instance != null) {
            // Check if name is assigned
            if (instance.getAttributeValue(ReactomeJavaConstants.name) == null) {
                Element graphicElm = entryElm.getChild("graphics");
                String displayName = graphicElm.getAttributeValue("name");
                instance.addAttributeValue(ReactomeJavaConstants.name, displayName);
                instance.setDisplayName(displayName);
            }
            if (name != null && name.length() > 0)
                keggIdToInstance.put(name, instance);
        }
        return instance;
    }
    
    private Document loadKeggPathway(String fileName) throws Exception {
        SAXBuilder builder = new SAXBuilder();
        Document doc = builder.build(fileName);
        return doc;
    }
    
    public XMLFileAdaptor getReactomeModel() {
        return fileAdaptor;
    }
    
    @Before
    public void setUpTest() throws Exception {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    @Test
    public void testBatchConvert() throws Exception {
        KeggAnalyzer analyzer = new KeggAnalyzer();
        List<String> pathwayNames = analyzer.getNonMetabolismPathways();
        //String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/KEGG/hsa/";
        String dirName = FIConfiguration.getConfiguration().get("KEGG_HSA_KGML_DIR");
        List<String> fileNames = new ArrayList<String>();
        for (String pName : pathwayNames) {
            System.out.println("file name: " + pName);
            fileNames.add(dirName + pName + ".xml");
        }
        dbAdaptor = new MySQLAdaptor("localhost",
                                     "reactome_39_plus_i",
                                     "root",
                                     "macmysql01",
                                     3306);
//        dbAdaptor = new MySQLAdaptor("localhost",
//                                     "gk_central_031309",
//                                     "root",
//                                     "macmysql01",
//                                     3306);
//        dbAdaptor = new MySQLAdaptor("localhost",
//                                     "reactome_28_plus_i_myisam",
//                                     "root",
//                                     "macmysql01",
//                                     3306);
        convert(fileNames);
        //fileAdaptor.save(FIConfiguration.getConfiguration().get("KEGG_DIR + "kegg.rtpj");
//        fileAdaptor.save(FIConfiguration.getConfiguration().get("KEGG_DIR + "kegg101110.rtpj");
        fileAdaptor.save(FIConfiguration.getConfiguration().get("KEGG_CONVERTED_FILE"));
    }
    
    @Test
    public void runBatchConvert() throws Exception {
        // These three variables should be modified for each construction
        String dirName = FIConfiguration.getConfiguration().get("KEGG_HSA_KGML_DIR");
        String outFileName = FIConfiguration.getConfiguration().get("KEGG_CONVERTED_FILE");
        dbAdaptor = new MySQLAdaptor("localhost",
                                     FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                     FIConfiguration.getConfiguration().get("DB_USER"),
                                     FIConfiguration.getConfiguration().get("DB_PWD"));
        
        List<String> fileNames = new ArrayList<String>();
        File dir = new File(dirName);
        for (File file : dir.listFiles()) {
            String fileName = file.getName();
            if(fileName.endsWith(".xml"))
                fileNames.add(file.getAbsolutePath());
        }
        
        convert(fileNames);
        fileAdaptor.save(outFileName);
    }
    
    @Test
    public void testConvert() throws Exception {
        String dir = "/Users/wgm/Documents/caBIG_R3/datasets/KEGG/031209/hsa/";
        String fileName = dir + "hsa05010.xml";
        KeggToReactomeConverter converter = new KeggToReactomeConverter();
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
                                                  "reactome_28_plus_i",
                                                  "root",
                                                  "macmysql01",
                                                  3306);
        converter.setMySQLAdaptor(dbAdaptor);
        converter.convert(fileName);
        XMLFileAdaptor fileAdaptor = converter.getReactomeModel();
        String outFileName = dir + "hsa05010.rtpj";
        fileAdaptor.save(outFileName);
    }
    
    @Test
    public void fixHasEventProblemInKEGG() throws Exception {
        String srcFileName = FIConfiguration.getConfiguration().get("KEGG_DIR") + "kegg101110_fixed.rtpj";
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        fileAdaptor.setSource(srcFileName);
//        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
//                                                   "reactome_28_plus_i",
//                                                   "root",
//                                                   "macmysql01",
//                                                   3306);
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("brie8.cshl.edu",
                                                  "test_reactome_28_plus_i",
                                                  "authortool",
                                                  "T001test",
                                                  3306);
        PersistenceManager.getManager().setActiveFileAdaptor(fileAdaptor);
        PersistenceManager.getManager().setActiveMySQLAdaptor(dbAdaptor);
        SynchronizationManager manager = SynchronizationManager.getManager();
        // Update small molecules
//        GKInstance keggDbInst = dbAdaptor.fetchInstance(486224L);
//        Collection<?> dbInsts = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.SimpleEntity,
//                                                                   ReactomeJavaConstants.dataSource,
//                                                                   "=",
//                                                                   keggDbInst);
//        Collection<?> localInsts = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.SimpleEntity);
//        for (Iterator<?> it = localInsts.iterator(); it.hasNext();) {
//            GKInstance inst = (GKInstance) it.next();
//            String name = (String) inst.getAttributeValue(ReactomeJavaConstants.name);
//            for (Iterator<?> it1 = dbInsts.iterator(); it1.hasNext();) {
//                GKInstance dbInst = (GKInstance) it1.next();
//                String dbName = (String) dbInst.getAttributeValue(ReactomeJavaConstants.name); 
//                if (dbName.equals(name)) {
//                    PostProcessHelper.updateFromDB(inst, dbInst, manager);
//                    it1.remove();
//                    break;
//                }
//            }
//        }
//        updateInstancesOnDisplayName(ReactomeJavaConstants.DatabaseIdentifier, 
//                                     fileAdaptor,
//                                     dbAdaptor,
//                                     manager);
//        updateInstancesOnCrossReference(ReactomeJavaConstants.EntityWithAccessionedSequence, 
//                                        fileAdaptor,
//                                        dbAdaptor, 
//                                        manager);
//        updateInstancesInClass(ReactomeJavaConstants.DefinedSet,
//                               fileAdaptor,
//                               dbAdaptor, 
//                               manager);
//        updateInstancesOnDisplayName(ReactomeJavaConstants.Complex,
//                               fileAdaptor,
//                               dbAdaptor, 
//                               manager);
//        // A special case for complex
//        localInsts = fileAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
//                                                          ReactomeJavaConstants._displayName,
//                                                          "=",
//                                                          "unknown");
//        GKInstance inst = (GKInstance) localInsts.iterator().next();
//        localInsts = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
//                                                          ReactomeJavaConstants._displayName,
//                                                          "=",
//                                                          "unknown");
//        GKInstance dbInst = (GKInstance) localInsts.iterator().next();
//        PostProcessHelper.updateFromDB(inst, dbInst, manager);
//        updateInstancesInClass(ReactomeJavaConstants.Interaction,
//                               fileAdaptor,
//                               dbAdaptor, 
//                               manager);
//        // Three special cases
//        updateInstanceOnDBID(-10112L, fileAdaptor, 478071L, dbAdaptor, manager);
//        updateInstanceOnDBID(-10111L, fileAdaptor, 478070L, dbAdaptor, manager);
//        updateInstanceOnDBID(-10102L, fileAdaptor, 478061L, dbAdaptor, manager);
        // Check Pathways finally
        Collection<?> localInsts = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Pathway);
        for (Iterator<?> it = localInsts.iterator(); it.hasNext();) {
            GKInstance inst = (GKInstance) it.next();
            GKInstance crossRef = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.crossReference);
            Collection c = dbAdaptor.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway,
                                                              ReactomeJavaConstants.crossReference,
                                                              "=",
                                                              crossRef);
            if (c == null || c.size() != 1) {
                System.out.println(inst + " cannot be mapped correctly!");
                continue;
            }
            GKInstance dbInst = (GKInstance) c.iterator().next();
            // Check hasEvents
            List<?> localValues = inst.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
            List<?> dbValues = dbInst.getAttributeValuesList(ReactomeJavaConstants.hasEvent);
            if (localValues.size() == dbValues.size()) {
                // Update
                PostProcessHelper.updateFromDB(inst, dbInst, manager);
            }
            else
                inst.setDBID(dbInst.getDBID());
        }
        fileAdaptor.save(FIConfiguration.getConfiguration().get("KEGG_DIR") + "kegg101110_fixed_1.rtpj");
    }
    
    private void updateInstancesOnDisplayName(String clsName,
                                              XMLFileAdaptor fileAdaptor,
                                              MySQLAdaptor dbAdaptor,
                                              SynchronizationManager manager) throws Exception {
        Collection<?> localInsts = fileAdaptor.fetchInstancesByClass(clsName);
        for (Iterator<?> it = localInsts.iterator(); it.hasNext();) {
            GKInstance inst = (GKInstance) it.next();
            if (inst.getDBID() > 0)
                continue;
            Collection<?> c = dbAdaptor.fetchInstanceByAttribute(clsName,
                                                                 ReactomeJavaConstants._displayName,
                                                                 "=",
                                                                 inst.getDisplayName());
            if (c == null || c.size() == 0) {
                System.out.println("Cannot be matched: " + inst);
                continue;
            }
            if (c.size() > 1)
                System.out.println("More than one mappings: " + inst);
            for (Iterator<?> it1 = c.iterator(); it1.hasNext();) {
                GKInstance dbInst = (GKInstance) it1.next();
                GKInstance dataSource = (GKInstance) dbInst.getAttributeValue(ReactomeJavaConstants.dataSource);
                if (dataSource != null && dataSource.getDisplayName().equals("KEGG")) {
                    PostProcessHelper.updateFromDB(inst, dbInst, manager);
                }
            }
        }
    }
    
    private void updateInstancesOnCrossReference(String clsName,
                                                 XMLFileAdaptor fileAdaptor,
                                                 MySQLAdaptor dbAdaptor,
                                                 SynchronizationManager manager) throws Exception {
        Collection<?> localInsts = fileAdaptor.fetchInstancesByClass(clsName);
        for (Iterator<?> it = localInsts.iterator(); it.hasNext();) {
            GKInstance inst = (GKInstance) it.next();
            if (inst.getDBID() > 0)
                continue;
            GKInstance crossRef = (GKInstance) inst.getAttributeValue(ReactomeJavaConstants.crossReference);
            Collection<?> c = dbAdaptor.fetchInstanceByAttribute(clsName,
                                                                 ReactomeJavaConstants.crossReference,
                                                                 "=",
                                                                 crossRef);
            if (c == null || c.size() == 0) {
                System.out.println("Cannot be matched: " + inst);
                continue;
            }
            for (Iterator<?> it1 = c.iterator(); it1.hasNext();) {
                GKInstance dbInst = (GKInstance) it1.next();
                GKInstance dataSource = (GKInstance) dbInst.getAttributeValue(ReactomeJavaConstants.dataSource);
                if (dataSource != null && dataSource.getDisplayName().equals("KEGG")) {
                    PostProcessHelper.updateFromDB(inst, dbInst, manager);
                }
            }
        }
    }
    
    private void updateInstanceOnDBID(Long localDbId,
                                      XMLFileAdaptor fileAdaptor,
                                      Long dbId,
                                      MySQLAdaptor dbAdaptor,
                                      SynchronizationManager manager) throws Exception {
        GKInstance inst = fileAdaptor.fetchInstance(localDbId);
        GKInstance dbInst = dbAdaptor.fetchInstance(dbId);
        PostProcessHelper.updateFromDB(inst, dbInst, manager);
    }

    private void updateInstancesInClass(String clsName,
                                        XMLFileAdaptor fileAdaptor,
                                        MySQLAdaptor dbAdaptor,
                                        SynchronizationManager manager) throws Exception, InvalidAttributeException {
        Collection<?> localInsts;
        localInsts = fileAdaptor.fetchInstancesByClass(clsName);
        for (Iterator<?> it = localInsts.iterator(); it.hasNext();) {
            GKInstance inst = (GKInstance) it.next();
            Set<?> matched = dbAdaptor.fetchIdenticalInstances(inst);
            if (matched == null || matched.size() == 0) {
                System.out.println("Cannot be matched: " + inst);
            }
            else if (matched.size() > 0) {
                for (Iterator<?> it1 = matched.iterator(); it1.hasNext();) {
                    GKInstance dbInst = (GKInstance) it1.next();
                    GKInstance dataSource = (GKInstance) dbInst.getAttributeValue(ReactomeJavaConstants.dataSource);
                    if (dataSource != null && dataSource.getDisplayName().equals("KEGG")) {
                        PostProcessHelper.updateFromDB(inst, dbInst, manager);
                    }
                }
            }
            else {
                GKInstance dbInst = (GKInstance) matched.iterator().next();
                PostProcessHelper.updateFromDB(inst, dbInst, manager);
            }
        }
    }
}