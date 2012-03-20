/*
 * Created on Feb 27, 2006
 *
 */
package org.reactome.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.gk.model.GKInstance;
import org.gk.model.PersistenceAdaptor;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to do some analysis about the panther dataset.
 * @author guanming
 */
public class PantherAnalyzer extends ReactomeAnalyzer {
    private final String DATASET_DIR = "datasets/";
    //private final String MAPPING_FILE_NAME = DATASET_DIR + "Panther/Version1.3/SequenceAssociationPathway1.3";
//    private final String MAPPING_FILE_NAME = DATASET_DIR + "Panther/Version2.5/SequenceAssociationPathway2.5.txt";
    
    public PantherAnalyzer() {
        //dataSourceId = 206011L; for version 1
        //dataSourceId = 210689L; // for version 2
        dataSourceId = 191282L; // version 3.0
    }
    
//    protected PersistenceAdaptor getMySQLAdaptor() throws Exception {
//        if (dba == null) {
////            dba = new MySQLAdaptor("localhost",
////                                   "reactome_plus_i_v2",
////                                   "root",
////                                   "macmysql01",
////                                   3306);
////            dba = new MySQLAdaptor("localhost",
////                                   "manuel_gk_central",
////                                   "root",
////                                   "macmysql01",
////                                   3306);
////            XMLFileAdaptor fileAdpator = new XMLFileAdaptor();
////            fileAdpator.setSource(DATASET_DIR + "Panther/Version1.3/Panther.rtpj");
////            dba = fileAdpator;
//        }
//        return dba;
//    }
    
    protected Collection prepareComplexes() throws Exception {
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) {
            MySQLAdaptor dba = (MySQLAdaptor) adaptor;
            // GKInstance for dataSource pantherdb
            GKInstance pantherdb = getDataSource();
            Collection complexes = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Complex,
                                                                ReactomeJavaConstants.dataSource,
                                                                "=",
                                                                pantherdb);
            SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Complex);
            SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.hasComponent);
            dba.loadInstanceAttributeValues(complexes, att);
            return complexes;
        }
        else {
            XMLFileAdaptor fileadaptor = (XMLFileAdaptor) adaptor;
            return fileadaptor.fetchInstancesByClass(ReactomeJavaConstants.Complex);
        }
    }
    
    protected Collection prepareReactions() throws Exception {
        PersistenceAdaptor adaptor = getMySQLAdaptor();
        if (adaptor instanceof MySQLAdaptor) {
            // Load necessary attributes
            MySQLAdaptor dba = (MySQLAdaptor) adaptor;
            GKInstance pantherdb = getDataSource();
            Collection reactions = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Reaction,
                                                                ReactomeJavaConstants.dataSource,
                                                                "=",
                                                                pantherdb);
            Collection cas = dba.fetchInstancesByClass(ReactomeJavaConstants.CatalystActivity);
            Collection regulations = dba.fetchInstancesByClass(ReactomeJavaConstants.Regulation);
            Collection entities = dba.fetchInstanceByAttribute(ReactomeJavaConstants.EntityWithAccessionedSequence,
                                                               ReactomeJavaConstants.dataSource,
                                                               "=",
                                                               pantherdb);
            // Load precedingEvent values
            SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Reaction);
            SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.input);
            dba.loadInstanceAttributeValues(reactions, att);
            att = cls.getAttribute(ReactomeJavaConstants.output);
            dba.loadInstanceAttributeValues(reactions, att);
            att = cls.getAttribute(ReactomeJavaConstants.catalystActivity);
            dba.loadInstanceAttributeValues(reactions, att);
            cls = dba.getSchema().getClassByName(ReactomeJavaConstants.CatalystActivity);
            att = cls.getAttribute(ReactomeJavaConstants.physicalEntity);
            dba.loadInstanceAttributeValues(cas, att);
            cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Regulation);
            att = cls.getAttribute(ReactomeJavaConstants.regulatedEntity);
            dba.loadInstanceAttributeValues(regulations, att);
            att = cls.getAttribute(ReactomeJavaConstants.regulator);
            dba.loadInstanceAttributeValues(regulations, att);
            cls = dba.getSchema().getClassByName(ReactomeJavaConstants.EntityWithAccessionedSequence);
            att = cls.getAttribute(ReactomeJavaConstants.referenceEntity);
            dba.loadInstanceAttributeValues(entities, att);
            return reactions;
        }
        else {
            XMLFileAdaptor fileAdaptor = (XMLFileAdaptor) adaptor;
            return fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Reaction);
        }
    }
    
    protected List<GKInstance> getTopics() throws Exception {
        PersistenceAdaptor dba = getMySQLAdaptor();
        GKInstance dataSource = getDataSource();
        Collection collection = dba.fetchInstanceByAttribute(ReactomeJavaConstants.Pathway, 
                                                             ReactomeJavaConstants.dataSource,
                                                             "=",
                                                             dataSource);
        List<GKInstance> topics = new ArrayList<GKInstance>();
        for (Iterator it = collection.iterator(); it.hasNext();) {
            GKInstance gkInstance = (GKInstance) it.next();
            topics.add(gkInstance);
        }
        return topics;
    }
    
    public void countPositionToCompartment() throws Exception {
        Set<String> values = new HashSet<String>();
        String fileName = DATASET_DIR + File.separator + "Panther" + File.separator + "PositionToCompartment.txt";
        FileReader reader = new FileReader(fileName);
        BufferedReader bufferedRead = new BufferedReader(reader);
        String line = null;
        Pattern pattern = Pattern.compile(">(\\w*)</");
        Matcher matcher = null;
        while ((line = bufferedRead.readLine()) != null) {
            // >inside</
            matcher = pattern.matcher(line);
            if (matcher.find()) {
                String value = matcher.group(1);
                values.add(value);
            }
        }
        System.out.println("positionToCompartment: " + values);
    }
    
    /**
     * Count how many human proteins are in the Panther pathway files.
     */
    public void countHumanProteins() throws Exception {
        // Get the human proteins accession numbers in UniProts
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> uniAccIDsMap = uniAnalyzer.loadUniProtIDsMap();
        // TREMBL 
        //fileName = DATASET_DIR + "UniProt" + File.separator + "uniprot_trembl_human.dat";
        //processUniProtIds(fileName, uniProtIds);
        System.out.println("Total UniProt: " + uniAccIDsMap.size());
        // Panther map file: Pathway Component to UniProt IDs
        String fileName = DATASET_DIR + "Panther" + File.separator + "SequenceAssociationPathway1.13";
        FileReader fileReader = new FileReader(fileName);
        BufferedReader bufferedReader = new BufferedReader(fileReader);
        String line = null;
        String[] tokens = null;
        Set<String> pantherIds = new HashSet<String>();
        Set<String> nonredundantIds = new HashSet<String>();
        String uniID = null;
        while ((line = bufferedReader.readLine()) != null) {
            tokens = line.split("\t");
            uniID = tokens[4];
            if (uniAccIDsMap.containsKey(uniID)) {
                pantherIds.add(uniID);
                nonredundantIds.add(uniAccIDsMap.get(uniID));
            }
        }
        System.out.println("UniProt in Panther: " + pantherIds.size());
        System.out.println("UnitProt in Panther (nonredundant): " + nonredundantIds.size());
        // Check how many panther ids have been in the Reactome already
//        Set<String> reactomeIds = loadReactomeIds();
//        int c = 0;
//        for (String id : pantherIds) {
//            if (reactomeIds.contains(id))
//                c++;
//        }
//        System.out.println("Panther ID in Reactome: " + c);
    }

    public void compareDataWithDavids() throws Exception {
        String davidFileName = "results/interaction/PantherFromDavidInteractions090606.txt";
        String pantherFileName = "results/interaction/PantherInteractions090706.txt";
        FileUtility fu = new FileUtility();
        Set<String> david = fu.loadInteractions(davidFileName);
        Set<String> newPanther = fu.loadInteractions(pantherFileName);
        System.out.println("Total in David: " + david.size());
        System.out.println("Total in New: " + newPanther.size());
        david.removeAll(newPanther);
        System.out.println("David only: " + david.size());
//      Need to replace the Uniprot Ids
        UniProtAnalyzer uniAnalyzer = new UniProtAnalyzer();
        Map<String, String> idMap = uniAnalyzer.loadUniProtIDsMap();
        int index = 0;
        Set<String> convertedDavid = new HashSet<String>();
        String id1, id2;
        int compare = 0;
        for (String pair : david) {
            index = pair.indexOf(" ");
            id1 = pair.substring(0, index);
            id2 = pair.substring(index + 1);
            if (id1.equals(idMap.get(id1)) &&
                id2.equals(idMap.get(id2)))
                convertedDavid.add(pair);
//            id1 = idMap.get(id1);
//            id2 = idMap.get(id2);
//            compare = id1.compareTo(id2);
//            if (compare < 0)
//                convertedDavid.add(id1 + " " + id2);
//            else
//                convertedDavid.add(id2 + " " + id1);
        }
        System.out.println("Converted: " + convertedDavid.size());
        // Want to print 100
        List<String> list = new ArrayList<String>(convertedDavid);
        for (int i = 0; i < 100; i++)
            System.out.println(list.get(i));
    }
}
