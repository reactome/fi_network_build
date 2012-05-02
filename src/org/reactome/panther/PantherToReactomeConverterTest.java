/*
 * Created on Feb 28, 2006
 *
 */
package org.reactome.panther;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.log4j.PropertyConfigurator;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.junit.Test;
import org.reactome.convert.common.Converter;
import org.reactome.convert.common.ConverterHandler;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

public class PantherToReactomeConverterTest {
    //private final String PANTHER_DIR = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/Version1.3/SBML_1.3/SBML/";
    //private final String PANTHER_DIR = "/Users/wgm/Documents/gkteam/Arabidopsis/";
    private PantherToReactomeConverter converter;
    
    public PantherToReactomeConverterTest() {
        try {
            setUp();
        }
        catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    protected void setUp() throws Exception {
        PropertyConfigurator.configure("resources/log4j.properties");
//        AuthorToolAppletUtilities.setSchemaFileName("schema_extended");
        XMLFileAdaptor reactomeAdaptor = new XMLFileAdaptor();
        converter = new PantherToReactomeConverter();
        converter.setFileAdaptor(reactomeAdaptor);
//        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
//                                                  "gk_central_101606",
//                                                  "root",
//                                                  "macmysql01",
//                                                  3306);
        MySQLAdaptor dbAdaptor = getDBA();
        converter.setDatabaseAdaptor(dbAdaptor);
//        super.setUp();
    }
    
    private MySQLAdaptor getDBA() throws Exception {
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_release_28_human_fis",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "gk_central_031309",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                            FIConfiguration.getConfiguration().get("DB_USER"),
                                            FIConfiguration.getConfiguration().get("DB_PWD"),
                                            3306);
        return dba;
    }
    
    /**
     * This method is used to test Andreash's new implementations based on model
     * versioning.
     * @throws Exception
     */
    @Test
    public void testNewConverter() throws Exception {
        MySQLAdaptor dbAdaptor = getDBA();
        ConverterHandler handler = ConverterHandler.getInstance();
        String pantherModelName = "Angiotensin_II-stimulated_signaling_through_G_proteins_and_beta-arrestin.xml";
        pantherModelName = FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR") + pantherModelName;
        Converter cdConverter = handler.getConverter(pantherModelName);
        XMLFileAdaptor reactomeAdaptor = new XMLFileAdaptor();
        cdConverter.setFileAdaptor(reactomeAdaptor);
        cdConverter.setDatabaseAdaptor(dbAdaptor);
        GKInstance pathway = cdConverter.convert(pantherModelName);
        GKInstance summation = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.summation);
        if (summation != null) {
            String text = (String) summation.getAttributeValue(ReactomeJavaConstants.text);
            System.out.println(pathway.getDisplayName() + ": " + text);
        }
        String destFileName = FIConfiguration.getConfiguration().get("PANTHER_DIR") + "Angiotensin_II.rtpj";
        cdConverter.save(destFileName);
    }
    
    /**
     * This test method is using Andreash's new implementation based on model versioning.
     * Note only one XMLFileAdaptor is used so that all converted instances can be in the
     * same Reactome Curator Tool project file.
     * @throws Exception
     */
    @Test
    public void testNewBatchConverter() throws Exception {
        MySQLAdaptor dbAdaptor = getDBA();
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        SchemaClass refGeneProd = fileAdaptor.getSchema().getClassByName(ReactomeJavaConstants.ReferenceGeneProduct);
        System.out.println("ReferenceGeneProduct: " + refGeneProd.getName());
        String dirName = FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR");
        System.out.println(dirName);
        File[] sbmlFiles = new File(dirName).listFiles();
        System.out.println("Total files: " + sbmlFiles.length);
        List<String> fileNames = new ArrayList<String>();
        // Load all model files
        for (File file : sbmlFiles) {
            String fileName = file.getName();
            if (fileName.endsWith(".xml"))
                fileNames.add(file.getAbsolutePath());
        }
        // for test
//        fileNames.clear();
//        fileNames.add(FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR") + "Inflammation_mediated_by_chemokine_and_cytokine_signaling_pathway.xml");
        // Add these statements so that this pathway can be handled first.
        // Otherwise, some instances in this pathway may be messed up with others.
        String tmp = FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR") + "Inflammation_mediated_by_chemokine_and_cytokine_signaling_pathway.xml";
        fileNames.remove(tmp);
        fileNames.add(0, tmp);
        ConverterHandler handler = ConverterHandler.getInstance();
        for (String fileName : fileNames) {
            System.out.println("Converting " + fileName + "...");
            Converter converter = handler.getConverter(fileName);
            converter.setDatabaseAdaptor(dbAdaptor);
            converter.setFileAdaptor(fileAdaptor);
            GKInstance pathway = converter.convertBeforePost(fileName);
            System.out.println("Done: " + pathway.getDisplayName());
        }
        // Do post-processing.
        PantherPostProcessor postProcessor = new PantherPostProcessor();
        postProcessor.postProcess(dbAdaptor, fileAdaptor);
        postProcessor.otherProcesses(dbAdaptor, fileAdaptor);
        fileAdaptor.save(FIConfiguration.getConfiguration().get("PANTHER_CONVERTED_FILE"));
    }
    
    public void attachSpecies() throws Exception {
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("brie8.cshl.edu",
                                                  "test_converted_pathways",
                                                  "authortool",
                                                  "T001test",
                                                  3306);
        // This is pantherdb data sources
        GKInstance pantherDS = dbAdaptor.fetchInstance(new Long(210683L));
        assert pantherDS != null;
        GKInstance humanSpecies = dbAdaptor.fetchInstance(new Long(48887L));
        assert humanSpecies != null;
        String[] clses = new String[] {
                ReactomeJavaConstants.Event,
                ReactomeJavaConstants.Complex,
                ReactomeJavaConstants.EntitySet,
                ReactomeJavaConstants.GenomeEncodedEntity,
                ReactomeJavaConstants.Polymer
        };
        try {
            dbAdaptor.startTransaction();
            for (String cls : clses) {
                System.out.println("handling " + cls + "...");
                Collection c = dbAdaptor.fetchInstanceByAttribute(cls, 
                                                                  ReactomeJavaConstants.dataSource, 
                                                                  "=", 
                                                                  pantherDS);
                if (c == null || c.size() == 0)
                    continue;
                SchemaClass schemaCls = dbAdaptor.getSchema().getClassByName(cls);
                SchemaAttribute att = schemaCls.getAttribute(ReactomeJavaConstants.species);
                dbAdaptor.loadInstanceAttributeValues(c, att);
                for (Iterator it = c.iterator(); it.hasNext();) {
                    GKInstance instance = (GKInstance) it.next();
                    GKInstance species = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.species);
                    if (species != null)
                        continue;
                    instance.setAttributeValue(ReactomeJavaConstants.species, humanSpecies);
                    dbAdaptor.updateInstanceAttribute(instance, 
                                                      ReactomeJavaConstants.species);
                }
            }
            dbAdaptor.commit();
        }
        catch (Exception e){
            dbAdaptor.rollback();
            throw e;
        }
    }
    
    /**
     * Use ReactomeProjectDumper for dumping the database.
     * @throws Exception
     */
//    public void testDumpToDB() throws Exception {
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "converted_pathways",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
//        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/Version1.3/Panther_1_3.rtpj";
//        PantherPostProcessor postProcessor = new PantherPostProcessor();
//        postProcessor.dumpToDB(fileName, dba);
//    }
    
    public void testBatchConvert() throws Exception {
        //String destDir = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/Version1.3/";
        //String destFileName = destDir + "Panther_1_3.rtpj";
        String destFileName = FIConfiguration.getConfiguration().get("PANTHER_CONVERTED_FILE");
        File[] sbmlFiles = new File(FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR")).listFiles();
        List<String> fileNames = new ArrayList<String>();
        for (File file : sbmlFiles) {
            String fileName = file.getName();
            if (fileName.endsWith(".xml"))
                fileNames.add(file.getAbsolutePath());
        }
        converter.convert(fileNames);
        converter.getFileAdaptor().save(destFileName);
    }
    
    public void batchConvertForIndividualFiles() throws Exception {
        String destDir = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/Version1.3/converted/";
        File[] sbmlFiles = new File(FIConfiguration.getConfiguration().get("PANTHER_DIR")).listFiles();
        for (File file : sbmlFiles) {
            String fileName = file.getName();
            if (fileName.endsWith(".xml")) {
                XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
                converter.setFileAdaptor(fileAdaptor);
                GKInstance pathway = converter.convert(file.getAbsolutePath());
                String destFile = destDir + fileName.substring(0, fileName.length() - 4) + ".rtpj";
                converter.getFileAdaptor().save(destFile);
            }
        }
    }
    
    public void testConvert() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        MySQLAdaptor dbAdaptor = getDBA();
        // Have to generate a JDOM Document object
        //String pantherModelName = "Androgen_estrogen_progesterone_biosynthesis.xml";
        //String pantherModelName = "Ascorbate_degradation.xml";
        //String pantherModelName = "Nicotinic_acetylcholone_receptor_signaling_pathway.xml";
        //String pantherModelName = "Adrenaline_synthesis.xml";
        //String pantherModelName = "Wnt_signaling_pathway.xml";
        //String pantherModelName = "Ionotropic_glutamate_receptor_pathway.xml";
        //String pantherModelName = "Angiogenesis.xml";
        //String pantherModelName = "Plasminogen_activating_cascade.xml";
        //String pantherModelName = "TGF_beta_signaling_pathway.xml";
        //String pantherModelName = "Heterotrimeric_G_protein_signaling_pathway_Gi_and_Gs_mediated_pathway.xml";
        //String pantherModelName = "Cholesterol_biosynthesis.xml";
        //String pantherModelName = "Alzheimer_disease_presenilin_pathway.xml";
        //String pantherModelName = "ath00010_L2.xml"; // From Arabidopsis from KEGG
        //String pantherModelName = "Angiotensin_II-stimulated_signaling_through_G_proteins_and_beta-arrestin.xml";
        String pantherModelName = "GABA-B_receptor_II_pathway.xml";
        pantherModelName = FIConfiguration.getConfiguration().get("PANTHER_FILES_DIR") + pantherModelName;
        //pantherModelName = "/Users/wgm/Documents/gkteam/Lisa/actin dendritic nucleation - fixed.xml";
        ConverterHandler handler = ConverterHandler.getInstance();
        Converter converter = handler.getConverter(pantherModelName);
        converter.setFileAdaptor(fileAdaptor);
        converter.setDatabaseAdaptor(dbAdaptor);
        
        GKInstance pathway = converter.convertBeforePost(pantherModelName);
        System.out.println("Pathway Name: " + pathway.getDisplayName());
        GKInstance summation = (GKInstance) pathway.getAttributeValue(ReactomeJavaConstants.summation);
        if (summation != null)
            System.out.println("Summation: " + summation.getAttributeValue(ReactomeJavaConstants.text));
        //String destFileName = "/Users/wgm/Documents/gkteam/Lisa/" + pathway.getDisplayName() + ".rtpj";
        //String destFileName = "/Users/wgm/Documents/gkteam/bernard/msb4100057-s1.rtpj";
        // Do post-processing.
        PantherPostProcessor postProcessor = new PantherPostProcessor();
        postProcessor.postProcess(dbAdaptor, fileAdaptor);
        postProcessor.otherProcesses(dbAdaptor, fileAdaptor);
        
        String destFileName = FIConfiguration.getConfiguration().get("PANTHER_DIR") + "SinglePathwayTest.rtpj";
        fileAdaptor.save(destFileName);
    }
    
    public void testUniAnnotations() throws Exception {
        //String rtpjFileName = PANTHER_DIR + "converted" + File.separator + "Wnt signaling pathway.rtpj";
        String rtpjFileName = "/Users/wgm/Documents/caBIG_R3/datasets/Panther/converted/AllPathways.rtpj";
        PantherPostProcessor annotator = new PantherPostProcessor();
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
                                                  "gk_central_072706",
                                                  "root",
                                                  "macmysql01",
                                                  3306);
        annotator.annotate(rtpjFileName, dbAdaptor);
    }
    
    public void processDavidChEBIMapFiles() throws IOException {
        // Load all ChEBI ids from obo files
        List<List<String>> ids = new ArrayList<List<String>>();
        String fileName = "/Users/wgm/Documents/caBIG_R3/datasets/ChEBI/chebi.obo";
//        [Term]
//        id: CHEBI:16042
//        name: halide anions
//        alt_id: CHEBI:5605
//        alt_id: CHEBI:14384
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        List<String> list = null;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("[Term]")) {
                list = new ArrayList<String>();
                ids.add(list);
            }
            else if (line.startsWith("id")) {
                String id = line.substring(10);
                list.add(id);
            }
            else if (line.startsWith("alt_id")) {
                String id = line.substring(14);
                list.add(id);
            }
        }
        fu.close();
        // Add these lines after ids
        String dir = "/Users/wgm/Documents/gkteam/david/pantherMap/";
        String inFileName = dir + "species_to_chebi_id.txt";
        fu.setInput(inFileName);
        String outFileName = dir + "SpeciesToChEBIId.txt";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFileName);
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            //10-formyl-THF   698
            String[] tokens = line.split("\t");
            // Find ChEBI id
            for (List<String> l : ids) {
                if (l.contains(tokens[1])) {
                    builder.append(tokens[0]).append("\t");
                    for (Iterator<String> it = l.iterator(); it.hasNext();) {
                        builder.append(it.next());
                        if (it.hasNext())
                            builder.append("\t");
                    }
                    break;
                }
            }
            if (builder.length() > 0)
                outFu.printLine(builder.toString());
            else
                outFu.printLine(line); // Some ids might be not in the obo file
            builder.setLength(0);
        }
        fu.close();
        outFu.close();
    }
}
