/*
 * Created on Jul 13, 2012
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.SchemaAttribute;
import org.gk.schema.SchemaClass;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to handle files needed by Cytoscape Plug-in and its associated web application.
 * @author gwu
 *
 */
public class CytoscapePlugInFileHandler {
    private final static Logger logger = Logger.getLogger(CytoscapePlugInFileHandler.class);
    
    public CytoscapePlugInFileHandler() {
    }
    
    @Test
    public void generatePlugInFiles() throws Exception {
        validateInteractionTypes();
        copyFiles();
        generateConfigurationFile();
        copyFIHibernateConfigFile();
    }
    
    /**
     * Make sure all interaction types have been registered for annotations. Otherwise, null exception 
     * will be thrown during FI annotating.
     * @throws Exception
     */
    @Test
    public void validateInteractionTypes() throws Exception {
        // Query interaction types used in the database
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                            FIConfiguration.getConfiguration().get("DB_USER"),
                                            FIConfiguration.getConfiguration().get("DB_PWD"));
        Collection<?> interactions = dba.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
        SchemaClass cls = dba.getSchema().getClassByName(ReactomeJavaConstants.Interaction);
        SchemaAttribute att = cls.getAttribute(ReactomeJavaConstants.interactionType);
        dba.loadInstanceAttributeValues(interactions, att);
        Set<String> types = new HashSet<String>();
        for (Iterator<?> it = interactions.iterator(); it.hasNext();) {
            GKInstance interaction = (GKInstance) it.next();
            String type = (String) interaction.getAttributeValue(ReactomeJavaConstants.interactionType);
            if (type != null)
                types.add(type);
        }
        logger.info("Total interaction types: " + types.size());
        // Load listed InteractionTypes in the XML file
        String fileName = "resources/InteractionTypeMapper.xml";
        SAXBuilder builder = new SAXBuilder();
        org.jdom.Document document = builder.build(new File(fileName));
        Element root = document.getRootElement();
        List<?> children = root.getChildren("type");
        Set<String> listedTypes = new HashSet<String>();
        for (Iterator<?> it = children.iterator(); it.hasNext();) {
            Element typeElm = (Element) it.next();
            listedTypes.add(typeElm.getAttributeValue("name"));
        }
        logger.info("Listed interaction types: " + listedTypes.size());
        types.removeAll(listedTypes);
        if (types.size() > 0) {
            logger.error("Interaction types have not been listed: " + types.size());
            System.out.println("Add the following types into file resources/InteractionTypeMapper.xml " +
                                "and make sure reverse types and directions are correct:");
            for (String type : types) {
                // Print out in the XML element format for easy integation
                System.out.println("<type name=\"" + type + "\" reverse=\"" + type + "\" direction=\"-\" />");
            }
            throw new IllegalStateException("Interaction types have not been listed: " + types.size());
        }
    }
    
    @Test
    public void testValidateInteractionTypes() throws Exception {
        PropertyConfigurator.configure("resources/log4j.properties");
        validateInteractionTypes();
    }
    
    private void copyFIHibernateConfigFile() throws IOException {
        String srcFileName = "resources/funcIntHibernate.cfg.xml";
        String destFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + File.separator + "funcIntHibernate.cfg.xml";
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        fu.setOutput(destFileName);
        String line = null;
        while ((line = fu.readLine()) != null) {
            if (line.contains("connection.url")) {
                String dbName = FIConfiguration.getConfiguration().get("FI_DB_NAME");
                line = "<property name=\"connection.url\">jdbc:mysql://localhost:3306/test_" + dbName + "</property>";
            }
            else if (line.contains("connection.username")) {
                line = "<property name=\"connection.username\">authortool</property>";
            }
            else if (line.contains("connection.password")) {
                line = "<property name=\"connection.password\">T001test</property>";
            }
            fu.printLine(line);
        }
        fu.close();
    }
    
    /**
     * Generate a server configuration file based on the output files. Each
     * release will have output files named different because of time-stamps.
     */
    private void generateConfigurationFile() throws IOException {
        // Use a template file to generate the configuration to be used in the sever side
        String srcFileName = "resources/servlet.config.prop";
        String destFileName = FIConfiguration.getConfiguration().get("SERVLET_CONFIG_FILE");
        FileUtility fu = new FileUtility();
        fu.setInput(srcFileName);
        fu.setOutput(destFileName);
        String line = null;
        Map<String, String> keyToFileName = getKeyToFileNameMap();
        while ((line = fu.readLine()) != null) {
            for (String key : keyToFileName.keySet()) {
                if (line.startsWith(key)) {
                    String value = FIConfiguration.getConfiguration().get(keyToFileName.get(key));
                    // Just need the file only
                    File file = new File(value);
                    if (!file.exists())
                        throw new IllegalStateException("Cannot find file: " + value);
                    line = assignValue(line, file.getName());
                    break;
                }
            }
            // Two special cases
            if (line.startsWith("Reactome.src.dbName") || line.startsWith("elv.dbName")) {
                String value = FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME");
                int index = line.indexOf("=");
                line = line.substring(0, index + 1) + "test_" + value; // This name pattern should be followed always
            }
            String year = FIConfiguration.getConfiguration().get("YEAR");
            line = line.replaceAll("caBigR3WebApp", "caBigR3WebApp" + year);
            fu.printLine(line);
        }
        fu.close();
    }
    
    private Map<String, String> getKeyToFileNameMap() {
        Map<String, String> map = new HashMap<String, String>();
        map.put("fi.gene.file", "GENE_FI_FILE_NAME");
        map.put("fi.bigComp.gene.file", "GENE_FI_BIG_COMP_FILE_NAME");
        map.put("fi.pathway.gene.file", "GENE_FI_PATHWAY_FILE_NAME");
        map.put("fi.predicted.gene.file", "GENE_FI_PREDICTED_FILE_NAME");
        map.put("name.to.pathways", "GENE_TO_TOPIC");
        map.put("protein.acc.to.name", "PROTEIN_ACCESSION_TO_NAME_FILE");
        map.put("name.to.reactome.pathways", "GENE_TO_REACTOME_PATHWAYS");
        return map;
    }
    
    private String assignValue(String line, String value) {
        int index = line.indexOf("=");
        return line.substring(0, index + 1) + value; // Add 1 to include "="
    }
    
    private void copyFiles() throws IOException {
        FIConfiguration config = FIConfiguration.getConfiguration();
        String[] srcFileNames = new String[] {
                config.get("GO_DIR") + "gene_association.goa_human",
                config.get("GO_DIR") + "GO.terms_and_ids.txt",
                config.get("KEGG_DIR") + "map_title.tab",
                config.get("KEGG_HSA_KGML_DIR") + "hsa.list",
                "resources" + File.separator + "InteractionTypeMapper.xml",
        };
        String[] targetNames = new String[] {
                "gene_association.goa_human",
                "GO.terms_and_ids.txt",
                "kegg_map_title.tab",
                "kegg_hsa.list",
                "InteractionTypeMapper.xml"
        };
        FileUtility fu = new FileUtility();
        logger.info("Copying files to the results directory ...");
        for (int i = 0; i < srcFileNames.length; i++) {
            String src = srcFileNames[i];
            String dest = config.get("RESULT_DIR") + File.separator + targetNames[i];
            logger.info("Copying " + src);
            fu.copy(new File(src), 
                    new File(dest));
        }
    }
    
}
