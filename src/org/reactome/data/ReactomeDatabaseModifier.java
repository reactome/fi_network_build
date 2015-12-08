/*
 * Created on Dec 5, 2013
 *
 */
package org.reactome.data;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to make some changes to the Reactome database.
 * @author gwu
 *
 */
public class ReactomeDatabaseModifier {
    private static Logger logger = Logger.getLogger(ReactomeDatabaseModifier.class);
    private FIConfiguration config = FIConfiguration.getConfiguration();
    
    /**
     * Default constructor.
     */
    public ReactomeDatabaseModifier() {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    /**
     * Since the slice database contains human UniProt used by pathways only. However,
     * for the FI network, we want to have all human UniProt in the database for imported
     * interactions. This method is used to copy these missed human UniProts to the 
     * Reactome plus i database.
     * @throws Exception
     */
    @Test
    public void copyHumanReferenceGeneProducts() throws Exception {
        MySQLAdaptor srcDBA = new MySQLAdaptor("localhost", 
                                               config.get("REACTOME_GK_CENTRAL_DB_NAME"),
                                               config.get("DB_USER"),
                                               config.get("DB_PWD"));
        
        MySQLAdaptor targetDBA = new MySQLAdaptor("localhost", // always 
                                                  config.get("REACTOME_SOURCE_DB_NAME"),
                                                  config.get("DB_USER"), 
                                                  config.get("DB_PWD"));
        
        GKInstance human = srcDBA.fetchInstance(48887L);
        @SuppressWarnings("unchecked")
        Collection<GKInstance> humanRefGeneProducts = srcDBA.fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceGeneProduct,
                                                                                      ReactomeJavaConstants.species,
                                                                                      "=",
                                                                                      human);
        logger.info("Total human ReferenceGeneProduct in the source database: " + humanRefGeneProducts.size());
        // Check if any instance is in the targetDBA
        try {
            targetDBA.startTransaction();
            int total = 0;
            for (GKInstance inst : humanRefGeneProducts) {
                GKInstance targetInst = targetDBA.fetchInstance(inst.getDBID());
                if (targetInst != null) {
                    if (targetInst.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct))
                        continue;
                    logger.error(inst + " has different schema in the target database!");
                    continue;
                }
                logger.info("Copy " + inst);
                targetDBA.storeInstance(inst, true);    
                total ++;
            }
            targetDBA.commit();
            logger.info("Total copied: " + total);
        }
        catch(Exception e) {
            logger.error(e);
            targetDBA.rollback();
        }
    }
    
    /**
     * Use this method to change myisam to innodb.
     * @throws Exception
     */
    @Test
    public void changeMyISAMToInnodb() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost", // Alwasy 
                                            config.get("REACTOME_SOURCE_DB_NAME"),
                                            config.get("DB_USER"), 
                                            config.get("DB_PWD"));
        Connection conn = dba.getConnection();
        List<String> tables = getTables(conn);
        dropFullTextIndex(tables, conn);
        changeToInnoDB(tables, conn);
    }
    
    private List<String> getTables(Connection conn) throws Exception {
        List<String> tables = new ArrayList<String>();
        Statement statement = conn.createStatement();
        ResultSet resultset = statement.executeQuery("Show Tables");
        while (resultset.next()) {
            String name = resultset.getString(1);
            tables.add(name);
        }
        resultset.close();
        statement.close();
        return tables;
    }
    
    /**
     * Since full text index is not supported in InnoDB. Existing Indices should be dropped first.
     * @param tables
     * @param conn
     * @throws Exception
     */
    private void dropFullTextIndex(List<String> tables, Connection conn) throws Exception {
        Statement statement = conn.createStatement();
        for (String table : tables) {
            ResultSet result = statement.executeQuery("SHOW INDEX FROM " + table);
            while(result.next()) {
                String indexType = result.getString(11);
                if (indexType != null && indexType.equals("FULLTEXT")) {
                    String keyName = result.getString(3);
                    logger.info("Drop full text index: " + table + "." + keyName);
                    conn.createStatement().execute("ALTER TABLE " + table + " DROP INDEX " + keyName);
                }
            }
            result.close();
        }
        statement.close();
    }
    
    private void changeToInnoDB(List<String> tables, Connection conn) throws Exception {
        Statement statement = conn.createStatement();
        for (String table : tables) {
            logger.info("Alter Table to InnoDB: " + table);
            statement.execute("ALTER TABLE " + table + " ENGINE=InnoDB");
        }
        statement.close();
    }
}
