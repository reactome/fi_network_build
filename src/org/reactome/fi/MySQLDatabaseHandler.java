/*
 * Created on Jul 10, 2012
 *
 */
package org.reactome.fi;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import org.apache.log4j.Logger;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to dump mysql database and convert the dumped file from Innodb to MyIsam.
 * The conversion from Innodb to MyISAM cannot be done for the FI database because of foreign keys
 * constraints imposed by hibernate.
 * @author gwu
 *
 */
public class MySQLDatabaseHandler {
    private static final Logger logger = Logger.getLogger(MySQLDatabaseHandler.class);
    
    public MySQLDatabaseHandler() {
    }
    
    /**
     * Convert original INNODB FI database into myisam.
     * @throws Exception
     */
    @Test
    public void dumpFIDatabase() throws Exception {
        String fiDbName = FIConfiguration.getConfiguration().get("FI_DB_NAME");
        dump(fiDbName);
    }
    
    public void dumpReactomeSourceDatabase() throws Exception {
        String fiDbName = FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME");
        dump(fiDbName);
    }
    
    private void dump(String dbName) throws Exception {
        logger.info("Dump database " + dbName + "...");
        String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + File.separator + dbName + ".sql";
        String command = FIConfiguration.getConfiguration().get("MYSQLDUMP") + 
                         " -u" + FIConfiguration.getConfiguration().get("DB_USER") + 
                         " -p" + FIConfiguration.getConfiguration().get("DB_PWD") + 
                         " " + dbName;
        Process process = Runtime.getRuntime().exec(command);
        InputStream output = process.getInputStream();
        InputStreamReader isr = new InputStreamReader(output);
        BufferedReader bis = new BufferedReader(isr);
        String line = null;
        PrintWriter pr = new PrintWriter(fileName);
        while ((line = bis.readLine()) != null) {
            if (line.contains("ENGINE=InnoDB")) {
//                System.out.println(line);
                line = line.replaceAll("InnoDB", "MyISAM");
            }
            pr.println(line);
        }
        bis.close();
        isr.close();
        output.close();
        pr.close();
    }
    
    @Test
    public void dumpDatabases() throws Exception {
        dumpFIDatabase();
        dumpReactomeSourceDatabase();
    }
    
}
