/*
 * Created on Mar 1, 2012
 *
 */
package org.reactome.fi;

import java.io.File;

import org.apache.log4j.PropertyConfigurator;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.convert.common.ReactomeProjectDumper;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to save converted pathway databases in the rtpj format 
 * into the Reactome source database.
 * @author gwu
 * @TODO: TO be moved to a better place.
 */
public class ConvertedPathwayDbDumper {
    
    public ConvertedPathwayDbDumper() {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    @Test
    public void dump() throws Exception {
        ReactomeProjectDumper dumper = new ReactomeProjectDumper();
        dumper.setMySQLAdaptor(getDBA());
        String[] fileNames = getFileNames();
        dumper.dumpToDB(fileNames);
    }
    
    private String[] getFileNames() {
        // Keep the order
        String[] fileNames = new String[] {
                FIConfiguration.getConfiguration().get("KEGG_CONVERTED_FILE"),
                FIConfiguration.getConfiguration().get("NATURE_PID_CURATED_CONVERTED"),
                FIConfiguration.getConfiguration().get("NATURE_PID_BIOCARTA_CONVERTED"),
                FIConfiguration.getConfiguration().get("PANTHER_CONVERTED_FILE"),
                FIConfiguration.getConfiguration().get("TRED_CONVERTED_FILE"),
                FIConfiguration.getConfiguration().get("ENCODE_TFF_CONVERTED_FILE")
        };
        // Make sure all file existing
        for (String fileName : fileNames) {
            File file = new File(fileName);
            if (!file.exists())
                throw new IllegalStateException(fileName + " doesn't exist!");
        }
        return fileNames;
    }
    
    private MySQLAdaptor getDBA() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost", 
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"), 
                                            FIConfiguration.getConfiguration().get("DB_USER"), 
                                            FIConfiguration.getConfiguration().get("DB_PWD"));
        return dba;
    }
    
}
