/*
 * Created on Jan 16, 2012
 *
 */
package org.reactome.b2rPostProcessor;

import java.io.File;

import org.apache.log4j.PropertyConfigurator;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.junit.Test;
import org.reactome.biopax.BioPAXToReactomeConverter;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This simple class is used to run NCI-PID converting from BioPAX level 2 to the Reactome data format.
 * @author gwu
 *
 */
public class NciPIDConverterRunner {
    
    public NciPIDConverterRunner() {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    private MySQLAdaptor getDBA() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                            FIConfiguration.getConfiguration().get("DB_USER"), 
                                            FIConfiguration.getConfiguration().get("DB_PWD"));
        return dba;
    }
    
    @Test
    public void runConvertOfCurated() throws Exception {
        convert(FIConfiguration.getConfiguration().get("NATURE_PID_CURATED"),
                FIConfiguration.getConfiguration().get("NATURE_PID_CURATED_CONVERTED"),
                false);
    }
    
    @Test
    public void runConvertOfBioCarta() throws Exception {
        // A GO term to id file should be created for this covnerting
        FileUtility fu = new FileUtility();
        String srcName = FIConfiguration.getConfiguration().get("GO_DIR") + "GO.terms_and_ids.txt";
        File file = new File(srcName);
        if (!file.exists())
            throw new IllegalStateException("GO.terms_and_ids.txt doesn't exist in the GO_DIR! Please download it from http://www.geneontology.org/doc/GO.terms_and_ids");
        File dest = new File("resources" + File.separator + "GO.terms_and_ids.txt");
        fu.copy(file, dest);
        convert(FIConfiguration.getConfiguration().get("NATURE_PID_BIOCARTA"),
                FIConfiguration.getConfiguration().get("NATURE_PID_BIOCARTA_CONVERTED"),
                true);
    }
    
    private void convert(String srcFileName, 
                         String destFileName,
                         boolean isBioCarta) throws Exception {
        BioPAXToReactomeConverter converter = new BioPAXToReactomeConverter();
        MySQLAdaptor dba = getDBA();
        converter.setReactomeDB(dba);
        if (isBioCarta) {
            NciPIDBToRPostProcessor postProcessor = (NciPIDBToRPostProcessor) converter.getMapperFactory().getPostProcessor();
            postProcessor.setDataSourceName("BioCarta - Imported by PID");
        }
        long time1 = System.currentTimeMillis();
        converter.convert(srcFileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for converting: " + (time2 - time1));
        XMLFileAdaptor reactomeAdaptor = converter.getReactomeModel();
        reactomeAdaptor.save(destFileName);
    }
    
}
