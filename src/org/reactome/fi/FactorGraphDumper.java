/*
 * Created on Feb 11, 2015
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Test;
import org.reactome.factorgraph.FactorGraph;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.pathway.factorgraph.PathwayToFactorGraphConverter;
import org.reactome.pathway.factorgraph.ReactomePathwayFGRunner;
import org.reactome.r3.util.FileUtility;
import org.reactome.r3.util.JAXBBindableList;

/**
 * This simple class is used to convert all pathway diagrams, which have entities drawn,
 * into factor graphs and them dump them into a zipped XML file by using JAXB.
 * @author gwu
 *
 */
public class FactorGraphDumper {
    private static Logger logger = Logger.getLogger(FactorGraphDumper.class);
    
    /**
     * Default constructor.
     */
    public FactorGraphDumper() {
        FileUtility.initializeLogging();
    }
    
    /**
     * This is the actual method to perform dumping.
     * @throws Exception
     */
    @Test
    public void dump() throws Exception {
        MySQLAdaptor dba = new MySQLAdaptor("localhost", // Always assume it is at the localhost
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                            FIConfiguration.getConfiguration().get("DB_USER"), 
                                            FIConfiguration.getConfiguration().get("DB_PWD"));
        ReactomePathwayFGRunner runner = new ReactomePathwayFGRunner();
        runner.setAdaptor(dba);
        List<GKInstance> pathways = runner.getPathwayList();
        logger.info("Total pathways for converting into factor graphs: " + pathways.size());
        PathwayToFactorGraphConverter converter = new PathwayToFactorGraphConverter();
        long time1 = System.currentTimeMillis();
        int count = 1;
        JAXBBindableList<FactorGraph> fgList = new JAXBBindableList<FactorGraph>();
        for (GKInstance pathway : pathways) {
            logger.info(count +": " + pathway);
            FactorGraph fg = converter.convertPathway(pathway);
            fgList.getList().add(fg);
            count ++;
//            if (count == 4)
//                break;
        }
        long time2 = System.currentTimeMillis();
        logger.info("Total time: " + (time2 - time1) + " ms");
        // Export into a file
        // Have to list both classes here.
        JAXBContext context = JAXBContext.newInstance(JAXBBindableList.class, FactorGraph.class);
        Marshaller marshaller = context.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        String fileName = FIConfiguration.getConfiguration().get("FACTOR_GRAPH_FILE");
        FileOutputStream fos = new FileOutputStream(fileName);
        ZipOutputStream zos = new ZipOutputStream(fos);
        // The last part should be zip in the file name
        int index = fileName.lastIndexOf(".");
        fileName = fileName.substring(0, index);
        File file = new File(fileName);
        ZipEntry entry = new ZipEntry(file.getName());
        zos.putNextEntry(entry);
        marshaller.marshal(fgList, zos);
        zos.close();
    }
    
    @Test
    public void testRead() throws Exception {
        JAXBContext context = JAXBContext.newInstance(JAXBBindableList.class, FactorGraph.class);
        String fileName = FIConfiguration.getConfiguration().get("FACTOR_GRAPH_FILE");
        File file = new File(fileName);
        Unmarshaller unmarshaller = context.createUnmarshaller();
        FileInputStream fis = new FileInputStream(file);
        ZipInputStream zis = new ZipInputStream(fis);
        ZipEntry entry = zis.getNextEntry(); // Have to call this method
        JAXBBindableList<FactorGraph> list = (JAXBBindableList<FactorGraph>) unmarshaller.unmarshal(zis);
        logger.info("Total factor graphs: " + list.getList().size());
        for (FactorGraph fg : list.getList()) {
            logger.info(fg.getName());
        }
    }
    
}
