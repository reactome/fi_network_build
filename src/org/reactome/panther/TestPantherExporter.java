package org.reactome.panther;

import java.io.File;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.junit.Test;
import org.reactome.convert.common.Converter;
import org.reactome.convert.common.ConverterHandler;
import org.reactome.convert.common.Pair;

public class TestPantherExporter {
    static Logger logger = Logger.getLogger(TestPantherExporter.class);
    
    /**
     * @param args
     */
    public static void main(String[] args) {
        if (true)
            convertSingleFile();
        else
            convertBatch();
    }
    
    private static String getTime() {
        Calendar myCal = Calendar.getInstance();
        return myCal.get(Calendar.HOUR_OF_DAY) + ":" + myCal.get(Calendar.MINUTE) + ":" + myCal.get(Calendar.SECOND);
    }
    
    private static void log(String logText) {
        System.out.print("["+getTime()+"] "+logText);
    }
    
    private static void convertBatch() {
        try {
            PropertyConfigurator.configure("resources/log4j.properties");
            
            String destDir = "E:\\panther\\results\\";
            String PANTHER_DIR = "E:\\panther\\SBML\\";
            
            ConverterHandler handler = ConverterHandler.getInstance();
            
            XMLFileAdaptor reactomeAdaptor;
            MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost", "test_reactome_23_pathway_diagram", "myUsr",
                                                      "myUsr", 3306);
            
            File[] sbmlFiles = new File(PANTHER_DIR).listFiles();
            List<File> files = new ArrayList<File>();
            
            for (File myFile : sbmlFiles) {
                if (myFile.getName().endsWith(".xml"))
                    files.add(myFile);
            }
            
            for (int i = 0; i < files.size(); i++) {
                File file = files.get(i);
                String fileName = file.getName();
                try {
                    log((i+1)+"/"+files.size()+": Converting "+fileName+"... ");
                    Pair convPair = (Pair)handler.autoDetect(file.getAbsolutePath());
                    Converter cdConverter = handler.getConverter(convPair);
                    reactomeAdaptor = new XMLFileAdaptor();
                    cdConverter.setFileAdaptor(reactomeAdaptor);
                    cdConverter.setDatabaseAdaptor(dbAdaptor);
                    GKInstance pathway = cdConverter.convert(file.getAbsolutePath());
                    String destFileName = destDir + fileName.substring(0, fileName.lastIndexOf(".")) + ".rtpj";
                    cdConverter.save(destFileName);
                    System.out.println("done.");
                } catch (Exception ex) {
                    System.out.println("failed ("+ex.getClass().getName().toString()+")");
                }	
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    private static void convertSingleFile() {
        try {
            PropertyConfigurator.configure("resources/log4j.properties");
            
            //String PANTHER_DIR = "E:\\panther\\";
            String PANTHER_DIR = "E:\\panther\\SBML\\";
            String pantherModelName = "Synaptic_vesicle_trafficking.xml";
            String fileName = pantherModelName.replace(".xml", "");
            pantherModelName = PANTHER_DIR + pantherModelName;
            
            ConverterHandler handler = ConverterHandler.getInstance();
            
            Pair<String, String> convPair = handler.autoDetect(pantherModelName);
            
            Converter cdConverter = handler.getConverter(convPair);
            
            XMLFileAdaptor reactomeAdaptor = new XMLFileAdaptor();
            cdConverter.setFileAdaptor(reactomeAdaptor);
            
            MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
                                                      "test_reactome_23_pathway_diagram", 
                                                      "myUsr",
                                                      "myUsr", 
                                                      3306);
            cdConverter.setDatabaseAdaptor(dbAdaptor);
            
            GKInstance pathway = cdConverter.convert(pantherModelName);
            
            String destFileName = "E:\\Panther\\results\\" + fileName + ".rtpj";
            cdConverter.save(destFileName);
            
            System.out.println("done");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}

