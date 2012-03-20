/*
 * Created on Aug 21, 2006
 *
 */
package org.reactome.psi;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.PropertyConfigurator;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.HibernateUtil;

public class PsiMiToReactomeUnitTest {
    private MySQLAdaptor dba = null;
    
    public PsiMiToReactomeUnitTest() {
        String log4jFileName = "resources" + File.separator + "log4j.properties";
        PropertyConfigurator.configure(log4jFileName);
    }
    
    private MySQLAdaptor getMySQLAdaptor() throws Exception {
        if (dba != null)
            return dba;
        dba = new MySQLAdaptor("localhost",
                               FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                               FIConfiguration.getConfiguration().get("DB_USER"),
                               FIConfiguration.getConfiguration().get("DB_PWD"),
                               3306);
        return dba;
    }
    
    public void postProcessBIND() throws Exception {
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        String fileName = "/Users/wgm/Documents/tmp/BIND.rtpj";
        fileAdaptor.setSource(fileName);
        MySQLAdaptor dbAdaptor = getMySQLAdaptor();
        BINDPsiMiToReactomePostProcessor processor = new BINDPsiMiToReactomePostProcessor();
        long time1 = System.currentTimeMillis();
        processor.postProcess(dbAdaptor, fileAdaptor);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for post-processing: " + (time2 - time1));
        fileAdaptor.save("/Users/wgm/Documents/tmp/BIND_POST.rtpj");
    }
    
    public void convertBIND() throws Exception {
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/BIND/";
        String fileName = dirName + "psi/taxid9606.1.psi.xml";
        List<String> fileNameList = new ArrayList<String>(1);
        fileNameList.add(fileName);
        String outFileName = dirName + "BIND.rtpj";
        // For BIND, want to do postProcess after all entries are converted since
        // there are a lot of entries in one file.
        PsiMiToReactomePostProcessor postProcessor = new BINDPsiMiToReactomePostProcessor();
        convert(fileNameList,
                outFileName,
                postProcessor);
    }
    
    private void convert(List<String> fileNameList,
                         String outputFileName,
                         PsiMiToReactomePostProcessor postProcessor) throws Exception {
        PsiMiToReactomConverter converter = new PsiMiToReactomConverter();
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        converter.setReactomeAdaptor(fileAdaptor);
        MySQLAdaptor dbAdaptor = getMySQLAdaptor();
        converter.setMySQLAdaptor(dbAdaptor);
        PsiMiLoader loader = new PsiMiLoader();
        loader.setPostProcessor(converter);
        for (String fileName : fileNameList) {
            long startTime = System.currentTimeMillis();
            System.out.println("Starting parsing file: " + fileName);
            loader.parse(fileName);
            long endTime = System.currentTimeMillis();
            System.out.println("Total parsing time: " + (endTime - startTime));
        }
        System.out.println("Starting post-processing...");
        long time1 = System.currentTimeMillis();
        postProcessor.postProcess(dbAdaptor, fileAdaptor);
        long time2 = System.currentTimeMillis();
        System.out.println("Total post-processing time: " + (time2 - time1));
        fileAdaptor.save(outputFileName);           
    }
    
    public void convertHPRD() throws Exception {
        //String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/HPRD/PSI-MI/";
        //String fileName = dirName + "HPRD_SINGLE_PSIMI_060106_FIXED.xml";
        String dirName = FIConfiguration.getConfiguration().get("HPRD_DIR");
        String fileName = dirName + "HPRD_SINGLE_PSIMI_090107.xml";
        //String fileName = dirName + "HPRD_SINGLE_PSIMI_010107.xml";
        List<String> fileNameList = new ArrayList<String>(1);
        fileNameList.add(fileName);
        //String outputFileName = dirName + "HPRD060106.rtpj";
        String outputFileName = dirName + "HPRD.rtpj";
        PsiMiToReactomePostProcessor postProcessor = new HPRDPsiMiToReactomePostProcessor();
        convert(fileNameList,
                outputFileName,
                postProcessor);
    }
    
    /**
     * This method is used to convert a PSI-MI file to a Reactome curator tool project.
     * @throws Exception
     */
    public void convertBioGrid() throws Exception {
        String dirName = FIConfiguration.getConfiguration().get("BIOGRID_DIR");
        String fileName = dirName + "BIOGRID-ORGANISM-Homo_sapiens-2.0.50.psi25.xml";
        List<String> fileNameList = new ArrayList<String>();
        fileNameList.add(fileName);
        String outFileName = dirName + "BioGrid.rtpj";
        PsiMiToReactomePostProcessor postProcessor = new BioGridPsiMiToReactomePostProcessor();
        convert(fileNameList,
                outFileName,
                postProcessor);
    }
    
    @Test
    public void convertIntAct() throws Exception {
        String dirName = FIConfiguration.getConfiguration().get("INTACT_HUMAN_DIR");
        List<String> fileList = new ArrayList<String>();
        File[] files = new File(dirName).listFiles();
        for (File file : files) {
            String name = file.getName();
            if (name.endsWith(".xml") &&
                !name.endsWith("negative.xml") && // Don't want to consider negative
                name.startsWith("human")) // Use human PPIs only
                fileList.add(file.getAbsolutePath());
//            if (fileList.size() > 0)
//                break;
        }
//        fileList.add(dirName + "human_small-01.xml");
//        fileList.add(dirName + "human_small-06.xml");
        System.out.println("Total files: " + fileList.size());
        String outFileName = dirName + "IntAct.rtpj";
        PsiMiToReactomePostProcessor postProcessor = new PsiMiToReactomePostProcessor();
        postProcessor.setDataSourceName("IntAct");
        postProcessor.setDataSourceUrl("http://www.ebi.ac.uk/intact/site/");
        
//        // Just for test
//        fileList.clear();
//        fileList.add(dirName + "human_18.xml");
        
        convert(fileList,
                outFileName,
                postProcessor);
    }
    
//    public void removeDatabseIdentifierForManuel() throws Exception {
//        String dirName = "/Users/wgm/Documents/Manuel/";
//        String source = dirName + "Manuel_1.rtpj";
//        String dest = dirName + "Manuel_2.rtpj";
//        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
//        fileAdaptor.setSource(source);
//        // Want to delete all DatabaseIdentifiers, which cannot provide any more information 
//        // for Interactions.
//        try {
//            Collection databaseIdentifiers = fileAdaptor.fetchInstancesByClass(ReactomeJavaConstants.Interaction);
//            for (Iterator it = databaseIdentifiers.iterator(); it.hasNext();) {
//                GKInstance instance = (GKInstance) it.next();
//                GKInstance crossRef = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.crossReference);
//                if (crossRef != null) {
//                    instance.removeAttributeValueNoCheck(ReactomeJavaConstants.crossReference, crossRef);
//                    fileAdaptor.removeFromClassMap(crossRef);
//                }
//            }
//            fileAdaptor.save(dest);
//        }
//        catch(Exception e) {
//            e.printStackTrace();
//        }
//    }
    
    public void testConvertFromFile() throws Exception {
        PsiMiLoader loader = new PsiMiLoader();
        PsiMiToReactomConverter converter = new PsiMiToReactomConverter();
        XMLFileAdaptor fileAdaptor = new XMLFileAdaptor();
        converter.setReactomeAdaptor(fileAdaptor);
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
                                                  "gk_central_072706",
                                                  "root",
                                                  "macmysql01",
                                                  3306);
        converter.setMySQLAdaptor(dbAdaptor);
        loader.setPostProcessor(converter);
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/";
        String[] filesNames = new String[] {
                //"/Users/wgm/Documents/caBIG_R3/datasets/BIND/psi/SmallSample.xml"
                //"/Users/wgm/Documents/caBIG_R3/datasets/IntAct/human_small-10_negative.xml",
                //dirName + "IntAct/human_rual-2005-2_01.xml"
                //dirName + "IntAct/human_small-05.xml"
                //"/Users/wgm/Documents/caBIG_R3/datasets/IntAct"
                //dirName + "HPRD/PSI-MI/HPRD_SINGLE_PSIMI_060106_FIXED.xml"
        };
        List<String> fileList = new ArrayList<String>();
        File file = null;
        for (String name : filesNames) {
            file = new File(name);
            if (file.isDirectory()) {
                File[] files = file.listFiles();
                for (File f : files) {
                    if (f.getName().endsWith(".xml")) {
                        fileList.add(f.getAbsolutePath());
                    }
                }
            }
            else
                fileList.add(name);
        }
        for (String fileName : fileList) {
            long startTime = System.currentTimeMillis();
            System.out.println("Starting parse file: " + fileName);
            loader.parse(fileName);
            //Entry entry = loader.getLastEntry();
            //printProperties(entry, "");
            long endTime = System.currentTimeMillis();
            System.out.println("Total parsing time: " + (endTime - startTime));
        }
        // For BIND, want to do postProcess after all entries are converted since
        // there are a lot of entries in one file.
        PsiMiToReactomePostProcessor postProcessor = new PsiMiToReactomePostProcessor();
        System.out.println("Starting post-processing...");
        long time1 = System.currentTimeMillis();
        // This is for HPRD
        postProcessor.setDataSourceName("HPRD");
        postProcessor.setDataSourceUrl("http://www.hprd.org");
        postProcessor.postProcess(dbAdaptor, fileAdaptor);
        long time2 = System.currentTimeMillis();
        System.out.println("Total post-processing time: " + (time2 - time1));
        fileAdaptor.save("/Users/wgm/Documents/tmp/HPRD.rtpj");
    }
    
    public void convert() throws Exception {
        long time1 = System.currentTimeMillis();
        String configFileName = "resources/hibernate.cfg.xml";
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(new File(configFileName));
        PsiMiToReactomConverter converter = new PsiMiToReactomConverter();
        converter.setPostProcessor(new PsiMiToReactomePostProcessor());
        converter.setReactomeAdaptor(new XMLFileAdaptor());
        MySQLAdaptor dbAdaptor = new MySQLAdaptor("localhost",
                                                  "gk_central_072706",
                                                  "root",
                                                  "macmysql01",
                                                  3306);
        converter.setMySQLAdaptor(dbAdaptor);
        converter.setSessionFactory(sessionFactory);  
        converter.convert();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for converting: " + (time2 - time1));
        XMLFileAdaptor fileAdaptor = converter.getReactomeAdaptor();
        fileAdaptor.save("/Users/wgm/Documents/tmp/BIND_IntAct.rtpj");
    }
//    
//    public static void main(String[] args) {
//        try {
//            PsiMiToReactomeUnitTest test = new PsiMiToReactomeUnitTest();
//            test.testConvertFromFile();
//        }
//        catch(Exception e) {
//            e.printStackTrace();
//        }
//    }
}
