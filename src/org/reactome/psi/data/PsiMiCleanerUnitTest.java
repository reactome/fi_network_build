/*
 * Created on Apr 18, 2006
 *
 */
package org.reactome.psi.data;

import java.io.File;
import java.lang.reflect.Method;

import junit.framework.TestCase;

import org.apache.log4j.PropertyConfigurator;
import org.hibernate.SessionFactory;
import org.reactome.fi.util.HibernateUtil;
import org.reactome.psi.PsiMiCleaner;

public class PsiMiCleanerUnitTest extends TestCase{
    private String configFileName = "resources/hibernate.cfg.xml";
    private PsiMiCleaner cleaner;
    
    public PsiMiCleanerUnitTest() {   
    }
    
    protected void setUp() throws Exception {
        String log4jFileName = "resources" + File.separator + "log4j.prop";
        PropertyConfigurator.configure(log4jFileName);
        File configFile = new File(configFileName);
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(configFile);
        cleaner = new PsiMiCleaner();
        cleaner.setSessionFactory(sessionFactory);
        super.setUp();
    }
    
    public void testDbReference() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.generateDbReferenceMergeList("dbReference.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time: " + (time2 - time1));
    }
    
    public void testUpdateXref() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.updateXref("dbReference.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time: " + (time2 - time1));
    }
    
    public void testUpdateDbReference() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.cleanUpDbReference("dbReference.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for updateDbReference: " + (time2 - time1));
    }
    
    public void testGenerateXrefList() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.generateXrefMergeList("xref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateXrefList: " + (time2 - time1));
    }
    
    public void testUpdateXrefReferrers() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.updateXrefReferrers("xref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for updateXrefReferrers: " + (time2 - time1));
    }
    
    public void testCleanUpXref() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.cleanUpXref("xref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for cleanUpXref: " + (time2 - time1));
    }
    
    public void testGenerateBibrefMergeList() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.generateBibrefMergeList("bibref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateBibrefMergeList: " + (time2 - time1));
    }
    
    public void testUpdateBibrefReferrers() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.updateBibrefReferrers("bibref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for updateBibrefReferrers: " + (time2 - time1));
    }
    
    public void testCleanUpBibref() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.cleanUpBibref("bibref.txt");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for cleanUpBibref: " + (time2 - time1));
    }
    
    public void testGenerateOpenCVMergeList() throws Exception {
        long time1 = System.currentTimeMillis();
        cleaner.generateOpenCVMergeList();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for generateOpenCVMergeList: " + (time2 - time1));
    }
    
    public void testUpdateOpenCVReferrers() throws Exception {
        testMethod("updateOpenCVReferrers");
    }
    
    public void testCleanUpOpenCV() throws Exception {
        testMethod("cleanUpOpenCV");
    }
    
    public void testGenerateBioSourceMergeList() throws Exception {
        testMethod("generateBioSourceMergeList");
    }
    
    public void testUpdateBioSourceReferrers() throws Exception {
        testMethod("updateBioSourceReferrers");
    }
    
    public void testCleanUpBioSource() throws Exception {
        testMethod("cleanUpBioSource");
    }
    
    public void testgenerateSourceMergeList() throws Exception {
        testMethod("generateSourceMergeList");
    }
    
    public void testUpdateSourceReferrers() throws Exception {
        testMethod("updateSourceReferrers");
    }
    
    public void testCleanUpSource() throws Exception {
        testMethod("cleanUpSource");
    }
    
    private void testMethod(String methodName) throws Exception {
        // Find the method first
        Method[] methods = cleaner.getClass().getMethods();
        Method targetMethod = null;
        for (Method m : methods) {
            if (m.getName().equals(methodName)) {
                targetMethod = m;
                break;
            }
        }
        long time1 = System.currentTimeMillis();
        targetMethod.invoke(cleaner, new Object[]{});
        long time2 = System.currentTimeMillis();
        System.out.printf("Time for %s: %d%n", methodName, (time2 - time1));
    }
    
    public void testAttachUniProt() throws Exception {
        long time1 = System.currentTimeMillis();
        String mapFileName = "/Users/wgm/Documents/caBIG_R3/datasets/BIND/GI2UniMap.txt";
        cleaner.attachUniProtId(mapFileName);
        long time2 = System.currentTimeMillis();
        System.out.println("Time for attachUniProtId: " + (time2 - time1));
    }
}
