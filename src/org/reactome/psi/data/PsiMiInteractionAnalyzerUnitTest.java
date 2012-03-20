/*
 * Created on Apr 27, 2006
 *
 */
package org.reactome.psi.data;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;

import org.apache.log4j.PropertyConfigurator;
import org.hibernate.Query;
import org.hibernate.ScrollableResults;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.fi.util.HibernateUtil;
import org.reactome.psi.*;

public class PsiMiInteractionAnalyzerUnitTest extends TestCase {
    private String configFileName = "resources/hibernate.cfg.xml";
    //private String psiSchemaFileName = "/Users/wgm/Documents/caBIG_R3/datasets/PSI-MI/MIF25.xsd";
    private File configFile = new File(configFileName);
    private PsiMiInteractionAnalyzer analyzer;
    
    public PsiMiInteractionAnalyzerUnitTest() {   
    }
    
    protected void setUp() throws Exception {
        String log4jFileName = "resources" + File.separator + "log4j.properties";
        PropertyConfigurator.configure(log4jFileName);
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(configFile);
        analyzer = new PsiMiInteractionAnalyzer(sessionFactory);
        super.setUp();
    }
    
    public void testExtractInteractions() throws Exception {
        long time1 = System.currentTimeMillis();
        analyzer.extractInteractions("results/interaction/HPRDInteractions.xml");
        long time2 = System.currentTimeMillis();
        System.out.println("Time for extractInteractions: " + (time2 - time1));
    }

    public void testAdd() throws Exception {
        BaseLocation baseLocation1 = new BaseLocation();
        BaseLocation.Interval interval1 = new BaseLocation.Interval();
        interval1.setBegin(100);
        interval1.setEnd(110);
        BaseLocation.Interval interval2 = new BaseLocation.Interval();
        interval2.setBegin(500);
        interval2.setEnd(505);
        baseLocation1.setBeginInterval(interval1);
        baseLocation1.setEndInterval(interval2);
        
        PsiMiDBAdaptor dbAdaptor = new PsiMiDBAdaptor();
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(configFile);
        dbAdaptor.setSessionFactory(sessionFactory);
        dbAdaptor.save(baseLocation1, null);
        
        sessionFactory.close();
    }
    
    public void testList() throws Exception {
        Session session = prepareSession();
        Query query = session.createQuery("from " + DbReference.class.getName());
        ScrollableResults scrollable = query.scroll();
        long time1 = System.currentTimeMillis();
        StringBuilder builder = new StringBuilder();
        int totalMerge = 0;
        Query query1 = session.createQuery("from org.reactome.psi.DbReference r where " +
                "r.db = :db and " +
                "r.id = :id and " +
                "r.secondary = :secondary and " +
                "r.version = :version and " +
                "r.refType = :refType");
        int c = 0;
        while (!scrollable.isLast()) {
            scrollable.next();
            DbReference dbReference = (DbReference) scrollable.get(0);
            c ++;
            query1.setProperties(dbReference);
            System.out.println("Query dbreference: " + c);
            List list = query1.list();
            totalMerge += list.size();
        }
        long time2 = System.currentTimeMillis();
        System.out.printf("Total merge: %d%n", totalMerge);
        System.out.println("Total Time: " + (time2 - time1));
        FileOutputStream fos = new FileOutputStream("entrez.txt");
        PrintStream ps = new PrintStream(fos);
        ps.print(builder.toString());
        ps.close();
        fos.close();
        //List list = session.createQuery("from org.reactome.psi.Names").list();
        //System.out.println("Size: " + list.size());
        //for (Iterator it = list.iterator(); it.hasNext();) {
          //  Names entry = (Names) it.next();
          //  System.out.printf("%s:%s:%s%n", entry.getShortLabel(), entry.getFullName(), entry.getAlias());
        //}
        session.getTransaction().commit();
    }

    private Session prepareSession() {
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(configFile);
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        return session;
    }
    
    public void testViewInteraction() throws Exception {
        long time1 = System.currentTimeMillis();
        Session session = prepareSession();
        Query query = session.createQuery("from " + Interaction.class.getName());
        int start = 0;
        int max = 2000;
        query.setFirstResult(start);
        query.setMaxResults(max);
        List list = query.list();
        Interaction interaction = null;
        FileWriter fileWriter = new FileWriter("results/interaction/interactions.txt");
        PrintWriter printWriter = new PrintWriter(fileWriter);
        StringBuilder stringBuilder = new StringBuilder();
        List<DbReference> uniXrefList = new ArrayList<DbReference>();
        Set<String> expNames = new HashSet<String>();
        Set<String> interactorTypeNames = new HashSet<String>();
        while (list.size() > 0) {
            for (Iterator it = list.iterator(); it.hasNext();) {
                interaction = (Interaction) it.next();
                stringBuilder.append(interaction.getDbId());
                List<Participant> participants = interaction.getParticipantList();
                for (Participant p : participants) {
                    Interactor interactor = p.getInteractor();
                    // Get xref
                    Xref xref = interactor.getXref();
                    if (xref == null) {
                        // Most of them are small molecules
                        //System.out.println("Interactor has not xref: " + 
                          //      interactor.getDbId() + " " + 
                            //    interactor.getNames().getShortLabel());
                        continue;
                    }
                    DbReference uniRef = getUniProtKBRef(xref);
                    if (uniRef != null) {  
                        uniXrefList.add(uniRef);
                    }
                    // Get InteractorType
                    OpenCV type = interactor.getInteractorType();
                    if (type != null)
                        interactorTypeNames.add(type.getNames().getShortLabel());
                }
                if (uniXrefList.size() < 2)
                    continue; // Nothing to output
                for (DbReference dbRef : uniXrefList) {
                    stringBuilder.append("\t").append(dbRef.getId());
                }
                // Check experiments
                List<Experiment> experiments = interaction.getExperimentList();
                if (experiments != null) {
                    for (Experiment exp : experiments) {
                        OpenCV openCV = exp.getInteractionDetectionMethod();
                        if (openCV != null) {
                            expNames.add(openCV.getNames().getShortLabel());
                        }
                    }
                    for (String name : expNames)
                        stringBuilder.append("\t").append(name);
                }
                //Check InteractionType
                for (String type : interactorTypeNames)
                    stringBuilder.append("\ttype:").append(type);
                printWriter.println(stringBuilder.toString());
                stringBuilder.setLength(0);
                uniXrefList.clear();
                interactorTypeNames.clear();
            }
            session.clear(); // To empty cache to keep the memory usage small.
            start += list.size();
            query.setFirstResult(start);
            list = query.list();
        }
        printWriter.close();
        fileWriter.close();
        long time2 = System.currentTimeMillis();
        System.out.println("Time for view: " + (time2 - time1));
    }
    
    public void testPsiMiLoader() throws Exception {
        PsiMiLoader loader = new PsiMiLoader();
        PsiMiDBAdaptor dbAdaptor = new PsiMiDBAdaptor();
        SessionFactory sessionFactory = HibernateUtil.getSessionFactory(configFile);
        dbAdaptor.setSessionFactory(sessionFactory);
        loader.setPostProcessor(dbAdaptor);
        String dirName = "/Users/wgm/Documents/caBIG_R3/datasets/";
        String[] filesNames = new String[] {
                //"/Users/wgm/Documents/caBIG_R3/datasets/BIND/psi/SmallSample.xml"
                //"/Users/wgm/Documents/caBIG_R3/datasets/BIND/psi/taxid9606.1.psi.xml",
                //"/Users/wgm/Documents/caBIG_R3/datasets/IntAct/human_small-10_negative.xml",
                //dirName + "IntAct/human_rual-2005-2_01.xml"
                //"/Users/wgm/Documents/caBIG_R3/datasets/IntAct"
                dirName + "HPRD/PSI-MI/HPRD_SINGLE_PSIMI_060106_FIXED.xml"
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
    }
    
    public void testPsiMiLoaderParse() throws Exception {
        PsiMiLoader loader = new PsiMiLoader();
        String[] filesNames = new String[] {
                "/Users/wgm/Documents/caBIG_R3/datasets/BIND/psi/SmallSample.xml",
                "/Users/wgm/Documents/caBIG_R3/datasets/BIND/psi/taxid9606.1.psi.xml",
                "/Users/wgm/Documents/caBIG_R3/datasets/IntAct/human_small-10_negative.xml",
                "/Users/wgm/Documents/caBIG_R3/datasets/IntAct/human_small-10.xml"
        };
        for (String fileName : filesNames) {
            long startTime = System.currentTimeMillis();
            loader.parse(fileName);
            //printProperties(entry, "");
            long endTime = System.currentTimeMillis();
            System.out.printf("Total parsing time for %s: %s%n", fileName, (endTime - startTime));
        }
    }
    
    private void printProperties(Object entry, String indent) throws Exception {
        printProperties(entry, indent, System.out);
    }
    
    private void printProperties(Object entry, String indent, PrintStream os) throws Exception {
        if (entry == null)
            return;
        BeanInfo beanInfo = Introspector.getBeanInfo(entry.getClass());
        PropertyDescriptor[] propertyDesc = beanInfo.getPropertyDescriptors();
        os.print(indent);
        os.println(entry.toString());
        for (PropertyDescriptor pd : propertyDesc) {
            Method method = pd.getReadMethod();
            if (method == null)
                continue; // Might be write-only, e.g., names in OpenCVExperimentalWrapper.
            Object value = method.invoke(entry, new Object[]{});
            if (value == null)
                continue;
            os.print(indent);
            os.print(pd.getName() + ": ");
            String fullName = value.getClass().getName();
            if (fullName.startsWith("org.reactome.psi."))
                printProperties(value, indent + "\t", os);
            else if (value instanceof List) {
                List list = (List) value;
                if (list.size() == 0) {
                    os.print("\n");
                    continue;
                }
                for (Iterator it = list.iterator(); it.hasNext();) {
                    printProperties(it.next(), indent + "\t", os);
                }
            }
            else
                os.println(value.toString());
        }
    }
    
    public void testPsiMiModel() throws Exception {
        PsiMiModel constants = new PsiMiModel();
        //constants.setPsiMiSchema(psiSchemaFileName);
        Map<String, Class> elementTypeMap = constants.getElementTypeMap();
        Set<String> keys = elementTypeMap.keySet();
        for (String key : keys) {
            Class type = elementTypeMap.get(key);
            System.out.printf("%s->%s%n", key, type);
        }
    }
    
    public void testMerge() throws Exception {
        Session session = prepareSession();
        long time1 = System.currentTimeMillis();
        Query query = session.createQuery("from org.reactome.psi.DbReference");
        Connection connection = session.connection();
        PreparedStatement statement = connection.prepareStatement("SELECT dbId FROM Xref where primaryRef=? or secondaryRef=?");
        int first = 0;
        int max = 20000;
        query.setFirstResult(first);
        query.setMaxResults(max);
        List list = query.list();
        int total = 0;
        Set<DbReference> set = new HashSet<DbReference>();
        int repeat = 0;
        int xrefRepeat = 0;
        while (list.size() > 0) {
            total += list.size();
            first += list.size();
            System.out.println("First: " + first);
            for (Iterator it = list.iterator(); it.hasNext();) {
                DbReference ref = (DbReference) it.next();
                if (set.contains(ref)) {
                    repeat ++;
                    statement.setLong(1, ref.getDbId());
                    statement.setLong(2, ref.getDbId());
                    ResultSet result = statement.executeQuery();
                    while (result.next())
                        xrefRepeat ++;
                }
                else
                    set.add(ref);
            }
            query.setFirstResult(first);
            list = query.list();
            session.clear();
            set.clear();
        }
        long time2 = System.currentTimeMillis();
        System.out.println("Total: " + total);
        System.out.println("Total Time: " + (time2 - time1));
        System.out.println("Repeat: " + repeat);
        System.out.println("XRefRepeat: " + xrefRepeat);
    }
    
    private DbReference getUniProtKBRef(Xref xref) {
        DbReference dbRef = xref.getPrimaryRef();
        if (dbRef.getDb().startsWith("uniprot"))
            return dbRef;
        List<DbReference> dbRefList = xref.getSecondaryRefList();
        if (dbRefList == null || dbRefList.size() == 0)
            return null;
        for (DbReference ref : dbRefList) {
            if (ref.getDb().startsWith("uniprot"))
                return ref;
        }
        return null;
    }
    
    public void grepPrimaryDBReferenceInInteractor() throws Exception {
        Map<String, Set<String>> gi2UniMap = new HumanPsiMiInteractionAnalyzer().generateGI2UniMap();
        Session session = prepareSession();
        Query query = session.createQuery("SELECT distinct xref from org.reactome.psi.Interactor");
        ScrollableResults scrollable = query.scroll();
        int totalEntrezProtein = 0;
        int uniprot = 0;
        int missedUniProt = 0;
        boolean isFound = false;
        int c = 0;
        while (scrollable.next()) {
            Xref xref = (Xref) scrollable.get(0);
            DbReference primaryRef = xref.getPrimaryRef();
            if (primaryRef.getDb().equals("Entrez Protein")) {
                totalEntrezProtein ++;
                isFound = false;
                List<DbReference> secondaryRefList = xref.getSecondaryRefList();
                //if (secondaryRefList == null || secondaryRefList.size() == 0)
                //    continue;
                for (DbReference ref : secondaryRefList) {
                    if (ref.getDb().equals("uniprotkb")) {
                        uniprot ++;
                        isFound = true;
                        break;
                    }
                }
                if (!isFound) {
                    if (c++ < 100)
                        System.out.println(primaryRef.getId());
                    Set<String> set = gi2UniMap.get(primaryRef.getId());
                    if (set != null && set.size() > 0)
                        missedUniProt ++;
                }
            }
        }
        System.out.println("Total Entrez Protein: " + totalEntrezProtein);
        System.out.println("Total UniProt: " + uniprot);
        System.out.println("Missed UniProt Mappings: " + missedUniProt);
    }
    
    public void testViewXrefInInteractor() throws Exception {
        Session session = prepareSession();
        Query query = session.createQuery("SELECT distinct xref from org.reactome.psi.Interactor");
        ScrollableResults scrollable = query.scroll();
        int c = 0;
        Set<String> primaryDBs = new HashSet<String>();
        Set<String> secondaryDBs = new HashSet<String>();
        while (scrollable.next()) {
            Xref xref = (Xref) scrollable.get(0);
            DbReference primaryRef = xref.getPrimaryRef();
//            if (primaryRef != null && primaryRef.getDb() != null) {
//                primaryDBs.add(primaryRef.getDb());
//            }
//            List<DbReference> secondaryRefList = xref.getSecondaryRefList();
//            if (secondaryRefList != null) {
//                for (DbReference ref : secondaryRefList)
//                    secondaryDBs.add(ref.getDb());
//            }
            if (primaryRef != null && primaryRef.getDb().equals("uniprotkb")) 
                continue;
//            DbReference secondaryRef = xref.getSecondaryRefList();
//            if (secondaryRef != null && secondaryRef.getDb().equals("uniprotkb"))
//                continue;
//            if (primaryRef == null && secondaryRef == null)
//                continue;
            if (primaryRef.getDb().equals("Entrez Protein"))// &&
                //secondaryRef == null)
                c ++;
        }
        System.out.println("Total number of Xref needing fix: " + c);
        //System.out.println("Primary DB: " + primaryDBs);
        //System.out.println("Secondary DB: " + secondaryDBs);
    }
}
