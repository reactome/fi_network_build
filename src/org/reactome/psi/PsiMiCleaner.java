/*
 * Created on Apr 18, 2006
 *
 */
package org.reactome.psi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hibernate.Query;
import org.hibernate.ScrollableResults;
import org.hibernate.Session;
import org.hibernate.SessionFactory;



/**
 * This class is used to do some cleaning after Psi-Mi data is loaded: make these instances
 * non-reduant: DbReference, Xref, Bibref, OpenCV, BioSource, Source.
 * @author guanming
 *
 */
public class PsiMiCleaner {
    private SessionFactory sessionFactory = null;
    
    public PsiMiCleaner() {
        
    }
    
    public void setSessionFactory(SessionFactory factory) {
        this.sessionFactory = factory;
    }
    
    /**
     * This method is used to add UniProt DbReference to Xref that are used by
     * Interactor. Interactors from BIND don't have UniProt IDs. A mapping file
     * is needed and passed from argument.
     * @param mapFileName
     * @throws Exception
     */
    public void attachUniProtId(String mapFileName) throws Exception {
        Map<String, String> gi2uid = new HashMap<String, String>();
        Map<String, String> uid2uacc = new HashMap<String, String>();
        loadGi2UniProtMap(mapFileName, gi2uid, uid2uacc);
        // Prepare session
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        // Load all DbReference for UniProtKB to increase peformace
        List dbReferences = session.createQuery("FROM org.reactome.psi.DbReference ref WHERE ref.db=\'uniprotkb\'").list();
        System.out.println("DbReferneces: " + dbReferences.size());
        // Create a cache 
        Map<String, DbReference> dbReferenceMap = new HashMap<String, DbReference>();
        for (Iterator it = dbReferences.iterator(); it.hasNext();) {
            DbReference dbReference = (DbReference) it.next();
            dbReferenceMap.put(dbReference.getSecondary(), dbReference);
        }
        Query query = session.createQuery("SELECT distinct xref from org.reactome.psi.Interactor");
        ScrollableResults scrollable = query.scroll();
        int c = 0;
        boolean isFound = false;
        while (scrollable.next()) {
            Xref xref = (Xref) scrollable.get(0);
            DbReference primaryRef = xref.getPrimaryRef();
            if (primaryRef != null && primaryRef.getDb().equals("uniprotkb")) 
                continue;
            List<DbReference> secondaryRefList = xref.getSecondaryRefList();
            if (secondaryRefList != null) {
                isFound = false;
                for (DbReference ref : secondaryRefList)
                    if (ref.getDb().equals("uniprotkb")) {
                        isFound = true;
                        break;
                    }
                if (isFound)
                    continue;
            }
            if (primaryRef == null && secondaryRefList == null)
                continue;
            if (primaryRef.getDb().startsWith("Entrez") &&
                secondaryRefList == null) {
                if(attachUniProtId(xref, gi2uid, uid2uacc, dbReferenceMap, session))
                    c++;
            }
        }
        System.out.println("Total number of Xref Fix: " + c);
        session.getTransaction().commit();
    }
    
    private boolean attachUniProtId(Xref xref, 
                                    Map<String, String> gi2uid,
                                    Map<String, String> uid2uacc, 
                                    Map<String, DbReference> dbReferenceMap,
                                    Session session) throws Exception {
        DbReference primaryRef = xref.getPrimaryRef();
        String gid = primaryRef.getId();
        String uid = gi2uid.get(gid);
        if (uid == null)
            return false; // Cannot find a match for gid
        // Try to find if a DbReference has already existed
        DbReference uniXref = dbReferenceMap.get(uid);
        if (uniXref == null) {
            uniXref = createUniXref();
            uniXref.setSecondary(uid);
            uniXref.setId(uid2uacc.get(uid));
            session.persist(uniXref);
            dbReferenceMap.put(uid, uniXref);
        }
        xref.addSecondaryRef(uniXref);
        session.update(xref);
        return true;
    }
    
    private DbReference createUniXref() {
        //uniprotkb | MI:0486 | Q96HA0 | q96ha0_human | TrEMBL_26 | identity | MI:0356   
        DbReference dbref = new DbReference();
        dbref.setDb("uniprotkb");
        dbref.setDbAc("MI:0486");
        dbref.setRefType("identity");
        dbref.setRefTypeAc("MI:0356");
        return dbref;
    }
    
    private void loadGi2UniProtMap(String mapFileName, 
                                   Map<String, String> gi2uid, 
                                   Map<String, String> uid2uacc) throws IOException {
        FileReader fileReader = new FileReader(mapFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = reader.readLine(); // Escape the first line
        while ((line = reader.readLine()) != null) {
            String[] tokens = line.split("\t"); 
            // Lower case is used in intact database.
            String uid = tokens[1].toLowerCase();
            String uaccs = tokens[2];
            int index = uaccs.indexOf(";");
            String uacc = null;
            if (index < 0)
                uacc = uaccs.trim();
            else
                uacc = uaccs.substring(0, index).trim();
            uid2uacc.put(uid, uacc);
            String gids = tokens[7];
            String[] gidTokens = gids.split("(;|,)");
            for (String gid : gidTokens) {
                gi2uid.put(gid.trim(), uid);
            }
        }
    }
    
    public void cleanUpSource() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("source.txt");
        cleanTable("Source", getStatement(), idMap);
    }
    
    public void updateSourceReferrers() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("source.txt");
        updateTable("Entry", "source", getStatement(), idMap);
    }
    
    public void generateSourceMergeList() throws Exception {
        Map<String, List<Long>> mergedMap = new HashMap<String, List<Long>>();
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from org.reactome.psi.Source");
        List list = query.list();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Source source = (Source) it.next();
            Names names = source.getNames();
            String nameKey = null;
            if (names != null) {
                nameKey = names.getShortLabel() + "|source|" + names.getFullName();
            }
            Bibref bibref = source.getBibref();
            String bibrefKey = null;
            if (bibref != null)
                bibrefKey = bibref.getDbId() + "";
            Xref xref = source.getXref();
            String xrefKey = null;
            if (xref != null)
                xrefKey = xref.getDbId() + "";
            String key = source.getRelease() + "|source|" + 
                         source.getReleaseDate() + "|source|" + 
                         nameKey + "|source|" + 
                         bibrefKey + "|source|" + 
                         xrefKey;
            List<Long> ids = mergedMap.get(key);
            if (ids == null) {
                ids = new ArrayList<Long>();
                mergedMap.put(key, ids);
            }
            ids.add(source.getDbId());
        }
        Map<Long, List<Long>> idMap = generateIdMap(mergedMap);
        outputMap(idMap, "source.txt");
    }
    
    public void cleanUpBioSource() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("biosource.txt");
        cleanTable("BioSource", getStatement(), idMap);
    }
    
    public void updateBioSourceReferrers() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("biosource.txt");
        Statement stat = getStatement();
        // Create a list of table.columns for merging
        String[][] tableColNames = new String[][] {
                {"ExperimentHostOrganism", "biosourceDbId"}, // Table name followed by column names
                {"Interactor", "organism"},
                {"ParticipantHostOrganism", "bioSourceId"}
        };
        for (String[] tableCol : tableColNames) {
            String tableName = tableCol[0];
            for (int i = 1; i < tableCol.length; i++) {
                String colName = tableCol[i];
                updateTable(tableName, colName, stat, idMap);
            }
            System.out.println("Done: " + tableName);
        }
    }
    
    public void generateBioSourceMergeList() throws Exception {
        Map<String, List<Long>> mergedMap = new HashMap<String, List<Long>>();
        Statement stat = getStatement();
        String query = "SELECT dbId, ncbiTaxId, cellType, compartment, tissue FROM BioSource";
        ResultSet resultSet = stat.executeQuery(query);
        while (resultSet.next()) {
            Long dbId = resultSet.getLong(1);
            Integer taxId = resultSet.getInt(2);
            Long cellType = resultSet.getLong(3);
            Long compartment = resultSet.getLong(4);
            Long tissue = resultSet.getLong(5);
            String key = taxId + "|biosource|" + 
                         cellType + "|biosource|" + 
                         compartment + "|biosource|" + 
                         tissue;
            List<Long> ids = mergedMap.get(key);
            if (ids == null) {
                ids = new ArrayList<Long>();
                mergedMap.put(key, ids);
            }
            ids.add(dbId);
        }
        stat.close();
        Map<Long, List<Long>> idMap = generateIdMap(mergedMap);
        outputMap(idMap, "biosource.txt");
    }
    
    public void cleanUpOpenCV() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("opencv.txt");
        cleanTable("opencv", getStatement(), idMap);
    }
    
    public void updateOpenCVReferrers() throws Exception {
        Map<Long, List<Long>> idMap = inputMap("opencv.txt");
        Statement stat = getStatement();
        // Create a list of table.columns for merging
        String[][] tableColNames = new String[][] {
                {"BaseLocation", "endStatus", "startStatus"}, // Table name followed by column names
                {"BioSource", "cellType", "compartment", "tissue"},
                {"Confidence", "unit"},
                {"Experiment", "featureDetectionMethod", "interactionDetectionMethod", "participantIdentificationMethod"},
                {"Feature", "featureDetectionMethod", "featureType"},
                {"Interaction", "interactionType"},
                {"Interactor", "interactorType"},
                {"OpenCVExperimentalWrapper", "openCV"},
                {"Participant", "biologicalRole"}
        };
        for (String[] tableCol : tableColNames) {
            String tableName = tableCol[0];
            for (int i = 1; i < tableCol.length; i++) {
                String colName = tableCol[i];
                updateTable(tableName, colName, stat, idMap);
            }
            System.out.println("Done: " + tableName);
        }
    }
    
    public void generateOpenCVMergeList() throws Exception {
        String outputFileName = "opencv.txt";
        Map<String, List<Long>> mergedMap = new HashMap<String, List<Long>>();
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from org.reactome.psi.OpenCV");
        int first = 0;
        int max = 10000;
        query.setFirstResult(first);
        query.setMaxResults(max);
        List list = query.list();
        int cycle = 0;
        while (list.size() > 0) {
            cycle ++;
            for (Iterator it = list.iterator(); it.hasNext();) {
                OpenCV cv = (OpenCV) it.next();
                Names names = cv.getNames();
                Xref xref = cv.getXref();
                String xrefKey = null;
                if (xref != null)
                    xrefKey = xref.getDbId() + "";
                String key = names.getShortLabel() + "!opencv!" + 
                             names.getFullName() + "!opencv!" + 
                             xrefKey;
                List<Long> ids = mergedMap.get(key);
                if (ids == null) {
                    ids = new ArrayList<Long>();
                    mergedMap.put(key, ids);
                }
                ids.add(cv.getDbId());
            }
            first += list.size();
            query.setFirstResult(first);
            list = query.list();
            session.clear();
            System.out.printf("Cycle %d: %d%n", cycle, first);
        }
        Map<Long, List<Long>> idMap = generateIdMap(mergedMap);
        outputMap(idMap, outputFileName);
    }
    
    public void generateBibrefMergeList(String outputFileName) throws Exception {
        Map<String, List<Long>> mergedMap = new HashMap<String, List<Long>>();
        Statement stat = getStatement();
        ResultSet resultSet = stat.executeQuery("SELECT dbId, xref FROM Bibref");
        while (resultSet.next()) {
            Long dbId = resultSet.getLong(1);
            String xref = resultSet.getString(2);
            List<Long> ids = mergedMap.get(xref);
            if (ids == null) {
                ids = new ArrayList<Long>();
                mergedMap.put(xref, ids);
            }
            ids.add(dbId);
        }
        stat.close();
        Map<Long, List<Long>> idMap = generateIdMap(mergedMap);
        outputMap(idMap, outputFileName);
    }
    
    public void generateXrefMergeList(String outputFileName) throws Exception {
        Map<String, List<Long>> mergedMap = new HashMap<String, List<Long>>();
        Statement stat = getStatement();
        ResultSet resultSet = stat.executeQuery("SELECT dbId, primaryRef, secondaryRef " +
                "FROM Xref");
        while (resultSet.next()) {
            Long dbId = resultSet.getLong(1);
            Long primaryRef = resultSet.getLong(2);
            Long secondardRef = resultSet.getLong(3);
            String key = primaryRef + "-" + secondardRef;
            List<Long> ids = mergedMap.get(key);
            if (ids == null) {
                ids = new ArrayList<Long>();
                mergedMap.put(key, ids);
            }
            ids.add(dbId);
        }
        stat.close();
        Map<Long, List<Long>> idMap = generateIdMap(mergedMap);
        outputMap(idMap, outputFileName);
    }
    
    private Map<Long, List<Long>> generateIdMap(Map<String, List<Long>> map) {
        Map<Long, List<Long>> idMap = new HashMap<Long, List<Long>>();
        Set<String> keys = map.keySet();
        for (String key : keys) {
            List<Long> ids = map.get(key);
            if (ids.size() < 2)
                continue;
            idMap.put(ids.remove(0), ids);
        }
        return idMap;
    }
    
    public void generateDbReferenceMergeList(String outputFileName) throws Exception {
        Map<DbReference, List<Long>> mergedMap = new HashMap<DbReference, List<Long>>();
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Query query = session.createQuery("from org.reactome.psi.DbReference");
        int first = 0;
        int max = 20000;
        query.setFirstResult(first);
        query.setMaxResults(max);
        List list = query.list();
        int cycle = 0;
        while (list.size() > 0) {
            cycle ++;
            for (Iterator it = list.iterator(); it.hasNext();) {
                DbReference reference = (DbReference) it.next();
                if (mergedMap.containsKey(reference)) {
                    List<Long> dbIds = mergedMap.get(reference);
                    dbIds.add(reference.getDbId());
                }
                else {
                    List<Long> dbIds = new ArrayList<Long>();
                    mergedMap.put(reference, dbIds);
                }
            }
            first += list.size();
            query.setFirstResult(first);
            list = query.list();
            session.clear();
            System.out.printf("Cycle %d: %d%n", cycle, first);
        }
        session.getTransaction().commit();
        Map<Long, List<Long>> idMap = new HashMap<Long, List<Long>>();
        Set<DbReference> references = mergedMap.keySet();
        for (DbReference ref : references) {
            List<Long> ids = mergedMap.get(ref);
            if (ids.size() == 0)
                continue;
            idMap.put(ref.getDbId(), mergedMap.get(ref));
        }
        outputMap(idMap, outputFileName);
    }
    
    private void outputMap(Map<Long, List<Long>> map, String fileName) throws IOException {
        StringBuilder builder = new StringBuilder();
        Set<Long> keys = map.keySet();
        for (Long key : keys) {
            List<Long> ids = map.get(key);
            builder.append(key).append("\t");
            for (Long id : ids) {
                builder.append(id).append(",");
            }
            builder.append("\n");
        }
        FileWriter fileWriter = new FileWriter(fileName);
        BufferedWriter writer = new BufferedWriter(fileWriter);
        writer.write(builder.toString());
        writer.close();
        fileWriter.close();
    }
    
    private Map<Long, List<Long>> inputMap(String mapFileName) throws IOException {
        Map<Long, List<Long>> idMap = new HashMap<Long, List<Long>>();
        FileReader fileReader = new FileReader(mapFileName);
        BufferedReader reader = new BufferedReader(fileReader);
        String line = null;
        int index = 0;
        while ((line = reader.readLine()) != null) {
            index = line.indexOf("\t");
            String key = line.substring(0, index);
            String ids = line.substring(index + 1);
            String[] ids1 = ids.split(",");
            List<Long> idList = new ArrayList<Long>(ids1.length);
            for (String tmp : ids1) {
                idList.add(new Long(tmp));
            }
            idMap.put(new Long(key), idList);
        }
        return idMap;
    }
    
    /**
     * Remove redundent DbReference objects.
     * @param dbRefFileName
     * @throws Exception
     */
    public void cleanUpDbReference(String dbRefFileName) throws Exception {
        Map<Long, List<Long>> idMap = inputMap(dbRefFileName);
        cleanTable("DbReference", getStatement(), idMap);
    }
    
    public void cleanUpXref(String xrefFileName) throws Exception {
        Map<Long, List<Long>> idMap = inputMap(xrefFileName);
        cleanTable("Xref", getStatement(), idMap);
    }
    
    public void cleanUpBibref(String bibrefFileName) throws Exception {
        Map<Long, List<Long>> idMap = inputMap(bibrefFileName);
        cleanTable("Bibref", getStatement(), idMap);
    }
    
    private void cleanTable(String tableName, Statement stat, Map<Long, List<Long>> idMap) throws Exception {
        StringBuilder queryBuilder = new StringBuilder();
        Set<Long> keys = idMap.keySet();
        try {
            for (Long key : keys) {
                List<Long> ids = idMap.get(key);
                queryBuilder.setLength(0);
                queryBuilder.append("DELETE FROM ");
                queryBuilder.append(tableName);
                queryBuilder.append(" WHERE dbId IN (");
                for (Iterator<Long> it = ids.iterator(); it.hasNext();) {
                    queryBuilder.append(it.next());
                    if (it.hasNext())
                        queryBuilder.append(",");
                }
                queryBuilder.append(")");
                stat.executeUpdate(queryBuilder.toString());
            }
        }
        catch(SQLException e) {
            System.out.println("Error in: " + queryBuilder.toString());
            throw e;
        }
    }
    
    public void updateXref(String dbRefFileName) throws Exception {
        Map<Long, List<Long>> idMap = inputMap(dbRefFileName);
        Statement stat = getStatement();
        // First for primaryRef
        updateTable("Xref", "primaryRef", stat, idMap);
        updateTable("Xref", "secondaryRef", stat, idMap);
    }
    
    public void updateXrefReferrers(String xrefFileName) throws Exception {
        String[] tableNames = new String[]{
                "Source",
                "Participant",
                "OpenCV",
                "Interactor",
                "Interaction",
                "Feature",
                "Experiment",
                "Bibref"
        };
        Map<Long, List<Long>> idMap = inputMap(xrefFileName);
        Statement stat = getStatement();
        // To create update statemenet
        for (String tableName : tableNames) {
            updateTable(tableName, "xref", stat, idMap);
            System.out.println("Done: " + tableName);
        }
    }
    
    public void updateBibrefReferrers(String bibrefFileName) throws Exception {
        String[] tableNames = new String[]{
                "Source",
                "Experiment"
        };
        Map<Long, List<Long>> idMap = inputMap(bibrefFileName);
        Statement stat = getStatement();
        for (String tableName : tableNames) {
            updateTable(tableName, "bibref", stat, idMap);
            System.out.println("Done: " + tableName);
        }
    }
    
    private void updateTable(String tableName, 
                             String colName, 
                             Statement stat, 
                             Map<Long, List<Long>> idMap) throws Exception {
        Set<Long> keys = idMap.keySet();
        StringBuilder queryBuilder = new StringBuilder();
        for (Long key : keys) {
            List<Long> ids = idMap.get(key);
            queryBuilder.setLength(0);
            queryBuilder.append("UPDATE ").append(tableName);
            queryBuilder.append(" SET ").append(colName).append("=").append(key);
            queryBuilder.append(" WHERE ").append(colName).append(" in (");
            for (Iterator<Long> it = ids.iterator(); it.hasNext();) {
                queryBuilder.append(it.next());
                if (it.hasNext())
                    queryBuilder.append(",");
            }
            queryBuilder.append(")");
            stat.executeUpdate(queryBuilder.toString());
        }
    }
    
    private Statement getStatement() throws Exception {
        Session session = sessionFactory.getCurrentSession();
        session.beginTransaction();
        Connection connection = session.connection();
        connection.setAutoCommit(true);
        return connection.createStatement();
    }
    
}
