/*
 * Created on Nov 28, 2006
 *
 */
package org.reactome.data;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import org.ensembl.compara.driver.ComparaDriver;
import org.ensembl.compara.driver.ComparaDriverFactory;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to query ensembl MySQL database.
 * @author guanming
 *
 */
public class EnsemblAnalyzer {
    private ComparaDriver driver;
    
    protected void setUp() throws Exception {
        driver = ComparaDriverFactory.createComparaDriver("ensembldb.ensembl.org", 
                                                          3306, 
                                                          "ensembl_compara_41", 
                                                          "anonymous", 
                                                           "");
    }
    
    public void simpleTest() throws Exception {
        List<String> ids = new ArrayList<String>();
        ids.add("P26599-2");
        Map<String, Integer> counts = getGeneCopyNumbers(ids);
        for (String id : ids)
            System.out.printf("%s - %d%n", id, counts.get(id));
    }
    
    public int getGeneCopyNumber(String uniProtId) throws Exception {
        List<String> ids = new ArrayList<String>();
        ids.add(uniProtId);
        Map<String, Integer> counts = getGeneCopyNumbers(ids);
        return counts.get(uniProtId);
    }
    
    public Map<String, Integer> getGeneCopyNumbers(List<String> uniProtIds) throws Exception {
        if (driver == null)
            setUp();
        Map<String, Integer> countMap = new HashMap<String, Integer>();
        Connection connection = driver.getConnection();
        String sql = "select f.family_id from family_member f, member m " +
                     "where m.stable_id = ? and m.member_id = f.member_id";
        PreparedStatement stat = connection.prepareStatement(sql);
        // Query for members
        String sql1 = "select count(m.member_id) from member m, family_member f " +
                      "where m.member_id = f.member_id and f.family_id = ? " +
                      "and m.taxon_id = 9606 and m.stable_id like 'ENSG%'";
        PreparedStatement stat1 = connection.prepareStatement(sql1);
        String id = null;
        for (String id1 : uniProtIds) {
            id = preProcessId(id1);
            stat.setString(1, id);
            ResultSet result = stat.executeQuery();
            Integer familyId = null;
            while (result.next())
                familyId = result.getInt(1);
            result.close();
            //System.out.println("Family ID for " + id + ": " + familyId);
            if (familyId == null) {
                countMap.put(id1, 0);
                continue;
            }
            stat1.setInt(1, familyId);
            int count = 0;
            result = stat1.executeQuery();
            while (result.next())
                count = result.getInt(1);
            //System.out.println("Total Gene Number: " + count);
            countMap.put(id1, count);
        }
        stat.close();
        stat1.close();
        // Need to close connection so that it can be reused.
        connection.close();
        return countMap;
    }
    
    private String preProcessId(String id) {
        // Don't use variant: it is not in the ensembl database
        int index = id.indexOf("-");
        if (index < 0)
            return id;
        return id.substring(0, index);
    }
    
    private Connection getConnection() throws Exception {
        Class.forName("com.mysql.jdbc.Driver");
        String url = "jdbc:mysql://localhost:3306/" + FIConfiguration.getConfiguration().get("ENSEMBL_COMPARA_DATABASE");
        Properties info = new Properties();
        info.setProperty("user", FIConfiguration.getConfiguration().get("DB_USER"));
        info.setProperty("password", FIConfiguration.getConfiguration().get("DB_PWD"));
        Connection connection = DriverManager.getConnection(url, info);
        return connection;
    }
    
    /**
     * This method is used to dump protein families for four organisms.
     * @throws Exception
     */
    @Test
    public void dumpProteinFamilies() throws Exception {
        // Pick uniprot identifiers only
        String queryString = "select f.family_id, m.stable_id from family_member f," +
        		" member m where f.member_id = m.member_id and m.taxon_id = ? " +
        		"and m.source_name like 'UniProt%'";
        Connection connection = getConnection();
        PreparedStatement stat = connection.prepareStatement(queryString);
        List<Integer> neededTaxonIds = getNeededTaxonIds();
        Map<String, Set<String>> familyToProteins = new HashMap<String, Set<String>>();
        for (Integer taxonId : neededTaxonIds) {
            stat.setInt(1, taxonId);
            ResultSet resultset = stat.executeQuery();
            while (resultset.next()) {
                String family = resultset.getString(1);
                String uniProtId = resultset.getString(2);
                Set<String> set = familyToProteins.get(family);
                if (set == null) {
                    set = new HashSet<String>();
                    familyToProteins.put(family, set);
                }
                set.add(taxonId + ":" + uniProtId);
            }
            resultset.close();
//            System.out.println("Finish taxon: " + taxonId);
        }
        stat.close();
        connection.close();
        System.out.println("Size of familyToProteins(): " + familyToProteins.size());
        filterFamilies(familyToProteins);
        FileUtility fu = new FileUtility();
        fu.saveSetMapInSort(familyToProteins, 
                            FIConfiguration.getConfiguration().get("ENSEMBL_PROTEIN_FAMILIES"));
    }
    
    public Map<String, Set<String>> loadYeastToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("4932");
    }
    
    public Map<String, Set<String>> loadWormToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("6239");
    }
    
    public Map<String, Set<String>> loadFlyToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("7227");
    }
    
    public Map<String, Set<String>> loadMouseToHumanMapInUniProt() throws IOException {
        return loadToHumanMapInUniProt("10090");
    }
    
    private Map<String, Set<String>> loadToHumanMapInUniProt(String taxonId) throws IOException {
        FileUtility fu = new FileUtility();
        String fileName = FIConfiguration.getConfiguration().get("ENSEMBL_PROTEIN_FAMILIES");
        Map<String, Set<String>> familyToProteins = fu.loadSetMap(fileName);
        Set<String> humanIds = new HashSet<String>();
        Set<String> otherIds = new HashSet<String>();
        // To be returned
        Map<String, Set<String>> map = new HashMap<String, Set<String>>();
        for (String family : familyToProteins.keySet()) {
            Set<String> proteins = familyToProteins.get(family);
            humanIds.clear();
            otherIds.clear();
            splitIds(proteins, humanIds, otherIds, taxonId);
            for (String otherId : otherIds) {
                Set<String> humanSet = map.get(otherId);
                if (humanSet == null) {
                    humanSet = new HashSet<String>();
                    map.put(otherId, humanSet);
                }
                humanSet.addAll(humanIds);
            }
        }
        return map;
    }
    
    private void splitIds(Set<String> proteins,
                          Set<String> humanIds,
                          Set<String> otherIds,
                          String taxonId) {
        for (String protein : proteins) {
            if (protein.startsWith("9606:")) {
                // This is a human protein
                humanIds.add(protein.substring(5));
            }
            else if (protein.startsWith(taxonId)) {
                otherIds.add(protein.substring(taxonId.length() + 1)); // 1 for ":".
            }
        }
    }
     
    
    /**
     * Filter protein families so that a family containing at least two species and one of them
     * should be homo sapiens.
     * @param familyToProteins
     */
    private void filterFamilies(Map<String, Set<String>> familyToProteins) {
        System.out.println("Total families before filtering: " + familyToProteins.size());
        for (Iterator<String> it = familyToProteins.keySet().iterator(); it.hasNext();) {
            String family = it.next();
            Set<String> proteins = familyToProteins.get(family);
            Set<String> species = extractSpecies(proteins);
            if (species.size() == 1 ||
                !species.contains("9606")) {
                it.remove();
            }
        }
        System.out.println("Total families after filtering: " + familyToProteins.size());
    }
    
    private Set<String> extractSpecies(Set<String> proteins) {
        Set<String> species = new HashSet<String>();
        int index = 0;
        for (String protein : proteins) {
            index = protein.indexOf(":");
            species.add(protein.substring(0, index));
        }
        return species;
    }
    
    /**
     * This method is used to filter a big sequence.txt file downloaded from ensembl
     * compara to sequence files used for four species only.
     * @throws Exception
     */
    public void filterSequences() throws Exception {
        List<Long> sequenceIds = getSequenceIdsForNeededTaxons();
        Set<Long> idSet = new HashSet<Long>(sequenceIds); // Should be fast
        System.out.println("Total ids: " + idSet.size());
        String inFile = FIConfiguration.getConfiguration().get("ENSEMBL_DIR") + "sequence.txt";
        String outFile = FIConfiguration.getConfiguration().get("ENSEMBL_DIR") + "sequence_filtered.txt";
        FileUtility inFu = new FileUtility();
        inFu.setInput(inFile);
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outFile);
        String line = null;
        int index = 0;
        String id = null;
        while ((line = inFu.readLine()) != null) {
            index = line.indexOf("\t");
            id = line.substring(0, index);
            if(idSet.contains(Long.parseLong(id)))
                outFu.printLine(line);
        }
        inFu.close();
        outFu.close();
    }
    
    private List<Long> getSequenceIdsForNeededTaxons() throws Exception {
        Connection connection = getConnection();
        List<Integer> taxonIds = getNeededTaxonIds();
        List<Long> sequenceIds = new ArrayList<Long>();
        String query = "SELECT sequence_id FROM member WHERE taxon_id = ?";
        PreparedStatement stat = connection.prepareStatement(query);
        for (Integer id : taxonIds) {
            stat.setInt(1, id);
            ResultSet resultSet = stat.executeQuery();
            while (resultSet.next()) {
                sequenceIds.add(resultSet.getLong(1));
            }
            resultSet.close();
        }
        stat.close();
        connection.close();
        return sequenceIds;
    }
    
    private List<Integer> getNeededTaxonIds() {
        List<Integer> rtns = new ArrayList<Integer>();
        rtns.add(9606); // homo sapiens
        rtns.add(4932); // S. cerevisiae
        rtns.add(6239); // C. elegans
        rtns.add(7227); // D. melanogaster
        rtns.add(10090); // Mus musculus
        return rtns;
    }
    
}
