/*
 * Created on Apr 8, 2008
 *
 */
package org.reactome.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.junit.Test;
import org.reactome.fi.util.ChecksumCalculator;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.funcInt.Protein;

/**
 * This class is used to handle protein sequence
 * @author wgm
 *
 */
public class ProteinSequenceHandler {
    private final String UNIPROT_DIR = FIConfiguration.getConfiguration().get("UNIPROT_DIR");
    private final String BIND_DIR = "../../caBIG_R3/datasets/BIND/";
    private final String NCBI_PROTEIN_URL = "http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?db=protein&retmode=text&id=";
    private final String HPRD_DIR = "../../caBIG_R3/datasets/HPRD/";
    // Used to cache some maps
    private Map<String, Sequence> idToSeq;
    private Map<String, String> checksumToDbAcc;
    // To calcualte checksum
    private ChecksumCalculator checksum;
    
    public ProteinSequenceHandler() {
    }
    
    /**
     * Use this method to consolidate a list of protein ids based on sequence checksums.
     * Some of protein ids may point to the same sequences.
     * @param ids
     * @return
     * @throws Exception
     */
    public Set<String> consolidateProteinIds(Set<String> ids) throws Exception {
        System.out.println("Before sequence consolidating: " + ids.size());
        if (idToSeq == null) {
            idToSeq = getAllSequences();
        }
        if (checksumToDbAcc == null) {
            checksumToDbAcc = generateUniqueChecksumToAccessionMap();
        }
        Set<String> idsInChecksum = new HashSet<String>();
        Sequence seq;
        int compare = 0;
        for (String id : ids) {
            // Need to test
            seq = getSequence(id, idToSeq);
            if (seq == null) {// Cannot do anything
                idsInChecksum.add(id);
                continue;
            }
            String checksum = seq.getChecksum();
            String dbAcc = checksumToDbAcc.get(checksum);
            id = getAccessionFromDbAcc(dbAcc);
            idsInChecksum.add(id);
        }
        System.out.println("After sequence consolidating: " + idsInChecksum.size());
        return idsInChecksum;
    }
    
    /**
     * Use this method to consolidate a list of protein pairs (FIs) based on sequence checksums.
     * Some of ids may point to the same proteins based on sequences.
     * @param interactions
     * @return
     * @throws Exception
     */
    public Set<String> consolidateInteractionsUseChecksum(Set<String> interactions) throws Exception {
        System.out.println("Before sequence consolidating: " + interactions.size());
        if (idToSeq == null) {
            idToSeq = getAllSequences();
        }
        if (checksumToDbAcc == null) {
            checksumToDbAcc = generateUniqueChecksumToAccessionMap();
        }
        Set<String> intInChecksum = new HashSet<String>();
        int index = 0;
        String id1, id2;
        Sequence seq1, seq2;
        int compare = 0;
        for (String pair : interactions) {
            index = pair.indexOf("\t");
            id1 = pair.substring(0, index);
            id2 = pair.substring(index + 1);
            // Need to test
            seq1 = getSequence(id1, idToSeq);
            if (seq1 == null)
                continue;
            seq2 = getSequence(id2, idToSeq);
            if (seq2 == null)
                continue;
            String checksum1 = seq1.getChecksum();
            String dbAcc1 = checksumToDbAcc.get(checksum1);
            id1 = getAccessionFromDbAcc(dbAcc1);
            String checksum2 = seq2.getChecksum();
            String dbAcc2 = checksumToDbAcc.get(checksum2);
            id2 = getAccessionFromDbAcc(dbAcc2);
            compare = id1.compareTo(id2);
            if (compare < 0)
                intInChecksum.add(id1 + "\t" + id2);
            else if (compare > 0)
                intInChecksum.add(id2 + "\t" + id1);
        }
        System.out.println("After sequence consolidating: " + intInChecksum.size());
        return intInChecksum;
    }
    
    private String getAccessionFromDbAcc(String dbAcc) {
        int index = dbAcc.indexOf(":");
        return dbAcc.substring(index + 1);
    }
    
    /**
     * The passed id should have db name information. e.g.: UniProt:P00533.
     * @param id
     * @return
     * @throws Exception
     */
    public Sequence getSequence(String dbAccession) throws Exception {
        if (idToSeq == null)
            idToSeq = getAllSequences();
        Sequence seq = idToSeq.get(dbAccession);
        if (seq == null) {
            // in case a splice form. Make another try.
            if (dbAccession.startsWith("UniProt") &&
                dbAccession.contains("-")) {
                // fall back to the parent sequence
                int index = dbAccession.indexOf("-");
                String tmp = dbAccession.substring(0, index);
                seq = idToSeq.get(tmp);
                if (seq != null)
                    System.out.println("Original Sequence is used for: " + dbAccession);

            }
        }
        return seq;
    }
    
    /**
     * From a sequence checksum to get a unique protein database accession.
     * @param checksum
     * @return
     * @throws Exception
     */
    public String getDbAccFromChecksum(String checksum) throws Exception {
        if (checksumToDbAcc == null)
            checksumToDbAcc = generateUniqueChecksumToAccessionMap();
        return checksumToDbAcc.get(checksum);
    }
    
    private Sequence getSequence(String id, 
                                 Map<String, Sequence> idToSequence) {
        Sequence sequence = idToSequence.get("UniProt:" + id);
        if (sequence != null)
            return sequence;
        if (id.contains("-")) {
            // UniProt splice form
            int index = id.indexOf("-");
            id = id.substring(0, index);
            sequence = idToSequence.get("UniProt:" + id);
            //if (sequence == null)
            //    System.out.println("Cannot get sequence: " + id);
            return sequence;
        }
        sequence = idToSequence.get("HPRD:" + id);
        if (sequence != null)
            return sequence;
        sequence = idToSequence.get("Entrez Protein:" + id);
        //if (sequence == null)
        //    System.out.println("Cannot get sequence: " + id);
        return sequence;
    }
    
    /**
     * Use this method to get a map from a sequence checksum to a unique database
     * accession number. The database accession is picked up in the order of this:
     * SwissProt, Trembl, HPRD, and Entrez Protein. The client should call this method
     * to consolidate database accessions based on amino acid sequences.
     * @return
     * @throws Exception
     */
    public Map<String, String> generateUniqueChecksumToAccessionMap() throws Exception {
        Map<String, Sequence> swissProtSeqs = loadSwissProtSequences();
        Map<String, Sequence> swissProtVarSeqs = loadSwissProtVarSequences();
        Map<String, Sequence> tremblSeqs = loadTremblSequences();
//        Map<String, Sequence> hprdSeqs = loadHprdSequences();
//        Map<String, Sequence> entrezSeqs = loadEntrezProteinSequences();
        List<Map<String, Sequence>> sequenceMaps = new ArrayList<Map<String, Sequence>>();
        // The following adding sequence is important: to make sure 
        // accession numbers from SwissProt is used first, then Trembl.
        sequenceMaps.add(swissProtSeqs);
        sequenceMaps.add(swissProtVarSeqs);
        sequenceMaps.add(tremblSeqs);
//        sequenceMaps.add(hprdSeqs);
//        sequenceMaps.add(entrezSeqs);
        Map<String, String> checksumTodbAcc = new HashMap<String, String>();
        for (Map<String, Sequence> map : sequenceMaps) {
            for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
                String dbAcc = it.next();
                Sequence sequence = map.get(dbAcc);
                String checksum = sequence.getChecksum();
                if (checksumTodbAcc.containsKey(checksum))
                    continue;
                checksumTodbAcc.put(checksum,
                                    dbAcc);
            }
        }
        return checksumTodbAcc;
    }
    
    /**
     * This method is used to check some properties for all sequences.
     * @throws Exception
     */
    @Test
    public void checkSequences() throws Exception {
        Map<String, Sequence> swissProtSeqs = loadSwissProtSequences();
        Map<String, Sequence> swissProtVarSeqs = loadSwissProtVarSequences();
        Map<String, Sequence> tremblSeqs = loadTremblSequences();
        Map<String, Sequence> entrezSeqs = loadEntrezProteinSequences();
        Map<String, Sequence> hprdSeqs = loadHprdSequences();
        // Check how many sequences we have
        Set<String> totalSequences = new HashSet<String>();
        Set<String> totalChecksums = new HashSet<String>();
        List<Map<String, Sequence>> sequenceList = new ArrayList<Map<String, Sequence>>();
        sequenceList.add(swissProtSeqs);
        sequenceList.add(swissProtVarSeqs);
        sequenceList.add(tremblSeqs);
        sequenceList.add(entrezSeqs);
        sequenceList.add(hprdSeqs);
        for (Map<String, Sequence> map : sequenceList) {
            for (Iterator<String> it = map.keySet().iterator(); it.hasNext();) {
                String id = it.next();
                Sequence sequence = map.get(id);
                totalSequences.add(sequence.getSequence());
                totalChecksums.add(sequence.getChecksum());
            }
        }
        System.out.println("Total sequences: " + totalSequences.size());
        System.out.println("Total checksum: " + totalChecksums.size());
    }
    
    /**
     * This method is used to return all sequences from UniProt (including splice forms),
     * HPRD and some Entrez proteins (only used by the FI database).
     * @return
     * @throws Exception
     */
    public Map<String, Sequence> getAllSequences() throws Exception {
        if (idToSeq == null) {
            Map<String, Sequence> allSequences = new HashMap<String, Sequence>();
            Map<String, Sequence> swissProtSeqs = loadSwissProtSequences();
            Map<String, Sequence> swissProtVarSeqs = loadSwissProtVarSequences();
            Map<String, Sequence> tremblSeqs = loadTremblSequences();
//            Map<String, Sequence> entrezSeqs = loadEntrezProteinSequences();
//            Map<String, Sequence> hprdSeqs = loadHprdSequences();
            allSequences.putAll(swissProtSeqs);
            allSequences.putAll(swissProtVarSeqs);
            allSequences.putAll(tremblSeqs);
//            allSequences.putAll(entrezSeqs);
//            allSequences.putAll(hprdSeqs);
            idToSeq = allSequences;
        }
        return idToSeq;
    }
    
    /**
     * Get the map from sequence checksums to uniprot ids.
     * @return
     * @throws Exception
     */
    public Map<String, String> getChecksumToUniprotId() throws IOException {
        Map<String, Sequence> allSequences = new HashMap<String, Sequence>();
        Map<String, Sequence> swissProtSeqs = loadSwissProtSequences();
//        Map<String, Sequence> swissProtVarSeqs = loadSwissProtVarSequences();
        Map<String, Sequence> tremblSeqs = loadTremblSequences();
        allSequences.putAll(swissProtSeqs);
//        allSequences.putAll(swissProtVarSeqs);
        allSequences.putAll(tremblSeqs);
        Map<String, String> checksumToId = new HashMap<String, String>();
        int index = 0;
        for (String id : allSequences.keySet()) {
            Sequence seq = allSequences.get(id);
            // remove database label
            index = id.indexOf(":");
            id = id.substring(index + 1);           
            checksumToId.put(seq.getChecksum(), id);
        }
        return checksumToId;
    }
    
    public Map<String, String> getUniProtIdToChecksum() throws IOException {
        Map<String, Sequence> allSequences = new HashMap<String, Sequence>();
        Map<String, Sequence> swissProtSeqs = loadSwissProtSequences();
        Map<String, Sequence> swissProtVarSeqs = loadSwissProtVarSequences();
        Map<String, Sequence> tremblSeqs = loadTremblSequences();
        allSequences.putAll(swissProtSeqs);
        allSequences.putAll(swissProtVarSeqs);
        allSequences.putAll(tremblSeqs);
        Map<String, String> idToChecksum = new HashMap<String, String>();
        int index = 0;
        for (String id : allSequences.keySet()) {
            Sequence seq = allSequences.get(id);
            // remove database label
            index = id.indexOf(":");
            id = id.substring(index + 1);           
            idToChecksum.put(id, seq.getChecksum());
        }
        return idToChecksum;
    }
    
    /**
     * This method is used to load protein sequences for HPRD proteins.
     * @return
     * @throws Exception
     */
    public Map<String, Sequence> loadHprdSequences() throws Exception {
        String fileName = HPRD_DIR + "HprdProteinSequences.fasta";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        String ac = null;
        StringBuilder seqBuilder = new StringBuilder();
        Map<String, Sequence> acToSeq = new HashMap<String, Sequence>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith(">")) {
                // We are interested in HUMAN proteins only
                ac = line.substring(1);
                if (seqBuilder.length() > 0) {
                    // Need to create a new Sequence object
                    Sequence sequence = new Sequence();
                    String seq = seqBuilder.toString();
                    sequence.setSequence(seq);
                    sequence.setChecksum(calculateChecksum(seq));
                    acToSeq.put(ac, sequence);
                    seqBuilder.setLength(0);
                }
            }
            else {
                // There are only two types of lines: Ac definition and Sequence
                seqBuilder.append(line);
            }
        }
        fu.close();
        return acToSeq;
    }
    
    /**
     * This method is used to generate a sequence file for HPRD ids from a list
     * HPRD XML files that contain information for HPRD entries.
     * @throws Exception
     */
    @Test
    public void generateHPRDSequenceFile() throws Exception {
        String dirName = HPRD_DIR + "HPRD_XML_060106/";
        String outputFileName = HPRD_DIR + "HprdProteinSequences.fasta";
        FileUtility outFu = new FileUtility();
        outFu.setOutput(outputFileName);
        File dir = new File(dirName);
        File[] files = dir.listFiles();
        FileUtility fu = new FileUtility();
        String fileName = null;
        int index = 0;
        String id = null;
        String line = null;
        boolean inSequence = false;
        for (File file : files) {
            fileName = file.getName();
            if (!fileName.endsWith(".xml"))
                continue;
            // Extrac HPRD id from file name
            index = fileName.indexOf(".");
            id = fileName.substring(0, index);
            // Start read the file and extract aa sequence from file
            fu.setInput(file.getAbsolutePath());
            while ((line = fu.readLine()) != null) {
                line = line.trim();
                if (line.startsWith("<protein_sequence")) {
                    String aaSeq = getSequenceFromHPRDFile(fu);
                    // Need to split aaSeq into 60 AA in one line as a standard fast file does
                    outputSequenceInFasta(aaSeq,
                                     "HPRD:" + id,
                                     outFu);
                    break;
                }
            }
            fu.close();
        }
        outFu.close();
    }
    
    private void outputSequenceInFasta(String seq,
                                       String id,
                                       FileUtility outFu) throws IOException {
        outFu.printLine(">" + id);
        // Need to split sequence as sixty AAs per line
        char[] aas = seq.toCharArray();
        int length = 60;
        int cycle = aas.length / 60; 
        if (cycle * length < aas.length)
            cycle ++;
        int currentLength = length;
        for (int i = 0; i < cycle; i++) {
            if (i * length + length > aas.length)
                currentLength = aas.length - i * length;
            outFu.printLine(new String(aas, i * length, currentLength));
        }
    }
    
    private String getSequenceFromHPRDFile(FileUtility fu) throws IOException {
        StringBuilder seqBuilder = new StringBuilder();
        String line = null;
        while ((line = fu.readLine()) != null) {
            line = line.trim();
            if (line.equals("</protein_sequence>"))
                break;
            seqBuilder.append(line);
        }
        String seq = seqBuilder.toString();
        seq = seq.replace(" ", "");
        return seq.toUpperCase();
    }
    
    /**
     * This method is used to check the total number of proteins in our FI database
     * based on sequence.
     * @throws Exception
     */
    @Test
    public void analyzeProteinsBasedOnSequences() throws Exception {
        // Cannot use hibernate API to query a database different from
        // set in the configure file though we can still create another
        // configure file. However, then we need to use another mapping
        // file since we have changed some protein properties
        Class.forName("com.mysql.jdbc.Driver");
        String url = "jdbc:mysql://localhost:3306/FunctionalInteractions_v3";
        Connection connection = DriverManager.getConnection(url,
                                                            "root",
                                                            "macmysql01");
        Statement statement = connection.createStatement();
        String query = "SELECT dbName, accession FROM protein";
        List<Protein> proteins = new ArrayList<Protein>();
        ResultSet resultSet = statement.executeQuery(query);
        while (resultSet.next()) {
            String dbName = resultSet.getString(1);
            String accession = resultSet.getString(2);
            Protein protein = new Protein();
            protein.setPrimaryAccession(accession);
            protein.setPrimaryDbName(dbName);
            proteins.add(protein);
        }
        resultSet.close();
        statement.close();
        connection.close();
        System.out.println("Total proteins: " + proteins.size());
        // We used four sources for protein identifiers: SwissProt, Trembl, HPRD
        // and Entrez Protein. 
        Map<String, Sequence> idToSequence = getAllSequences();
        Set<String> usedSequences = new HashSet<String>();
        Set<String> usedChecksums = new HashSet<String>();
        int totalUnknown = 0;
        int index = 0;
        for (Protein protein : proteins) {
            String id = protein.getPrimaryAccession();
            // The first isoform from UniProt should use the main sequence
            if (id.endsWith("-1")) {
                index = id.lastIndexOf("-");
                id = id.substring(0, index);
            }
            Sequence sequence = idToSequence.get(protein.getPrimaryDbName() + ":" +
                                                 protein.getPrimaryAccession());
            if (sequence == null) {
                System.out.println(id + ": no sequence");
                totalUnknown ++;
                continue;
            }
            usedSequences.add(sequence.getSequence());
            usedChecksums.add(sequence.getChecksum());
        }
        System.out.println("Total proteins: " + proteins.size());
        System.out.println("Total used sequence: " + usedSequences.size());
        System.out.println("Total used checksum: " + usedChecksums.size());
        System.out.println("Total unknown: " + totalUnknown);
        // Calculate what ids have been merged
        Map<String, List<String>> checksumToIds = new HashMap<String, List<String>>();
        for (Protein protein : proteins) {
            String id = protein.getPrimaryAccession();
            if (id.endsWith("-1")) {
                index = id.lastIndexOf("-");
                id = id.substring(0, index);
            }
            Sequence sequence = idToSequence.get(protein.getPrimaryDbName() + ":" + id);
            if (sequence == null) {
                System.out.println(id + ": no sequence");
                totalUnknown ++;
                continue;
            }
            List<String> list = checksumToIds.get(sequence.getChecksum());
            if (list == null) {
                list = new ArrayList<String>();
                checksumToIds.put(sequence.getChecksum(),
                                  list);
            }
            list.add(id);
        }
        // Check what id list has more than one
        int c = 0;
        for (Iterator<String> it = checksumToIds.keySet().iterator(); it.hasNext();) {
            String checksum = it.next();
            List<String> idList = checksumToIds.get(checksum);
            if (idList.size() > 1) {
                System.out.println(checksum + ": " + idList);
                c ++;
            }
        }
        System.out.println("Ids should be merged: " + c);
    }
    
    /**
     * This method is used to load sequence information from the file got from
     * Entrez protein database server.
     * @return
     * @throws Exception
     */
    public Map<String, Sequence> loadEntrezProteinSequences() throws Exception {
        // Actual work
        String fileName = BIND_DIR + "EntrezProteinFromNBCI.txt";
        Map<String, Sequence> entrezToSeq = new HashMap<String, Sequence>();
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        String ac = null;
        String source;
        int index = 0;
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("GI:")) {
                ac = line.substring(3); // After GI:
            }
            else if (line.startsWith("SOURCE")) {
                // Check if it is from human
                index = line.indexOf(" ");
                source = line.substring(index + 1).trim();
                if (!source.equals("Homo sapiens (human)")) {
                    ac = null; // Used as a flag
                }
            }
            else if (ac != null && line.startsWith("ORIGIN")) {
                // Need to get sequence information
                String aaSeq = getSequenceFromNCBI(fu);
                String checksumValue = calculateChecksum(aaSeq);
                Sequence sequence = new Sequence();
                sequence.setSequence(aaSeq);
                sequence.setChecksum(checksumValue);
                entrezToSeq.put("Entrez Protein:" + ac, sequence);
            }
        }
        fu.close();
        return entrezToSeq;
    }
    
    private String getSequenceFromNCBI(FileUtility fu) throws IOException {
        StringBuilder builder = new StringBuilder();
        String line = null;
        String line1 = null;
        int index;
        while ((line = fu.readLine()) != null) {
            if (line.equals("//"))
                break;
            line1 = line.trim();
            // 1 msrqtatalp tgtskcppsq rvpaltgtta snndlaslfe cpvcfdyvlp pilqcqsghl
            index = line1.indexOf(" ");
            line1 = line1.substring(index + 1);
            builder.append(line1);
        }
        String aa = builder.toString();
        aa = aa.replaceAll(" ", ""); // Need to remove any empty
        return aa.toUpperCase();
    }
    
    /**
     * This method is used to get the full entries for a list of Entrez Protein
     * used by BIND and cannot be mapped to SwissProt.
     * @throws Exception
     */
    @Test
    public void fetchEntrezProteinsFromNCBI() throws Exception {
        List<String> ids = generateEntrezProteinList();
        // These two ids are wrong:
        ids.remove("1"); // Wrong id: doesn't exist in Entrez Protein
        ids.remove("2295431"); // This is a nucleotide id
        // Want to output these query results
        String outputFileName = "EntrezProteinFromNBCI.txt";
        FileUtility fu = new FileUtility();
        fu.setOutput(BIND_DIR + outputFileName);
        int c = 0;
        for (String id : ids) {
            System.out.println("Protein: " + c++);
            fu.printLine("GI:" + id);
            URL url = new URL(NCBI_PROTEIN_URL + id);
            InputStream is = url.openStream();
            InputStreamReader reader = new InputStreamReader(is);
            BufferedReader bufferedReader = new BufferedReader(reader);
            String line = null;
            while ((line = bufferedReader.readLine()) != null) {
                fu.printLine(line);
            }
            bufferedReader.close();
            reader.close();
            is.close();
        }
        fu.close();
    }
    
    public List<String> generateEntrezProteinList() throws Exception {
        String clsName = "com.mysql.jdbc.Driver";
        Class.forName(clsName).newInstance();
        String url = "jdbc:mysql://localhost:3306/FunctionalInteractions_v3";
        Connection connect = DriverManager.getConnection(url, 
                                                         "root", 
                                                         "macmysql01");
        String query = "SELECT accession FROM Protein WHERE dbName = 'Entrez Protein'";
        Statement stat = connect.createStatement();
        ResultSet result = stat.executeQuery(query);
        List<String> entrezProteins = new ArrayList<String>();
        while (result.next()) {
           // System.out.println(result.getString(1));
            entrezProteins.add(result.getString(1));
        }
        result.close();
        stat.close();
        connect.close();
        return entrezProteins;
    }
    
    /**
     * This method is used to find what ids have been missed in the Entrez Protein database.
     * The query results are got by using Batch Query at Entrez:
     * http://www.ncbi.nlm.nih.gov/sites/batchentrez?db=Nucleotide
     * @throws Exception
     */
    @Test
    public void checkEntrezQueryResults() throws Exception { 
        List<String> entrezProteins = generateEntrezProteinList();
        // Need to load the query results
        String fileName = BIND_DIR + "EntrezProteinSearchResults.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        // Used to check if a line starting with number:
        String ac = null;
        List<String> queriedIds = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("//")) // for comment
                continue;
            if (line.startsWith("gi")) {
                // gi|4506065|ref|NP_002727.1|[4506065]
                int index1 = line.indexOf("|");
                int index2 = line.indexOf("|", index1 + 1);
                ac = line.substring(index1 + 1,
                                    index2);
                queriedIds.add(ac);
                continue;
            }
        }
        fu.close();
        // There are two ids that cannot be queried. Use the following
        // statements to find these two ids
        entrezProteins.removeAll(queriedIds);
        System.out.println("Queried Ids: " + entrezProteins);
    }
    
    /**
     * For some unknown reason, some Entrez proteins cannot generate summary text.
     * This method is used to get a list of Entrez proteins that cannot generate
     * summary text.
     * @throws Exception
     */
    @Test
    public void validateQueryResults() throws Exception {
        String fileName = BIND_DIR + "EntrezProteinSearchResults.txt";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        // Used to check if a line starting with number:
        line = "1: CAA43659";
        String regExp = "^(\\d+):";
        Pattern pattern = Pattern.compile(regExp);
        while ((line = fu.readLine()) != null) {
            Matcher matcher = pattern.matcher(line);
            if (matcher.find()) {
                if (line.contains("cannot get document summary")) {
                    // Eg: 1183:  UID 25952131: cannot get document summary
                    int index1 = line.indexOf("UID");
                    int index2 = line.indexOf(":", index1);
                    String id = line.substring(index1 + 3,
                                               index2).trim();
                    System.out.println(id);
                }
            }
        }
        fu.close();
    }
    
    @Test
    public void testLoadSequences() throws Exception {
        long time1 = System.currentTimeMillis();
        long mem1 = Runtime.getRuntime().totalMemory();
        System.out.println("Before loading: " + mem1/1024);
        //Map<String, Sequence> acToSequence = loadSwissProtVarSequences();
        Map<String, Sequence> acToSequence = loadEntrezProteinSequences();
        System.out.println("Total sequences: " + acToSequence.size());
        // Just check the first five
        int c = 0;
        for (Iterator<String> it = acToSequence.keySet().iterator(); it.hasNext();) {
            String ac = it.next();
            Sequence seq = acToSequence.get(ac);
            System.out.println(ac + " " + seq.getChecksum());
            System.out.println(seq.getSequence());
            System.out.println();
            c ++;
            if (c > 5)
                break;
        }
        long time2 = System.currentTimeMillis();
        long mem2 = Runtime.getRuntime().totalMemory();
        System.out.println("After loading: " + mem2 / 1024);
        System.out.println("Total time: " + (time2 - time1));
    }
    
    public Map<String, Sequence> loadSwissProtVarSequences() throws IOException {
        String fileName = UNIPROT_DIR + "uniprot_sprot_varsplic.fasta";
        FileUtility fu = new FileUtility();
        fu.setInput(fileName);
        String line = null;
        String ac = null;
        StringBuilder seqBuilder = new StringBuilder();
        Map<String, Sequence> acToSeq = new HashMap<String, Sequence>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith(">") && line.contains("_HUMAN")) {
                // We are interested in HUMAN proteins only
                ac = getAcFromFasta(line);
                if (seqBuilder.length() > 0) {
                    // Need to create a new Sequence object
                    Sequence sequence = new Sequence();
                    String seq = seqBuilder.toString();
                    sequence.setSequence(seq);
                    sequence.setChecksum(calculateChecksum(seq));
                    acToSeq.put("UniProt:" + ac, sequence);
                    seqBuilder.setLength(0);
                }
            }
            else {
                // There are only two types of lines: Ac definition and Sequence
                seqBuilder.append(line);
            }
        }
        fu.close();
        return acToSeq;
    }

    private String calculateChecksum(String sequence) {
        if (checksum == null)
            checksum = new ChecksumCalculator();
        return checksum.calculateChecksum(sequence);
    }
    
    private String getAcFromFasta(String line) {
        // >sp|P31946-2|1433B_HUMAN Isoform Short of P31946 - Homo sapiens (Human)
        int index = line.indexOf("|") + 1;
        return line.substring(index, line.indexOf("|", index));
    }
    
    public Map<String, Sequence> loadSwissProtSequences() throws IOException {
        String fileName = UNIPROT_DIR + "uniprot_sprot_human.dat";
        return loadUniProtSequences(fileName);
    }
    
    public Map<String, Sequence> loadTremblSequences() throws IOException {
        String fileName = UNIPROT_DIR + "uniprot_trembl_human.dat";
        return loadUniProtSequences(fileName);
    }
    
    private Map<String, Sequence> loadUniProtSequences(String fileName) throws IOException {
        Map<String, Sequence> accToSequence = new HashMap<String, Sequence>();
        String ac = null, sequence, checksum;
        int index;
        FileUtility fu = new FileUtility();
        String line = null;
        fu.setInput(fileName);
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("AC") && ac == null) {
                // It is possible AC has been split into several lines
                index = line.indexOf(";");
                ac = line.substring(5, index);
            }
            else if (line.startsWith("SQ")) {
                checksum = getChecksum(line);
                sequence = getSequence(fu);
                Sequence seq = new Sequence();
                seq.setChecksum(checksum);
                seq.setSequence(sequence);
                accToSequence.put("UniProt:" + ac, seq);
                // Need to reinitialization
                ac = checksum = sequence = null;
            }
        }
        fu.close();
        return accToSequence;
    }
    
    private String getSequence(FileUtility fu) throws IOException {
        String line = null;
        StringBuilder builder = new StringBuilder();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("//")) {
                // Start a new entry
                break;
            }
            builder.append(line.trim());
        }
        // Need to remove any space
        String rtn = builder.toString();
        return rtn.replaceAll(" ", "");
    }
    
    private String getChecksum(String line) {
        //SQ   SEQUENCE   245 AA;  27951 MW;  0BCA59BF97595485 CRC64;
        int index = line.lastIndexOf(";", line.length() - 2); // Pick the second last ";"
        String sub = line.substring(index + 1, line.length() - 1).trim();
        index = sub.indexOf(" ");
        return sub.substring(0, index).trim();
    }
    
    
    class Sequence {
        private String sequence;
        private String checksum;
        
        public Sequence() {
            
        }

        public String getSequence() {
            return sequence;
        }

        public void setSequence(String sequence) {
            this.sequence = sequence;
        }

        public String getChecksum() {
            return checksum;
        }

        public void setChecksum(String checksum) {
            this.checksum = checksum;
        }
        
    }
    
    
}
