/*
 * Created on Jan 9, 2009
 *
 */
package org.reactome.tred;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.fi.EntrezGeneAnalyzer;
import org.reactome.fi.util.FileUtility;
import org.reactome.funcInt.Interaction;
import org.reactome.hibernate.HibernateFIReader;

/**
 * This class is used to analyze TF targets.
 * @author wgm
 *
 */
public class TFTargetAnalyzer {
    private String SOURCE_DIR = "datasets/RegulatoryNetwork/";
    private FileUtility fu = new FileUtility();
    
    public TFTargetAnalyzer() {
    }
 
    /**
     * This method is used to compare the targets for TP53 downloaded from three 
     * TF databases.
     * @throws Exception
     */
    @Test
    public void compareP53TargetsFromDBs() throws Exception {
        String fileName = SOURCE_DIR + "TP53TargetsFromDBs.txt";
        fu.setInput(fileName);
        String line = null;
        List<String> transfacTargets = new ArrayList<String>();
        List<String> itfpTargets = new ArrayList<String>();
        List<String> tredTargets = new ArrayList<String>();
        List<String> knownTredTargets = new ArrayList<String>();
        while ((line = fu.readLine()) != null) {
            if (line.startsWith("#"))
                continue;
            if (line.startsWith(">TRANSFAC")) {
                // Process transfact targets
                while ((line = fu.readLine()) != null) {
                    if (line.length() == 0)
                        break;
                    //IN   T05038 ADA3; human, Homo sapiens.
                    String[] tokens = line.split(" +");
                    // Gene name is the third
                    String gene = tokens[2];
                    gene = gene.substring(0, gene.length() - 1);
                    transfacTargets.add(gene);
                }
            }
            else if (line.startsWith(">ITFP")) {
                // Process ITFP targets
                while ((line = fu.readLine()) != null) {
                    if (line.length() == 0)
                        break;
                    //TP53    ACLY    0.0444109   normal gene
                    String[] tokens = line.split("\t");
                    itfpTargets.add(tokens[1]);
                }
            }
            else if (line.startsWith(">TRED")) {
                while ((line = fu.readLine()) != null) {
                    if (line.length() == 0)
                        break;
                    if (line.startsWith("#"))
                        continue;
                    String[] tokens = line.split("\t");
                    // Want to get the known binding only
                    if (tokens[5].equals("known"))
                        knownTredTargets.add(tokens[0]);
                    tredTargets.add(tokens[0]);
                }
            }
        }
        EntrezGeneAnalyzer entrezAnalyzer = new EntrezGeneAnalyzer();
        Set<String> nTranTargets = entrezAnalyzer.normalizeGeneNames(transfacTargets);
        Set<String> nItfpTargets = entrezAnalyzer.normalizeGeneNames(itfpTargets);
        Set<String> nTredTargets = entrezAnalyzer.normalizeGeneNames(tredTargets);
        Set<String> nKnownTredTarSet = entrezAnalyzer.normalizeGeneNames(knownTredTargets);
        // Print out targets
        printout(nTranTargets, ">Transfac " + nTranTargets.size());
        printout(nItfpTargets, ">ITFP " + nItfpTargets.size());
        printout(nTredTargets, ">TRED " + nTredTargets.size());
        printout(nKnownTredTarSet, ">Known TRED " + nKnownTredTarSet.size());
        // Check sharing
        // Transfac and TRED
        Set<String> copyOfTran = new HashSet<String>(nTranTargets);
        copyOfTran.retainAll(nTredTargets);
        System.out.println("Shared between transfac and tred: " + copyOfTran.size());
        copyOfTran.retainAll(nItfpTargets);
        System.out.println("Shared among all three: " + copyOfTran.size());
        copyOfTran = new HashSet<String>(nTranTargets);
        copyOfTran.retainAll(nItfpTargets);
        System.out.println("Shared between transfac and itfp: " + copyOfTran.size());
        Set<String> copyOfItfp = new HashSet<String>(nItfpTargets);
        copyOfItfp.retainAll(nTredTargets);
        System.out.println("Shared between ITFP and TRED: " + copyOfItfp.size());
        Set<String> fiPartners = fetchTP53PartnersFromFIDB();
        printout(fiPartners, ">In our FI database " + fiPartners.size());
        System.out.println("Sharing between FI DB and Transfac:");
        checkSharing(fiPartners, nTranTargets);
        System.out.println("Sharing between FI DB and TRED:");
        checkSharing(fiPartners, nTredTargets);
        System.out.println("Sharing between FI DB and ITFG:");
        checkSharing(fiPartners, nItfpTargets);
        System.out.println("Sharing between FI DB and Known TRED:");
        checkSharing(fiPartners, nKnownTredTarSet);
    }
    
    private void checkSharing(Set<String> source,
                              Set<String> target) {
        Set<String> copy = new HashSet<String>(source);
        copy.retainAll(target);
        System.out.println("Shared genes: " + copy.size());
        for (String gene : copy)
            System.out.println("\t" + gene);
    }
    
    private void printout(Collection<String> list, String title) {
        System.out.println(title);
        //for (String gene : list)
        //    System.out.println(gene);
        //System.out.println();
    }
    
    @Test
    public void checkTP53Sequences() throws IOException {
        String fileName = SOURCE_DIR + "TwoTP53.txt";
        fu.setInput(fileName);
        String seq1 = fu.readLine();
        char[] bp1 = seq1.toCharArray();
        String seq2 = fu.readLine();
        char[] bp2 = seq2.toCharArray();
        for (int i = 0; i < bp1.length; i++) {
            char c1 = bp1[i];
            char c2 = bp2[i];
            if (c1 != c2)
                System.out.println(i + ": " + c1 + " -> " + c2);
        }
        fu.close();
    }
    
    private Set<String> fetchTP53PartnersFromFIDB() throws Exception {
        HibernateFIReader reader = new HibernateFIReader();
        SessionFactory sf = reader.initSession();
        Session session = sf.openSession();
        List<Interaction> interactions = reader.queryFIsForName("TP53", session);
        System.out.println("Total interactions: " + interactions.size());
        // Take proteins only from Human PPI and known ones
        Set<String> proteins = new HashSet<String>();
        for (Interaction i : interactions) {
            if (i.getEvidence() == null) {
                proteins.add(i.getFirstProtein().getShortName());
                proteins.add(i.getSecondProtein().getShortName());
            }
            else if (i.getEvidence().getHumanInteraction() != null &&
                     i.getEvidence().getHumanInteraction()) {
                proteins.add(i.getFirstProtein().getShortName());
                proteins.add(i.getSecondProtein().getShortName());
            }
        }
        session.close();
        proteins.remove("TP53");
        System.out.println("Total proteins: " + proteins.size());
        return proteins;
    }
    
}
