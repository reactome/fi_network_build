/*
 * Created on Apr 10, 2009
 *
 */
package org.reactome.tred;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.PropertyConfigurator;
import org.hibernate.Query;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FileUtility;
import org.reactome.fi.util.HibernateUtil;

/**
 * Use a Hibernate API to query the database.
 * @author wgm
 *
 */
public class TREDHiberanteReader {
    private SessionFactory sf;
    private FileUtility fu;
    
    public TREDHiberanteReader() {
        fu = new FileUtility();
    }
    
    private void setup() throws Exception {
        // For logging for Java
        PropertyConfigurator.configure("resources/log4j.properties");
        String configFileName = "resources/TREDHibernate.cfg.xml";
        File configFile = new File(configFileName);
        sf = HibernateUtil.getSessionFactory(configFile);
    }
    
    public Session getSession() throws Exception {
        if (sf == null)
            setup();
        return sf.openSession();
    }
    
    /**
     * This method is used to dump human TF and target interactions.
     * @throws Exception
     */
    @Test
    public void dumpHumanTFTargetInterctions() throws Exception {
        setup();
        Session session = sf.openSession();
        List<FactorPromoter> fps = fetchHumanFactorPromoters(session);
        filterHumanFactorPromoters(fps);
        // Hold TF to targets 
        Map<String, Set<String>> factorToGenes = new HashMap<String, Set<String>>();
        for (FactorPromoter fp : fps) {
            Factor factor = fp.getFactor();
            Set<String> targetGenes = factorToGenes.get(factor.getPrimaryName());
            if (targetGenes == null) {
                targetGenes = new HashSet<String>();
                factorToGenes.put(factor.getPrimaryName(), targetGenes);
            }
            Gene gene = fp.getPromoter().getGene();
            targetGenes.add(gene.getPrimaryName());
        }
        // Generate a file
        // For high-quality TF/Target interactions: FactorPromoterQuality 1, 2 and 
        // PromoterQuality 1, 2
        String outFileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "TREDTFTargetInteractions.txt";
        // For all TF/Target interactions
        //String outFileName = R3Constants.RESULT_DIR + "TREDTFTargetInteractions_All.txt";
        fu.saveSetMap(factorToGenes, outFileName);
    }
    
    /**
     * Filter a list of human FactorPromoter based on qualities.
     * @param fps
     */
    public void filterHumanFactorPromoters(List<FactorPromoter> fps) {
        for (Iterator<FactorPromoter> it = fps.iterator(); it.hasNext();) {
            FactorPromoter fp = it.next();
            if (!acceptFactorPromoter(fp)) {
                it.remove();
                continue;
            }
            Promoter promoter = fp.getPromoter();
            if (!acceptPromoter(promoter)) {
                it.remove();
                continue;
            }
            Gene gene = promoter.getGene();
            // Check species
            if (!gene.getSpecies().equals("human, Homo sapiens")) {
                it.remove();
                continue; // This is very strange. It should not be checked. It is more like a bug in the database.
            }
        }
    }
    
    /**
     * Fetch the list of human FactorPromoters, which contain information about TF/Target interactions.
     * @param session
     * @return
     * @throws Exception
     */
    public List<FactorPromoter> fetchHumanFactorPromoters(Session session) throws Exception {
        List<FactorPromoter> rtn = new ArrayList<FactorPromoter>();
        String queryString = "FROM " + Factor.class.getName() + " f WHERE f.species = ?";
        // Query for human Factors
        Query query = session.createQuery(queryString).setString(0, "human, Homo sapiens");
        List list = query.list();
        // Query for factor promoters
        queryString = "FROM " + FactorPromoter.class.getName() + " fp WHERE fp.factor = ?";
        // Hold TF to targets 
        Map<String, Set<String>> factorToGenes = new HashMap<String, Set<String>>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Factor factor = (Factor) it.next();
            query = session.createQuery(queryString).setParameter(0, factor);
            List fps = query.list(); 
            if (fps == null || fps.size() == 0)
                continue;
            rtn.addAll(fps);
        }
        return rtn;
    }
    
    private boolean acceptPromoter(Promoter promoter) {
        PromoterQuality promoterQuality = promoter.getQuality();
        if (promoterQuality.getId() == 1 ||
            promoterQuality.getId() == 2) {
            return true;
        }
        return false;
//        return true;
    }
   
    private boolean acceptFactorPromoter(FactorPromoter fp) {
        FactorPromoterQuality quality = fp.getQuality();
        String qualityName = quality.getName();
        if (qualityName.equals("known") ||
            qualityName.equals("likely")) {
            return true;
        }
        return false;
//        return true;
    }
    
    @Test
    public void testQuery() throws Exception {
        setup();
        Session session = sf.openSession();
        String queryString = "FROM " + Factor.class.getName() + " f WHERE f.species = ?";
        Query query = session.createQuery(queryString).setString(0, "human, Homo sapiens");
        List list = query.list();
        int index = 1;
        Set<String> factorNames = new HashSet<String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Factor factor = (Factor) it.next();
            factorNames.add(factor.getPrimaryName());
            //System.out.println(index + ": " + factor.getPrimaryName());
            index ++;
        }
        System.out.println("Total human factors in the DB: " + list.size());
        System.out.println("Total human factor names in the DB: " + factorNames.size());
        // Get its promoters
        queryString = "FROM " + FactorPromoter.class.getName() + " fp WHERE fp.factor = ?";
        Factor testFactor = null;
        for (Iterator it = list.iterator(); it.hasNext();) {
            Factor factor = (Factor) it.next();
            query = session.createQuery(queryString).setParameter(0, factor);
            List list1 = query.list();
            System.out.println(factor.getPrimaryName() + ": " + list1.size());
            if (factor.getPrimaryName().equals("E2F-4"))
                testFactor = factor;
        }
        query = session.createQuery(queryString).setParameter(0, testFactor);
        list = query.list();
        index = 0;
        System.out.println("Test factor: " + testFactor.getPrimaryName());
        Set<String> genes = new HashSet<String>();
        for (Iterator it = list.iterator(); it.hasNext();) {
            FactorPromoter fp = (FactorPromoter) it.next();
            if (acceptFactorPromoter(fp)) {
                Promoter promoter = fp.getPromoter();
                PromoterQuality promoterQuality = promoter.getQuality();
                if (promoterQuality.getId() == 1 ||
                    promoterQuality.getId() == 2) {
                    Gene gene = promoter.getGene();
                    // Check gene species
                    if (gene.getSpecies().equals("human, Homo sapiens")) {
                        index ++;
                        System.out.println(index + ": " + gene.getPrimaryName() + ", " + promoter.getId());
                        genes.add(gene.getPrimaryName());
                    }
                }
                //System.out.println(index + ": " + promoter.getName() + ", " + gene.getPrimaryName());
            }
        }
        System.out.println("Total factor promoters with known or likely qualities: " + index);
        System.out.println("Total genes: " + genes.size());
        session.close();
    }
    
    public Map<String, Set<String>> loadPrimaryNameToAllNames() throws Exception {
        setup();
        Session session = sf.openSession();
        Map<String, Set<String>> primaryNameToAllNames = new HashMap<String, Set<String>>();
        String queryString = "SELECT f.primaryName, f.allNames FROM " + Factor.class.getName() + " f WHERE f.species = ?";
        // Query for human Factors
        Query query = session.createQuery(queryString).setString(0, "human, Homo sapiens");
        List list = query.list();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Object[] names = (Object[]) it.next();
            String primaryName = names[0].toString();
            String allNames = names[1].toString();
            String[] tokens = allNames.split(", ");
            Set<String> set = new HashSet<String>();
            for (String t : tokens)
                set.add(t);
            primaryNameToAllNames.put(primaryName, set);
        }
        queryString = "SELECT g.primaryName, g.allNames FROM " + Gene.class.getName() + " g WHERE g.species = 'human, Homo sapiens' AND g.primaryName != 'n/a'";
        list = query.list();
        for (Iterator it = list.iterator(); it.hasNext();) {
            Object[] names = (Object[]) it.next();
            String primaryName = names[0].toString();
            String allNames = names[1].toString();
            String[] tokens = allNames.split(", ");
            Set<String> set = new HashSet<String>();
            for (String t : tokens)
                set.add(t);
            primaryNameToAllNames.put(primaryName, set);
        }
        session.close();
        return primaryNameToAllNames;
    }
}
