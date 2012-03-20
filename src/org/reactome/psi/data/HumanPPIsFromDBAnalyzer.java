/*
 * Created on Apr 1, 2009
 *
 */
package org.reactome.psi.data;

import java.util.Set;

import org.gk.model.GKInstance;
import org.junit.Test;
import org.reactome.data.CPathAnalyzer;
import org.reactome.data.ProteinIdFilters;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.fi.util.FeatureChecker;
import org.reactome.fi.util.FileUtility;

/**
 * This class is used to handle PPIs that have been dumped into the Reactome database already.
 * @author wgm
 *
 */
public class HumanPPIsFromDBAnalyzer {
    
    public HumanPPIsFromDBAnalyzer() {
    }
    
    /**
     * Source IDs for human PPI databases.
     * @return
     */
    private Long[] getDataSourceIds() {
        Long[] sourceIds = new Long[] {
                // Based on reactome_release_28_plus_i.
                616255L, // IntAct
                //816565L, // HPRD
                //724811L  // BioGrid
        };
        return sourceIds;
    }
    
    /**
     * This method is used to check human PPIs from the PPI databases
     * used as feature in NBC based on odds ratio.
     * @throws Exception
     */
    @Test
    public void testPPIAsFeature() throws Exception {
        String[] fileNames = new String[] {
                "HumanPPIs_Less2_intact.txt",
                "HumanPPIs_HPRD.txt",
                "HumanPPIs_BioGrid.txt"
        };
        FeatureChecker checker = new FeatureChecker();
        FileUtility fu = new FileUtility();
        for (String fileName : fileNames) {
            Set<String> ppis = fu.loadInteractions(FIConfiguration.getConfiguration().get("RESULT_DIR") + fileName);
            System.out.println("File: " + fileName);
            System.out.println("Total PPIs: " + ppis.size());
            checker.checkFeatureOddsRatio(ppis);
        }
    }
    
    /**
     * This method is used to dump human PPIs in UniProt accession
     * numbers into files.
     * @throws Exception
     */
    @Test
    public void dumpPPIsInUniProts() throws Exception {
        // Used to filer to human UniProt accession numbers.
        ProteinIdFilters idFilter = new ProteinIdFilters();
        FileUtility fu = new FileUtility();
        Long[] sourceIds = getDataSourceIds();
        for (Long sourceId : sourceIds) {
            CPathAnalyzer analyzer = new CPathAnalyzer();
            analyzer.setDataSourceId(sourceId);
            GKInstance dataSource = analyzer.getDataSource();
            Set<String> ppis = analyzer.extractInteractionSet();
            Set<String> filtered = idFilter.cleanUpVsUniProt(ppis);
            System.out.println("Total PPIs in " + dataSource.getDisplayName() + ": " + filtered.size());
            String fileName = FIConfiguration.getConfiguration().get("RESULT_DIR") + "HumanPPIs_Less2_" + dataSource.getDisplayName() + ".txt";
            fu.saveInteractions(filtered, fileName);
        }
    }
}
