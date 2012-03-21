/*
 * Created on Mar 20, 2012
 *
 */
package org.reactome.fi;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.junit.Test;
import org.reactome.b2rPostProcessor.NciPIDConverterRunner;
import org.reactome.data.EnsemblAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.kegg.KeggToReactomeConverter;
import org.reactome.panther.PantherToReactomeConverterTest;
import org.reactome.tred.TREDToReactomeConverter;

/**
 * This class groups a list of methods that should be run for building a new version of FI network.
 * Here is the sequence to run these methods:
 * 1). prepareMappingFiles()
 * 2). convertPathwayDBs()
 * 3). dumpPathwayDBs()
 * 4). prepareNBCFeatures()
 * 5). trainNBC()
 * 6). predictFIs()
 * 7). buildFIDb()
 * @author gwu
 *
 */
public class FINetworkBuilder {
    private static Logger logger = Logger.getLogger(FINetworkBuilder.class);
    
    public FINetworkBuilder() {
        PropertyConfigurator.configure("resources/log4j.properties");
    }
    
    /**
     * Running this method to prepare several files that are used to do mapping.
     * This method should be run first.
     * @throws Exception
     */
    @Test
    public void prepareMappingFiles() throws Exception {
        UniProtAnalyzer uniProtAnalyzer = new UniProtAnalyzer();
        logger.info("Running UniProtAnalyzer.generateUniProtIDsMap()...");
        uniProtAnalyzer.generateUniProtIDsMap();
        logger.info("Running UniProtAnalyzer.generateEntrezGeneToUniProt()...");
        uniProtAnalyzer.generateEntrezGeneToUniProt();
        logger.info("Running UniProtAnalyzer.generateUniToPfamMap()...");
        uniProtAnalyzer.generateUniToPfamMap();
        EnsemblAnalyzer ensemblAnalyzer = new EnsemblAnalyzer();
        logger.info("Running EnsemblAnalyzer.dumpProteinFamilies()...");
        ensemblAnalyzer.dumpProteinFamilies();
    }
    
    /**
     * Convert non-Reactome pathway databases into a local Reactome curator tool
     * projects before dumping into the merged database.
     * @throws Exception
     */
    @Test
    public void convertPathwayDBs() throws Exception {
        KeggToReactomeConverter keggConverter = new KeggToReactomeConverter();
        logger.info("Running KeggToReactomeConverter.runBatchConver()...");
        keggConverter.runBatchConvert();
        NciPIDConverterRunner nciPidConverter = new NciPIDConverterRunner();
        logger.info("Running NciPIDConverrerRunner.runConvertOfCurated()...");
        nciPidConverter.runConvertOfCurated();
        logger.info("Running NciPIDConverterRunner.runConvertOfBiocarta()...");
        nciPidConverter.runConvertOfBioCarta();
        PantherToReactomeConverterTest pantherConverter = new PantherToReactomeConverterTest();
        logger.info("Running PantherToReactomeConverterTest.testNewBatchConverter()...");
        pantherConverter.testNewBatchConverter();
        TREDToReactomeConverter tredConverter = new TREDToReactomeConverter();
        logger.info("Running TREDToReactomeConverter.doConvert()...");
        tredConverter.doConvert();
        // Check protein coverage
        ProteinAndInteractionCount count = new ProteinAndInteractionCount();
        logger.info("Count proteins...");
        count.checkUniProtNumbersInConvertedDBs();
    }
    
    public void dumpPathwayDBs() throws Exception {
        
    }
    
    public void prepareNBCFeatures() throws Exception {
        
    }
    
    public void trainNBC() throws Exception {
        
    }
    
    public void predictFIs() throws Exception {
        
    }
    
    public void buildFIDb() throws Exception {
        
    }
}
