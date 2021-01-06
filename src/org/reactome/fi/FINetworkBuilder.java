/*
 * Created on Mar 20, 2012
 *
 */
package org.reactome.fi;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.junit.Test;
import org.reactome.b2rPostProcessor.NciPIDConverterRunner;
import org.reactome.data.EncodeTFTargetToReactomeConverter;
import org.reactome.data.EnsemblAnalyzer;
import org.reactome.data.GODataAnalyzerV2;
import org.reactome.data.GOTermLoader;
import org.reactome.data.IRefIndexMITTabAnalyzer;
import org.reactome.data.MicroarrayDataAnalyzer;
import org.reactome.data.PfamAnalyzer;
import org.reactome.data.ReactomeAnalyzer;
import org.reactome.data.UniProtAnalyzer;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.hibernate.HibernateFIReader;
import org.reactome.kegg.KeggToReactomeConverter;
import org.reactome.panther.PantherToReactomeConverterTest;
import org.reactome.psi.data.PsiMiOrthologyAnalyzer;
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
 * 8). generateCytoscapePlugInFiles()
 * 9). compareFilesToPreviousVersion()
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
        // For GO annotation. As of 2016, go.obo is downloaded and terms_and_ids.txt
        // should be generated before other procedures.
        GOTermLoader loader = new GOTermLoader();
        logger.info("Generating GO.terms_and_ids.txt from go.obo...");
        loader.generateGOTermsToIdsFromObo();
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
        
        TREDToReactomeConverter tredConverter = new TREDToReactomeConverter();
        logger.info("Running TREDToReactomeConverter.doConvert()...");
        tredConverter.doConvert();

        // Gene expressions have been used as support evidences for TF/Target interactions
        // from ENCODE. So have to call these methods first in order to get co-expression
        // values
        // These two methods will be called again in method prepareNBCFeatures(). The results
        // should be the same!
        MicroarrayDataAnalyzer microarrayAnalyzer = new MicroarrayDataAnalyzer();
        logger.info("Running MicroarrayDataAnalyzer.normalizeLeeGeneExp()...");
        microarrayAnalyzer.normalizeLeeGeneExp();
        logger.info("Running MicroarrayDataAnalyzer.generatePrietoCarlosGeneExpFile()...");
        microarrayAnalyzer.generatePrietoCarlosGeneExpFile();
        
        logger.info("Running EncodeTFTargetToReactomeConverter.convert()...");
        new EncodeTFTargetToReactomeConverter().convert();

        PantherToReactomeConverterTest pantherConverter = new PantherToReactomeConverterTest();
        logger.info("Running PantherToReactomeConverterTest.testNewBatchConverter()...");
        pantherConverter.testNewBatchConverter();

        // Check protein coverage
        // Some files have not be generated in the following call. This method should not be called.
        ProteinAndInteractionCount count = new ProteinAndInteractionCount();
        logger.info("Count proteins...");
        count.checkUniProtNumbersInConvertedDBs();

    }
    
    /**
     * Dump converted pathways from other databases into the extended Reactome database.
     * Before running method, try to take a look at the converted projects to make sure
     * nothing weird there.
     * @throws Exception
     */
    @Test
    public void dumpPathwayDBs() throws Exception {
        ConvertedPathwayDbDumper dumper = new ConvertedPathwayDbDumper();
        logger.info("Running ConvertedPathwayDbDumper.dumpPathwayDBs()...");
        dumper.dump();
    }
    
    /**
     * Extract Pathway FIs from Reactome and other imported pathway dbs, and dump
     * them into a local file. Before running this method, make sure you have set
     * dataSourceIds correctly in method ReactomeAnalyzer.getPathwayDBAnalyzers().
     * @throws Exception
     */
    @Test
    public void dumpPathwayFIs() throws Exception {
        logger.info("Running FIFileAnalyzer.dumpPathwayFIs()...");
        FIFileAnalyzer fiFileAnalyzer = new FIFileAnalyzer();
        fiFileAnalyzer.dumpPathwayFIs();
    }
    
    /**
     * Run this method to prepare features for NBC training.
     * @throws Exception
     */
    @Test
    public void prepareNBCFeatures() throws Exception {
        // The following statements are used for PPI based features
        EnsemblAnalyzer ensemblAnalyzer = new EnsemblAnalyzer();
        logger.info("Running EnsemblAnalyzer.dumpProteinFamilies()...");
        ensemblAnalyzer.dumpProteinFamilies();
        // Extract PPIs from PSITAB files
        IRefIndexMITTabAnalyzer iRefIndexAnalyzer = new IRefIndexMITTabAnalyzer();
        logger.info("Running IRefIndexMITTabAnalyzer.loadHumanPPIs()...");
        iRefIndexAnalyzer.loadHumanPPIs();
        logger.info("Running IRefIndexMITTabAnalyzer.loadFlyPPIs()...");
        iRefIndexAnalyzer.loadFlyPPIs();
        logger.info("Running IRefIndexMITTabAnalyzer.loadWormPPIs()...");
        iRefIndexAnalyzer.loadWormPPIs();
        logger.info("Running IRefIndexMITTabAnalyzer.loadYeastPPIs()...");
        iRefIndexAnalyzer.loadYeastPPIs();
        logger.info("Running IRefIndexMITTabAnalyzer.loadMousePPIs()...");
        iRefIndexAnalyzer.loadMousePPIs();
        // Map PPIs in other species to human PPIs via ensembl compara
        PsiMiOrthologyAnalyzer psimiAnalyzer = new PsiMiOrthologyAnalyzer();
        logger.info("Running PsiMiOrthologyAnalyzer.generateHumanPPIsFromYeastInUniProt()...");
        psimiAnalyzer.generateHumanPPIsFromYeastInUniProt();
        logger.info("Running PsiMiOrthologyAnalyzer.generateHumanPPIsFromWormInUniProt()...");
        psimiAnalyzer.generateHumanPPIsFromWormInUniProt();
        logger.info("Running PsiMiOrthologyAnalyzer.generateHumanPPIsFromFlyInUniProt()...");
        psimiAnalyzer.generateHumanPPIsFromFlyInUniProt();
        logger.info("Running PsiMiOrthologyAnalyzer.generateHumanPPIsFromMouseInUniProt()...");
        psimiAnalyzer.generateHumanPPIsFromMouseInUniProt();
        // Check PPIs' odds ratioes
        logger.info("Checking odds ratio...");
        psimiAnalyzer.testPPIsFeatures();
        
        // The following statements are used for gene expression based features: two
        // data sets have been used here.
        MicroarrayDataAnalyzer microarrayAnalyzer = new MicroarrayDataAnalyzer();
        logger.info("Running MicroarrayDataAnalyzer.normalizeLeeGeneExp()...");
        microarrayAnalyzer.normalizeLeeGeneExp();
        logger.info("Checking its odds ratio...");
        microarrayAnalyzer.checkCoExpFromPavlidis();
        logger.info("Running MicroarrayDataAnalyzer.generatePrietoCarlosGeneExpFile()...");
        microarrayAnalyzer.generatePrietoCarlosGeneExpFile();
        logger.info("Check its odds ratio...");
        microarrayAnalyzer.checkCoExpFromPrietoCarlos();
        
        // Domain-domain interactions
        PfamAnalyzer pfamAnalyzer = new PfamAnalyzer();
        logger.info("Running PfamAnalyzer.convertIntToPfamIDs()...");
        pfamAnalyzer.convertIntToPfamIDs();
        logger.info("Checking its odds ratio...");
        pfamAnalyzer.testPfamFeature();
        
        // Nothing needs to be generated. But we want to check GO features
        GODataAnalyzerV2 goAnalyzer = new GODataAnalyzerV2();
        logger.info("Checking GO features...");
        goAnalyzer.testGOFeatures();
    }
    
    /**
     * Train NBC, create data for ROC curve, check protein coverge and generate
     * a predicted FI file. The data file generated from calcualteROCPoints() should
     * be used in R for generating ROC curve and calculate AUC.
     * @throws Exception
     */
    @Test
    public void trainNBC() throws Exception {
        NBCAnalyzer nbcAnalyzer = new NBCAnalyzer();
        logger.info("Running NBCAnalyzer.calculateNBCBasedOnReactome()...");
        nbcAnalyzer.calculateNBCBasedOnReactome();
        logger.info("Running NBCAnalyzer.calculateROCPoints()...");
        nbcAnalyzer.calculateROCPoints();
        // Generate protein pairs having domain-domain interactions 
        logger.info("Running NBCAnalyzer.checkSharedBPPairAndDomainPairs()...");
        nbcAnalyzer.checkSharedBPPairAndDomainPair();
        logger.info("Running checkCutoffValueForPredictedFIs()...");
        nbcAnalyzer.checkCutoffValueForPredictedFIs();
        // At this stage, NBCGUITest() should be run to check contributions of each
        // feature and see if anything is weird!
        // This method should be run separately
        // NBCGUITest.main()
    }
    
    /**
     * After investigating the results from method trainNBC(), choose an appropriate 
     * cutoff value, and set in the configuration file: results/configuration.prop.
     * Usually 0.50 or above should be chosen unless you have a strong reason not so.
     * @throws Exception
     */
    @Test
    public void predictFIs() throws Exception {
        NBCAnalyzer nbcAnalyzer = new NBCAnalyzer();
        logger.info("Running NBCAnalyzer.generatePredictedFIs()...");
        nbcAnalyzer.generatePredictedFIs();
    }
    
    /**
     * Generate a hibernate based FI database and dump FIs in gene names into 
     * files.
     * @throws Exception
     */
    @Test
    public void buildFIDb() throws Exception {
        // Need to generate database schema first. To do there, an
        // empty database should be available first.
        FIDBBuilder dbBuilder = new FIDBBuilder();
        logger.info("Running FIDBBuilder.generateSchema()...");
        dbBuilder.generateSchema();
        // Dump extracted pathway FIs into the database
        logger.info("Running FIDBBuilder.dump()..."); 
        dbBuilder.dump();
        // Dump predicted FIs into the database
        logger.info("Running FIDBBuilder.dumpPredicted()...");
        dbBuilder.dumpPredictedFIs();
    }
    
    /**
     * Use this method to generate a list of static files to be used by Cytoscape plug-in.
     * @throws Exception
     */
    @Test
    public void generateCytoscapePlugInFiles() throws Exception {
        // These two files are used for reaction-based enrichment analysis
        ReactomeAnalyzer reactomeAnalyzer = new ReactomeAnalyzer();
        logger.info("Running ReactomeAnalyzer.generateGenesToReactionsMap()...");
        reactomeAnalyzer.generateGenesToReactionsMap();
        logger.info("Running ReactomeAnalyzer.generateFIsToReactionsMap()...");
        reactomeAnalyzer.generateFIsToReactionsMap();
//        if (true)
//            return;
        // Generate a file containing FIs in gene names so that it can be used in 
        // plug-in and data analysis
        HibernateFIReader hibernateReader = new HibernateFIReader();
        logger.info("Running HibernateFIReader.generateFIFileInGeneInHiberante()...");
        hibernateReader.generateFIFileInGeneInHibernate();
//        if (true)
//            return;
        // Create a map from pathway fis to their original sources
        logger.info("Running HibernateFIReader.generateFIInGeneSourceFile()...");
        hibernateReader.generateFIInGeneSourceFile();

        // Generate the largest graph component
        FIGraphAnalyzer graphAnalyzer = new FIGraphAnalyzer();
        logger.info("Running FIGraphAnalyzer.analyzeComponents()...");
        graphAnalyzer.analyzeComponents();
        // Mapping from protein accessions to names
        // This file is not used here right now. But it is used in the plug-in.
        logger.info("Running HiberanteFIReader.generateAccessionToProteinNameMap()...");
        hibernateReader.generateAccessionToProteinNameMap();
        PathwayGeneSetGenerator genesetGenerator = new PathwayGeneSetGenerator();
        logger.info("Running PathwayGeneSetGenerator.generateReactomePathwayListBasedOnDiagrams()...");
        // Generate a flattened list of pathways from the Reactome database. Each Reactome pathway should have
        // an associated pathway diagram fully laid-out.
        genesetGenerator.generateReactomePathwayListBasedOnDiagrams();
        // Generate Pathways to Genes mappings for enrichment analysis
        logger.info("Running PathwayGeneSetGenerator.generateProteinNameToPathwayMap()...");
        genesetGenerator.generateProteinNameToPathwayMap();
        // Generate a list of all pathways in the Reactome database. This list of pathways is used for 
        // pathway enrichment analysis in a hierarchical way (aka not based on flatenned list)
        // Generate GMT file for GSEA analysis
        genesetGenerator.generateHumanFiles();
        // For mouse files
        genesetGenerator.generateMouseFiles();
        // Dump MySQL databases to MyISAM types so that they can be loaded into the deployment machine
        // or be backed up
        MySQLDatabaseHandler dbHandler = new MySQLDatabaseHandler();
        dbHandler.dumpDatabases();
        // Copy files needed by the plug-in web application
        CytoscapePlugInFileHandler fileHandler = new CytoscapePlugInFileHandler();
        fileHandler.generatePlugInFiles();
        // The following statements are used to generate the matrix file for the biggest FI network component.
        // The generated matrix should be moved to the OICR to calculate the heat kernel for the HotNet implementation
        // using R (using script runHeatKernel_R.sh, which should be modified for the correct file names). 
        // The generated matrix text file should NOT be copied into the WEB-INF folder.
//        HotNetMatrixCalculator hotnetMatrixCalculator = new HotNetMatrixCalculator();
//        hotnetMatrixCalculator.testCalculateHeatKernel();
        // The following method should be called after the kernel file was generated from using R
//        HotNetMatrixCalculator hotnetMatrixCalculator = new HotNetMatrixCalculator();
//        hotnetMatrixCalculator.generateSerializedMatrixFile();
        
        // Generate factor graphs for pathway diagrams to be used by the Cytoscape app
        // Consider running this with its own step to avoid a long logging.
        FactorGraphDumper fgDumper = new FactorGraphDumper();
        fgDumper.dump();
        // Don't forget to generate FIsInGene_xxxxxx_with_annotations.txt file by using
        // the method in the fiviz_corews project: 
        // The method should be run is in class org.reactome.r3.fi.InteractionAnnotator.annotateAllFIs()
        // in project FIVizWS_corews. This should be moved to here in the future.
//        logger.info("Running InteractionAnnotator.annotateAllFIs()...");
//        InteractionAnnotator annotator = new InteractionAnnotator();
//        annotator.annoateAllFIs();
    }
    
    /**
     * This method is used to generate a comparison file to compare files
     * generated for this version to a previous version. The results are better
     * to be viewed in Excel, and copied to the combined_logging.txt file.
     * @throws IOException
     */
    @Test
    public void compareFilesToPreviousVersion() throws IOException {
        FIConfiguration config = FIConfiguration.getConfiguration();
        String year = config.get("YEAR");
        String preYear = (new Integer(year) - 1) + "";
        File dir = new File(config.get("RESULT_DIR"));
        Map<String, File> nameToFile = getNameToFile(dir);
        
        // After migrating to git. For 2017 only.
//        String preDirName = "/Users/wug/Documents/eclipse_workspace/FINetworkBuild_CVS/results/";
//        File preDir = new File(preDirName + preYear);
        
        File preDir = new File("results/" + preYear);
        Map<String, File> preNameToFile = getNameToFile(preDir);
        Set<String> allNames = new HashSet<String>(nameToFile.keySet());
        allNames.addAll(preNameToFile.keySet());
        List<String> nameList = new ArrayList<String>(allNames);
        Collections.sort(nameList);
        // Just print out files in the system out for easy view
        System.out.println("\t" + year + "\t" + preYear);
        for (String name : nameList) {
            if (name.startsWith("."))
                continue; // Don't show hidden file
            File file = nameToFile.get(name);
            File preFile = preNameToFile.get(name);
            System.out.println(name + "\t" + 
                               (file == null ? "NA" : file.length()) + "\t" + 
                               (preFile == null ? "NA" : preFile.length()));
        }
    }
    
    private Map<String, File> getNameToFile(File dir) {
        Map<String, File> nameToFile = new HashMap<String, File>();
        File[] files = dir.listFiles();
        for (File file : files) 
            nameToFile.put(file.getName(), file);
        return nameToFile;
    }
    
}
