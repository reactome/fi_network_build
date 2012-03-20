/*
 * Created on Apr 15, 2009
 *
 */
package org.reactome.tred;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.persistence.XMLFileAdaptor;
import org.hibernate.Session;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;

/**
 * This class is used to convert data in TRED database to a Reactome Curator Tool project.
 * @author wgm
 *
 */
public class TREDToReactomeConverter {
    private static final Logger logger = Logger.getLogger(TREDToReactomeConverter.class);
    private XMLFileAdaptor fileAdaptor;
    // Make sure each name has one EWAS
    private Map<String, GKInstance> nameToEntity;
    // Want to consolidate TF/Target based on names. Key interaction based on factor/target key
    private Map<String, GKInstance> interactionMap;
    // For literature reference
    private Map<String, GKInstance> pubmedIdToInstance;
    // To consolidate summation based on source_desc and pubmedid. It is assumed that
    // only one summation can be generated from one pubmed id. However, this may not be true.
    private Map<String, GKInstance> sourceExtractToSummation;
    // Gene id to DatabaseIdentifier to create cross link to gene pages in TRED
    private Map<Integer, GKInstance> geneIdToDatabaseIdentifier;
    private GKInstance tredDataSource;
    // For human
    private GKInstance humanSpecies;
    // For post-process
    private TREDToReactomePostProcessor postProcessor;
    // For database link
    private MySQLAdaptor dba;
    
    public TREDToReactomeConverter() {
    }
    
    private void init() throws Exception {
        fileAdaptor = new XMLFileAdaptor();
        nameToEntity = new HashMap<String, GKInstance>();
        interactionMap = new HashMap<String, GKInstance>();
        pubmedIdToInstance = new HashMap<String, GKInstance>();
        sourceExtractToSummation = new HashMap<String, GKInstance>();
        geneIdToDatabaseIdentifier = new HashMap<Integer, GKInstance>();
        tredDataSource = fileAdaptor.createNewInstance(ReactomeJavaConstants.ReferenceDatabase);
        tredDataSource.setAttributeValue(ReactomeJavaConstants.name, "TRED");
        tredDataSource.setAttributeValue(ReactomeJavaConstants.url, 
                                         "http://rulai.cshl.edu/cgi-bin/TRED/tred.cgi?process=home");
        tredDataSource.setDisplayName("TRED");
        humanSpecies = fileAdaptor.createNewInstance(ReactomeJavaConstants.Species);
        humanSpecies.setAttributeValue(ReactomeJavaConstants.name, "Homo sapiens");
    }
    
    public void setPostProcessor(TREDToReactomePostProcessor postProcessor) {
        this.postProcessor = postProcessor;
    }
    
    public void setMySQLAdaptor(MySQLAdaptor dba) {
        this.dba = dba;
    }
    
    /**
     * The main method to do converting.
     * @param destFileName
     * @throws Exception
     */
    public void convert(String destFileName) throws Exception {
        init();
        TREDHiberanteReader hibernateReader = new TREDHiberanteReader();
        Session session = hibernateReader.getSession();
        logger.info("Starting fetching FactorPromoter...");
        long time1 = System.currentTimeMillis();
        List<FactorPromoter> fps = hibernateReader.fetchHumanFactorPromoters(session);
        hibernateReader.filterHumanFactorPromoters(fps);
        long time2 = System.currentTimeMillis();
        logger.info("Done fetching FactorPromoter: " + (time2 - time1));
        for (FactorPromoter fp : fps) {
            GKInstance interction = convertFactorPromoter(fp);
        }
        session.close();
        // For post processing
        if (postProcessor != null && dba != null) {
            postProcessor.postProcess(dba, fileAdaptor);
        }
        // Do a saving
        fileAdaptor.save(destFileName);
    }
    
    /**
     * Create a TargettedInteraction from FactorPromoter.
     * @param fp
     * @return
     * @throws Exception
     */
    private GKInstance convertFactorPromoter(FactorPromoter fp) throws Exception {
        // Check attributes
        Factor factor = fp.getFactor();
        GKInstance factorInstance = createEWASBasedOnName(factor.getPrimaryName(),
                                                          factor.getAllNames());
        Promoter promoter = fp.getPromoter();
        Gene gene = promoter.getGene();
        GKInstance geneInstance = createEWASBasedOnName(gene.getPrimaryName(),
                                                        gene.getAllNames());
        String interactionKey = factorInstance.getDBID() + "-" + geneInstance.getDBID();
        GKInstance interaction = interactionMap.get(interactionKey);
        if (interaction == null) {
            // Create interaction and add attributes
            interaction = fileAdaptor.createNewInstance(ReactomeJavaConstants.TargettedInteraction);
            interaction.setAttributeValue(ReactomeJavaConstants.factor, 
                                          factorInstance);
            interaction.setAttributeValue(ReactomeJavaConstants.target,
                                          geneInstance);
            // Used as a note
            interaction.setAttributeValue(ReactomeJavaConstants.definition, 
                                          "TRED FactorPromoter id: " + fp.getId());
            interaction.setDisplayName(factorInstance.getDisplayName() + "-" + 
                                       geneInstance.getDisplayName());
            // Add a cross link to TRED web gene id
            GKInstance crossReference = createCrossReferenceForGene(gene.getId());
            interaction.addAttributeValue(ReactomeJavaConstants.crossReference, crossReference);
            interaction.setAttributeValue(ReactomeJavaConstants.species, humanSpecies);
            interactionMap.put(interactionKey, interaction);
        }
        else {
            // Add FactorPromoter Id
            String definition = (String) interaction.getAttributeValue(ReactomeJavaConstants.definition);
            definition = definition + ", " + fp.getId();
            interaction.setAttributeValue(ReactomeJavaConstants.definition, definition);
        }
        addEvidence(interaction, fp);
        return interaction;
    }
    
    private GKInstance createCrossReferenceForGene(Integer geneId) throws Exception {
        GKInstance instance = geneIdToDatabaseIdentifier.get(geneId);
        if (instance != null)
            return instance;
        instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.DatabaseIdentifier);
        instance.setAttributeValue(ReactomeJavaConstants.referenceDatabase, tredDataSource);
        instance.setAttributeValue(ReactomeJavaConstants.identifier, "Gene:" + geneId);
        geneIdToDatabaseIdentifier.put(geneId, instance);
        return instance;
    }
    
    /**
     * Extract FactorPromoterEvidence to interaction instance.
     * @param interaction
     * @param fp
     * @throws Exception
     */
    private void addEvidence(GKInstance interaction,
                             FactorPromoter fp) throws Exception {
        Set<FactorPromoterEvidence> evidences = fp.getEvidences();
        if (evidences == null) {
            logger.warn("FactorPromoter has no evidence: " + fp.getId());
            return;
        }
        for (FactorPromoterEvidence evidence : evidences) {
            // Make sure all evidences are from pubmed
            Source source = evidence.getSource();
            if (!source.getName().equals("PubMed")) {
                continue;
            }
            String pubmedId = evidence.getSourceAccession();
            // For literature reference. It will be easier to code if literture reference is attached
            // to summation. However, it is more Reactome like if this literature reference is attached
            // to interaction directly.
            GKInstance literature = createLiteratureReference(pubmedId);
            List values = interaction.getAttributeValuesList(ReactomeJavaConstants.literatureReference);
            if (values == null || !values.contains(literature))
                interaction.addAttributeValue(ReactomeJavaConstants.literatureReference, literature);
            // For summation
            GKInstance summation = createSummation(evidence);
            if (summation != null) {
                values = interaction.getAttributeValuesList(ReactomeJavaConstants.summation);
                if (values == null || !values.contains(summation))
                    interaction.addAttributeValue(ReactomeJavaConstants.summation, summation);
            }
        }
    }
    
    private GKInstance createSummation(FactorPromoterEvidence evidence) throws Exception {
        String extract = evidence.getSourceExtract();
        if (extract == null || extract.length() == 0)
            return null;
        GKInstance instance = sourceExtractToSummation.get(extract);
        if (instance != null)
            return instance;
        instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.Summation);
        sourceExtractToSummation.put(evidence.getSourceExtract(), instance);
        instance.setAttributeValue(ReactomeJavaConstants.text, extract);
        // Note: source_desc is basically artile titles plus journal information. Will not
        // bother to parse these information.
        return instance;
    }
    
    private GKInstance createLiteratureReference(String pubmedId) throws Exception {
        GKInstance instance = pubmedIdToInstance.get(pubmedId);
        if (instance != null)
            return instance;
        instance = fileAdaptor.createNewInstance(ReactomeJavaConstants.LiteratureReference);
        instance.setAttributeValue(ReactomeJavaConstants.pubMedIdentifier, 
                                   Integer.parseInt(pubmedId));
        pubmedIdToInstance.put(pubmedId, instance);
        return instance;
    }
    
    /**
     * Create an EWAS based on primary name.
     * @param primaryName
     * @param allNames
     * @return
     * @throws Exception
     */
    private GKInstance createEWASBasedOnName(String primaryName,
                                             String allNames) throws Exception {
        GKInstance ewas = nameToEntity.get(primaryName);
        if (ewas != null)
            return ewas;
        ewas = fileAdaptor.createNewInstance(ReactomeJavaConstants.EntityWithAccessionedSequence);
        ewas.addAttributeValue(ReactomeJavaConstants.name, primaryName);
        // Did a parse for all names
        if (allNames != null) {
            String[] tokens = allNames.split(", ");
            for (String token : tokens) {
               if (token.equals(primaryName))
                   continue;
               ewas.addAttributeValue(ReactomeJavaConstants.name, token);
            }
        }
        ewas.setAttributeValue(ReactomeJavaConstants.species, humanSpecies);
        InstanceDisplayNameGenerator.setDisplayName(ewas); // Need to be used by interction.
        nameToEntity.put(primaryName, ewas);
        return ewas;
    }
    
    /**
     * Use this method to fire the convert method.
     * @throws Exception
     */
    @Test
    public void doConvert() throws Exception {
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_28_plus_i",
//                                            "root",
//                                            "macmysql01",
//                                            3306);
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            FIConfiguration.getConfiguration().get("REACTOME_SOURCE_DB_NAME"),
                                            FIConfiguration.getConfiguration().get("DB_USER"),
                                            FIConfiguration.getConfiguration().get("DB_PWD"),
                                            3306);
        setMySQLAdaptor(dba);
        TREDToReactomePostProcessor postProcessor = new TREDToReactomePostProcessor();
        setPostProcessor(postProcessor);
        convert(FIConfiguration.getConfiguration().get("TRED_CONVERTED_FILE"));
    }
    
}
