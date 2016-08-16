/*
 * Created on Jun 17, 2010
 * Moved from caBigR3Web project to the ReactomeFINetwork build project for generating the annotated
 * FI file.
 */
package org.reactome.hibernate;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.hibernate.Query;
import org.hibernate.Session;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.junit.Test;
import org.reactome.fi.util.FIConfiguration;
import org.reactome.funcInt.FIAnnotation;
import org.reactome.funcInt.Interaction;
import org.reactome.funcInt.ReactomeSource;
import org.reactome.funcInt.ReactomeSourceType;
import org.reactome.r3.util.FileUtility;

/**
 * This class is used to annotate functional interactions.
 * @author wgm
 *
 */
public class InteractionAnnotator extends HibernateFIReader {
    private static final Logger logger = Logger.getLogger(InteractionAnnotator.class);
    private MySQLAdaptor dba;
    // Load types
    private Map<String, FIAnnotation> nameToType;
    // We use this as a map to search of type 
    private Map<String, String> reverseNameToName;
    // For return reactome sources for pathway fis
    private Map<String, Set<Long>> fiToSources;
    // A flag to limit annotations to the Reactome data source only
    private boolean useReactomeSourceOnly = false; // Default should be false
    
    public InteractionAnnotator() throws Exception {
        initInteractionTypes();
    }
    
    /**
     * Query Reactome sources DB IDs for a set of FIs.
     * @param fis
     * @return
     * @throws IOException
     */
    public Map<String, Set<Long>> queryPathwayFIsSources(String[] fis) {
        if (fiToSources == null)
            return null;
        Map<String, Set<Long>> rtnFiToSource = new HashMap<String, Set<Long>>();
        for (String fi : fis) {
            Set<Long> sources = fiToSources.get(fi);
            if (sources == null) {
                // Do a flip in case proteins in FIs are not sorted
                String tmpFI = flipFI(fi);
                sources = fiToSources.get(tmpFI);
                if (sources == null)
                    continue;
            }
            rtnFiToSource.put(fi, sources);
        }
        return rtnFiToSource;
    }
    
    private String flipFI(String fi) {
        int index = fi.indexOf("\t");
        String protein1 = fi.substring(0, index);
        String protein2 = fi.substring(index + 1);
        return protein2 + "\t" + protein1;
    }
    
    public void setSourceDBA(MySQLAdaptor dba) {
        this.dba = dba;
        dba.initDumbThreadForConnection();
    }
    
    /**
     * Query the type for the specified two proteins. There should be only one type returned.
     * So this method has merged types for interactions name1:name2 and name2:name1.
     * @param name1
     * @param name2
     * @param session
     * @return
     * @throws Exception
     */
    public FIAnnotation annotate(String name1,
                                 String name2,
                                 Query query) throws Exception {
        if (nameToType == null)
            initInteractionTypes();
        List<Interaction> interactions = queryFIsBasedOnGeneNames(name1, 
                                                                  name2,
                                                                  query);
        if (interactions == null || interactions.size() == 0) {
            logger.error(name1 + ", " + name2 + " cannot find an interaction!");
            return null;
        }
        return annotate(name1, 
                        name2,
                        interactions);
    }
    
    /**
     * Use this method to control the use of query for faster performance.
     * @param name1
     * @param name2
     * @param query
     * @return
     * @throws Exception
     */
    private List<Interaction> queryFIsBasedOnGeneNames(String name1, 
                                                       String name2,
                                                       Query query) throws Exception {
        List<Interaction> rtn = new ArrayList<Interaction>();
        query.setString(0, name1);
        query.setString(1, name2);
        List interactions = query.list();
        if (interactions != null && interactions.size() > 0)
            rtn.addAll(interactions);
        // Do a reverse search
        query.setString(0, name2);
        query.setString(1, name1);
        interactions = query.list();
        if (interactions != null && interactions.size() > 0)
            rtn.addAll(interactions);
        return rtn;
    }
    
    /**
     * Annotate an Interaction object.
     * @param interaction
     * @return
     * @throws Exception
     */
    public FIAnnotation annotate(Interaction interaction) throws Exception {
        Set<FIAnnotation> annotations = queryType(interaction,
                                                  interaction.getFirstProtein().getShortName(),
                                                  interaction.getSecondProtein().getShortName());
        return mergeTypes(annotations);
    }

    private FIAnnotation annotate(String name1, 
                                  String name2,
                                  List<Interaction> interactions) throws Exception {
        Set<FIAnnotation> rtnTypes = new HashSet<FIAnnotation>();
        for (Interaction interaction : interactions) {
            Set<FIAnnotation> fiTypes = queryType(interaction, name1, name2);
            rtnTypes.addAll(fiTypes);
        }
        FIAnnotation rtn = mergeTypes(rtnTypes);
        // To get the score for it based on interactions. Extracted FIs have score 1.0.
        // It is possible that several interactions have been mixed together. 
        // See: ReactomeR3CytoscacpePlugin.displayInteraction
        double score = 0;
        // Find the highest score in this loop.
        for (Interaction i : interactions) {
            if (i.getEvidence() != null) {
                if (i.getEvidence().getProbability() > score)
                    score = i.getEvidence().getProbability();
            }
            else { // Extracted FI
                score = 1.0d;
                break;
            }
        }
        rtn.setScore(score);
        return rtn;
    }
    
    public String generateType(Set<String> types) {
        if (types == null || types.size() == 0)
            return "unknown";
        List<String> list = new ArrayList<String>(types);
        Collections.sort(list);
        StringBuilder builder = new StringBuilder();
        for (Iterator<String> it = list.iterator(); it.hasNext();) {
            builder.append(it.next());
            if (it.hasNext())
                builder.append("; ");
        }
        return builder.toString();
    }
    
    private FIAnnotation mergeTypes(Set<FIAnnotation> types) {
        if (types.size()  == 1) {
            FIAnnotation original = types.iterator().next();
            return original.cloneType();
        }
        FIAnnotation merged = new FIAnnotation();
        Set<String> typeNames = new HashSet<String>();
        Set<String> firstLetter = new HashSet<String>(1);
        Set<String> thirdLetter = new HashSet<String>(1);
        int index = 0;
        for (FIAnnotation type : types) {
            typeNames.add(type.getAnnotation());
            if (type.getDirection() != null && type.getDirection().length() > 1) {
                index = type.getDirection().indexOf("-");
                if (index == 0)
                    thirdLetter.add(type.getDirection().charAt(1) + "");
                else if (index == 1)
                    firstLetter.add(type.getDirection().charAt(0) + "");
            }
        }
        merged.setAnnotation(generateType(typeNames));
        // Need to generate direction
        String direction = ((firstLetter.size() == 1) ? firstLetter.iterator().next() : "") + 
                           "-" + // Always
                           ((thirdLetter.size() == 1) ? thirdLetter.iterator().next() : "");
        merged.setDirection(direction);
        return merged;
    }
    
    private Set<FIAnnotation> queryType(Interaction interaction, 
                                        String name1,
                                        String name2) throws Exception {
        Set<FIAnnotation> fiTypes = new HashSet<FIAnnotation>();
        if (interaction == null)
            return fiTypes;
        if (interaction.getEvidence() != null) {
            fiTypes.add(nameToType.get("predicted"));
            return fiTypes;
        }
        // Need to query the Reactome Source
        Set<ReactomeSource> sources = interaction.getReactomeSources();
        boolean needReverse = false;
        // Need to see if a reverse type should be used.
        if (interaction.getFirstProtein().getShortName().equals(name2))
            needReverse = true;
        FIAnnotation fiType = null;
        for (ReactomeSource src : sources) {
            // Check if the annotation should be limited to Reactome only
            if (useReactomeSourceOnly && !src.getDataSource().equals("Reactome"))
                continue;
            if (src.getSourceType() == ReactomeSourceType.COMPLEX)
                fiType = nameToType.get("complex");
            else if (src.getSourceType() == ReactomeSourceType.INTERACTION) {
                fiType = extractInteractionType(src, 
                                                interaction,
                                                needReverse);
            }
            else if (src.getSourceType() == ReactomeSourceType.REACTION) {
                String type = extractTypeFromReaction(src, interaction);
                fiType = nameToType.get(type);
                if (fiType != null && needReverse)
                    fiType = fiType.generateReverseType();
            }
            else if (src.getSourceType() == ReactomeSourceType.TARGETED_INTERACTION) {
                String type = extractTypeFromTargetedInteraction(src, interaction);
                fiType = nameToType.get(type);
                if (fiType != null && needReverse)
                    fiType = fiType.generateReverseType();
            }
            if (fiType == null) {
                logger.error(name1 + " " + name2 + " has no type!" + " Interaction: " + interaction.getDbId());
                fiType = nameToType.get("unknown"); // Used as an error mark
            }
            fiTypes.add(fiType);
        }
        return fiTypes;
    }
    
    /**
     * This method is used to extract type based on TargetedInteraction.
     * @param src
     * @return
     * @throws Exception
     */
    private String extractTypeFromTargetedInteraction(ReactomeSource src,
                                                      Interaction interaction) throws Exception {
//        if (src.getDataSource().equals("TRED")) { // Right now covers TRED TF/Target interactions only.
        // First check based on ids
            Long dbId = src.getReactomeId();
            GKInstance instance = dba.fetchInstance(dbId);
            // Get factor
            GKInstance factor = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.factor);
            Set<GKInstance> refPepSeqs = InstanceUtilities.grepRefPepSeqsFromPhysicalEntity(factor);
            Set<String> factorIds = new HashSet<String>();
            Set<String> factorNames = new HashSet<String>();
            for (GKInstance refPepSeq : refPepSeqs) {
                String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                factorIds.add(id);
                if (refPepSeq.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
                    String geneName = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.geneName);
                    if (geneName != null)
                        factorNames.add(geneName);
                }
            }
            // Get target
            GKInstance target = (GKInstance) instance.getAttributeValue(ReactomeJavaConstants.target);
            refPepSeqs = InstanceUtilities.grepRefPepSeqsFromPhysicalEntity(target);
            Set<String> targetIds = new HashSet<String>();
            Set<String> targetNames = new HashSet<String>();
            for (GKInstance refPepSeq : refPepSeqs) {
                String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                targetIds.add(id);
                if (refPepSeq.getSchemClass().isValidAttribute(ReactomeJavaConstants.geneName)) {
                    String geneName = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.geneName);
                    if (geneName != null)
                        targetNames.add(geneName);
                }
            }
            String proteinId1 = interaction.getFirstProtein().getPrimaryAccession();
            String protein1 = interaction.getFirstProtein().getShortName();
            String proteinId2 = interaction.getSecondProtein().getPrimaryAccession();
            String protein2 = interaction.getSecondProtein().getShortName();
            if (factorIds.contains(proteinId1) && targetIds.contains(proteinId2)) {
                return "expression regulates";
            }
            if (factorIds.contains(proteinId2) && targetIds.contains(proteinId1)) {
                return "expression regulated by";
            }
        // If we cannot find mapping based on ids, we may try using gene names since all TF/target interactions
        // are loaded based on gene names and may choose a ids that are not the same as in the database.
            if (factorNames.contains(protein1) && targetNames.contains(protein2)) {
                return "expression regulates";
            }
            if (factorNames.contains(protein2) && targetNames.contains(protein1)) {
                return "expression regulated by";
            }
//        }
        return null;
    }
    
    /**
     * Extract the FI type from a reaction. There are four types for a FI extracted from a
     * reaction: Input, Catalyze (Or Catalyzed), Inhibit (Inhibited), Activate (Activated).
     * @param src
     * @param interaction
     * @return
     * @throws Exception
     */
    private String extractTypeFromReaction(ReactomeSource src,
                                           Interaction interaction) throws Exception {
        String firstProteinId = interaction.getFirstProtein().getPrimaryAccession();
        String secondProteinId = interaction.getSecondProtein().getPrimaryAccession();
        // Need to map back to the reaction
        GKInstance reaction = dba.fetchInstance(src.getReactomeId());
        // List types first
        List list = reaction.getAttributeValuesList(ReactomeJavaConstants.input);
        Set<GKInstance> inputs = new HashSet<GKInstance>();
        if (list != null) {
            for (Object obj : list) {
                GKInstance input = (GKInstance) obj;
                inputs.add(input);
            }
        }
        Set<String> inputIds = extractIds(inputs);
        Set<GKInstance> catalysts = getCatalysts(reaction);
        Set<String> catalystIds = extractIds(catalysts);
        // Check regulators
        Set<GKInstance> activators = getRegulators(reaction, 
                                                   ReactomeJavaConstants.PositiveRegulation);
        Set<String> activatorIds = extractIds(activators);
        Set<GKInstance> inhibitors = getRegulators(reaction, 
                                                   ReactomeJavaConstants.NegativeRegulation);
        Set<String> inhibitorIds = extractIds(inhibitors);
        // Generate the types
        if (activatorIds.contains(firstProteinId) && !activatorIds.contains(secondProteinId)) {
            return "activate";
        }
        if (inhibitorIds.contains(firstProteinId) && !inhibitorIds.contains(secondProteinId))
            return "inhibit";
        if (catalystIds.contains(firstProteinId) && inputIds.contains(secondProteinId))
            return "catalyze";
        // Check the other direction
        if (activatorIds.contains(secondProteinId) && !activatorIds.contains(firstProteinId))
            return "activated by";
        if (inhibitorIds.contains(secondProteinId) && !inhibitorIds.contains(firstProteinId))
            return "inhibited by";
        if (catalystIds.contains(secondProteinId) && inputIds.contains(firstProteinId))
            return "catalyzed by";
        if (inputIds.contains(firstProteinId) && inputIds.contains(secondProteinId))
            return "input";
        return "reaction"; // Cannot see the difference
    }
    
    /**
     * Extract ids from the set of GKInstances.
     */
    private Set<String> extractIds(Set<GKInstance> instances) throws Exception {
        Set<String> ids = new HashSet<String>();
        for (GKInstance instance : instances) {
            Set<GKInstance> refPepSeqs = InstanceUtilities.grepRefPepSeqsFromPhysicalEntity(instance);
            for (GKInstance refPepSeq : refPepSeqs) {
                String id = (String) refPepSeq.getAttributeValue(ReactomeJavaConstants.identifier);
                if (id != null)
                    ids.add(id);
            }
        }
        return ids;
    }
    
    private Set<GKInstance> getCatalysts(GKInstance reaction) throws Exception {
        Set<GKInstance> catalysts = new HashSet<GKInstance>();
        List cas = reaction.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
        if (cas != null) {
            for (Iterator it1 = cas.iterator(); it1.hasNext();) {
                GKInstance ca = (GKInstance) it1.next();
                List list = ca.getAttributeValuesList(ReactomeJavaConstants.physicalEntity);
                if (list != null)
                    catalysts.addAll(list);
            }
        }
        return catalysts;
    }
    
    private Set<GKInstance> getRegulators(GKInstance reaction,
                                          String clsName) throws Exception {
        Set<GKInstance> regulators = new HashSet<GKInstance>();
        Collection regulations = reaction.getReferers(ReactomeJavaConstants.regulatedEntity);
        if (regulations != null) {
            for (Iterator it1 = regulations.iterator(); it1.hasNext();) {
                GKInstance regulation = (GKInstance) it1.next();
                if (regulation.getSchemClass().isa(clsName)) {
                    List list = regulation.getAttributeValuesList(ReactomeJavaConstants.regulator);
                    for (Iterator it2 = list.iterator(); it2.hasNext();) {
                        GKInstance tmp = (GKInstance) it2.next();
                        if (tmp.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity))
                            regulators.add(tmp);
                    }
                }
            }
        }
        return regulators;
    }
    
    /**
     * A helper method to extract interaction type from a ReactomeSource object.
     * @param src
     * @return
     * @throws Exception
     */
    private FIAnnotation extractInteractionType(ReactomeSource src,
                                          Interaction interaction,
                                          boolean needReverse) throws Exception {
        GKInstance instance = dba.fetchInstance(src.getReactomeId());
        String type = (String) instance.getAttributeValue(ReactomeJavaConstants.interactionType);
        if (type == null)
            return nameToType.get("interaction");
        // Do a very simple reverse mapping
        List<?> interactors = instance.getAttributeValuesList(ReactomeJavaConstants.interactor);
        // Check the first type only
        GKInstance firstValue = (GKInstance) interactors.get(0);
        Set<GKInstance> tmp = new HashSet<GKInstance>();
        tmp.add(firstValue);
        Set<String> firstIdentifiers = extractIds(tmp);
        if (needReverse) {
            if (firstIdentifiers.contains(interaction.getFirstProtein().getPrimaryAccession())) {
                FIAnnotation rtn = nameToType.get(type);
                return rtn.generateReverseType();
            }
            if (firstIdentifiers.contains(interaction.getSecondProtein().getPrimaryAccession())) {
                FIAnnotation rtn = nameToType.get(type);
                return rtn;
            }
        }
        else {
            if (firstIdentifiers.contains(interaction.getFirstProtein().getPrimaryAccession())) {
                FIAnnotation rtn = nameToType.get(type);
                return rtn;
            }
            if (firstIdentifiers.contains(interaction.getSecondProtein().getPrimaryAccession())) {
                FIAnnotation rtn = nameToType.get(type);
                return rtn.generateReverseType();
            }
        }
        return nameToType.get("interaction"); // Too complicated!!!
    }
    
    private String extractInteractionTypeInString(ReactomeSource src) throws Exception {
        GKInstance instance = dba.fetchInstance(src.getReactomeId());
        String type = (String) instance.getAttributeValue(ReactomeJavaConstants.interactionType);
        if (type == null)
            return "interaction";
        return type;
    }
    
    /**
     * Load FI types from a pre-generated file and add some types known already (e.g. predicted).
     * @return
     * @throws Exception
     */
    private void initInteractionTypes() throws Exception {
        nameToType= new HashMap<String, FIAnnotation>();
        reverseNameToName = new HashMap<String, String>();
        SAXBuilder builder = new SAXBuilder();
        Document document = builder.build("resources/InteractionTypeMapper.xml");
        List list = document.getDocument().getRootElement().getChildren("type");
        for (Iterator it = list.iterator(); it.hasNext();) {
            Element elm = (Element) it.next();
            String name = elm.getAttributeValue("name");
            String reverse = elm.getAttributeValue("reverse");
            String direction = elm.getAttributeValue("direction");
            FIAnnotation type = new FIAnnotation();
            type.setAnnotation(name);
            type.setDirection(direction);
            type.setReverseAnnotation(reverse);
            nameToType.put(name, type);
            reverseNameToName.put(reverse, name);
        }
        // Some specific type
        FIAnnotation predictedType = new FIAnnotation();
        predictedType.setAnnotation("predicted");
        predictedType.setReverseAnnotation("predicted");
        predictedType.setDirection("-");
        nameToType.put(predictedType.getAnnotation(), predictedType);
        // Some specific type
        FIAnnotation unknownType = new FIAnnotation();
        unknownType.setAnnotation("unknown");
        nameToType.put(unknownType.getAnnotation(), unknownType);
        loadReactionTypes(nameToType);
    }
    
    private void loadReactionTypes(Map<String, FIAnnotation> nameToType) {
        String[] types = new String[] {
                "complex",
                "activate",
                "inhibit",
                "catalyze",
                "activated by",
                "inhibited by",
                "catalyzed by",
                "input",
                "reaction", // Cannot see the difference
                // Following for TF/Target interactions
                "expression regulates",
                "expression regulated by",
                "interaction"
        };
        String[] reverseTypes = new String[] {
                "complex",
                "activated by",
                "inhibited by",
                "catalyzed by",
                "activate",
                "inhibite",
                "catalyze",
                "input",
                "reaction", // Cannot see the difference
                "expression regulated by",
                "expression regulates",
                "interaction"
        };
        String[] directions = new String[] {
                "-",
                "->",
                "-|",
                "->",
                "<-",
                "|-",
                "<-",
                "-",
                "-",
                "->",
                "<-",
                "-"
        };
        for (int i = 0; i < types.length; i++) {
            FIAnnotation type = new FIAnnotation();
            type.setAnnotation(types[i]);
            type.setReverseAnnotation(reverseTypes[i]);
            type.setDirection(directions[i]);
            nameToType.put(type.getAnnotation(), type);
        }
    }
    
    /**
     * Used to check a single FI in case it is wrong.
     * @throws Exception
     */
    @Test
    public void checkFIAnnotation() throws Exception {
        FIConfiguration config = FIConfiguration.getConfiguration();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            config.get("REACTOME_SOURCE_DB_NAME"),
                                            config.get("DB_USER"),
                                            config.get("DB_PWD"));
        setSourceDBA(dba);
        // Objects related to hibernate session
        initSession();
        Session session = sessionFactory.openSession();
        Query query = session.createQuery("FROM Interaction as i WHERE i.firstProtein.shortName = ? "
                + "AND i.secondProtein.shortName = ?");
//        String fi = "TRAPPC2P1\tZBTB33";
//        String fi = "BRCA1\tTRAPPC2B";
        // Because HSPA1A and HSAP1B have the same aa sequences, another ENCODE interaction
        // EP300 HSAP1B is merged into the following interaction, which makes an ENCODE annotation
        // not work any more. Needs to manually annotate it.
//        String fi = "EP300\tHSPA1A";
        // A list of FIs that cannot be annotated because of sequence-based merging
        // In 2015 version
        String[] fis = new String[] {
//                "BRCA1 TRAPPC2B", //has no type! Interaction: 51029
                "EP300 HSPA1A", // has no type! Interaction: 36573
                "HSF1 HSPA1A", // has no type! Interaction: 131579
//                "HSPA1A JUN", // has no type! Interaction: 17546
//                "HSPA1A PPARGC1A", // has no type! Interaction: 93702
//                "TRAPPC2B ZBTB33" // has no type! Interaction: 158201
        };
        for (String fi : fis) {
            String[] genes = fi.split(" ");
            FIAnnotation fiAnnot = annotate(genes[0], 
                                            genes[1],
                                            query);
            String text = String.format("%s\t%s\t%s\t%.2f",
                                        fi,
                                        fiAnnot.getAnnotation(),
                                        fiAnnot.getDirection(),
                                        fiAnnot.getScore());
            System.out.println(text);
        }
        session.close();
    }
    
    /**
     * Use this method to generate a FI network containing FIs from Reactome only.
     * The annotations will be performed for Reactome source only.
     * @throws Exception
     */
    @Test
    public void generateAnnotatedReactomeFIs() throws Exception {
        FIConfiguration config = FIConfiguration.getConfiguration();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            config.get("REACTOME_SOURCE_DB_NAME"),
                                            config.get("DB_USER"),
                                            config.get("DB_PWD"));
        setSourceDBA(dba);
        
        // Auto-generate the FI file name with annotations
        String fiFileName = config.get("GENE_FI_FILE_NAME");
        int index = fiFileName.lastIndexOf(".");
        String outFileName = fiFileName.substring(0, index) + "_Reactome_Annotated" + fiFileName.substring(index);
        // Intiailize hibernate
        initSession();
        Session session = sessionFactory.openSession();
        List<?> annotatedFIs = fetchAnnotatedFIs(session);
        List<Interaction> reactomeFIs = new ArrayList<Interaction>();
        // Filter to Reactome FIs
        // For output
        Set<String> fisInGenes = new HashSet<String>();
        for (Object obj : annotatedFIs) {
            Interaction fi = (Interaction) obj;
            Set<ReactomeSource> sources = fi.getReactomeSources();
            for (ReactomeSource src : sources) {
                if (src.getDataSource().equals("Reactome")) {
                    reactomeFIs.add(fi);
                    String fiInGene = getFIInName(fi);
                    if (fiInGene != null)
                        fisInGenes.add(fiInGene);
                    break;
                }
            }
        }
        logger.info("Total Reactome FIs: " + reactomeFIs.size());
        useReactomeSourceOnly = true;
        boolean hasUnknown = annotate(reactomeFIs, outFileName);
        session.close();
        if (hasUnknown) {
            // Signaling to stop for manual fix
            throw new IllegalStateException("Unknown annotation has been found. Need a manual fix. Check errors in logging!");
        }
        // Generate a FI file without annotation
        outFileName = fiFileName.substring(0, index) + "_Reactome" + fiFileName.substring(index);
        List<String> fiList = new ArrayList<String>(fisInGenes);
        Collections.sort(fiList);
        FileUtility fu = new FileUtility();
        fu.saveCollection(fiList, outFileName);
    }
    
    /**
     * This method is used to annotate FIs so that they can be placed in the download
     * folder and used for quick performance.
     * @throws Exception
     */
    @Test
    public void annoateAllFIs() throws Exception {
        // Parameters for the 2009 version
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_28_plus_i_myisam",
//                                            "root",
//                                            "macmysql01");
//        String fiFileName = "WebContent/WEB-INF/FIsInGene_041709.txt";
//        String outFileName = "WebContent/WEB-INF/FIsInGene_041709_with_annotations.txt";
//        File configFile = new File("WebContent/WEB-INF/funcIntHibernate.cfg.xml");
        
        // Parameters for the 2012 version
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_41_plus_i",
//                                            "root",
//                                            "macmysql01");
//        String dirName = "/Users/gwu/Documents/EclipseWorkspace/FINetworkBuild/results/2012/";
//        String fiFileName = dirName + "FIsInGene_071012.txt";
//        String outFileName = dirName + "FIsInGene_071012_with_annotations.txt";
//        
//        File configFile = new File("WebContent/WEB-INF/funcIntHibernate.cfg.xml");
        
        // Parameters for the 2013 version
//        MySQLAdaptor dba = new MySQLAdaptor("localhost",
//                                            "reactome_47_plus_i",
//                                            "root",
//                                            "macmysql01");
//        String dirName = "/Users/gwu/Documents/EclipseWorkspace/FINetworkBuild/results/2013/";
//        String fiFileName = dirName + "FIsInGene_121013.txt";
//        String outFileName = dirName + "FIsInGene_121013_with_annotations.txt";
//        
//        File hibernateConfig = new File("WebContent/WEB-INF/funcIntHibernate.cfg.xml");
        
        // The following running are based on configuration in the FINetworkBuild project, which
        // should work for all versions.
        FIConfiguration config = FIConfiguration.getConfiguration();
        MySQLAdaptor dba = new MySQLAdaptor("localhost",
                                            config.get("REACTOME_SOURCE_DB_NAME"),
                                            config.get("DB_USER"),
                                            config.get("DB_PWD"));
        setSourceDBA(dba);
        
        // Auto-generate the FI file name with annotations
        String fiFileName = config.get("GENE_FI_FILE_NAME");
        int index = fiFileName.lastIndexOf(".");
        String outFileName = fiFileName.substring(0, index) + "_with_annotations" + fiFileName.substring(index);
        
        // Intiailize hibernate
        initSession();
        Session session = sessionFactory.openSession();
        // This is a special case occurred in 2014 version of the FI network: because of the same amino acid 
        // sequence of TRAPPC2 and TRAPPC2P1, the following FI cannot be mapped and a manual annotation has 
        // to be performed. 
        // The following FI doesn't appear in 2015. So the following code is commented out! However, there is another
        // FI, BRCA1\tTRAPPC2B, which is assumed to have the same problem, and need a manual annotation
        // after the annotation and add the annotation back to the annotation file.
//        String fi = "TRAPPC2P1\tZBTB33";
//        String[] genes = fi.split("\t");
//        FIAnnotation fiAnnot = annotate(genes[0], genes[1]);
//        String text = String.format("%s\t%s\t%s\t%.2f",
//                                       fi,
//                                       fiAnnot.getAnnotation(),
//                                       fiAnnot.getDirection(),
//                                       fiAnnot.getScore());
//        System.out.println(text);
//        if (true)
//            return;
        // We will start from querying FIs for faster performance
        List<?> predictedFIs = fetchPredictedFIs(session);
        List<?> annotatedFIs = fetchAnnotatedFIs(session);
        List allFIs = new ArrayList(predictedFIs);
        allFIs.addAll(annotatedFIs);
        boolean hasUnknown = annotate(allFIs, outFileName);
        session.close();
        if (hasUnknown) {
            // Signaling to stop for manual fix
            throw new IllegalStateException("Unknown annotation has been found. Need a manual fix. Check errors in logging!");
        }
    }
    
    private boolean annotate(List<?> allFIs,
                             String outFileName) throws IOException, Exception {
        Map<String, List<Interaction>> fiInGeneToFIs = new HashMap<String, List<Interaction>>();
        for (Object obj : allFIs) {
            Interaction i = (Interaction) obj;
            String fiInGene = getFIInName(i);
            if (fiInGene == null)
                continue; // There are many cases like this. See comments in other places calling this method.
            List<Interaction> list = fiInGeneToFIs.get(fiInGene);
            if (list == null) {
                list = new ArrayList<Interaction>();
                fiInGeneToFIs.put(fiInGene, list);
            }
            list.add(i);
        }
        logger.info("Total FIs for annotations: " + fiInGeneToFIs.size());
        // Sort it for easy comparison with other files
        List<String> fiInGeneList = new ArrayList<String>(fiInGeneToFIs.keySet());
        Collections.sort(fiInGeneList);
        
        FileUtility fu = new FileUtility();
        fu.setOutput(outFileName);
        fu.printLine("Gene1\tGene2\tAnnotation\tDirection\tScore");
        int count = 0;
        boolean hasUnknown = false;
        long time1 = System.currentTimeMillis();
        for (String fiInGene : fiInGeneList) {
//            if (!(line.contains("ATF2") && line.contains("CDKN1A")))
//                continue;
            List<Interaction> interactions = fiInGeneToFIs.get(fiInGene);
            String[] tokens = fiInGene.split("\t");
            FIAnnotation annotation = annotate(tokens[0],
                                               tokens[1],
                                               interactions);
            if (annotation.getAnnotation().equals("unknown"))
                hasUnknown = true;
            String outLine = String.format("%s\t%s\t%s\t%.2f",
                                           fiInGene,
                                           annotation.getAnnotation(),
                                           annotation.getDirection(),
                                           annotation.getScore());
            fu.printLine(outLine);
            count ++;
//            if (count == 20000)
//                break;
//            if (count % 5000 == 0)
//                logger.info("Total annotated FIs: " + count);
        }
        fu.close();
        logger.info("Total interactions: " + count);
        long time2 = System.currentTimeMillis();
        logger.info("Total time: " + (time2 - time1));
        return hasUnknown;
    }
}
